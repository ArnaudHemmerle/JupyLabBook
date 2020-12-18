# custom libraries
from lib.extraction.common import PyNexus as PN

import numpy as np
import lmfit as lm
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import time
import subprocess
import sys

def Treat(nxs_filename, recording_dir,
          channel0, thetazfactor, wavelength, thetac,
          binsize, computeqz, 
          absorbers='', logx=False, logy=False, logz=False, nblevels=50, cmap='jet',
          working_dir='', moytocreate=[10, 20, 40],
          show_data_stamps=False, plot=False, save=False, verbose=False):
    
    '''
    Call functions for extracting, plotting, and saving a GIXD scan.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    channel0 : float
        vertical channel corresponding to the Vineyard's peak
    thetazfactor : float
        factor for conversion from channel to radian in the vertical direction (rad/channel)
    wavelength : float
        wavelength in nm
    thetac : float
        critical angle in rad
    binsize : int
        size in pixels of the vertical binning (along qz)
    computeqz : bool
        switch from pixels to qz in the vertical direction
    absorbers : str, optional
        text to display indicating which absorber was used
    logx : bool, optional
        log on the x axis of the integrated profile
    logy : bool, optional
        log on the y axis of the integrated profile
    logz : bool, optional
        log on the image
    nblevels : int, optional
        number of color levels for the image display
    cmap : str, optional
        colormap of the image
    working_dir : str, optional
        directory where the treated files will be stored
    moytocreate : array_like, optional
        binsizes to be saved
    show_data_stamps : bool, optional
        print the list of sensors from the nexus file
    plot : bool, optional
        plot the GIXD
    save : bool, optional
        save the GIXD
    verbose : bool, optional
        verbose mode


    Returns
    -------
    array_like
        x, an array containing either qxy (nm^-1) values, if qxy is available in the list of sensors,
        or delta (deg) if not
    array_like
        y, an array containing either qz (nm^-1) values, if computeqz is True,
        or vertical channels if not        
    array_like
        daty, rods integrated over the whole vertical axis of the detector            
    array_like
        datyTop, rods integrated over the top half vertical axis of the detector            
    array_like
        datyBottom, rods integrated over the bottom half vertical axis of the detector     
    array_like
        datyFirstQuarter, rods integrated over the bottom quarter vertical axis of the detector                        
    array_like
        mat, each line corresponds to a position of delta.
        It is the matrix corresponding to the image displayed in plot            
    array_like
        mat_binned, original matrix after binning            
    array_like
        ch_binned, array of channels after binning            
    float or None
        mean_pi, the average of the surface pressure (mN/m) over the scan (None if pressure not found)            
    float or None
        mean_area, the average of the area per molecule (nm^2) over the scan (None if area not found)        
    float or None
        mean_gamma, the average of gamma (deg) over the scan (None if gamma not found)


    Raises
    ------
    SystemExit('Pilatus not found')
        when Pilatus is not found
    SystemExit('No sensor found')
        when delta and qxy are not found
    SystemExit('gamma not found')
        when gamma not found and computeqz is True
    '''

    x, y, xlabel, ylabel, column_x,\
    daty, datyTop, datyBottom, datyFirstQuarter,\
    mat, mat_binned, ch_binned, \
    mean_pi, mean_area, mean_gamma = \
    Extract(nxs_filename, recording_dir,
            channel0, thetazfactor, wavelength, thetac,
            binsize, computeqz,
            show_data_stamps, verbose)

    if plot:
        Plot(x, y, xlabel, ylabel, 
             daty, datyTop, datyBottom, datyFirstQuarter,
             mat_binned, ch_binned,
             mean_pi, mean_area, mean_gamma,
             nxs_filename, absorbers, logx, logy, logz, nblevels, cmap)

    if save:
        Save(x, daty, datyTop, datyBottom, datyFirstQuarter, mat, moytocreate, mean_gamma, column_x,
             channel0, thetazfactor, wavelength, thetac,
             nxs_filename, recording_dir, working_dir,
             computeqz, verbose)


    return x, y,\
           daty, datyTop, datyBottom, datyFirstQuarter,\
           mat, mat_binned, ch_binned, \
           mean_pi, mean_area, mean_gamma


def Extract(nxs_filename, recording_dir,
            channel0, thetazfactor, wavelength, thetac,
            binsize, computeqz,
            show_data_stamps, verbose):

    '''
    Extract the nexus scan and return useful quantities for GIXD.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    channel0 : float
        vertical channel corresponding to the Vineyard's peak
    thetazfactor : float
        factor for conversion from channel to radian in the vertical direction (rad/channel)
    wavelength : float
        wavelength in nm
    thetac : float
        critical angle in rad
    binsize : int
        size in pixels of the vertical binning (along qz)
    computeqz : bool
        switch from pixels to qz in the vertical direction
    show_data_stamps : bool
        print the list of sensors from the nexus file
    verbose : bool
        verbose mode


    Returns
    -------
    array_like
        x, an array containing either qxy (nm^-1) values, if qxy is available in the list of sensors,
        or delta (deg) if not
    array_like
        y, an array containing either qz (nm^-1) values, if computeqz is True,
        or vertical channels if not        
    str
        xlabel, indicating whether x corresponds to qxy or delta (useful for plot)            
    str
        ylabel, indicating whether y corresponds to qz or channels (useful for plot)           
    int
        column_x, the column corresponding to the x values in stamps0D (useful for save)            
    array_like
        daty, rods integrated over the whole vertical axis of the detector            
    array_like
        datyTop, rods integrated over the top half vertical axis of the detector            
    array_like
        datyBottom, rods integrated over the bottom half vertical axis of the detector     
    array_like
        datyFirstQuarter, rods integrated over the bottom quarter vertical axis of the detector                        
    array_like
        mat, each line corresponds to a position of delta.
        It is the matrix corresponding to the image displayed in plot            
    array_like
        mat_binned, original matrix after binning            
    array_like
        ch_binned, array of channels after binning            
    float or None
        mean_pi, the average of the surface pressure (mN/m) over the scan (None if pressure not found)            
    float or None
        mean_area, the average of the area per molecule (nm^2) over the scan (None if area not found)        
    float or None
        mean_gamma, the average of gamma (deg) over the scan (None if gamma not found)


    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    SystemExit('Pilatus not found')
        when Pilatus is not found
    SystemExit('No sensor found')
        when delta and qxy are not found
    SystemExit('gamma not found')
        when gamma not found and computeqz is True
    '''
        
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
        sys.exit('Nexus not found')
        
    else:
        
        # Extract sensors from the Pilatus
        column_z=None
        column_qxy=None
        column_delta=None
        column_pi=None
        column_area=None
        column_gamma=None
        
        if verbose: print(PN._BLUE+" - Open Nexus Data File :"+ PN._RESET)
        if verbose: print('\t'+nxs_path)
        try:
            nexus=PN.PyNexusFile(nxs_path, fast=True)
        except:
            print(PN._RED,'\t Nexus file seems not to exist or is not correct',PN._RESET)
            sys.exit('Nexus not found')
        
        nbpts=np.int(nexus.get_nbpts())
        if verbose: print("\t. Number of data points: ", nbpts)
        
        # Get stamps
        stamps=nexus.extractStamps()
        if show_data_stamps : print("\t. Available Counters:")
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                if show_data_stamps : print("\t\t", i, ' -------> ', stamps[i][1])
                if stamps[i][1].lower()=='pilatus':
                    column_z=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])

                    
        # Extract 0D data
        stamps0D, data=nexus.extractData('0D')
        
        # Find columns
        for i in range(len(stamps0D)):
            if stamps0D[i][1] is not None:
                if stamps0D[i][1].lower()=='delta':
                    column_delta=i
                if stamps0D[i][1].lower()=='qxy':
                    column_qxy=i
                if stamps0D[i][1].lower()=='surfacepressure':
                    column_pi=i
                if stamps0D[i][1].lower()=='areapermolecule':
                    column_area=i
                if stamps0D[i][1].lower()=='gamma':
                    column_gamma=i
                    
        
        # Check that Pilatus data are present (images)
        if column_z is not None:
            if verbose: print('\t. Pilatus data found, (column %d, alias %s)'%(column_z, stamps[column_z][1]))
        else:
            print(PN._RED,'\t. No pilatus data found', PN._RESET)
            nexus.close()
            sys.exit('Pilatus not found')

        # Check what is the x input (qxy or delta)
        if column_qxy is None:
            if column_delta is None:
                print(PN._RED,'\t. No usual actuator for GIXD found, stop here', PN._RESET)
                sys.exit('No sensor found')
            else:
                column_x=column_delta
                if verbose: print('\t. delta data found, (column %d, alias %s)'%(column_x, stamps[column_x][1]))
        else:
            column_x=column_qxy
            if verbose: print('\t. qxy data found, (column %d, alias %s)'%(column_x, stamps[column_x][1]))

        # Find start and stop without Nan
        istart=0
        istop=nbpts
        i=0
        while np.isnan(data[0][i]):
            istart=i
            i=i+1
        i=nbpts-1
        while np.isnan(data[0][i]):
            istop=i+1
            i=i-1
        istop=istop-1
        if verbose: print('\t. Valid data between points %d and %d'%(istart, istop))

        # Compute and print mean values
        if column_pi is not None:
            mean_pi=data[column_pi][istart:istop-1].mean()
            if verbose: print('\t. Surface pressure data found, mean value %3.4g ± %3.4g mN/m'
                              %(mean_pi,data[column_pi][istart:istop-1].std() ))
        else:
            print('\t. No surface pressure data found')
            mean_pi = None
            
        if column_area is not None:
            mean_area=data[column_area][istart:istop-1].mean()
            if verbose: print('\t. Area per molecule data found, mean value %3.4g ± %3.4g nm2 per molecule'
                              %(mean_area, data[column_area][istart:istop-1].std()))
        else:
            print('\t. No area per molecule data found')
            mean_area = None
            
        if column_gamma is not None:
            mean_gamma=data[column_gamma][istart:istop-1].mean()
            if verbose: print('\t. Gamma motor data found, mean value %3.4g deg'%(mean_gamma))
        else:
            print('\t. No gamma motor data found')
            mean_gamma = None
            if computeqz:
                print('\t. gamma is required to compute qz. Add gamma to the sensors.')
                nexus.close()
                sys.exit('gamma not found')

        # Load images
        stamps, images=nexus.extractData('2D')        
        nexus.close()
        
        # Get positions of the dead pixels
        try:
            dead_pixels = np.genfromtxt('lib/extraction/common/dead_pixels.dat', dtype = 'uint16', delimiter = ', ')
        except:
            print('Careful: the file lib/extraction/common/dead_pixels.dat was not found. Taking no dead pixels.')
            dead_pixels = []
            
        
        daty=[]
        datyTop=[]
        datyBottom=[]
        datyFirstQuarter=[]
        mat=[]

        for i in range(istart, istop, 1):
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i+1, nbpts))
            sys.stdout.flush()
            
            image=images[0][i]
            
            for i,j in dead_pixels:
                image[i,j]=0.0
                
            # Keep only the ROI corresponding to the Soller
            ROI=[510, 350, 130, 692]
            image=image[ROI[1]:ROI[1]+ROI[3], ROI[0]:ROI[0]+ROI[2]]
            
            # Replace the intensity of the dead zones with a value of 0
            image=np.where(image<0., 0., image)
            
            # Rod (integration along the horizontal axis)
            rod=image.sum(axis=1)
            
            # Final image with rods stacked (dim: nb_angle x vertical_dim_ROI)
            mat.append(rod)
            
            # Integrated rods
            daty.append(rod.sum())
            
            # Integrate the rod on different parts of the detector only
            datyTop.append(rod[0:np.int(ROI[3]/2)].sum())
            datyBottom.append(rod[np.int(ROI[3]/2):ROI[3]].sum())
            datyFirstQuarter.append(rod[np.int(3*ROI[3]/4):ROI[3]].sum())
            
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()

        # Convert in numpy array
        daty = np.array(daty)[istart:istop]
        datyTop = np.array(datyTop)[istart:istop]
        datyBottom = np.array(datyBottom)[istart:istop]
        datyFirstQuarter = np.array(datyFirstQuarter)[istart:istop]
        mat = np.array(mat)            
        
        # Extract x values (qxy or delta)
        x = data[column_x]
        x = x[istart:istop]

        # Bin the matrix, return binned vertical channels and binned matrix
        ch_binned, mat_binned = Groupe(mat, binsize=binsize)

        # Extract y values (qz or vertical channels)
        if computeqz==True:
            # Compute and return qz
            thetaz=thetac+(mean_gamma*np.pi/180.0)+(channel0-ch_binned)*thetazfactor
            y=2.0*np.pi*np.sin(thetaz)/wavelength
        else:
            # Return the vertical channels 
            y=ch_binned
            
        # Nature and unit of the axis
        if column_qxy is None:
            # x axis corresponds to delta
            xlabel = str(stamps0D[column_x][1])
        else:
            # x axis corresponds to qxy
            xlabel = str(stamps0D[column_x][1]+' ('+stamps0D[column_x][2]+')')
        
        if computeqz:
            ylabel = 'qz (nm-1)'
        else:
            ylabel = r'$vertical\ channels$'
            
        return x, y, xlabel, ylabel, column_x,\
               daty, datyTop, datyBottom, datyFirstQuarter,\
               mat, mat_binned, ch_binned, \
               mean_pi, mean_area, mean_gamma
    
def Plot(x, y, xlabel, ylabel, 
         daty, datyTop, datyBottom, datyFirstQuarter,
         mat_binned, ch_binned,
         mean_pi, mean_area, mean_gamma,
         nxs_filename, absorbers, logx, logy, logz, nblevels, cmap):
    '''
    Plot GIXD data.

    Parameters
    ----------
    x : array_like
        either qxy (nm^-1) or delta (deg) values
    y : array_like
        either qz (nm^-1) values or vertical channels               
    xlabel : str
        label for qxy or delta            
    ylabel : str
        label for qz or channels
    daty : array_like
        rods integrated over the whole vertical axis of the detector            
    datyTop : array_like
        rods integrated over the top half vertical axis of the detector            
    datyBottom : array_like
        rods integrated over the bottom half vertical axis of the detector     
    datyFirstQuarter : array_like
        rods integrated over the bottom quarter vertical axis of the detector                        
    mat_binned : array_like
        binned matrix corresponding to the GIXD image to be displayed          
    ch_binned : array_like
        array of channels after binning            
    mean_pi : float or None
        the average of the surface pressure (mN/m) over the scan      
    mean_area : float or None
        the average of the area per molecule (nm^2) over the scan       
    mean_gamma : float or None
        the average of gamma (deg) over the scan       
    nxs_filename : str
        nexus filename        
    absorbers : str, optional
        text to display indicating which absorber was used
    logx : bool, optional
        log on the x axis of the integrated profile
    logy : bool, optional
        log on the y axis of the integrated profile
    logz : bool, optional
        log on the image
    nblevels : int, optional
        number of color levels for the image display
    cmap : str, optional
        colormap of the image
    '''

    # Print absorbers
    if absorbers != '':
        print("\t. Absorbers:", str(absorbers))

    # Error bars along y
    error=np.sqrt(daty)

    # Create Graph  
    fig=plt.figure(nxs_filename, figsize=(12,5))
    fig.subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.16)
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)

    # Plot the integrated spectra

    if (not logx) and (not logy):
        ax1.errorbar(x, daty, yerr=error, fmt='ro', label="full")
        ax1.plot(x, datyTop, 'k-', label='top')
        ax1.plot(x, datyBottom, 'b-', label='bottom')
        ax1.plot(x, datyFirstQuarter, 'r-', label='bottom quarter')
        ax1.legend()
    elif logx and (not logy):
        ax1.semilogx(x, daty, 'ro')
    elif (not logx) and (logy):
        ax1.semilogy(x, daty, 'ro')
    elif logx and (logy):
        ax1.loglog(x, daty, 'ro')
    ax1.set_xlabel(xlabel, labelpad=13, fontsize='large')
    ax1.set_ylabel("qz integrated intensity", labelpad=13, fontsize='large')

    # Plot the matrix
    if logz:
        ax2.contourf(x, y, np.log(mat_binned.transpose()))
    else:
        zmax=mat_binned.max()
        zmin=mat_binned.min()
        ax2.contourf(x, y, (mat_binned.transpose()),levels=np.linspace(zmin, zmax, nblevels), cmap=cmap)

    ax2.set_ylabel(ylabel, fontsize='large')
    ax2.set_xlabel(xlabel, labelpad=13, fontsize='large')

    if mean_pi is not None:
        fig.text(.04, .05, r'$\pi = %3.4gmN.m^{-1}$'%(mean_pi),
                 fontsize='large', color='red')
    if mean_gamma is not None:
        fig.text(.96, .05, r'$\gamma = %3.4g deg$'%(mean_gamma),
                 fontsize='large', color='red', horizontalalignment='right')

    S_suptitle=nxs_filename[nxs_filename.rfind('/')+1:]
    fig.suptitle(S_suptitle, fontsize='x-large')
    
    plt.show()

        
def Save(x, daty, datyTop, datyBottom, datyFirstQuarter, mat, moytocreate, mean_gamma, column_x,
         channel0, thetazfactor, wavelength, thetac,
         nxs_filename, recording_dir, working_dir,
         computeqz, verbose):   
    '''
    Save GIXD data.
    XXX_1D.mat: each line corresponds to a position of delta.
                It is the matrix corresponding to the image displayed.
    XXX_1D.dat: each line corresponds to a position of delta.
                It contains the value of each sensor, and integration along qz:
                - QzIntegrated : over the whole detector,
                - QzIntegratedTop : over its top half,
                - QzIntegratedBottom : over its bottom half,
                - QzIntegratedBottomQuarter : over its bottom quarter.
    Binned data:
        - XXX_1D.matNN : binning of the matrix, with NN the number of points per bin.
        - XXX_1D_qz.datNN : to convert bin number to qz in XXX_1D.matNN.
        - XXX_1D.moyNN : a more convenient way to represent the binned matrices
                         with a 3 columns (qxy, qz, intensity) display.

    Parameters
    ----------
    x : array_like
        either qxy (nm^-1) or delta (deg) values
    daty : array_like
        rods integrated over the whole vertical axis of the detector            
    datyTop : array_like
        rods integrated over the top half vertical axis of the detector            
    datyBottom : array_like
        rods integrated over the bottom half vertical axis of the detector     
    datyFirstQuarter : array_like
        rods integrated over the bottom quarter vertical axis of the detector                        
    mat : array_like
        Original matrix corresponding to the GIXD image
    moytocreate : array_like, optional
        binsizes to be saved    
    mean_gamma : float or None
        the average of gamma (deg) over the scan     
    column_x : int
        the column corresponding to the x values in stamps0D
    channel0 : float
        vertical channel corresponding to the Vineyard's peak
    thetazfactor : float
        factor for conversion from channel to radian in the vertical direction (rad/channel)
    wavelength : float
        wavelength in nm
    thetac : float
        critical angle in rad  
    nxs_filename : str
        nexus filename        
    recording_dir : str
        directory where the nexus file is stored
    working_dir : str, optional
        directory where the treated files will be stored            
    computeqz : bool
        switch from pixels to qz in the vertical direction
    verbose : bool, optional
        verbose mode            
    '''

    # Extract 0D data
    nxs_path = recording_dir+nxs_filename
    nexus = PN.PyNexusFile(nxs_path)
    stamps0D, data = nexus.extractData('0D')
    nbpts = np.int(nexus.get_nbpts())
    
    # Find start and stop without Nan
    istart=0
    istop=nbpts
    i=0
    while np.isnan(data[0][i]):
        istart=i
        i=i+1
    i=nbpts-1
    while np.isnan(data[0][i]):
        istop=i+1
        i=i-1
    istop=istop-1

    # Create Save Name
    savename=working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]+'_1D'

    # Save the original matrix
    np.savetxt(savename+'.mat', mat)
    if verbose: print('\t. Original, non binned, matrix saved in:')
    if verbose: print("\t", savename+'.mat')

    # Take care of scalar data
    tosave=np.zeros((daty.shape[0], 4), float)
    tosave[:,0]=daty
    tosave[:,1]=datyTop
    tosave[:,2]=datyBottom
    tosave[:,3]=datyFirstQuarter

    # Concatenate scalar data
    data=np.array(data).transpose()
    data=data[istart:istop,:]
    data_tosave=np.concatenate((data,tosave), axis=1)
    
    f=open(savename+'.dat', "w")
    
    # Create and save header
    s=""
    for i in range(len(stamps0D)):
        if stamps0D[i][1] is not None:
            s=s+stamps0D[i][1]+'\t'
        else:
            s=s+stamps0D[i][0]+'\t'
    s=s+'QzIntegrated \t QzIntegratedTop \t QzIntegratedBottom \tQzIntegratedBottomQuarter\n'
    f.write(s)
    
    # Save data
    for i in range(data_tosave.shape[0]):
        s=""
        for j in range(data_tosave.shape[1]):
            s=s+str(data_tosave[i,j])+'\t'
        f.write(s+'\n')
    f.close()
    
    if verbose: print('\t. Scalar data saved in:')
    if verbose: print("\t", savename+'.dat')

    # Save as moy
    for binsize in moytocreate:
        
        # Bin the matrix
        ch_binned, mat_binned= Groupe(mat, binsize=binsize)

        # Extract y values (qz or vertical channels)
        if computeqz:
            # Compute and return qz
            thetaz=thetac+(mean_gamma*np.pi/180.0)+(channel0-ch_binned)*thetazfactor
            y=2.0*np.pi*np.sin(thetaz)/wavelength
        else:
            # Return the vertical channels 
            y=ch_binned
            
        # Create the moy images
        moy=np.zeros((nbpts*y.shape[0],3), float)
        index=0
        for i in range(data.shape[0]):
            for j in range(y.shape[0]):
                moy[index,0]=x[i]
                moy[index,1]=y[j]
                moy[index,2]=mat_binned[i,j]
                index=index+1
                
        f=open(savename+'.moy'+str(binsize), 'w')
        f.write(stamps0D[column_x][1]+'\t'+'Qz \t Intensity\n')    
        for i in range(moy.shape[0]):
            f.write(str(moy[i,0])+'\t'+str(moy[i,1])+'\t'+str(moy[i,2])+'\n')          
        f.close()
        
        # Save the matrix
        np.savetxt(savename+'.mat'+str(binsize), mat_binned)
        
        # Save the Qz
        if computeqz:
            np.savetxt(savename+'_qz.dat'+str(binsize), y)
            if verbose: print('\t. qz values saved in:')
            if verbose: print('\t'+savename+'_qz.dat'+str(binsize))
                
        if verbose:
            print('\t. Binned matrix saved in:')
            print("\t", savename+'.mat'+str(binsize))

            print('\t. XYZ data saved in:')
            print("\t", savename+'.moy'+str(binsize))

    plt.show()

 
def Extract_channel0(nxs_filename='SIRIUS_test.nxs', recording_dir='',
                     binsize=10, logx=False, logy=False, logz=False, nblevels=50, cmap='jet',
                     show_data_stamps = True, verbose=True):
    '''
    Extract the channel of the Vineyard's peak.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    binsize : int, optional
        size in pixels of the vertical binning (along qz)
    logx : bool, optional
        log on the x axis of the integrated profile
    logy : bool, optional
        log on the y axis of the integrated profile
    logz : bool, optional
        log on the image
    nblevels : int, optional
        number of color levels for the image display
    cmap : str, optional
        colormap of the image
    show_data_stamps : bool, optional
        print the list of sensors from the nexus file
    verbose : bool, optional
        verbose mode


    Returns
    -------
    int
        channel0, the vertical channel corresponding to the Vineyard's peak

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file is not found
    SystemExit('Pilatus not found')
        when Pilatus is not found
    SystemExit('No sensor found')
        when delta and qxy are not found
    SystemExit('gamma not found')
        when gamma not found and computeqz is True
    '''

        
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
        sys.exit('Nexus not found')
        
    else:
        
        # Extract sensors from the Pilatus
        column_z=None
        column_qxy=None
        column_delta=None
        column_pi=None
        column_area=None
        column_gamma=None
        
        if verbose: print(PN._BLUE+" - Open Nexus Data File :"+ PN._RESET)
        if verbose: print('\t'+nxs_path)
        try:
            nexus=PN.PyNexusFile(nxs_path, fast=True)
        except:
            print(PN._RED,'\t Nexus file seems not to exist or is not correct',PN._RESET)
            sys.exit('Nexus not found')
        
        nbpts=np.int(nexus.get_nbpts())
        if verbose: print("\t. Number of data points: ", nbpts)
        
        # Get stamps
        stamps=nexus.extractStamps()
        if show_data_stamps : print("\t. Available Counters:")
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                if show_data_stamps : print("\t\t", i, ' -------> ', stamps[i][1])
                if stamps[i][1].lower()=='pilatus':
                    column_z=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])

                    
        # Extract 0D data
        stamps0D, data=nexus.extractData('0D')
        
        # Find columns
        for i in range(len(stamps0D)):
            if stamps0D[i][1] is not None:
                if stamps0D[i][1].lower()=='delta':
                    column_delta=i
                if stamps0D[i][1].lower()=='qxy':
                    column_qxy=i
                if stamps0D[i][1].lower()=='surfacepressure':
                    column_pi=i
                if stamps0D[i][1].lower()=='areapermolecule':
                    column_area=i
                if stamps0D[i][1].lower()=='gamma':
                    column_gamma=i
                    
        
        # Check that Pilatus data are present (images)
        if column_z is not None:
            if verbose: print('\t. Pilatus data found, (column %d, alias %s)'%(column_z, stamps[column_z][1]))
        else:
            print(PN._RED,'\t. No pilatus data found', PN._RESET)
            nexus.close()
            sys.exit('Pilatus not found')

        # Check what is the x input (qxy or delta)
        if column_qxy is None:
            if column_delta is None:
                print(PN._RED,'\t. No usual actuator for GIXD found, stop here', PN._RESET)
                sys.exit('No sensor found')
            else:
                column_x=column_delta
                if verbose: print('\t. delta data found, (column %d, alias %s)'%(column_x, stamps[column_x][1]))
        else:
            column_x=column_qxy
            if verbose: print('\t. qxy data found, (column %d, alias %s)'%(column_x, stamps[column_x][1]))

        # Find start and stop without Nan
        istart=0
        istop=nbpts
        i=0
        while np.isnan(data[0][i]):
            istart=i
            i=i+1
        i=nbpts-1
        while np.isnan(data[0][i]):
            istop=i+1
            i=i-1
        istop=istop-1
        if verbose: print('\t. Valid data between points %d and %d'%(istart, istop))

        # Compute and print mean values
        if column_pi is not None:
            mean_pi=data[column_pi][istart:istop-1].mean()
            if verbose: print('\t. Surface pressure data found, mean value %3.4g ± %3.4g mN/m'
                              %(mean_pi,data[column_pi][istart:istop-1].std() ))
        else:
            print('\t. No surface pressure data found')
            mean_pi = None
            
        if column_area is not None:
            mean_area=data[column_area][istart:istop-1].mean()
            if verbose: print('\t. Area per molecule data found, mean value %3.4g ± %3.4g nm2 per molecule'
                              %(mean_area, data[column_area][istart:istop-1].std()))
        else:
            print('\t. No area per molecule data found')
            mean_area = None
            
        if column_gamma is not None:
            mean_gamma=data[column_gamma][istart:istop-1].mean()
            if verbose: print('\t. Gamma motor data found, mean value %3.4g deg'%(mean_gamma))
        else:
            print('\t. No gamma motor data found')
            mean_gamma = None
            if computeqz:
                print('\t. gamma is required to compute qz. Add gamma to the sensors.')
                nexus.close()
                sys.exit('gamma not found')

        # Load images
        stamps, images=nexus.extractData('2D')
        
        # Get positions of the dead pixels
        try:
            dead_pixels = np.genfromtxt('lib/extraction/common/dead_pixels.dat', dtype = 'uint16', delimiter = ', ')
        except:
            print('Careful: the file lib/extraction/common/dead_pixels.dat was not found. Taking no dead pixels.')
            dead_pixels = []
            
        
        daty=[]
        datyTop=[]
        datyBottom=[]
        datyFirstQuarter=[]
        mat=[]

        for i in range(istart, istop, 1):
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i+1, nbpts))
            sys.stdout.flush()
            
            image=images[0][i]
            
            for i,j in dead_pixels:
                image[i,j]=0.0
                
            # Keep only the ROI corresponding to the Soller
            ROI=[510, 350, 130, 692]
            image=image[ROI[1]:ROI[1]+ROI[3], ROI[0]:ROI[0]+ROI[2]]
            
            # Replace the intensity of the dead zones with a value of 0
            image=np.where(image<0., 0., image)
            
            # Rod (integration along the horizontal axis)
            rod=image.sum(axis=1)
            
            # Final image with rods stacked (dim: nb_angle x vertical_dim_ROI)
            mat.append(rod)
            
            # Integrated rods
            daty.append(rod.sum())
            
            # Integrate the rod on different parts of the detector only
            datyTop.append(rod[0:np.int(ROI[3]/2)].sum())
            datyBottom.append(rod[np.int(ROI[3]/2):ROI[3]].sum())
            datyFirstQuarter.append(rod[np.int(3*ROI[3]/4):ROI[3]].sum())
            
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()

        # Convert in numpy array
        daty = np.array(daty)[istart:istop]
        datyTop = np.array(datyTop)[istart:istop]
        datyBottom = np.array(datyBottom)[istart:istop]
        datyFirstQuarter = np.array(datyFirstQuarter)[istart:istop]
        mat = np.array(mat)        
        
        # Extract x values (qxy or delta)
        x = data[column_x]
        x = x[istart:istop]       
        
        # Bin the matrix, return binned vertical channels and binned matrix
        ch_binned, mat_binned = Groupe(mat, binsize=binsize)
            
        # Create Graph
        # The graph has to be splitted into two parts for good rendering in PDF
        fig=plt.figure(1, figsize=(12,5))
        fig.subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.16)
        ax1=fig.add_subplot(121)
        ax2=fig.add_subplot(122)

        # Plot the integrated spectra
        error=np.sqrt(daty)
        
        if (not logx) and (not logy):
            ax1.errorbar(x, daty, yerr=error, fmt='ro', label="full")
            ax1.plot(x, datyTop, 'k-', label='top')
            ax1.plot(x, datyBottom, 'b-', label='bottom')
            ax1.plot(x, datyFirstQuarter, 'r-', label='bottom quarter')
            ax1.legend()
        elif logx and (not logy):
            ax1.semilogx(x, daty, 'ro')
        elif (not logx) and (logy):
            ax1.semilogy(x, daty, 'ro')
        elif logx and (logy):
            ax1.loglog(x, daty, 'ro')
        ax1.set_xlabel(stamps0D[column_x][1]+' ('+stamps0D[column_x][2]+')', labelpad=13, fontsize='large')
        ax1.set_ylabel("Qz integrated intensity", labelpad=13, fontsize='large')

        # Plot the matrix
        if logz:
            ax2.contourf(x, ch_binned, np.log(mat_binned.transpose()))
        else:
            zmax=mat_binned.max()
            zmin=mat_binned.min()
            ax2.contourf(x, ch_binned, (mat_binned.transpose()),levels=np.linspace(zmin, zmax, nblevels), cmap=cmap)

        ax2.set_ylabel(r'$vertical\ channels$', fontsize='large')
        ax2.set_xlabel(stamps0D[column_x][1]+' ('+stamps0D[column_x][2]+')', labelpad=13, fontsize='large')
                
        if column_pi is not None:
            fig.text(.04, .05, r'$\pi = %3.4gmN.m^{-1}$'%(mean_pi), fontsize='large', color='red')
        if column_gamma is not None:
            fig.text(.96, .05, r'$\gamma = %3.4g deg$'%(mean_gamma), fontsize='large', color='red', horizontalalignment='right')

        rod=mat.sum(axis=0)
        fig.suptitle(nxs_filename[nxs_filename.rfind('/')+1:], fontsize='x-large')
        plt.show()
        
        fig=plt.figure(1, figsize=(12,8))
        ax3=fig.add_axes([0.125, 0.13, 0.775, 0.37])        
        ax3.plot(rod)
        
        # Extract the channel of the Vineyard peak
        channel0=rod.argmax()
        
        ax3.text(channel0*0.95, rod[channel0] , 'Channel of Vineyard Peak ($\mathregular{\\theta_c}$): %d'%(int(channel0)), 
                 fontsize='x-large', horizontalalignment='right', color='red')
        plt.plot((channel0, channel0), (rod[channel0]*1.1, rod[channel0]*1.3), 'r-', lw=2)
        ax3.set_xlabel('channels', fontsize='large')    
        ax3.set_ylabel('Q$\mathregular{_{xy}}$ - Integrated Intensity', fontsize='large')
        plt.show()
        
        print(PN._RED, 'Data not saved. To save data, run a GIXD on the scan.', PN._RESET)
        print(PN._RED, 'Channel0: %g'%channel0, PN._RESET)
       
        
    return channel0


    
def Groupe(mat, binsize=10):
    '''
    Bin a matrix along the vertical axis.

    Parameters
    ----------
    mat : array_like
        matrix to bin
    binsize : int, optional
        size in pixels of the vertical binning

    Returns
    -------
    tupple of arrays
        (ch_binned, mat_binned), array of channels after binning and matrix after binning
    '''
    mat_binned=[]
    
    for i in range(mat.shape[0]):
        ch_binned=[]
        z=[]
        j=0
        while j+binsize<mat.shape[1]:
            ch_binned.append(j+float(binsize)/2.0)
            z.append(mat[i, j:j+binsize].sum())
            j=j+binsize
        mat_binned.append(z)
        
    return (np.array(ch_binned), np.array(mat_binned))


def Calib_thetaz(calib_thetaz_data):
    '''
    Extract the factor for calibration of thetaz from a linear fit on gamma vs channel. Plot the calibration.

    Parameters
    ----------
    calib_thetaz_data : array_like
        numpy array containing the values of each gamma and corresponding channel.
        For example : np.array([
                      [0,      970],
                      [-1,     899],
                      [-2,     827]])


    Returns
    -------
    float
        thetazfactor, factor for conversion from channel to radian in the vertical direction (rad/channel)
    '''
    
    fig=plt.figure(1, figsize=(10,5))
    ax=fig.add_subplot(111)
    xfit=calib_thetaz_data[:,0]
    yfit=calib_thetaz_data[:,1]
    ax.plot(xfit,yfit, 'bo')
    A,B = Linear_fit(xfit,yfit)
    x=np.linspace(xfit.min(), xfit.max(), 100)
    ax.plot(x, Linear(x, A, B), 'r-', lw=2)
    
    ax.set_title('1D detector channel calibration')
    ax.set_xlabel('Gamma (deg)')
    ax.set_ylabel('Channels')

    fig.text(0.2, .8, "%3.5g channels per degree"%(B))
    fig.text(0.2, .75, "%3.5g degrees per channel"%(1.0/B))
    fig.text(0.2, .7, "%3.5g radians per channel"%((np.pi/B)/180.0))

    plt.show()
    
    thetazfactor = (np.pi/B)/180.0
    
    return thetazfactor


######################## FITTING FUNCTIONS #####################   

def Linear(x, A, B):
    '''Returns A+B*x'''
    return A+B*x


def Linear_fit(xfit, yfit, verbose=False):
    '''
    Linear fit with LMFIT : yfit = Coeff*xfit+Cste

    Parameters
    ----------
    xfit : array_like
        x array
    yfit : array_like
        y array
    verbose : bool, optional
        verbose mode

    Returns
    -------
    tupple of float
        (Cste, Coeff), the coeff and cste returned by the fit.
    '''

    def Residuals_Linear(params, x, y):
        return y-(params['Cste']+params['Coeff']*x)
    
    fitparams=lm.Parameters()
    nbpts=xfit.shape[0]
    B=(yfit[nbpts-1]-yfit[0])/(xfit[nbpts-1]-xfit[0])
    A=yfit[nbpts-1]-B*xfit[nbpts-1]
    fitparams.add_many(('Cste', A, True, -np.inf, yfit.max()*1.0, None),
                           ('Coeff', B, True, -10*B, 10*B, None),
                       )
    fitter = lm.Minimizer(Residuals_Linear, fitparams, fcn_args=(xfit, yfit))
    result=fitter.minimize()
        # Print result if asked via verbose
    if verbose:
        print(lm.fit_report(result))
        
    return (result.params['Cste'], result.params['Coeff'])

