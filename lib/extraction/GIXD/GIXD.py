# custom libraries
from lib.extraction.common import PyNexus as PN


import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import os
import time
import subprocess
import sys
import lmfit as L
from scipy.special import erf
from PIL import Image
import io
from contextlib import redirect_stdout


__version__ = '1.0.3'

"""
Here are defined the custom functions used for analysis of data in the JupyLabBook.
Please observe a few rules for easier debugging:
- use the same keywords for modules as the one in the beginning of the file (np for numpy, plt for pyplot, ...)
- use the filename as the main file identifier (not the scan number, not the full path); i.e. your function should use
the arguments : nxs_filename='SIRIUS_test.nxs' and  recording_dir=''. 
To access the nexus file construct the path within the function (for ex. nxs_path = recording_dir+nxs_filename)
(The reason is that the scan_index may interfere with the day, the month or the year).
"""


##########################################################################################
###################################### GIXD ##############################################
##########################################################################################

def Extract_channel_Qc(nxs_filename='SIRIUS_test.nxs', working_dir='', recording_dir='',
                       logx=False, logy=False, logz=False):
    """
    Extract and return the channel corresponding to Qc from the position of the Vineyard's peak.
    """
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        column_z=None
        column_qxy=None
        column_delta=None
        column_pi=None
        column_area=None
        column_gamma=None
        print(PN._BLUE+" - Open Nexus Data File :"+ PN._RESET)
        print('\t'+nxs_path)
        try:
            nexus=PN.PyNexusFile(nxs_path)
        except:
            print(PN._RED,'\t Nexus file seems not to exist or is not correct',PN._RESET)
            return
        nbpts=np.int(nexus.get_nbpts())
        print("\t. Number of data points: ", nbpts)
        # Get stamps
        stamps=nexus.extractStamps()
        print("\t. Available Counters:")
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                print("\t\t", i, ' -------> ', stamps[i][1])
                if stamps[i][1].lower()=='pilatus':
                    column_z=i
            else:
                print("\t\t",i, ' -------> ', stamps[i][0])
        # Extract 0D data
        sys.stdout.write('Extracting 0D data\r')
        sys.stdout.flush()
        stamps0D, data=nexus.extractData('0D')
        sys.stdout.write('                  \r')
        sys.stdout.flush()
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
            print('\t. Pilatus data found, (column %d, alias %s)'%(column_z, stamps[column_z][1]))
        else:
            print(PN._RED,'\t. No pilatus data found', PN._RESET)
            nexus.close()
            return
        
        if column_qxy is None:
            if column_delta is None:
                print(PN._RED,'\t. No usual actuator for GIXD found, stop there', PN._RESET)
                nexus.close()
                return
            else:
                column_x=column_delta
                print('\t. delta data found, (column %d, alias %s)'%(column_x, stamps[column_x][1]))
        else:
            column_x=column_qxy
            print('\t. qxy data found, (column %d, alias %s)'%(column_x, stamps[column_x][1]))
        
        # find start and stop without Nan
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
            
        # Compute and print mean values
        if column_pi is not None:
            mean_pi=data[column_pi][istart:istop].mean()
            print('\t. Surface pressure data found, mean value %3.4g ± %3.4g mN/m'%(mean_pi,data[column_pi][istart:istop].std() ))
        else:
            print('\t. No surface pressure data found')
        if column_area is not None:
            mean_area=data[column_area][istart:istop].mean()
            print('\t. Area per molecule data found, mean value %3.4g ± %3.4g nm2 per molecule'%(mean_area, data[column_area][istart:istop].std()))
        else:
            print('\t. No area per molecule data found')
        if column_gamma is not None:
            mean_gamma=data[column_gamma][istart:istop].mean()
            print('\t. Gamma motor data found, mean value %3.4g deg'%(mean_gamma))
        else:
            print('\t. No gamma motor data found')
        
        # Load images
        stamps, images=nexus.extractData('2D')
        nexus.close()
        
        # Create the Qxy Qz map from images
        daty=[]
        datytop=[]
        datybottom=[]
        datyFirstQuarter=[]
        mat=[]
        
        for i in range(istart, istop, 1):
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i+1, nbpts))
            sys.stdout.flush()
            image=images[0][i]
            # dead pixels
            for i,j in dead_pixels:
                image[i,j]=0.0           
            # Keep only the ROI
            ROI=[510, 350, 130, 692]
            image=image[ROI[1]:ROI[1]+ROI[3], ROI[0]:ROI[0]+ROI[2]]
            image=np.where(image==-2, 0, image)
            mat.append(image.sum(axis=1))
            rod=image.sum(axis=1)
            daty.append(rod.sum())
            datytop.append(rod[0:np.int(ROI[3]/2)].sum())
            datybottom.append(rod[np.int(ROI[3]/2):ROI[3]].sum())
            datyFirstQuarter.append(rod[np.int(3*ROI[3]/4):ROI[3]].sum())
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()
            
        # Create Graph
        # The graph has to be splitted into two parts for good rendering in PDF
        fig=plt.figure(1, figsize=(12,5))
        fig.subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.16)
        ax1=fig.add_subplot(121)
        ax2=fig.add_subplot(122)
                        
        # convert new data as numpy array
        daty=np.array(daty)
        datytop=np.array(datytop)
        datybottom=np.array(datybottom)
        original_mat=np.array(mat)

        # extract x values (qxy or delta)
        x=data[column_x]

        # Plot the integrated spectra
        error=np.sqrt(daty)
        if (not logx) and (not logy):
            ax1.errorbar(x, daty, yerr=error, fmt='ro', label="full")
            ax1.plot(x, datytop, 'k-', label='top')
            ax1.plot(x, datybottom, 'b-', label='bottom')
            ax1.plot(x, datyFirstQuarter, 'r-', label='Bottom Quarter')
        elif logx and (not logy):
            ax1.semilogx(x, daty, fmt='ro')
        elif (not logx) and (logy):
            ax1.semilogy(x, daty, fmt='ro')
        elif logx and (logy):
            ax1.loglog(x, daty, fmt='ro')
        ax1.set_xlabel(stamps0D[column_x][1]+' ('+stamps0D[column_x][2]+')', labelpad=13, fontsize='large')
        ax1.set_ylabel("Qz integrated intensity", labelpad=13, fontsize='large')
                
        # Bin the matrix
        ch, mat= Groupe(original_mat, binsize=10)

        # Plot the matrix
        if logz:
            ax2.contourf(x[istart:istop], ch, np.log(mat.transpose()))
        else:
            ax2.contourf(x[istart:istop], ch, (mat.transpose()), cmap='jet')

        ax2.set_ylabel(r'$vertical\ channels$', fontsize='large')
        ax2.set_xlabel(stamps0D[column_x][1]+' ('+stamps0D[column_x][2]+')', labelpad=13, fontsize='large')
                
        if column_pi is not None:
            fig.text(.04, .05, r'$\pi = %3.4gmN.m^{-1}$'%(mean_pi), fontsize='large', color='red')
        if column_gamma is not None:
            fig.text(.96, .05, r'$\gamma = %3.4g deg$'%(mean_gamma), fontsize='large', color='red', horizontalalignment='right')

        original_mat=np.where(original_mat<=-2, 0, original_mat)
        rod=original_mat.sum(axis=0)
        fig.suptitle(nxs_filename[nxs_filename.rfind('/')+1:], fontsize='x-large')
        plt.show()
        
        fig=plt.figure(1, figsize=(12,8))
        ax3=fig.add_axes([0.125, 0.13, 0.775, 0.37])        
        ax3.plot(rod)
        i_max=rod.argmax()
        
        ax3.text(i_max*0.95, rod[i_max] , 'Channel of Vineyard Peak ($\mathregular{\\theta_c}$): %d'%(int(i_max)), 
                 fontsize='x-large', horizontalalignment='right', color='red')
        plt.plot((i_max, i_max), (rod[i_max]*1.1, rod[i_max]*1.3), 'r-', lw=2)
        ax3.set_xlabel('channels', fontsize='large')    
        ax3.set_ylabel('Q$\mathregular{_{xy}}$ - Integrated Intensity', fontsize='large')
        plt.show()
        
        print(PN._RED, 'Data not saved. To save data, run a GIXD on the scan.', PN._RESET)
        print(PN._RED, 'Channel0: %g'%i_max, PN._RESET)
       
        
    return i_max


def Groupe(mat, binsize=10):
    tmp=[]
    
    for i in range(mat.shape[0]):
        x=[]
        y=[]
        z=[]
        j=0
        while j+binsize<mat.shape[1]:
            y.append(j+float(binsize)/2.0)
            z.append(mat[i, j:j+binsize].sum())
            j=j+binsize
        tmp.append(z)
    return (np.array(y), np.array(tmp))


def Extract(nxs_filename, recording_dir,
            channel0, thetazfactor, wavelength, thetac,
            binsize, computeqz,
            show_data_stamps, verbose, absorbers):
    
    """
    Extract the GIXD scan. 
    """
    
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
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
            return
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
        
        # Get absorbers
        if absorbers != '':
            print("\t. Absorbers:", str(absorbers))
                    
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
            return

        if column_qxy is None:
            if column_delta is None:
                print(PN._RED,'\t. No usual actuator for GIXD found, stop there', PN._RESET)
                return
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
            if verbose: print('\t. Surface pressure data found, mean value %3.4g ± %3.4g mN/m'%(mean_pi,data[column_pi][istart:istop-1].std() ))
        else:
            print('\t. No surface pressure data found')
            mean_pi = None
        if column_area is not None:
            mean_area=data[column_area][istart:istop-1].mean()
            if verbose: print('\t. Area per molecule data found, mean value %3.4g ± %3.4g nm2 per molecule'%(mean_area, data[column_area][istart:istop-1].std()))
        else:
            print('\t. No area per molecule data found')
            mean_area = None
        if column_gamma is not None:
            mean_gamma=data[column_gamma][istart:istop-1].mean()
            if verbose: print('\t. Gamma motor data found, mean value %3.4g deg'%(mean_gamma))
        else:
            print('\t. No gamma motor data found')
            mean_gamma = None

        # Load images
        stamps, images=nexus.extractData('2D')
        
        # Get positions of the dead pixels
        dead_pixels = np.genfromtxt('lib/extraction/common/dead_pixels.dat', dtype = 'uint16', delimiter = ', ')

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
                
            # Keep only the ROI corresponding to the Sollers
            ROI=[510, 350, 130, 692]
            image=image[ROI[1]:ROI[1]+ROI[3], ROI[0]:ROI[0]+ROI[2]]
            
            # Replace the dead zones with a value of -2 with a value of 0
            image=np.where(image<0, 0, image)
           
            # rod (integration along the horizontal axis)
            rod=image.sum(axis=1)
            
            # final image with rods stacked (dim: nb_angle x vertical_dim_ROI)
            mat.append(rod)
            
            # integrated rod (for 2D plots)
            daty.append(rod.sum())
            
            # integrate the rod on different parts of the detector only
            datyTop.append(rod[0:np.int(ROI[3]/2)].sum())
            datyBottom.append(rod[np.int(ROI[3]/2):ROI[3]].sum())
            datyFirstQuarter.append(rod[np.int(3*ROI[3]/4):ROI[3]].sum())
            
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()

        # convert in numpy array
        daty = np.array(daty)[istart:istop]
        datyTop = np.array(datyTop)[istart:istop]
        datyBottom = np.array(datyBottom)[istart:istop]
        datyFirstQuarter = np.array(datyFirstQuarter)[istart:istop]
        mat = np.array(mat)            
        
        # extract x values (qxy or delta)
        x = data[column_x]
        x = x[istart:istop]

        # Bin the matrix, return vertical channels and matrix
        ch_binned, mat_binned = Groupe(mat, binsize=binsize)

        # extract y values (qz or vertical channels)
        if computeqz==True:
            # Compute and return qz
            thetaz=thetac+(mean_gamma*np.pi/180.0)+(channel0-ch_binned)*thetazfactor
            y=2.0*np.pi*np.sin(thetaz)/wavelength
        else:
            # Return the vertical channels 
            y=ch_binned
            
        # Nature and unit of the axis
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
         mat, mat_binned, ch_binned,
         mean_pi, mean_area, mean_gamma,
         nxs_filename, logx, logy, logz,
         nblevels, absorbers, cmap):

        
        # error bar along y
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

        
def Save(x, daty, datyTop, datyBottom, datyFirstQuarter, mat, moytocreate, mean_gamma, column_x,
         channel0, thetazfactor, wavelength, thetac,
         nxs_filename, recording_dir, working_dir,
         computeqz, verbose):   
    
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

        # extract y values (qz or vertical channels)
        if computeqz:
            # Compute and return qz
            thetaz=thetac+(mean_gamma*np.pi/180.0)+(channel0-ch_binned)*thetazfactor
            y=2.0*np.pi*np.sin(thetaz)/wavelength
        else:
            # Return the vertical channels 
            y=ch_binned
            
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

