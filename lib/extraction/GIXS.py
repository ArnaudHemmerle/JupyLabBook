# custom libraries
from lib.extraction.common import PyNexus as PN

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from PIL import Image
import os
import sys

def Extract(nxs_filename, recording_dir,
            wavelength, thetai, distance,
            pixel_PONI_x, pixel_PONI_y, pixel_size,
            force_gamma_delta, fgamma, fdelta,
            show_data_stamps=False, verbose=False):

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
    show_data_stamps : bool, optional
        print the list of sensors from the nexus file
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
        
        i_pilatus=None
        
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
                    i_pilatus=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])    

        # Check that Pilatus data are present (images)
        if i_pilatus is not None:
            if verbose: print('\t. Pilatus data found, (column %d, alias %s)'%(i_pilatus, stamps[i_pilatus][1]))
        else:
            print(PN._RED,'\t. No pilatus data found', PN._RESET)
            nexus.close()
            sys.exit('Pilatus not found')
                    
                
        images = np.zeros([nbpts, 1043, 981])

        for i in range(nbpts):
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i+1, nbpts))
            sys.stdout.flush()

            #Extract the images from the nexus file
            stamp, image = nexus.extract_one_data_point(stamps[i_pilatus][0], i, verbose = False)

            
            # Get positions of the dead pixels
            try:
                dead_pixels = np.genfromtxt('lib/extraction/common/dead_pixels.dat', dtype = 'uint16', delimiter = ', ')
            except:
                print('Careful: the file lib/extraction/common/dead_pixels.dat was not found. Taking no dead pixels.')
                dead_pixels = []
            
            #Remove the dead pixels
            for ii,jj in dead_pixels:
                image[ii,jj]=0.0   

            images[i,:] = image    
            
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()

        #Extract the values of each elements of the nxs
        s, data = nexus.extractData('0D')
        nexus.close()

        i_gamma = None
        i_delta = None
        for i in range(len(stamps)):
            if (stamps[i][1] != None and stamps[i][1].lower()=='delta'):
                i_delta = i
            if (stamps[i][1] != None and stamps[i][1].lower()=='gamma'):
                i_gamma = i    

        if (i_gamma != None) and (force_gamma_delta==False):
            gamma = np.mean(data[i_gamma])
            print(PN._RED,'\t. Gamma found. gamma = %3.4g deg'%(gamma), PN._RESET)
        else:
            gamma = fgamma
            print("")
            print(PN._RED,'\t. No gamma found! gamma = %g'%gamma, PN._RESET)
        if (i_delta != None) and (force_gamma_delta==False):
            delta = np.mean(data[i_delta])
            print(PN._RED,'\t. Delta found. delta = %3.4g deg'%(delta), PN._RESET)
        else:
            delta = fdelta
            print(PN._RED,'\t. No delta found! delta = %g'%delta, PN._RESET)

        if verbose: 
            print('\t. For more details on the geometry, see:')
            print('\t \t -Fig.2 in doi:10.1107/S0909049512022017')
            print('\t \t -Slide 4 in http://gisaxs.com/files/Strzalka.pdf')            

        # Sum the images over time
        images_sum = images.sum(axis=0)
        
        # Replace dead zones with an intensity of -2.
        images_sum = np.where(images_sum<0., -2., images_sum)
        
        pixels_x = np.arange(0,np.shape(images_sum)[1],1)
        pixels_y = np.arange(0,np.shape(images_sum)[0],1)

        xx, yy = np.meshgrid(pixels_x, pixels_y)

        # alphai (incident angle)
        alphai = thetai    

        pixel_direct_x = pixel_PONI_x-distance/pixel_size*np.tan(delta*np.pi/180.)
        pixel_direct_y = pixel_PONI_y+distance/pixel_size*np.tan(gamma*np.pi/180.)

        # 2*theta in rad
        twotheta = np.arctan(pixel_size*(xx-pixel_direct_x)/distance)
        
        # alpha_f in rad
        deltay0 = distance*np.tan(alphai*np.pi/180.)
        alphaf = np.arctan( (pixel_size*(pixel_direct_y-yy)-deltay0)/distance)
        
        # q in nm^-1
        k0 = 2*np.pi/wavelength
        qz = k0*(np.sin(alphaf)+np.sin(alphai))
        qxy_approx = 2*k0*np.sin(twotheta/2.)

        # True values of q (not used here)
        #qx = k0*(np.cos(alphaf)*np.cos(twotheta)-np.cos(alphai))
        #qy = k0*np.cos(alphaf)*np.sin(twotheta)
        #qxy = np.sqrt(np.square(qx)+np.square(qy))
        #q = np.sqrt(np.square(qxy)+np.square(qz))
        
        # Profiles
        # Put all the negative pixels to zero before integration            
        integrated_qxy = np.where(images_sum>0, images_sum, 0.).sum(axis=1)
        integrated_qz = np.where(images_sum>0, images_sum, 0.).sum(axis=0)
        
        pixel_x_array = np.arange(0,len(integrated_qz))
        pixel_y_array = np.arange(0,len(integrated_qxy))
        
        alphaf_array = np.arctan( (pixel_size*(pixel_direct_y-pixel_y_array)-deltay0)/distance)
        qz_array = 2*np.pi/wavelength*(np.sin(alphaf_array)+np.sin(alphai))
        twotheta_array = np.arctan(pixel_size*(pixel_x_array-pixel_direct_x)/distance)
        # Careful, approx. qxy=4*pi/lambda*sin(2*theta/2)
        qxy_array = 4*np.pi/wavelength*np.sin(twotheta_array/2.)
        
        return images_sum, qxy_approx, qz,\
               integrated_qxy, integrated_qz, qxy_array, qz_array
    
    
def Plot(images_sum, qxy_approx, qz,
         integrated_qxy, integrated_qz, qxy_array, qz_array,
         nxs_filename, absorbers, logz, cmap,
         qxymin, qxymax, qzmin, qzmax):
    '''
    Plot the integrated image and the corresponding profiles.
    
    Parameters
    ----------
    images_sum : array_like
        the 2D image to plot
    integrated_x : array_like
        the profile integrated along the horizontal axis  
    integrated_y : array_like
        the profile integrated along the vertical axis  
    nxs_filename : str
        nexus filename
    absorbers : str
        text to display indicating which absorber was used
    logz : bool
        log on the image
    cmap : str
        colormap of the image
    xmin : float
        min limit of the vertical profile plot (integrated over the horizontal axis)
    xmax : float
        max limit of the vertical profile plot (integrated over the horizontal axis)
    ymin : float
        min limit of the horizontal profile plot (integrated over the vertical axis)
    ymax : float
        max limit of the horizontal profile plot (integrated over the vertical axis)
    '''    
    # Print absorbers
    if absorbers != '':
        print("\t. Absorbers:", str(absorbers))                 
   
    fig = plt.figure(figsize=(15,15))
    fig.subplots_adjust(top=0.95)
    fig.suptitle(nxs_filename.split('\\')[-1], fontsize='x-large')

    # Divide the grid in 2x2
    outer = gridspec.GridSpec(2, 2, wspace=0.2)

    # Divide the left row in 2x1
    inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                    subplot_spec=outer[0], hspace=0.5)


    #############################################
    # Plot a profile along y (integrated over x)
    ax1 = fig.add_subplot(inner[0])


    if logz: ax1.set_yscale('log')
    ax1.set(xlabel = 'qz (nm^-1)', ylabel = 'integration along qxy')

    try:
        # Compute best limits
        temp = integrated_qxy[(qz_array<qzmax) & (qzmin<qz_array)]
        qzmin_plot = np.min(temp[temp>0])   
        qzmax_plot = np.max(integrated_qxy[(qz_array<qzmax) & (qzmin<qz_array)])
        ax1.set_xlim(qzmin,qzmax)
        ax1.set_ylim(0.8*qzmin_plot,1.2*qzmax_plot) 
        ax1.plot(qz_array, integrated_qxy)
    except:
        # Go back to automatic limits if bad limits given
        ax1.plot(qz_array, integrated_qxy)

    ###################################################
    # Plot a profile along x (integrated over y)
    ax2 = fig.add_subplot(inner[1])

    if logz: ax2.set_yscale('log')
    ax2.set(xlabel = 'qxy (nm^-1)', ylabel = 'integration along qz')

    try:
        # Compute best limits
        temp = integrated_qz[(qxy_array<qxymax) & (qxymin<qxy_array)]
        qzmin_plot = np.min(temp[temp>0])   
        qzmax_plot = np.max(integrated_qz[(qxy_array<qxymax) & (qxymin<qxy_array)])  
        ax2.set_xlim(qxymin,qxymax) 
        ax2.set_ylim(0.8*qzmin_plot,1.2*qzmax_plot)    
        ax2.plot(qxy_array, integrated_qz)            
    except:
        # Go back to automatic limits if bad limits given
        ax2.plot(qxy_array, integrated_qz) 

    #Divide the right row in 1x1
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[1], wspace=0.1, hspace=0.1)

    #####################################################
    # Show the full image integrated over the scan
    ax2 = fig.add_subplot(inner[0])

    ax2.pcolormesh(qxy_approx, qz, images_sum, cmap = cmap, shading = 'auto',
                   norm = colors.LogNorm(), rasterized=True)
    ax2.set_xlabel('qxy (nm^-1)', fontsize='large')
    ax2.set_ylabel('qz (nm^-1)', fontsize='large')       

    plt.show()   