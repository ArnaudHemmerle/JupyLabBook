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
    
    
def Treat(nxs_filename, recording_dir,
          xmin=0., xmax=980., ymin=0., ymax=1042.,
          absorbers='', logz=True, cmap='viridis',
          working_dir='', show_data_stamps=False, plot=False, save=False, verbose=False):
        
    '''
    Call functions for extracting, plotting, and saving a sum over the Pilatus images in a scan.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    xmin : float, optional
        min limit of the vertical profile plot (integrated over the horizontal axis)
    xmax : float, optional
        max limit of the vertical profile plot (integrated over the horizontal axis)
    ymin : float, optional
        min limit of the horizontal profile plot (integrated over the vertical axis)
    ymax : float, optional
        max limit of the horizontal profile plot (integrated over the vertical axis)
    absorbers : str, optional
        text to display indicating which absorber was used
    logz : bool, optional
        log on the image
    cmap : str, optional
        colormap of the image
    working_dir : str, optional
        directory where the treated files will be stored
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
        images_sum, the Pilatus image integrated over the scan (usually time).
    array_like
        integrated_x, an array containing the profile integrated along the horizontal axis  
    array_like
        integrated_y, an array containing the profile integrated along the vertical axis        


    Raises
    ------
    SystemExit('Pilatus not found')
        when Pilatus is not found
    '''

    images_sum, integrated_x, integrated_y = \
    Extract(nxs_filename, recording_dir, show_data_stamps, verbose)

    if plot:
        Plot(images_sum, integrated_x, integrated_y,
             nxs_filename, absorbers, logz, cmap,
             xmin, xmax, ymin, ymax)

    if save:
        Save(images_sum, integrated_x, integrated_y, nxs_filename,
             working_dir, verbose)
        
    return images_sum, integrated_x, integrated_y
         

def Extract(nxs_filename, recording_dir,
            show_data_stamps=False, verbose=False):
    '''
    Extract the Pilatus images from the scan and return the sum.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    show_data_stamps : bool
        print the list of sensors from the nexus file
    verbose : bool
        verbose mode


    Returns
    -------
    array_like
        images_sum, the Pilatus image integrated over the scan (usually time).
    array_like
        integrated_x, an array containing the profile integrated along the horizontal axis  
    array_like
        integrated_y, an array containing the profile integrated along the vertical axis        


    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    SystemExit('Pilatus not found')
        when Pilatus is not found
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

        nexus.close()  
        
        # Sum the images over time
        images_sum = images.sum(axis=0)
        
        # Replace dead zones with an intensity of -2.
        images_sum = np.where(images_sum<0., -2., images_sum)
        
        # Integrate over the horizontal axis
        # Put the dead zones to 0. for the integration
        integrated_x = np.where(images_sum>0, images_sum, 0.).sum(axis=1)
        
        # Integrate over the vertical axis
        # Put the dead zones to 0. for the integration
        integrated_y = np.where(images_sum>0, images_sum, 0.).sum(axis=0)
        
        return images_sum, integrated_x, integrated_y


def Plot(images_sum, integrated_x, integrated_y,
         nxs_filename, absorbers, logz, cmap,
         xmin, xmax, ymin, ymax):
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
   
    
    pixels_x = np.arange(0,np.shape(images_sum)[1],1)
    pixels_y = np.arange(0,np.shape(images_sum)[0],1)

    xx, yy = np.meshgrid(pixels_x, pixels_y)

    fig = plt.figure(figsize=(15,15))

    # Divide the grid in 2x2
    outer = gridspec.GridSpec(2, 2, wspace=0.2)

    # Divide the left row in 2x1
    inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                    subplot_spec=outer[0], hspace=0.5)

    ############################################
    # Plot a profile along y (integrated over x)
    ax1 = fig.add_subplot(inner[0])
    
    if logz: ax1.set_yscale('log')
    ax1.set(xlabel = 'vertical pixel (y)', ylabel = 'integration along x')

    try:
        # Compute best limits
        ax1.set_xlim(ymin,ymax)
        temp = integrated_x[int(ymin):int(ymax)]
        ax1.set_ylim(np.min(temp[temp>0])*0.8,np.max(integrated_x[int(ymin):int(ymax)])*1.2)
        ax1.plot(integrated_x)
    except:
        # Go back to automatic limits if bad limits given
        ax1.plot(integrated_x)

    ############################################
    # Plot a profile along x (integrated over y)
    ax2 = fig.add_subplot(inner[1])
    
    if logz: ax2.set_yscale('log')
    ax2.set(xlabel = 'horizontal pixel (x)', ylabel = 'integration along y')

    try:
        # Compute best limits
        ax2.set_xlim(xmin,xmax)
        temp = integrated_y[int(xmin):int(xmax)]
        ax2.set_ylim(np.min(temp[temp>0])*0.8,np.max(integrated_y[int(xmin):int(xmax)])*1.2)
        ax2.plot(integrated_y)
    except:
        # Go back to automatic limits if bad limits given
        ax2.plot(integrated_y)

    #Divide the right row in 1x1
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[1], wspace=0.1, hspace=0.1)

    ############################################
    # Show the full image integrated over the scan
    ax0 = fig.add_subplot(inner[0])
    if logz:
        ax0.pcolormesh(xx, yy, images_sum, cmap = cmap, shading = 'auto',
                       norm = colors.LogNorm(), rasterized=True)
    else:
        ax0.pcolormesh(xx, yy, images_sum, cmap = cmap, shading = 'auto',
                       rasterized=True)
    ax0.set(xlabel = 'horizontal pixel (x)', ylabel ='vertical pixel (y)')
    ax0.invert_yaxis()
    fig.subplots_adjust(top=0.95)
    fig.suptitle(nxs_filename.split('\\')[-1], fontsize=16)
    plt.show()


def Save(images_sum, integrated_x, integrated_y, nxs_filename, working_dir, verbose):
    '''
    Save.
    XXX_pilatus_sum.mat: the matrix corresponding to the image displayed, in ascii.
    XXX_pilatus_sum.tiff: the matrix corresponding to the image displayed, in tiff.

    XXX_integrated_x.dat: the profile integrated along the horizontal axis.
    XXX_integrated_y.dat: the profile integrated along the vertical axis.

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
    working_dir : str
        directory where the treated files will be stored
    verbose : bool
        verbose mode
    '''

    # Create Save Name
    savename=working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]

    # profile y
    np.savetxt(savename+'_integrated_y.dat', integrated_y,
               delimiter = '\t', comments='', header ='YIntegrated')
    # profile x
    np.savetxt(savename+'_integrated_x.dat', integrated_x,
               delimiter = '\t', comments='', header ='XIntegrated')

    # 2D matrix
    np.savetxt(savename+'_pilatus_sum.mat', images_sum)

    # tiff image
    im = Image.fromarray(images_sum)
    im.save(savename+'_pilatus_sum.tiff')

    if verbose: 
        print('\t. Original matrix saved in:')
        print("\t", savename+'.mat')
        print(" ")
        print('\t. Tiff saved in:')
        print("\t", savename+'.tiff')
        print(" ")
        print('\t. Profiles saved in:')
        print("\t", savename+'_integrated_x.dat')
        print("\t", savename+'_integrated_y.dat')
        print(" ")
    


