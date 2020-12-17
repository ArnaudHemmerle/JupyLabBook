# custom libraries
from lib.extraction.common import PyNexus as PN

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import os
import sys
import io
from contextlib import redirect_stdout


def Treat(nxs_filename, recording_dir, working_dir='', fast=True, show_data_stamps=False, plot=False, save=False, verbose=False):
    '''
    Call functions for extracting, plotting, and saving an Isotherm.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    working_dir : str, optional
        directory where the treated files will be stored
    fast : bool, optional
        triger fast extract of the nexus file
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
        area, an array containing the list of area values
    array_like
        pressure, an array containing the list of pressure values 
    array_like
        time, an array containing the list of time values      

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    '''
        
    area, pressure, time = Extract(nxs_filename, recording_dir, fast, show_data_stamps, verbose)
    
    if plot:
        Plot(area, pressure, time, nxs_filename)
        
    if save:
        Save(nxs_filename, recording_dir, fast, working_dir, verbose)
        
    return area, pressure, time
        
        
def Extract(nxs_filename, recording_dir, fast, show_data_stamps, verbose):
    '''
    Extract the isotherm.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    fast : bool
        triger fast extract of the nexus file
    show_data_stamps : bool
        print the list of sensors from the nexus file
    verbose : bool
        verbose mode

    Returns
    -------
    array_like
        area, an array containing the list of area values
    array_like
        pressure, an array containing the list of pressure values 
    array_like
        time, an array containing the list of time values      

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    '''

    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
        sys.exit('Nexus not found')
         
    else:
        column_pi=None
        column_area=None
        column_time=None
        if verbose: print(PN._BLUE+" - Open Nexus Data File :"+ PN._RESET)
        if verbose: print('\t'+nxs_path)
        try:
            nexus=PN.PyNexusFile(nxs_path, fast=fast)
        except:
            print(PN._RED,'\t Nexus file seems not to exist or is not correct',PN._RESET)
            sys.exit('Nexus not found')
            
        nbpts=np.int(nexus.get_nbpts())
        if verbose: print("\t. Number of data points: ", nbpts)
            
        # Get stamps and Data
        stamps, data = nexus.extractData('0D')
        nexus.close()
        
        # Explore stamps
        if show_data_stamps : print("\t. Available Counters:")
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                if show_data_stamps : print("\t\t", i, ' -------> ', stamps[i][1])
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])

        # Find columns
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                if stamps[i][1].lower()=='surfacepressure':
                    column_pi=i
                if stamps[i][1].lower()=='areapermolecule':
                    column_area=i
            if stamps[i][0].find('sensorsRelTimestamps')>-1:
                column_time=i

        if column_time is not None:
            if column_area is not None:
                if column_pi is not None:
                    cont=True
        if cont:
            if verbose:
                print('\t. Area per molecule found column %d'%(column_area))
                print('\t. Surface pressure per molecule found column %d'%(column_pi))
                print('\t. Time per molecule found column %d'%(column_time))
            
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
                istop=istop-1

                print('\t. Valid data between points %d and %d'%(istart, istop))
                
        area = data[column_area]
        pressure = data[column_pi]
        time = data[column_time]
        
        return area, pressure, time
        
 
def Plot(area, pressure, time, nxs_filename):
    '''
    Plot the isotherm.

    Parameters
    ----------
    area : array_like
        list of area values
    pressure : array_like
        list of pressure values
    time : array_like
        list of time values
    nxs_filename : str
        nexus filename
    '''    
    fig=plt.figure(1, figsize=(12,5))
    fig.subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.16)
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    ax3=ax2.twinx()

    ax1.plot(area, pressure, 'r-', lw=2)
    ax2.plot(time, area, 'k-', lw=2)
    ax3.plot(time, pressure, 'b-', lw=2)

    S=nxs_filename[nxs_filename.rfind('/')+1:nxs_filename.rfind('.')]

    fig.suptitle(S, fontsize=20) 

    ax1.set_xlabel('Area per Molecule (nm$\mathregular{^2}$)', fontsize=14)
    ax1.set_ylabel('Surface pressure (mN/m)', fontsize=14)

    ax2.set_xlabel('Time (sec)', fontsize=14)
    ax2.set_ylabel('Area per Molecule (nm$\mathregular{^2}$)', fontsize=14)
    ax3.set_ylabel('Surface pressure (mN/m)', fontsize=14, color='b')

    plt.show()
    
    
def Save(nxs_filename, recording_dir, fast, working_dir, verbose):
    '''
    Use the PyNexus library to convert the Nexus into a .dat file.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    fast : bool
        triger fast extract of the nexus file
    working_dir : str
        directory where the treated files will be stored
    verbose : bool
        verbose mode
    '''
   
    savename = working_dir + nxs_filename[:nxs_filename.rfind('.nxs')]

    # We assume extraction was already checked with Extract
    nxs_path = recording_dir+nxs_filename
    nexus=PN.PyNexusFile(nxs_path, fast=fast)
    
    # Get stamps and Data
    stamps, data = nexus.extractData('0D')
    
    f = io.StringIO()
    # Avoid printing sensors in the notebook
    with redirect_stdout(f):
        old_nexus_filename = nexus.filename
        # Save in working dir
        nexus.filename = working_dir+nxs_filename
        nexus.savePointExtractedData((stamps, data))
        nexus.filename = old_nexus_filename
    out = f.getvalue()        

    nexus.close()  
    
    if verbose:
        print('\t. Nexus converted to .dat file:')
        print("\t", savename+'.dat')
        print(" ")        
    
    