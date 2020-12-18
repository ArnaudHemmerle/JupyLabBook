# custom libraries
from lib.extraction.common import PyNexus as PN

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import sys
import io
from contextlib import redirect_stdout
import lmfit as lm
from scipy.special import erf

def Treat(nxs_filename, recording_dir, xLabel, yLabel, working_dir='', show_data_stamps=False,
          plot=False, save=False, verbose=False):
    
    '''
    Call functions for extracting, plotting, and saving a 1D scan, i.e. extraction of two sensors from a scan.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    xLabel : str
        exact name of the x sensor, as it appears in the data stamps
    yLabel : str
        exact name of the y sensor, as it appears in the data stamps        
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
        xData, an array containing the list of x values
    array_like
        yData, an array containing the list of y values

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    SystemExit('Sensor not found')
        when sensor not found
    '''
        
    xData, yData = Extract(nxs_filename, recording_dir, xLabel, yLabel, show_data_stamps, verbose)
    
    if plot:
        Plot(xData, yData, xLabel, yLabel)
        
    if save:
        # For 1D plot we consider fast extraction only
        fast = True
        Save(nxs_filename, recording_dir, fast, working_dir, verbose)
        
    return  xData, yData

def Extract(nxs_filename, recording_dir, xLabel, yLabel, show_data_stamps, verbose):
    '''
    Extract the sensors.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    xLabel : str
        exact name of the x sensor, as it appears in the data stamps
    yLabel : str
        exact name of the y sensor, as it appears in the data stamps        
    show_data_stamps : bool, optional
        print the list of sensors from the nexus file
    verbose : bool, optional
        verbose mode

    Returns
    -------
    array_like
        xData, an array containing the list of x values
    array_like
        yData, an array containing the list of y values

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    SystemExit('Sensor not found')
        when sensor not found
    '''
    
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
        sys.exit('Nexus not found')
         
    else:
        if verbose: print(PN._BLUE+" - Open Nexus Data File :"+ PN._RESET)
        if verbose: print('\t'+nxs_path)
        try:
            nexus=PN.PyNexusFile(nxs_path, fast=True)
        except:
            print(PN._RED,'\t Nexus file seems not to exist or is not correct',PN._RESET)
            sys.exit('Nexus not found')
            
        nbpts=np.int(nexus.get_nbpts())
        if verbose: print("\t. Number of data points: ", nbpts)
            
        # Get stamps and Data
        stamps0D, data0D = nexus.extractData('0D')
        nexus.close()    
        
        # Explore stamps
        if show_data_stamps : print("\t. Available Counters:")
        for i in range(len(stamps0D)):
            if stamps0D[i][1] is not None:
                if show_data_stamps : print("\t\t", i, ' -------> ', stamps0D[i][1])
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps0D[i][0])
        
    
        sensor_list = [stamps0D[i][0] if stamps0D[i][1]== None else stamps0D[i][1] for i in range(len(stamps0D))]

        if xLabel in sensor_list:
            xArg = sensor_list.index(xLabel)
        else:
            print(PN._RED,'\t Sensor %s is not in the sensor list'%(xLabel),PN._RESET)
            sys.exit('Sensor not found')           
            
        if yLabel in sensor_list:
            yArg = sensor_list.index(yLabel)
        else:
            print(PN._RED,'\t Sensor %s is not in the sensor list'%(yLabel),PN._RESET)
            sys.exit('Sensor not found')  
          
        # Select the data to fit               
        istart=0
        i=0
        while np.isnan(data0D[xArg][i]):
            istart=i
            i=i+1
        while not np.isnan(data0D[yArg][i]) and i<data0D[yArg].shape[0]-1:
            istop=i+1
            i=i+1
            
        xData = data0D[xArg][istart:istop+1]
        yData= data0D[yArg][istart:istop+1]        
        
        return xData, yData
    
def Plot(xData, yData, xLabel, yLabel):
    '''
    1D Plot of the sensors.

    Parameters
    ----------
    xData : array_like
        list of x values
    yData : array_like
        list of y values        
    xLabel : str
        exact name of the x sensor, as it appears in the data stamps
    yLabel : str
        exact name of the y sensor, as it appears in the data stamps        

    '''    
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(xData, yData, 'o-')
    ax.set_xlabel(xLabel, fontsize=16)
    ax.set_ylabel(yLabel, fontsize=16)
    plt.show()

    
def Save(nxs_filename, recording_dir, fast, working_dir, verbose):
    '''
    Use the PyNexus library to convert the Nexus file into a .dat file.

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
    
    
######################## FITTING FUNCTIONS #####################    

def NormalFunction(x, mu, sigma, A, B, C):
    retour=(x-mu)/(sigma*np.sqrt(2.0))
    retour=(1/(sigma*np.sqrt(2.0*np.pi)))*np.exp(-retour*retour)
    return A+B*x+C*retour

def ResidualsNormFunction(params, x, y):
    # Unpack parameters
    A=params['Linear_Cste']
    B=params['Linear_Coeff']
    C=params['Amplitude']
    sigma=params['sigma']
    mu=params['mu']
    return y-NormalFunction(x, mu, sigma, A,B, C)

def GaussianFit(nxs_filename, recording_dir,
                xLabel='', yLabel='', verbose=False):
    
    '''
    Fit sensor_x vs sensor_y with a Gaussian.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    xLabel : str
        exact name of the x sensor, as it appears in the data stamps
    yLabel : str
        exact name of the y sensor, as it appears in the data stamps        
    verbose : bool, optional
        verbose mode

    Returns
    -------
    array_like
        xData, an array containing the list of x values
    array_like
        yData, an array containing the list of y values
    array_like
        yFit, an array containing the list of fitted y values

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    SystemExit('Sensor not found')
        when sensor not found
    '''
    
    xData, yData = Extract(nxs_filename, recording_dir, xLabel, yLabel, show_data_stamps=False, verbose=verbose)

    # Create the graph
    fig=plt.figure(1)
    ax=fig.add_subplot(111)

    # Fit the data using LMFIT
    # Define the parameters and first guess
    fitparams=lm.Parameters()
    nbpts=xData.shape[0]
    B=(yData[nbpts-1]-yData[0])/(xData[nbpts-1]-xData[0])
    A=yData[nbpts-1]-B*xData[nbpts-1]
    fitparams.add_many(('Linear_Cste', A, True, -np.inf, yData.max()*1.0, None),
                       ('Linear_Coeff', B, True, -10*B, 10*B, None),
                   ('Amplitude', yData.max(), True, 0.0, yData.max()*1.1, None),
                   ('sigma',np.abs(xData.max()-xData.min())/3., True, 0.0, xData.max()-xData.min(), None),
                   ('mu', (xData.min()+xData.max())/2., True, xData.min(), xData.max(), None),
                   )
    
    # Fit initialisation and fit
    fitter = lm.Minimizer(ResidualsNormFunction, fitparams, fcn_args=(xData, yData))
    result=fitter.minimize()   

    # Print result if asked via verbose
    if verbose:
        print(lm.fit_report(result))    

    yFit = NormalFunction(xData,result.params['mu'], result.params['sigma'],
                          result.params['Linear_Cste'], result.params['Linear_Coeff'], result.params['Amplitude'])
    
                
    # Plot first guess
    #ax.plot(xData,NormalFunction(xData,fitparams['mu'], fitparams['sigma'], fitparams['Linear_Cste'], fitparams['Linear_Coeff'], fitparams['Amplitude']), 'k--', lw=1)
        
    # plot the fitted data
    ax.plot(xData, yData, 'o-', label=nxs_filename[nxs_filename.rfind('_')+1:nxs_filename.rfind('.')])
    # plot the fit result
    ax.plot(xData, yFit, 'r-', lw=2)
    ax.legend(fontsize=16)
    ax.set_xlabel(xLabel, fontsize=16)
    ax.set_ylabel(yLabel, fontsize=16)
    ax.text(xData.min()*1.05, yData.max()*0.75,'Center %3.4g'%( result.params['mu']), fontsize=12)
    ax.text(xData.min()*1.05, yData.max()*0.65,'FWHM %3.4g'%(2.0*np.sqrt(2.0*np.log(2.))*result.params['sigma']), fontsize=12)
    ax.set_title('Gaussian Fit', fontsize=14)
    
    plt.show()

    return xData, yData, yFit
    

def NormRepFunction(x, mu, sigma, Amp, Cst, sens):
    """Normal repartition function"""
    retour=sens*(x-mu)/(sigma*np.sqrt(2.))
    retour=Cst+0.5*Amp*(1.+erf(retour))
    return retour

       
def ResidualsNormRepFunction(params, x, y, sens):
    # Unpack parameters
    C=params['Constant_Coeff']
    A=params['Amplitude']
    sigma=params['sigma']
    mu=params['mu']
    return y-NormRepFunction(x, mu, sigma, A, C, sens)
    
    
def GaussianRepartitionFit(nxs_filename, recording_dir,
                            xLabel='', yLabel='', verbose=False):
    '''
    Fit sensor_x vs sensor_y with an Erf function.

    Parameters
    ----------
    nxs_filename : str
        nexus filename
    recording_dir : str
        directory where the nexus file is stored
    xLabel : str
        exact name of the x sensor, as it appears in the data stamps
    yLabel : str
        exact name of the y sensor, as it appears in the data stamps        
    verbose : bool, optional
        verbose mode

    Returns
    -------
    array_like
        xData, an array containing the list of x values
    array_like
        yData, an array containing the list of y values
    array_like
        yFit, an array containing the list of fitted y values

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file not found
    SystemExit('Sensor not found')
        when sensor not found
    '''
    
    xData, yData = Extract(nxs_filename, recording_dir, xLabel, yLabel, show_data_stamps=False, verbose=verbose)
        
    # Find the decrease direction
    if xData[0]<xData[xData.shape[0]-1]:
        if yData[0]<yData[xData.shape[0]-1]:
            sens=1.
        else:
            sens=-1.
    else:
        if yData[0]>yData[xData.shape[0]-1]:
            sens=-1.
        else:
            sens=1.

    # Fit the data using LMFIT
    # Define the parameters and first guess
    fitparams=lm.Parameters()
    fitparams.add_many(('Constant_Coeff', yData.min(), True, yData.min()*0.9, yData.max()*1.1, None),
                   ('Amplitude', yData.max()-yData.min(), True, 0.0, yData.max()*1.1, None),
                   ('sigma',np.abs(xData.max()-xData.min())/3., True, 0.0, xData.max()-xData.min(), None),
                   ('mu', (xData.min()+xData.max())/2., True, xData.min(), xData.max(), None),
                   )
    
    # Fit initialisation and fit
    fitter = lm.Minimizer(ResidualsNormRepFunction, fitparams, fcn_args=(xData, yData, sens))
    result=fitter.minimize()
    # Print result if asked via verbose
    if verbose:
        print(lm.fit_report(result))    
        
    yFit = NormRepFunction(xData,result.params['mu'], result.params['sigma'],
                           result.params['Amplitude'], result.params['Constant_Coeff'], sens)
    

    # Create the graph
    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    
    # Plot first guess
    #ax.plot(xData,NormRepFunction(xData,fitparams['mu'], fitparams['sigma'], fitparams['Amplitude'], fitparams['Constant_Coeff'], sens), 'k--', lw=1)

    # Plot the fitted data
    ax.plot(xData, yData, 'o-', label=nxs_filename[nxs_filename.rfind('_')+1:nxs_filename.rfind('.')])
    # Plot the fit result
    ax.plot(xData, yFit, 'r-', lw=2)

    # Plot the associated gaussian function
    ax2=ax.twinx()
    ax2.plot(xData,NormalFunction(xData,result.params['mu'], result.params['sigma'], 0.0, 0.0, result.params['Amplitude']),
             'b-', lw=1)
    ax.legend(fontsize=16)
    ax.set_xlabel(xLabel, fontsize=16)
    ax.set_ylabel(yLabel, fontsize=16)
    if sens==1:
        fig.text(0.2, 0.65,'Center %3.4g'%( result.params['mu']), fontsize=12)
        fig.text(0.2, 0.55,'FWHM %3.4g'%(2.0*np.sqrt(2.0*np.log(2.))*result.params['sigma']), fontsize=12)
    else:
        fig.text(0.7, 0.65,'Center %3.4g'%( result.params['mu']), fontsize=12)
        fig.text(0.7, 0.55,'FWHM %3.4g'%(2.0*np.sqrt(2.0*np.log(2.))*result.params['sigma']), fontsize=12)
    ax.set_title('Normal Repartition Function Fit', fontsize=14)
    
    plt.show()
    
    return xData, yData, yFit
    