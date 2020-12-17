import PyNexus as PN
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

# Enter here the values of the dead pixels on the pilatus
dead_pixels = [
(877,528),
(1018,881),
(922,847),
(382,432),
(640,859),
(640,860) 
]

##########################################################################################
###################################### ALIGNMENT #########################################
##########################################################################################

def Gaussian_fit(nxs_filename='SIRIUS_test.nxs', recording_dir='',
                 xLabel='', yLabel='', xlog=False, ylog=False, verbose=False):
    
    """ 
    Fit the scan by a Gaussian function.
    
    Input parameters:
        nxs_filename : scan to fit, SIRIUS_year_month_day_index.nxs
        recording_dir : recording directory
        x : alias of the x-axis to use
        y : alias of the y-axis to use
        xlog : x-axis in log scale
        ylog : y-axis in log scale
        verbose : give details of the fit result
        
    Output :
        None
    """
    nxs_path = recording_dir+nxs_filename
    
    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        # Create the graph
        fig=plt.figure(1)
        ax=fig.add_subplot(111)
        
        # Open Nexus file and Read data
        try:
            nexus=PN.PyNexusFile(nxs_path)
            stamps, data= nexus.extractData('0D')
            nexus.close()
        except:
            print(PN._RED+'File %s seems not to be a regular nexus file'%(nxs_filename)+PN._RESET )
    
        # Find the sensors and actuators to plot
        ix=None
        iy=None
        for i in range(len(stamps)):
                if stamps[i][1] is not None:
                    if stamps[i][1].lower().find(xLabel.lower())>-1:
                        ix=i
                    if stamps[i][1].lower().find(yLabel.lower())>-1:
                        iy=i
        if ix is None:
            print(PN._RED+"%s is not in the list of registered sensors or actuators, stop there"%(xLabel)+PN._RESET)
            exit()
        if iy is None:
            print(PN._RED+"%s is not in the list of registered sensors or actuators, stop there"%(yLabel)+PN._RESET)
            exit()    
        # Select the data to fit               
        istart=0
        i=0
        while np.isnan(data[ix][i]):
            istart=i
            i=i+1
        while not np.isnan(data[iy][i]) and i<data[iy].shape[0]-1:
            istop=i+1
            i=i+1
        xfit=data[ix][istart:istop+1]
        yfit=data[iy][istart:istop+1]           
        
        # Fit the data using Larch module
        # Define the parameters and fist guess
        fitparams=L.Parameters()
        nbpts=xfit.shape[0]
        B=(yfit[nbpts-1]-yfit[0])/(xfit[nbpts-1]-xfit[0])
        A=yfit[nbpts-1]-B*xfit[nbpts-1]
        fitparams.add_many(('Linear_Cste', A, True, -np.inf, yfit.max()*1.0, None),
                           ('Linear_Coeff', B, True, -10*B, 10*B, None),
                       ('Amplitude', yfit.max(), True, 0.0, yfit.max()*1.1, None),
                       ('sigma',np.abs(xfit.max()-xfit.min())/3., True, 0.0, xfit.max()-xfit.min(), None),
                       ('mu', (xfit.min()+xfit.max())/2., True, xfit.min(), xfit.max(), None),
                       )
        # Plot first guess
        #ax.plot(xfit,NormalFunction(xfit,fitparams['mu'], fitparams['sigma'], fitparams['Linear_Cste'], fitparams['Linear_Coeff'], fitparams['Amplitude']), 'k--', lw=1)
        # Fit initialisation and fit
        fitter = L.Minimizer(ResidualsNormFunction, fitparams, fcn_args=(xfit, yfit))
        result=fitter.minimize()
        # Print result if asked via verbose
        if verbose:
            print(L.fit_report(result))
        # plot the fitted data
        ax.plot(data[ix], data[iy], 'o-', label=nxs_filename[nxs_filename.rfind('_')+1:nxs_filename.rfind('.')])
        # plot the fit result
        ax.plot(xfit,NormalFunction(xfit,result.params['mu'], result.params['sigma'], result.params['Linear_Cste'],
                                    result.params['Linear_Coeff'], result.params['Amplitude']), 'r-', lw=2)
        ax.legend(fontsize=16)
        ax.set_xlabel(stamps[ix][1], fontsize=16)
        ax.set_ylabel(stamps[iy][1], fontsize=16)
        ax.text(xfit.min()*1.05, yfit.max()*0.75,'Center %3.4g'%( result.params['mu']), fontsize=12)
        ax.text(xfit.min()*1.05, yfit.max()*0.65,'FWHM %3.4g'%(2.0*np.sqrt(2.0*np.log(2.))*result.params['sigma']), fontsize=12)
        ax.set_title('Gaussian Fit', fontsize=14)



def GaussianRepartition_fit(nxs_filename='SIRIUS_test.nxs', recording_dir='',
                            xLabel='', yLabel='', xlog=False, ylog=False, verbose=False):
    """ 
    Fit the scan by the repartition function of a Normal distribution (Gaussian Beam).
    
    Input parameters:
        nxs_filename : scan to fit, SIRIUS_year_month_day_index.nxs
        recording_dir : recording directory
        x : alias of the x-axis to use
        y : alias of the y-axis to use
        xlog : x-axis in log scale
        ylog : y-axis in log scale
        verbose : give details of the fit result
        
    Output :
        None
    """
    nxs_path = recording_dir+nxs_filename
    
    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        # Create the graph
        fig=plt.figure(1)
        ax=fig.add_subplot(111)
        
        # Open Nexus file and Read data
        try:
            nexus=PN.PyNexusFile(nxs_path)
            stamps, data= nexus.extractData('0D')
            nexus.close()
        except:
            print(PN._RED+'File %s seems not to be a regular nexus file'%(nxs_filename)+PN._RESET )
            
            
        # Find the sensors to extract and plot
        ix=None
        iy=None
        for i in range(len(stamps)):
                if stamps[i][1] is not None:
                    if stamps[i][1].lower().find(xLabel.lower())>-1:
                        ix=i
                    if stamps[i][1].lower().find(yLabel.lower())>-1:
                        iy=i
            
        if ix is None:
            print(PN._RED+"%s is not in the list of registered sensors or actuators, stop there"%(xLabel)+PN._RESET)
            exit()
        if iy is None:
            print(PN._RED+"%s is not in the list of registered sensors or actuators, stop there"%(yLabel)+PN._RESET)
            exit() 
    
        # Select the data to fit               
        istart=0
        i=0
        while np.isnan(data[ix][i]):
            istart=i
            i=i+1
        while not np.isnan(data[iy][i]) and i<data[iy].shape[0]-1:
            istop=i+1
            i=i+1
        xfit=data[ix][istart:istop+1]
        yfit=data[iy][istart:istop+1]           

        
        # Find the decrease direction
        if xfit[0]<xfit[xfit.shape[0]-1]:
            if yfit[0]<yfit[xfit.shape[0]-1]:
                sens=1.
            else:
                sens=-1.
        else:
            if yfit[0]>yfit[xfit.shape[0]-1]:
                sens=-1.
            else:
                sens=1.
        
        # Fit the data using LMFIT
        # Define the parameters and first guess
        fitparams=L.Parameters()
        fitparams.add_many(('Constant_Coeff', yfit.min(), True, yfit.min()*0.9, yfit.max()*1.1, None),
                       ('Amplitude', yfit.max()-yfit.min(), True, 0.0, yfit.max()*1.1, None),
                       ('sigma',np.abs(xfit.max()-xfit.min())/3., True, 0.0, xfit.max()-xfit.min(), None),
                       ('mu', (xfit.min()+xfit.max())/2., True, xfit.min(), xfit.max(), None),
                       )
        # Plot first guess
        #ax.plot(xfit,NormRepFunction(xfit,fitparams['mu'], fitparams['sigma'], fitparams['Amplitude'], fitparams['Constant_Coeff'], sens), 'k--', lw=1)
        # Fit initialisation and fit
        fitter = L.Minimizer(ResidualsNormRepFunction, fitparams, fcn_args=(xfit, yfit, sens))
        result=fitter.minimize()
        # Print result if asked via verbose
        if verbose:
            print(L.fit_report(result))
        # Plot the fitted data
        ax.plot(data[ix], data[iy], 'o-', label=nxs_filename[nxs_filename.rfind('_')+1:nxs_filename.rfind('.')])
        # Plot the fit result
        ax.plot(xfit,NormRepFunction(xfit,result.params['mu'], result.params['sigma'], result.params['Amplitude'],
                                     result.params['Constant_Coeff'], sens), 'r-', lw=2)

        # Plot the associated gaussian function
        ax2=ax.twinx()
        ax2.plot(xfit,NormalFunction(xfit,result.params['mu'], result.params['sigma'], 0.0, 0.0, result.params['Amplitude']),
                 'b-', lw=1)
        ax.legend(fontsize=16)
        ax.set_xlabel(stamps[ix][1], fontsize=16)
        ax.set_ylabel(stamps[iy][1], fontsize=16)
        if sens==1:
            fig.text(0.2, 0.65,'Center %3.4g'%( result.params['mu']), fontsize=12)
            fig.text(0.2, 0.55,'FWHM %3.4g'%(2.0*np.sqrt(2.0*np.log(2.))*result.params['sigma']), fontsize=12)
        else:
            fig.text(0.7, 0.65,'Center %3.4g'%( result.params['mu']), fontsize=12)
            fig.text(0.7, 0.55,'FWHM %3.4g'%(2.0*np.sqrt(2.0*np.log(2.))*result.params['sigma']), fontsize=12)
        ax.set_title('Normal Repartition Function Fit', fontsize=14)

        
def ResidualsNormRepFunction(params, x, y, sens):
    # Unpack parameters
    C=params['Constant_Coeff']
    A=params['Amplitude']
    sigma=params['sigma']
    mu=params['mu']
    return y-NormRepFunction(x, mu, sigma, A, C, sens)


def ResidualsNormFunction(params, x, y):
    # Unpack parameters
    A=params['Linear_Cste']
    B=params['Linear_Coeff']
    C=params['Amplitude']
    sigma=params['sigma']
    mu=params['mu']
    return y-NormalFunction(x, mu, sigma, A,B, C)

def NormRepFunction(x, mu, sigma, Amp, Cst, sens):
    """
    Normal repartition function
    """
    retour=sens*(x-mu)/(sigma*np.sqrt(2.))
    retour=Cst+0.5*Amp*(1.+erf(retour))
    return retour

def NormalFunction(x, mu, sigma, A, B, C):
    retour=(x-mu)/(sigma*np.sqrt(2.0))
    retour=(1/(sigma*np.sqrt(2.0*np.pi)))*np.exp(-retour*retour)
    return A+B*x+C*retour


def Linear(x, A, B):
    return A+B*x

def Residuals_Linear(params, x, y):
    return y-(params['Cste']+params['Coeff']*x)

def Linear_fit(xfit,yfit, verbose=False):
    fitparams=L.Parameters()
    nbpts=xfit.shape[0]
    B=(yfit[nbpts-1]-yfit[0])/(xfit[nbpts-1]-xfit[0])
    A=yfit[nbpts-1]-B*xfit[nbpts-1]
    fitparams.add_many(('Cste', A, True, -np.inf, yfit.max()*1.0, None),
                           ('Coeff', B, True, -10*B, 10*B, None),
                       )
    fitter = L.Minimizer(Residuals_Linear, fitparams, fcn_args=(xfit, yfit))
    result=fitter.minimize()
        # Print result if asked via verbose
    if verbose:
        print(L.fit_report(result))
    return (result.params['Cste'], result.params['Coeff'])

  
##########################################################################################
###################################### XRF ###############################################
##########################################################################################

def Extract_XRF(nxs_filename='SIRIUS_test.nxs', working_dir='', recording_dir='',
                logz=False, list_elems=[0,1,2,3], first_channel=0, last_channel=2048,
                use_eV=False, gain=10., eV0=0., arr_peaks=None,
                show_data_stamps=False, verbose=False, absorbers='', fast=True,
                plot_spectrogram=False, plot_first_last=False, plot_sum=False):
    """
    Extract, correct with ICR/OCR, and plot the fluo spectrum. 
    """
        
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        if verbose: print(PN._BLUE+" - Open Nexus Data File :"+ PN._RESET)
        if verbose: print('\t'+nxs_path)
        try:            
            nexus=PN.PyNexusFile(nxs_path, fast=fast)
        except:
            print(PN._RED,'\t Nexus file seems not to exist or is not correct',PN._RESET)
            return
        nbpts=np.int(nexus.get_nbpts())
        if verbose: print("\t. Number of data points: ", nbpts)
            
        # Get stamps
        stamps, data= nexus.extractData()

        # Save data
        f = io.StringIO()
        # Avoid printing sensors in the notebook
        with redirect_stdout(f):
            old_nexus_filename = nexus.filename
            # Save in working dir
            nexus.filename = working_dir+nxs_filename
            nexus.savePointExtractedData((stamps, data))
            nexus.saveOneDExtractedData((stamps, data))
            nexus.filename = old_nexus_filename
        out = f.getvalue()
        
        if verbose: 
            savename = working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]
            print('\t. 0D data saved in:')
            print("\t", savename+'.dat')
            print('\t. Spectrum(s) saved in:')
            print("\t", savename+'_fluospectrum*.mat')
                

        nexus.close()
        if show_data_stamps : print("\t. Available Counters:")
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                if show_data_stamps : print("\t\t", i, ' -------> ', stamps[i][1])
                if stamps[i][1].lower()=='pilatus':
                    columnz=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])
         
        # Get absorbers
        if absorbers != '':
            print("\t. Absorbers:", str(absorbers))        
        
    def extract_and_correct(ind_spectrum):
        """Extract the requested fluospectrum from the nexus file and correct it with ICR/OCR"""

        is_icr_found = False
        is_ocr_found = False
        for i in range(len(stamps)):
            if (stamps[i][1] != None and stamps[i][1].lower() == "fluoicr0"+ind_spectrum):
                fluoicr = data[i]
                is_icr_found = True
            if (stamps[i][1] != None and stamps[i][1].lower() == "fluoocr0"+ind_spectrum):
                fluoocr = data[i]
                is_ocr_found = True
            if (stamps[i][1] != None and stamps[i][1].lower() == "fluospectrum0"+ind_spectrum):
                fluospectrum = data[i]
            if (stamps[i][1] == None and stamps[i][0].lower() == "integration_time"):
                integration_time = data[i]
                
        if is_icr_found:
            ICR = fluoicr
            if is_ocr_found:
                OCR = fluoocr
            else:
                print(PN._RED+"OCR not found in data. Taking OCR = spectrum_intensity/counting_time."+PN._RESET)
                OCR = np.array([np.sum(fluospectrum[n])/integration_time[n] for n in range(len(fluospectrum))])
                
            ratio = np.array([ICR[n]/OCR[n] if (~np.isclose(OCR[n],0.) & ~np.isnan(OCR[n]) & ~np.isnan(ICR[n]))
                              else 0. for n in range(len(ICR))])
            spectrums_corr = np.array([fluospectrum[n]*ratio[n] for n in range(len(ratio))])
            is_extract_ok = True
            return is_extract_ok, spectrums_corr
                
        else:
            print(PN._RED+"ICR not found in data. Check if the box \'Elements\' is right."+PN._RESET)
            print(PN._RED+"Try to put 4 in the box \'Elements\' for the single-element detector."+PN._RESET)
            print(PN._RED+"Try to put 0, 1, 2, 3 in the box \'Elements\' for the four-elements detector."+PN._RESET)

            is_extract_ok = False
            return is_extract_ok, None
           
    # Correct each chosen element with ICR/OCR and sum them
    allspectrums_corr = np.zeros((nbpts, 2048))

    for i in list_elems:
        is_extract_ok, allspectrums_corr_i = extract_and_correct(str(i))
        if is_extract_ok:
            allspectrums_corr  += allspectrums_corr_i
        else:
            break
    
    if is_extract_ok:
        ind_non_zero_spectrums = np.where(np.sum(allspectrums_corr, axis = 1)>10.)[0]
        list_ranges = np.split(ind_non_zero_spectrums, np.where(np.diff(ind_non_zero_spectrums) != 1)[0]+1)
        first_non_zero_spectrum = ind_non_zero_spectrums[0]
        last_non_zero_spectrum = ind_non_zero_spectrums[-1]

        channels = np.arange(int(first_channel), int(last_channel+1))
        eVs = channels*gain+eV0
        spectrums = allspectrums_corr[0:last_non_zero_spectrum+1,
                                      int(first_channel):int(last_channel+1)]

        if plot_spectrogram:

            fig = plt.figure(figsize=(12,4.6))
            ax1 = fig.add_subplot(111)
            ax1.set_title(nxs_filename.split('\\')[-1], fontsize='x-large')
            ax1.set_xlabel('spectrum index', fontsize='large')
            ax1.set_xlim(left = 0, right = last_non_zero_spectrum)

            if use_eV:
                xx, yy = np.meshgrid(np.arange(0,last_non_zero_spectrum+1), eVs)
                ax1.set_ylabel('eV', fontsize='large')
            else:
                xx, yy = np.meshgrid(np.arange(0,last_non_zero_spectrum+1), channels)
                ax1.set_ylabel('channel', fontsize='large')          

            if logz:
                ax1.pcolormesh(xx, yy, spectrums.transpose(), cmap='viridis',  shading = 'auto',
                               norm = colors.LogNorm(), rasterized=True)
            else:
                ax1.pcolormesh(xx, yy, spectrums.transpose(), cmap='viridis',  shading = 'auto',
                               rasterized=True)

            plt.show()

        if plot_sum:
            fig = plt.figure(figsize=(12,4.5))
            ax1 = fig.add_subplot(111)
            ax1.set_ylabel('counts', fontsize='large')
            if logz: ax1.set_yscale('log')
            if use_eV:
                ax1.set_xlabel('eV', fontsize='large')
                line1, = ax1.plot(eVs, np.sum(spectrums, axis = 0), 'b.-', label='Sum of spectrums')   
            else:
                ax1.set_xlabel('channel', fontsize='large')
                line1, = ax1.plot(channels, np.sum(spectrums, axis = 0), 'b.-', label='Sum of spectrums') 

            if arr_peaks[0][0]!= None :  
        
                # Plot the peak positions
               
                # Prepare a list of colors and linestyles
                colors_axv = iter(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
                                   '#7f7f7f', '#bcbd22', '#17becf']*20)   
                linestyles_axv = iter(['--', '-.', '-', ':']*40)
             
                # Rearrange the peaks to plot them by increasing energy
                arr_peaks = np.array(arr_peaks)
                arg_position_peaks = np.argsort([float(elem[1]) for elem in arr_peaks])
                val_position_peaks = arr_peaks[arg_position_peaks][:,1]
                labels_peaks = arr_peaks[arg_position_peaks][:,0]
                
                axvlines = []
                for i in range(len(arr_peaks)):
                    axvlines.append(ax1.axvline(float(val_position_peaks[i]), label = str(labels_peaks[i]),
                                                color = next(colors_axv), linestyle = next(linestyles_axv)))
                
                axvlegends = ax1.legend(handles=axvlines, fontsize=10, ncol = len(arr_peaks)//16+1,
                                        bbox_to_anchor=(1.01, 1.), loc='upper left',  borderaxespad=0.) 
                plt.gca().add_artist(axvlegends)  

            ax1.legend(handles=[line1], fontsize='large', loc='upper left')     
            plt.show()

        if plot_first_last:    
            #Plot the selected channel range
            fig = plt.figure(figsize=(12,4.5))
            ax1 = fig.add_subplot(111)
            ax1.set_ylabel('counts', fontsize='large')
            if logz: ax1.set_yscale('log')
            if use_eV:
                ax1.set_xlabel('eV', fontsize='large')        
                line1, = ax1.plot(eVs, spectrums[first_non_zero_spectrum], 'b.-', label='First spectrum')
                line2, = ax1.plot(eVs, spectrums[-1], 'r.-', label='Last spectrum')            
            else:
                ax1.set_xlabel('channel', fontsize='large')        
                line1, = ax1.plot(channels, spectrums[first_non_zero_spectrum], 'b.-', label='First spectrum')
                line2, = ax1.plot(channels, spectrums[-1], 'r.-', label='Last spectrum')
            
            if arr_peaks[0][0]!= None :  
                
                # Rearrange the peaks to plot them by increasing energy
                arr_peaks = np.array(arr_peaks)
                arg_position_peaks = np.argsort([float(elem[1]) for elem in arr_peaks])
                val_position_peaks = arr_peaks[arg_position_peaks][:,1]
                labels_peaks = arr_peaks[arg_position_peaks][:,0]
                
                axvlines = []
                for i in range(len(arr_peaks)):
                    axvlines.append(ax1.axvline(float(val_position_peaks[i]), label = str(labels_peaks[i]),
                                                color = next(colors_axv), linestyle = next(linestyles_axv)))
                
                axvlegends = ax1.legend(handles=axvlines, fontsize=10, ncol = len(arr_peaks)//16+1,
                                        bbox_to_anchor=(1.01, 1.), loc='upper left',  borderaxespad=0.)
                plt.gca().add_artist(axvlegends)  
            
            ax1.legend(handles=[line1, line2], fontsize='large', loc='upper left')    
            plt.show()


##########################################################################################
###################################### GENERAL ###########################################
##########################################################################################

def Plot_1D(nxs_filename='SIRIUS_test.nxs', working_dir='', recording_dir='',
            xLabel='zs', yLabel='alphax'):
    """
    Simple 1D plot. When called by FrontendFunctions by clicking on Treat Scan, 
    it will plot the current selection in the interactive 1D plot.
    It is necessary to call this function to have the plot in the pdf.
    """
    nexus = PN.PyNexusFile(recording_dir+nxs_filename, fast=True)
    stamps0D, data0D = nexus.extractData('0D')
    
    # Save data
    f = io.StringIO()
    # Avoid printing sensors in the notebook
    with redirect_stdout(f):
        old_nexus_filename = nexus.filename
        # Save in working dir
        nexus.filename = working_dir+nxs_filename
        nexus.savePointExtractedData((stamps0D, data0D))
        nexus.filename = old_nexus_filename
    out = f.getvalue()
   
    nexus.close()
    sensor_list = [stamps0D[i][0] if stamps0D[i][1]== None else stamps0D[i][1] for i in range(len(stamps0D))]

    xArg = sensor_list.index(xLabel)
    yArg = sensor_list.index(yLabel)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(data0D[xArg], data0D[yArg], 'o-')
    ax.set_xlabel(xLabel, fontsize=16)
    ax.set_ylabel(yLabel, fontsize=16)
    plt.show()

    
