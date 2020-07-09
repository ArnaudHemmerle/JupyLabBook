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


__version__ = '0.17'

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

def Calib_thetaz(data):
    """
    Function to do the channel-thetaz calibration.
    Return the calibration factor channel to thetaz in rad/chan.
    """
    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    xfit=data[:,0]
    yfit=data[:,1]
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

    return (np.pi/B)/180.0
    
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
        columnz=None
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
                    columnz=i
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
        if columnz is not None:
            print('\t. Pilatus data found, (column %d, alias %s)'%(columnz, stamps[columnz][1]))
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
                columnx=column_delta
                print('\t. delta data found, (column %d, alias %s)'%(columnx, stamps[columnx][1]))
        else:
            columnx=column_qxy
            print('\t. qxy data found, (column %d, alias %s)'%(columnx, stamps[columnx][1]))
        
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
        datx=data[columnx]

        # Plot the integrated spectra
        error=np.sqrt(daty)
        if (not logx) and (not logy):
            ax1.errorbar(datx, daty, yerr=error, fmt='ro', label="full")
            ax1.plot(datx, datytop, 'k-', label='top')
            ax1.plot(datx, datybottom, 'b-', label='bottom')
            ax1.plot(datx, datyFirstQuarter, 'r-', label='Bottom Quarter')
        elif logx and (not logy):
            ax1.semilogx(datx, daty, fmt='ro')
        elif (not logx) and (logy):
            ax1.semilogy(datx, daty, fmt='ro')
        elif logx and (logy):
            ax1.loglog(datx, daty, fmt='ro')
        ax1.set_xlabel(stamps0D[columnx][1]+' ('+stamps0D[columnx][2]+')', labelpad=13, fontsize='large')
        ax1.set_ylabel("Qz integrated intensity", labelpad=13, fontsize='large')
                
        # Bin the matrix
        ch, mat= Groupe(original_mat, binsize=10)

        # Plot the matrix
        if logz:
            ax2.contourf(datx[istart:istop], ch, np.log(mat.transpose()))
        else:
            ax2.contourf(datx[istart:istop], ch, (mat.transpose()), cmap='jet')

        ax2.set_ylabel(r'$vertical\ channels$', fontsize='large')
        ax2.set_xlabel(stamps0D[columnx][1]+' ('+stamps0D[columnx][2]+')', labelpad=13, fontsize='large')
                
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


def Extract_GIXD(nxs_filename='SIRIUS_test.nxs', working_dir='', recording_dir='',
                 logx=False, logy=False, logz=False,
                 channel0=600, thetazfactor=0.01, wavelength=0.155, thetac=0.0028, thetai=0.002,
                 binsize=10, computeqz=True, nblevels=50, moytocreate=(10, 20, 40),
                 show_data_stamps=False, verbose=False, absorbers='', cmap='jet', plot_true_GIXD=False):
    
    """
    Extract, plot, and save the GIXD scan. 
    """
    
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        columnz=None
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
                    columnz=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])
        
        # Get absorbers
        if absorbers != '':
            print("\t. Absorbers:", str(absorbers))
                    
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
        if columnz is not None:
            if verbose: print('\t. Pilatus data found, (column %d, alias %s)'%(columnz, stamps[columnz][1]))
        else:
            print(PN._RED,'\t. No pilatus data found', PN._RESET)
            nexus.close()
            return

        if column_qxy is None:
            if column_delta is None:
                print(PN._RED,'\t. No usual actuator for GIXD found, stop there', PN._RESET)
                return
            else:
                columnx=column_delta
                if verbose: print('\t. delta data found, (column %d, alias %s)'%(columnx, stamps[columnx][1]))
        else:
            columnx=column_qxy
            if verbose: print('\t. qxy data found, (column %d, alias %s)'%(columnx, stamps[columnx][1]))

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
        if column_area is not None:
            mean_area=data[column_area][istart:istop-1].mean()
            if verbose: print('\t. Area per molecule data found, mean value %3.4g ± %3.4g nm2 per molecule'%(mean_area, data[column_area][istart:istop-1].std()))
        else:
            print('\t. No area per molecule data found')
        if column_gamma is not None:
            mean_gamma=data[column_gamma][istart:istop-1].mean()
            if verbose: print('\t. Gamma motor data found, mean value %3.4g deg'%(mean_gamma))
        else:
            print('\t. No gamma motor data found')

        # Load images
        stamps, images=nexus.extractData('2D')

        # Create the Qxy Qz map from images
        daty=[]
        datyTop=[]
        datyBottom=[]
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
            image=np.where(image<0, 0, image)
            mat.append(image.sum(axis=1))

            #Create rods
            rod=image.sum(axis=1)
            daty.append(rod.sum())
            datyTop.append(rod[0:np.int(ROI[3]/2)].sum())
            datyBottom.append(rod[np.int(ROI[3]/2):ROI[3]].sum())
            datyFirstQuarter.append(rod[np.int(3*ROI[3]/4):ROI[3]].sum())
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()

        # Create Graph
        fig=plt.figure(nxs_filename, figsize=(12,5))
        fig.subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.16)
        ax1=fig.add_subplot(121)
        ax2=fig.add_subplot(122)

        # convert new data as numpy array
        daty=np.array(daty)
        datyTop=np.array(datyTop)
        datyBottom=np.array(datyBottom)
        datyFirstQuarter=np.array(datyFirstQuarter)
        original_mat=np.array(mat)

        # extract x values (qxy or delta)
        datx=data[columnx]

        # Plot the integrated spectra
        error=np.sqrt(daty)
        if (not logx) and (not logy):
            ax1.errorbar(datx[istart:istop], daty[istart:istop], yerr=error, fmt='ro', label="full")
            ax1.plot(datx[istart:istop], datyTop[istart:istop], 'k-', label='top')
            ax1.plot(datx[istart:istop], datyBottom[istart:istop], 'b-', label='bottom')
            ax1.plot(datx[istart:istop], datyFirstQuarter[istart:istop], 'r-', label='bottom quarter')
            ax1.legend()
        elif logx and (not logy):
            ax1.semilogx(datx[istart:istop], daty[istart:istop], 'ro')
        elif (not logx) and (logy):
            ax1.semilogy(datx[istart:istop], daty[istart:istop], 'ro')
        elif logx and (logy):
            ax1.loglog(datx[istart:istop], daty[istart:istop], 'ro')
        ax1.set_xlabel(stamps0D[columnx][1]+' ('+stamps0D[columnx][2]+')', labelpad=13, fontsize='large')
        ax1.set_ylabel("Qz integrated intensity", labelpad=13, fontsize='large')

        # Bin the matrix
        ch, mat= Groupe(original_mat, binsize=binsize)

        # Compute qz
        if computeqz==True:
            thetaz=thetac+(mean_gamma*np.pi/180.0)+(channel0-ch)*thetazfactor
            qz=2.0*np.pi*np.sin(thetaz)/wavelength
        else:
            qz=ch

        # Plot the matrix
        if logz:
            ax2.contourf(datx[istart:istop], qz, np.log(mat.transpose()))
        else:
            zmax=mat.max()
            zmin=mat.min()
            ax2.contourf(datx[istart:istop], qz, (mat.transpose()),levels=np.linspace(zmin, zmax, nblevels), cmap=cmap)

        if computeqz:
            ax2.set_ylabel('Qz (nm-1)', fontsize='large')
        else:
            ax2.set_ylabel(r'$vertical\ channels$', fontsize='large')
        ax2.set_xlabel(stamps0D[columnx][1]+' ('+stamps0D[columnx][2]+')', labelpad=13, fontsize='large')
        
        if column_pi is not None:
            fig.text(.04, .05, r'$\pi = %3.4gmN.m^{-1}$'%(mean_pi),
                     fontsize='large', color='red')
        if column_gamma is not None:
            fig.text(.96, .05, r'$\gamma = %3.4g deg$'%(mean_gamma),
                     fontsize='large', color='red', horizontalalignment='right')
        
        original_mat=np.where(original_mat<=-2, 0, original_mat)

        S_suptitle=nxs_filename[nxs_filename.rfind('/')+1:]
        fig.suptitle(S_suptitle, fontsize='x-large')


        if plot_true_GIXD:
            """
            Plot the true Qx, Qy, and Qz without the approx. of alphai~0 and qx~0.
            The change from theta to alpha seems redondant, but allows a better understanding of the GISAXS-like geometry.
            Note that the data saved in the .moy/.mat files are NOT the ones used in the true GIXD plots.
            Most of the time the qy~0 is valid, but one should check that it is the case once per geometry used.
            """
            # Plot true GIXD
            twotheta = np.array(data[column_delta][istart:istop])*np.pi/180.
            xx, yy = np.meshgrid(twotheta, ch)

            # Conversion factor alphaf - channel (in rad/chan)
            alpha_factor = thetazfactor
            # Critical angle in rad
            alphac = thetac
            # Channel of Vineyard's peak
            channelc = channel0
            # alphaf (exit angle towards the detector)
            alphaf = alphac + mean_gamma*np.pi/180. + (channelc-yy)*alpha_factor
            # alphai (incident angle)
            alphai = thetai

            # True qx, qy, qz in nm^-1
            k0 = 2*np.pi/wavelength
            qx = k0*(np.cos(alphaf)*np.cos(twotheta)-np.cos(alphai))
            qy = k0*np.cos(alphaf)*np.sin(twotheta)
            qz = k0*(np.sin(alphaf)+np.sin(alphai))
            qxy = np.sqrt(np.square(qx)+np.square(qy))
            #q = np.sqrt(np.square(qxy)+np.square(qz))

            fig=plt.figure(figsize=(12,5))
            fig.subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.16)
            ax1=fig.add_subplot(121)
            ax1.pcolormesh(twotheta*180./np.pi, alphaf*180./np.pi, mat.transpose(), cmap=cmap, rasterized=True)
            ax1.set_xlabel('2 theta (deg)', fontsize='large')
            ax1.set_ylabel('alpha_f (deg)', fontsize='large')

            ax2=fig.add_subplot(122)
            ax2.pcolormesh(qxy, qz, mat.transpose(), cmap=cmap, rasterized=True)
            ax2.set_xlabel('qxy (nm^-1)', fontsize='large')
            ax2.set_ylabel('qz (nm^-1)', fontsize='large')
            fig.suptitle('True GIXD', fontsize='x-large')

            if column_pi is not None:
                fig.text(.04, .05, r'$\pi = %3.4gmN.m^{-1}$'%(mean_pi),
                         fontsize='large', color='red')
            if column_gamma is not None:
                fig.text(.96, .05, r'$\gamma = %3.4g deg$'%(mean_gamma),
                         fontsize='large', color='red', horizontalalignment='right')


            if verbose:
                print('\t. For more details on the geometry, see:')
                print('\t \t -Fig.2 in doi:10.1107/S0909049512022017')
                print('\t \t -Slide 4 in http://gisaxs.com/files/Strzalka.pdf')

        # Create Save Name
        savename=working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]+'_1D'

        # Save the original matrix
        np.savetxt(savename+'.mat', original_mat)
        if verbose: print('\t. Original, non binned matrix saved in:')
        if verbose: print("\t", savename+'.mat')

        # Take care of scalar data
        tosave=np.zeros((daty.shape[0], 4), float)
        tosave[:,0]=daty
        tosave[:,1]=datyTop
        tosave[:,2]=datyBottom
        tosave[:,3]=datyFirstQuarter
        #print(np.array(data[istart:istop]).transpose().shape, tosave.shape)
        data=np.array(data).transpose()
        data=data[istart:istop,:]
        new_data=np.concatenate((data,tosave), axis=1)
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
        for i in range(new_data.shape[0]):
            s=""
            for j in range(new_data.shape[1]):
                s=s+str(new_data[i,j])+'\t'
            f.write(s+'\n')
        f.close()
        if verbose: print('\t. Scalar data saved in:')
        if verbose: print("\t", savename+'.dat')

        # Save as moy
        for binsize in moytocreate:
            # Bin the matrix
            ch, mat= Groupe(original_mat, binsize=binsize)

            # Compute qz
            if computeqz==True:
                thetaz=thetac+mean_gamma*np.pi/180.0+(channel0-ch)*thetazfactor
                qz=2.0*np.pi*np.sin(thetaz)/wavelength
            else:
                qz=ch
            moy=np.zeros((nbpts*qz.shape[0],3), float)
            index=0
            for i in range(data.shape[0]):
                for j in range(qz.shape[0]):
                    moy[index,0]=datx[i]
                    moy[index,1]=qz[j]
                    moy[index,2]=mat[i,j]
                    index=index+1
            f=open(savename+'.moy'+str(binsize), 'w')
            f.write(stamps0D[columnx][1]+'\t'+'Qz \t Intensity\n')
            for i in range(moy.shape[0]):
                f.write(str(moy[i,0])+'\t'+str(moy[i,1])+'\t'+str(moy[i,2])+'\n')
            f.close()
            # Save the matrix
            np.savetxt(savename+'.mat'+str(binsize), mat)
            # Save the Qz
            if computeqz:
                np.savetxt(savename+'_qz'+str(binsize)+'.dat'+str(binsize), qz)
                if verbose: print('\t. Qz values saved in:')
                if verbose: print('\t'+savename+'_qz'+str(binsize)+'.dat')
            if verbose:
                print('\t. Binned matrix saved in:')
                print("\t", savename+'.mat'+str(binsize))

                print('\t. XYZ    data saved in:')
                print("\t", savename+'.moy'+str(binsize))
            
        plt.show()
        
        
def Plot_isotherm(nxs_filename='SIRIUS_test.nxs', working_dir = '', recording_dir='',
                  show_data_stamps=False, verbose=False, fast=True):
    
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
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
            return
        nbpts=np.int(nexus.get_nbpts())
        if verbose: print("\t. Number of data points: ", nbpts)
            
        # Get stamps and Data
        stamps, data=nexus.extractData('0D')
        
        # Save data
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
            savename = working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]
            print('\t. Data saved in:')
            print("\t", savename+'.dat')
                        
        
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
            if verbose: print('\t. Area per molecule found column %d'%(column_area))
            if verbose: print('\t. Surface pressure per molecule found column %d'%(column_pi))
            if verbose: print('\t. Time per molecule found column %d'%(column_time))
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

            if verbose: print('\t. Valid data between points %d and %d'%(istart, istop))

            # Create Graph

            fig=plt.figure(1, figsize=(12,5))
            fig.subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.16)
            ax1=fig.add_subplot(121)
            ax2=fig.add_subplot(122)
            ax3=ax2.twinx()

            ax1.plot(data[column_area], data[column_pi], 'r-', lw=2)
            ax2.plot(data[column_time], data[column_area], 'k-', lw=2)
            ax3.plot(data[column_time], data[column_pi], 'b-', lw=2)

            
            S=nxs_filename[nxs_filename.rfind('/')+1:nxs_filename.rfind('.')]
            
            fig.suptitle(S, fontsize=20) 

            ax1.set_xlabel('Area per Molecule (nm$\mathregular{^2}$)', fontsize=14)
            ax1.set_ylabel('Surface pressure (mN/m)', fontsize=14)

            ax2.set_xlabel('Time (sec)', fontsize=14)
            ax2.set_ylabel('Area per Molecule (nm$\mathregular{^2}$)', fontsize=14)
            ax3.set_ylabel('Surface pressure (mN/m)', fontsize=14, color='b')
        
        

            
##########################################################################################
####################################### GIXS #############################################
##########################################################################################

            
def Extract_GIXS(nxs_filename='SIRIUS_test.nxs', working_dir='', recording_dir='',
                 logz=True, wavelength=0.155, thetai=0.002, distance=2722,
                 pixel_PONI_x=490, pixel_PONI_y=975, pixel_size=0.172,
                 number_bins_x=10, number_bins_y=10, xmin=0., xmax=1., ymin=0., ymax=1.,
                 show_data_stamps=False, force_gamma_delta=False, fgamma=0., fdelta=0.,
                 verbose=False, absorbers='', cmap='viridis',
                 plot_twotheta_alphaf=False, plot_qxy_qz=False, plot_qxy_q=False):
    
    """
    Extract, plot, and save the GIXS scan. 
    """
    
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        i_pilatus=None
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
                    i_pilatus=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])
        
        # Get absorbers
        if absorbers != '':
            print("\t. Absorbers:", str(absorbers))
                    
        # Check that Pilatus data are present (images)
        if i_pilatus is not None:
            if verbose: print('\t. Pilatus data found, (column %d, alias %s)'%(i_pilatus, stamps[i_pilatus][1]))
        else:
            print(PN._RED,'\t. No pilatus data found', PN._RESET)
            nexus.close()
            return

        images = np.zeros([nbpts, 1043, 981])

        for i in range(nbpts):
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i+1, nbpts))
            sys.stdout.flush()

            #Extract the images from the nexus file
            stamp, image = nexus.extract_one_data_point(stamps[i_pilatus][0], i, verbose = False)

            #Remove the dead pixels
            for ii,jj in dead_pixels:
                image[ii,jj]=0.0   

            images[i,:] = image    
            
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()

            
        images_sum = images.sum(axis=0)

        #Replace the dead zone (spacing between chips) on the detector with -2. 
        images_sum = np.where(images_sum>0, images_sum, -2.)
        
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
            if verbose: print('\t. Gamma motor data found, mean value %3.4g deg'%(gamma))
        else:
            gamma = fgamma
            print(PN._RED,'\t. No gamma found! gamma = %g'%gamma, PN._RESET)
        if (i_delta != None) and (force_gamma_delta==False):
            delta = np.mean(data[i_delta])
            if verbose: print('\t. Delta motor data found, mean value %3.4g deg'%(delta))
        else:
            delta = fdelta
            print(PN._RED,'\t. No delta found! delta = %g'%delta, PN._RESET)

        if verbose: 
            print('\t. For more details on the geometry, see:')
            print('\t \t -Fig.2 in doi:10.1107/S0909049512022017')
            print('\t \t -Slide 4 in http://gisaxs.com/files/Strzalka.pdf')

        pixels_x = np.arange(0,np.shape(images_sum)[1],1)
        pixels_y = np.arange(0,np.shape(images_sum)[0],1)

        xx, yy = np.meshgrid(pixels_x, pixels_y)

        # alphai (incident angle)
        alphai = thetai    

        pixel_direct_x = pixel_PONI_x-distance/pixel_size*np.tan(delta*np.pi/180.)
        pixel_direct_y = pixel_PONI_y-distance/pixel_size*np.tan(gamma*np.pi/180.)

        # 2*theta in rad
        twotheta = np.arctan(pixel_size*(xx-pixel_direct_x)/distance)

        # alpha_f in rad
        deltay0 = distance*np.tan(alphai*np.pi/180.)
        alphaf = np.arctan( (pixel_size*(pixel_direct_y-yy)-deltay0)/distance)
        
        # True qx, qy, qz in nm^-1
        k0 = 2*np.pi/wavelength
        qx = k0*(np.cos(alphaf)*np.cos(twotheta)-np.cos(alphai))
        qy = k0*np.cos(alphaf)*np.sin(twotheta)
        qz = k0*(np.sin(alphaf)+np.sin(alphai))
        qxy = np.sqrt(np.square(qx)+np.square(qy))
        q = np.sqrt(np.square(qxy)+np.square(qz))
              
        fig = plt.figure(figsize=(15,15))
        fig.subplots_adjust(top=0.95)
        fig.suptitle(nxs_filename.split('\\')[-1], fontsize='x-large')
        
        #Divide the grid in 2x2
        outer = gridspec.GridSpec(2, 2, wspace=0.2)

        #Divide the left row in 2x1
        inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[0], hspace=0.5)
        
        if plot_twotheta_alphaf:
            
            #Plot a profile along y (integrated over x)
            ax0 = fig.add_subplot(inner[0])   
            alphaf_bins = np.linspace(0.,np.max(alphaf),number_bins_y)

            twotheta_binned = np.zeros(len(alphaf_bins)-1)
            alphaf_binned = np.zeros(len(alphaf_bins)-1)
            alphaf_flat = np.ravel(alphaf)
            twotheta_flat = np.ravel(images_sum)
            for n in range(len(alphaf_bins)-1):
                alphaf_bin_inf = alphaf_bins[n]
                alphaf_bin_sup = alphaf_bins[n+1]
                twotheta_binned[n] = np.sum(twotheta_flat[(alphaf_flat<alphaf_bin_sup) & (alphaf_bin_inf<alphaf_flat)])
                alphaf_binned[n] = (alphaf_bin_sup+alphaf_bin_inf)/2.                 
                
            if logz: ax0.set_yscale('log') 
                
            # Define lims of the plot
            ax0.set_xlim(ymin,ymax)         
            temp = twotheta_binned[(alphaf_binned*180./np.pi<ymax) & (ymin<alphaf_binned*180./np.pi)]
            ymin_plot = np.min(temp[temp>0])   
            ymax_plot = np.max(twotheta_binned[(alphaf_binned*180./np.pi<ymax) & (ymin<alphaf_binned*180./np.pi)])
            ax0.set_ylim(0.8*ymin_plot*180./np.pi,1.2*ymax_plot*180./np.pi)       
            
            ax0.set_xlabel('alpha_f (deg)', fontsize='large')
            ax0.set_ylabel('2 theta (deg)', fontsize='large')
            ax0.plot(alphaf_binned*180./np.pi, twotheta_binned*180./np.pi)            
            
            
            #Plot a profile along x (integrated over y)
            ax1 = fig.add_subplot(inner[1])
       
            twotheta_bins = np.linspace(0.,np.max(twotheta),number_bins_x)

            alphaf_binned = np.zeros(len(twotheta_bins)-1)
            twotheta_binned = np.zeros(len(twotheta_bins)-1)
            twotheta_flat = np.ravel(twotheta)
            alphaf_flat = np.ravel(images_sum)

            for n in range(len(twotheta_bins)-1):
                twotheta_bin_inf = twotheta_bins[n]
                twotheta_bin_sup = twotheta_bins[n+1]
                alphaf_binned[n] = np.sum(alphaf_flat[(twotheta_flat<twotheta_bin_sup) & (twotheta_bin_inf<twotheta_flat)])
                twotheta_binned[n] = (twotheta_bin_sup+twotheta_bin_inf)/2.  
                
            if logz: ax1.set_yscale('log')
                
            # Define lims of the plot
            ax1.set_xlim(xmin,xmax)          
            temp = alphaf_binned[(twotheta_binned*180./np.pi<xmax) & (xmin<twotheta_binned*180./np.pi)]
            ymin_plot = np.min(temp[temp>0])   
            ymax_plot = np.max(alphaf_binned[(twotheta_binned*180./np.pi<xmax) & (xmin<twotheta_binned*180./np.pi)])
            ax1.set_ylim(0.8*ymin_plot*180./np.pi,1.2*ymax_plot*180./np.pi)                       
                                          
                
            ax1.set_xlabel('2 theta (deg)', fontsize='large')
            ax1.set_ylabel('alpha_f (deg)', fontsize='large')
            ax1.plot(twotheta_binned*180./np.pi, alphaf_binned*180./np.pi)

            #Divide the right row in 1x1
            inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                            subplot_spec=outer[1], wspace=0.1, hspace=0.1)

            #Show the full image integrated over the scan
            ax2 = fig.add_subplot(inner[0])

            ax2.pcolormesh(twotheta*180./np.pi, alphaf*180./np.pi, images_sum, cmap = cmap,
                           norm = colors.LogNorm(), rasterized=True)
            ax2.set_xlabel('2 theta (deg)', fontsize='large')
            ax2.set_ylabel('alpha_f (deg)', fontsize='large')


        if plot_qxy_qz:
            
            #Plot a profile along y (integrated over x)
            ax0 = fig.add_subplot(inner[0])
            
            qz_bins = np.linspace(0.,np.max(qz),number_bins_y)

            qxy_binned = np.zeros(len(qz_bins)-1)
            qz_binned = np.zeros(len(qz_bins)-1)
            qz_flat = np.ravel(qz)
            qxy_flat = np.ravel(images_sum)
            for n in range(len(qz_bins)-1):
                qz_bin_inf = qz_bins[n]
                qz_bin_sup = qz_bins[n+1]
                qxy_binned[n] = np.sum(qxy_flat[(qz_flat<qz_bin_sup) & (qz_bin_inf<qz_flat)])
                qz_binned[n] = (qz_bin_sup+qz_bin_inf)/2.                 
                
            if logz: ax0.set_yscale('log') 

            # Define lims of the plot
            ax0.set_xlim(ymin,ymax)
            
            temp = qxy_binned[(qz_binned<ymax) & (ymin<qz_binned)]
            ymin_plot = np.min(temp[temp>0])   
            ymax_plot = np.max(qxy_binned[(qz_binned<ymax) & (ymin<qz_binned)])
            ax0.set_ylim(0.8*ymin_plot,1.2*ymax_plot)       
            
            ax0.set_xlabel('qz (nm^-1)', fontsize='large')
            ax0.set_ylabel('qxy (nm^-1)', fontsize='large')
            ax0.plot(qz_binned, qxy_binned)            
            
            #Plot a profile along x (integrated over y)
            ax1 = fig.add_subplot(inner[1])
     
            qxy_bins = np.linspace(0.,np.max(qxy),number_bins_x)

            qz_binned = np.zeros(len(qxy_bins)-1)
            qxy_binned = np.zeros(len(qxy_bins)-1)
            qxy_flat = np.ravel(qxy)
            qz_flat = np.ravel(images_sum)

            for n in range(len(qxy_bins)-1):
                qxy_bin_inf = qxy_bins[n]
                qxy_bin_sup = qxy_bins[n+1]
                qz_binned[n] = np.sum(qz_flat[(qxy_flat<qxy_bin_sup) & (qxy_bin_inf<qxy_flat)])
                qxy_binned[n] = (qxy_bin_sup+qxy_bin_inf)/2.  
                
            if logz: ax1.set_yscale('log')

            # Define lims of the plot
            ax1.set_xlim(xmin,xmax)          
            temp = qz_binned[(qxy_binned<xmax) & (xmin<qxy_binned)]
            ymin_plot = np.min(temp[temp>0])   
            ymax_plot = np.max(qz_binned[(qxy_binned<xmax) & (xmin<qxy_binned)])
            ax1.set_ylim(0.8*ymin_plot,1.2*ymax_plot)                              
 
            ax1.set_xlabel('qxy (nm^-1)', fontsize='large')
            ax1.set_ylabel('qz (nm^-1)', fontsize='large')
            ax1.plot(qxy_binned, qz_binned)

            #Divide the right row in 1x1
            inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                            subplot_spec=outer[1], wspace=0.1, hspace=0.1)

            #Show the full image integrated over the scan
            ax2 = fig.add_subplot(inner[0])

            ax2.pcolormesh(qxy, qz, images_sum, cmap = cmap,
                           norm = colors.LogNorm(), rasterized=True)
            ax2.set_xlabel('qxy (nm^-1)', fontsize='large')
            ax2.set_ylabel('qz (nm^-1)', fontsize='large')

 
        if plot_qxy_q:
            
            #Plot a profile along y (integrated over x)
            ax0 = fig.add_subplot(inner[0])

            q_bins = np.linspace(0.,np.max(q),number_bins_y)

            qxy_binned = np.zeros(len(q_bins)-1)
            q_binned = np.zeros(len(q_bins)-1)
            q_flat = np.ravel(q)
            qxy_flat = np.ravel(images_sum)
            for n in range(len(q_bins)-1):
                q_bin_inf = q_bins[n]
                q_bin_sup = q_bins[n+1]
                qxy_binned[n] = np.sum(qxy_flat[(q_flat<q_bin_sup) & (q_bin_inf<q_flat)])
                q_binned[n] = (q_bin_sup+q_bin_inf)/2.                 
                
            if logz: ax0.set_yscale('log') 

            # Define lims of the plot
            ax0.set_xlim(ymin,ymax)
            
            temp = qxy_binned[(q_binned<ymax) & (ymin<q_binned)]
            ymin_plot = np.min(temp[temp>0])   
            ymax_plot = np.max(qxy_binned[(q_binned<ymax) & (ymin<q_binned)])
            ax0.set_ylim(0.8*ymin_plot,1.2*ymax_plot)                                 
             
            ax0.set_xlabel('q (nm^-1)', fontsize='large')
            ax0.set_ylabel('qxy (nm^-1)', fontsize='large')
            ax0.plot(q_binned, qxy_binned)            
            
            #Plot a profile along x (integrated over y)
            ax1 = fig.add_subplot(inner[1])
       
            qxy_bins = np.linspace(0.,np.max(qxy),number_bins_x)

            q_binned = np.zeros(len(qxy_bins)-1)
            qxy_binned = np.zeros(len(qxy_bins)-1)
            qxy_flat = np.ravel(qxy)
            q_flat = np.ravel(images_sum)
            
            for n in range(len(qxy_bins)-1):
                qxy_bin_inf = qxy_bins[n]
                qxy_bin_sup = qxy_bins[n+1]
                q_binned[n] = np.sum(q_flat[(qxy_flat<qxy_bin_sup) & (qxy_bin_inf<qxy_flat)])
                qxy_binned[n] = (qxy_bin_sup+qxy_bin_inf)/2.  
                
            if logz: ax1.set_yscale('log')
                
            # Define lims of the plot
            ax1.set_xlim(xmin,xmax)          
            temp = q_binned[(qxy_binned<xmax) & (xmin<qxy_binned)]
            ymin_plot = np.min(temp[temp>0])   
            ymax_plot = np.max(q_binned[(qxy_binned<xmax) & (xmin<qxy_binned)])
            ax1.set_ylim(0.8*ymin_plot,1.2*ymax_plot)             
                
            ax1.set_xlabel('qxy (nm^-1)', fontsize='large')
            ax1.set_ylabel('q (nm^-1)', fontsize='large')
            ax1.plot(qxy_binned, q_binned)

            #Divide the right row in 1x1
            inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                            subplot_spec=outer[1], wspace=0.1, hspace=0.1)

            #Show the full image integrated over the scan
            ax2 = fig.add_subplot(inner[0])

            ax2.pcolormesh(qxy, q, images_sum, cmap = cmap,
                           norm = colors.LogNorm(), rasterized=True)
            ax2.set_xlabel('qxy (nm^-1)', fontsize='large')
            ax2.set_ylabel('q (nm^-1)', fontsize='large')


        plt.show()
        plt.close()    
        
        # Create Save Name
        savename=working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]
            
        np.savetxt(savename+'_pilatus_sum.mat', images_sum)

        im = Image.fromarray(images_sum)
        im.save(savename+'_pilatus_sum.tiff')

        if verbose: 
            print('\t. Original matrix saved in:')
            print("\t", savename+'.mat')
            print(" ")
            print('\t. Tiff saved in:')
            print("\t", savename+'.tiff')
            print(" ")
            
       
  
##########################################################################################
###################################### XRF ###############################################
##########################################################################################

def Extract_XRF(nxs_filename='SIRIUS_test.nxs', working_dir='', recording_dir='',
                logz=False, list_elems=[0,1,2,3], first_channel=0, last_channel=2048,
                use_eV=False, gain=10., eV0=0.,
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
                ax1.pcolormesh(xx, yy, spectrums.transpose(), cmap='viridis', norm = colors.LogNorm(), rasterized=True)
            else:
                ax1.pcolormesh(xx, yy, spectrums.transpose(), cmap='viridis', rasterized=True)

            plt.show()

        if plot_sum:
            fig = plt.figure(figsize=(12,4.5))
            ax1 = fig.add_subplot(111)
            ax1.set_ylabel('counts', fontsize='large')
            if logz: ax1.set_yscale('log')
            if use_eV:
                ax1.set_xlabel('eV', fontsize='large')
                ax1.plot(eVs, np.sum(spectrums, axis = 0), 'b.-', label='Sum of spectrums')
            else:
                ax1.set_xlabel('channel', fontsize='large')
                ax1.plot(channels, np.sum(spectrums, axis = 0), 'b.-', label='Sum of spectrums')  
            ax1.legend(fontsize='large')
            plt.show()

        if plot_first_last:    
            #Plot the selected channel range
            fig = plt.figure(figsize=(12,4.5))
            ax1 = fig.add_subplot(111)
            ax1.set_ylabel('counts', fontsize='large')
            if logz: ax1.set_yscale('log')
            if use_eV:
                ax1.set_xlabel('eV', fontsize='large')        
                ax1.plot(eVs, spectrums[first_non_zero_spectrum], 'b.-', label='First spectrum')
                ax1.plot(eVs, spectrums[-1], 'r.-', label='Last spectrum')            
            else:
                ax1.set_xlabel('channel', fontsize='large')        
                ax1.plot(channels, spectrums[first_non_zero_spectrum], 'b.-', label='First spectrum')
                ax1.plot(channels, spectrums[-1], 'r.-', label='Last spectrum')
            ax1.legend(fontsize='large')
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

    
def Extract_pilatus_sum(nxs_filename='SIRIUS_test.nxs', working_dir='', recording_dir='', logz=True,
                        xmin=0, xmax=980, ymin=0, ymax=1042,
                        show_data_stamps=False, verbose=False, absorbers='', fast=True, cmap='viridis'):
    
    """
    Extract, plot, and save the sum of all the pilatus images in a scan. 
    """
    
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        i_pilatus=None
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
        stamps=nexus.extractStamps()
        if show_data_stamps : print("\t. Available Counters:")
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                if show_data_stamps : print("\t\t", i, ' -------> ', stamps[i][1])
                if stamps[i][1].lower()=='pilatus':
                    i_pilatus=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])

        # Get absorbers
        if absorbers != '':
            print("\t. Absorbers:", str(absorbers))                    
                       
        # Check that Pilatus data are present (images)
        if i_pilatus is not None:
            if verbose: print('\t. Pilatus data found, (column %d, alias %s)'%(i_pilatus, stamps[i_pilatus][1]))
        else:
            print(PN._RED,'\t. No pilatus data found', PN._RESET)
            nexus.close()
            return

        images = np.zeros([nbpts, 1043, 981])

        for i in range(nbpts):
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i+1, nbpts))
            sys.stdout.flush()

            #Extract the images from the nexus file
            stamp, image = nexus.extract_one_data_point(stamps[i_pilatus][0], i, verbose = False)

            #Remove the dead pixels
            for ii,jj in dead_pixels:
                image[ii,jj]=0.0   

            images[i,:] = image    
            
            sys.stdout.write('                                                                                  \r')
            sys.stdout.flush()

        nexus.close()  
        
        images_sum = images.sum(axis=0)

        #Replace the dead zone (spacing between chips) on the detector with -2. 
        images_sum = np.where(images_sum>0, images_sum, -2.)

        pixels_x = np.arange(0,np.shape(images_sum)[1],1)
        pixels_y = np.arange(0,np.shape(images_sum)[0],1)

        xx, yy = np.meshgrid(pixels_x, pixels_y)

              
        fig = plt.figure(figsize=(15,15))

        #Divide the grid in 2x2
        outer = gridspec.GridSpec(2, 2, wspace=0.2)

        #Divide the left row in 2x1
        inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[0], hspace=0.5)

        #Plot a profile along y (integrated over x)
        ax1 = fig.add_subplot(inner[0])
        profile_y = images_sum.sum(axis=1)
        ax1.set_xlim(ymin,ymax)
        
        temp = profile_y[int(ymin):int(ymax)]
        ax1.set_ylim(np.min(temp[temp>0])*0.8,np.max(profile_y[int(ymin):int(ymax)])*1.2)
        if logz: ax1.set_yscale('log')
        ax1.set(xlabel = 'vertical pixel (y)', ylabel = 'integration along x')
        ax1.plot(profile_y)

        #Plot a profile along x (integrated over y)
        ax2 = fig.add_subplot(inner[1])
        profile_x = images_sum.sum(axis=0)
        ax2.set_xlim(xmin,xmax)
        temp = profile_x[int(xmin):int(xmax)]
        ax2.set_ylim(np.min(temp[temp>0])*0.8,np.max(profile_x[int(xmin):int(xmax)])*1.2)
        if logz: ax2.set_yscale('log')
        ax2.set(xlabel = 'horizontal pixel (x)', ylabel = 'integration along y')
        ax2.plot(profile_x)

        #Divide the right row in 1x1
        inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                        subplot_spec=outer[1], wspace=0.1, hspace=0.1)

        #Show the full image integrated over the scan
        ax0 = fig.add_subplot(inner[0])
        if logz:
            ax0.pcolormesh(xx, yy, images_sum, cmap = cmap, norm = colors.LogNorm(), rasterized=True)
        else:
            ax0.pcolormesh(xx, yy, images_sum, cmap = cmap, rasterized=True)
        ax0.set(xlabel = 'horizontal pixel (x)', ylabel ='vertical pixel (y)')
        ax0.invert_yaxis()
        fig.subplots_adjust(top=0.95)
        fig.suptitle(nxs_filename.split('\\')[-1], fontsize=16)
        plt.show()
        plt.close()

        # Create Save Name
        savename=working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]
        
        # profile y
        np.savetxt(savename+'_profile_y.dat', profile_y)
        
        # profile x
        np.savetxt(savename+'_profile_x.dat', profile_x)
        
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
        