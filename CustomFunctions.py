import PyNexus as PN
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import time
import subprocess
import sys
import lmfit as L
from scipy.special import erf


__version__ = '0.1'

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

        # Create the Qxy Qz map from images
        daty=[]
        datytop=[]
        datybottom=[]
        datyFirstQuarter=[]
        mat=[]
        
        for i in range(istart, istop, 1):
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i, nbpts))
            sys.stdout.flush()
            image=images[0][i]
            # dead pixels
            image[877,528]=0.0
            image[922,847]=0.0
            image[1018,881]=0.0
            image[382,432]=0.0
            image[640,859]=0.0
            image[640,860]=0.0
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
        # The graph has to be splitted into for good rendering in PDF
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
        
        print('Data not saved. To save data, run a GIXD on the scan.')
        print('Channel0: %g'%i_max)
        
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
                 show_data_stamps=False, show_saved_data=False, plot_true_GIXD=False):
    
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
        if show_data_stamps : print("\t. Available Counters:")
        for i in range(len(stamps)):
            if stamps[i][1] is not None:
                if show_data_stamps : print("\t\t", i, ' -------> ', stamps[i][1])
                if stamps[i][1].lower()=='pilatus':
                    columnz=i
            else:
                if show_data_stamps : print("\t\t",i, ' -------> ', stamps[i][0])
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
                return
            else:
                columnx=column_delta
                print('\t. delta data found, (column %d, alias %s)'%(columnx, stamps[columnx][1]))
        else:
            columnx=column_qxy
            print('\t. qxy data found, (column %d, alias %s)'%(columnx, stamps[columnx][1]))

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

        print('\t. Valuable data between points %d and %d'%(istart, istop))

        # Compute and print mean values
        if column_pi is not None:
            mean_pi=data[column_pi][istart:istop-1].mean()
            print('\t. Surface pressure data found, mean value %3.4g ± %3.4g mN/m'%(mean_pi,data[column_pi][istart:istop-1].std() ))
        else:
            print('\t. No surface pressure data found')
        if column_area is not None:
            mean_area=data[column_area][istart:istop-1].mean()
            print('\t. Area per molecule data found, mean value %3.4g ± %3.4g nm2 per molecule'%(mean_area, data[column_area][istart:istop-1].std()))
        else:
            print('\t. No area per molecule data found')
        if column_gamma is not None:
            mean_gamma=data[column_gamma][istart:istop-1].mean()
            print('\t. Gamma motor data found, mean value %3.4g deg'%(mean_gamma))
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
            sys.stdout.write('Treat image %d/%d                                                                 \r'%(i, nbpts))
            sys.stdout.flush()
            image=images[0][i]
            # dead pixels
            image[877,528]=0.0
            image[922,847]=0.0
            image[1018,881]=0.0
            image[382,432]=0.0
            image[640,859]=0.0
            image[640,860]=0.0
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
            ax2.contourf(datx[istart:istop], qz, (mat.transpose()),levels=np.linspace(zmin, zmax, nblevels), cmap='jet')

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
            ax1.pcolormesh(twotheta*180./np.pi, alphaf*180./np.pi, mat.transpose(), cmap = 'jet')
            ax1.set_xlabel('2 theta (deg)', fontsize='large')
            ax1.set_ylabel('alpha_f (deg)', fontsize='large')

            ax2=fig.add_subplot(122)
            ax2.pcolormesh(qxy, qz, mat.transpose(), cmap = 'jet')
            ax2.set_xlabel('qxy (nm^-1)', fontsize='large')
            ax2.set_ylabel('qz (nm^-1)', fontsize='large')
            fig.suptitle('True GIXD', fontsize='x-large')

            if column_pi is not None:
                fig.text(.04, .05, r'$\pi = %3.4gmN.m^{-1}$'%(mean_pi),
                         fontsize='large', color='red')
            if column_gamma is not None:
                fig.text(.96, .05, r'$\gamma = %3.4g deg$'%(mean_gamma),
                         fontsize='large', color='red', horizontalalignment='right')


            print('\t. For more details on the geometry, see:')
            print('\t \t -Fig.2 in doi:10.1107/S0909049512022017')
            print('\t \t -Slide 4 in http://gisaxs.com/files/Strzalka.pdf')

        # Create Save Name
        savename=working_dir+nxs_filename[:nxs_filename.rfind('.nxs')]+'_1D'

        # Save the original matrix
        np.savetxt(savename+'.mat', original_mat)
        if show_saved_data: print('\t. Original, non binned matrix saved in:')
        if show_saved_data: print("\t\t", savename+'.mat')

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
        if show_saved_data: print('\t. Scalar data saved in:')
        if show_saved_data: print("\t\t", savename+'.dat')

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
                if show_saved_data: print('\t. Qz values saved in:')
                if show_saved_data: print('\t\t'+savename+'_qz'+str(binsize)+'.dat')
            if show_saved_data:
                print('\t. Binned matrix saved in:')
                print("\t\t", savename+'.mat'+str(binsize))

                print('\t. XYZ    data saved in:')
                print("\t\t", savename+'.moy'+str(binsize))
        if not show_saved_data:
            print('\t. Data saved in text format')
            
        plt.show()
        
        
def Plot_isotherm(nxs_filename='SIRIUS_test.nxs', recording_dir='', show_data_stamps=False):
    
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(PN._RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+PN._RESET)
        print(('\t\t recording directory : '+recording_dir))
    else:
        column_pi=None
        column_area=None
        column_time=None
        print(PN._BLUE+" - Open Nexus Data File :"+ PN._RESET)
        print('\t'+nxs_path)
        try:
            nexus=PN.PyNexusFile(nxs_path)
        except:
            print(PN._RED,'\t Nexus file seems not to exist or is not correct',PN._RESET)
            return
        nbpts=np.int(nexus.get_nbpts())
        print("\t. Number of data points: ", nbpts)
        # Get stamps and Data
        stamps, data=nexus.extractData('0D')
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

            print('\t. Valuable data between points %d and %d'%(istart, istop))

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
###################################### GENERAL ###########################################
##########################################################################################

def Plot_1D(nxs_filename='SIRIUS_test.nxs', recording_dir='',
            xLabel='zs', yLabel='alphax'):
    """
    Simple 1D plot. When called by FrontendFunctions, will plot the current selection in the interactive 1D plot.
    It is necessary to call this function to have the plot in the pdf.
    """
    nexus = PN.PyNexusFile(recording_dir+nxs_filename, fast=True)
    stamps0D, data0D = nexus.extractData('0D')
    sensor_list = [stamps0D[i][0] if stamps0D[i][1]== None else stamps0D[i][1] for i in range(len(stamps0D))]

    xArg = sensor_list.index(xLabel)
    yArg = sensor_list.index(yLabel)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(data0D[xArg], data0D[yArg], 'o-')
    ax.set_xlabel(stamps0D[xArg][1], fontsize=16)
    ax.set_ylabel(stamps0D[yArg][1], fontsize=16)
    plt.show()
