import numpy as np
import lmfit as lm

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