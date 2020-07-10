# JupyLabBook
A Jupyter Notebook used as a lab book on the beamline SIRIUS (SOLEIL Synchrotron).

# Current status
JupyLabBook is presently in development & testing.  

### Last versions of modules:  
FrontendFunctions.py: 0.18  
CustomFunctions.py: 0.18  


# What is JupyLabBook?
JupyLabBook is an interactive lab book to be used on SIRIUS. It is meant to be as close as possible as the paper-like lab book that you would use on the beamline (except for the ruler and the tape).  
The user can easily plot raw data, do data reduction, and write sample description without typing a single line of code. A PDF with a linear presentation of the notebook can be generated at any time of the experiment.  
Expert users can define their own functions to perform more involved data analysis or customed presentation of the data, if needed.

# What is not JupyLabBook?
JupyLabBook was not designed to be a notebook for **analysis** but for **data reduction** only.  
Think of it as the traditional lab book, that should remain untouched after the end of the experiment. Notebooks for specific analysis (GIXS, GIXD, XRF ...) can be provided as well on demand.
JupyLabBook is not made of paper: it can crash, it can be corrupted, it can be deleted by mistake ... We advice users to always keep a written paper-like lab book in parallel with any interactive notebook. Especially for all the info concerning the samples used, that cannot be retrieved afterwards.

# How to use JupyLabBook?
1) Contrary to its name, JupyLabBook cannot be used with Jupyter Lab. It has to be used with Jupyter Notebook. If you are running Jupyter Lab, use the option "Launch Classic Notebook" in the Help menu.  

2) After short introduction by the beamline staff, you should be able to generate the plots by clicking on the right buttons.
 
# Installation notes
Recommandations for installation on Windows 10:
1) Install Anaconda Python >= 3.7 (conda -version >= 4.8.3)  
2) Install MikTex >= 2.9.7. Restart after installation.  
3) Install the package lmfit by typing in the conda prompt:  
```conda install -c conda-forge lmfit```



 
