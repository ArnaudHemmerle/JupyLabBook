# JupyLabBook
A Jupyter Notebook used as an interactive lab book on the beamline SIRIUS (SOLEIL Synchrotron).

### Last versions of modules:  
FrontendFunctions.py: 1.0.1    
CustomFunctions.py: 1.0.3    

# User manual

See the file JupyLabBook_User_Manual.pdf in the repository.


# What is JupyLabBook?
JupyLabBook is an interactive lab book to be used on SIRIUS. It is meant to be as close as possible as the paper-like lab book that you would use on the beamline (except for the ruler and the tape).  
The user can easily plot raw data, do data reduction, and write sample description without typing a single line of code. A PDF with a linear presentation of the notebook can be generated at any time of the experiment.  
Expert users can define their own functions to perform more involved data analysis or customed presentation of the data, if needed.

# What is not JupyLabBook?
JupyLabBook was not designed to be a notebook for **analysis**, but for **data reduction** only.  
Think of it as the traditional lab book, that should remain untouched after the end of the experiment. Notebooks for specific analysis (GIXS, GIXD, XRF ...) can be provided as well on demand.
JupyLabBook is not made of paper: it can crash, it can be corrupted, it can be deleted by mistake ... In the present state of development, we advise users to keep a written paper-like lab book in parallel with any interactive notebook. Especially for all the info concerning the samples used, that cannot be retrieved afterwards.
