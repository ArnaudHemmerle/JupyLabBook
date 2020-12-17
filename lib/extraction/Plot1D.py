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

#REFAIRE BIEN

def Treat(nxs_filename, recording_dir, working_dir, 
          xLabel, yLabel):
    """
    Simple 1D plot. When called by FrontendFunctions by clicking on Treat Scan, 
    it will plot the current selection in the interactive 1D plot.
    It is necessary to call this function to have the plot in the pdf.
    """
    
    # Extract
    nexus = PN.PyNexusFile(recording_dir+nxs_filename, fast=True)
    stamps0D, data0D = nexus.extractData('0D')
   
    sensor_list = [stamps0D[i][0] if stamps0D[i][1]== None else stamps0D[i][1] for i in range(len(stamps0D))]
    
    xArg = sensor_list.index(xLabel)
    yArg = sensor_list.index(yLabel)

    # Plot
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(data0D[xArg], data0D[yArg], 'o-')
    ax.set_xlabel(xLabel, fontsize=16)
    ax.set_ylabel(yLabel, fontsize=16)
    plt.show()

    # Save
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