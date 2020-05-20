import ipywidgets as widgets
import PyNexus as PN
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
import base64
from IPython.display import clear_output, Javascript, display, HTML
import time
import subprocess
import sys
import nbformat as nbf
import CustomFunctions as CF

__version__ = '0.3'

"""
-Here are defined all the functions relevant to the front end of JupyLabBook,
i.e. the widgets (buttons, interactive plots), the self-generation of cells, conversion and saving of files ...
-Arguments are passed from the Notebook through objects belonging to the classes Experiment and Scan only.
-Modify the function Choose_treatment() to add your custom functions defined in CustomFunctions.py
"""

class Scan:
    """
    Class Scan is used to pass arguments concerning the current scan only.
    """
    def __init__(self):
        pass


def Print_version():
    print("Versions of modules used:")
    print("CustomFunctions: %s"%CF.__version__)
    print("FrontendFunctions: %s"%__version__)
    print("PyNexus: %s"%PN.__version__)
    print("Check that you are using the last versions of the modules and read the manual on: \n%s"%"https://github.com/ArnaudHemmerle/JupyLabBook"+'\n')

def Check_files(expt):

    Print_version()
    
    if not os.path.exists(expt.working_dir):
        print(PN._RED+"The following folder does not exist and should be created:"+PN._RESET)
        print(expt.working_dir)
        print("Data analysis will be saved in this folder.")
        print("")
    if not os.path.exists(expt.recording_dir):
        print(PN._RED+"The following folder does not exist:"+PN._RESET)
        print(expt.recording_dir)
        print("This folder should contain the nexus files.")
        print("")
    if not os.path.exists(expt.notebook_name):
        print(PN._RED+"The following file does not exist:"+PN._RESET)
        print(expt.notebook_name)
        print("Check that you have assigned the correct notebook name to expt.notebook_name")
        print("")
    if not os.path.exists('latex_template.tplx'):
        print(PN._RED+"The following file does not exist:"+PN._RESET)
        print('latex_template.tplx')
        print("This file contains the template for generating PDF and should be placed in the same folder as the notebook.")
        print("") 


def Define_scan_identifiers(scan, expt):
    """
    Create a series of identifiers for the current scan.
    """
    # For example:
    # scan.nxs = 'SIRIUS_2017_12_11_08042.nxs'
    # scan.path = '/Users/arnaudhemmerle/recording/SIRIUS_2017_12_11_08042.nxs'
    # scan.id = 'SIRIUS_2017_12_11_08042'
    # scan.number = 8042
    
    scan.path = expt.recording_dir+scan.nxs
    scan.id = scan.nxs[:-4]
    split_name = scan.nxs.split('.')[0].split('_')
    scan.number = int(scan.nxs.split('.')[0].split('_')[-1])


def Find_command_in_logs(scan, expt):
    """
    Find the command corresponding to the scan in the logs.
    """
    command_found = False
    scan.command = 'No command found'
    

    for log_file in expt.list_logs_files:
        with open(expt.logs_dir+log_file) as f:
            for line in f:
                if "#" not in line: temp = line
                if scan.id in line:
                    scan.command = temp
                    command_found = True
            f.close()

        if command_found:
            break
    
    
def Choose_action(expt):
    """
    Choose the next action to do. Help the choice by printing the command corresponding to the scan displayed in the list.
    Take an object from the class Experiment.
    Return an object from the class Scan with info related to the scan which will be treated. 
    """
    
    # Create an object for the current scan
    scan = Scan()

    expt.list_nxs_files = [file for file in sorted(os.listdir(expt.recording_dir)) if 'nxs' in file][::-1]
    if expt.list_nxs_files == []:
        print(PN._RED+'There is no nexus file in the recording folder.'+PN._RESET)
        print(PN._RED+'Recording folder: %s'%expt.recording_dir+PN._RESET)
        expt.list_nxs_files = ['SIRIUS_NoFileFound_00_00_00.nxs']
    
    expt.list_logs_files = [file for file in sorted(os.listdir(expt.logs_dir)) if 'log' in file][::-1]
    
    def selection_scan(nxs_file = expt.list_nxs_files):
        """
        Called by the widget to select the scan to be treated.
        """
    
        # Generate several identifiers for the scan
        scan.nxs = nxs_file
        Define_scan_identifiers(scan, expt)
        
        # Find the scan in the log files and extract/display the corresponding command
        Find_command_in_logs(scan, expt)

        if scan.command == 'No command found':
            print('Log file not found.')
        else:
            print("Corresponding command: ", scan.command)


    def on_button_treat_clicked(b):
        """
        Generate and execute cells corresponding to the chosen scan.
        """
        clear_output(wait=False)
        
        Create_cell(code='FF.Choose_treatment(scan, expt)',
                    position ='below', celltype='code', is_print=False)

        Create_cell(code='### '+scan.id+': '+scan.command,
                    position ='below', celltype='markdown', is_print=True)
        
       
    def on_button_list_clicked(b):
        """
        Selection of multiple scans.
        """
        clear_output(wait=False)
        Choose_list_scans(expt)
        
    def on_button_refresh_clicked(b):
        """
        Re-execute the cell to have the list of files updated.
        """
        # Create a unique id based on epoch time
        display_id = int(time.time()*1e9)

        display(Javascript("""IPython.notebook.execute_selected_cells();"""),display_id=display_id)

        # Necessary hack to avoid self-execution of cells at notebook re-opening
        # See http://tiny.cc/fnf3nz
        display(Javascript(""" """), display_id=display_id, update=True)
        
        
    def on_button_calibthetaz_clicked(b):
        """
        Create a cell for the calibration of thetaz.
        """
        clear_output(wait=False)
        
        Create_cell(code=
                    'calib_thetaz_data=np.array(['+'\n'+
                    '# gamma  channel'+'\n'+     
                    '[0,      970],'+'\n'+
                    '[-1,     899],'+'\n'+
                    '[-2,     827],'+'\n'+
                    '[-3,     755],'+'\n'+
                    '[-4,     683],'+'\n'+
                    '[-5,     611],'+'\n'+
                    '[-6,     539],'+'\n'+
                    '])'+'\n'+
                    'expt.thetazfactor=CF.Calib_thetaz(calib_thetaz_data)',
                    position='below', celltype='code', is_print = True)
        
        Create_cell(code='## Calibration thetaz', position='below', celltype='markdown', is_print = True)
        

        Create_cell(code='scan = FF.Choose_action(expt)', position ='at_bottom', celltype='code', is_print=False)        
        
        
    def on_button_form_clicked(b):
        """
        Create a form.
        """
        clear_output(wait=False)
        Create_cell(code='FF.Create_form()',position ='below', celltype='code', is_print=False)
        Create_cell(code='scan = FF.Choose_action(expt)', position ='at_bottom', celltype='code', is_print=False)        
     
            
    def on_button_export_clicked(b):
        """
        Export the notebook to PDF.
        """
        
        print('Export in progress...')
        
        export_done = Export_nb_to_pdf(expt.notebook_name)
        
        if export_done:
            print('Notebook exported to %s.pdf'%expt.notebook_name.split('.')[0])
        else:
            print("There was something wrong with the export to pdf. Please try again.")
            

            
    def on_button_HR_export_clicked(b):
        """
        Export the notebook to PDF in high-resolution.
        """

        print('Export in progress (it may take a long time)...')
        
        # Do first an export to low resolution PDF.
        export_LR_done = Export_nb_to_pdf(expt.notebook_name)
                   
        # Dupplicate the notebook with an extra cell at the beginning setting the 
        # images to pdf (High Resolution)
        nb = nbf.read(expt.notebook_name, as_version=4)
        code='is_HR = True'
        new_cell = nbf.v4.new_code_cell(code)
        nb['cells'].insert(0, new_cell)              
        nbf.write(nb, expt.notebook_name[:-6]+'_HR.ipynb')
        
        # Execute all the cells in the HR notebook
        command = 'jupyter nbconvert --to notebook --inplace --execute --allow-errors --ExecutePreprocessor.timeout=-1 '
        command+= expt.notebook_name[:-6]+'_HR.ipynb'
        rc = subprocess.call(command,shell=True)
        if rc>1:
            export_ipynb_done = False
        else:
            export_ipynb_done = True
        
        # Export the HR ipynb to PDF
        export_HR_done = Export_nb_to_pdf(expt.notebook_name[:-6]+'_HR.ipynb')
        
        # Remove the HR ipynb
        os.remove(expt.notebook_name[:-6]+'_HR.ipynb')
        
        if export_LR_done:
            print('Notebook exported in low resolution to %s.pdf'%expt.notebook_name.split('.')[0])
            
        if (export_LR_done and export_ipynb_done and export_HR_done) :
            print('Notebook exported in high resolution to %s_HR.pdf'%expt.notebook_name.split('.')[0])
        else:
            print("There was something wrong with the export to pdf. Please try again.")
            
    
    # Display the widgets
   
    # Click to treat a single scan
    button_treat = widgets.Button(description="Treat scan")
    button_treat.on_click(on_button_treat_clicked)
    
    # Click to treat a list of scans
    button_list = widgets.Button(description="Treat a list of scans")
    button_list.on_click(on_button_list_clicked)
    
    # Click to refresh the list of files
    button_refresh = widgets.Button(description="Refresh")
    button_refresh.on_click(on_button_refresh_clicked)
    
    # Click to do the calibthetaz
    button_calibthetaz = widgets.Button(description="Calib. theta z")
    button_calibthetaz.on_click(on_button_calibthetaz_clicked)
        
    # Click to fill the form
    button_form = widgets.Button(description="Fill form")
    button_form.on_click(on_button_form_clicked)
    
    # Click to export to pdf
    button_export = widgets.Button(description="Export to PDF")
    button_export.on_click(on_button_export_clicked)
    
    # Click to export to pdf in high resolution
    button_HR_export = widgets.Button(description="Export to PDF in HR")
    button_HR_export.on_click(on_button_HR_export_clicked)

    buttons0 = widgets.HBox([button_form, button_calibthetaz])
    display(buttons0)
    
    buttons1 = widgets.HBox([button_export, button_HR_export])
    display(buttons1)
    
    # Select a single scan and print the corresponding command
    w_print_scan = widgets.interact(selection_scan)
    w_print_scan.widget.children[0].description = 'Next scan:'
    w_print_scan.widget.children[0].layout = {'width': '400px'}
    
    buttons2 = widgets.HBox([button_treat, button_list, button_refresh])
    display(buttons2)
    
    return scan


def Choose_list_scans(expt):
    """
    Allow to choose multiple scans obtained in the same conditions.
    Take an object from the class Experiment.
    """

       
    scan_list = widgets.SelectMultiple(
    options=expt.list_nxs_files,
    rows=10,
    layout= {'width': '800px'},
    style={'description_width': 'initial'},
    description='Select scans (hold shift or ctrl):'
    )
    
    display(scan_list)
    
    
    def on_button_GIXD_clicked(b):     
   
        for nxs_file in scan_list.value:
            
            scan = Scan()
            
            # Generate several identifiers for the scan
            scan.nxs = nxs_file
            Define_scan_identifiers(scan, expt)

            # Find the scan in the log files and extract/display the corresponding command
            Find_command_in_logs(scan, expt)
            
            Create_cell(code='CF.Extract_GIXD(nxs_filename=\''+scan.nxs+'\','+
                        'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                        'logx=False, logy=False, logz=False, '+
                        'channel0=expt.channel0, thetazfactor=expt.thetazfactor, '+
                        'wavelength=expt.wavelength, thetac=expt.thetac, thetai=expt.thetai, '+
                        'binsize=expt.binsize, computeqz=True, nblevels=expt.nblevels, moytocreate=expt.moytocreate, '+
                        'show_data_stamps=False, show_saved_data=False, plot_true_GIXD=False)',
                        position='below', celltype='code', is_print = True)
            
            Create_cell(code='### '+scan.id+': '+scan.command,
                    position ='below', celltype='markdown', is_print=True)
            
            
    def on_button_true_GIXD_clicked(b):

        for nxs_file in scan_list.value:

            scan = Scan()

            # Generate several identifiers for the scan
            scan.nxs = nxs_file
            Define_scan_identifiers(scan, expt)

            # Find the scan in the log files and extract/display the corresponding command
            Find_command_in_logs(scan, expt)
        
            Create_cell(code='CF.Extract_GIXD(nxs_filename=\''+scan.nxs+'\','+
                        'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                        'logx=False, logy=False, logz=False, '+
                        'channel0=expt.channel0, thetazfactor=expt.thetazfactor, '+
                        'wavelength=expt.wavelength, thetac=expt.thetac, thetai=expt.thetai, '+
                        'binsize=expt.binsize, computeqz=True, nblevels=expt.nblevels, moytocreate=expt.moytocreate, '+
                        'show_data_stamps=False, show_saved_data=False, plot_true_GIXD=True)',
                        position='below', celltype='code', is_print = True)
            
            Create_cell(code='### '+scan.id+': '+scan.command,
                    position ='below', celltype='markdown', is_print=True)
            
            
    def on_button_plot_isotherm_clicked(b):

        for nxs_file in scan_list.value:
            
            scan = Scan()

            # Generate several identifiers for the scan
            scan.nxs = nxs_file
            Define_scan_identifiers(scan, expt)

            # Find the scan in the log files and extract/display the corresponding command
            Find_command_in_logs(scan, expt)
            Create_cell(code='CF.Plot_isotherm(nxs_filename=\''+scan.nxs+'\', recording_dir=expt.recording_dir,'+
                        'show_data_stamps=False)',
                        position='below', celltype='code', is_print = True)
            
            Create_cell(code='### '+scan.id+': '+scan.command,
                    position ='below', celltype='markdown', is_print=True)
        
        
    def on_button_next_clicked(b):
        clear_output(wait=False)
        Create_cell(code='scan = FF.Choose_action(expt)',
                    position ='at_bottom', celltype='code', is_print=False)
        
    
    button_GIXD = widgets.Button(description="Plot GIXD")
    button_GIXD.on_click(on_button_GIXD_clicked)
    
    button_true_GIXD = widgets.Button(description="Plot True GIXD")
    button_true_GIXD.on_click(on_button_true_GIXD_clicked)
    
    button_plot_isotherm = widgets.Button(description="Plot isotherm")
    button_plot_isotherm.on_click(on_button_plot_isotherm_clicked)
    
    button_next = widgets.Button(description="Next action")
    button_next.on_click(on_button_next_clicked)

    buttons0 = widgets.HBox([button_GIXD, button_true_GIXD, button_plot_isotherm, button_next])
    display(buttons0)
   

def Export_nb_to_pdf(nb_name):
    """
    Save the current state of the notebook (including the widgets).
    Export the notebook to pdf using a command line through the OS.
    Return a boolean which is True if the export suceeded without error/warning.
    """
    display(HTML('<script>Jupyter.menubar.actions._actions["widgets:save-with-widgets"].handler()</script>') )
    t0 = time.time()
    rc = 1
    while rc>0:
        if (time.time()-t0) > 100:
            # Timeout before PDF export is considered as failed
            export_done = False
            break
        else:
            time.sleep(3)
            command = 'jupyter nbconvert '
            command+= nb_name
            command+= ' --to pdf '
            command+= ' --TagRemovePreprocessor.remove_cell_tags=\"[\'notPrint\']\" ' # Remove the widgets from the PDF
            command+= ' --no-input ' # Remove the code cells
            command+= '--template latex_template.tplx' # Custom template
            rc = subprocess.call(command,shell=True)
            if rc==0: export_done = True
                
    return export_done

def Delete_current_cell():
    """Delete the cell which called this function in the IPython Notebook."""
    
    display(Javascript(
        """
        var index = IPython.notebook.get_selected_cells_indices();
        IPython.notebook.delete_cell(index);
        """
    ))

def Create_cell(code='', position='below', celltype='markdown', is_print = False):
    """Create a cell in the IPython Notebook.
    code: unicode, Code to fill the new cell with.
    celltype: unicode, Type of cells "code" or "markdown".
    position: unicode, Where to put the cell "below" or "at_bottom"
    is_print: boolean, To decide if the cell is printed in the pdf report
    The code requires direct use of Javascript and is thus not usable in Jupyter Lab.
    """

    encoded_code = (base64.b64encode(code.encode())).decode()

    # Create a unique id based on epoch time
    display_id = int(time.time()*1e9)

    if is_print:
        display(Javascript("""
        var cell = IPython.notebook.insert_cell_{0}("{1}");
        cell.set_text(atob("{2}"));
        cell.execute();
        """.format(position, celltype, encoded_code)),display_id=display_id)

    else:
        display(Javascript("""
        var cell = IPython.notebook.insert_cell_{0}("{1}");
        cell.set_text(atob("{2}"));
        cell.metadata.tags = ['notPrint']
        cell.execute();
        """.format(position, celltype, encoded_code)),display_id=display_id)

    # Necessary hack to avoid self-execution of cells at notebook re-opening
    # See http://tiny.cc/fnf3nz
    display(Javascript(""" """), display_id=display_id, update=True)


def Set_interactive_1D(scan):
    """
    Take an object from the class Scan.
    1) Extract the sensors from the nxs file.
    2) Set an interactive 1D plot.
    """

    nexus = PN.PyNexusFile(scan.path, fast=True)
    stamps0D, data0D = nexus.extractData('0D')
    sensor_list = [stamps0D[i][0] if stamps0D[i][1]== None else stamps0D[i][1] for i in range(len(stamps0D))]


    def plot_interactive_1D(xLabel, yLabel):
        xArg = sensor_list.index(xLabel)
        yArg = sensor_list.index(yLabel)

        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(data0D[xArg], data0D[yArg], 'o-')
        ax.set_xlabel(xLabel, fontsize=16)
        ax.set_ylabel(yLabel, fontsize=16)
        plt.show()

        scan.xLabel = xLabel
        scan.yLabel = yLabel

    widgets.interact(plot_interactive_1D, xLabel=sensor_list, yLabel=sensor_list)
    


def Choose_treatment(scan, expt):
    """
    Take objects from the class Scan and Experiment.
    1) Allow the user to plot sensors from the scan using the interactive widget.
    2) Display the buttons for choosing the next action.
    """
    
    
    # Define the function called when clicking the button
    # DEFINE HERE A FUNCTION TO CREATE A CELL CALLING YOUR CUSTOM FUNCTION
    def on_button_1D_clicked(b):
        Create_cell(code='CF.Plot_1D(nxs_filename=\''+scan.nxs+'\', recording_dir=expt.recording_dir,'+
                    'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)

    def on_button_GIXD_clicked(b):
        Create_cell(code='CF.Extract_GIXD(nxs_filename=\''+scan.nxs+'\','+
                    'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                    'logx=False, logy=False, logz=False, '+
                    'channel0=expt.channel0, thetazfactor=expt.thetazfactor, '+
                    'wavelength=expt.wavelength, thetac=expt.thetac, thetai=expt.thetai, '+
                    'binsize=expt.binsize, computeqz=True, nblevels=expt.nblevels, moytocreate=expt.moytocreate, '+
                    'show_data_stamps=False, show_saved_data=False, plot_true_GIXD=False)',
                    position='below', celltype='code', is_print = True)

    def on_button_true_GIXD_clicked(b):
        Create_cell(code='CF.Extract_GIXD(nxs_filename=\''+scan.nxs+'\','+
                    'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                    'logx=False, logy=False, logz=False, '+
                    'channel0=expt.channel0, thetazfactor=expt.thetazfactor, '+
                    'wavelength=expt.wavelength, thetac=expt.thetac, thetai=expt.thetai, '+
                    'binsize=expt.binsize, computeqz=True, nblevels=expt.nblevels, moytocreate=expt.moytocreate, '+
                    'show_data_stamps=False, show_saved_data=False, plot_true_GIXD=True)',
                    position='below', celltype='code', is_print = True)

    def on_button_next_clicked(b):
        #clear_output(wait=False)
        
        Delete_current_cell()
        
        Create_cell(code='scan = FF.Choose_action(expt)',
                    position ='at_bottom', celltype='code', is_print=False)

    def on_button_fit_erf_clicked(b):
        Create_cell(code='CF.GaussianRepartition_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)  
        
    def on_button_fit_gau_clicked(b):
        Create_cell(code='CF.Gaussian_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)   
        
    def on_button_vineyard_clicked(b):
        Create_cell(code='expt.channel0 = CF.Extract_channel_Qc(nxs_filename=\''+scan.nxs+'\','+
                    'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                    'logx=False, logy=False, logz=False)',
                    position='below', celltype='code', is_print = True)
        
    def on_button_plot_isotherm_clicked(b):
        Create_cell(code='CF.Plot_isotherm(nxs_filename=\''+scan.nxs+'\', recording_dir=expt.recording_dir,'+
                    'show_data_stamps=False)',
                    position='below', celltype='code', is_print = True)

    # Display the widgets    
    # ADD HERE A CALL TO YOUR BUTTON
    button_fit_gau = widgets.Button(description="Fit with gaussian")
    button_fit_gau.on_click(on_button_fit_gau_clicked)
    
    button_fit_erf = widgets.Button(description="Fit with erf")
    button_fit_erf.on_click(on_button_fit_erf_clicked)
    
    # Click to extract the Vineyard
    button_vineyard = widgets.Button(description="Extract Vineyard")
    button_vineyard.on_click(on_button_vineyard_clicked)
    
    button_1D = widgets.Button(description="Add plot to report")
    button_1D.on_click(on_button_1D_clicked)
    
    button_GIXD = widgets.Button(description="Plot GIXD")
    button_GIXD.on_click(on_button_GIXD_clicked)
    
    button_true_GIXD = widgets.Button(description="Plot True GIXD")
    button_true_GIXD.on_click(on_button_true_GIXD_clicked)
    
    button_plot_isotherm = widgets.Button(description="Plot isotherm")
    button_plot_isotherm.on_click(on_button_plot_isotherm_clicked)
    
    button_next = widgets.Button(description="Next action")
    button_next.on_click(on_button_next_clicked)
   
    buttons0 = widgets.HBox([button_fit_gau, button_fit_erf, button_1D])
    display(buttons0)
    
        
    # Set up an interactive 1D plot
    Set_interactive_1D(scan)
    
    
    buttons1 = widgets.HBox([button_GIXD, button_true_GIXD, button_plot_isotherm])
    display(buttons1)

    buttons2 = widgets.HBox([button_vineyard, button_next])
    display(buttons2)

def Create_form():

    short_layout = widgets.Layout(width='200px', height='40px')
    medium_layout = widgets.Layout(width='300px', height='40px')
    long_layout = widgets.Layout(width='800px', height='40px')
    style = {'description_width': 'initial'}

    title = widgets.Text(
            placeholder='Title',
            description='Experiment title:',
            layout=long_layout,
            style=style)

    number = widgets.Text(
            placeholder='Number',
            description='Experiment number:',
            layout=medium_layout,
            style=style)

    ttype = widgets.Text(
            placeholder='Proposal, Commissioning, ...',
            description='Type:',
            layout=medium_layout,
            style=style)

    safety = widgets.Dropdown(
            options=['Green', 'Yellow', 'Red'],
            value='Yellow',
            description='Safety:',
            layout=widgets.Layout(width='150px'),
            style=style)

    date  = widgets.Text(
            placeholder='DD/MM/YYYY - DD/MM/YYYY',
            description='Date:',
            layout=medium_layout,
            style=style)

    proposer = widgets.Text(
               placeholder='Name',
               description='Main proposer:',
               layout=short_layout,
               style=style)

    contact = widgets.Text(
               placeholder='Name',
               description='Local contact:',
               layout=medium_layout,
               style=style)

    users = widgets.Text(
            placeholder='Names',
            description='Users (on site):',
            layout=long_layout,
            style=style)

    recording = widgets.Text(
               placeholder='Full path to the recording directory',
               description='Recording directory:',
               layout=long_layout,
               style=style)


    current = widgets.Text(
               placeholder='450 mA, 500 mA, ...',
               description='Current:',
               layout=medium_layout,
               style=style)

    mode     = widgets.Text(
               placeholder='Hybrid, Top-up, ...',
               description='Mode:',
               layout=medium_layout,
               style=style)

    dcm   = widgets.Dropdown(
            options=['Si111', 'InSb', 'Not Used'],
            value='Si111',
            description='DCM:',
            layout=widgets.Layout(width='150px'),
            style=style)

    mgm   = widgets.Dropdown(
            options=['In use', 'Not used'],
            value='Not used',
            description='MGM:',
            layout=widgets.Layout(width='150px'),
            style=style)

    energy  = widgets.Text(
               placeholder='Value(s)',
               description='Energy (keV):',
               layout=medium_layout,
               style=style)

    energy_type = widgets.Dropdown(
                  options=['Fixed', 'Variable'],
                  value='Fixed',
                  description='Fixed/Variable energy:',
                  layout=widgets.Layout(width='300px'),
                  style=style)

    wavelength = widgets.Text(
                 placeholder='Value(s)',
                 description='Wavelength (nm):',
                 layout=medium_layout,
                 style=style)

    harmonic  =  widgets.Text(
                 placeholder='Value(s)',
                 description='Harmonic:',
                 layout=short_layout,
                 style=style)

    polarisation = widgets.Text(
                 placeholder='Value(s)',
                 value='LH',
                 description='Polarisation:',
                 layout=short_layout,
                 style=style)

    phase = widgets.Text(
                 placeholder='Value(s)',
                 value='0',
                 description='Phase (deg):',
                 layout=short_layout,
                 style=style)

    m1  = widgets.Dropdown(
          options=['M1-A Pt Track', 'M1-A B4C Track', 'M1-B (B4C)', 'No M1'],
          value='M1-A Pt Track',
          description='M1:',
          layout=widgets.Layout(width='200px'),
          style=style)

    m2  = widgets.Dropdown(
          options=['M2 Pt Track', 'M2 B4C Track', 'No M2'],
          value='M2 Pt Track',
          description='M2:',
          layout=widgets.Layout(width='200px'),
          style=style)

    m3  = widgets.Dropdown(
          options=['M3 Pt Track', 'M3 B4C Track', 'No M3'],
          value='No M3',
          description='M3:',
          layout=widgets.Layout(width='200px'),
          style=style)

    m4  = widgets.Dropdown(
          options=['M4 Pt Track', 'M4 Si Track', 'No M4'],
          value='M4 Pt Track',
          description='M4:',
          layout=widgets.Layout(width='200px'),
          style=style)

    horizontal_foc = widgets.Checkbox(
                     value=True,
                     description='Horizontal focalisation',
                     layout=long_layout,
                     style=style)

    vertical_foc = widgets.Checkbox(
                   value=True,
                   description='Vertical focalisation',
                   layout=long_layout,
                   style=style)


    horizontal_size = widgets.Text(
                      placeholder='Value(s)',
                      description='Horizontal beamsize (mm):',
                      layout=medium_layout,
                      style=style)

    vertical_size   = widgets.Text(
                      placeholder='Value(s)',
                      description='Vertical beamsize (mm):',
                      layout=medium_layout,
                      style=style)

    mon1  = widgets.Text(
            placeholder='Empty or type of mon',
            description='mon1:',
            layout=short_layout,
            style=style)

    mon2  = widgets.Text(
            placeholder='Empty or type of mon',
            description='mon2:',
            layout=short_layout,
            style=style)

    mon3  = widgets.Text(
            placeholder='Empty or type of mon',
            description='mon3:',
            layout=short_layout,
            style=style)

    mon4  = widgets.Text(
            placeholder='Empty or type of mon',
            description='mon4:',
            layout=short_layout,
            style=style)

    detectors  = widgets.Textarea(
                 placeholder='Type of detectors',
                 description='Detectors:',
                 layout=long_layout,
                 style=style)

    remarks     = widgets.Textarea(
                 placeholder='Remarks',
                 description='Remarks:',
                 layout=long_layout,
                 style=style)

    display(title)
    display(widgets.HBox([date, number]))
    display(widgets.HBox([ttype, safety]))
    print('-'*100)
    display(widgets.HBox([proposer, contact]))
    display(users)
    print('-'*100)
    display(recording)
    print('\033[1m'+'Machine:')
    display(widgets.HBox([current, mode]))
    print('\033[1m'+'Optics:')
    display(widgets.HBox([dcm, mgm]))
    display(widgets.HBox([m1, m2, m3, m4]))
    print('\033[1m'+'Beam:')
    display(widgets.HBox([energy_type, energy, wavelength]))
    display(widgets.HBox([harmonic, polarisation, phase]))
    display(widgets.HBox([horizontal_foc, vertical_foc]))
    display(widgets.HBox([horizontal_size, vertical_size]))
    print('\033[1m'+'Monitors and XBPM:')
    display(widgets.HBox([mon1, mon2, mon3, mon4]))
    display(detectors)
    print('\033[1m'+'Remarks:')
    display(remarks)
    
    def on_button_clicked(b):
        
        txt = []
        txt.append('$\LARGE \\textbf{SIRIUS Beamline}:\\textbf{Experiment %s}$'%number.value)

        #Add spaces
        ptitle = ('\ '.join(title.value.split(' ')))
        txt.append('$\Large \\color{red}{\\bf %s}$'%ptitle)

        txt.append('* %s %s'%(ttype.description,ttype.value)+'\n'
                    +'* %s %s'%(safety.description,safety.value)+'\n'
                    +'* %s %s'%(date.description,date.value))

        txt.append('* %s %s'%(proposer.description,proposer.value)+'\n'
                    +'* %s %s'%(contact.description,contact.value)+'\n'
                    +'* %s %s'%(users.description,users.value)+'\n'
                    +'* %s %s'%(recording.description,recording.value))

        txt.append('* Machine:'+'\n'
                    +'\t * %s %s'%(current.description,current.value)+'\n'
                    +'\t * %s %s'%(mode.description,mode.value))

        txt.append('* Optics:'+'\n'
                    +'\t * %s %s'%(dcm.description,dcm.value)+'\n'
                    +'\t * %s %s'%(mgm.description,mgm.value)+'\n'
                    +'\t * %s %s'%(m1.description,m1.value)+'\n'
                    +'\t * %s %s'%(m2.description,m2.value)+'\n'
                    +'\t * %s %s'%(m3.description,m3.value)+'\n'
                    +'\t * %s %s'%(m4.description,m4.value))

        txt.append('* Beam:'+'\n'
                    +'\t * %s %s'%(energy_type.description,energy_type.value)+'\n'
                    +'\t * %s %s'%(energy.description,energy.value)+'\n'
                    +'\t * %s %s'%(wavelength.description,wavelength.value)+'\n'
                    +'\t * %s %s'%(harmonic.description,harmonic.value)+'\n'
                    +'\t * %s %s'%(polarisation.description,polarisation.value)+'\n'
                    +'\t * %s %s'%(phase.description,phase.value)+'\n'
                    +'\t * %s %s'%(horizontal_foc.description,horizontal_foc.value)+'\n'
                    +'\t * %s %s'%(vertical_foc.description,vertical_foc.value)+'\n'
                    +'\t * %s %s'%(horizontal_size.description,vertical_size.value))

        txt.append('* Monitors and XBPM:'+'\n'
                    +'\t * %s %s'%(mon1.description,mon1.value)+'\n'
                    +'\t * %s %s'%(mon2.description,mon2.value)+'\n'
                    +'\t * %s %s'%(mon3.description,mon3.value)+'\n'
                    +'\t * %s %s'%(mon4.description,mon4.value)+'\n'
                    +'\t * %s %s'%(detectors.description,detectors.value))

        txt.append('* %s %s'%(remarks.description,remarks.value))

        txt.reverse()
        
        for elem in txt:
            Create_cell(code=elem, position ='below', celltype='markdown', is_print=True)
        
        Create_cell(code='# Form', position ='below', celltype='markdown', is_print=True)
        
    button = widgets.Button(description="Print form")
    out = widgets.Output()
    display(button,out)
    button.on_click(on_button_clicked)
