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

__version__ = '0.7'

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
    if not os.path.exists(expt.logs_dir):
        print(PN._RED+"The following folder does not exist:"+PN._RESET)
        print(expt.logs_dir)
        print("This folder should contain the log files.")
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
                    # Remove the line jump
                    scan.command = temp.replace('\n','')
                    command_found = True
            f.close()

        if command_found:
            break
    
    
def Choose_action(expt):
    """
    Choose the next action to do. Help the choice by printing the command corresponding to the scans displayed in the list.
    Take an object from the class Experiment.
    Return a list of objects from the class Scan with info related to the scans which will be treated. 
    """
    
    # Define the list of nxs files in the recording directory
    expt.list_nxs_files = [file for file in sorted(os.listdir(expt.recording_dir)) if 'nxs' in file][::-1]
    if expt.list_nxs_files == []:
        print(PN._RED+'There is no nexus file in the recording folder.'+PN._RESET)
        print(PN._RED+'Recording folder: %s'%expt.recording_dir+PN._RESET)
        expt.list_nxs_files = ['SIRIUS_NoFileFound_00_00_00.nxs']
    
    # Define the list of log files in the log directory
    expt.list_logs_files = [file for file in sorted(os.listdir(expt.logs_dir)) if 'log' in file][::-1]
    
    def selection_scans(nxs_files):
        """
        Called by the widget to select the scans to be treated.
        """
        
        expt.scans = []
        
        for nxs_file in nxs_files:
            
            # Create an object scan
            scan = Scan()
        
            # Generate several identifiers for the scan
            scan.nxs = nxs_file
            Define_scan_identifiers(scan, expt)
        
            # Find the scan in the log files and extract/display the corresponding command
            Find_command_in_logs(scan, expt)

            if scan.command == 'No command found':
                print('%s: log file not found.'%scan.nxs)
            else:
                print('%s: %s'%(scan.nxs,scan.command))
                
            expt.scans.append(scan)
           

    def on_button_treat_clicked(b):
        """
        Generate and execute cells corresponding to the chosen scans.
        """
        clear_output(wait=False)
        
        Create_cell(code='FF.Choose_treatment(expt)',
                    position ='below', celltype='code', is_print=False)

        if len(expt.scans)==1:
            Create_cell(code='### '+expt.scans[0].id+': '+expt.scans[0].command,
                        position ='below', celltype='markdown', is_print=True)

        
    def on_button_refresh_clicked(b):
        """
        Re-execute the cell to update it.
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
        

        Create_cell(code='FF.Choose_action(expt)', position ='at_bottom', celltype='code', is_print=False)        
        
        
    def on_button_form_clicked(b):
        """
        Create a form.
        """
        clear_output(wait=False)
        Create_cell(code='FF.Create_form()',position ='below', celltype='code', is_print=False)
        Create_cell(code='FF.Choose_action(expt)', position ='at_bottom', celltype='code', is_print=False)        
     
            
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

                
    def on_button_markdown_clicked(b):
        """
        Insert a markdown cell above the current cell.
        """ 
        Create_cell(code='**Insert your comment (double-click here and replace current text).**',
                    position ='above', celltype='markdown', is_print=True)
    
                    
    def on_button_script_clicked(b):
        """
        Insert a script as a markdown cell.
        """ 

        # Check if there is already a path for scripts directories
        # If not, use the recording directory
        try:
            path_to_dir_default = expt.path_to_dir
        except:
            path_to_dir_default = expt.recording_dir

        path_to_dir = widgets.Text(
                value=path_to_dir_default,
                description='Scripts directory:',
                layout=widgets.Layout(width='900px', height='40px'),
                style={'description_width': 'initial'})
        display(path_to_dir)

        script_name = widgets.Text(
                placeholder='script_name.ipy',
                description='Script name:',
                layout=widgets.Layout(width='400px', height='40px'),
                style={'description_width': 'initial'})
        display(script_name)

        def on_button_display_clicked(b):
            """
            Display a script in a markdown cell.
            """ 
            
            on_button_refresh_clicked(b)
            
            # Pass the current value of the directory to the default one
            expt.path_to_dir = path_to_dir.value

            text_file = open(path_to_dir.value+script_name.value)
            file_content = text_file.read()
            text_file.close()
            code = '```python'+file_content+'```'

            Create_cell(code='### '+ script_name.value, position ='above', celltype='markdown', is_print=True)
            Create_cell(code=code, position ='above', celltype='markdown', is_print=True)
            

        button_display = widgets.Button(description="Add script to report")
        button_display.on_click(on_button_display_clicked)

        display(button_display)
    
    # Display the widgets
   
    # Click to treat a single scan
    button_treat = widgets.Button(description="Treat scan(s)")
    button_treat.on_click(on_button_treat_clicked)
    
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
    
    # Click to insert a markdown cell
    button_markdown = widgets.Button(description="Insert comment")
    button_markdown.on_click(on_button_markdown_clicked)
    
    # Click to insert a script
    button_script = widgets.Button(description="Insert script")
    button_script.on_click(on_button_script_clicked)

    buttons0 = widgets.HBox([button_treat, button_refresh])
    display(buttons0)
      
    # Widget for selection of multiple scans
    w_print_scans = widgets.interact(selection_scans, nxs_files = widgets.SelectMultiple(options=expt.list_nxs_files))
    w_print_scans.widget.children[0].description = 'Next scan(s):'
    w_print_scans.widget.children[0].layout = {'width': '400px'}
    
    buttons1 = widgets.HBox([button_form, button_calibthetaz])
    display(buttons1)

    buttons2 = widgets.HBox([button_markdown, button_script, button_export])
    display(buttons2)
    

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
    nexus.close()
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
    


def Choose_treatment(expt):
    """
    Take object from the class Experiment.
    1) Allow the user to plot sensors from the scan using the interactive widget (if only one scan was selected).
    2) Display the buttons for choosing the next action.
    """
    
    # Styling options for widgets
    style = {'description_width': 'initial'}
    tiny_layout = widgets.Layout(width='150px', height='40px')
    short_layout = widgets.Layout(width='200px', height='40px')
    
    # Define the function called when clicking the button
    # DEFINE HERE A FUNCTION TO CREATE A CELL CALLING YOUR CUSTOM FUNCTION

    def on_button_GIXD_clicked(b):
        for scan in expt.scans:
            
            Create_cell(code='CF.Extract_GIXD(nxs_filename=\''+scan.nxs+'\','+
                        'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                        'logx=False, logy=False, logz=False, '+
                        'channel0=expt.channel0, thetazfactor=expt.thetazfactor, '+
                        'wavelength=expt.wavelength, thetac=expt.thetac, thetai=expt.thetai, '+
                        'binsize=expt.binsize, computeqz=True, nblevels=expt.nblevels, moytocreate=expt.moytocreate, '+
                        'show_data_stamps='+str(w_stamps.value)+', verbose='+str(w_verbose.value)+', plot_true_GIXD=False)',
                        position='below', celltype='code', is_print = True)
            
            if len(expt.scans)>1:
                Create_cell(code='### '+scan.id+': '+scan.command,
                        position ='below', celltype='markdown', is_print=True)

    def on_button_true_GIXD_clicked(b):
        for scan in expt.scans:
       
            Create_cell(code='CF.Extract_GIXD(nxs_filename=\''+scan.nxs+'\','+
                        'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                        'logx=False, logy=False, logz=False, '+
                        'channel0=expt.channel0, thetazfactor=expt.thetazfactor, '+
                        'wavelength=expt.wavelength, thetac=expt.thetac, thetai=expt.thetai, '+
                        'binsize=expt.binsize, computeqz=True, nblevels=expt.nblevels, moytocreate=expt.moytocreate, '+
                        'show_data_stamps='+str(w_stamps.value)+', verbose='+str(w_verbose.value)+', plot_true_GIXD=True)',
                        position='below', celltype='code', is_print = True)
        
            if len(expt.scans)>1:
                Create_cell(code='### '+scan.id+': '+scan.command,
                        position ='below', celltype='markdown', is_print=True)

    def on_button_pilatus_clicked(b):
        for scan in expt.scans:

            Create_cell(code='CF.Extract_pilatus_sum(nxs_filename=\''+scan.nxs+'\','+ 
                       'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                       'logz=True, show_data_stamps='+str(w_stamps.value)+', verbose='+str(w_verbose.value)+', cmap=expt.cmap)',
                        position='below', celltype='code', is_print = True)
            
            if len(expt.scans)>1:
                Create_cell(code='### '+scan.id+': '+scan.command,
                        position ='below', celltype='markdown', is_print=True)     
        
    def on_button_GIXS_angles_clicked(b):
        for scan in expt.scans:
            
            Create_cell(code='CF.Extract_GIXS(nxs_filename=\''+scan.nxs+'\','+ 
                       'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                       'logz=True, wavelength=expt.wavelength, thetai=expt.thetai, distance=expt.distance, '+
                       'pixel_PONI_x=expt.pixel_PONI_x, pixel_PONI_y=expt.pixel_PONI_y, '+
                       'show_data_stamps='+str(w_stamps.value)+', verbose='+str(w_verbose.value)+' ,cmap=expt.cmap, '+
                        'plot_twotheta_alphaf=True)',
                        position='below', celltype='code', is_print = True)

            if len(expt.scans)>1:
                Create_cell(code='### '+scan.id+': '+scan.command,
                        position ='below', celltype='markdown', is_print=True)
               
    def on_button_GIXS_qxy_qz_clicked(b):
        for scan in expt.scans:
            
            Create_cell(code='CF.Extract_GIXS(nxs_filename=\''+scan.nxs+'\','+ 
                       'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                       'logz=True, wavelength=expt.wavelength, thetai=expt.thetai, distance=expt.distance, '+
                       'pixel_PONI_x=expt.pixel_PONI_x, pixel_PONI_y=expt.pixel_PONI_y,  '+
                       'show_data_stamps='+str(w_stamps.value)+', verbose='+str(w_verbose.value)+' ,cmap=expt.cmap, '+
                        'plot_qxy_qz=True)',
                        position='below', celltype='code', is_print = True)
            
            if len(expt.scans)>1:
                Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)
        
    def on_button_GIXS_qxy_q_clicked(b):
        for scan in expt.scans:
            
            Create_cell(code='CF.Extract_GIXS(nxs_filename=\''+scan.nxs+'\','+ 
                   'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                   'logz=True, wavelength=expt.wavelength, thetai=expt.thetai, distance=expt.distance, '+
                   'pixel_PONI_x=expt.pixel_PONI_x, pixel_PONI_y=expt.pixel_PONI_y, '+
                   'show_data_stamps='+str(w_stamps.value)+', verbose='+str(w_verbose.value)+' ,cmap=expt.cmap, '+
                   'plot_qxy_q=True)',
                        position='below', celltype='code', is_print = True)
            
            if len(expt.scans)>1:
                Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)
                
                
    def on_button_GIXS_clicked(b):
        
        # Checkboxes for options       

        # logz
        try: value = expt.GIXS_logz
        except: value = True
        w_GIXS_logz = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='log z')
        
        # wavelength
        try: value = expt.wavelength
        except: value = 0.155
        w_wavelength = widgets.FloatText(value=value, style=style, layout=short_layout, description='wavelength (nm)')

        # thetai
        try: value = expt.thetai
        except: value = 0.002
        w_thetai = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='thetai (rad)')
        
        # distance 
        try: value = expt.distance
        except: value = 1000.
        w_distance = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='distance (mm)')

        # pixel_PONI_x 
        try: value = expt.pixel_PONI_x
        except: value = 0.
        w_pixel_PONI_x = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='PONIx (pix)')
        
        # pixel_PONI_y 
        try: value = expt.pixel_PONI_y
        except: value = 0.
        w_pixel_PONI_y = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='PONIy (pix)')
        
        # pixel_size 
        try: value = expt.pixel_size
        except: value = 0.172
        w_pixel_size = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='Pixel size (um)')
            
        # show_data_stamps
        try: value = expt.show_data_stamps
        except: value = False
        w_show_data_stamps = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print sensors')

        # verbose
        try: value = expt.verbose
        except: value = False
        w_verbose = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print scan info')
          
        # cmap
        try: value = expt.cmap
        except: value = 'viridis'
        w_cmap = widgets.Text(value=value, style=style, layout=tiny_layout, description='cmap')
        
        # plot_twotheta_alphaf/plot_qxy_qz/plot_qxy_q
        try: value = expt.GIXS_plot_type
        except: value = 'pixels only'
        w_GIXS_plot_type = widgets.Select(value=value, style=style, rows=4,
                                          options=['pixels only', 'angles', 'qxy/qz', 'qxy/q'], description='Plot type')
            
        display(widgets.HBox([w_show_data_stamps, w_verbose, w_GIXS_logz, w_cmap, w_pixel_size]))        
        display(widgets.HBox([w_wavelength, w_distance, w_thetai, w_pixel_PONI_x, w_pixel_PONI_y]))
        

        
        def on_button_plot_clicked(b):
            
            # Pass current values as default values
            expt.GIXS_logz = w_GIXS_logz.value
            expt.wavelength = w_wavelength.value
            expt.thetai = w_thetai.value
            expt.distance = w_distance.value
            expt.pixel_PONI_x = w_pixel_PONI_x.value
            expt.pixel_PONI_y = w_pixel_PONI_y.value
            expt.pixel_size = w_pixel_size.value
            expt.show_data_stamps = w_show_data_stamps.value
            expt.verbose = w_verbose.value
            expt.cmap = w_cmap.value
            expt.GIXS_plot_type = w_GIXS_plot_type.value
            
            # Pass plot type to params    
            expt.plot_twotheta_alphaf = True if w_GIXS_plot_type.value == 'angles' else False
            expt.plot_qxy_qz = True if w_GIXS_plot_type.value == 'qxy/qz' else False
            expt.plot_qxy_q = True if w_GIXS_plot_type.value == 'qxy/q' else False

            
            for scan in expt.scans:

                Create_cell(code='CF.Extract_GIXS(nxs_filename=\''+scan.nxs+'\','+ 
                       'working_dir=expt.working_dir, recording_dir=expt.recording_dir,'+
                        'logz='+str(expt.GIXS_logz)+','+
                        'wavelength='+str(expt.wavelength)+','+
                        'thetai='+str(expt.thetai)+','+
                        'distance='+str(expt.distance)+','+
                        'pixel_PONI_x='+str(expt.pixel_PONI_x)+','+
                        'pixel_PONI_y='+str(expt.pixel_PONI_y)+','+
                        'pixel_size='+str(expt.pixel_size)+','+  
                        'show_data_stamps='+str(expt.show_data_stamps)+','+
                        'verbose='+str(expt.verbose)+','+
                        'cmap=\''+str(expt.cmap)+'\','+
                        'plot_twotheta_alphaf='+str(expt.plot_twotheta_alphaf)+','+
                        'plot_qxy_qz='+str(expt.plot_qxy_qz)+','+
                        'plot_qxy_q='+str(expt.plot_qxy_q)+')',
                        position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    Create_cell(code='### '+scan.id+': '+scan.command,
                                position ='below', celltype='markdown', is_print=True)            
            
        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        display(widgets.HBox([w_GIXS_plot_type,button_plot]))

    def on_button_fluo_clicked(b):
        
        # Checkboxes for options       
       
        # show_data_stamps
        try: value = expt.show_data_stamps
        except: value = False
        w_stamps = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Print sensors')

        # verbose
        try: value = expt.verbose
        except: value = False
        w_verbose = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Print scan info')
  
        # plot_spectrogram
        try: value = expt.is_plot_spectrogram
        except: value = True
        w_spectrogram = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Plot spectrogram')
        
        # plot_first_last
        try: value = expt.is_plot_first_last
        except: value = True
        w_first_last = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Plot first&last spectrum')
        
        # plot_sum
        try: value = expt.is_plot_sum
        except: value = True
        w_sum = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Plot sum of spectrums')
        
        # logz
        try: value = expt.is_XRF_logz
        except: value = True
        w_XRF_logz = widgets.Checkbox(value=value, style=style, layout = short_layout, description='log z')
        
        # list_elems
        try: value = expt.elems_str
        except: value = '0, 1, 2, 3'
        w_elems = widgets.Text(value=value, description='Elements', style=style, layout = short_layout)

        # first_channel
        try: value = expt.first_channel
        except: value = 0
        w_first_channel = widgets.IntText(value=value, description='First channel', style=style, layout = short_layout)
        
        # last_channel
        try: value = expt.last_channel
        except: value = 2047
        w_last_channel = widgets.IntText(value=value, description='Last channel', style=style, layout = short_layout)

        
        display(widgets.HBox([w_stamps, w_verbose]))
        display(widgets.HBox([w_spectrogram, w_first_last, w_sum]))
        display(widgets.HBox([w_XRF_logz, w_elems, w_first_channel, w_last_channel]))
            
        def on_button_plot_clicked(b):

            # Pass current values as default values
            expt.is_plot_spectrogram = w_spectrogram.value
            expt.is_plot_first_last = w_first_last.value
            expt.is_plot_sum = w_sum.value
            expt.is_XRF_logz = w_XRF_logz.value
            expt.elems_str = w_elems.value
            expt.first_channel = w_first_channel.value
            expt.last_channel = w_last_channel.value
            expt.show_data_stamps = w_stamps.value
            expt.verbose = w_verbose.value
            
            # Convert elems_str into a list
            list_elems = [int(expt.elems_str.split(',')[i]) for i in range(len(expt.elems_str.split(',')))]
            
            for scan in expt.scans:

                Create_cell(code='CF.Extract_fluo('+
                           'nxs_filename=\''+scan.nxs+'\','+ 
                           'working_dir=expt.working_dir,recording_dir=expt.recording_dir,'+
                           'logz='+str(w_XRF_logz.value)+','+
                           'list_elems='+str(list_elems)+','+
                           'first_channel='+str(w_first_channel.value)+','+
                           'last_channel='+str(w_last_channel.value)+','+
                           'show_data_stamps='+str(w_stamps.value)+','+
                           'verbose='+str(w_verbose.value)+','+
                           'plot_spectrogram='+str(w_spectrogram.value)+','+
                           'plot_first_last='+str(w_first_last.value)+','+
                           'plot_sum='+str(w_sum.value)+
                           ')',
                           position='below', celltype='code', is_print = True)  

                if len(expt.scans)>1:
                    Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)  
                    
        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        
        display(button_plot)

    def on_button_isotherm_clicked(b):
        for scan in expt.scans:
            
            Create_cell(code='CF.Plot_isotherm(nxs_filename=\''+scan.nxs+'\', recording_dir=expt.recording_dir, '+
                       'show_data_stamps='+str(w_stamps.value)+', verbose='+str(w_verbose.value)+')',
                        position='below', celltype='code', is_print = True)
            
            if len(expt.scans)>1:
                Create_cell(code='### '+scan.id+': '+scan.command,
                        position ='below', celltype='markdown', is_print=True)

    # Actions relevant for single scan analysis only
    def on_button_1D_clicked(b):
        scan = expt.scans[0]
        Create_cell(code='CF.Plot_1D(nxs_filename=\''+scan.nxs+'\', recording_dir=expt.recording_dir,'+
                    'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)
        
    def on_button_fit_erf_clicked(b):
        scan = expt.scans[0]
        Create_cell(code='CF.GaussianRepartition_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)  
        
    def on_button_fit_gau_clicked(b):
        scan = expt.scans[0]
        Create_cell(code='CF.Gaussian_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)   
        
    def on_button_vineyard_clicked(b):
        scan = expt.scans[0]
        Create_cell(code='expt.channel0 = CF.Extract_channel_Qc(nxs_filename=\''+scan.nxs+'\','+
                    'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                    'logx=False, logy=False, logz=False)',
                    position='below', celltype='code', is_print = True)
        
    # Next action   
    def on_button_next_clicked(b):
        #clear_output(wait=False)
        
        # Save the default value of booleans
        
        Delete_current_cell()
        
        Create_cell(code='FF.Choose_action(expt)',
                    position ='at_bottom', celltype='code', is_print=False)        
        
    def on_button_markdown_clicked(b):
        """
        Insert a markdown cell above the current cell.
        """ 
        Create_cell(code='**Insert your comment (double-click here and replace current text).**',
                    position ='above', celltype='markdown', is_print=True)
        
       
    # Display the widgets    
    # ADD HERE A CALL TO YOUR BUTTON
    button_fit_gau = widgets.Button(description="Fit with gaussian")
    button_fit_gau.on_click(on_button_fit_gau_clicked)
    
    button_fit_erf = widgets.Button(description="Fit with erf")
    button_fit_erf.on_click(on_button_fit_erf_clicked)
    
    button_vineyard = widgets.Button(description="Extract Vineyard")
    button_vineyard.on_click(on_button_vineyard_clicked)
    
    button_1D = widgets.Button(description="Add plot to report")
    button_1D.on_click(on_button_1D_clicked)
    
    button_GIXD = widgets.Button(description="Plot GIXD")
    button_GIXD.on_click(on_button_GIXD_clicked)
    
    button_true_GIXD = widgets.Button(description="Plot true GIXD")
    button_true_GIXD.on_click(on_button_true_GIXD_clicked)
    
    button_fluo = widgets.Button(description="Plot fluo")
    button_fluo.on_click(on_button_fluo_clicked)
    
    button_isotherm = widgets.Button(description="Plot isotherm")
    button_isotherm.on_click(on_button_isotherm_clicked)
    
    button_pilatus = widgets.Button(description="Plot pilatus")
    button_pilatus.on_click(on_button_pilatus_clicked)

    button_GIXS = widgets.Button(description="Plot GIXS")
    button_GIXS.on_click(on_button_GIXS_clicked)
        
    button_next = widgets.Button(description="Next action")
    button_next.on_click(on_button_next_clicked)
    
    button_markdown = widgets.Button(description="Insert comment")
    button_markdown.on_click(on_button_markdown_clicked)
      
    if len(expt.scans)==1:
        # Options for single scan analysis only
        
        # Buttons for general treatment
        buttons0 = widgets.HBox([button_fit_gau, button_fit_erf, button_1D])
        display(buttons0)

        # Set up an interactive 1D plot
        Set_interactive_1D(expt.scans[0])

    else:
        print("Selected scans:")
        for scan in expt.scans:
              print('%s: %s'%(scan.nxs,scan.command))
        
    # Buttons for specific treatment
    buttons1 = widgets.HBox([button_GIXD, button_true_GIXD, button_fluo, button_isotherm])
    display(buttons1)
    
    buttons2 = widgets.HBox([button_pilatus, button_GIXS])
    display(buttons2)

    buttons3 = widgets.HBox([button_vineyard, button_markdown, button_next])
    display(buttons3)
    
    
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
                     description='Horizontal focalisation:',
                     layout=long_layout,
                     style=style)

    vertical_foc = widgets.Checkbox(
                   value=True,
                   description='Vertical focalisation:',
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
                    +'\t * %s %s'%(horizontal_size.description,horizontal_size.value)+'\n'
                    +'\t * %s %s'%(vertical_size.description,vertical_size.value))

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
        
        Create_cell(code='# Experimental setup', position ='below', celltype='markdown', is_print=True)
        
    button = widgets.Button(description="Print form")
    out = widgets.Output()
    display(button,out)
    button.on_click(on_button_clicked)
