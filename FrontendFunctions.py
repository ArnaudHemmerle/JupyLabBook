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
import math
import ipysheet

__version__ = '1.0'

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
    """
    Take an object from the class Experiment.
    1) Print the versions of custom modules
    2) Check if files and folders exist
    3) Create the first cell
    """
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
        
    Create_cell(code='FF.Choose_action(expt)', position ='at_bottom', celltype='code', is_print=False)   


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
    scan_found = False
    scan.command = 'No command found'
    
    for log_file in expt.list_logs_files:
        with open(expt.logs_dir+log_file) as f:
            for line in f:
                if "#" not in line: temp = line
                if scan.id in line:
                    # Remove the line jump
                    scan.command = temp.replace('\n','')
                    scan_found = True
            f.close()

        if scan_found:
            break
 

def Find_absorbers_in_logs(scan, expt):
    """
    Find the absorbers of the scan in the logs.
    """
    scan_found = False
    absorbers = 'No abs found' 
    
    for log_file in expt.list_logs_files:
        with open(expt.logs_dir+log_file) as f:
            for line in f:
                if "Aborbers" in line: 
                    temp = line.split(': ')[-1]
                if "Absorbers" in line: 
                    temp = line.split('- ')[-1]
                if scan.id in line:
                    # Remove the line jump
                    absorbers = temp.replace('\n','')
                    scan_found = True
            f.close()

        if scan_found:
            break   
            
    return absorbers
    
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
        #clear_output(wait=False)
        
        Create_cell(code='FF.Choose_treatment(expt)',
                    position ='below', celltype='code', is_print=False)

        if len(expt.scans)==1:
            Create_cell(code='### '+expt.scans[0].id+': '+expt.scans[0].command,
                        position ='below', celltype='markdown', is_print=True)
            
        Delete_current_cell()

        
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
        
    def on_button_convert_logs_clicked(b):
        """
        Convert the logs in human-readable files. 
        """
        Convert_logs(expt)
        
    def on_button_calibthetaz_clicked(b):
        """
        Create a cell for the calibration of thetaz.
        """
        #clear_output(wait=False)
        
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
        
        Delete_current_cell()
        
    def on_button_form_clicked(b):
        """
        Create a form.
        """
        #clear_output(wait=False)
        Create_cell(code='FF.Create_form()',position ='below', celltype='code', is_print=False)
        Create_cell(code='FF.Choose_action(expt)', position ='at_bottom', celltype='code', is_print=False)        
         
        Delete_current_cell()
            
    def on_button_wm_clicked(b):
        """
        Print the motors positions from a log file.
        """
        Print_wm(expt)

    def on_button_commands_clicked(b):
        """
        Print list of commands from a log file.
        """
        Print_commands(expt)        
        
    def on_button_export_clicked(b):
        """
        Export the notebook to PDF.
        """
        
        print('Export in progress...')
        
        export_done = Export_nb_to_pdf(expt.notebook_name)
        
        if export_done:
            print('Notebook exported to %s.pdf'%expt.notebook_name.split('.')[0])
        else:
            print("There was something wrong with the export to pdf.")
            print("Did you rename the Notebook? If yes:")
            print("1) Change the value of expt.notebook_name in the first cell (top of the Notebook).")
            print("2) Re-execute the first cell.")
            print("3) Try to export the pdf again in the last cell (bottom of the Notebook).")

                
    def on_button_markdown_clicked(b):
        """
        Insert a markdown cell below the current cell.
        """ 
        
        Delete_current_cell()
        
        Create_cell(code='', position ='below', celltype='markdown', is_print=True, is_execute=False)
    
        Create_cell(code='FF.Choose_action(expt)', position ='at_bottom', celltype='code', is_print=False)
        
    def on_button_script_clicked(b):
        """
        Insert a script as a markdown cell.
        """ 
        Print_script(expt)
    
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

    # Click to export human-readable logs
    button_convert_logs = widgets.Button(description="Convert logs")
    button_convert_logs.on_click(on_button_convert_logs_clicked)
      
    # Click to fill the form
    button_form = widgets.Button(description="Fill form")
    button_form.on_click(on_button_form_clicked)
    
    # Click to print motors
    button_wm = widgets.Button(description="Insert positions")
    button_wm.on_click(on_button_wm_clicked)
    
    # Click to print a list of commands
    button_commands = widgets.Button(description="Insert commands")
    button_commands.on_click(on_button_commands_clicked)
    
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
    w_print_scans = widgets.interact(selection_scans,
                                     nxs_files = widgets.SelectMultiple(
                                                 options=expt.list_nxs_files,
                                                 rows=10)
                                    )
    w_print_scans.widget.children[0].description = 'Next scan(s):'
    w_print_scans.widget.children[0].layout = {'width': '400px'}
    
    buttons1 = widgets.HBox([button_form, button_calibthetaz, button_convert_logs, button_export])
    display(buttons1)

    buttons2 = widgets.HBox([button_markdown, button_script, button_wm, button_commands])
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

def Create_cell(code='', position='below', celltype='markdown', is_print = False, is_execute = True):
    """Create a cell in the IPython Notebook.
    code: unicode, Code to fill the new cell with.
    celltype: unicode, Type of cells "code" or "markdown".
    position: unicode, Where to put the cell "below" or "at_bottom"
    is_print: boolean, To decide if the cell is printed in the pdf report
    The code requires direct use of Javascript and is thus not usable in Jupyter Lab.
    """

    # Delay to ensure unique id
    time.sleep(0.1)
    
    encoded_code = (base64.b64encode(code.encode())).decode()

    # Create a unique id based on epoch time
    display_id = int(time.time()*1e9)

    js_code = """var cell = IPython.notebook.insert_cell_{0}("{1}");
              cell.set_text(atob("{2}"));
              """
    if not is_print: js_code += """cell.metadata.tags = ['notPrint']
                                """
        
    if is_execute: js_code += """cell.execute();
                                """
    display(Javascript(js_code.format(position, celltype, encoded_code)),display_id=display_id)

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
    Allow the user to choose the treatment to be applied, via widgets.
    """
    
    # Styling options for widgets
    style = {'description_width': 'initial'}
    tiny_layout = widgets.Layout(width='150px', height='40px')
    short_layout = widgets.Layout(width='200px', height='40px')
    medium_layout = widgets.Layout(width='250px', height='40px')
    large_layout = widgets.Layout(width='300px', height='40px')
    
    # Define the function called when clicking the button
    # DEFINE HERE A FUNCTION TO CREATE A CELL CALLING YOUR CUSTOM FUNCTION

    def on_button_GIXD_clicked(b):

        # logx
        try: value = expt.GIXD_logx
        except: value = False
        w_GIXD_logx = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='log x')
        
        # logy
        try: value = expt.GIXD_logy
        except: value = False
        w_GIXD_logy = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='log y')
        
        # GIXD_logz
        try: value = expt.GIXD_logz
        except: value = False
        w_GIXD_logz = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='log z')
        
        # channel0
        try: value = expt.channel0
        except: value = 640
        w_channel0 = widgets.FloatText(value=value, style=style, layout=short_layout, description='Vineyard (chan)')

        # thetazfactor
        try: value = expt.thetazfactor
        except: value = 0.000243
        w_thetazfactor = widgets.FloatText(value=value, style=style, layout=large_layout,
                                           description='thetazfactor (rad/chan)')
        
        # wavelength
        try: value = expt.wavelength
        except: value = 0.155
        w_wavelength = widgets.FloatText(value=value, style=style, layout=short_layout, description='wavelength (nm)')

        # thetac
        try: value = expt.thetac
        except: value = 0.0028
        w_thetac = widgets.FloatText(value=value, style=style, layout=short_layout, description='thetac (rad)')

        # binsize
        try: value = expt.binsize
        except: value = 10
        w_binsize = widgets.IntText(value=value, style=style, layout=tiny_layout, description='binsize')

        # computeqz
        try: value = expt.computeqz
        except: value = True
        w_computeqz = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Compute qz')
        
        # nblevels
        try: value = expt.nblevels
        except: value = 50
        w_nblevels = widgets.IntText(value=value, style=style, layout=tiny_layout, description='nblevels')
        
        # moytocreate
        try: value = expt.moytocreate_str
        except: value = '10, 20, 40'
        w_moytocreate_str = widgets.Text(value=value, description='moy to create', style=style, layout = short_layout)
            
        # show_data_stamps
        try: value = expt.show_data_stamps
        except: value = False
        w_show_data_stamps = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print sensors')

        # verbose
        try: value = expt.verbose
        except: value = False
        w_verbose = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print scan info')

        # show_absorbers
        try: value = expt.show_absorbers
        except: value = False
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print abs')        
        
        # GIXD_cmap
        try: value = expt.GIXD_cmap
        except: value = 'jet'
        w_GIXD_cmap = widgets.Select(value=value, style=style, rows=5, description='cmap',
                                options=['viridis', 'jet', 'Greys', 'cividis', 'hot'])        
            
        display(widgets.HBox([w_show_data_stamps, w_verbose, w_show_absorbers,
                              w_GIXD_logx, w_GIXD_logy, w_GIXD_logz, w_GIXD_cmap]))        
        display(widgets.HBox([w_binsize, w_nblevels, w_moytocreate_str, w_channel0, w_computeqz]))
        display(widgets.HBox([w_wavelength, w_thetac, w_thetazfactor]))
                
        
        def on_button_plot_clicked(b):

            # Pass current values as default values
            expt.GIXD_logx = w_GIXD_logx.value
            expt.GIXD_logy = w_GIXD_logy.value
            expt.GIXD_logz = w_GIXD_logz.value
            expt.channel0 = w_channel0.value
            expt.thetazfactor = w_thetazfactor.value
            expt.wavelength = w_wavelength.value
            expt.thetac = w_thetac.value
            expt.binsize = w_binsize.value
            expt.computeqz = w_computeqz.value
            expt.nblevels = w_nblevels.value
            expt.moytocreate_str = w_moytocreate_str.value
            expt.show_data_stamps = w_show_data_stamps.value
            expt.verbose = w_verbose.value
            expt.show_absorbers = w_show_absorbers.value
            expt.GIXD_cmap = w_GIXD_cmap.value
            
            # Convert moytocreate_str into a list
            list_moytocreate = [int(expt.moytocreate_str.split(',')[i]) for i in range(len(expt.moytocreate_str.split(',')))]
                        
            for scan in expt.scans:
                
                # Extract absorbers
                if expt.show_absorbers:
                    absorbers = Find_absorbers_in_logs(scan, expt)
                else:
                    absorbers = ''
                             
                Create_cell(code='CF.Extract_GIXD(nxs_filename=\''+scan.nxs+'\','+
                            'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                            'logx='+str(expt.GIXD_logx)+','+
                            'logy='+str(expt.GIXD_logy)+','+
                            'logz='+str(expt.GIXD_logz)+','+
                            'channel0='+str(expt.channel0)+','+
                            'thetazfactor='+str(expt.thetazfactor)+','+
                            'wavelength='+str(expt.wavelength)+','+
                            'thetac='+str(expt.thetac)+','+
                            'binsize='+str(expt.binsize)+','+
                            'computeqz='+str(expt.computeqz)+','+
                            'nblevels='+str(expt.nblevels)+','+
                            'moytocreate='+str(list_moytocreate)+','+
                            'show_data_stamps='+str(expt.show_data_stamps)+','+
                            'verbose='+str(expt.verbose)+','+
                            'absorbers='+'\''+str(absorbers)+'\''+','+
                            'cmap=\''+str(expt.GIXD_cmap)+'\')',
                            position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)

            # Do as if the button next was clicked
            on_button_next_clicked(b)

        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        display(button_plot)

        

    def on_button_pilatus_clicked(b):
        
        # pilatus_logz
        try: value = expt.pilatus_logz
        except: value = True
        w_pilatus_logz = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='log z')
      
        # show_data_stamps
        try: value = expt.show_data_stamps
        except: value = False
        w_show_data_stamps = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print sensors')

        # verbose
        try: value = expt.verbose
        except: value = False
        w_verbose = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print scan info')

        # show_absorbers
        try: value = expt.show_absorbers
        except: value = False
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print abs')
        
        # pilatus_cmap
        try: value = expt.pilatus_cmap
        except: value = 'Greys'
        w_pilatus_cmap = widgets.Select(value=value, style=style, rows=5, description='cmap',
                                        options=['viridis', 'jet', 'Greys', 'cividis', 'hot'])
        
        # xmin
        try: value = expt.xmin
        except: value = 0
        w_xmin = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='x min (pix)')        

        # xmax
        try: value = expt.xmax
        except: value = 980
        w_xmax = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='x max (pix)')         
 
        # ymin
        try: value = expt.ymin
        except: value = 0
        w_ymin = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='y min (pix)')        

        # ymax
        try: value = expt.ymax
        except: value = 1042
        w_ymax = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='y max (pix)')   

        display(widgets.HBox([w_show_data_stamps, w_verbose, w_show_absorbers, w_pilatus_logz, w_pilatus_cmap]))  
        display(widgets.HBox([w_xmin, w_xmax])) 
        display(widgets.HBox([w_ymin, w_ymax])) 

        def on_button_plot_clicked(b):

            # Pass current values as default values
            expt.pilatus_logz = w_pilatus_logz.value
            expt.show_data_stamps = w_show_data_stamps.value
            expt.verbose = w_verbose.value
            expt.show_absorbers = w_show_absorbers.value
            expt.pilatus_cmap = w_pilatus_cmap.value
            expt.xmin = w_xmin.value
            expt.xmax = w_xmax.value
            expt.ymin = w_ymin.value
            expt.ymax = w_ymax.value

            for scan in expt.scans:
                
                # Extract absorbers
                if expt.show_absorbers:
                    absorbers = Find_absorbers_in_logs(scan, expt)
                else:
                    absorbers = ''                

                Create_cell(code='CF.Extract_pilatus_sum(nxs_filename=\''+scan.nxs+'\','+ 
                           'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                            'logz='+str(expt.pilatus_logz)+','+
                            'xmin='+str(expt.xmin)+','+
                            'xmax='+str(expt.xmax)+','+
                            'ymin='+str(expt.ymin)+','+
                            'ymax='+str(expt.ymax)+','+                            
                            'show_data_stamps='+str(expt.show_data_stamps)+','+
                            'verbose='+str(expt.verbose)+','+
                            'absorbers='+'\''+str(absorbers)+'\''+','+
                            'cmap=\''+str(expt.pilatus_cmap)+'\''+')',
                            position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)     

            # Do as if the button next was clicked
            on_button_next_clicked(b)                    
                    
        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        display(button_plot)  
       
    def on_button_GIXS_clicked(b):
        
        # Checkboxes for options       

        # GIXS_logz
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
        w_thetai = widgets.FloatText(value=value, style=style, layout=short_layout, description='thetai (rad)')
        
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

        # force_gamma_delta 
        try: value = expt.force_gamma_delta
        except: value = False
        w_force_gamma_delta = widgets.Checkbox(value=value, style=style, layout=short_layout, description='Force Gamma&Delta')        

        # fgamma
        try: value = expt.fgamma
        except: value = 0.
        w_fgamma = widgets.FloatText(value=value, style=style, layout=short_layout, description='Forced Gamma')        

        # fdelta
        try: value = expt.fdelta
        except: value = 0.
        w_fdelta = widgets.FloatText(value=value, style=style, layout=short_layout, description='Forced Delta')         
        
        # pixel_size 
        try: value = expt.pixel_size
        except: value = 0.172
        w_pixel_size = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='Pixel size (um)')

        # xmin
        try: value = expt.xmin
        except: value = 0.
        w_xmin = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='x min')   
 
        # xmax
        try: value = expt.xmax
        except: value = 1.
        w_xmax = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='x max')  
        
        # ymin
        try: value = expt.ymin
        except: value = 0.
        w_ymin = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='y min')   
 
        # ymax
        try: value = expt.ymax
        except: value = 1.
        w_ymax = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='y max')          
        
        # show_data_stamps
        try: value = expt.show_data_stamps
        except: value = False
        w_show_data_stamps = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print sensors')
       
        # verbose
        try: value = expt.verbose
        except: value = False
        w_verbose = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print scan info')

        # show_absorbers
        try: value = expt.show_absorbers
        except: value = False
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print abs')
        
        # GIXS_cmap
        try: value = expt.GIXS_cmap
        except: value = 'viridis'
        w_GIXS_cmap = widgets.Select(value=value, style=style, rows=5, description='cmap',
                                options=['viridis', 'jet', 'Greys', 'cividis', 'hot'])
     
                 
        display(widgets.HBox([w_show_data_stamps, w_verbose, w_show_absorbers, w_GIXS_logz, w_GIXS_cmap, w_pixel_size]))        
        display(widgets.HBox([w_wavelength, w_distance, w_thetai, w_pixel_PONI_x, w_pixel_PONI_y]))
        display(widgets.HBox([w_force_gamma_delta, w_fgamma, w_fdelta]))
        display(widgets.HBox([w_xmin, w_xmax])) 
        display(widgets.HBox([w_ymin, w_ymax])) 
        
        def on_button_plot_clicked(b):
            
            # Pass current values as default values
            expt.GIXS_logz = w_GIXS_logz.value
            expt.wavelength = w_wavelength.value
            expt.thetai = w_thetai.value
            expt.distance = w_distance.value
            expt.pixel_PONI_x = w_pixel_PONI_x.value
            expt.pixel_PONI_y = w_pixel_PONI_y.value
            expt.pixel_size = w_pixel_size.value
            expt.xmin = w_xmin.value
            expt.xmax = w_xmax.value
            expt.ymin = w_ymin.value
            expt.ymax = w_ymax.value
            expt.force_gamma_delta = w_force_gamma_delta.value
            expt.fgamma = w_fgamma.value
            expt.fdelta = w_fdelta.value
            expt.show_data_stamps = w_show_data_stamps.value
            expt.verbose = w_verbose.value
            expt.show_absorbers = w_show_absorbers.value
            expt.GIXS_cmap = w_GIXS_cmap.value

            
            for scan in expt.scans:

                # Extract absorbers
                if expt.show_absorbers:
                    absorbers = Find_absorbers_in_logs(scan, expt)
                else:
                    absorbers = ''                
                
                Create_cell(code='CF.Extract_GIXS(nxs_filename=\''+scan.nxs+'\','+ 
                       'working_dir=expt.working_dir, recording_dir=expt.recording_dir,'+
                        'logz='+str(expt.GIXS_logz)+','+
                        'wavelength='+str(expt.wavelength)+','+
                        'thetai='+str(expt.thetai)+','+
                        'distance='+str(expt.distance)+','+
                        'pixel_PONI_x='+str(expt.pixel_PONI_x)+','+
                        'pixel_PONI_y='+str(expt.pixel_PONI_y)+','+
                        'pixel_size='+str(expt.pixel_size)+','+  
                        'xmin='+str(expt.xmin)+','+
                        'xmax='+str(expt.xmax)+','+
                        'ymin='+str(expt.ymin)+','+
                        'ymax='+str(expt.ymax)+','+
                        'show_data_stamps='+str(expt.show_data_stamps)+','+
                        'force_gamma_delta='+str(expt.force_gamma_delta)+','+
                        'fgamma='+str(expt.fgamma)+','+
                        'fdelta='+str(expt.fdelta)+','+
                        'verbose='+str(expt.verbose)+','+ 
                        'absorbers='+'\''+str(absorbers)+'\''+','+
                        'cmap=\''+str(expt.GIXS_cmap)+'\''+')',
                        position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    Create_cell(code='### '+scan.id+': '+scan.command,
                                position ='below', celltype='markdown', is_print=True)            

            # Do as if the button next was clicked
            on_button_next_clicked(b)
            
        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        display(button_plot)

    def on_button_XRF_clicked(b):
        
        # Checkboxes for options       
       
        # show_data_stamps
        try: value = expt.show_data_stamps
        except: value = False
        w_show_data_stamps = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Print sensors')

        # verbose
        try: value = expt.verbose
        except: value = False
        w_verbose = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Print scan info')

        # show_absorbers
        try: value = expt.show_absorbers
        except: value = False
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Print abs')
        
        # fastextract
        try: value = expt.fastextract
        except: value = True
        w_fastextract = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Fast extract')

        # plot_spectrogram
        try: value = expt.plot_spectrogram
        except: value = True
        w_plot_spectrogram = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Plot spectrogram')
        
        # plot_first_last
        try: value = expt.plot_first_last
        except: value = True
        w_plot_first_last = widgets.Checkbox(value=value, style=style,
                                             layout = short_layout, description='Plot first&last spectrums')
        
        # plot_sum
        try: value = expt.plot_sum
        except: value = True
        w_plot_sum = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Plot sum of spectrums')
        
        # logz
        try: value = expt.XRF_logz
        except: value = True
        w_XRF_logz = widgets.Checkbox(value=value, style=style, layout = short_layout, description='log z')
        
        # list_elems
        try: value = expt.elems_str
        except: value = '4'
        w_elems_str = widgets.Text(value=value, description='Elements', style=style, layout = short_layout)

        # first_channel
        try: value = expt.first_channel
        except: value = 0
        w_first_channel = widgets.IntText(value=value, description='First channel', style=style, layout = short_layout)
        
        # last_channel
        try: value = expt.last_channel
        except: value = 2047
        w_last_channel = widgets.IntText(value=value, description='Last channel', style=style, layout = short_layout)

        # use_eV
        try: value = expt.use_eV
        except: value = False
        w_use_eV = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Use eV')

        # gain
        try: value = expt.gain
        except: value = 10.
        w_gain = widgets.FloatText(value=value, description='Gain', style=style, layout = short_layout)        

        # eV0
        try: value = expt.eV0
        except: value = 0.
        w_eV0 = widgets.FloatText(value=value, description='eV0', style=style, layout = short_layout)
        
        # is_identify_peaks
        expt.is_identify_peaks = False
        
        display(widgets.HBox([w_show_data_stamps, w_verbose, w_show_absorbers, w_fastextract]))
        display(widgets.HBox([w_plot_spectrogram, w_plot_first_last, w_plot_sum]))
        display(widgets.HBox([w_XRF_logz, w_elems_str, w_first_channel, w_last_channel]))
        display(widgets.HBox([w_use_eV, w_gain, w_eV0]))

        def on_button_identify_peaks_clicked(b):
            
            sheet = ipysheet.easy.sheet(columns=3, rows=20 ,column_headers = ['Name','Position','Use?(y/n)'])
            
            # Fill the sheet with previous values
            try: to_fill = expt.arr_peaks_full
            except: to_fill = np.array([[None,None,None] for i in range(20)])
                
            # ipysheet does not work correctly with None entries
            # It is necessary to fill first the cells with something
            for i in range(3):
                ipysheet.easy.column(i,  to_fill[:,i])

            def on_button_validate_clicked(b):
                
                expt.is_identify_peaks = True
                
                # Get values from the sheet
                expt.arr_peaks_full = ipysheet.numpy_loader.to_array(ipysheet.easy.current())
                
                # Remove the empty lines and make array of tuples (name, eV)
                expt.arr_peaks = [(elem[0],elem[1],elem[2]) for elem in expt.arr_peaks_full if elem[0]!=None]
                
                # Send only the lines with y in the third column
                expt.arr_peaks = [elem[0:2] for elem in expt.arr_peaks if elem[2]=='y']
                
                for scan in expt.scans:                                                              
                
                    print("Peaks on scan %s"%scan.nxs)
                    # Extract and plot the XRF with the peaks when validate is clicked
                    CF.Extract_XRF(nxs_filename=scan.nxs,
                                   working_dir=expt.working_dir,
                                   recording_dir=expt.recording_dir,
                                   logz=w_XRF_logz.value,
                                   list_elems=w_elems_str.value,
                                   first_channel=w_first_channel.value,
                                   last_channel=w_last_channel.value,
                                   use_eV=w_use_eV.value,
                                   gain=w_gain.value,
                                   eV0=w_eV0.value,
                                   arr_peaks=expt.arr_peaks,
                                   show_data_stamps=False,
                                   verbose=False,
                                   absorbers='',
                                   fast=w_fastextract.value,
                                   plot_spectrogram=False,
                                   plot_first_last=False,
                                   plot_sum=True)
                          
                
            button_validate = widgets.Button(description="Validate peaks")
            button_validate.on_click(on_button_validate_clicked)

            display(button_validate)
            display(sheet)
            
        def on_button_plot_clicked(b):
            
            if not expt.is_identify_peaks: expt.arr_peaks = [(None,None)]
            
            # Pass current values as default values
            expt.XRF_logz = w_XRF_logz.value
            expt.elems_str = w_elems_str.value
            expt.first_channel = w_first_channel.value
            expt.last_channel = w_last_channel.value
            expt.use_eV = w_use_eV.value
            expt.gain = w_gain.value
            expt.eV0 = w_eV0.value
            expt.show_data_stamps = w_show_data_stamps.value
            expt.verbose = w_verbose.value
            expt.show_absorbers = w_show_absorbers.value
            expt.fastextract = w_fastextract.value
            expt.plot_spectrogram = w_plot_spectrogram.value
            expt.plot_first_last = w_plot_first_last.value
            expt.plot_sum = w_plot_sum.value
           
            # Convert elems_str into a list
            list_elems = [int(expt.elems_str.split(',')[i]) for i in range(len(expt.elems_str.split(',')))]
            
            for scan in expt.scans:
                
                # Extract absorbers
                if expt.show_absorbers:
                    absorbers = Find_absorbers_in_logs(scan, expt)
                else:
                    absorbers = ''
                                            
                Create_cell(code='CF.Extract_XRF('+
                           'nxs_filename=\''+scan.nxs+'\','+ 
                           'working_dir=expt.working_dir,recording_dir=expt.recording_dir,'+
                           'logz='+str(expt.XRF_logz)+','+
                           'list_elems='+str(list_elems)+','+
                           'first_channel='+str(expt.first_channel)+','+
                           'last_channel='+str(expt.last_channel)+','+
                           'use_eV='+str(expt.use_eV)+','+
                           'gain='+str(expt.gain)+','+
                           'eV0='+str(expt.eV0)+','+
                           'arr_peaks='+str(expt.arr_peaks)+','+
                           'show_data_stamps='+str(expt.show_data_stamps)+','+
                           'verbose='+str(expt.verbose)+','+
                           'absorbers='+'\''+str(absorbers)+'\''+','+
                           'fast='+str(expt.fastextract)+','+
                           'plot_spectrogram='+str(expt.plot_spectrogram)+','+
                           'plot_first_last='+str(expt.plot_first_last)+','+
                           'plot_sum='+str(expt.plot_sum)+
                           ')',
                           position='below', celltype='code', is_print = True)  

                if len(expt.scans)>1:
                    Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)  

            # Do as if the button next was clicked
            on_button_next_clicked(b)                    

        button_identify_peaks = widgets.Button(description="Identify peaks")
        button_identify_peaks.on_click(on_button_identify_peaks_clicked)            
            
        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        
        display(widgets.HBox([button_identify_peaks, button_plot]))
    
    def on_button_isotherm_clicked(b):
            
        # show_data_stamps
        try: value = expt.show_data_stamps
        except: value = False
        w_show_data_stamps = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print sensors')

        # verbose
        try: value = expt.verbose
        except: value = False
        w_verbose = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print scan info')
        
        # fastextract
        try: value = expt.fastextract
        except: value = True
        w_fastextract = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Fast extract')

        
        display(widgets.HBox([w_show_data_stamps, w_verbose, w_fastextract]))   
        
        def on_button_plot_clicked(b):

            # Pass current values as default values
            expt.show_data_stamps = w_show_data_stamps.value
            expt.verbose = w_verbose.value
            expt.fastextract = w_fastextract.value

            for scan in expt.scans:

                Create_cell(code='CF.Plot_isotherm(nxs_filename=\''+scan.nxs+'\','+
                            'working_dir=expt.working_dir,recording_dir=expt.recording_dir,'+
                            'show_data_stamps='+str(w_show_data_stamps.value)+
                            ', verbose='+str(w_verbose.value)+', '+
                            'fast='+str(w_fastextract.value)+')',
                            position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)

            # Do as if the button next was clicked
            on_button_next_clicked(b)                    
                    
                    
        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        
        display(button_plot)                    
                    
    # Actions relevant for single scan analysis only
    def on_button_1D_clicked(b):            
        scan = expt.scans[0]
        Create_cell(code='CF.Plot_1D(nxs_filename=\''+scan.nxs+'\','+
                    'working_dir=expt.working_dir,recording_dir=expt.recording_dir,'+
                    'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)
        
        # Do as if the button next was clicked
        on_button_next_clicked(b) 
        
    def on_button_fit_erf_clicked(b):
        scan = expt.scans[0]
        Create_cell(code='CF.GaussianRepartition_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)  

        # Do as if the button next was clicked
        on_button_next_clicked(b)        
        
    def on_button_fit_gau_clicked(b):
        scan = expt.scans[0]
        Create_cell(code='CF.Gaussian_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)   
        
        # Do as if the button next was clicked
        on_button_next_clicked(b)  
        
    def on_button_vineyard_clicked(b):
        scan = expt.scans[0]
        Create_cell(code='expt.channel0 = CF.Extract_channel_Qc(nxs_filename=\''+scan.nxs+'\','+
                    'working_dir=expt.working_dir, recording_dir=expt.recording_dir, '+
                    'logx=False, logy=False, logz=False)',
                    position='below', celltype='code', is_print = True)
 
        # Do as if the button next was clicked
        on_button_next_clicked(b)  

    # Next action   
    def on_button_next_clicked(b):
        #clear_output(wait=False)
        
        Delete_current_cell()
        
        Create_cell(code='FF.Choose_action(expt)',
                    position ='at_bottom', celltype='code', is_print=False)        
        
    def on_button_markdown_clicked(b):
        """
        Insert a markdown cell below the current cell.
        """ 
        Delete_current_cell()
        
        Create_cell(code='', position ='below', celltype='markdown', is_print=True, is_execute=False)
    
        Create_cell(code='FF.Choose_treatment(expt)', position ='at_bottom', celltype='code', is_print=False)
       
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
    
    button_XRF = widgets.Button(description="Plot XRF")
    button_XRF.on_click(on_button_XRF_clicked)
    
    button_isotherm = widgets.Button(description="Plot isotherm")
    button_isotherm.on_click(on_button_isotherm_clicked)
    
    button_pilatus = widgets.Button(description="Plot Pilatus")
    button_pilatus.on_click(on_button_pilatus_clicked)

    button_GIXS = widgets.Button(description="Plot GIXS")
    button_GIXS.on_click(on_button_GIXS_clicked)
        
    button_next = widgets.Button(description="Next action")
    button_next.on_click(on_button_next_clicked)
    
    button_markdown = widgets.Button(description="Insert comment")
    button_markdown.on_click(on_button_markdown_clicked)

    # Buttons for general treatment
    buttons0 = widgets.HBox([button_markdown, button_next])
    display(buttons0)    
    
    if len(expt.scans)==1:
        # Option for single scan analysis only

        # Set up an interactive 1D plot
        Set_interactive_1D(expt.scans[0])

    else:
        print("Selected scans:")
        for scan in expt.scans:
              print('%s: %s'%(scan.nxs,scan.command))
        
    # Buttons for specific treatment
    buttons1 = widgets.HBox([button_vineyard, button_fit_gau, button_fit_erf, button_1D])
    display(buttons1)
    
    buttons2 = widgets.HBox([button_GIXD, button_XRF, button_isotherm, button_pilatus, button_GIXS])
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
        
        # Remove the widget when done
        Delete_current_cell()
        
    button = widgets.Button(description="Print form")
    out = widgets.Output()
    display(button,out)
    button.on_click(on_button_clicked)

    
def Print_wm(expt):
    """
    Print the positions of the motors/sensors from extracting the command wm in a log file.
    """

    try:
        default_value = expt.default_wm_log_file
    except:
        default_value = expt.list_logs_files[0]
    
    w_select_log = widgets.Dropdown(
                 options=expt.list_logs_files,
                 value=default_value,
                 layout=widgets.Layout(width='400px'),
                 style={'description_width': 'initial'})

    def on_button_select_clicked(b):
        
        # Keep history of the chosen file
        expt.default_wm_log_file = w_select_log.value
        
        # Construct the list of wm in the log
        list_wm = []
        log_file =  w_select_log.value
        with open(expt.logs_dir+log_file) as f:
            for line in f:
                if "# " in line: temp = line
                if ("wm " in line and "pwm" not in line and "ERROR" not in line):
                    list_wm.append(temp.replace('\n','')+'; '+line.replace('\n',''))
        f.close()

        if list_wm == []:
            print(PN._RED+"No positions in the log file: "+log_file+PN._RESET)
            print(PN._RED+"Choose another log file."+PN._RESET)
            return

        w_select_wm = widgets.Dropdown(
                     options=list_wm,
                     value=list_wm[-1],
                     layout=widgets.Layout(width='700px'),
                     style={'description_width': 'initial'})

        display(w_select_wm)

        def on_button_print_clicked(b):
             
            # Extract the date of the wm (used as an identifier to be found in the log)
            date_wm = w_select_wm.value.split('; ')[0]

            # Find the right lines in the log
            list_sensors = []
            is_date_found = False
            count = 0
            with open(expt.logs_dir+log_file) as f:
                for line in f:
                    if date_wm in line: 
                        is_date_found = True
                    if (is_date_found and '-------' in line):
                        count += 1    
                    if (count == 2 and '-------' not in line):
                        # append name, value, unit
                        list_sensors.append(line.split()[2:5])

                        # if unit is not present
                        if ('[' in line.split()[4]):  
                            list_sensors[-1][-1] = ' '

                        # if unit is 'No unit'
                        if ('No' in line.split()[4]):  
                            list_sensors[-1][-1] = ' '

            f.close()

            # Display as a markdown table
            max_elem_per_line = 6
            nb_elem = len(list_sensors)
            if nb_elem >= max_elem_per_line:
                nb_line = math.ceil(nb_elem/max_elem_per_line)
                nb_per_line = math.ceil(nb_elem/nb_line)
            else:
                nb_per_line = max_elem_per_line

            list_sensors_str = []
            for i in range(0,nb_elem,nb_per_line):
                list_sensors_cut = list_sensors[i:i+nb_per_line]

                tbw_str = '' 
                for sensor in list_sensors_cut:
                    tbw_str += '|'+str(sensor[0])
                tbw_str += '|'+' \n'

                for sensor in list_sensors_cut:
                    tbw_str += '|'+':-:'
                tbw_str += '|'+' \n'

                for sensor in list_sensors_cut:
                    tbw_str += '|'+str(sensor[1])
                tbw_str += '|'+' \n'

                for sensor in list_sensors_cut:
                    tbw_str += '|'+str(sensor[2])
                tbw_str += '|'

                list_sensors_str.append(tbw_str)

            for sensor_str in list_sensors_str[::-1]:
                # Create markdown cells with the tables
                Create_cell(code=sensor_str, position ='below', celltype='markdown', is_print=True)

            # Put title    
            Create_cell(code='### '+w_select_wm.value.split('; ')[1], position ='below', celltype='markdown', is_print=True)
            
            Delete_current_cell()
        
            Create_cell(code='FF.Choose_action(expt)',
                        position ='at_bottom', celltype='code', is_print=False)  
            
        button = widgets.Button(description="Insert positions")
        display(button)
        button.on_click(on_button_print_clicked)


    button = widgets.Button(description="Select log")
    button.on_click(on_button_select_clicked)

    display(widgets.HBox([w_select_log, button]))  
    
def Print_script(expt):
    """
    Collect info from a wm command in a log and add it to the report.
    """

    # Check if there is already a path for scripts directories
    # If not, use the recording directory
    try:
        path_to_dir_default = expt.path_to_dir
    except:
        path_to_dir_default = expt.recording_dir

    # Widget to write path
    w_path_to_dir = widgets.Text(
            value=path_to_dir_default,
            description='Scripts directory:',
            layout=widgets.Layout(width='800px', height='40px'),
            style={'description_width': 'initial'})

    def on_button_validate_path_clicked(b):
        """
        Validate the path of the script folder. Open selection for the script.
        """

        if not os.path.exists(w_path_to_dir.value):
            print(PN._RED+"Wrong folder name."+PN._RESET)
            print("")
            return

        # Pass the current value of the directory to the default one
        expt.path_to_dir = w_path_to_dir.value

        # Define the list of log files in the log directory
        list_scripts_files = [file for file in sorted(os.listdir(expt.path_to_dir)) if '.ipy' in file][::-1]

        if len(list_scripts_files) < 1:
            print(PN._RED+"There is no script in this folder."+PN._RESET)
            print("")
            return                

        # Widget to select script
        w_select_script = widgets.Dropdown(
             options=list_scripts_files,
             value=list_scripts_files[-1],
             layout=widgets.Layout(width='300px'),
             style={'description_width': 'Script:'})

        def on_button_insert_script_clicked(b):
            """
            Display a script in a markdown cell.
            Add scan numbers if asked.
            """ 

            # Get and insert the script
            path_to_dir = w_path_to_dir.value
            script_name = w_select_script.value
            text_file = open(path_to_dir+script_name)
            
            # Import original script
            script_init = text_file.read()
            
            # Split the different lines
            script_split = script_init.split('\n')
            
            
            # Check if the value of first script is a number
            # If not, do not add the counting
            
            if np.char.isnumeric(w_index_first_scan.value):

                
                ###################################
                # Add scan number
                ###################################

                # List of names which triggers an increment of the scan number
                list_trig = ['cont_regh', 'scan']

                count = int(w_index_first_scan.value)

                is_loop = False
                for i in range(len(script_split)):

                    if is_loop:

                        # Fist run on the loop to collect info
                        loop_length = 0
                        nb_scan_in_loop = 0
                        is_first_run = True
                        for k in range(i,len(script_split)):
                            if is_first_run: 
                                if ('    ' in script_split[k]) or ('\t' in script_split[k]):
                                    loop_length +=1
                                    if any(elem in script_split[k] for elem in list_trig):
                                        nb_scan_in_loop+=1
                                else:
                                    is_first_run = False

                        # Attribute the scan numbers
                        for j in range(repet_nb-1):
                            for k in range(i,i+loop_length):
                                if any(elem in script_split[k] for elem in list_trig):
                                        script_split[k] = script_split[k]+' #'+str(count)
                                        count+=1

                        is_loop = False

                    # Detect a loop
                    if ('in range' in script_split[i]):
                        repet_nb = int(''.join([str(s) for s in script_split[i] if s.isdigit()]))
                        is_loop = True

                    # Regular scan not in a loop
                    if any(elem in script_split[i] for elem in list_trig):
                        script_split[i] = script_split[i]+' #'+str(count)
                        count+=1

                ##################################
            
            
            # Re-create the script with added scan numbers
            script_modif ='\n'.join(script_split)
            
            text_file.close()
            
            
            code = '```python\n'+script_modif+'```'

            Create_cell(code='### '+ script_name, position ='above', celltype='markdown', is_print=True)
            Create_cell(code=code, position ='above', celltype='markdown', is_print=True)
            
            Delete_current_cell()
        
            Create_cell(code='FF.Choose_action(expt)',
                        position ='at_bottom', celltype='code', is_print=False)  

        
        # Get number of the first scan
        w_index_first_scan = widgets.Text(value='',
                                             style = {'description_width': 'initial'},
                                             layout = widgets.Layout(width='200px'),
                                             description='Index first scan')
        

        button_insert_script = widgets.Button(description="Insert script")
        button_insert_script.on_click(on_button_insert_script_clicked)

        display(widgets.HBox([w_select_script, w_index_first_scan, button_insert_script]))

    button_validate_path = widgets.Button(description="Validate path")
    button_validate_path.on_click(on_button_validate_path_clicked)
    display(widgets.HBox([w_path_to_dir, button_validate_path]))

    
############################################################
################## PRINTING COMMANDS #######################    
############################################################

    
def Print_commands(expt):
    """
    Print a list of commands from a log file. Add the scan numbers at the end of each line.
    """

    try:
        default_value = expt.default_wm_log_file
    except:
        default_value = expt.list_logs_files[0]

    w_select_log = widgets.Dropdown(
                 options=expt.list_logs_files,
                 value=default_value,
                 layout=widgets.Layout(width='400px'),
                 style={'description_width': 'initial'})

    def on_button_select_clicked(b):

        # Keep history of the chosen file
        expt.default_wm_log_file = w_select_log.value

        pathToFile = expt.logs_dir+w_select_log.value

        # Extract the formatted commands
        rlog_lines = Extract_commands(pathToFile)
        
        w_select_commands = widgets.SelectMultiple(options=rlog_lines,rows=30,
                                                   layout=widgets.Layout(width='800px'),
                                                   value=[rlog_lines[-1]])

        display(w_select_commands)

        def on_button_insert_commands_clicked(b):
            """
            Display selected commands and scan numbers in a markdown cell.
            """ 

            code = '```python\n'+''.join(w_select_commands.value)+'```'

            Create_cell(code=code, position ='above', celltype='markdown', is_print=True)

            Delete_current_cell()

            Create_cell(code='FF.Choose_action(expt)',
                        position ='at_bottom', celltype='code', is_print=False)  

        # Click to insert a script
        button_insert_commands = widgets.Button(description="Insert commands")
        button_insert_commands.on_click(on_button_insert_commands_clicked)

        display(button_insert_commands)

    button = widgets.Button(description="Select log")
    button.on_click(on_button_select_clicked)

    display(widgets.HBox([w_select_log, button]))  
 

def Convert_logs(expt):
    """
    Convert all the logs in a human-readable version with date and commands.
    """
    
    rlogs_dir = expt.working_dir+'readable_logs/'
    
    if not os.path.exists(rlogs_dir):
            os.mkdir(rlogs_dir)   
    
    print("Conversion of logs into a human-readable format in the folder:") 
    print(rlogs_dir)
    print(" ")
    
    # Convert all the logs
    for file in expt.list_logs_files:
        rlog_lines = Extract_commands(expt.logs_dir+file) 
        
        wfile = open(rlogs_dir+'_r'+file[1:],'w')
        for k in range(len(rlog_lines)):
            wfile.write(rlog_lines[k])
        wfile.close()
        
        print('_r'+file[1:]+': Done.')
        
    print(" ")    
    print('Finished!')
        
def Extract_commands(pathToFile):
    """
    Extract commands from a log file. Returns an array of all the formatted commands.
    """

    # Define the elements to be removed of the original log file
    remove_elem_0 = ['#', ']', '\n']
    remove_elem = [
                  'DevError', 'desc =', 'origin =', 'reason =', 'severity =',
                  'connection'
                   ]

    # List of possible dates
    date_list = ['# Mon', '# Tue', '# Wed', '# Thu', '# Fri', '# Sat', '# Sun']
    date = ''

    with open(pathToFile, 'r') as f:
        log_lines = f.readlines()

    rlog_lines = np.empty([], dtype ='<U1000')

    for i in range(len(log_lines)):
        if any(date in log_lines[i] for date in date_list):
            date = log_lines[i].replace('\n',' ')[2:]
        elif 'Scan File Name' in log_lines[i]:
            #scan_number = ''.join(filter(str.isdigit, log_lines[i]))
            scan_number = log_lines[i].split('_')[-1].replace('\n','')
            while scan_number[0]=='0':
                # Delete trailing zero
                scan_number = scan_number[1:]
            rlog_lines[-1] = rlog_lines[-1][:-1] +  ' #' + scan_number +'\n' 
        elif '*** ABORT ***' in log_lines[i]:
            rlog_lines = np.append(rlog_lines, 'ABORTED :'+rlog_lines[-1])
        elif 'Scan aborted' in log_lines[i]:
            rlog_lines = np.append(rlog_lines, 'SCAN ABORTED: '+rlog_lines[-1])

        elif log_lines[i][0] in remove_elem_0:
            pass
        elif any([elem in log_lines[i] for elem in remove_elem]):
            pass
        else:
            rlog_lines = np.append(rlog_lines, date+log_lines[i])

    # Remove first empty line
    rlog_lines = rlog_lines[1:]    

    return rlog_lines
