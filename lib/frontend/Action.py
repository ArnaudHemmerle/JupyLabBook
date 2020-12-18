from . import Utils
from lib.extraction.common import PyNexus as PN

import os
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
import time
import subprocess
import math

"""Frontend library for all the widgets concerning the actions in the notebook."""


class Scan:
    """ Class Scan is used to pass arguments concerning the current scan only."""
    def __init__(self):
        pass

    
def Check_and_init(expt):
    '''
    Check if files and folders exist, then create the first cell.

    Parameters
    ----------
    expt : object
        object from the class Experiment
    '''

    print("Data reduction will be saved in the folder:")
    if not os.path.exists(expt.working_dir):
        print(PN._RED+"Careful, the following folder does not exist and should be created:"+PN._RESET)
    print(expt.working_dir)
    print("")
        
    print("The original nexus files should be in the folder:")    
    if not os.path.exists(expt.recording_dir):
        print(PN._RED+"Careful, the following folder does not exist:"+PN._RESET)
    print(expt.recording_dir)
    print("")
    
    print("The log files should be in the folder:")
    if not os.path.exists(expt.logs_dir):
        print(PN._RED+"Careful, the following folder does not exist:"+PN._RESET)
    print(expt.logs_dir)
    print("")        
    
    if not os.path.exists(expt.notebook_name):
        print(PN._RED+"Careful, assign the correct notebook name to expt.notebook_name."+PN._RESET)
        print("")
        
    if not os.path.exists('latex_template.tplx'):
        print(PN._RED+"The following file does not exist:"+PN._RESET)
        print('latex_template.tplx')
        print("This file contains the template for generating PDF and should be placed in the same folder as the notebook.")
        print("") 
        
    Utils.Create_cell(code='FE.Action.Choose(expt)', position ='at_bottom', celltype='code', is_print=False)    

def Choose(expt):
    '''
    Choose the next action to do. Help the choice by printing the command corresponding to the scans displayed in the list.
    Create a list of objects from the class Scan, with info related to the scans which will be treated. 

    Parameters
    ----------
    expt : object
        object from the class Experiment
    '''    

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
        
        Utils.Create_cell(code='FE.Treatment.Choose(expt)',
                    position ='below', celltype='code', is_print=False)

        if len(expt.scans)==1:
            Utils.Create_cell(code='### '+expt.scans[0].id+': '+expt.scans[0].command,
                        position ='below', celltype='markdown', is_print=True)
            
        Utils.Delete_current_cell()

        
    def on_button_refresh_clicked(b):
        """
        Re-execute the cell to update it.
        """
        
        Utils.Refresh_current_cell()
        
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
        
        Utils.Create_cell(code=
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
                    'expt.thetazfactor=GIXD.Calib_thetaz(calib_thetaz_data)',
                    position='below', celltype='code', is_print = True)
        
        Utils.Create_cell(code='## Calibration thetaz', position='below', celltype='markdown', is_print = True)
        

        Utils.Create_cell(code='FE.Action.Choose(expt)', position ='at_bottom', celltype='code', is_print=False)        

        Utils.Delete_current_cell()
        
    def on_button_form_clicked(b):
        """
        Create a form.
        """
        #clear_output(wait=False)
        Utils.Create_cell(code='FE.Action.Create_form()',position ='below', celltype='code', is_print=False)
        Utils.Create_cell(code='FE.Action.Choose(expt)', position ='at_bottom', celltype='code', is_print=False)        
         
        Utils.Delete_current_cell()
            
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
        
        Utils.Delete_current_cell()
        
        Utils.Create_cell(code='', position ='below', celltype='markdown', is_print=True, is_execute=False)
    
        Utils.Create_cell(code='FE.Action.Choose(expt)', position ='at_bottom', celltype='code', is_print=False)
        
    def on_button_script_clicked(b):
        """
        Insert a script as a markdown cell.
        """ 
        Print_script(expt)
        
    def on_button_image_clicked(b):
        """
        Insert an image in a markdown cell.
        """ 
        Insert_image(expt)        
        
    
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
    
    # Click to insert an image
    button_image = widgets.Button(description="Insert image")
    button_image.on_click(on_button_image_clicked)

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

    buttons2 = widgets.HBox([button_markdown, button_script, button_wm, button_commands, button_image])
    display(buttons2)
    

def Export_nb_to_pdf(nb_name):
    '''
    Export the notebook to pdf using a command line through the OS.

    Parameters
    ----------
    nb_name : str
        full name of the notebook. Ex: 'JupyLabBook.ipynb'

    Returns
    -------
    bool
        export_done, True if the export suceeded without error/warning
    '''
    
    # Save the current state of the notebook (including the widgets)
    Utils.Save_nb()
    
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


    
def Create_form():
    '''Display the widget for creation and printing of the experiment form.'''    
    
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

        ########################################
        # Prepare the title in several parts
        # Avoid having a line larger than the page

        # Max number of characters allowed by line
        max_length = 70

       
        title_split = title.value.split(' ')
        title_blocks = [] 

        j=0
        for k in range(0,len(title_split)):
            title_part = ''
            for i in range(j,len(title_split)): 
                if len(title_part)<max_length:
                    title_part += title_split[i]+' '
                    j=j+1
                else:
                    break
            title_blocks.append(title_part)        

        for title_block in title_blocks:
            if title_block != '':        
                txt.append('$\Large \\color{red}{\\bf %s}$'%('\ '.join(title_block.split(' '))))

        
        # Former (simpler version), but allows line longer than the page width
        #ptitle = ('\ '.join(title.value.split(' ')))
        #txt.append('$\Large \\color{red}{\\bf %s}$'%ptitle)
        ########################################
        
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
            Utils.Create_cell(code=elem, position ='below', celltype='markdown', is_print=True)
        
        Utils.Create_cell(code='# Experimental setup', position ='below', celltype='markdown', is_print=True)
        
        # Remove the widget when done
        Utils.Delete_current_cell()
        
    button = widgets.Button(description="Print form")
    out = widgets.Output()
    display(button,out)
    button.on_click(on_button_clicked)

    
def Print_wm(expt):
    '''
    Print the positions of the motors/sensors from extraction of the command wm in a log file.

    Parameters
    ----------
    expt : object
        object from the class Expt
    '''    

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
        with open(expt.logs_dir+log_file, encoding="utf8", errors='ignore') as f:
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
                Utils.Create_cell(code=sensor_str, position ='below', celltype='markdown', is_print=True)

            # Put title    
            Utils.Create_cell(code='### '+w_select_wm.value.split('; ')[1], position ='below', celltype='markdown', is_print=True)
            
            Utils.Delete_current_cell()
        
            Utils.Create_cell(code='FE.Action.Choose(expt)',
                        position ='at_bottom', celltype='code', is_print=False)  
            
        button = widgets.Button(description="Insert positions")
        display(button)
        button.on_click(on_button_print_clicked)


    button = widgets.Button(description="Select log")
    button.on_click(on_button_select_clicked)

    display(widgets.HBox([w_select_log, button]))  
    
def Print_script(expt):
    '''
    Collect info from a wm command in a log and add it to the report.

    Parameters
    ----------
    expt : object
        object from the class Experiment
    '''

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

            Utils.Create_cell(code='### '+ script_name, position ='above', celltype='markdown', is_print=True)
            Utils.Create_cell(code=code, position ='above', celltype='markdown', is_print=True)
            
            Utils.Delete_current_cell()
        
            Utils.Create_cell(code='FE.Action.Choose(expt)',
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

    
    
    
def Insert_image(expt):
    '''
    Insert an image in the notebook.

    Parameters
    ----------
    expt : object
        object from the class Expt
    '''  

    # Check if there is already a path for images directories
    # If not, use the recording directory
    try:
        path_to_img_default = expt.path_to_img
    except:
        path_to_img_default = expt.recording_dir

    # Widget to write path
    w_path_to_img = widgets.Text(
            value=path_to_img_default,
            description='Images directory:',
            layout=widgets.Layout(width='800px', height='40px'),
            style={'description_width': 'initial'})

    def on_button_validate_path_clicked(b):
        """
        Validate the path of the images folder. Open selection for the image.
        """

        if not os.path.exists(w_path_to_img.value):
            print(PN._RED+"Wrong folder name."+PN._RESET)
            print("")
            return

        # Pass the current value of the directory to the default one
        expt.path_to_img = w_path_to_img.value

        # Define the list of img files in the directory
        list_img_files = [file for file in sorted(os.listdir(expt.path_to_img))][::-1]

        if len(list_img_files) < 1:
            print(PN._RED+"There is no image in this folder."+PN._RESET)
            print("")
            return                

        # Widget to select image
        w_select_img = widgets.Dropdown(
             options=list_img_files,
             value=list_img_files[-1],
             layout=widgets.Layout(width='300px'),
             style={'description_width': 'Image:'})

        def on_button_insert_image_clicked(b):
            """
            Insert an image in a markdown cell.
            """ 

            # Get and insert the image
            path_to_img = w_path_to_img.value+w_select_img.value

            Utils.Create_cell(code='![]('+ path_to_img+')', position ='above', celltype='markdown', is_print=True)
            
            Utils.Delete_current_cell()
        
            Utils.Create_cell(code='FE.Action.Choose(expt)',
                        position ='at_bottom', celltype='code', is_print=False)  

        button_insert_image = widgets.Button(description="Insert image")
        button_insert_image.on_click(on_button_insert_image_clicked)

        display(widgets.HBox([w_select_img, button_insert_image]))

    button_validate_path = widgets.Button(description="Validate path")
    button_validate_path.on_click(on_button_validate_path_clicked)
    display(widgets.HBox([w_path_to_img, button_validate_path]))
    

    
def Print_commands(expt):
    '''
    Print a list of commands from a log file. Add the scan numbers at the end of each line.

    Parameters
    ----------
    expt : object
        object from the class Experiment
    '''


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
        
        # Convert numpy.ndarray into a simple array
        # (ipywidget does not like numpy array as options)
        
        rlog_lines = [str(rlog_line) for rlog_line in rlog_lines]
        
        w_select_commands = widgets.SelectMultiple(options=rlog_lines,rows=30,
                                                   layout=widgets.Layout(width='800px'),
                                                   value=[str(rlog_lines[-1])])

        display(w_select_commands)

        def on_button_insert_commands_clicked(b):
            """
            Display selected commands and scan numbers in a markdown cell.
            """ 

            code = '```python\n'+''.join(w_select_commands.value)+'```'

            Utils.Create_cell(code=code, position ='above', celltype='markdown', is_print=True)

            Utils.Delete_current_cell()

            Utils.Create_cell(code='FE.Action.Choose(expt)',
                        position ='at_bottom', celltype='code', is_print=False)  

        # Click to insert a script
        button_insert_commands = widgets.Button(description="Insert commands")
        button_insert_commands.on_click(on_button_insert_commands_clicked)

        display(button_insert_commands)

    button = widgets.Button(description="Select log")
    button.on_click(on_button_select_clicked)

    display(widgets.HBox([w_select_log, button]))  
 

def Convert_logs(expt):
    '''
    Convert all the logs in a human-readable version with date and commands.

    Parameters
    ----------
    expt : object
        object from the class Experiment
    '''

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
    '''
    Extract commands from a log file.

    Parameters
    ----------
    expt : object
        object from the class Experiment
        
    Returns
    -------
    array_like
        rlog_lines, an array containing all the formatted commands
    '''

    # Define the elements to be removed of the original log file
    remove_elem_0 = ['#', ']', '\n']
    remove_elem = [
                  'DevError', 'desc =', 'origin =', 'reason =', 'severity =',
                  'connection'
                   ]

    # List of possible dates
    date_list = ['# Mon', '# Tue', '# Wed', '# Thu', '# Fri', '# Sat', '# Sun']
    date = ''

    with open(pathToFile, 'r', encoding="utf8", errors='ignore') as f:
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


def Find_command_in_logs(scan, expt):
    '''
    Find the command corresponding to the scan in the logs.

    Parameters
    ----------
    scan : object
        object from the class Scan

    expt : object
        object from the class Experiment
    '''

    scan_found = False
    scan.command = 'No command found'
    
    for log_file in expt.list_logs_files:
        with open(expt.logs_dir+log_file, encoding="utf8", errors='ignore') as f:
            for line in f:
                if "#" not in line: temp = line
                if scan.id in line:
                    # Remove the line jump
                    scan.command = temp.replace('\n','')
                    scan_found = True
            f.close()

        if scan_found:
            break
 



def Define_scan_identifiers(scan, expt):
    '''
    Create a series of identifiers for the current scan.

    Parameters
    ----------
    scan : object
        object from the class Scan

    expt : object
        object from the class Experiment
    '''

    # For example:
    # scan.nxs = 'SIRIUS_2017_12_11_08042.nxs'
    # scan.path = '/Users/arnaudhemmerle/recording/SIRIUS_2017_12_11_08042.nxs'
    # scan.id = 'SIRIUS_2017_12_11_08042'
    # scan.number = 8042
    
    scan.path = expt.recording_dir+scan.nxs
    scan.id = scan.nxs[:-4]
    split_name = scan.nxs.split('.')[0].split('_')
    scan.number = int(scan.nxs.split('.')[0].split('_')[-1])

