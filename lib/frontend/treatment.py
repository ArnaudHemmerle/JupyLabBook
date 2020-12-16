from . import utils
from lib.extraction.common import PyNexus as PN

import ipywidgets as widgets
import ipysheet
import numpy as np
import matplotlib.pyplot as plt


"""Frontend library for all the widgets concerning the treatments in the notebook."""

def Choose(expt):
    '''
    Choose the treatment to be applied to the selected scan.

    Parameters
    ----------
    expt : object
        object from the class Experiment
    '''
    
    # Styling options for widgets
    style = {'description_width': 'initial'}
    tiny_layout = widgets.Layout(width='150px', height='40px')
    short_layout = widgets.Layout(width='200px', height='40px')
    medium_layout = widgets.Layout(width='250px', height='40px')
    large_layout = widgets.Layout(width='300px', height='40px')
    
    # Define the function called when clicking the button
    # DEFINE HERE A FUNCTION TO CREATE A CELL CALLING YOUR CUSTOM FUNCTION

    def on_button_GIXD_clicked(b):

        # plot
        try: value = expt.plot
        except: value = True
        w_plot = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Plot')
        
        # save
        try: value = expt.save
        except: value = True
        w_save = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Save')
        
        # GIXD_logx
        try: value = expt.GIXD_logx
        except: value = False
        w_GIXD_logx = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='log x')
        
        # GIXD_logy
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
        w_channel0 = widgets.FloatText(value=value, style=style, layout=short_layout, description='Vineyard (channel)')

        # thetazfactor
        try: value = expt.thetazfactor
        except: value = 0.000243
        w_thetazfactor = widgets.FloatText(value=value, style=style, layout=large_layout,
                                           description='thetazfactor (rad/channel)')
        
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
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print absorbers')        
        
        # GIXD_cmap
        try: value = expt.GIXD_cmap
        except: value = 'jet'
        w_GIXD_cmap = widgets.Select(value=value, style=style, rows=5, description='cmap',
                                options=['viridis', 'jet', 'Greys', 'cividis', 'hot'])        
            
        display(widgets.HBox([w_plot, w_save, w_show_data_stamps, w_verbose, w_show_absorbers,
                              w_GIXD_logx, w_GIXD_logy, w_GIXD_logz, w_GIXD_cmap]))        
        display(widgets.HBox([w_binsize, w_nblevels, w_moytocreate_str, w_channel0, w_computeqz]))
        display(widgets.HBox([w_wavelength, w_thetac, w_thetazfactor]))
                
        
        def on_button_plot_clicked(b):

            # Pass current values as default values
            expt.plot = w_plot.value
            expt.save = w_save.value
            expt.GIXD_logy = w_GIXD_logy.value
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
                             
                utils.Create_cell(code=
                            'x, y, '+
                            'daty, datyTop, datyBottom, datyFirstQuarter, '+
                            'mat, mat_binned, ch_binned, '+
                            'mean_pi, mean_area, mean_gamma =\\\n'+
                            'GIXD.Treat(nxs_filename=\''+scan.nxs+'\', '+
                            'recording_dir=expt.recording_dir, '+
                            'channel0='+str(expt.channel0)+', '+
                            'thetazfactor='+str(expt.thetazfactor)+', '+
                            'wavelength='+str(expt.wavelength)+', '+
                            'thetac='+str(expt.thetac)+', '+
                            'binsize='+str(expt.binsize)+', '+
                            'computeqz='+str(expt.computeqz)+', '+
                            'absorbers='+'\''+str(absorbers)+'\''+', '+
                            'logx='+str(expt.GIXD_logx)+', '+
                            'logy='+str(expt.GIXD_logy)+', '+
                            'logz='+str(expt.GIXD_logz)+', '+
                            'nblevels='+str(expt.nblevels)+', '+
                            'cmap=\''+str(expt.GIXD_cmap)+'\''+', '+
                            'working_dir=expt.working_dir, '+ 
                            'moytocreate='+str(list_moytocreate)+', '+
                            'show_data_stamps='+str(expt.show_data_stamps)+', '+
                            'plot='+str(expt.plot)+', '+
                            'save='+str(expt.save)+', '+
                            'verbose='+str(expt.verbose)+')',
                            position='below', celltype='code', is_print = True)

     
                if len(expt.scans)>1:
                    utils.Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)

            # Do as if the button next was clicked
            on_button_next_clicked(b)

        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        display(button_plot)
        

    def on_button_pilatus_clicked(b):

        # plot
        try: value = expt.plot
        except: value = True
        w_plot = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Plot')
        
        # save
        try: value = expt.save
        except: value = True
        w_save = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Save')
        
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
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print absorbers')
        
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

        display(widgets.HBox([w_plot, w_save, w_show_data_stamps, w_verbose, w_show_absorbers,
                              w_pilatus_logz, w_pilatus_cmap]))  
        display(widgets.HBox([w_xmin, w_xmax])) 
        display(widgets.HBox([w_ymin, w_ymax])) 

        def on_button_plot_clicked(b):

            # Pass current values as default values
            expt.plot = w_plot.value
            expt.save = w_save.value
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
                    
            
                utils.Create_cell(code='images_sum, integrated_x, integrated_y =\\\n'+
                            'PilatusSum.Treat(nxs_filename=\''+scan.nxs+'\' ,'+ 
                            'recording_dir=expt.recording_dir, '+
                            'xmin='+str(expt.xmin)+', '+
                            'xmax='+str(expt.xmax)+', '+
                            'ymin='+str(expt.ymin)+', '+
                            'ymax='+str(expt.ymax)+', '+                                      
                            'absorbers='+'\''+str(absorbers)+'\''+', '+
                            'logz='+str(expt.pilatus_logz)+', '+  
                            'cmap=\''+str(expt.pilatus_cmap)+'\''+', '+
                            'working_dir=expt.working_dir, '+     
                            'show_data_stamps='+str(expt.show_data_stamps)+', '+
                            'plot='+str(expt.plot)+', '+
                            'save='+str(expt.save)+', '+
                            'verbose='+str(expt.verbose)+')',
                            position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    utils.Create_cell(code='### '+scan.id+': '+scan.command,
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
        w_force_gamma_delta = widgets.Checkbox(value=value, style=style, layout=short_layout, description='Force gamma&delta')        

        # fgamma
        try: value = expt.fgamma
        except: value = 0.
        w_fgamma = widgets.FloatText(value=value, style=style, layout=short_layout, description='Forced gamma (deg)')        

        # fdelta
        try: value = expt.fdelta
        except: value = 0.
        w_fdelta = widgets.FloatText(value=value, style=style, layout=short_layout, description='Forced delta (deg)')         
        
        # pixel_size 
        try: value = expt.pixel_size
        except: value = 0.172
        w_pixel_size = widgets.FloatText(value=value, style=style, layout=tiny_layout, description='Pixel size (um)')

        # qxymin
        try: value = expt.qxymin
        except: value = 0.
        w_qxymin = widgets.FloatText(value=value, style=style, layout=short_layout, description='qxy min (nm-1)')   
 
        # qxymax
        try: value = expt.qxymax
        except: value = 1.
        w_qxymax = widgets.FloatText(value=value, style=style, layout=short_layout, description='qxy max (nm-1)')  
        
        # qzmin
        try: value = expt.qzmin
        except: value = 0.
        w_qzmin = widgets.FloatText(value=value, style=style, layout=short_layout, description='qz min (nm-1)')   
 
        # qzmax
        try: value = expt.qzmax
        except: value = 1.
        w_qzmax = widgets.FloatText(value=value, style=style, layout=short_layout, description='qz max (nm-1)')          
        
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
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout=tiny_layout, description='Print absorbers')
        
        # GIXS_cmap
        try: value = expt.GIXS_cmap
        except: value = 'viridis'
        w_GIXS_cmap = widgets.Select(value=value, style=style, rows=5, description='cmap',
                                options=['viridis', 'jet', 'Greys', 'cividis', 'hot'])
     
                 
        display(widgets.HBox([w_show_data_stamps, w_verbose, w_show_absorbers, w_GIXS_logz, w_GIXS_cmap, w_pixel_size]))        
        display(widgets.HBox([w_wavelength, w_distance, w_thetai, w_pixel_PONI_x, w_pixel_PONI_y]))
        display(widgets.HBox([w_force_gamma_delta, w_fgamma, w_fdelta]))
        display(widgets.HBox([w_qxymin, w_qxymax])) 
        display(widgets.HBox([w_qzmin, w_qzmax])) 
        
        def on_button_plot_clicked(b):
            
            # Pass current values as default values
            expt.GIXS_logz = w_GIXS_logz.value
            expt.wavelength = w_wavelength.value
            expt.thetai = w_thetai.value
            expt.distance = w_distance.value
            expt.pixel_PONI_x = w_pixel_PONI_x.value
            expt.pixel_PONI_y = w_pixel_PONI_y.value
            expt.pixel_size = w_pixel_size.value
            expt.qxymin = w_qxymin.value
            expt.qxymax = w_qxymax.value
            expt.qzmin = w_qzmin.value
            expt.qzmax = w_qzmax.value
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
                
                utils.Create_cell(code='CF.Extract_GIXS(nxs_filename=\''+scan.nxs+'\','+ 
                       'working_dir=expt.working_dir, recording_dir=expt.recording_dir,'+
                        'logz='+str(expt.GIXS_logz)+','+
                        'wavelength='+str(expt.wavelength)+','+
                        'thetai='+str(expt.thetai)+','+
                        'distance='+str(expt.distance)+','+
                        'pixel_PONI_x='+str(expt.pixel_PONI_x)+','+
                        'pixel_PONI_y='+str(expt.pixel_PONI_y)+','+
                        'pixel_size='+str(expt.pixel_size)+','+  
                        'qxymin='+str(expt.qxymin)+','+
                        'qxymax='+str(expt.qxymax)+','+
                        'qzmin='+str(expt.qzmin)+','+
                        'qzmax='+str(expt.qzmax)+','+
                        'show_data_stamps='+str(expt.show_data_stamps)+','+
                        'force_gamma_delta='+str(expt.force_gamma_delta)+','+
                        'fgamma='+str(expt.fgamma)+','+
                        'fdelta='+str(expt.fdelta)+','+
                        'verbose='+str(expt.verbose)+','+ 
                        'absorbers='+'\''+str(absorbers)+'\''+','+
                        'cmap=\''+str(expt.GIXS_cmap)+'\''+')',
                        position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    utils.Create_cell(code='### '+scan.id+': '+scan.command,
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
        w_show_absorbers = widgets.Checkbox(value=value, style=style, layout = short_layout, description='Print absorbers')
        
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

            # Fill the sheet with previous values or None entries
            try: to_fill = expt.arr_peaks_full
            except: to_fill = np.array([[None,None,None] for i in range(20)])
 
           # Determine the number of rows to dynamically add some empty rows
            nb_filled_rows = len([elem for elem in to_fill if (elem[0]!=None and elem[0]!='')])
            nb_empty_rows = len([elem for elem in to_fill if (elem[0]==None or elem[0]=='')])
            if nb_empty_rows<15:
                to_fill = np.append([elem for elem in to_fill if (elem[0]!=None and elem[0]!='')],
                                    np.array([[None,None,None] for i in range(15)]), axis = 0)
                
            sheet = ipysheet.easy.sheet(columns=3, rows=len(to_fill) ,column_headers = ['Name','Position','Use?(y/n)'])
            

  
            # ipysheet does not work correctly with no entries
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

                # Convert elems_str into a list
                list_elems = [int(w_elems_str.value.split(',')[i]) for i in range(len(w_elems_str.value.split(',')))]
                
                
                for scan in expt.scans:                                                              
                
                    print("Peaks on scan %s"%scan.nxs)
                    # Extract and plot the XRF with the peaks when validate is clicked
                    CF.Extract_XRF(nxs_filename=scan.nxs,
                                   working_dir=expt.working_dir,
                                   recording_dir=expt.recording_dir,
                                   logz=w_XRF_logz.value,
                                   list_elems=list_elems,
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
                                            
                utils.Create_cell(code='CF.Extract_XRF('+
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
                    utils.Create_cell(code='### '+scan.id+': '+scan.command,
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

                utils.Create_cell(code='CF.Plot_isotherm(nxs_filename=\''+scan.nxs+'\','+
                            'working_dir=expt.working_dir,recording_dir=expt.recording_dir,'+
                            'show_data_stamps='+str(w_show_data_stamps.value)+
                            ', verbose='+str(w_verbose.value)+', '+
                            'fast='+str(w_fastextract.value)+')',
                            position='below', celltype='code', is_print = True)

                if len(expt.scans)>1:
                    utils.Create_cell(code='### '+scan.id+': '+scan.command,
                            position ='below', celltype='markdown', is_print=True)

            # Do as if the button next was clicked
            on_button_next_clicked(b)                    
                    
                    
        button_plot = widgets.Button(description="Plot")
        button_plot.on_click(on_button_plot_clicked)
        
        display(button_plot)                    
                    
    # Actions relevant for single scan analysis only
    def on_button_1D_clicked(b):            
        scan = expt.scans[0]
        utils.Create_cell(code='CF.Plot_1D(nxs_filename=\''+scan.nxs+'\','+
                    'working_dir=expt.working_dir,recording_dir=expt.recording_dir,'+
                    'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)
        
        # Do as if the button next was clicked
        on_button_next_clicked(b) 
        
    def on_button_fit_erf_clicked(b):
        scan = expt.scans[0]
        utils.Create_cell(code='CF.GaussianRepartition_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)  

        # Do as if the button next was clicked
        on_button_next_clicked(b)        
        
    def on_button_fit_gau_clicked(b):
        scan = expt.scans[0]
        utils.Create_cell(code='CF.Gaussian_fit(nxs_filename=\''+scan.nxs+'\', recording_dir = expt.recording_dir,'+
                        'xLabel=\''+scan.xLabel+'\', yLabel=\''+scan.yLabel+'\')',
                    position='below', celltype='code', is_print = True)   
        
        # Do as if the button next was clicked
        on_button_next_clicked(b)  
        
    def on_button_vineyard_clicked(b):
        scan = expt.scans[0]
        utils.Create_cell(code='expt.channel0 = GIXD.Extract_channel0(nxs_filename=\''+scan.nxs+'\','+
                    'recording_dir=expt.recording_dir, binsize=10, logx=False, '+
                    'logy=False, logz=False, nblevels=50, cmap=\'jet\', '+
                    'show_data_stamps = True, verbose=True)',
                    position='below', celltype='code', is_print = True)
 
    
        # Do as if the button next was clicked
        on_button_next_clicked(b)  

    # Next action   
    def on_button_next_clicked(b):
        #clear_output(wait=False)
        
        utils.Delete_current_cell()
        
        utils.Create_cell(code='FE.action.Choose(expt)',
                    position ='at_bottom', celltype='code', is_print=False)        
        
    def on_button_markdown_clicked(b):
        """
        Insert a markdown cell below the current cell.
        """ 
        utils.Delete_current_cell()
        
        utils.Create_cell(code='', position ='below', celltype='markdown', is_print=True, is_execute=False)
    
        utils.Create_cell(code='FE.treatment.Choose(expt)', position ='at_bottom', celltype='code', is_print=False)
       
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
    
    
    

def Find_absorbers_in_logs(scan, expt):
    '''
    Find absorbers in the logs.

    Parameters
    ----------
    scan : object
        object from the class Scan
    expt : object
        object from the class Experiment
        
    Returns
    -------
    str
        absorbers       
    '''
    scan_found = False
    absorbers = 'No absorber found' 
    
    for log_file in expt.list_logs_files:
        with open(expt.logs_dir+log_file, encoding="utf8", errors='ignore') as f:
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

def Set_interactive_1D(scan):
    '''
    Extract the sensors from the nxs file and set an interactive 1D plot.

    Parameters
    ----------
    scan : object
        object from the class Scan
    '''    
    
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
    

