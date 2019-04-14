import numpy as np
import pandas as pd
import colorlover as cl
import sys
sys.path.append('../')
import colors
import io
import base64
'''
This code converts reads input file and output a dictionary which includes the following:
     theta_values for various circular plots,
     plot_color(convert color into rgb values if needed), 
     plot_opacity, 
     plot_width, 
     plot_line_type,
     offset_theta (annotation),
'''



def chr_info(input_file_path, 
             sep='\t', header='infer', 
             custom_label=False,
             custom_spacing=False,
             custom_color=False,
             dash_dict=None
             ):
   ### ONGOING
    try:
        chr_info_pd = pd.read_csv(input_file_path, sep=sep, header=header, engine='python').fillna(0)
    except Exception:
        try:
            chr_info_pd = pd.read_csv(input_file_path, sep=sep)
        except Exception:
            try:
                chr_info_pd = pd.read_csv(input_file_path, sep='\t')
            except Exception:
                try:
                    chr_info_pd = pd.read_csv(input_file_path, sep=' ')
                except Exception:
                    try:
                        chr_info_pd = pd.read_csv(input_file_path, sep=',')
                    except Exception:
                        print('unable to read ideogram, printing input file path below')
                        print(input_file_path)
        

    if dash_dict is not None:
        for i in range(len(chr_info_pd)):
            if custom_label is False:
                
                chr_info_pd = chr_info_pd[chr_info_pd.iloc[:,0].isin(dash_dict['chromosome_checklist'])]

            else:
                chr_info_pd = chr_info_pd[chr_info_pd.iloc[:,2].isin(dash_dict['chromosome_checklist'])]


    chr_info=np.array(chr_info_pd.iloc[:])
    # responsive to dashapp checklist
    # print(chr_info)


    assert len(chr_info) > 0
    assert chr_info.shape[1] >= 2
    config_dict = {}

    chr_name = chr_info[:,0]
    chr_size=chr_info[:,1]
    default_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrx','chry','chr23','chr24','chr0','chrM','chrUn','chrNA']


    if custom_label:
        try: 
            chr_label = chr_info[:,2]
        except IndexError:
            print ('3rd column is missing for chromosome custom label, using chromosome names as label')
            chr_label = chr_name
    else:
        chr_label = chr_name
    
    if custom_spacing:
        try:
            chr_spacing = chr_info[:,4]
        except IndexError:
            print ('5th column is missing for inter-chromosome custom spacing, using default spacing values')
            chr_spacing = np.ones(chr_size.shape)*(0.05*sum(chr_size)/len(chr_size))
    else:
        chr_spacing = np.ones(chr_size.shape)*(0.05*sum(chr_size)/len(chr_size))
    
        
    if custom_color:
        try:
            chr_color = colors.to_rgb(chr_info[:,3])
        except IndexError:
            print ('4th column is missing for custom chromosome color, using default method to assign chromosome colors')
            if len(chr_info) <= 30:
                chr_color = colors.to_rgb(default_list)[:len(chr_info)]
            else: 
                assert len(chr_info) <= 710
                chr_color = colors.random_rgb(len(chr_info))
    else: 
        if len(chr_info) <= 30:
                chr_color = colors.to_rgb(default_list)[:len(chr_info)]
        else: 
            assert len(chr_info) <= 710
            chr_color = colors.random_rgb(len(chr_info))
    
    config_dict['chr_color_dict'] = dict(zip(chr_label, chr_color))
    config_dict['chr_fillcolor'] = chr_color

    # I assume the start is half of the first chr_spacing

    if len(chr_info) == 1:
        start = chr_spacing[0]/2
        ideogram_bin = [[start, chr_size[0] + start]]
        # it would be like: [[10,110]] for example
    else:
        ideogram_bin=[]
        start = chr_spacing[0]/2
        for i in range(len(chr_info)-1):
            end = start + chr_size[i]
            ideogram_bin.append([start,end])
            start = end + chr_spacing[i+1]
        end=start + chr_size[-1]
        ideogram_bin.append([start,end])

    config_dict['chr_name'] = chr_name
    config_dict['chr_size'] = chr_size
    config_dict['chr_label'] = chr_label
    config_dict['chr_label_dict'] = dict(zip(chr_name, chr_label))
    
    config_dict['chr_spacing'] = chr_spacing
    config_dict['ideogram_bin'] =ideogram_bin

    return config_dict


def ideogram_coord_config(chr_info, npoints=1000,
                          show_major_tick=False, major_tick_spacing=40000000,
                          show_minor_tick=False, minor_tick_spacing=20000000,
                          show_tick_label=False, tick_label_spacing=40000000, 
                         ):
    '''This function helps measuring accumulative fine coordinates for ideogram, major, minor ticks and tick labels, default is only ideogram
    '''
    # ideogram_accum_coord_list is used for ideogram patch, so a reversed array is needed, 2X
    config_dict={}

    config_dict['chr_label'] = chr_info['chr_label']

    SUM = sum(chr_info['chr_size']) + sum(chr_info['chr_spacing'])

    config_dict['SUM'] = SUM

    step = SUM/npoints
    ideogram_bin = chr_info['ideogram_bin']
    ideogram_accum_coord_list = []
    for i in range(len(ideogram_bin)):
        ideogram_accum_coord_list.append(np.concatenate((np.arange(ideogram_bin[i][0], ideogram_bin[i][1]+1, step), np.arange(ideogram_bin[i][0], ideogram_bin[i][1]+1, step)[::-1])))
    config_dict['ideogram_accum_coord_list'] = ideogram_accum_coord_list

    if show_major_tick:
        # 0,0,1,1,2,2 ...
        major_tick_accum_coord_list = []
        for i in range(len(ideogram_bin)):
            major_tick_accum_coord_list.append(np.repeat(np.arange(ideogram_bin[i][0], ideogram_bin[i][1]+1, major_tick_spacing), 2))
        config_dict['major_tick_accum_coord_list'] = major_tick_accum_coord_list

    if show_minor_tick:
        minor_tick_accum_coord_list=[]
        for i in range(len(ideogram_bin)):
            minor_tick_accum_coord_list.append(np.repeat(np.arange(ideogram_bin[i][0], ideogram_bin[i][1]+1, minor_tick_spacing), 2))
        config_dict['minor_tick_accum_coord_list'] = minor_tick_accum_coord_list

    if show_tick_label:
        tick_label_accum_coord_list = []
        tick_label_non_accum_list = []
        for i in range(len(ideogram_bin)):
            tick_label_accum_coord_list.append(np.arange(ideogram_bin[i][0], ideogram_bin[i][1]+1, tick_label_spacing))
            tick_label_non_accum_list.append(np.arange(0, chr_info['chr_size'][i] + 1, tick_label_spacing))
        
        config_dict['tick_label_non_accum_list'] = tick_label_non_accum_list
        config_dict['tick_label_accum_coord_list'] = tick_label_accum_coord_list

    return config_dict


def read_data(input_file_path, category, chr_info, sep='\t', header='infer'):
    
    pd_data = pd.read_csv(input_file_path, sep=sep, header=header, engine='python').fillna('NA')
    

    ####ONGOING 
    pd_data = pd_data[pd_data.iloc[:,0].isin(chr_info['chr_name'])]
    for i in range(len(pd_data)):
        pd_data.iloc[i,0] = chr_info['chr_label_dict'][pd_data.iloc[i,0]]

    if category in ['link', 'ribbon', 'twistedribbon']:
        pd_data = pd_data[pd_data.iloc[:,3].isin(chr_info['chr_name'])]
        for i in range(len(pd_data)):
            pd_data.iloc[i,3] = chr_info['chr_label_dict'][pd_data.iloc[i,3]]

    return np.array(pd_data.iloc[:])

    
def data_array(input_file_path, category, chr_info, sep='\t', header='infer', 
               colorcolumn=None, to_rgb=True, sortbycolor=False
              ):
    '''This function helps:
         calculating accumulative fine coordinates from the data
         convert chr_name to chr_label
         possibly convert color string to rgb colors (for cytoband, heatmap)
         possibly sort data by colorcolumn (for cytoband, heatmap, annotation)
         Note that you should not sort data for other plots'''


    '''notice that from ring data doesn't need to be processed'''
    # this function is used for copmlex_config module as well as the hovertext module

    # accum_coord_config contains: category, data_chr data_input_coord,data_accum_coord,
    ### The need for sorting is because plotly layout can only take in one color at a time, 
    # so I have to append different colors each time, if data is sorted by color, we would reduce the number of append and save time
    ### for cytoband and heatmap, the default is to sortbycolor, for other plots, unless color is specified by a column, the default is that they have the same color

    supported_categories = ['cytoband', 'histogram', 'line', 'area', 'scatter', 'tile', 'heatmap', 'link', 'ribbon', 'twistedribbon', 'connector', 'highlight', 'annotation']
    ## for ring data, we don't need to use this function
    if category not in supported_categories:
        raise ValueError('Please choose a supported category, please note that ideogram, ring data does not go through here')

    data_dict = {}
    input_pd = pd.read_csv(input_file_path, sep=sep, header=header, engine='python').fillna(0)
 
    chr_ideogram_bin = chr_info['ideogram_bin']

    # convert chr_name to chr_label, e.g. hs1 => chr1

    input_pd = input_pd[input_pd.iloc[:,0].isin(chr_info['chr_name'])]
    for i in range(len(input_pd)):
        input_pd.iloc[i,0] = chr_info['chr_label_dict'][input_pd.iloc[i,0]]

    if category in ['link', 'ribbon', 'twistedribbon']:
        input_pd = input_pd[input_pd.iloc[:,3].isin(chr_info['chr_name'])]
        for i in range(len(input_pd)):
            input_pd.iloc[i,3] = chr_info['chr_label_dict'][input_pd.iloc[i,3]]


    input_array = np.array(input_pd.iloc[:])
    # can't make this assertion in dash app since if all data is on the chromosome which we uncheck input_array will be non-existent
    # assert len(input_array) > 0

    data_chr = input_array[:,0]
    
    


    if category in ['link', 'ribbon', 'twistedribbon']:
        try:
            data_chr_1 = input_array[:,3]
        except IndexError:
            print ('error, printing input_array')
            print (input_array) 

        '''
        for i in range(len(input_array)):
            try:
                data_chr_1[i] = chr_label_dict[data_chr_1[i]]
            except KeyError:
                np.delete(input_array, i, 0)
        '''

        coord = input_array[:,1:3]
        coord_1 = input_array[:,4:6]

        for i in range(len(chr_ideogram_bin)): 
            for j in range(len(input_array)):
                if data_chr[j] == chr_info['chr_label'][i]:
                    coord[j,:] += chr_ideogram_bin[i][0]
                if data_chr_1[j] == chr_info['chr_label'][i]:
                    coord_1[j,:] += chr_ideogram_bin[i][0]

    elif category in ['cytoband', 'histogram', 'heatmap', 'highlight', 'tile', 'connector']:
        coord = input_array[:,1:3]
        for i in range(len(chr_ideogram_bin)): 
            for j in range(len(input_array)):
                if data_chr[j] == chr_info['chr_label'][i]:
                    coord[j,:] += chr_ideogram_bin[i][0]
            
    elif category in ['scatter', 'line', 'area', 'annotation']:
        pos = input_array[:,1]
        for i in range(len(chr_ideogram_bin)): 
            for j in range(len(input_array)):
                if data_chr[j] == chr_info['chr_label'][i]:
                    pos[j] += chr_ideogram_bin[i][0]

    if isinstance(colorcolumn, int):
        if to_rgb:
            input_array[:,colorcolumn] = colors.to_rgb(input_array[:,colorcolumn])

    if sortbycolor == True:
        if not isinstance(colorcolumn, int):
            data_dict['sortindex'] = list(range(len(input_array)))
        else:
            indices = input_array[:,colorcolumn].argsort()
            input_array = input_array[indices]

            # DEBUGGING:
            #print ('the real sortindex for {} is: \t'.format(category))
            #print (indices)

            data_dict['sortindex'] = indices
    else:
        data_dict['sortindex'] = list(range(len(input_array)))
    data_dict['data_array'] = input_array
    return data_dict
    