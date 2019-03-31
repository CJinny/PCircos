'''
This code helps parse the json config file and convert it into a dictionary
font_family reference:   https://www.w3schools.com/cssref/playit.asp?filename=playcss_font-family
'''
import sys
sys.path.append('../')
import pandas as pd
import json

def nested_eval(X):
    '''
    this function converts any stringed boolean (e.g. "True"), stringed nontype ("None"), stringed float ("1.5") stringed list (e.g. "[1,2,3]", ["1","1.2","1.3"]) to boolean, nontype and list variable,
    non-stringed units won't be affected. Any stringed variable nested inside list or dict will be reviewed.
    Note that ordinary string (e.g. "ABC" is not recognized and eval(string) will result in NameError, therefore they won't be converted in this function
    Reason for this function is that json doesn't store True, it needs "True"
    '''

    if isinstance(X, str):
        try:
            X = eval(X)
        except (NameError, SyntaxError):
            pass
    elif isinstance(X, dict):
        for key in X.keys():
            X[key] = nested_eval(X[key])
    elif isinstance(X, list):
        X = [*map(lambda x: nested_eval(x), X)]
    else:
        pass
    return X

def json_dict(json_config_path):
    with open(json_config_path, 'r') as f:
        json_read = f.read()
    f.close()

    # JSONDecodeError will be thrown if json file is invalid
    json_config = json.loads(json_read)
    config_dict = nested_eval(json_config)
    return config_dict


def nested_fill_dict(input_dict, default_dict):
    
    '''Assigning default values if key not found
    For nested dictionary or dict nested within list, search thoroughly and append default values if key not found!
    '''
    ## USAGE: nested_fill_dict(my_json_config, default_json_config)
    ## ONGOING

    for key in default_dict.keys():
        
        
        if key not in input_dict:
            if isinstance(default_dict[key], (str, int, float)):
                input_dict[key] = default_dict[key]
            
            elif key in ['xaxis', 'yaxis']:
                input_dict[key] = default_dict[key]
            else:
                continue
        
        elif isinstance(input_dict[key], dict):
            ## In cases where a plot is missing, for example no histogram is to be plotted, I want to make sure not to add default_dict['Category']['histogram] into the data
            if key != 'Category':
                nested_fill_dict(input_dict[key], default_dict[key])
            else: 
                nested_fill_dict(input_dict['Category'], default_dict['Category'])
                #for subkey in input_dict['Category'].keys():
                    #nested_fill_dict(input_dict['Category'][subkey], default_dict['Category'][subkey])

        elif isinstance(input_dict[key], list):
            # only dictionary list will be further investigated, stringed list and numeric list wont
            for element in input_dict[key]:
                if isinstance(element, dict):                    
                    nested_fill_dict(element, default_dict[key])
                else:
                    pass
        else:
            # if value is string, bool, int or anything else
            pass
            


def json2dict(input_json_path, dash_dict=None):
    '''Note that after assignment of default values, some additional information has to be provided: scatterplot default colors etc, I'll do this in layout and trace module '''
    
    input_json_dict = json_dict(input_json_path)

    try:
        defaultdict = json_dict("default_params.json")
    except FileNotFoundError:
        try:
            defaultdict = json_dict("config/default_params.json")
        except FileNotFoundError:
            try:
                defaultdict = json_dict("/home/jincui/anaconda3/bin/PCircos/config/default_params.json")
            except FileNotFoundError:
                print ('default_params json file not found')

    nested_fill_dict(input_json_dict, defaultdict)

    if dash_dict:

        input_json_dict['Category']['ideogram']['ideogram'].update(dict(degreerange=dash_dict['degreerange']))
        #input_json_dict['Category']['ideogram']['majortick'].update(dict(spacing=dash_dict['majortick_spacing']))
        #input_json_dict['Category']['ideogram']['minortick'].update(dict(spacing=dash_dict['minortick_spacing']))
        input_json_dict['Category']['ideogram']['ticklabel'].update(dict(textformat=dash_dict['tick_format']))

        '''
        for key in input_json_dict['Category'].keys():
            if key in ['ideogram', 'annotation', 'highlight']:
                continue

            if not isinstance(input_json_dict['Category'][key], list):
                
                if key != 'cytoband':
                    
                    radius_tmp = dash_dict['radius_args'].pop(0)
                    input_json_dict['Category'][key]['radius'].update(dict(R0=radius_tmp[0], R1=radius_tmp[1]))

                    if key not in ['tile', 'connector', 'link', 'heatmap', 'scatter', 'line']:
                        input_json_dict['Category'][key]['layout']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key]['layout']['fillcolor'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue 

                    elif key in ['tile', 'connector', 'link']:
                        input_json_dict['Category'][key]['layout']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key]['layout']['line']['color'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue 

                    elif key == 'scatter':
                        input_json_dict['Category'][key]['trace']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key]['trace']['marker']['color'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue 

                    elif key == 'line':
                        input_json_dict['Category'][key]['trace']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key]['trace']['line']['color'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue 

                    elif key == 'heatmap':
                        
                        input_json_dict['Category'][key]['layout']['opacity'] = dash_dict['opacity_args'].pop(0)

                        if dash_dict['color_args'][0] in ['RdYlBu', 'Spectral', 'RdYlGn', 'PiYG', 'PuOr', 'PRGn', 'BrBG', 'RdBu', 'RdGy']:
                            input_json_dict['Category'][key]['palatte'].update(dict(ncolor=11, scale='div', palatte=dash_dict['color_args'].pop(0)))
                        elif dash_dict['color_args'][0] in ['Reds', 'YlOrRd', 'RdPu', 'YlOrBr', 'Greens', 'YlGnBu', 'GnBu', 'BuPu', 'Greys', 
                                                            'Oranges', 'OrRd', 'BuGn', 'PuBu', 'PuRd', 'Blues', 'PuBuGn', 'YlGn', 'Purples']:
                            input_json_dict['Category'][key]['palatte'].update(dict(ncolor=9, scale='seq', palatte=dash_dict['color_args'].pop(0)))
                        else:
                            dash_dict['color_args'].pop(0)
                                
                else:
                    input_json_dict['Category'][key]['layout']['opacity'] = dash_dict['opacity_args'].pop(0)

            else:
                assert key != 'cytoband'

                for i in range(len(input_json_dict['Category'][key])):
                    #radius_tmp = dash_dict['radius_args'].pop(0)
                    #input_json_dict['Category'][key][i]['radius'].update(dict(R0=radius_tmp[0], R1=radius_tmp[1]))

                    
                    if key not in ['tile', 'connector', 'link', 'heatmap', 'scatter', 'line']:
                        input_json_dict['Category'][key][i]['layout']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key][i]['layout']['fillcolor'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue 

                    elif key in ['tile', 'connector', 'link']:
                        input_json_dict['Category'][key][i]['layout']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key][i]['layout']['line']['color'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue

                    elif key == 'scatter':
                        input_json_dict['Category'][key][i]['trace']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key][i]['trace']['marker']['color'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue


                    elif key == 'line':
                        input_json_dict['Category'][key][i]['trace']['opacity'] = dash_dict['opacity_args'].pop(0)
                        if dash_dict['color_args'][0] and len(dash_dict['color_args'][0]) >= 3:
                            input_json_dict['Category'][key][i]['trace']['line']['color'] = dash_dict['color_args'].pop(0)
                        else:
                            dash_dict['color_args'].pop(0)
                            continue

                    elif key == 'heatmap':
                        
                        input_json_dict['Category'][key][i]['layout']['opacity'] = dash_dict['opacity_args'].pop(0)

                        if dash_dict['color_args'][0] in ['RdYlBu', 'Spectral', 'RdYlGn', 'PiYG', 'PuOr', 'PRGn', 'BrBG', 'RdBu', 'RdGy']:
                            input_json_dict['Category'][key][i]['palatte'].update(dict(ncolor=11, scale='div', palatte=dash_dict['color_args'].pop(0)))
                        elif dash_dict['color_args'][0] in ['Reds', 'YlOrRd', 'RdPu', 'YlOrBr', 'Greens', 'YlGnBu', 'GnBu', 'BuPu', 'Greys', 
                                                            'Oranges', 'OrRd', 'BuGn', 'PuBu', 'PuRd', 'Blues', 'PuBuGn', 'YlGn', 'Purples']:
                            input_json_dict['Category'][key][i]['palatte'].update(dict(ncolor=9, scale='seq', palatte=dash_dict['color_args'].pop(0)))
                        else:
                            dash_dict['color_args'].pop(0)
            '''

    return input_json_dict



### FOR DEBUGGING ONLY
'''
def default_dict():
    try:
        default_json_dict = json_dict("default_params.json")
    except FileNotFoundError:
        try:
            default_json_dict = json_dict("config/default_params.json")
        except FileNotFoundError:
            try:
                default_json_dict = json_dict("/home/jincui/anaconda3/bin/PCircos/config/default_params.json")
            except FileNotFoundError:
                print ('default_params json file not found')
    assert isinstance(default_json_dict, dict)
    return default_json_dict'''