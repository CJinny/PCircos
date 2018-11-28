'''
TO DO:

  add margin in demo

  dash:
    chr checklist, 
    plot color states is not working,
    give states to rangeslider, mark original value
    option for anglelimit (90, 180, 360), angleoffset(0, -90)

    data log2 transform

'''
import sys
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import plotly.graph_objs as go
import numpy as np
import pandas as pd
import maths
import colors
import colorlover as cl
from Complex import Complex
from fig import Figure
from time import time

from config import json_config, coord_config

__author__ = 'Jin Cui'
__version__ = '1.0.0'
__date__ = 'August 27 2018'

if (sys.version_info[0]!=3):
    raise Exception('PCircos requires Python 2, your current Python version is {}.{}.{}'.
                    format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))



external_css = ['pcircos.css', 'jquery.dataTables.min.css', 'jquery-ui.css']
external_js = [ "https://code.jquery.com/jquery-3.3.1.js",
                "https://code.jquery.com/ui/1.10.4/jquery-ui.js",
                "https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js",
                "https://code.jquery.com/ui/1.10.4/jquery-ui.js",
                "https://raw.githubusercontent.com/CJinnny/plotly_circos/master/pcircos.js"
                ]

app = dash.Dash(__name__)

'''
app = dash.Dash(__name__, external_stylesheets=external_css,
                          external_scripts=external_js  )
'''



fig_instance = Figure(input_json_path=sys.argv[1])
SUM = fig_instance.SUM
x_range = fig_instance.config_dict['General']['xaxis']['range']

try:
    assert sum(x_range) == 0
except AssertionError:
    print ('Please make sure xaxis plus and minus limits have the same absolute value')

y_range = fig_instance.config_dict['General']['yaxis']['range']
try:
    assert sum(y_range) == 0
except AssertionError:
    print ('Please make sure yaxis plus and minus limits have the same absolute value')

try:
    assert x_range == y_range
except AssertionError:
    print ('Warning, please make sure xaxis and yaxis range is the same')

rlimit = (x_range[1] + y_range[1])/2.0

# extract default value from user input json file

degreerange = fig_instance.config_dict['Category']['ideogram']['ideogram']['degreerange'].copy()

majortick_spacing = fig_instance.config_dict['Category']['ideogram']['majortick']['spacing']
minortick_spacing = fig_instance.config_dict['Category']['ideogram']['minortick']['spacing']

checkbox_options = [*map(lambda x: {'label': str(x), 'value': str(x)}, fig_instance.get_chr_info()['chr_label'])].copy()

checkbox_values = fig_instance.get_chr_info()['chr_label'].tolist()
#checkbox_values = fig_instance.get_chr_info()['chr_label'].tolist().copy()



id_dict = {}

'''
Notice that scatter and line plot uses trace whereas other plots use layout

e.g. {'histogram': {'hist_id_1': config_dict['Category']['histogram'][0]},
                    'hist_id_2': config_dict['Category']['histogram'][1]},
                    },
      'scatter': {}
      
      
      }

'''

for key in fig_instance.config_dict['Category'].keys():

    # ideogram and annotation is not customized
    if key in ['ideogram', 'annotation', 'highlight']:
        continue

    if not isinstance(fig_instance.config_dict['Category'][key], list):
        
        if 'id' in fig_instance.config_dict['Category'][key]:
            id_dict[key] = {fig_instance.config_dict['Category'][key]['id']: fig_instance.config_dict['Category'][key]}
        else:
            raise ValueError('id is missing for {}'.format(key))

    else:
        id_dict[key] = {fig_instance.config_dict['Category'][key][0]['id']: fig_instance.config_dict['Category'][key][0]}
        
        for i in range(1, len(fig_instance.config_dict['Category'][key])):
            
            if 'id' in fig_instance.config_dict['Category'][key][i]:
                id_dict[key].update({fig_instance.config_dict['Category'][key][i]['id']: fig_instance.config_dict['Category'][key][i]})
            else:
                raise ValueError('id is missing for {} [{}]'.format(key, i))




id_list = []
for key in id_dict.keys():
    for subkey in id_dict[key].keys():
        id_list.append(subkey)
if not len(set(id_list)) == len(id_list):
    raise ValueError('Please make sure to input unique id value for each plot')





def plot_radius():
    # html.Div([])
    res = []
    for key in id_dict.keys():
        if key == 'cytoband':
            continue
        for subkey in id_dict[key].keys():
            res.extend([html.P(subkey.capitalize()),
                        dcc.RangeSlider(id='radius of {}'.format(subkey),
                                        min=0,
                                        max=rlimit,
                                        step=0.01,
                                        value=[id_dict[key][subkey]['radius']['R0'], id_dict[key][subkey]['radius']['R1']]
                                        )
            ])


    return html.Details([
                    html.Summary('Select plot radius range'),
                    html.Div(res)
                    ]) 

def plot_color():
    res = []
    for key in id_dict.keys():
        if key == 'cytoband':
            continue
        for subkey in id_dict[key].keys():
            res.extend([html.P(subkey.capitalize()),
                        dcc.Input(id='color of {}'.format(subkey),
                                    type='text',
                                    value='',
                                ),
                        html.Button(id='button of {}'.format(subkey), children='Submit'),            
            ])
    return html.Details([
                    html.Summary('Input plot color'),
                    html.Div(res)
                    ])

def plot_opacity():
    res = []
    for key in id_dict.keys():
        if key in ['scatter', 'line']:
            for subkey in id_dict[key].keys():
                res.extend([html.P(subkey.capitalize()),
                        dcc.Slider(id='opacity of {}'.format(subkey),
                                    min=0,
                                    max=1,
                                    step=0.05,
                                    value=id_dict[key][subkey]['trace']['opacity']
                                    )                
                ])
        else:
            for subkey in id_dict[key].keys():
                res.extend([html.P(subkey.capitalize()),
                        dcc.Slider(id='opacity of {}'.format(subkey),
                                    min=0,
                                    max=1,
                                    step=0.05,
                                    value=id_dict[key][subkey]['layout']['opacity']
                                    )                
                ])
    return html.Details([
                    html.Summary('Select plot opacity'),
                    html.Div(res)
                    ])

app.layout = html.Div([
                html.Div([
                html.Details([
                    html.Summary('Ideogram'),
                    html.Div([
                        html.P('Degree range'),
                        dcc.RangeSlider(
                            id='degreerange',
                            min=0,
                            max=360,
                            step=10,
                            value=degreerange
                        ),
                        html.P('Select chromosomes'),
                        dcc.Checklist(
                            id='chromosome_checklist',
                            options=checkbox_options,
                            values=checkbox_values,
                            labelStyle={'display': 'inline-block'}
                        ),
                        html.P('Select majortick spacing'),
                        dcc.Slider(
                            id='majortick_spacing',
                            min=min(1000*(SUM//2000000), majortick_spacing),
                            max=max(1000*(SUM//1000000), majortick_spacing),
                            step=0.1*(max(1000*(SUM//1000000), majortick_spacing) - min(1000*(SUM//2001000), majortick_spacing)),
                            value=majortick_spacing
                        ),
                        html.P('Select minortick spacing'),
                        dcc.Slider(
                            id='minortick_spacing',
                            min=min(1000*(SUM//8000000), minortick_spacing),
                            max=max(1000*(SUM//4000000), minortick_spacing),
                            step=0.1*(max(1000*(SUM//1000000), minortick_spacing) - min(1000*(SUM//2001000), minortick_spacing)),
                            value=minortick_spacing
                        ),
                        html.P('Select tick label format'),
                        dcc.RadioItems(
                            id='tick_format',
                            options=[{'label': i, 'value': i} for i in ['Mb', 'Kb', 'Gb']],
                            value='Mb',
                            labelStyle={'display': 'inline-block'}
                        ),
                    ])
                ]),
                plot_radius(),
                plot_color(),
                plot_opacity(),
                ], className='sidebar'),
                #style={'width': '20%', 'display': 'inline-block', 'float': 'left', 'paddingRight': '10px', 'fontFamily': 'Times New Roman Times Serif'}),
                
                html.Div([
                    dcc.Graph(id='PCircos_update')
                ], className='graph')
                # style={'width': '80%', 'float': 'right', 'display': 'inline-block' })
            ], className="total")
                


def dash_input():

    dash_inputs = [Input('degreerange', 'value'),
                   Input('chromosome_checklist', 'values'),
                   Input('majortick_spacing', 'value'),
                   Input('minortick_spacing', 'value'),
                   Input('tick_format', 'value'),
                   ]

    for key in id_dict.keys():
        if key != 'cytoband':
            for subkey in id_dict[key].keys():
                dash_inputs.append(Input('radius of {}'.format(subkey), 'value'))
    
    for key in id_dict.keys():
        if key != 'cytoband':
            for subkey in id_dict[key].keys():
                dash_inputs.append(Input('button of {}'.format(subkey), 'value'))
    
    for key in id_dict.keys():
        for subkey in id_dict[key].keys():
            dash_inputs.append(Input('opacity of {}'.format(subkey), 'value'))

    return dash_inputs


def dash_state():
    dash_states = []
    for key in id_dict.keys():
        if key != 'cytoband':
            for subkey in id_dict[key].keys():
                dash_states.append(State('color of {}'.format(subkey), 'value'))
    return dash_states


@app.callback(
    Output('PCircos_update', 'figure'),
    dash_input(),
    dash_state()
    )

def update_output(degreerange_value,
                  chromosome_checklist_value,
                  majortick_spacing_value,
                  minortick_spacing_value, 
                  tick_format_value, 
                  *args):
    n = 0
    for key in id_dict.keys():
        if key == 'cytoband':
            continue
        for subkey in id_dict[key].keys():
            n += 1

    radius_args = args[0:n]
    color_args = args[n:2*n]
    opacity_args = args[2*n:]



    dash_dict = dict(degreerange=degreerange_value,
                     chromosome_checklist=chromosome_checklist_value,
                     majortick_spacing=majortick_spacing_value,
                     minortick_spacing=minortick_spacing_value,
                     tick_format=tick_format_value,
                     radius_args=list(radius_args),
                     color_args=list(color_args),
                     opacity_args=list(opacity_args))
    ### SELECT CHROMOSOME
    #print ('figure instance is:')
    #print (Figure(input_json_path=sys.argv[1], dash_dict=dash_dict).fig())
    return Figure(input_json_path=sys.argv[1], dash_dict=dash_dict).fig()

if __name__ == '__main__':
    #print (app.layout)
    app.run_server()