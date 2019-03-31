

## dependencies: 
    # dash
    # dash_core_components
    # dash_html_components
    # Input, Output, State

# colormode == 'custom' is the only choice for heatmap
# colormode == 'custom' is disableed for line, area, link, ribbon, twistedribbon
# colormode == 'ideogram' is disabled for link, ribbon, twistedribbon

# colormode is disabled for connector

import sys
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_colorscales
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
from colorpicker_box import ColorPickerBox
import dash_daq as daq
from dashapp_functions import *
import io
import base64
import json

__author__ = 'Jin Cui'
__version__ = '1.0.0'
__date__ = 'August 27 2018'

if (sys.version_info[0]!=3):
    raise Exception('PCircos requires Python 2, your current Python version is {}.{}.{}'.
                    format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

UPLOADBOX_STYLES = {'display': 'block', 'position': 'relative', 'width': '95%', 'height': '40px',
                    'line-height': '40px', 'border-width': '1px', 'border-style': 'dashed',
                    'border-radius': '5px', 'text-align': 'center', 'marginTop': '1.2em', 
                    'marginBottom': '0.8em', 'overflow': 'auto'
                    }
FS_STYLES = {'marginTop': '0.6em', 'marginBottom': '0.4em'}

'''
external_css = ['pcircos.css', 'jquery.dataTables.min.css', 'jquery-ui.css', 'https://codepen.io/chriddyp/pen/bWLwgP.css']
external_js = [ "https://code.jquery.com/jquery-3.3.1.js",
                "pcircos.js",
                "https://code.jquery.com/ui/1.10.4/jquery-ui.js",
                "https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js",
                "https://code.jquery.com/ui/1.10.4/jquery-ui.js",
                ]
'''

external_css = [
        'https://codepen.io/chriddyp/pen/bWLwgP.css',
        {
            'href': 'https://use.fontawesome.com/releases/v5.7.2/css/all.css',
            'rel': 'stylesheet',
            'integrity': 'sha384-fnmOCqbTlWIlj8LyTjo7mOUStjsKC4pOpQbqyi7RrhN7udi9RwhKkMHpvLbHG9Sr',
            'crossorigin': 'anonymous'
        }]

app = dash.Dash(__name__, external_stylesheets=external_css)

app.config.supress_callback_exceptions = True


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

checkbox_options = [*map(lambda x: {'label': str(x), 'value': str(x)}, fig_instance.get_chr_info()['chr_label'])].copy()

checkbox_values = fig_instance.get_chr_info()['chr_label'].tolist()
#checkbox_values = fig_instance.get_chr_info()['chr_label'].tolist().copy()




app.layout = html.Div([
                html.Div([
                    html.Details([
                        html.Summary('Ideogram'),
                        html.Div([

                            html.Details([
                                html.Summary('Ideogram upload', className='summary-secondary'),
                                html.Div([
                                    dcc.Upload(
                                        id='ideogram-upload',
                                        children=html.Div([
                                            'Upload ideogram here'
                                            ]),
                                        style=UPLOADBOX_STYLES,
                                    ),
                                    html.P('Choose a field separator', style=FS_STYLES),
                                    dcc.RadioItems(
                                        id='ideogram-fs',
                                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                                        value='Tab',
                                        labelStyle={'dispaly': 'inline-block'}
                                    ),
                                    html.P('Select degree range'),
                                    html.Div([
                                        dcc.Input(
                                            id='ideogram-degreerange-min', 
                                            type='number', 
                                            value=12, 
                                            step=1, 
                                            min=0, 
                                            max=179, 
                                            style={'width': '45%'}
                                        ),
                                        html.P('-', style={'display': 'inline-block'}),
                                        dcc.Input(
                                            id='ideogram-degreerange-max', 
                                            type='number', 
                                            value=348, 
                                            step=1, 
                                            min=180, 
                                            max=360, 
                                            style={'width': '45%'} 
                                        ),
                                    ], style={'display': 'inline-block'}),

                                    html.P('Ideogram fill opacity'),
                                    dcc.Input(
                                        id='ideogram-opacity', 
                                        type='number', 
                                        value=0.6, 
                                        step=0.05, 
                                        min=0, 
                                        max=1,
                                        style={'width': '45%', 'marginBottom': '1.2em'}
                                    ),

                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Ideogram annotation', className='summary-secondary'),
                                html.Div([
                                    html.P('Enable chromosome annotations', className='booleanswitch'),
                                    daq.BooleanSwitch(
                                        id='ideogram-chrannotation-enabler',
                                        on=True,
                                        style={'display':'inline-block'}
                                    ),

                                    html.Div([
                                        html.P('Input annotation radius'),
                                        dcc.Input(
                                            id='ideogram-chrannotation-radius', 
                                            type='number', 
                                            value=1.25, 
                                            min=1, 
                                            max=2, 
                                            step=0.05,
                                            style={'width': '45%'}
                                        ),
                                        html.P('Select annotation font size (px)'),
                                        dcc.Input(
                                            id='ideogram-chrannotation-fontsize', 
                                            type='number', 
                                            value=15, 
                                            min=5, 
                                            max=40, 
                                            step=1
                                        ),
                                        html.P('Select annotation font type'),
                                        dcc.RadioItems(
                                            id='ideogram-chrannotation-fonttype',
                                            options=[{'label': i, 'value': i} for i in ['normal', 'bold', 'italic', 'bold+italic']],
                                            value='bold',
                                            labelStyle={'display': 'inline-block'}
                                        ),
                                        html.P('Select annotation font color'),
                                        ColorPickerBox(id='ideogram-chrannotation-fontcolor', value={'rgb': {'r': 0, 'g': 0, 'b': 0, 'a': 1}}),
                                        html.P('Input annotation angle offset'),
                                        dcc.Input(
                                            id='ideogram-chrannotation-angleoffset',
                                            type='number',
                                            value=-90,
                                            min=-90,
                                            max=90,
                                            step=90,
                                            style={'width': '45%'}
                                        ),
                                        html.P('Input annotation angle limit'),
                                        dcc.Input(
                                            id='ideogram-chrannotation-anglelimit',
                                            type='number', 
                                            value=360,
                                            min=180,
                                            max=360,
                                            step=180,
                                            style={'width': '45%'}
                                        ),

                                    ], id='ideogram-chrannotation-font'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Ideogram filter', className='summary-secondary'),
                                html.Div([
                                    html.P('Select chromosomes'),
                                    dcc.Checklist(id='chromosome_checklist',
                                                options=checkbox_options,
                                                values=checkbox_values,
                                                labelStyle={'display': 'inline-block'}
                                    ),
                                    html.Div([
                                        html.P('Enable custom color', className='booleanswitch'),
                                        daq.BooleanSwitch(
                                            id='ideogram-customcolor-enabler',
                                            on=False,
                                            style={'display':'inline-block'}
                                        ),
                                    ]),
                                    html.Div(id='ideogram-customcolorlist'),
                                ], className='indent')
                            ]),
                            
                            html.Details([
                                html.Summary('Ideogram ticks', className='summary-secondary'),
                                html.Div([
                                    html.P('Enable ideogram ticks', className='booleanswitch'),
                                    daq.BooleanSwitch(
                                        id='ideogram-ticks-enabler',
                                        on=True,
                                        style={'display': 'inline-block'}
                                    ),

                                    html.Div([
                                        html.P('Input majortick spacing'),
                                    
                                        dcc.Input(
                                            id='ideogram-majortick-spacing',
                                            type='number',
                                            value=20000000,
                                            step=1000000,
                                        ),
                                        html.P('Input minortick spacing'),
                                        dcc.Input(
                                            id='ideogram-minortick-spacing',
                                            type='number',
                                            value=5000000,
                                            step=250000,
                                        ),

                                        html.P('Select tick label format'),
                                        dcc.RadioItems( 
                                            id='ideogram-tick-format',
                                            options=[{'label': i, 'value': i} for i in ['Mb', 'Kb', 'Gb']],
                                            value='Mb',
                                            labelStyle={'display': 'inline-block'},
                                        ),
                                        ], id='ideogram-tick-controls'
                                    ),
                                ], className='indent')
                            ]),

                        ], className='indent'),
                    ]),
                    html.Details([
                        html.Summary('Backgrounds'),
                        
                        html.Div([

                            html.Details([
                                html.Summary('Cytoband', className='summary-secondary'),
                                html.Div([
                                    dcc.Upload(
                                        id='cytoband-upload',
                                        children=html.Div(['Upload cytoband here']),
                                        style=UPLOADBOX_STYLES,
                                    ),
                                    html.P('Choose a field separator', style=FS_STYLES),
                                    dcc.RadioItems(
                                        id='cytoband-fs',
                                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                                        value='Tab',
                                        labelStyle={'dispaly': 'inline-block'}
                                    ),
                                    html.P('Select opacity value'),

                                    dcc.Input(id='cytoband-opacity', type='number', value=0.9, min=0, max=1, step=0.05),

                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Ring', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of ring element(s)'),

                                    html.Div([
                                        dcc.Input(id='ring-number', value=0, min=0, max=200, type='number', style={'display': 'inline-block'}),
                                        daq.BooleanSwitch(
                                            id='ring-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    
                                    html.Div(id='ring-expand'),
                                  
                                    ], style={'paddingLeft': '0.2em', 'marginBottom': '1.2em'}
                                ),
                            ]),

                            html.Details([
                                html.Summary('Highlight', className='summary-secondary'),
                                html.Div([
                                    dcc.Upload(
                                        id='highlight-upload',
                                        children=html.Div(['Upload highlight here']),
                                        style=UPLOADBOX_STYLES,
                                    ),
                                    
                                    html.P('Choose a field separator'),
                                    dcc.RadioItems(
                                        id='highlight-fs',
                                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                                        value='Tab',
                                        labelStyle={'dispaly': 'inline-block'}
                                    ),
                                ], className='indent'),
                            ]),

                            html.Details([
                                html.Summary('Annotation', className='summary-secondary'),
                                html.Div([
                                    dcc.Upload(
                                        id='annotation-upload',
                                        children=html.Div(['Upload annotation here']),
                                        style=UPLOADBOX_STYLES,
                                    ),
                                    html.P('Choose a field separator'),
                                    dcc.RadioItems(
                                        id='annotation-fs',
                                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                                        value='Tab',
                                        labelStyle={'dispaly': 'inline-block'}
                                    ),
                                ], className='indent')
                            ]),

                        ], className='indent'),
  
                        ]),

                    html.Details([
                        html.Summary('Plots'),
                        html.Div([

                            html.Details([
                                html.Summary('Histogram', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of histogram(s)'),
                                    html.Div([
                                        dcc.Input(id='histogram-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='histogram-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='histogram-expand'),
                                ], className='indent'),
                            ]),

                            html.Details([
                                html.Summary('Scatter', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of scatter(s)'),
                                    html.Div([
                                        dcc.Input(id='scatter-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='scatter-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='scatter-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Line', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of line(s)'),
                                    html.Div([
                                        dcc.Input(id='line-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='line-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='line-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Area', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of area(s)'),
                                    html.Div([
                                        dcc.Input(id='area-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='area-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='area-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Tile', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of tile(s)'),
                                    html.Div([
                                        dcc.Input(id='tile-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='tile-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='tile-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Heatmap', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of heatmap(s)'),
                                    html.Div([
                                        dcc.Input(id='heatmap-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='heatmap-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='heatmap-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Connector', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of connector(s)'),
                                    html.Div([
                                        dcc.Input(id='connector-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='connector-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='connector-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Link', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of link(s)'),
                                    html.Div([
                                        dcc.Input(id='link-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='link-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='link-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Ribbon', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of ribbon(s)'),
                                    html.Div([
                                        dcc.Input(id='ribbon-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='ribbon-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='ribbon-expand'),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Twisted ribbon', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of twisted ribbon(s)'),
                                    html.Div([
                                        dcc.Input(id='twistedribbon-number', value=0, min=0, max=20, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='twistedribbon-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    html.Div(id='twistedribbon-expand'),
                                ], className='indent')
                            ]),





                         
                        ])
                    ]),

                html.Button(id='submit-button', children='Submit'),

                html.Div([
                    dcc.Store(id='ideogram-output'),
                    dcc.Store(id='cytoband-output'),
                    dcc.Store(id='ring-output'),
                    dcc.Store(id='highlight-output'),
                    dcc.Store(id='annotation-output'),

                    dcc.Store(id='histogram-output'),
                    dcc.Store(id='scatter-output'),
                    dcc.Store(id='line-output'),
                    dcc.Store(id='area-output'),
                    dcc.Store(id='tile-output'),
                    dcc.Store(id='heatmap-output'),
                    dcc.Store(id='connector-output'),
                    dcc.Store(id='link-output'),
                    dcc.Store(id='ribbon-output'),
                    dcc.Store(id='twistedribbon-output'),


                    dcc.Store(id='merge_all')

                ]),

                ], className='sidebar'),
                #style={'width': '20%', 'display': 'inline-block', 'float': 'left', 'paddingRight': '10px', 'fontFamily': 'Times New Roman Times Serif'}),
                
                html.Div([
                    dcc.Graph(id='PCircos_update')
                ], className='graph')
                # style={'width': '80%', 'float': 'right', 'display': 'inline-block' })
            ], className="total")
                


def dash_input():
    return [Input('submit-button', 'n_clicks')]



def dash_state():
    dash_states = [State('ideogram-degreerange-min', 'value'),
                   State('ideogram-degreerange-max', 'value'),
                   State('chromosome_checklist', 'values'),
                   State('ideogram-tick-format', 'value'),
                   ]
    return dash_states




@app.callback(
    Output('ideogram-chrannotation-font', 'style'),
    [
        Input('ideogram-chrannotation-enabler', 'on')
    ]
)
def enable_ideogramchrannotation(bool_value):
    if bool_value is True:
        return {'display': 'block', 'paddingLeft': '0.2em'}
    else:
        return {'display': 'none'}



## display or hide custom colors
@app.callback(
    Output('ideogram-customcolorlist', 'style'),
    [
        Input('ideogram-customcolor-enabler', 'on')
    ]
)
def enable_ideogramcustomcolor(bool_value):
    if bool_value is True:
        return {'display': 'block', 'marginTop': '0.8em'}
    else:
        return {'display': 'none'}

## enable selection of custom colors for each chromsoome
@app.callback(
    Output('ideogram-customcolorlist', 'children'),
    [
        Input('chromosome_checklist', 'values')
    ]
)

def ideogram_colorbox(chromosome_checklist):
    return expand_chromosome_color(chromosome_checklist)



@app.callback(
    Output('ideogram-tick-controls', 'style'),
    [
        Input('ideogram-ticks-enabler', 'on')
    ]
)
def ideogram_tick_controls(bool_value):
    if bool_value is True:
        return {'display': 'block', 'paddingLeft': '0.2em'}
    else:
        return {'display': 'none'}


## primary output into hidden divs

@app.callback(
    Output('ideogram-output', 'data'),
    [
        Input('ideogram-upload', 'contents'),
        Input('ideogram-fs', 'value'),
        Input('ideogram-opacity', 'value'),
        Input('ideogram-chrannotation-enabler', 'on'),
        Input('ideogram-chrannotation-radius', 'value'),
        Input('ideogram-chrannotation-fontsize', 'value'),
        Input('ideogram-chrannotation-fonttype', 'value'),
        Input('ideogram-chrannotation-fontcolor', 'value'),
        Input('ideogram-chrannotation-angleoffset', 'value'),
        Input('ideogram-chrannotation-anglelimit', 'value'),
        Input('ideogram-degreerange-min', 'value'),
        Input('ideogram-degreerange-max', 'value'),
        Input('chromosome_checklist', 'values'),
        Input('ideogram-ticks-enabler', 'on'),
        Input('ideogram-majortick-spacing', 'value'),
        Input('ideogram-minortick-spacing', 'value'),
        Input('ideogram-tick-format', 'value'),

    ],
)

def ideogram_callback(contents, fs, opacity, chrannotation,
                      chrannotation_radius, chrannotation_fontsize, chrannotation_fonttype,
                      chrannotation_fontcolor, chrannotation_angleoffset, chrannotation_anglelimit,
                      degreerange_min, degreerange_max, chromosome_checklist, 
                      ticks_enabler, majortick_spacing, minortick_spacing, tick_format): 
    
    content_string = interp_contents(contents)
    sep = interp_fs(fs)

    degreerange = [degreerange_min, degreerange_max]

    ## need to add decoding function before reading content_string by pandas
    ideogram_dict = {
        'file': {
            'content_string': content_string, 'header': 'infer', 'sep': sep
        },
        'degreerange': degreerange,
        'showfillcolor': True,
        'layout': {'opacity': opacity},
        'chrannotation': {
            'show': chrannotation, 
            'radius': {'R': chrannotation_radius},
            'fonttype': chrannotation_fonttype,
            'textangle': {
                'angleoffset': chrannotation_angleoffset,
                'anglelimit': chrannotation_anglelimit
            },
            'layout': {
                'font': {
                    'size': chrannotation_fontsize,
                    'color': revert_rgb(chrannotation_fontcolor)
                }
            },                   
        },
        'majortick': {
            'show': ticks_enabler,
            'spacing': majortick_spacing
        },
        'minortick': {
            'show': ticks_enabler,
            'spacing': minortick_spacing
        },
        'ticklabel': {
            'show': ticks_enabler,
            'spacing': majortick_spacing,
            'textformat': tick_format,
            'textangle': {
                'angleoffset': chrannotation_angleoffset,
                'anglelimit': chrannotation_anglelimit
            }
        }
    }
    #print(ideogram_dict)
    return json.dumps(ideogram_dict)


@app.callback(
    Output('cytoband-output', 'data'),
    [
        Input('cytoband-upload', 'contents'),
        Input('cytoband-fs', 'value'),
        Input('cytoband-opacity', 'value'),
    ],
)

def cytoband_callback(contents, fs, opacity):

    content_string = interp_contents(contents)
    sep = interp_fs(fs)

    if content_string == None:
        cytoband_dict = {'show': False}
    else:
        cytoband_dict = {
            'file': {
                'content_string': content_string, 'header': 'infer', 'sep': sep
            },
            'layout': {'opacity': opacity}
        }
    print(cytoband_dict)
    return json.dumps(cytoband_dict)


@app.callback(
    Output('ring-number', 'disabled'),
    [
        Input('ring-number-lock', 'on')
    ]
)
def ring_number_lock(bool_value):
    return bool_value



@app.callback(
    Output('ring-expand', 'children'),
    [
        Input('ring-number', 'value')
    ],
)

def n_ring(number):
    return expand_ring(number)




@app.callback(
    Output('highlight-output', 'data'),
    [
        Input('highlight-upload', 'contents'),
        Input('highlight-fs', 'value')
    ],
)

def highlight_callback(contents, fs):
    content_string = interp_contents(contents)
    sep = interp_fs(fs)
    if content_string == None:
        highlight_dict = {'show': False}
    else:
        highlight_dict = {
            'file': {
                'content_string': content_string, 'header': 'infer', 'sep': sep
            }
        }
    #print(highlight_dict)
    return json.dumps(highlight_dict)
            

@app.callback(
    Output('annotation-output', 'data'),
    [
        Input('annotation-upload', 'contents'),
        Input('annotation-fs', 'value')
    ],
)

def annotation_callback(contents, fs):
    content_string = interp_contents(contents)
    sep = interp_fs(fs)
    if content_string == None:
        annotation_dict = {'show': False}
    else:
        annotation_dict = {
            'file': {
                'content_string': content_string, 'header': 'infer', 'sep': sep
            }
        }
    return json.dumps(annotation_dict)


@app.callback(
    Output('histogram-number', 'disabled'),
    [
        Input('histogram-number-lock', 'on')
    ]
)
def histogram_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('histogram-expand', 'children'),
    [
        Input('histogram-number', 'value')
    ]
)
def n_histogram(number):
    return expand_histogram(number)

  


@app.callback(
    Output('scatter-number', 'disabled'),
    [
        Input('scatter-number-lock', 'on')
    ]
)

def scatter_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('scatter-expand', 'children'),
    [
        Input('scatter-number', 'value')
    ]
)
def n_scatter(number):
    return expand_scatter(number)



@app.callback(
    Output('line-number', 'disabled'),
    [
        Input('line-number-lock', 'on')
    ]
)
def line_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('line-expand', 'children'),
    [
        Input('line-number', 'value')
    ]
)
def n_line(number):
    return expand_line(number)


@app.callback(
    Output('area-number', 'disabled'),
    [
        Input('area-number-lock', 'on')
    ]
)
def area_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('area-expand', 'children'),
    [
        Input('area-number', 'value')
    ]
)
def n_area(number):
    return expand_area(number)

@app.callback(
    Output('tile-number', 'disabled'),
    [
        Input('tile-number-lock', 'on')
    ]
)
def tile_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('tile-expand', 'children'),
    [
        Input('tile-number', 'value')
    ]
)
def n_tile(number):
    return expand_tile(number)



@app.callback(
    Output('heatmap-number', 'disabled'),
    [
        Input('heatmap-number-lock', 'on')
    ]
)
def heatmap_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('heatmap-expand', 'children'),
    [
        Input('heatmap-number', 'value')
    ]
)
def n_heatmap(number):
    return expand_heatmap(number)



@app.callback(
    Output('connector-number', 'disabled'),
    [
        Input('connector-number-lock', 'on')
    ]
)
def connector_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('connector-expand', 'children'),
    [
        Input('connector-number', 'value')
    ]
)
def n_connector(number):
    return expand_connector(number)


@app.callback(
    Output('link-number', 'disabled'),
    [
        Input('link-number-lock', 'on')
    ]
)
def link_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('link-expand', 'children'),
    [
        Input('link-number', 'value')
    ]
)
def n_link(number):
    return expand_link(number)


@app.callback(
    Output('ribbon-number', 'disabled'),
    [
        Input('ribbon-number-lock', 'on')
    ]
)
def ribbon_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('ribbon-expand', 'children'),
    [
        Input('ribbon-number', 'value')
    ]
)
def n_ribbon(number):
    return expand_ribbon(number)


@app.callback(
    Output('twistedribbon-number', 'disabled'),
    [
        Input('twistedribbon-number-lock', 'on')
    ]
)
def twistedribbon_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('twistedribbon-expand', 'children'),
    [
        Input('twistedribbon-number', 'value')
    ]
)
def n_twistedribbon(number):
    return expand_twistedribbon(number)


## ONGOING
## Difficult part, need to parse n_plots
'''
@app.callback(
    Output('histogram-colormode-choice_0', 'children'),
    [
        Input('histogram-colormode_0', 'value')
    ]
)
def histogram_colormode_choice(choice):
    if choice == 'Mono':
        return html.Div([
            html.P('Select color:'),
            ColorPickerBox(id='histogram_mono_0'.format(i), 
                            value='red', label='', 
                            style={'display': 'inline-block','margin': '0.2em'}
                            ),
            ], style={'display': 'inline-block'} )
    elif choice == 'Custom':
        return html.Div([
            html.P('Select which column contains the color'),
            dcc.Input(id='histogram_mono_0'.format(i),
                        type='number', 
                        value=12, 
                        step=1, 
                        min=0, 
                        max=179, 
                        style={'width': '45%'}
                        ),
            ], style={'display': 'inline-block'} )
'''

'''
for i in range(20):
    try:
        @app.callback(
            Output('histogram-colormode-choice_{}'.format(i), 'children'),
            [
                Input('histogram-colormode_{}'.format(i), 'value')
            ]
        )
        def histogram_colormode_choice(choice):
            if choice == 'Mono':
                return html.Div([
                    html.P('Select color:'),
                    ColorPickerBox(id='histogram_mono_{}'.format(i), 
                                   value='red', label='', 
                                   style={'display': 'inline-block','margin': '0.2em'}
                                  ),
                    ], style={'display': 'inline-block'} )
            elif choice == 'Custom':
                return html.Div([
                    html.P('Select which column contains the color'),
                    dcc.Input(id='histogram_mono_{}'.format(i),
                              type='number', 
                              value=12, 
                              step=1, 
                              min=0, 
                              max=179, 
                              style={'width': '45%'}
                             ),
                    ], style={'display': 'inline-block'} )

    except Exception:
        pass

'''
#for val in range(app.layout['histogram-number'].value):





@app.callback(
    Output('PCircos_update', 'figure'),
    dash_input(),
    dash_state()
    )

def update_output(_,
                  degreerange_min,
                  degreerange_max,
                  chromosome_checklist_value,
                  tick_format_value, 
                  
                  #*args
                  ):
    
    degreerange_value = [degreerange_min, degreerange_max]
 
    dash_dict = dict(degreerange=degreerange_value,
                     chromosome_checklist=chromosome_checklist_value,
                     #majortick_spacing=majortick_spacing_value,
                     #minortick_spacing=minortick_spacing_value,
                     tick_format=tick_format_value,
                     )
    ### SELECT CHROMOSOME
    #print ('figure instance is:')
    #print (Figure(input_json_path=sys.argv[1], dash_dict=dash_dict).fig())
    return Figure(input_json_path=sys.argv[1], dash_dict=dash_dict).fig()

if __name__ == '__main__':
    #print (app.layout)
    app.run_server(debug=True)


# python3 dashapp.py demo_data/demo_params.json