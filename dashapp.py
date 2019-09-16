

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
from dash_dict import *
from dash.exceptions import PreventUpdate


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
P_NELEMENT_STYLES = {'width': '10.5em'}
SWITCH_STYLES = {'display': 'inline-block', 'marginTop': '0.3em' }
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
server = app.server

#app.config.supress_callback_exceptions = True


#fig_instance = Figure(input_json_path=sys.argv[1])
#checkbox_options = [*map(lambda x: {'label': str(x), 'value': str(x)}, fig_instance.get_chr_info()['chr_label'])].copy()
#checkbox_values = fig_instance.get_chr_info()['chr_label'].tolist()




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
                                            value=0, 
                                            step=10, 
                                            min=0, 
                                            max=179, 
                                            style={'width': '45%'}
                                        ),
                                        html.P('-', style={'display': 'inline-block'}),
                                        dcc.Input(
                                            id='ideogram-degreerange-max', 
                                            type='number', 
                                            value=360, 
                                            step=10, 
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
                                    
                                    expand_ring(10),
                                  
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
                                        dcc.Input(id='histogram-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='histogram-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='histogram-expand'),
                                    html.Div([expand_histogram()]),

                                ], className='indent'),
                            ]),

                            html.Details([
                                html.Summary('Scatter', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of scatter(s)'),
                                    html.Div([
                                        dcc.Input(id='scatter-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='scatter-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='scatter-expand'),
                                    html.Div([expand_scatter()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Line', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of line(s)'),
                                    html.Div([
                                        dcc.Input(id='line-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='line-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='line-expand'),
                                    html.Div([expand_line()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Area', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of area(s)'),
                                    html.Div([
                                        dcc.Input(id='area-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='area-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='area-expand'),
                                    html.Div([expand_area()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Tile', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of tile(s)'),
                                    html.Div([
                                        dcc.Input(id='tile-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='tile-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='tile-expand'),
                                    html.Div([expand_tile()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Heatmap', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of heatmap(s)'),
                                    html.Div([
                                        dcc.Input(id='heatmap-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='heatmap-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='heatmap-expand'),
                                    html.Div([expand_heatmap()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Connector', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of connector(s)'),
                                    html.Div([
                                        dcc.Input(id='connector-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='connector-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='connector-expand'),
                                    html.Div([expand_connector()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Link', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of link(s)'),
                                    html.Div([
                                        dcc.Input(id='link-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='link-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='link-expand'),
                                    html.Div([expand_link()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Ribbon', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of ribbon(s)'),
                                    html.Div([
                                        dcc.Input(id='ribbon-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='ribbon-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='ribbon-expand'),
                                    html.Div([expand_ribbon()]),
                                ], className='indent')
                            ]),

                            html.Details([
                                html.Summary('Twisted ribbon', className='summary-secondary'),
                                html.Div([
                                    html.P('Input the number of twisted ribbon(s)'),
                                    html.Div([
                                        dcc.Input(id='twistedribbon-number', value=0, min=0, max=5, type='number', style={'display': 'inline-block', 'width': '45%'}),
                                        daq.BooleanSwitch(
                                            id='twistedribbon-number-lock',
                                            on=False,
                                            label='lock',
                                            labelPosition='bottom',
                                            style={'display': 'inline-block', 'float': 'right'}
                                        )
                                    ], style={'marginBottom': '1.2em'}),
                                    #html.Div(id='twistedribbon-expand'),
                                    html.Div([expand_twistedribbon()]),
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

                    dcc.Store(id='combined-output')

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
'''
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

@app.callback(
    Output('ideogram-customcolorlist', 'children'),
    [
        Input('chromosome-checklist', 'values')
    ]
)

def ideogram_colorbox(chromosome_checklist):
    return expand_chromosome_color(chromosome_checklist)
'''


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
        #Input('chromosome-checklist', 'values'),
        Input('ideogram-ticks-enabler', 'on'),
        Input('ideogram-majortick-spacing', 'value'),
        Input('ideogram-minortick-spacing', 'value'),
        Input('ideogram-tick-format', 'value'),

    ],
)

def ideogram_callback(contents, fs, opacity, chrannotation,
                      chrannotation_radius, chrannotation_fontsize, chrannotation_fonttype,
                      chrannotation_fontcolor, chrannotation_angleoffset, chrannotation_anglelimit,
                      degreerange_min, degreerange_max, 
                      #chromosome_checklist, 
                      ticks_enabler, majortick_spacing, minortick_spacing, tick_format
                     ): 

    degreerange = [degreerange_min, degreerange_max]

    try:
        content_string = interp_contents(contents)
        decoded = base64.b64decode(content_string)
        path = io.StringIO(decoded.decode('utf-8'))
        df = pd.read_csv(path, sep=interp_fs(fs), header='infer')
        print(df.head())
        print('ideogram upload successful')
        print(path)
    except Exception:
        print('unable to print ideogram dict')
        

    ideogram_dict = {
        'patch':{
            'file': {
                'path': content_string, 'sep': interp_fs(fs)
            },
            'radius': {'R0': 1.0, 'R1': 1.1},
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
            "customoptions": {
                "customlabel": "True",
                "customspacing": "False",
                "customcolor": "True"
            },
        },
        'majortick': {
            'show': ticks_enabler,
            'spacing': majortick_spacing,
            "radius": {"R0": 1.1, "R1": 1.12},
        },
        'minortick': {
            'show': ticks_enabler,
            'spacing': minortick_spacing,
            "radius": {"R0": 1.1, "R1": 1.108},
        },
        'ticklabel': {
            'show': ticks_enabler,
            'spacing': majortick_spacing,
            "radius": {"R": 1.16},
            'textformat': tick_format,
            'textangle': {
                'angleoffset': chrannotation_angleoffset,
                'anglelimit': chrannotation_anglelimit
            }
        }
    }
   

    return ideogram_dict


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
                'path': content_string, 'header': 'infer', 'sep': sep
            },
            'layout': {'opacity': opacity}
        }
    #print(cytoband_dict)
    return cytoband_dict

'''
@app.callback(
    Output('ring-number', 'disabled'),
    [
        Input('ring-number-lock', 'on')
    ]
)
def ring_number_lock(bool_value):
    return bool_value
'''

'''
@app.callback(
    Output('ring-expand', 'children'),
    [
        Input('ring-number', 'value')
    ],
)

 
#def n_ring(number):
    #return expand_ring(10)
'''


# ONGOING 
@app.callback(
    Output('ring-output', 'data'),
    ring_input_list(),
    
)

def store_ring(*args):
    res = []
    for i in range(len(args)//4):
        di = {
            'radius': {
                'R0': args[4*i],
                'R1': args[4*i+1]
            },
            'layout': {
                'opacity': args[4*i+2],
                'fillcolor': 'rgb({},{},{})'.format(args[4*i+3]['rgb']['r'], args[4*i+3]['rgb']['g'], args[4*i+3]['rgb']['b']),
                'layer': 'below',
                'line': {
                    'width': 0
                }
            }
        }
        res.append(di)
    #print('debugging store_ring 0')
    #print(args)
    #print('debugging store_ring res')
    #print(res)
    return res




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
                'path': content_string, 'header': 'infer', 'sep': sep
            }
        }
    #print(highlight_dict)
    return highlight_dict
            

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
                'path': content_string, 'header': 'infer', 'sep': sep,
            },
            'customcolor': False
        }
    return annotation_dict


## display or hide depending on the number of histogram user inputs
### I have to list all possible combinations since one of the limitations of dash is multiple outputs

### histogram list of callbacks
@app.callback(
    Output('histogram-idx_0', 'style'),
    [Input('histogram-number', 'value')]
)
def toggle_histogram_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('histogram-idx_1', 'style'),
    [Input('histogram-number', 'value')]
)
def toggle_histogram_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('histogram-idx_2', 'style'),
    [Input('histogram-number', 'value')]
)
def toggle_histogram_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('histogram-idx_3', 'style'),
    [Input('histogram-number', 'value')]
)
def toggle_histogram_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('histogram-idx_4', 'style'),
    [Input('histogram-number', 'value')]
)
def toggle_histogram_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('histogram-colormode-mono_0', 'style'),
    [Input('histogram-colormode_0', 'value')]
)

def histogram_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('histogram-colormode-custom_0', 'style'),
    [Input('histogram-colormode_0', 'value')]
)

def histogram_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('histogram-colormode-mono_1', 'style'),
    [Input('histogram-colormode_1', 'value')]
)

def histogram_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('histogram-colormode-custom_1', 'style'),
    [Input('histogram-colormode_1', 'value')]
)

def histogram_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('histogram-colormode-mono_2', 'style'),
    [Input('histogram-colormode_2', 'value')]
)

def histogram_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('histogram-colormode-custom_2', 'style'),
    [Input('histogram-colormode_2', 'value')]
)

def histogram_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('histogram-colormode-mono_3', 'style'),
    [Input('histogram-colormode_3', 'value')]
)

def histogram_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('histogram-colormode-custom_3', 'style'),
    [Input('histogram-colormode_3', 'value')]
)

def histogram_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('histogram-colormode-mono_4', 'style'),
    [Input('histogram-colormode_4', 'value')]
)

def histogram_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('histogram-colormode-custom_4', 'style'),
    [Input('histogram-colormode_4', 'value')]
)

def histogram_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('histogram-output', 'data'),
    histogram_input_list()
)

def store_histogram(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[9*i]), 'sep': interp_fs(args[9*i+1])
            },
            'hovertextformat': args[9*i+2],
            'radius': {'R0': args[9*i+3], 'R1': args[9*i+4]},
            'layout': {'opacity': args[9*i+5]}
        }
        if args[9*i+6] == 'By Chromosome':
            di['colorcolumn'] = 'ideogram'
        elif args[9*i+6] == 'Mono':
            #print('histogram args[9*i+7] is:')
            #print(args[9*i+7])
            try:
                di['layout']['fillcolor'] = args[9*i+7]['hex']
            except Exception:
                di['layout']['fillcolor'] = 'rgb(241,112,19)'
        elif args[9*i+6] == 'Custom':
            di['colorcolumn'] = args[9*i+8]
        res.append(di)
    #print(res)
    return res



#### scatter list of callbacks
@app.callback(
    Output('scatter-idx_0', 'style'),
    [Input('scatter-number', 'value')]
)
def toggle_scatter_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('scatter-idx_1', 'style'),
    [Input('scatter-number', 'value')]
)
def toggle_scatter_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('scatter-idx_2', 'style'),
    [Input('scatter-number', 'value')]
)
def toggle_scatter_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('scatter-idx_3', 'style'),
    [Input('scatter-number', 'value')]
)
def toggle_scatter_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('scatter-idx_4', 'style'),
    [Input('scatter-number', 'value')]
)
def toggle_scatter_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}




@app.callback(
    Output('scatter-colormode-mono_0', 'style'),
    [Input('scatter-colormode_0', 'value')]
)

def scatter_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('scatter-colormode-custom_0', 'style'),
    [Input('scatter-colormode_0', 'value')]
)

def scatter_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('scatter-colormode-mono_1', 'style'),
    [Input('scatter-colormode_1', 'value')]
)

def scatter_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('scatter-colormode-custom_1', 'style'),
    [Input('scatter-colormode_1', 'value')]
)

def scatter_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('scatter-colormode-mono_2', 'style'),
    [Input('scatter-colormode_2', 'value')]
)

def scatter_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('scatter-colormode-custom_2', 'style'),
    [Input('scatter-colormode_2', 'value')]
)

def scatter_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('scatter-colormode-mono_3', 'style'),
    [Input('scatter-colormode_3', 'value')]
)

def scatter_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('scatter-colormode-custom_3', 'style'),
    [Input('scatter-colormode_3', 'value')]
)

def scatter_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('scatter-colormode-mono_4', 'style'),
    [Input('scatter-colormode_4', 'value')]
)

def scatter_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('scatter-colormode-custom_4', 'style'),
    [Input('scatter-colormode_4', 'value')]
)

def scatter_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('scatter-output', 'data'),
    scatter_input_list()
)

def store_scatter(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[11*i]), 'sep': interp_fs(args[11*i+1])
            },
            'hovertextformat': args[11*i+2],
            'radius': {'R0': args[11*i+3], 'R1': args[11*i+4]},
            'trace': {'marker': {'size': args[11*i+6], 'opacity': args[11*i+5], 'symbol': args[11*i+7]}}
        }
        if args[11*i+8] == 'By Chromosome':
            di['colorcolumn'] = 'ideogram'
        elif args[11*i+8] == 'Mono':
            di['colorcolumn'] = 'None'
            print(args[11*i+9])
            di['trace']['marker']['color'] = args[11*i+9]['hex']
        elif args[11*i+8] == 'Custom':
            di['colorcolumn'] = args[11*i+10]
        res.append(di)
    #print('original dashapp store scatter res')
    #print(res)
    return res

###
#### line list of callbacks
@app.callback(
    Output('line-idx_0', 'style'),
    [Input('line-number', 'value')]
)
def toggle_line_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('line-idx_1', 'style'),
    [Input('line-number', 'value')]
)
def toggle_line_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('line-idx_2', 'style'),
    [Input('line-number', 'value')]
)
def toggle_line_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('line-idx_3', 'style'),
    [Input('line-number', 'value')]
)
def toggle_line_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('line-idx_4', 'style'),
    [Input('line-number', 'value')]
)
def toggle_line_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}



@app.callback(
    Output('line-colormode-mono_0', 'style'),
    [Input('line-colormode_0', 'value')]
)

def line_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('line-colormode-custom_0', 'style'),
    [Input('line-colormode_0', 'value')]
)

def line_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('line-colormode-mono_1', 'style'),
    [Input('line-colormode_1', 'value')]
)

def line_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('line-colormode-custom_1', 'style'),
    [Input('line-colormode_1', 'value')]
)

def line_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('line-colormode-mono_2', 'style'),
    [Input('line-colormode_2', 'value')]
)

def line_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('line-colormode-custom_2', 'style'),
    [Input('line-colormode_2', 'value')]
)

def line_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('line-colormode-mono_3', 'style'),
    [Input('line-colormode_3', 'value')]
)

def line_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('line-colormode-custom_3', 'style'),
    [Input('line-colormode_3', 'value')]
)

def line_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('line-colormode-mono_4', 'style'),
    [Input('line-colormode_4', 'value')]
)

def line_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('line-colormode-custom_4', 'style'),
    [Input('line-colormode_4', 'value')]
)

def line_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     


@app.callback(
    Output('line-output', 'data'),
    line_input_list()
)

def store_line(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[12*i]), 'sep': interp_fs(args[12*i+1])
            },
            'hovertextformat': args[12*i+2],
            'radius': {'R0': args[12*i+3], 'R1': args[12*i+4]},
            'trace': {'line': {'width': args[12*i+7], 'smoothing': args[12*i+8]}, 'opacity': args[12*i+5], 'marker': {'size': args[12*i+6]} }
        }
        if args[12*i+9] == 'By Chromosome':
            di['colorcolumn'] = 'ideogram'
        elif args[12*i+9] == 'Mono':
            di['colorcolumn'] = 'None'
            di['trace']['marker']['color'] = args[12*i+10]['hex']
            di['trace']['line']['color'] = args[12*i+10]['hex']
           
        elif args[12*i+9] == 'Custom':
            di['colorcolumn'] = args[12*i+11]
        res.append(di)
    #print(res)
    return res


#### area list of callbacks
@app.callback(
    Output('area-idx_0', 'style'),
    [Input('area-number', 'value')]
)
def toggle_area_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('area-idx_1', 'style'),
    [Input('area-number', 'value')]
)
def toggle_area_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('area-idx_2', 'style'),
    [Input('area-number', 'value')]
)
def toggle_area_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('area-idx_3', 'style'),
    [Input('area-number', 'value')]
)
def toggle_area_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('area-idx_4', 'style'),
    [Input('area-number', 'value')]
)
def toggle_area_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('area-colormode-mono_0', 'style'),
    [Input('area-colormode_0', 'value')]
)

def area_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('area-colormode-custom_0', 'style'),
    [Input('area-colormode_0', 'value')]
)

def area_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('area-colormode-mono_1', 'style'),
    [Input('area-colormode_1', 'value')]
)

def area_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('area-colormode-custom_1', 'style'),
    [Input('area-colormode_1', 'value')]
)

def area_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('area-colormode-mono_2', 'style'),
    [Input('area-colormode_2', 'value')]
)

def area_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('area-colormode-custom_2', 'style'),
    [Input('area-colormode_2', 'value')]
)

def area_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('area-colormode-mono_3', 'style'),
    [Input('area-colormode_3', 'value')]
)

def area_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('area-colormode-custom_3', 'style'),
    [Input('area-colormode_3', 'value')]
)

def area_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('area-colormode-mono_4', 'style'),
    [Input('area-colormode_4', 'value')]
)

def area_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('area-colormode-custom_4', 'style'),
    [Input('area-colormode_4', 'value')]
)

def area_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('area-output', 'data'),
    area_input_list()
)

def store_area(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[9*i]), 'sep': interp_fs(args[9*i+1])
            },
            'hovertextformat': args[9*i+2],
            'radius': {'R0': args[9*i+3], 'R1': args[9*i+4]},
            'layout': {'opacity': args[9*i+5]}
        }
        if args[9*i+6] == 'By Chromosome':
            di['colorcolumn'] = 'ideogram'
        elif args[9*i+6] == 'Mono':
            di['colorcolumn'] = 'None'
            #print('store_area arg')
            #print(args[9*i+7])
            try:
                di['layout']['fillcolor'] = args[9*i+7]['hex']
            except Exception:
                di['layout']['fillcolor'] = 'rgb(241,112,19)'
        elif args[9*i+6] == 'Custom':
            di['colorcolumn'] = args[9*i+8]
        res.append(di)
    #print(res)
    return res


#### tile list of callbacks
@app.callback(
    Output('tile-idx_0', 'style'),
    [Input('tile-number', 'value')]
)
def toggle_tile_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('tile-idx_1', 'style'),
    [Input('tile-number', 'value')]
)
def toggle_tile_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('tile-idx_2', 'style'),
    [Input('tile-number', 'value')]
)
def toggle_tile_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('tile-idx_3', 'style'),
    [Input('tile-number', 'value')]
)
def toggle_tile_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('tile-idx_4', 'style'),
    [Input('tile-number', 'value')]
)
def toggle_tile_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('tile-colormode-mono_0', 'style'),
    [Input('tile-colormode_0', 'value')]
)

def tile_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('tile-colormode-custom_0', 'style'),
    [Input('tile-colormode_0', 'value')]
)

def tile_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('tile-colormode-mono_1', 'style'),
    [Input('tile-colormode_1', 'value')]
)

def tile_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('tile-colormode-custom_1', 'style'),
    [Input('tile-colormode_1', 'value')]
)

def tile_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('tile-colormode-mono_2', 'style'),
    [Input('tile-colormode_2', 'value')]
)

def tile_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('tile-colormode-custom_2', 'style'),
    [Input('tile-colormode_2', 'value')]
)

def tile_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('tile-colormode-mono_3', 'style'),
    [Input('tile-colormode_3', 'value')]
)

def tile_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('tile-colormode-custom_3', 'style'),
    [Input('tile-colormode_3', 'value')]
)

def tile_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('tile-colormode-mono_4', 'style'),
    [Input('tile-colormode_4', 'value')]
)

def tile_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('tile-colormode-custom_4', 'style'),
    [Input('tile-colormode_4', 'value')]
)

def tile_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('tile-output', 'data'),
    tile_input_list()
)

def store_tile(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[10*i]), 'sep': interp_fs(args[10*i+1])
            },
            'hovertextformat': args[10*i+2],
            'radius': {'R0': args[10*i+3], 'R1': args[10*i+4]},
            'layout': {'opacity': args[10*i+5], 
                      'line': {'width': args[10*i+6]} }
        }
        if args[10*i+7] == 'By Chromosome':
            di['colorcolumn'] = 'ideogram'
        elif args[10*i+7] == 'Mono':
            try:
                di['layout']['fillcolor'] = args[10*i+8]['hex']
                di['layout']['line']['color'] = args[10*i+8]['hex']
            except Exception:
                di['layout']['fillcolor'] = 'rgb(241,112,19)'
                di['layout']['line']['color'] = 'rgb(241,112,19)'

        elif args[10*i+7] == 'Custom':
            di['colorcolumn'] = args[10*i+9]
        res.append(di)
    print(res)
    return res


#### heatmap list of callbacks
@app.callback(
    Output('heatmap-idx_0', 'style'),
    [Input('heatmap-number', 'value')]
)
def toggle_heatmap_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('heatmap-idx_1', 'style'),
    [Input('heatmap-number', 'value')]
)
def toggle_heatmap_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('heatmap-idx_2', 'style'),
    [Input('heatmap-number', 'value')]
)
def toggle_heatmap_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('heatmap-idx_3', 'style'),
    [Input('heatmap-number', 'value')]
)
def toggle_heatmap_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('heatmap-idx_4', 'style'),
    [Input('heatmap-number', 'value')]
)
def toggle_heatmap_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}

### ONGOING
### Dash ColorScales component colorscale property => palatte dict in config

@app.callback(
    Output('heatmap-output', 'data'),
    heatmap_input_list()
)
def store_heatmap(number, *args):
    #print('triggering heatmap')
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[9*i]), 'sep': interp_fs(args[9*i+1])
            },
            'hovertextformat': args[9*i+2],
            'radius': {'R0': args[9*i+3], 'R1': args[9*i+4]},
            'layout': {'opacity': args[9*i+5]},
            'palatte': { 'palatte': args[9*i+6], 
                         'reverse': args[9*i+7],
                         'scale': args[9*i+8],
                         }
        }
        res.append(di)
    #print(res)
    
    return res

#### connector list of callbacks
@app.callback(
    Output('connector-idx_0', 'style'),
    [Input('connector-number', 'value')]
)
def toggle_connector_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('connector-idx_1', 'style'),
    [Input('connector-number', 'value')]
)
def toggle_connector_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('connector-idx_2', 'style'),
    [Input('connector-number', 'value')]
)
def toggle_connector_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('connector-idx_3', 'style'),
    [Input('connector-number', 'value')]
)
def toggle_connector_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('connector-idx_4', 'style'),
    [Input('connector-number', 'value')]
)
def toggle_connector_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}

'''
@app.callback(
    Output('connector-colormode-mono_0', 'style'),
    [Input('connector-colormode_0', 'value')]
)

def connector_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('connector-colormode-custom_0', 'style'),
    [Input('connector-colormode_0', 'value')]
)

def connector_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('connector-colormode-mono_1', 'style'),
    [Input('connector-colormode_1', 'value')]
)

def connector_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('connector-colormode-custom_1', 'style'),
    [Input('connector-colormode_1', 'value')]
)

def connector_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('connector-colormode-mono_2', 'style'),
    [Input('connector-colormode_2', 'value')]
)

def connector_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('connector-colormode-custom_2', 'style'),
    [Input('connector-colormode_2', 'value')]
)

def connector_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('connector-colormode-mono_3', 'style'),
    [Input('connector-colormode_3', 'value')]
)

def connector_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('connector-colormode-custom_3', 'style'),
    [Input('connector-colormode_3', 'value')]
)

def connector_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('connector-colormode-mono_4', 'style'),
    [Input('connector-colormode_4', 'value')]
)

def connector_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('connector-colormode-custom_4', 'style'),
    [Input('connector-colormode_4', 'value')]
)

def connector_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
'''

@app.callback(
    Output('connector-output', 'data'),
    connector_input_list()
)
# 9 => 7
def store_connector(number, *args):
    res = []
    N = 7

    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[N*i]), 'sep': interp_fs(args[N*i+1])
            },
            'radius': {'R0': args[N*i+2], 'R1': args[N*i+3]},
            'layout': {'opacity': args[N*i+4], 'line': {'width': args[N*i+5]}}
        }
        try:
            di['layout']['line']['color'] = args[N*i+6]['hex']
        except Exception:
            di['layout']['line']['color'] = 'black'

        #if args[N*i+6] == 'Mono':
            #try:
                #di['layout']['line']['color'] = args[N*i+7]['hex']
            #except Exception:
                #di['layout']['line']['color'] = 'black'
        #elif args[9*i+6] == 'Custom':
            #di['colorcolumn'] = args[9*i+8]
        res.append(di)
    #print(res)
    return res


#### link list of callbacks
@app.callback(
    Output('link-idx_0', 'style'),
    [Input('link-number', 'value')]
)
def toggle_link_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('link-idx_1', 'style'),
    [Input('link-number', 'value')]
)
def toggle_link_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('link-idx_2', 'style'),
    [Input('link-number', 'value')]
)
def toggle_link_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('link-idx_3', 'style'),
    [Input('link-number', 'value')]
)
def toggle_link_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('link-idx_4', 'style'),
    [Input('link-number', 'value')]
)
def toggle_link_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('link-colormode-mono_0', 'style'),
    [Input('link-colormode_0', 'value')]
)

def link_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('link-colormode-custom_0', 'style'),
    [Input('link-colormode_0', 'value')]
)

def link_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('link-colormode-mono_1', 'style'),
    [Input('link-colormode_1', 'value')]
)

def link_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('link-colormode-custom_1', 'style'),
    [Input('link-colormode_1', 'value')]
)

def link_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('link-colormode-mono_2', 'style'),
    [Input('link-colormode_2', 'value')]
)

def link_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('link-colormode-custom_2', 'style'),
    [Input('link-colormode_2', 'value')]
)

def link_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('link-colormode-mono_3', 'style'),
    [Input('link-colormode_3', 'value')]
)

def link_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('link-colormode-custom_3', 'style'),
    [Input('link-colormode_3', 'value')]
)

def link_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('link-colormode-mono_4', 'style'),
    [Input('link-colormode_4', 'value')]
)

def link_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('link-colormode-custom_4', 'style'),
    [Input('link-colormode_4', 'value')]
)

def link_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('link-output', 'data'),
    link_input_list()
)

def store_link(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[11*i]), 'sep': interp_fs(args[11*i+1])
            },
            'hovertextformat': [args[11*i+2], args[11*i+3]],
            'radius': {'R0': args[11*i+4], 'R1': args[11*i+5]},
            'layout': {'opacity': args[11*i+6], 
                      'line': {'width': args[11*i+7]}}
        }
        if args[11*i+8] == 'Mono':
            try:
                di['layout']['line']['color'] = args[11*i+9]['hex']
            except Exception:
                di['layout']['line']['color'] = 'rgb(241,112,19)'
        elif args[11*i+8] == 'Custom':
            di['colorcolumn'] = args[11*i+10]
        res.append(di)
   #print(res)
    return res


#### ribbon list of callbacks
@app.callback(
    Output('ribbon-idx_0', 'style'),
    [Input('ribbon-number', 'value')]
)
def toggle_ribbon_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('ribbon-idx_1', 'style'),
    [Input('ribbon-number', 'value')]
)
def toggle_ribbon_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('ribbon-idx_2', 'style'),
    [Input('ribbon-number', 'value')]
)
def toggle_ribbon_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('ribbon-idx_3', 'style'),
    [Input('ribbon-number', 'value')]
)
def toggle_ribbon_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('ribbon-idx_4', 'style'),
    [Input('ribbon-number', 'value')]
)
def toggle_ribbon_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('ribbon-colormode-mono_0', 'style'),
    [Input('ribbon-colormode_0', 'value')]
)

def ribbon_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('ribbon-colormode-custom_0', 'style'),
    [Input('ribbon-colormode_0', 'value')]
)

def ribbon_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('ribbon-colormode-mono_1', 'style'),
    [Input('ribbon-colormode_1', 'value')]
)

def ribbon_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('ribbon-colormode-custom_1', 'style'),
    [Input('ribbon-colormode_1', 'value')]
)

def ribbon_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('ribbon-colormode-mono_2', 'style'),
    [Input('ribbon-colormode_2', 'value')]
)

def ribbon_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('ribbon-colormode-custom_2', 'style'),
    [Input('ribbon-colormode_2', 'value')]
)

def ribbon_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('ribbon-colormode-mono_3', 'style'),
    [Input('ribbon-colormode_3', 'value')]
)

def ribbon_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('ribbon-colormode-custom_3', 'style'),
    [Input('ribbon-colormode_3', 'value')]
)

def ribbon_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('ribbon-colormode-mono_4', 'style'),
    [Input('ribbon-colormode_4', 'value')]
)

def ribbon_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('ribbon-colormode-custom_4', 'style'),
    [Input('ribbon-colormode_4', 'value')]
)

def ribbon_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('ribbon-output', 'data'),
    ribbon_input_list()
)

def store_ribbon(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[10*i]), 'sep': interp_fs(args[10*i+1])
            },
            'hovertextformat': [args[10*i+2], args[10*i+3]],
            'radius': {'R0': args[10*i+4], 'R1': args[10*i+5]},
            'layout': {'opacity': args[10*i+6]}
        }
        
        if args[10*i+7] == 'Mono':
            try:
                di['layout']['fillcolor'] = args[10*i+8]['hex']
            except Exception:
                di['layout']['fillcolor'] = 'rgb(241,112,19)'
        elif args[10*i+7] == 'Custom':
            di['colorcolumn'] = args[10*i+9]
        res.append(di)
    #print(res)
    return res


#### twistedribbon list of callbacks
@app.callback(
    Output('twistedribbon-idx_0', 'style'),
    [Input('twistedribbon-number', 'value')]
)
def toggle_twistedribbon_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('twistedribbon-idx_1', 'style'),
    [Input('twistedribbon-number', 'value')]
)
def toggle_twistedribbon_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('twistedribbon-idx_2', 'style'),
    [Input('twistedribbon-number', 'value')]
)
def toggle_twistedribbon_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('twistedribbon-idx_3', 'style'),
    [Input('twistedribbon-number', 'value')]
)
def toggle_twistedribbon_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('twistedribbon-idx_4', 'style'),
    [Input('twistedribbon-number', 'value')]
)
def toggle_twistedribbon_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('twistedribbon-colormode-mono_0', 'style'),
    [Input('twistedribbon-colormode_0', 'value')]
)

def twistedribbon_colormode_mono_0(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('twistedribbon-colormode-custom_0', 'style'),
    [Input('twistedribbon-colormode_0', 'value')]
)

def twistedribbon_colormode_custom_0(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('twistedribbon-colormode-mono_1', 'style'),
    [Input('twistedribbon-colormode_1', 'value')]
)

def twistedribbon_colormode_mono_1(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('twistedribbon-colormode-custom_1', 'style'),
    [Input('twistedribbon-colormode_1', 'value')]
)

def twistedribbon_colormode_custom_1(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'} 

@app.callback(
    Output('twistedribbon-colormode-mono_2', 'style'),
    [Input('twistedribbon-colormode_2', 'value')]
)

def twistedribbon_colormode_mono_2(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
    
@app.callback(
    Output('twistedribbon-colormode-custom_2', 'style'),
    [Input('twistedribbon-colormode_2', 'value')]
)

def twistedribbon_colormode_custom_2(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('twistedribbon-colormode-mono_3', 'style'),
    [Input('twistedribbon-colormode_3', 'value')]
)

def twistedribbon_colormode_mono_3(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('twistedribbon-colormode-custom_3', 'style'),
    [Input('twistedribbon-colormode_3', 'value')]
)

def twistedribbon_colormode_custom_3(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}
     

@app.callback(
    Output('twistedribbon-colormode-mono_4', 'style'),
    [Input('twistedribbon-colormode_4', 'value')]
)

def twistedribbon_colormode_mono_4(value):
    if value == 'Mono':
        return {'display': 'block'}
    else: return {'display': 'none'}
     
@app.callback(
    Output('twistedribbon-colormode-custom_4', 'style'),
    [Input('twistedribbon-colormode_4', 'value')]
)

def twistedribbon_colormode_custom_4(value):
    if value == 'Custom':
        return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('twistedribbon-output', 'data'),
    twistedribbon_input_list()
)

def store_twistedribbon(number, *args):
    res = []
    for i in range(number):
        di = {
            'file': { 
                'path': interp_contents(args[10*i]), 'sep': interp_fs(args[10*i+1])
            },
            'hovertextformat': [args[10*i+2], args[10*i+3]],
            'radius': {'R0': args[10*i+4], 'R1': args[10*i+5]},
            'layout': {'opacity': args[10*i+6]}
        }
        
        if args[10*i+7] == 'Mono':
            try:
                di['layout']['fillcolor'] = args[10*i+8]['hex']
            except Exception:
                di['layout']['fillcolor'] = 'rgb(241,112,19)'
        elif args[10*i+7] == 'Custom':
            di['colorcolumn'] = args[10*i+9]
        res.append(di)
    #print(res)
    return res

'''
@app.callback(
    Output('combined-output', 'data'),
    [
        Input('ideogram-output', 'data'),
        Input('cytoband-output', 'data'),
        Input('ring-output', 'data'),
        Input('highlight-output', 'data'),
        Input('annotation-output', 'data'),
        Input('histogram-output', 'data'),
        Input('scatter-output', 'data'),
        Input('line-output', 'data'),
        Input('area-output', 'data'),
        Input('tile-output', 'data'),
        Input('heatmap-output', 'data'),
        Input('connector-output', 'data'),
        Input('link-output', 'data'),
        Input('ribbon-output', 'data'),
        Input('twistedribbon-output', 'data'),
    ]
)
def merge_all(*args):
    # it seems that I cannot remove keys with empty list otherwise a runtimeerror would occur saying dictionary changed size during iteration
    res = {
            'General': {'width': 1400, 'height': 1400},
            'Category': {
                'ideogram': args[0],
                'cytoband': args[1],
                'ring': args[2],
                'highlight': args[3],
                'annotation': args[4],
                'histogram': args[5],
                'scatter': args[6],
                'line': args[7],
                'area': args[8],
                'tile': args[9],
                'heatmap': args[10],
                'connector': args[11],
                'link': args[12],
                'ribbon': args[13],
                'twistedribbon': args[14]
                }
            }
    print('printing combined-all')
    print(res)
    
    return res    
'''

@app.callback(
    Output('histogram-number', 'disabled'),
    [
        Input('histogram-number-lock', 'on')
    ]
)
def histogram_number_lock(bool_value):
    return bool_value

@app.callback(
    Output('scatter-number', 'disabled'),
    [
        Input('scatter-number-lock', 'on')
    ]
)

def scatter_number_lock(bool_value):
    return bool_value


@app.callback(
    Output('line-number', 'disabled'),
    [
        Input('line-number-lock', 'on')
    ]
)
def line_number_lock(bool_value):
    return bool_value




@app.callback(
    Output('area-number', 'disabled'),
    [
        Input('area-number-lock', 'on')
    ]
)
def area_number_lock(bool_value):
    return bool_value


@app.callback(
    Output('tile-number', 'disabled'),
    [
        Input('tile-number-lock', 'on')
    ]
)
def tile_number_lock(bool_value):
    return bool_value



@app.callback(
    Output('heatmap-number', 'disabled'),
    [
        Input('heatmap-number-lock', 'on')
    ]
)
def heatmap_number_lock(bool_value):
    return bool_value


@app.callback(
    Output('connector-number', 'disabled'),
    [
        Input('connector-number-lock', 'on')
    ]
)
def connector_number_lock(bool_value):
    return bool_value


@app.callback(
    Output('link-number', 'disabled'),
    [
        Input('link-number-lock', 'on')
    ]
)
def link_number_lock(bool_value):
    return bool_value


@app.callback(
    Output('ribbon-number', 'disabled'),
    [
        Input('ribbon-number-lock', 'on')
    ]
)
def ribbon_number_lock(bool_value):
    return bool_value



@app.callback(
    Output('twistedribbon-number', 'disabled'),
    [
        Input('twistedribbon-number-lock', 'on')
    ]
)
def twistedribbon_number_lock(bool_value):
    return bool_value



@app.callback(
    Output('PCircos_update', 'figure'),
    dash_input(),
    [
        State('ideogram-output', 'data'),
        State('cytoband-output', 'data'),
        State('ring-output', 'data'),
        State('highlight-output', 'data'),
        State('annotation-output', 'data'),
        State('histogram-output', 'data'),
        State('scatter-output', 'data'),
        State('line-output', 'data'),
        State('area-output', 'data'),
        State('tile-output', 'data'),
        State('heatmap-output', 'data'),
        State('connector-output', 'data'),
        State('link-output', 'data'),
        State('ribbon-output', 'data'),
        State('twistedribbon-output', 'data'),
    ]
)
def merge_all(n_clicks, *args):
    # it seems that I cannot remove keys with empty list otherwise a runtimeerror would occur saying dictionary changed size during iteration
    if n_clicks is None:
        raise PreventUpdate

    res = {
            'General': {'width': 1400, 'height': 1400},
            'Category': {
                'ideogram': args[0],
                'cytoband': args[1],
                'ring': args[2],
                'highlight': args[3],
                'annotation': args[4],
                'histogram': args[5],
                'scatter': args[6],
                'line': args[7],
                'area': args[8],
                'tile': args[9],
                'heatmap': args[10],
                'connector': args[11],
                'link': args[12],
                'ribbon': args[13],
                'twistedribbon': args[14]
                }
            }
    
    #print('before updating dash_dict:')
    #print(res)
    #print('updating dash_dict state')
    

    remove_empty(res)
    convert_dict(res)
    #print('after convert_dict')
    #print(res)

    #print('debugging ring')
    #print(args[2])

    return Figure(dash_dict=res).fig()   




'''
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
'''


if __name__ == '__main__':
    #print (app.layout)
    try:
        app.run_server(debug=True, host='0.0.0.0', port=sys.argv[1])
    except Exception:
        app.run_server(debug=True, host='0.0.0.0', port=8050)

# python3 dashapp.py 