
import dash_core_components as dcc
import dash_html_components as html
import dash_colorscales
from colorpicker_box import ColorPickerBox
import io
import base64
import colors
import numpy as np
import dash
from dash.dependencies import Input, Output, State

UPLOADBOX_STYLES = {'display': 'block', 'position': 'relative', 'width': '95%', 'height': '40px',
                    'line-height': '40px', 'border-width': '1px', 'border-style': 'dashed',
                    'border-radius': '5px', 'text-align': 'center', 'marginTop': '1.2em', 
                    'marginBottom': '0.8em', 'overflow': 'auto'
                    }
FS_STYLES = {'marginTop': '0.6em', 'marginBottom': '0.4em'}
RdBu = ['#67001f', '#b2182b', '#d6604d', '#f4a582',  '#fddbc7', '#f7f7f7', '#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061']

def interp_contents(contents):
    if contents == None:
        return None
   
    else:
        _, content_string = contents.split(',')
        #decoded = base64.b64decode(content_string)
        #path = io.StringIO(decoded.decode('utf-8'))
        return content_string


def interp_fs(fs):
    assert fs in ['Tab', 'Blank', 'Comma']
    if fs == 'Tab':
        return '\t'
    elif fs == 'Blank':
        return ' '
    else:
        return ','

def interp_minmax(min, max):
    if min <= max:
        return min, max
    else:
        return max, min

# likely not needed
def interp_bool(bool_str):
    assert bool_str in ['True', 'False']
    if bool_str == 'True':
        return True
    else:
        return False

def convert_rgb(rgb_str):
    # convert 'rgb(255,0,255)' => {'rgb': {'r': 255, 'g': 0, 'b': 255, 'a': 1}}
    rgb_numeric = np.array(rgb_str.strip('rgb()').split(','))
    return {'rgb': {'r': rgb_numeric[0], 'g': rgb_numeric[1], 'b': rgb_numeric[2], 'a': 1}}

def revert_rgb(rgb_dict):
    # revert {'rgb': {'r': 255, 'g': 0, 'b': 255, 'a': 1}} => 'rgb(255,0,255)'
    assert 'rgb' in rgb_dict
    rgb_sub = rgb_dict['rgb']
    return 'rgb({},{},{})'.format(rgb_sub['r'], rgb_sub['g'], rgb_sub['b'])


def expand_chromosome_color(chromosome_checklist):

    default_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrx','chry','chr23','chr24','chr0','chrM','chrUn','chrNA']
    if len(chromosome_checklist) <= 30:
        initial_chr_color = colors.to_rgb(default_list)[:len(chromosome_checklist)]
    else:
        assert len(chromosome_checklist) <= 710
        initial_chr_color = colors.random_rgb(len(chromosome_checklist))
    
    def convert_rgb(rgb_str):
        # convert 'rgb(255,0,255)' => {'rgb': {'r': 255, 'g': 0, 'b': 255, 'a': 1}}
        rgb_numeric = np.array(rgb_str.strip('rgb()').split(','))
        return {'rgb': {'r': rgb_numeric[0], 'g': rgb_numeric[1], 'b': rgb_numeric[2], 'a': 1}}
    initial_chr_rgb = [*map(lambda a: convert_rgb(a), initial_chr_color)]


    # ColorPickerBox(id='fillcolor_ring_{}'.format(i), value={'rgb': {'r': 255, 'g': 0, 'b': 255, 'a': 1}}, style={'marginTop': '0.4em', 'marginBottom': '0.4em'}
    component_list = []
    for i in range(len(chromosome_checklist)):
        component_list.append(ColorPickerBox(id='chr_{}'.format(chromosome_checklist[i]), 
                                             value=initial_chr_rgb[i],
                                             label='{}'.format(chromosome_checklist[i]), 
                                             style={'display': 'inline-block','margin': '0.2em'})
                                            )

    return html.Div(component_list)
    


def expand_ring(n_clicks):
    # consider using Input number to replace rangeslider

    def single_element(i): 
        return  html.Div([
                    html.P('Select ring radius range for ring {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='ring-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='ring-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'} 
                        ),
                    ]),
                    html.P('Select opacity for ring {}'.format(i+1)),
                    dcc.Input(
                        id='ring-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.1,
                        style={'marginBottom': '1.2em'}
                    ),
                    html.Div([
                        html.Div('Select fillcolor for ring {}'.format(i+1)),
                        ColorPickerBox(id='ring-fillcolor_{}'.format(i), value={'rgb': {'r': 255, 'g': 0, 'b': 255, 'a': 1}}, style={'marginTop': '0.4em', 'marginBottom': '0.4em'}),
                        
                    ], style={'display': 'block'}
                    ),
                ], className='indent'
                )
    
    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])


# for createing inputs under callbacks in store_histogram
def histogram_input_list():
    input_list = [Input('histogram-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('histogram-upload_{}'.format(i), 'contents'),
            Input('histogram-fs_{}'.format(i), 'value'),
            Input('histogram-hovertext_{}'.format(i), 'value'),
            Input('histogram-range-min_{}'.format(i), 'value'),
            Input('histogram-range-max_{}'.format(i), 'value'),
            Input('histogram-opacity_{}'.format(i), 'value'),
            Input('histogram-colormode_{}'.format(i), 'value'),
            Input('histogram-colorpickerbox_{}'.format(i), 'value'),
            Input('histogram-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list

def scatter_input_list():
    input_list = [Input('scatter-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('scatter-upload_{}'.format(i), 'contents'),
            Input('scatter-fs_{}'.format(i), 'value'),
            Input('scatter-hovertext_{}'.format(i), 'value'),
            Input('scatter-range-min_{}'.format(i), 'value'),
            Input('scatter-range-max_{}'.format(i), 'value'),
            Input('scatter-opacity_{}'.format(i), 'value'),
            Input('scatter-markersize_{}'.format(i), 'value'),
            Input('scatter-markersymbol_{}'.format(i), 'value'),
            Input('scatter-colormode_{}'.format(i), 'value'),
            Input('scatter-colorpickerbox_{}'.format(i), 'value'),
            Input('scatter-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list


def line_input_list():
    input_list = [Input('line-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('line-upload_{}'.format(i), 'contents'),
            Input('line-fs_{}'.format(i), 'value'),
            Input('line-hovertext_{}'.format(i), 'value'),
            Input('line-range-min_{}'.format(i), 'value'),
            Input('line-range-max_{}'.format(i), 'value'),
            Input('line-opacity_{}'.format(i), 'value'),
            Input('line-markersize_{}'.format(i), 'value'),
            Input('line-linewidth_{}'.format(i), 'value'),
            Input('line-smoothing_{}'.format(i), 'value'),
            Input('line-colormode_{}'.format(i), 'value'),
            Input('line-colorpickerbox_{}'.format(i), 'value'),
            Input('line-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list


def area_input_list():
    input_list = [Input('area-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('area-upload_{}'.format(i), 'contents'),
            Input('area-fs_{}'.format(i), 'value'),
            Input('area-hovertext_{}'.format(i), 'value'),
            Input('area-range-min_{}'.format(i), 'value'),
            Input('area-range-max_{}'.format(i), 'value'),
            Input('area-opacity_{}'.format(i), 'value'),
            Input('area-colormode_{}'.format(i), 'value'),
            Input('area-colorpickerbox_{}'.format(i), 'value'),
            Input('area-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list

def tile_input_list():
    input_list = [Input('tile-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('tile-upload_{}'.format(i), 'contents'),
            Input('tile-fs_{}'.format(i), 'value'),
            Input('tile-hovertext_{}'.format(i), 'value'),
            Input('tile-range-min_{}'.format(i), 'value'),
            Input('tile-range-max_{}'.format(i), 'value'),
            Input('tile-opacity_{}'.format(i), 'value'),
            Input('tile-linewidth_{}'.format(i), 'value'),
            Input('tile-colormode_{}'.format(i), 'value'),
            Input('tile-colorpickerbox_{}'.format(i), 'value'),
            Input('tile-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list

def heatmap_input_list():
    input_list = [Input('heatmap-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('heatmap-upload_{}'.format(i), 'contents'),
            Input('heatmap-fs_{}'.format(i), 'value'),
            Input('heatmap-hovertext_{}'.format(i), 'value'),
            Input('heatmap-range-min_{}'.format(i), 'value'),
            Input('heatmap-range-max_{}'.format(i), 'value'),
            Input('heatmap-opacity_{}'.format(i), 'value'),
            Input('heatmap-palette_{}'.format(i), 'colorscale'),
            Input('heatmap-palatte-reverse_{}'.format(i), 'value'),
            Input('heatmap-palatte-scale_{}'.format(i), 'value'),
        ])
    return input_list

def connector_input_list():
    input_list = [Input('connector-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('connector-upload_{}'.format(i), 'contents'),
            Input('connector-fs_{}'.format(i), 'value'),
            Input('connector-range-min_{}'.format(i), 'value'),
            Input('connector-range-max_{}'.format(i), 'value'),
            Input('connector-opacity_{}'.format(i), 'value'),
            Input('connector-linewidth_{}'.format(i), 'value'),
            Input('connector-colormode_{}'.format(i), 'value'),
            Input('connector-colorpickerbox_{}'.format(i), 'value'),
            Input('connector-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list

def link_input_list():
    input_list = [Input('link-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('link-upload_{}'.format(i), 'contents'),
            Input('link-fs_{}'.format(i), 'value'),
            Input('link-hovertext_0_{}'.format(i), 'value'),
            Input('link-hovertext_1_{}'.format(i), 'value'),
            Input('link-range-min_{}'.format(i), 'value'),
            Input('link-range-max_{}'.format(i), 'value'),
            Input('link-opacity_{}'.format(i), 'value'),
            Input('link-linewidth_{}'.format(i), 'value'),
            Input('link-colormode_{}'.format(i), 'value'),
            Input('link-colorpickerbox_{}'.format(i), 'value'),
            Input('link-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list

def ribbon_input_list():
    input_list = [Input('ribbon-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('ribbon-upload_{}'.format(i), 'contents'),
            Input('ribbon-fs_{}'.format(i), 'value'),
            Input('ribbon-hovertext_0_{}'.format(i), 'value'),
            Input('ribbon-hovertext_1_{}'.format(i), 'value'),
            Input('ribbon-range-min_{}'.format(i), 'value'),
            Input('ribbon-range-max_{}'.format(i), 'value'),
            Input('ribbon-opacity_{}'.format(i), 'value'),
            Input('ribbon-colormode_{}'.format(i), 'value'),
            Input('ribbon-colorpickerbox_{}'.format(i), 'value'),
            Input('ribbon-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list

def twistedribbon_input_list():
    input_list = [Input('twistedribbon-number', 'value')]
    for i in range(5):
        input_list.extend([
            Input('twistedribbon-upload_{}'.format(i), 'contents'),
            Input('twistedribbon-fs_{}'.format(i), 'value'),
            Input('twistedribbon-hovertext_0_{}'.format(i), 'value'),
            Input('twistedribbon-hovertext_1_{}'.format(i), 'value'),
            Input('twistedribbon-range-min_{}'.format(i), 'value'),
            Input('twistedribbon-range-max_{}'.format(i), 'value'),
            Input('twistedribbon-opacity_{}'.format(i), 'value'),
            Input('twistedribbon-colormode_{}'.format(i), 'value'),
            Input('twistedribbon-colorpickerbox_{}'.format(i), 'value'),
            Input('twistedribbon-colorcolumn_{}'.format(i), 'value'),
        ])
    return input_list


#### State list
def histogram_State_list():
    State_list = [State('histogram-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('histogram-upload_{}'.format(i), 'contents'),
            State('histogram-fs_{}'.format(i), 'value'),
            State('histogram-hovertext_{}'.format(i), 'value'),
            State('histogram-range-min_{}'.format(i), 'value'),
            State('histogram-range-max_{}'.format(i), 'value'),
            State('histogram-opacity_{}'.format(i), 'value'),
            State('histogram-colormode_{}'.format(i), 'value'),
            State('histogram-colorpickerbox_{}'.format(i), 'value'),
            State('histogram-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def scatter_State_list():
    State_list = [State('scatter-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('scatter-upload_{}'.format(i), 'contents'),
            State('scatter-fs_{}'.format(i), 'value'),
            State('scatter-hovertext_{}'.format(i), 'value'),
            State('scatter-range-min_{}'.format(i), 'value'),
            State('scatter-range-max_{}'.format(i), 'value'),
            State('scatter-opacity_{}'.format(i), 'value'),
            State('scatter-markersize_{}'.format(i), 'value'),
            State('scatter-markersymbol_{}'.format(i), 'value'),
            State('scatter-colormode_{}'.format(i), 'value'),
            State('scatter-colorpickerbox_{}'.format(i), 'value'),
            State('scatter-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def line_State_list():
    State_list = [State('line-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('line-upload_{}'.format(i), 'contents'),
            State('line-fs_{}'.format(i), 'value'),
            State('line-hovertext_{}'.format(i), 'value'),
            State('line-range-min_{}'.format(i), 'value'),
            State('line-range-max_{}'.format(i), 'value'),
            State('line-opacity_{}'.format(i), 'value'),
            State('line-markersize_{}'.format(i), 'value'),
            State('line-linewidth_{}'.format(i), 'value'),
            State('line-smoothing_{}'.format(i), 'value'),
            State('line-colormode_{}'.format(i), 'value'),
            State('line-colorpickerbox_{}'.format(i), 'value'),
            State('line-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def area_State_list():
    State_list = [State('area-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('area-upload_{}'.format(i), 'contents'),
            State('area-fs_{}'.format(i), 'value'),
            State('area-hovertext_{}'.format(i), 'value'),
            State('area-range-min_{}'.format(i), 'value'),
            State('area-range-max_{}'.format(i), 'value'),
            State('area-opacity_{}'.format(i), 'value'),
            State('area-colormode_{}'.format(i), 'value'),
            State('area-colorpickerbox_{}'.format(i), 'value'),
            State('area-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def tile_State_list():
    State_list = [State('tile-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('tile-upload_{}'.format(i), 'contents'),
            State('tile-fs_{}'.format(i), 'value'),
            State('tile-hovertext_{}'.format(i), 'value'),
            State('tile-range-min_{}'.format(i), 'value'),
            State('tile-range-max_{}'.format(i), 'value'),
            State('tile-opacity_{}'.format(i), 'value'),
            State('tile-linewidth_{}'.format(i), 'value'),
            State('tile-colormode_{}'.format(i), 'value'),
            State('tile-colorpickerbox_{}'.format(i), 'value'),
            State('tile-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def heatmap_State_list():
    State_list = [State('heatmap-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('heatmap-upload_{}'.format(i), 'contents'),
            State('heatmap-fs_{}'.format(i), 'value'),
            State('heatmap-hovertext_{}'.format(i), 'value'),
            State('heatmap-range-min_{}'.format(i), 'value'),
            State('heatmap-range-max_{}'.format(i), 'value'),
            State('heatmap-opacity_{}'.format(i), 'value'),
            State('heatmap-palette_{}'.format(i), 'colorscale'),
            State('heatmap-palatte-reverse_{}'.format(i), 'value'),
            State('heatmap-palatte-scale_{}'.format(i), 'value'),
        ])
    return State_list

def connector_State_list():
    State_list = [State('connector-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('connector-upload_{}'.format(i), 'contents'),
            State('connector-fs_{}'.format(i), 'value'),
            State('connector-range-min_{}'.format(i), 'value'),
            State('connector-range-max_{}'.format(i), 'value'),
            State('connector-opacity_{}'.format(i), 'value'),
            State('connector-linewidth_{}'.format(i), 'value'),
            State('connector-colormode_{}'.format(i), 'value'),
            State('connector-colorpickerbox_{}'.format(i), 'value'),
            State('connector-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def link_State_list():
    State_list = [State('link-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('link-upload_{}'.format(i), 'contents'),
            State('link-fs_{}'.format(i), 'value'),
            State('link-hovertext_0_{}'.format(i), 'value'),
            State('link-hovertext_1_{}'.format(i), 'value'),
            State('link-range-min_{}'.format(i), 'value'),
            State('link-range-max_{}'.format(i), 'value'),
            State('link-opacity_{}'.format(i), 'value'),
            State('link-linewidth_{}'.format(i), 'value'),
            State('link-colormode_{}'.format(i), 'value'),
            State('link-colorpickerbox_{}'.format(i), 'value'),
            State('link-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def ribbon_State_list():
    State_list = [State('ribbon-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('ribbon-upload_{}'.format(i), 'contents'),
            State('ribbon-fs_{}'.format(i), 'value'),
            State('ribbon-hovertext_0_{}'.format(i), 'value'),
            State('ribbon-hovertext_1_{}'.format(i), 'value'),
            State('ribbon-range-min_{}'.format(i), 'value'),
            State('ribbon-range-max_{}'.format(i), 'value'),
            State('ribbon-opacity_{}'.format(i), 'value'),
            State('ribbon-colormode_{}'.format(i), 'value'),
            State('ribbon-colorpickerbox_{}'.format(i), 'value'),
            State('ribbon-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list

def twistedribbon_State_list():
    State_list = [State('twistedribbon-number', 'value')]
    for i in range(5):
        State_list.extend([
            State('twistedribbon-upload_{}'.format(i), 'contents'),
            State('twistedribbon-fs_{}'.format(i), 'value'),
            State('twistedribbon-hovertext_0_{}'.format(i), 'value'),
            State('twistedribbon-hovertext_1_{}'.format(i), 'value'),
            State('twistedribbon-range-min_{}'.format(i), 'value'),
            State('twistedribbon-range-max_{}'.format(i), 'value'),
            State('twistedribbon-opacity_{}'.format(i), 'value'),
            State('twistedribbon-colormode_{}'.format(i), 'value'),
            State('twistedribbon-colorpickerbox_{}'.format(i), 'value'),
            State('twistedribbon-colorcolumn_{}'.format(i), 'value'),
        ])
    return State_list



def expand_histogram():
    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='histogram-upload_{}'.format(i),
                        children=html.Div(['Upload histogram {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for histogram {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='histogram-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format for histogram {}'.format(i+1)),
                    dcc.Textarea(
                        id='histogram-hovertext_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Start: {}<br>End: {}<br>LogFC: {:.4f}\".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for histogram {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='histogram-range-min_{}'.format(i),
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='histogram-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for histogram {}'.format(i+1)),
                    dcc.Input(
                        id='histogram-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.1,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for histogram {}'.format(i+1)),
                    dcc.RadioItems(
                        id='histogram-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By Chromosome', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'}
                    ),

                    html.Div([
                        html.P('Select color for histogram {}'.format(i+1)),
                        ColorPickerBox(id='histogram-colorpickerbox_{}'.format(i), 
                                       label='', 
                                       style={'display': 'hidden','margin': '0.2em'}
                                       )
                        ], id='histogram-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for histogram {}'.format(i+1)),
                        dcc.Input(id='histogram-colorcolumn_{}'.format(i),
                                  type='number',
                                  min=4,
                                  value=4,
                                  step=1,
                                  style={'width': '45%'}
                                )
                        ], id='histogram-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),

        ], id='histogram-idx_{}'.format(i), style={'display': 'none'}, className='indent')
        
    return html.Div([single_element(i) for i in range(5)])


def expand_scatter():

    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='scatter-upload_{}'.format(i),
                        children=html.Div(['Upload scatter {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for scatter {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='scatter-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format for scatter {}'.format(i+1)),
                    dcc.Textarea(
                        id='scatter-hovertext_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Position: {}<br>Value: {:.2f}\".format(a[i,0], a[i,1], float(a[i,2])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for scatter {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='scatter-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='scatter-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for scatter {}'.format(i+1)),
                    dcc.Input(
                        id='scatter-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.9,
                        style={'width': '45%'}
                    ),
                    html.P('Select marker size (px) for scatter {}'.format(i+1)),
                    dcc.Input(
                        id='scatter-markersize_{}'.format(i),
                        type='number',
                        min=1,
                        max=10,
                        step=1,
                        value=5,
                        style={'width': '45%'}
                    ),
                    html.P('Select marker symbol for scatter {}'.format(i+1)),
                    dcc.Input(
                        id='scatter-markersymbol_{}'.format(i),
                        type='number',
                        min=0,
                        max=44,
                        step=1,
                        value=0,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for scatter {}'.format(i+1)),
                    dcc.RadioItems(
                        id='scatter-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By Chromosome', 'Custom']],
                        value='By Chromosome',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for scatter {}'.format(i+1)),
                        ColorPickerBox(id='scatter-colorpickerbox_{}'.format(i), 
                                        label='', 
                                        style={'display': 'hidden','margin': '0.2em'}
                                        )
                            ], id='scatter-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for scatter {}'.format(i+1)),
                        dcc.Input(id='scatter-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=4,
                                    value=4,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='scatter-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                ], id='scatter-idx_{}'.format(i), style={'display': 'none'}, className='indent')

    
    return html.Div([single_element(i) for i in range(5)])


def expand_line():

    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='line-upload_{}'.format(i),
                        children=html.Div(['Upload line {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for line {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='line-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format for line {}'.format(i+1)),
                    dcc.Textarea(
                        id='line-hovertext_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Position: {}<br>Value: {:.2f}\".format(a[i,0], a[i,1], float(a[i,2])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for line {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='line-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='line-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for line {}'.format(i+1)),
                    dcc.Input(
                        id='line-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.9,
                        style={'width': '45%'}
                    ),
                    html.P('Select marker size (px) for line {}'.format(i+1)),
                    dcc.Input(
                        id='line-markersize_{}'.format(i),
                        type='number',
                        min=1,
                        max=10,
                        step=1,
                        value=2,
                        style={'width': '45%'}
                    ),
                    html.P('Select line width (px) for line {}'.format(i+1)),
                    dcc.Input(
                        id='line-linewidth_{}'.format(i),
                        type='number',
                        min=0,
                        max=10,
                        step=1,
                        value=2,
                        style={'width': '45%'}
                    ),
                    html.P('Select smoothing for line {}'.format(i+1)),
                    dcc.Input(
                        id='line-smoothing_{}'.format(i),
                        type='number',
                        min=0,
                        max=1.3,
                        step=0.1,
                        value=0,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for line {}'.format(i+1)),
                    dcc.RadioItems(
                        id='line-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By Chromosome']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for line {}'.format(i+1)),
                        ColorPickerBox(id='line-colorpickerbox_{}'.format(i), 
                                        label='', 
                                        style={'display': 'hidden','margin': '0.2em'}
                                        )
                            ], id='line-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for line {}'.format(i+1)),
                        dcc.Input(id='line-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=4,
                                    value=4,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='line-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                ], id='line-idx_{}'.format(i), style={'display': 'none'}, className='indent')

    return html.Div([single_element(i) for i in range(5)])


def expand_area():

    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='area-upload_{}'.format(i),
                        children=html.Div(['Upload area {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for area {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='area-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format for area {}'.format(i+1)),
                    dcc.Textarea(
                        id='area-hovertext_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Position: {}<br>Value: {:.2f}\".format(a[i,0], a[i,1], float(a[i,2])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for area {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='area-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='area-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for area {}'.format(i+1)),
                    dcc.Input(
                        id='area-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.4,
                        style={'width': '45%'}
                    ),
                    
                    html.P('Select color mode for area {}'.format(i+1)),
                    dcc.RadioItems(
                        id='area-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By Chromosome', 'Custom']],
                        value='By Chromosome',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for area {}'.format(i+1)),
                        ColorPickerBox(id='area-colorpickerbox_{}'.format(i), 
                                        label='', 
                                        style={'display': 'hidden','margin': '0.2em'}
                                        )
                            ], id='area-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for area {}'.format(i+1)),
                        dcc.Input(id='area-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=4,
                                    value=4,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='area-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                ], id='area-idx_{}'.format(i), style={'display': 'none'}, className='indent')

    return html.Div([single_element(i) for i in range(5)])




def expand_tile():

    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='tile-upload_{}'.format(i),
                        children=html.Div(['Upload tile {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for tile {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='tile-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format for tile {}'.format(i+1)),
                    dcc.Textarea(
                        id='tile-hovertext_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Start: {}<br>End: {}<br>CNV: {:.4f}\".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for tile {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='tile-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='tile-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for tile {}'.format(i+1)),
                    dcc.Input(
                        id='tile-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.9,
                        style={'width': '45%'}
                    ),
                    html.P('Select line width (px) for line {}'.format(i+1)),
                    dcc.Input(
                        id='tile-linewidth_{}'.format(i),
                        type='number',
                        min=1,
                        max=10,
                        step=1,
                        value=3,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for tile {}'.format(i+1)),
                    dcc.RadioItems(
                        id='tile-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By Chromosome', 'Custom']],
                        value='By Chromosome',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for tile {}'.format(i+1)),
                        ColorPickerBox(id='tile-colorpickerbox_{}'.format(i), 
                                        label='', 
                                        style={'display': 'hidden','margin': '0.2em'}
                                        )
                            ], id='tile-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for tile {}'.format(i+1)),
                        dcc.Input(id='tile-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=4,
                                    value=4,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='tile-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                ], id='tile-idx_{}'.format(i), style={'display': 'none'}, className='indent')

    return html.Div([single_element(i) for i in range(5)])


def expand_heatmap():

    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='heatmap-upload_{}'.format(i),
                        children=html.Div(['Upload heatmap {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for heatmap {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='heatmap-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format for heatmap {}'.format(i+1)),
                    dcc.Textarea(
                        id='heatmap-hovertext_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}\".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for heatmap {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='heatmap-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='heatmap-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for heatmap {}'.format(i+1)),
                    dcc.Input(
                        id='heatmap-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.8,
                        style={'width': '45%'}
                    ),
                    html.P('Select colorscale for heatmap {}'.format(i+1)),
                    dash_colorscales.DashColorscales(
                        id='heatmap-palette_{}'.format(i),
                        colorscale=RdBu,
                        nSwatches=11,
                        fixSwatches=False
                    ),
                    html.P('Reverse colorscale for heatmap {}'.format(i+1)),
                    dcc.RadioItems(
                        id='heatmap-palatte-reverse_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['True', 'False']],
                        value='True',
                        labelStyle={'display': 'inline-block'}
                    ),
                    html.P('Choose type of colorscale for heatmap {}'.format(i+1)),
                    dcc.RadioItems(
                        id='heatmap-palatte-scale_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['div', 'seq']],
                        value='True',
                    )
                    
                ], id='heatmap-idx_{}'.format(i), style={'display': 'none'}, className='indent'
            )

    return html.Div([single_element(i) for i in range(5)])

def expand_connector():

    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='connector-upload_{}'.format(i),
                        children=html.Div(['Upload connector {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for connector {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='connector-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Select radius range for connector {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='connector-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='connector-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for connector {}'.format(i+1)),
                    dcc.Input(
                        id='connector-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=1,
                        style={'width': '45%'}
                    ),
                    html.P('Select line width for connector {}'.format(i+1)),
                    dcc.Input(
                        id='connector-linewidth_{}'.format(i),
                        type='number',
                        min=0.5,
                        max=10.0,
                        step=0.5,
                        value=1,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for connector {}'.format(i+1)),
                    dcc.RadioItems(
                        id='connector-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for connector {}'.format(i+1)),
                        ColorPickerBox(id='connector-colorpickerbox_{}'.format(i), 
                                       label='', 
                                       style={'display': 'hidden', 'margin': '0.2em'}
                                        )
                            ], id='connector-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for connector {}'.format(i+1)),
                        dcc.Input(id='connector-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=3,
                                    value=3,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='connector-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),

                    # I disabled custom segment radius ratio
                    
                ], id='connector-idx_{}'.format(i), style={'display': 'none'}, className='indent'
            )

   
    return html.Div([single_element(i) for i in range(5)])

def expand_link():
    # disable By chromosome color mode
    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='link-upload_{}'.format(i),
                        children=html.Div(['Upload link {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for link {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='link-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format (1st) for link {}'.format(i+1)),
                    dcc.Textarea(
                        id='link-hovertext_0_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}\".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Input hovertext format (2nd) for link {}'.format(i+1)),
                    dcc.Textarea(
                        id='link-hovertext_1_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}\".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for link {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='link-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='link-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for link {}'.format(i+1)),
                    dcc.Input(
                        id='link-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.8,
                        style={'width': '45%'}
                    ),
                    html.P('Select line width (px) for link {}'.format(i+1)),
                    dcc.Input(
                        id='link-linewidth_{}'.format(i),
                        type='number',
                        min=1,
                        max=10,
                        step=1,
                        value=3,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for link {}'.format(i+1)),
                    dcc.RadioItems(
                        id='link-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for link {}'.format(i+1)),
                        ColorPickerBox(id='link-colorpickerbox_{}'.format(i), 
                                       label='', 
                                       style={'display': 'hidden', 'margin': '0.2em'}
                                        )
                            ], id='link-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for link {}'.format(i+1)),
                        dcc.Input(id='link-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=6,
                                    value=6,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='link-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),

                ], id='link-idx_{}'.format(i), style={'display': 'none'}, className='indent'
            )

    return html.Div([single_element(i) for i in range(5)]) 
    

def expand_ribbon():
    # disable By chromosome color mode
    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='ribbon-upload_{}'.format(i),
                        children=html.Div(['Upload ribbon {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for ribbon {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='ribbon-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format (1st) for ribbon {}'.format(i+1)),
                    dcc.Textarea(
                        id='ribbon-hovertext_0_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}\".format(a[i,0], a[i,3], a[i,1], a[i,5], a[i,2], a[i,4]) ',
                        style={'width': '100%'}
                    ),
                    html.P('Input hovertext format (2nd) for ribbon {}'.format(i+1)),
                    dcc.Textarea(
                        id='ribbon-hovertext_1_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}\".format(a[i,3], a[i,0], a[i,5], a[i,1], a[i,4], a[i,2]) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for ribbon {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='ribbon-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='ribbon-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for ribbon {}'.format(i+1)),
                    dcc.Input(
                        id='ribbon-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.6,
                        style={'width': '45%'}
                    ),
                    
                    html.P('Select color mode for ribbon {}'.format(i+1)),
                    dcc.RadioItems(
                        id='ribbon-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for ribbon {}'.format(i+1)),
                        ColorPickerBox(id='ribbon-colorpickerbox_{}'.format(i), 
                                       label='', 
                                       style={'display': 'hidden', 'margin': '0.2em'}
                                        )
                            ], id='ribbon-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for ribbon {}'.format(i+1)),
                        dcc.Input(id='ribbon-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=6,
                                    value=6,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='ribbon-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),                    
                ], id='ribbon-idx_{}'.format(i), style={'display': 'none'}, className='indent')

    return html.Div([single_element(i) for i in range(5)]) 
    
def expand_twistedribbon():
    # disable By chromosome color mode
    def single_element(i):
        return html.Div([
                    html.Hr(),
                    dcc.Upload(
                        id='twistedribbon-upload_{}'.format(i),
                        children=html.Div(['Upload twisted ribbon {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose field separator for twisted ribbon {}'.format(i+1), style=FS_STYLES),
                    dcc.RadioItems(
                        id='twistedribbon-fs_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Tab', 'Blank', 'Comma']],
                        value='Tab',
                        labelStyle={'dispaly': 'inline-block'}
                    ),
                    html.P('Input hovertext format (1st) for twisted ribbon {}'.format(i+1)),
                    dcc.Textarea(
                        id='twistedribbon-hovertext_0_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}\".format(a[i,0], a[i,3], a[i,1], a[i,4], a[i,2], a[i,5]) ',
                        style={'width': '100%'}
                    ),
                    html.P('Input hovertext format (2nd) for twisted ribbon {}'.format(i+1)),
                    dcc.Textarea(
                        id='twistedribbon-hovertext_1_{}'.format(i),
                        placeholder='Enter hovertext format...',
                        value=' \"Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}\".format(a[i,3], a[i,0], a[i,4], a[i,1], a[i,5], a[i,2]) ',
                        style={'width': '100%'}
                    ),
                    html.P('Select radius range for twisted ribbon {}'.format(i+1)),
                    html.Div([
                        dcc.Input(
                            id='twistedribbon-range-min_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        ),
                        html.P('-', style={'display': 'inline-block'}),
                        dcc.Input(
                            id='twistedribbon-range-max_{}'.format(i), 
                            type='number', 
                            value=0, 
                            step=0.05, 
                            min=0, 
                            max=2, 
                            style={'width': '45%'}
                        )
                    ]),
                    html.P('Select opacity for twisted ribbon {}'.format(i+1)),
                    dcc.Input(
                        id='twistedribbon-opacity_{}'.format(i),
                        type='number',
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.6,
                        style={'width': '45%'}
                    ),
                    
                    html.P('Select color mode for twisted ribbon {}'.format(i+1)),
                    dcc.RadioItems(
                        id='twistedribbon-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div([
                        html.P('Select color for twistedribbon {}'.format(i+1)),
                        ColorPickerBox(id='twistedribbon-colorpickerbox_{}'.format(i), 
                                       label='', 
                                       style={'display': 'hidden', 'margin': '0.2em'}
                                        )
                            ], id='twistedribbon-colormode-mono_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),
                    html.Div([
                        html.P('Input color column for twisted ribbon {}'.format(i+1)),
                        dcc.Input(id='twistedribbon-colorcolumn_{}'.format(i),
                                    type='number',
                                    min=6,
                                    value=6,
                                    step=1,
                                    style={'width': '45%'}
                                )
                        ], id='twistedribbon-colormode-custom_{}'.format(i), style={'display': 'none'}, className='indent'
                    ),                    
                    
                ], id='twistedribbon-idx_{}'.format(i), style={'display': 'none'}, className='indent')

    return html.Div([single_element(i) for i in range(5)]) 





