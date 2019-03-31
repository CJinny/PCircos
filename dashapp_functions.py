
import dash_core_components as dcc
import dash_html_components as html
import dash_colorscales
from colorpicker_box import ColorPickerBox
import io
import base64
import colors
import numpy as np

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


'''
def generate_histogram_upload_id(i):
    return 'histogram-upload_{}'.format(i)

def generate_histogram_fs_id(i):
    return 'histogram-fs_{}'.format(i)

def generate_histogram_hovertext_id(i):
    return 'histogram-hovertext_{}'.format(i)

def generate_histogram_range_min_id(i):
    return 'histogram-range-min_{}'.format(i)

def generate_histogram_range_max_id(i):
    return 'histogram-range-max_{}'.format(i)

def generate_histogram_opacity_id(i):
    return 'histogram-opacity_{}'.format(i)

def generate_histogram_colormode_id(i):
    return 'histogram-colormode_{}'.format(i)
'''



def expand_histogram(n_clicks):
    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='histogram-upload_{}'.format(i),
                        children=html.Div(['Upload histogram {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for histogram {}'.format(i+1), style=FS_STYLES),
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
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By chromosome', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'}
                    ),
                    #html.Div(id='histogram-colormode-choice_{}'.format(i), className='indent')
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])


def expand_scatter(n_clicks):

    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='scatter-upload_{}'.format(i),
                        children=html.Div(['Upload scatter {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for scatter {}'.format(i+1), style=FS_STYLES),
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
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By chromosome', 'Custom']],
                        value='By chromosome',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div(id='scatter-colormode-choice_{}'.format(i), className='indent')
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])


def expand_line(n_clicks):

    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='line-upload_{}'.format(i),
                        children=html.Div(['Upload line {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for line {}'.format(i+1), style=FS_STYLES),
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
                        id='line-markersymbol_{}'.format(i),
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
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By chromosome']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div(id='line-colormode-choice_{}'.format(i), className='indent')
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])


def expand_area(n_clicks):

    def single_element(i):
        return html.Div([
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
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By chromosome', 'Custom']],
                        value='By chromosome',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div(id='area-colormode-choice_{}'.format(i), className='indent')
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])




def expand_tile(n_clicks):

    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='tile-upload_{}'.format(i),
                        children=html.Div(['Upload tile {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for tile {}'.format(i+1), style=FS_STYLES),
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
                        id='line-markersymbol_{}'.format(i),
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
                        options=[{'label': i, 'value': i} for i in ['Mono', 'By chromosome', 'Custom']],
                        value='By chromosome',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div(id='tile-colormode-choice_{}'.format(i), className='indent')
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])


def expand_heatmap(n_clicks):

    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='heatmap-upload_{}'.format(i),
                        children=html.Div(['Upload heatmap {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for heatmap {}'.format(i+1), style=FS_STYLES),
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
                    
                    html.P('Select palette for heatmap {}'.format(i+1)),
                    dash_colorscales.DashColorscales(
                        id='heatmap-palette_{}'.format(i),
                        colorscale=RdBu,
                        nSwatches=11,
                        fixSwatches=False
                    )
                    
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])

def expand_connector(n_clicks):

    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='connector-upload_{}'.format(i),
                        children=html.Div(['Upload connector {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for connector {}'.format(i+1), style=FS_STYLES),
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
                    # I disabled custom segment radius ratio
                    
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)])

def expand_link(n_clicks):
    # disable By chromosome color mode
    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='link-upload_{}'.format(i),
                        children=html.Div(['Upload link {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for link {}'.format(i+1), style=FS_STYLES),
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
                        id='line-markersymbol_{}'.format(i),
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
                    html.Div(id='link-colormode-choice_{}'.format(i), className='indent')
                    
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)]) 
    

def expand_ribbon(n_clicks):
    # disable By chromosome color mode
    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='ribbon-upload_{}'.format(i),
                        children=html.Div(['Upload ribbon {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for ribbon {}'.format(i+1), style=FS_STYLES),
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
                    html.P('Select line width (px) for ribbon {}'.format(i+1)),
                    dcc.Input(
                        id='line-markersymbol_{}'.format(i),
                        type='number',
                        min=1,
                        max=10,
                        step=1,
                        value=3,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for ribbon {}'.format(i+1)),
                    dcc.RadioItems(
                        id='ribbon-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div(id='ribbon-colormode-choice_{}'.format(i), className='indent')
                    
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)]) 
    
def expand_twistedribbon(n_clicks):
    # disable By chromosome color mode
    def single_element(i):
        return html.Div([
                    dcc.Upload(
                        id='twistedribbon-upload_{}'.format(i),
                        children=html.Div(['Upload twisted ribbon {} here'.format(i+1)]),
                        style=UPLOADBOX_STYLES,
                    ),
                    html.P('Choose a field separator for twisted ribbon {}'.format(i+1), style=FS_STYLES),
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
                    html.P('Select line width (px) for twisted ribbon {}'.format(i+1)),
                    dcc.Input(
                        id='line-markersymbol_{}'.format(i),
                        type='number',
                        min=1,
                        max=10,
                        step=1,
                        value=3,
                        style={'width': '45%'}
                    ),
                    html.P('Select color mode for twisted ribbon {}'.format(i+1)),
                    dcc.RadioItems(
                        id='twistedribbon-colormode_{}'.format(i),
                        options=[{'label': i, 'value': i} for i in ['Mono', 'Custom']],
                        value='Mono',
                        labelStyle={'display': 'inline-block'},
                    ),
                    html.Div(id='twistedribbon-colormode-choice_{}'.format(i), className='indent')
                    
                ], className='indent')

    if n_clicks is None:
        return html.Div([single_element(1)])
    else:    
        return html.Div([single_element(i) for i in range(n_clicks)]) 


