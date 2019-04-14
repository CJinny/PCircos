import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import itertools
import json
from dashapp_functions import *
import base64
import datetime
import io
import numpy as np
import pandas as pd

app = dash.Dash()
app.config.supress_callback_exceptions = True


def create_outputs():
    def single_output(i):
        return html.Div([
            
            html.P('output_{}'.format(i)),
            dcc.Input('input_{}'.format(i), type='number', min=0, max=3, step=0.5, value=0),
            dcc.RadioItems('radio_{}'.format(i), options=[{'label': i, 'value': i} for i in ['Mono', 'By chromosome', 'Custom']], value='Custom'),
            
        ], id='single_output_{}'.format(i), 
        style={'display': 'inline-block'}
        )

    return html.Div([single_output(i) for i in range(5)])



app.layout = html.Div([
    html.P('Type the number of components to show'),
    dcc.Input(id='number', type='number', value=0, step=1, min=0, max=5),
    html.Div([create_outputs()]),
    dcc.Store(id='output_print'),
    html.Div([expand_histogram()]),

    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ])),
    html.Div(id='output-upload')
])

def parse_contents(contents):
    try:
        _, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t', header='infer')
        print(df.head())
        
        return html.Div([str(df.iloc[0,0])])
    except Exception:
        return html.Div('Error occur')
    

@app.callback(
    Output('output-upload', 'children'),
    [Input('upload-data', 'contents')]
)
def upload_output(contents):
    _, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t', header='infer')
    print(df.tail())
    return parse_contents(contents)






@app.callback(
    Output('single_output_0', 'style'),
    [Input('number', 'value')]
)

def toggle_0(number):
    if number >= 1: return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('single_output_1', 'style'),
    [Input('number', 'value')]
)

def toggle_1(number):
    if number >= 2: return {'display': 'block'}
    else: return {'display': 'none'}

@app.callback(
    Output('single_output_2', 'style'),
    [Input('number', 'value')]
)

def toggle_2(number):
    if number >= 3: return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('single_output_3', 'style'),
    [Input('number', 'value')]
)

def toggle_3(number):
    if number >= 4: return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('single_output_4', 'style'),
    [Input('number', 'value')]
)

def toggle_4(number):
    if number >= 5: return {'display': 'block'}
    else: return {'display': 'none'}


@app.callback(
    Output('output_print', 'data'),
    [   
        Input('number', 'value'),
        Input('input_0', 'value'),
        Input('radio_0', 'value'),
        Input('input_1', 'value'),
        Input('radio_1', 'value'),
        Input('input_2', 'value'),
        Input('radio_2', 'value'),
        Input('input_3', 'value'),
        Input('radio_3', 'value'),
        Input('input_4', 'value'),
        Input('radio_4', 'value'),

    ]
)
def display_result(number, *arg):
    res = []
    for i in range(number):
        di = {'input': arg[2*i], 
              'radio': arg[2*i+1]
              }
        res.append(di)
        #res.append(arg[2*i+1])
        #res.append({'input': arg[i*2], 'radio': arg[i*2+1]})
    print(res)
    return res




'''
for i in range(20):
    @app.callback(
        Output('single_output_{}'.format(i), 'style'),
        [Input('number', 'value')]

    )
    def toggle_hide_i(number):
        
        if i <= number:
            
            return {'display': 'inline-block'}

        else:
            print(i, number)
            return {'display': 'block'}
'''

if __name__ == '__main__':
    app.run_server(debug=True, port=8051)
