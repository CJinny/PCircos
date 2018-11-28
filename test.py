import sys
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import plotly.graph_objs as go
import numpy as np
import pandas as pd

app = dash.Dash(__name__)
app.layout = html.Div([
    html.Div([
        html.P('Select checklist'),
        dcc.Checklist(
            id='checklist',
            options=[
                {'label': 'a', 'value': 'a'},
                {'label': 'b', 'value': 'b'},
                {'label': 'c', 'value': 'c'},
                {'label': 'd', 'value': 'd'},
            ],
            values=['a','b','d'],
            labelStyle={'display':'inline-block'}
        )
    ]),
    html.Div(id='result')
])

@app.callback(
    Output(component_id='result', component_property='children'),
    [Input(component_id='checklist', component_property='values')]
)
def update_result(input_value):
    return 'You have selected {}'.format(" ".join(input_value))

if __name__ == '__main__':
    app.run_server()