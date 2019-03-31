import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

app = dash.Dash()
server = app.server

app.layout = html.Div([
                       
                       html.Section(id="slideshow", children=[
                                                              html.Div(id="slideshow-container", children=[
                                                                                                           html.Div(id="image"),
                                                                                                           dcc.Interval(id='interval', interval=3000)
                                                                                                           ])
                                                              ])
                       
                       ])

@app.callback(Output('image', 'children'),
              [Input('interval', 'n_intervals')])
def display_image(n):
    if n == None or n % 3 == 1:
        img = html.Img(src="http://placeimg.com/625/225/any")
    elif n % 3 == 2:
        img = html.Img(src="http://placeimg.com/625/225/animals")
    elif n % 3 == 0:
        img = html.Img(src="http://placeimg.com/625/225/arch")
    else:
        img = "None"
    return img

if __name__ == '__main__':
    app.run_server(debug=True)
