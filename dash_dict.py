import numpy as np
import pandas as pd
import colorlover as cl
import dash_core_components as dcc
import dash_html_components as html
import dash_colorscales
import sys
sys.path.append('../')
import colors
from dashapp_functions import *
import io
import base64
from config import json_config, coord_config

'''
This code reads combined-output data and update Circos plot 
'''


'''
{'Category': {
    'ideogram': 
        {
            "file": {"content_string": null, 
                     "header": "infer", "sep": " "}, 
            "degreerange": [12, 348], 
            "showfillcolor": true, 
            "layout": {"opacity": 0.6}, 
            "chrannotation": {"show": true, 
                             "radius": {"R": 1.25}, 
                             "fonttype": "bold", 
                             "textangle": {"angleoffset": -90, "anglelimit": 360}, 
                             "layout": {"font": {"size": 15, "color": "rgb(0,0,0)"}}}, 
                             "majortick": {"show": true, "spacing": 20000000}, 
                             "minortick": {"show": true, "spacing": 5000000}, 
                             "ticklabel": {"show": true, "spacing": 20000000, 
                             "textformat": "Kb", "textangle": {"angleoffset": -90, "anglelimit": 360}}}', 
                             'cytoband': '{"show": false}', 
                             'ring': None, 
                             'highlight': '{"show": false}', 
                             'annotation': '{"show": false}', 
                             
                             'histogram': [
                                 {'file': {'content_string': None, 'sep': '\t'}, 
                                            'hovertextformat': ' "Chromosome: {}<br>Start: {}<br>End: {}<br>LogFC: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ', 
                                            'radius': {'R0': 0.05, 'R1': 0.2}, 'layout': {'opacity': 0.15, 'fillcolor': {'hex': '#583822', 'rgb': {'r': 88, 'g': 56, 'b': 34, 'a': 1}}}}], 
                                            
                            'scatter': [{'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': ' "Chromosome: {}<br>Position: {}<br>Value: {:.2f}".format(a[i,0], a[i,1], float(a[i,2])) ', 'radius': {'R0': 0, 'R1': 0}, 'trace': {'marker': {'size': 5, 'opacity': 0.9, 'symbol': 0}}, 'colorcolumn': 'ideogram'}, {'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': ' "Chromosome: {}<br>Position: {}<br>Value: {:.2f}".format(a[i,0], a[i,1], float(a[i,2])) ', 'radius': {'R0': 0, 'R1': 0}, 'trace': {'marker': {'size': 5, 'opacity': 0.9, 'symbol': 0}}, 'colorcolumn': 'ideogram'}], 
                             'line': [{'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': ' "Chromosome: {}<br>Position: {}<br>Value: {:.2f}".format(a[i,0], a[i,1], float(a[i,2])) ', 'radius': {'R0': 0, 'R1': 0}, 'marker': {'size': 2, 'color': {'hex': '#56723a', 'rgb': {'r': 86, 'g': 114, 'b': 58, 'a': 1}}}, 'line': {'width': 2, 'opacity': 0.9, 'smoothing': 0, 'color': {'hex': '#56723a', 'rgb': {'r': 86, 'g': 114, 'b': 58, 'a': 1}}}}, {'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': ' "Chromosome: {}<br>Position: {}<br>Value: {:.2f}".format(a[i,0], a[i,1], float(a[i,2])) ', 'radius': {'R0': 0, 'R1': 0}, 'marker': {'size': 2}, 'line': {'width': 2, 'opacity': 0.9, 'smoothing': 0}, 'colorcolumn': 'ideogram'}], 'area': [{'file': {'content_string': None, 'sep': ' '}, 'hovertextformat': ' "Chromosome: {}<br>Position: {}<br>Value: {:.2f}".format(a[i,0], a[i,1], float(a[i,2])) ', 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.4}, 'colorcolumn': 'ideogram'}, {'file': {'content_string': None, 'sep': ','}, 'hovertextformat': ' "Chromosome: {}<br>Position: {}<br>Value: {:.2f}".format(a[i,0], a[i,1], float(a[i,2])) ', 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.4}, 'colorcolumn': 'ideogram'}], 'tile': [{'file': {'content_string': None, 'sep': ','}, 'hovertextformat': ' "Chromosome: {}<br>Start: {}<br>End: {}<br>CNV: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ', 'radius': {'R0': 0.1, 'R1': 0.2}, 'layout': {'opacity': 0.9, 'line': {'width': 3}}, 'colorcolumn': 'ideogram'}, {'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': ' "Chromosome: {}<br>Start: {}<br>End: {}<br>CNV: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ', 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.9, 'line': {'width': 3}}, 'colorcolumn': 'ideogram'}], 'heatmap': [{'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ', 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.8}, 'palatte': {'palatte': ['#fff7ec', '#feebcf', '#fddcaf', '#fdca93', '#fdb27a', '#fc8d59', '#f26d4b', '#e1482f', '#c92113', '#a80001', '#7f0000'], 'reverse': 'True', 'scale': 'True'}}, {'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ', 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.8}, 'palatte': {'palatte': ['#8e0152', '#c51b7d', '#de77ae', '#f1b6da', '#fde0ef', '#f7f7f7', '#e6f5d0', '#b8e186', '#7fbc41', '#4d9221', '#276419'], 'reverse': 'True', 'scale': 'True'}}], 'connector': [{'file': {'content_string': None, 'sep': '\t'}, 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 1, 'line': {'width': 1, 'color': {'hex': '#452f49', 'rgb': {'r': 69, 'g': 47, 'b': 73, 'a': 1}}}}}, {'file': {'content_string': None, 'sep': '\t'}, 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 1, 'line': {'width': 1}}, 'colorcolumn': 5}], 'link': [{'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': [' "Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ', ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) '], 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.8, 'line': {'width': 3, 'color': None}}}, {'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': [' "Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ', ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Value: {:.4f}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) '], 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.8, 'line': {'width': 3, 'color': None}}}], 'ribbon': [{'file': {'content_string': None, 'sep': ','}, 'hovertextformat': [' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,0], a[i,3], a[i,1], a[i,5], a[i,2], a[i,4]) ', ' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,3], a[i,0], a[i,5], a[i,1], a[i,4], a[i,2]) '], 'radius': {'R0': 0.1, 'R1': 0.15}, 'layout': {'opacity': 0.6, 'fillcolor': None}}, {'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': [' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,0], a[i,3], a[i,1], a[i,5], a[i,2], a[i,4]) ', ' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,3], a[i,0], a[i,5], a[i,1], a[i,4], a[i,2]) '], 'radius': {'R0': 0.1, 'R1': 0.15}, 'layout': {'opacity': 0.6, 'fillcolor': None}}, {'file': {'content_string': None, 'sep': ' '}, 'hovertextformat': [' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,0], a[i,3], a[i,1], a[i,5], a[i,2], a[i,4]) ', ' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,3], a[i,0], a[i,5], a[i,1], a[i,4], a[i,2]) '], 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.6, 'fillcolor': None}}], 'twistedribbon': [{'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': [' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,0], a[i,3], a[i,1], a[i,4], a[i,2], a[i,5]) ', ' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,3], a[i,0], a[i,4], a[i,1], a[i,5], a[i,2]) '], 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.6}, 'colorcolumn': 9}, {'file': {'content_string': None, 'sep': '\t'}, 'hovertextformat': [' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,0], a[i,3], a[i,1], a[i,4], a[i,2], a[i,5]) ', ' "Chromosome: {} => {}<br>From: {} => {}<br>From: {} => {}".format(a[i,3], a[i,0], a[i,4], a[i,1], a[i,5], a[i,2]) '], 'radius': {'R0': 0, 'R1': 0}, 'layout': {'opacity': 0.6}, 'colorcolumn': 8}]}}
'''

def remove_empty(dash_dict):
    if dash_dict is None:
        pass
    else:
        if not isinstance(dash_dict, dict):
            pass
        else:
            for key in list(dash_dict.keys()):
                if dash_dict[key] in [None, []]:
                    del dash_dict[key] 

                elif key in ['cytoband', 'annotation', 'highlight']:
                    if 'file' not in dash_dict[key]:
                        del dash_dict[key]
                elif isinstance(dash_dict[key], dict):
                    remove_empty(dash_dict[key])
                    
                elif isinstance(dash_dict[key], list):
                    for i in range(len(dash_dict[key])):
                        remove_empty(dash_dict[key][i])
                        


def convert_dict(dash_dict):
    if dash_dict is None:
        pass
    else:

        if not isinstance(dash_dict, dict):
            pass
        else:

            for key in list(dash_dict.keys()):
                if key == 'path':
                    if dash_dict[key] is None:
                        pass
                    else:
                        decoded = base64.b64decode(dash_dict[key])
                        path = io.StringIO(decoded.decode('utf-8'))
                        dash_dict[key] = path
                    
                elif dash_dict[key] in [None, []]:
                    del dash_dict[key]
                   
                
                elif isinstance(dash_dict[key], dict):
                    convert_dict(dash_dict[key])
                    
                elif isinstance(dash_dict[key], list):
                    for i in range(len(dash_dict[key])):
                        convert_dict(dash_dict[key][i])
                
    

