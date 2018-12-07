# enable dash only in dash module, disable dash here



import numpy as np
import colorlover as cl

from config import coord_config, json_config
from IPython.display import HTML
from plotly.graph_objs import *
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
init_notebook_mode(connected=True)
import sys
import argparse
from argparse import ArgumentParser

from Complex import Complex
from fig import Figure
import maths
import colors
from time import time


__author__ = 'Jin Cui'
__version__ = '1.0.0'
__date__ = 'August 27 2018'

if (sys.version_info[0]!=3):
    raise Exception('PCircos requires Python 2, your current Python version is {}.{}.{}'.
                    format(sys.version_info[0],sys.version_info[1],sys.version_info[2]))


def run_PCircos():
    #print (json_config.default_dict())
    #print (sys.argv[1])
    #print ("input_json_dict: \t")
    #print (json_config.json2dict(sys.argv[1]))
    #print ("\t")



    fig_instance = Figure(input_json_path=sys.argv[1])


    ################
    # ATTEMPT TO MODIFY DEGREERANGE HERE AND SEE RESPONSE
    ##############

    #print (fig_instance.get_chr_info()['ideogram_bin'])
    #print ("\t")
    #print (fig_instance.get_data_array('histogram'))
    #print ("fig instance hovertext is: \t")
    #print (fig_instance.get_hovertext('histogram'))
    #print (fig_instance.get_data_complexes('scatter'))
    #print ("\t")
    #print ('trace is: \t')
    #print (fig_instance.trace())
    #print ('\t')
    #print ('layout is: \t')
    #print (fig_instance.layout())
    #print ('examine attribute categories: \t')
    #print (fig_instance.categories)
    #print ('examine ideogram_coord_config: \t')
    #print (fig_instance.ideogram_ideogram)
    #print ('\t')
    #print ('examine ideogram_coord_config: \t')
    #print (fig_instance.get_ideogram_coord_config())
    #print (fig_instance.layout())
    #print ('ideogram_path is: \t')
    #print (fig_instance.ideogram_path(fig_instance.get_ideogram_complex()))
    #print (fig_instance.get_data_complexes('link').reshape((2,6)))
    #print (fig_instance.get_paths_dict('link'))
    #print (fig_instance.get_hovertext('link'))


    #print ('histogram data complexes are: \t')
    #print (fig_instance.get_data_complexes('histogram'))
    #print ('histogram hovertexts are: \t')
    #print (fig_instance.get_hovertext('histogram'))

    #print ('histogram data complexes are: \t')
    #print (fig_instance.get_data_complexes('histogram'))
    
    #print ('histogram hovertexts are: \t')
    #print (fig_instance.get_hovertext('histogram'))

    #print ('histogram traces are: \t')
    #print (fig_instance.get_traces('histogram'))
    '''
    print ('get_data_array_dict(heatmap): \t')
    print (fig_instance.get_data_array_dict('heatmap'))
    print ('get_data_arrays for heatmap: \t')
    print (fig_instance.get_data_array('heatmap'))
    print ('get_read_data for heatmap is: \t')
    print (fig_instance.get_read_data('heatmap'))
    print ('get_hovertext for heatmap: \t')
    print (fig_instance.get_hovertext('heatmap'))

    print ('get_data_complex_array for heatmap: \t')
    print (fig_instance.get_data_complexes('heatmap'))'''

    #print ('get_data path for heatmap is: \t')
    #print (fig_instance.get_paths_dict('heatmap'))
    #print ('\t')
    #print ('\t')
    #print ('get_traces for heatmap: \t')
    #print (fig_instance.get_traces('heatmap'))

    #print ('fig_instance.config_dict is: \t')
    #print (fig_instance.config_dict)



    #print ('fig_instance categories is ')
    #print (fig_instance.categories['annotation'])

    

    fig = fig_instance.fig()
    try:
        filename = sys.argv[3]
        plot(fig, filename=filename)
    except IndexError:
        plot(fig)

if __name__ == "__main__":
    t=time()
    run_PCircos()
    print ('total run time:')
    print (time()-t)

    # python3 PCircos.py demo_data/demo_params.json > test.txt
