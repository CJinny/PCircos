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
__version__ = '1.1'
__date__ = 'April 18 2019'

if (sys.version_info[0]!=3):
    raise Exception('PCircos requires Python 2, your current Python version is {}.{}.{}'.
                    format(sys.version_info[0],sys.version_info[1],sys.version_info[2]))


def run_PCircos():

    fig_instance = Figure(input_json_path=sys.argv[1])


    fig = fig_instance.fig()
    try:
        filename = sys.argv[2]
        plot(fig, filename=filename)
    except IndexError:
        plot(fig)

if __name__ == "__main__":
    t=time()
    run_PCircos()
    print ('total run time:')
    print (time()-t)

    # python3 PCircos.py demo_data/demo_params.json > test.txt
