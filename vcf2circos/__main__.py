# enable dash only in dash module, disable dash here
import os
import sys
from vcf2circos.parseargs import Parseargs
from vcf2circos.datafactory import Datafactory
from pprint import pprint
from vcf2circos.utils import launch


sys.path.append(os.path.join(os.path.dirname(__file__), "."))
sys.path.append(os.path.abspath(os.path.join("../", "demo_data")))
# print(sys.path)

from itertools import count
import numpy as np
from collections import defaultdict

# import colorlover as cl

from IPython.display import HTML
from plotly.graph_objs import *
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

# init_notebook_mode(connected=True)

from Complex import Complex
from fig import Figure
from time import time

from os.path import join as osj
import plotly.io
import plotly
import pathlib
import json
import copy

import argparse


def main():
    t = time()
    args = Parseargs().parseargs()
    # Input
    input_file = args.input
    if os.path.exists(input_file):
        input_format = pathlib.Path(input_file).suffix.replace(".", "")
        print(f"[INFO] Input file: {input_file} (format '{input_format}')")
    else:
        print(f"[ERROR] No input file '{input_file}")
        exit()

    # Output
    output_file = args.output
    output_format = pathlib.Path(output_file).suffix.replace(".", "")

    if output_format not in [
        "html",
        "png",
        "jpg",
        "jpeg",
        "webp",
        "svg",
        "pdf",
        "eps",
        "json",
    ]:
        print(f"[ERROR] Output file format '{output_file}' not supported")
    print(f"[INFO] Output file: {output_file} (format '{output_format}')")

    # Export
    export_format = args.export
    if export_format:
        print(f"[INFO] Export format: {export_format}")
        output_export_file = ".".join(output_file.split(".")[:-1])+"."+export_format
        print(f"[INFO] Export file: {output_export_file}")
    # Options
    options_input = args.options or {}
    if options_input:
        if os.path.isfile(options_input):
            with open(options_input, "r") as f:
                options = json.loads(f.read())
            f.close()
        else:
            try:
                options = json.loads(options_input)
            except IndexError:
                options = {}
    else:
        options = {}
    if options:
        print(f"[INFO] Options provided.")
        options["File"] = options_input
        if args.assembly:
            assert args.assembly in os.listdir(osj(options["Static"], "Assembly")), (
                "ERROR genome assembly "
                + args.assembly
                + " is not available update \n your config or choose an available among this list "
                + ", ".join(
                    [
                        assemb
                        for assemb in os.listdir(osj(options["Static"], "Assembly"))
                    ]
                )
            )
            options["Assembly"] = args.assembly
    else:
        print(f"[INFO] Options not provided.")

    # Notebook mode
    # notebook_mode = args.notebook_mode
    # if notebook_mode:
    #    init_notebook_mode(connected=True)

    # Input

    if input_format in ["vcf", "gz"]:
        print("\n")
        js = Datafactory(input_file, options).plot_dict()
        fig_instance = Figure(dash_dict=js, options=options)
        # Export in vcf2circos JSON
        #if output_export_file:
        #    if not os.path.exists(os.path.dirname(os.path.abspath(output_export_file))):
        #        os.mkdir(os.path.dirname(output_export_file))
        #    with open(output_export_file, "w") as ef:
        #        ef.write(json.dumps(copy.deepcopy(js), indent=4))
        #        ef.close()

    elif input_format in ["json"]:
        fig_instance = Figure(input_json_path=input_file, options=options)
    else:

        print("[ERROR] input format not supported")

    fig = fig_instance.fig()
    fig["layout"]["showlegend"] = True
    fig.update_layout(legend_title="Legend")
    fig.update_layout(legend_xanchor="right")
    fig.update_layout(legend_x=1.1)
    # Order of the legend test
    dico = defaultdict(
        int,
        {
            k: v
            for k, v in {
                "GENES": 1,
                "MORBID_GENES": 2,
                "CNV": 3,
                "DEL": 4,
                "DUP": 5,
                "INS": 6,
                "INV": 7,
                "BND": 8,
                "SNV": 9,
            }.items()
        },
    )
    for scatter in fig["data"]:
        if scatter.showlegend is None and hasattr(scatter, "name"):
            if scatter.name is not None:
                scatter.legendrank = dico[scatter.name]
    try:
        if not output_format:
            plot(fig)
        if not os.path.exists(os.path.dirname(os.path.abspath(output_file))):
            os.mkdir(os.path.dirname(output_file))
        elif output_format in ["html"]:
            plot(fig, filename=output_file)
        else:
            print("[ERROR] output format not supported")
        if export_format in [
            "png",
            "jpg",
            "jpeg",
            "webp",
            "svg",
            "pdf",
            "eps",
            "json",
        ]:
            #plotly.io.write_image(fig, output_export_file, format="svg")
            fig.write_image(output_export_file)
    except IndexError:
        plot(fig)

    print("[INFO] total run time:")
    print("[INFO] " + str(time() - t) + " sec")


if __name__ == "__main__":
    main()
