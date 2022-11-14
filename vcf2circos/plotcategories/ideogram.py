from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import Colorpal
import pandas as pd
from os.path import join as osj
import os

# space, space between ring in option.example.json
# height hauteur du ring
# positon position from center

# just for test
# data = {
#    "chr_name": ["chr1", "chr2", "chr3"],
#    "chr_size": [249250621, 243199373, 198022430],
#    "chr_label": ["chr1", "chr2", "chr3"],
#    "chr_color": ["pink", "rosybrown", "firebrick"],
# }

list_graph_type = ["majortick", "minortick", "ticklabel"]


class Ideogram(Plotconfig):
    """
    "scatter":{
        "pattern":{
        ...
    },  "data":{
        ...
    }}
    """

    def __init__(self, plotconfig):
        self.plotconfig = plotconfig
        # assert os.path.exists(osj(self.options["Static"], "chr_size.txt"))
        self.chr_conf = pd.read_csv(
            osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                "chr." + self.options["Assembly"] + ".sorted.txt",
            ),
            sep="\t",
            header=0,
        )
        self.degreerange = [0, 360]
        self.showfillcolor = self.cast_bool(True)
        self.chrannotation = (
            {
                "show": "True",
                "radius": {"R": 1.25},
                "fonttype": "bold",
                "textangle": {"angleoffset": 0, "anglelimit": 360},
                "layout": {
                    "xref": "x",
                    "yref": "y",
                    "showarrow": False,
                    "font": {"size": 10, "color": "black"},
                },
            },
        )
        self.customoptions = {
            "customlabel": "True",
            "customspacing": "False",
            "customcolor": 3,
        }
        self.npoints = 1000
        self.radius = {"R0": 1.0, "R1": 1.1}
        self.layout = {
            "type": "path",
            "opacity": 0.9,
            "layer": "above",
            "line": {"color": "gray", "width": 2},
        }
        self.majortick = {
            "show": "True",
            "spacing": 30000000,
            "radius": {"R0": 1.1, "R1": 1.125},
            "layout": {
                "type": "path",
                "opacity": 0.9,
                "layer": "above",
                "line": {"color": "black", "width": 1},
            },
        }
        self.minortick = {
            "show": "True",
            "spacing": 5000000,
            "radius": {"R0": 1.1, "R1": 1.118},
            "layout": {
                "type": "path",
                "opacity": 0.9,
                "line": {"color": "black", "width": 0.5},
            },
        }
        self.ticklabel = {
            "show": "True",
            "spacing": 30000000,
            "radius": {"R": 1.16},
            "textformat": "Mb",
            "textangle": {"angleoffset": -90, "anglelimit": 360},
            "layout": {
                "xref": "x",
                "yref": "y",
                "showarrow": False,
                "font": {"family": "Times New Roman", "size": 8, "color": "black",},
            },
        }

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def data_ideogram(self):
        tmp = self.chr_conf.loc[
            self.chr_conf["chr_name"].isin(self.data["Chromosomes"])
        ]
        data = {
            "chr_name": tmp["chr_name"].to_list(),
            "chr_size": tmp["size"].to_list(),
            "chr_label": tmp["chr_name"].to_list(),
            "chr_color": list(Colorpal(len(tmp["chr_name"].to_list()))),
        }
        print(data["chr_name"])
        print(data["chr_color"])
        # exit()
        return data

    def merge_options(self):
        dico = {}
        dico["patch"] = {}
        # ideo = Ideogram()
        dico["patch"]["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": self.data_ideogram()},
        }
        dico["patch"]["show"] = self.show
        dico["patch"]["degreerange"] = self.degreerange
        dico["patch"]["showfillcolor"] = self.showfillcolor
        dico["patch"]["chrannotation"] = {
            "show": "True",
            "radius": {"R": 1.25},
            "fonttype": "bold",
            "textangle": {"angleoffset": 0, "anglelimit": 360},
            "layout": {
                "xref": "x",
                "yref": "y",
                "showarrow": False,
                "font": {"size": 10, "color": "black"},
            },
        }
        dico["patch"]["customoptions"] = {
            "customlabel": "True",
            "customspacing": "False",
            "customcolor": 3,
        }
        dico["patch"]["npoints"] = self.npoints
        dico["patch"]["radius"] = self.radius
        dico["patch"]["layout"] = self.layout
        dico["majortick"] = self.majortick
        dico["minortick"] = self.minortick
        dico["ticklabel"] = self.ticklabel
        return dico
        # Loopable as fuck
