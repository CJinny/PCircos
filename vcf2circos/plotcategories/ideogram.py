from vcf2circos.plotcategories.plotconfig import Plotconfig
import pandas as pd
from os.path import join as osj

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

    def __init__(
        self,
        filename,
        options,
        show,
        file,
        radius,
        sortbycolor,
        colorcolumn,
        hovertextformat,
        trace_car,
        data,
        layout,
    ):
        super().__init__(
            filename,
            options,
            show,
            file,
            radius,
            sortbycolor,
            colorcolumn,
            hovertextformat,
            trace_car,
            data,
            layout,
        )
        self.chr_conf = pd.read_csv(
            osj(self.options["Static"], "chr_size.txt"), sep="\t", header=0
        )
        self.data = self.process_vcf()
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
        self.customoptions = (
            {"customlabel": "True", "customspacing": "False", "customcolor": 3,},
        )
        self.npoints = (1000,)
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

    def adapt_data(self):
        data = {
            "chr_name": self.data["Chromosomes"],
            "chr_size": self.chr_conf.loc[
                self.chr_conf["chr_label"] == self.data["Chromosomes"]["chr_size"]
            ].to_list(),
            "chr_label": self.data["Chromosomes"],
            "chr_color": self.chr_conf.loc[
                self.chr_conf["chr_label"] == self.data["Chromosomes"]["chr_color"]
            ].to_list(),
        }
        data = self.data["Chromosomes"]
        return data

    def merge_options(self):
        dico = {}
        dico["patch"] = {}
        # ideo = Ideogram()
        dico["patch"]["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": self.adapt_data()},
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
