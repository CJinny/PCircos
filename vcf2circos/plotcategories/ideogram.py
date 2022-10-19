from vcf2circos.plotcategories.plotconfig import Plotconfig
import inspect

# space, space between ring in option.example.json
# height hauteur du ring
# positon position from center


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
        config_ring=None,
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

        self.ring = config_ring
        self.variants_position = self.ring["position"]
        self.variants_ring_space = self.ring["space"]
        self.variants_ring_height = self.ring["height"]

        self.ringval = [
            {
                # SNV
                "radius": {
                    "R0": self.variants_position
                    + (7 * self.variants_ring_space)
                    + (6 * self.variants_ring_height),
                    "R1": self.variants_position
                    + (7 * self.variants_ring_space)
                    + (7 * self.variants_ring_height),
                },
                "layout": {
                    "opacity": 0.1,
                    "fillcolor": "gray",
                    "layer": "below",
                    "line": {"color": "gray", "width": 1},
                },
            },
            {
                # level 5
                "radius": {
                    "R0": self.variants_position
                    + (6 * self.variants_ring_space)
                    + (5 * self.variants_ring_height),
                    "R1": self.variants_position
                    + (6 * self.variants_ring_space)
                    + (6 * self.variants_ring_height),
                },
                "layout": {
                    "opacity": 0.1,
                    "fillcolor": "lightgrey",
                    "layer": "below",
                    "line": {"color": "lightgrey", "width": 1},
                },
            },
            {
                # level 4
                "radius": {
                    "R0": self.variants_position
                    + (5 * self.variants_ring_space)
                    + (4 * self.variants_ring_height),
                    "R1": self.variants_position
                    + (5 * self.variants_ring_space)
                    + (5 * self.variants_ring_height),
                },
                "layout": {
                    "opacity": 0.1,
                    "fillcolor": "lightgrey",
                    "layer": "below",
                    "line": {"color": "lightgrey", "width": 1},
                },
            },
            {
                # level 3
                "radius": {
                    "R0": self.variants_position
                    + (4 * self.variants_ring_space)
                    + (3 * self.variants_ring_height),
                    "R1": self.variants_position
                    + (4 * self.variants_ring_space)
                    + (4 * self.variants_ring_height),
                },
                "layout": {
                    "opacity": 0.1,
                    "fillcolor": "lightgrey",
                    "layer": "below",
                    "line": {"color": "lightgrey", "width": 1},
                },
            },
            {
                # level 2
                "radius": {
                    "R0": self.variants_position
                    + (3 * self.variants_ring_space)
                    + (2 * self.variants_ring_height),
                    "R1": self.variants_position
                    + (3 * self.variants_ring_space)
                    + (3 * self.variants_ring_height),
                },
                "layout": {
                    "opacity": 0.1,
                    "fillcolor": "white",
                    "layer": "below",
                    "line": {"color": "white", "width": 1},
                },
            },
            {
                # level 1
                "radius": {
                    "R0": self.variants_position
                    + (2 * self.variants_ring_space)
                    + (1 * self.variants_ring_height),
                    "R1": self.variants_position
                    + (2 * self.variants_ring_space)
                    + (2 * self.variants_ring_height),
                },
                "layout": {
                    "opacity": 0.1,
                    "fillcolor": "lightgrey",
                    "layer": "below",
                    "line": {"color": "lightgrey", "width": 1},
                },
            },
            {
                # level 0
                "radius": {
                    "R0": self.variants_position
                    + (1 * self.variants_ring_space)
                    + (0 * self.variants_ring_height),
                    "R1": self.variants_position
                    + (1 * self.variants_ring_space)
                    + (1 * self.variants_ring_height),
                },
                "layout": {
                    "opacity": 0.1,
                    "fillcolor": "lightgrey",
                    "layer": "below",
                    "line": {"color": "lightgrey", "width": 1},
                },
            },
        ]

    def merge_options(self):
        dico = {}
        dico["patch"] = {}
        # ideo = Ideogram()
        dico["patch"]["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {
                "orient": "columns",
                "data": {
                    "chr_name": ["chr1", "chr2", "chr3"],
                    "chr_size": [249250621, 243199373, 198022430],
                    "chr_label": ["chr1", "chr2", "chr3"],
                    "chr_color": ["pink", "rosybrown", "firebrick"],
                },
            },
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
        # dico["ring"] = self.ringval

        #  self.options[attr] = self.attr
        return dico, self.ringval
        # Loopable as fuck

