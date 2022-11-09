from vcf2circos.plotcategories.plotconfig import Plotconfig
import numpy as np

# space, space between ring in option.example.json
# height hauteur du ring
# positon position from center


list_graph_type = ["majortick", "minortick", "ticklabel"]


class Ring(Plotconfig):
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
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def create_ring(self):
        rings_list = []
        # Copy number level
        for i, coeff in enumerate(self.rangescale):
            rings_list.append(
                {
                    "radius": {"R0": coeff, "R1": coeff + self.config_ring["height"],},
                    "layout": {
                        "opacity": 0.1,
                        "fillcolor": "lightgray",
                        "layer": "below",
                        "line": {"color": "lightgray", "width": 1},
                    },
                },
            )
        rings_list[2]["layout"]["fillcolor"] = "white"
        rings_list[2]["layout"]["line"]["color"] = "white"

        rings_list[6]["layout"]["fillcolor"] = "grey"
        rings_list[6]["layout"]["line"]["color"] = "grey"
        return rings_list
