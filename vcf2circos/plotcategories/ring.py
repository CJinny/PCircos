from vcf2circos.plotcategories.plotconfig import Plotconfig
import numpy as np
from pprint import pprint

# space, space between ring in option.example.json
# height hauteur du ring
# positon position from center


# list_graph_type = ["majortick", "minortick", "ticklabel"]


class Ring(Plotconfig):
    def __init__(self, plotconfig, ring_upper_var=None):
        self.plotconfig = plotconfig
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]
        self.ring_upper_var = ring_upper_var

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    # ideogram cyto internal at 1.0
    def ring(self, R0, R1):
        return {
            "radius": {"R0": R0, "R1": R1},
            "layout": {
                "opacity": 0.1,
                "fillcolor": "lightgray",
                "layer": "below",
                "line": {"color": "lightgray", "width": 1},
            },
        }

    def create_ring(self):
        rings_list = []
        # Copy number level
        for i, coeff in enumerate(self.rangescale):
            rings_list.append(self.ring(coeff, coeff + self.variants_ring_height))

        rings_list[2]["layout"]["fillcolor"] = "white"
        rings_list[2]["layout"]["line"]["color"] = "white"

        rings_list[6]["layout"]["fillcolor"] = "grey"
        rings_list[6]["layout"]["line"]["color"] = "grey"
        if self.ring_upper_var is not None:
            R1 = 1.03
            R0 = R1 - self.variants_ring_height
            for field in self.ring_upper_var:
                R1 -= self.variants_ring_space + self.variants_ring_height
                R0 -= self.variants_ring_space + self.variants_ring_height
                rings_list.append(self.ring(R0, R1))
        return rings_list
