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
        config_ring,
        rangescale,
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
        # In order to auto generate rings
        self.rangescale = rangescale
        self.config_ring = config_ring
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]

    def create_ring(self):
        rings_list = []
        # Copy number level
        for i, coeff in enumerate(self.rangescale):
            rings_list.append(
                {
                    "radius": {
                        "R0": coeff,
                        "R1": coeff + self.config_ring["height"],
                    },
                    "layout": {
                        "opacity": 0.1,
                        "fillcolor": "lightgray",
                        "layer": "below",
                        "line": {"color": "lightgray", "width": 1},
                    },
                },
            )
        rings_list[1]["layout"]["fillcolor"] = "white"
        rings_list[1]["layout"]["line"]["color"] = "white"

        rings_list[6]["layout"]["fillcolor"] = "grey"
        rings_list[6]["layout"]["line"]["color"] = "grey"
        return rings_list
