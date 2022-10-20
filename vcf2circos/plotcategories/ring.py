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
        min_l,
        max_l,
        nrings,
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
        # number of ring, min, max
        self.min_l = min_l
        self.max_l = max_l
        self.nrings = nrings
        self.config_ring = config_ring
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]

        coeff = np.linspace(self.min_l, self.max_l, num=self.nrings)
        for coeff in range_scale
            self.ringval = [
                {
                    "radius": {
                        "R0": self.variants_position
                        + (coeff * self.variants_ring_space)
                        + ((coeff - 1) * self.variants_ring_height),
                        "R1": self.variants_position
                        + (coeff * self.variants_ring_space)
                        + (coeff * self.variants_ring_height),
                    },
                    "layout": {
                        "opacity": 0.1,
                        "fillcolor": "gray",
                        "layer": "below",
                        "line": {"color": "gray", "width": 1},
                    },
                },
            ]

        # self.ringval = [
        #    {
        #        # SNV
        #        "radius": {
        #            "R0": self.variants_position
        #            + (7 * self.variants_ring_space)
        #            + (6 * self.variants_ring_height),
        #            "R1": self.variants_position
        #            + (7 * self.variants_ring_space)
        #            + (7 * self.variants_ring_height),
        #        },
        #        "layout": {
        #            "opacity": 0.1,
        #            "fillcolor": "gray",
        #            "layer": "below",
        #            "line": {"color": "gray", "width": 1},
        #        },
        #    },
        #    {
        #        # level 5
        #        "radius": {
        #            "R0": self.variants_position
        #            + (6 * self.variants_ring_space)
        #            + (5 * self.variants_ring_height),
        #            "R1": self.variants_position
        #            + (6 * self.variants_ring_space)
        #            + (6 * self.variants_ring_height),
        #        },
        #        "layout": {
        #            "opacity": 0.1,
        #            "fillcolor": "lightgrey",
        #            "layer": "below",
        #            "line": {"color": "lightgrey", "width": 1},
        #        },
        #    },
        #    {
        #        # level 4
        #        "radius": {
        #            "R0": self.variants_position
        #            + (5 * self.variants_ring_space)
        #            + (4 * self.variants_ring_height),
        #            "R1": self.variants_position
        #            + (5 * self.variants_ring_space)
        #            + (5 * self.variants_ring_height),
        #        },
        #        "layout": {
        #            "opacity": 0.1,
        #            "fillcolor": "lightgrey",
        #            "layer": "below",
        #            "line": {"color": "lightgrey", "width": 1},
        #        },
        #    },
        #    {
        #        # level 3
        #        "radius": {
        #            "R0": self.variants_position
        #            + (4 * self.variants_ring_space)
        #            + (3 * self.variants_ring_height),
        #            "R1": self.variants_position
        #            + (4 * self.variants_ring_space)
        #            + (4 * self.variants_ring_height),
        #        },
        #        "layout": {
        #            "opacity": 0.1,
        #            "fillcolor": "lightgrey",
        #            "layer": "below",
        #            "line": {"color": "lightgrey", "width": 1},
        #        },
        #    },
        #    {
        #        # level 2
        #        "radius": {
        #            "R0": self.variants_position
        #            + (3 * self.variants_ring_space)
        #            + (2 * self.variants_ring_height),
        #            "R1": self.variants_position
        #            + (3 * self.variants_ring_space)
        #            + (3 * self.variants_ring_height),
        #        },
        #        "layout": {
        #            "opacity": 0.1,
        #            "fillcolor": "white",
        #            "layer": "below",
        #            "line": {"color": "white", "width": 1},
        #        },
        #    },
        #    {
        #        # level 1
        #        "radius": {
        #            "R0": self.variants_position
        #            + (2 * self.variants_ring_space)
        #            + (1 * self.variants_ring_height),
        #            "R1": self.variants_position
        #            + (2 * self.variants_ring_space)
        #            + (2 * self.variants_ring_height),
        #        },
        #        "layout": {
        #            "opacity": 0.1,
        #            "fillcolor": "lightgrey",
        #            "layer": "below",
        #            "line": {"color": "lightgrey", "width": 1},
        #        },
        #    },
        #    {
        #        # level 0
        #        "radius": {
        #            "R0": self.variants_position
        #            + (1 * self.variants_ring_space)
        #            + (0 * self.variants_ring_height),
        #            "R1": self.variants_position
        #            + (1 * self.variants_ring_space)
        #            + (1 * self.variants_ring_height),
        #        },
        #        "layout": {
        #            "opacity": 0.1,
        #            "fillcolor": "lightgrey",
        #            "layer": "below",
        #            "line": {"color": "lightgrey", "width": 1},
        #        },
        #    },
        # ]
