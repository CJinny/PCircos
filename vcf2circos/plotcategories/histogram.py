from vcf2circos.plotcategories.plotconfig import Plotconfig
from os.path import join as osj
import pandas as pd
import os


class Histogram_(Plotconfig):
    """
    It need to create one histogram for each SV event FOR EACH SV height (from 0 copy number to 5), which will create the grey band between color dor
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
        self.config_ring = config_ring
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]
        self.rangescale = rangescale
        # corresponding to SNV InDel height 7th ring (after 0 to 5 copy number height)
        self.radius = {
            "R0": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 1) * self.variants_ring_height),
            "R1": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 2) * self.variants_ring_height),
        }
        print("#Range", self.rangescale)
        # self.radius = {"R0": 0.90, "R1": 0.92}
        self.file = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }
        self.hovertextformat = " \"<b>{}:{}-{}</b><br>ClinGen informations:<br>{}\".format(a[i,0], a[i,1], a[i,2], a[i,5].replace(';', '<br>').replace('%2C', '<br>   ')) "
        self.trace = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {"size": 5, "symbol": 0, "color": "gray", "opacity": 1},
        }
        self.layout = {
            "type": "path",
            "opacity": 1,
            "fillcolor": "gray",
            "line": {"color": "gray", "width": 5},
        }

    def data_histogram_variants(self) -> dict:
        # for each sv event regarding copy number
        file = self.file
        data = {
            "chr_name": [],
            "start": [],
            "end": [],
            "val": [],
            "color": [],
            "info": [],
        }
        for val in self.data["Variants"]:
            [data["chr_name"].append(chroms) for chroms in self.data["Chromosomes"]]
            data["start"].append(int(val["SV_start"]))
            data["end"].append(int(val["SV_end"]))
            data["val"].append(1)
            data["color"].append("grey")
            data["info"].append(
                ";".join([str(key) + "=" + str(value) for key, value in val.items()])
            )
        file["dataframe"]["data"] = data
        return file

    def merge_options(self) -> list:
        for cn in self.rangescale:
            data = {}
            data["show"] = self.show
            data["customfillcolor"] = "False"
            data["file"] = self.data_histogram_variants()
            data["sortbycolor"] = "False"
            data["colorcolumn"] = 4
            radius = (
                (
                    self.variants_position
                    + (cn * self.variants_ring_space)
                    + ((cn - 1) * self.variants_ring_height)
                )
                + (
                    self.variants_position
                    + (cn * self.variants_ring_space)
                    + (cn * self.variants_ring_height)
                )
            ) / 2
            data["radius"] = {
                "R0": radius,
                "R1": radius,
            }
            data["hovertextformat"] = self.hovertextformat
            data["trace"] = self.trace
            data["layout"] = self.layout
        return data
