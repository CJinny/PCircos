from vcf2circos.plotcategories.plotconfig import Plotconfig
from os.path import join as osj
import pandas as pd
import os


class Histogram_(Plotconfig):
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
        self.rangescale = rangescale
        self.radius = {
            "R0": self.options["Variants"]["rings"]["height"] * max(self.rangescale - 1)
            + 0.02,
            "R1": self.options["Variants"]["rings"]["height"] * max(self.rangescale)
            - 0.02,
        }
        print("#Range", self.rangescale)
        # self.radius = {"R0": 0.90, "R1": 0.92}
        self.file = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": self.data_histogram_variants()},
        }
        self.hovertextformat = " \"<b>{}:{}-{}</b><br>ClinGen informations:<br>{}\".format(a[i,0], a[i,1], a[i,2], a[i,5].replace(';', '<br>').replace('%2C', '<br>   ')) "
        self.trace = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {"size": 5, "symbol": 0, "color": "blue", "opacity": 1},
        }
        self.layout = {
            "type": "path",
            "opacity": 1,
            "fillcolor": "gray",
            "line": {"color": "gray", "width": 5},
        }

    def data_histogram_variants(self):
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
            data["color"].append("blue")
            data["info"].append(
                ";".join([str(key) + "=" + str(value) for key, value in val.items()])
            )
        return data

    def merge_options(self):
        data = {}
        data["show"] = self.show
        data["customfillcolor"] = "False"
        data["file"] = self.file
        data["sortbycolor"] = "False"
        data["colorcolumn"] = 4
        data["radius"] = self.radius
        data["hovertextformat"] = self.hovertextformat
        data["trace"] = self.trace
        data["layout"] = self.layout
        return data
