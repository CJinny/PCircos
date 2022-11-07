from vcf2circos.plotcategories.plotconfig import Plotconfig
from os.path import join as osj
import pandas as pd
import os
import numpy as np
from pprint import pprint


class Scatter_(Plotconfig):
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
        # self.df_data = pd.DataFrame.from_dict(self.data).astype(
        #    {
        #        "Chromosomes": str,
        #        "Genes": str,
        #        "Exons": str,
        #        "Variants": object,
        #        "Variants_type": str,
        #        "CopyNumber": int,
        #        "Color": str,
        #    }
        # )
        self.rangescale = rangescale
        self.config_ring = config_ring
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]
        # corresponding to SNV InDel height 7th ring (after 0 to 5 copy number height)
        self.radius = {
            "R0": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 1) * self.variants_ring_height),
            "R1": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 2) * self.variants_ring_height),
        }

    def scatter_variants(self):
        return self.adapt_data()

    def adapt_data(self):
        # no need end list compare to histogram cuz only dot
        # data = {
        #    "chr_name": [],
        #    "start": [],,
        #    "val": [],
        #    "ref": [],
        #    "alt": [],
        #    "type": [],
        #    "color": [],
        #    "hovertext": [],
        #    "symbol": [],
        #    "genes": [],
        #    "exons": [],
        # }
        #            "Chromosomes": [],
        #    "Genes": [],
        #    "Exons": [],
        #    "Record": [],
        #    "Variants": [],
        #    "Variants_type": [],
        #    "CopyNumber": [],
        #    "Color": []
        #
        # keys = ["Chromosome", "Genes", ]
        # data_f = {x:self.data[x] for x in keys}
        pass

    def merge_options(self, histo_data):
        final = []
        # list of dico from histo class
        for dico in histo_data:
            # key val in each dico
            if dico["trace"]["uid"].startswith("cnv_"):
                tmp = {}
                tmpe_se = []
                for key, val in dico["file"]["dataframe"]["data"].items():
                    if key not in ["start", "end"]:
                        tmp[key] = list(np.repeat(val, 2))
                    elif key == "color":
                        #utils color match en fonciton du type
                    #    tmp_se.append(val)
                tmp["start"] = []
                for s, e in zip(
                    dico["file"]["dataframe"]["data"]["start"],
                    dico["file"]["dataframe"]["data"]["end"],
                ):
                    tmp["start"].append(s)
                    tmp["start"].append(e)
                final.append(tmp)

        pprint(final)
        print(len(final[0]["chr_name"]))
        pass

    def scatter_cnv_level(self, cn):
        d = {}

        d["show"] = "True"
        d["customfillcolor"] = "False"
        d["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }

        d["sortbycolor"] = "False"
        d["colorcolumn"] = 7
        radius = (
            self.rangescale[cn]
            + self.rangescale[cn]
            + self.options["Variants"]["rings"]["height"]
        ) / 2
        d["radius"] = {
            "R0": radius,
            "R1": radius,
        }
        d["hovertextformat"] = self.hovertextformat
        d["trace"] = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 5,
                "symbol": d_file["dataframe"]["data"]["symbol"],
                "color": "gray",
                "opacity": 0.1,
            },
            "uid": "cnv_level_" + str(cn),
        }
        d["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.1,
            "fillcolor": "red",
            "line": {"color": "lightgray", "width": 5},
        }
        return d
