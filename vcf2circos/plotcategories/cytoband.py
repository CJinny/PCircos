from vcf2circos.plotcategories.plotconfig import Plotconfig
from os.path import join as osj
from itertools import repeat
import pandas as pd
import os


class Cytoband(Plotconfig):
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
        # Cytoband params
        # Need creation of dict in options attribute, regarding data input
        assert os.path.exists(
            osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                "cytoband_" + self.options["Assembly"] + "_chr_infos.txt.gz",
            )
        )
        self.cytoband_conf = pd.read_csv(
            osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                "cytoband_" + self.options["Assembly"] + "_chr_infos.txt.gz",
            ),
            sep="\t",
            header=0,
            compression="infer",
        )
        self.colorcolumn = 3
        self.sortbycolor = "True"
        self.hovertextformat = ' "<b>{}</b>".format(a[i,0])'
        self.trace = {
            "uid": "cytoband",
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {"size": 1, "symbol": 0, "color": None, "opacity": 1},  # 8
        }
        self.layout = {
            "type": "path",
            "layer": "below",
            "opacity": 1.0,
            "line": {"color": None, "width": 0},
        }

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def data_cytoband(self, chr_bnd):
        """
        histo band which will contains cytoband annotations, do not forget chromosomes carrying BND
        """
        chr_list = self.data["Chromosomes"]
        chr_list.extend([chrs_rec for chrs_rec in chr_bnd])
        chr_list = list(set(chr_list))
        tmp = self.cytoband_conf.loc[self.cytoband_conf["chr_name"].isin(chr_list)]
        data = {
            "chr_name": tmp["chr_name"].to_list(),
            "start": tmp["start"].tolist(),
            "end": tmp["end"].tolist(),
            "band_color": tmp["band_color"].tolist(),
            "band": tmp["band"].tolist(),
        }
        return data

    def merge_options(self, chr_bnd):
        """
        cytoband then histogram in list guess for the color
        """

        dico = {}
        dico
        dico["show"] = self.show
        dico["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": self.data_cytoband(chr_bnd)},
        }
        dico["sortbycolor"] = self.sortbycolor
        dico["colorcolumn"] = self.colorcolumn
        dico["hovertextformat"] = self.hovertextformat
        dico["trace"] = self.trace
        dico["layout"] = self.layout
        dico["trace"]["marker"]["color"] = list(
            repeat(
                self.options["Color"]["CYTOBAND"],
                len(dico["file"]["dataframe"]["data"]["chr_name"]),
            )
        )
        dico["layout"]["line"]["color"] = list(
            repeat(
                self.options["Color"]["CYTOBAND"],
                len(dico["file"]["dataframe"]["data"]["chr_name"]),
            )
        )
        return dico
