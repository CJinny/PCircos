from vcf2circos.plotcategories.plotconfig import Plotconfig
from os.path import join as osj
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
        self.file = file
        self.colorcolumn = 3
        self.sortbycolor = "True"
        self.hovertextformat = ' "<b>{}</b>".format(a[i,0])'
        self.trace = {
            "uid": "cytoband",
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {"size": 1, "symbol": 0, "color": "black", "opacity": 1,},  # 8
        }
        self.layout = {
            "type": "path",
            "layer": "below",
            "opacity": 1.0,
            "line": {"color": "black", "width": 0},
        }

    def data_cytoband(self):
        tmp = self.cytoband_conf.loc[
            self.cytoband_conf["chr_name"].isin(self.data["Chromosomes"])
        ]
        data = {
            "chr_name": tmp["chr_name"].to_list(),
            "start": tmp["start"].tolist(),
            "end": tmp["end"].tolist(),
            "band_color": tmp["band_color"].tolist(),
            "band": tmp["band"].tolist(),
        }
        return data

    def merge_options(self):
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
            "dataframe": {"orient": "columns", "data": self.data_cytoband()},
        }
        dico["sortbycolor"] = self.sortbycolor
        dico["colorcolumn"] = self.colorcolumn
        dico["hovertextformat"] = self.hovertextformat
        dico["trace"] = self.trace
        dico["layout"] = self.layout

        histo = []
        Cytoband_infos = {
            "show": "True",
            "file": {
                "path": "",
                "header": "infer",
                "sep": "\t",
                "dataframe": {"orient": "columns", "data": self.data_cytoband()},
            },
            "colorcolumn": 4,
            "radius": {"R0": 1, "R1": 1.1},
            "hovertextformat": " \"<b>{}:{}-{}<br>{}{}</b>\".format(a[i,0], a[i,1], a[i,2], a[i,0].replace('chr', ''), ''.join(a[i,5:]))",
            # "hovertextformat": " \"<b>{}</b>\".format(a[i,0])",
            "trace": {
                "uid": "cytoband_tile",
                "hoverinfo": "text",
                "mode": "markers",
                "marker": {
                    "size": 0,
                    "symbol": 0,  # 8
                    "color": self.data_cytoband()["band_color"],
                    "opacity": 0,
                },
                "hovertextformat": " \"<b>{}:{}-{}<br>{}{}</b>\".format(a[i,0], a[i,1], a[i,2], a[i,0].replace('chr', ''), ''.join(a[i,5:]))",
                "trace": {
                    "uid": "cytoband_tile",
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {
                        "size": 0,
                        "symbol": 0,  # 8
                        "color": self.data_cytoband()["band_color"],
                        "opacity": 0,
                    },
                },
            },
        }

        histo = [Cytoband_infos]
        return (dico, Cytoband_infos)
