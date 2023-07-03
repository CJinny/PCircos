from typing import Generator
from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import check_data_plot
from collections import OrderedDict
from os.path import join as osj
import numpy as np
from pprint import pprint
from itertools import repeat, chain


class Scatter_(Plotconfig):
    def __init__(self, plotconfig, data_histo):
        self.plotconfig = plotconfig
        self.hovertextformat = ' "<b>{}:{}</b> | {} > {}<br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,3], a[i,4], a[i,5], a[i,7])'
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
        self.data_histo = data_histo

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def adapt_genes(self, dico: dict) -> dict:
        tmp = {}
        od = OrderedDict()
        tmp["chr_name"] = list(
            chain(
                *zip(
                    dico["chr_name"],
                    dico["chr_name"],
                )
            )
        )
        tmp["start"] = []
        for s, e in zip(
            dico["start"],
            dico["end"],
        ):
            tmp["start"].append(s)
            tmp["start"].append(e)
        tmp["val"] = list(chain(*zip(dico["val"], dico["val"])))
        tmp["color"] = [var for var in dico["color"]]
        for key, val in tmp.items():
            if key != "color":
                od[key] = val
        for key, val in dico.items():
            if key not in od.keys() and key not in [
                "genes",
                "end",
                "exons",
                "hovertext",
                "symbol",
            ]:
                od[key] = list(chain(*zip(val, val)))
        od["color"] = list(
            self.morbid_genes(list(chain(*zip(dico["gene"], dico["gene"]))))
        )
        od["infos"] = list(repeat("", len(od["chr_name"])))
        od["hovertext"] = list(repeat("", len(od["chr_name"])))
        return od

    def val_data_col(self, dico):
        for key, val in dico["file"]["dataframe"]["data"].items():
            if key not in ["start", "end", "color", "genes", "exons"]:
                yield key, list(np.repeat(val, 2))

    def adapt_data(self) -> list:
        final = []
        # list of dico from histo class
        for dico in self.data_histo:
            if dico["trace"]["uid"].startswith("cnv_"):
                tmp = {}
                od = OrderedDict()
                tmp["chr_name"] = list(
                    chain(
                        *zip(
                            dico["file"]["dataframe"]["data"]["chr_name"],
                            dico["file"]["dataframe"]["data"]["chr_name"],
                        )
                    )
                )
                tmp["start"] = []
                for s, e in zip(
                    dico["file"]["dataframe"]["data"]["start"],
                    dico["file"]["dataframe"]["data"]["end"],
                ):
                    tmp["start"].append(s)
                    tmp["start"].append(e)
                tmp["color"] = [
                    self.options["Color"][var]
                    for var in dico["file"]["dataframe"]["data"]["type"]
                ]
                for key, val in tmp.items():
                    if key != "color":
                        od[key] = val

                for key, val in dico["file"]["dataframe"]["data"].items():
                    if key not in od.keys() and key not in [
                        "genes",
                        "end",
                        "exons",
                        "hovertext",
                        "symbol",
                    ]:
                        od[key] = list(chain(*zip(val, val)))
                od["color"] = list(chain(*zip(tmp["color"], tmp["color"])))
                od["hovertext"] = list(
                    chain(
                        *zip(
                            dico["file"]["dataframe"]["data"]["hovertext"],
                            dico["file"]["dataframe"]["data"]["hovertext"],
                        )
                    )
                )
                od["symbol"] = list(
                    chain(
                        *zip(
                            dico["file"]["dataframe"]["data"]["symbol"],
                            dico["file"]["dataframe"]["data"]["symbol"],
                        )
                    )
                )
                final.append([od, dico["radius"], dico["trace"]["uid"]])
            # Becarefull to not hoverride histogram data
            elif dico["trace"]["uid"] == "genes":
                tmp = self.adapt_genes(dico["file"]["dataframe"]["data"])
                tmp["color"] = list(self.morbid_genes(tmp["gene"]))
                final.append(
                    [
                        tmp,
                        dico["radius"],
                        dico["trace"]["uid"],
                    ]
                )
        return final

    def morbid_genes(self, genes: list) -> Generator:
        for g in genes:
            if g in self.df_morbid["genes"].to_list():
                yield self.options["Color"]["MORBID_GENES"]
            else:
                yield self.options["Color"]["GENES"]

    def merge_options(self):
        final = []
        data_list_list = self.adapt_data()
        for data, radius, level in data_list_list:
            # if no more mutations in copy number level remove dict
            if data["chr_name"]:
                if level == "genes":
                    colorcolumn = 3
                    symbol = 0
                    hovertextformat = ' "<b>{}:{}<br>Gene: {}</b><br>{}".format(a[i,0], a[i,1], a[i,4], a[i,6])'
                    final.append(
                        self.scatter_cnv_level(
                            data, radius, level, hovertextformat, symbol, colorcolumn
                        )
                    )
                else:
                    final.append(self.scatter_cnv_level(data, radius, level))
        for dico in final:
            check_data_plot(dico["file"]["dataframe"]["data"])
        return final

    def scatter_cnv_level(
        self, data, radius, level, hovertextformat=None, symbol=None, colorcolumn=None
    ):
        if hovertextformat is None:
            hovertextformat = self.hovertextformat
        if symbol is None:
            if "symbol" in data.keys():
                symbol = data["symbol"]
            else:
                symbol = data["val"]
        if colorcolumn is None:
            colorcolumn = 6
        d = {}
        d["show"] = "True"
        d["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }

        d["sortbycolor"] = "False"
        d["colorcolumn"] = colorcolumn
        d["radius"] = radius
        d["hovertextformat"] = hovertextformat
        d["trace"] = {
            "hoverinfo": "text",
            "mode": "markers",
            "opacity": 1,
            "marker": {
                "size": 5,
                "symbol": symbol,
                "color": data["color"],
                "opacity": 1,
            },
            "uid": level,
        }
        d["layout"] = {"showlegend": "True"}
        d["name"] = d["trace"]["uid"]
        for key, items in d["file"]["dataframe"]["data"].items():
            if key == "ref" or key == "alt":
                for j, var in enumerate(items):
                    if len(var) > 15:
                        items[j] = var[:15] + "..."
        return d
