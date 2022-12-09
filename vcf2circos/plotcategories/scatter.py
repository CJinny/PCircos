from typing import Generator
from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import variants_color, check_data_plot
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
        gene_scatter = {}
        gene_scatter["chr_name"] = dico["chr_name"] + dico["chr_name"]
        gene_scatter["start"] = dico["start"] + dico["end"]
        gene_scatter["val"] = dico["val"] + dico["val"]
        gene_scatter["color"] = dico["color"] + dico["color"]
        gene_scatter["gene"] = dico["gene"] + dico["gene"]
        gene_scatter["infos"] = list(repeat("", len(gene_scatter["chr_name"])))

        # pprint(gene_scatter)
        # gene_scatter["hovertext"] = dico["gene"] + dico["gene"]
        gene_scatter["hovertext"] = list(repeat("", len(gene_scatter["chr_name"])))
        # print(gene_scatter.keys())
        return gene_scatter

    def val_data_col(self, dico):
        for key, val in dico["file"]["dataframe"]["data"].items():
            if key not in ["start", "end", "color", "genes", "exons"]:
                # tmp[key] = list(np.repeat(val, 2))
                # utils color match en fonciton du type
                #    tmp_se.append(val)

                yield key, list(np.repeat(val, 2))

    def adapt_data(self) -> list:
        final = []
        # list of dico from histo class
        for dico in self.data_histo:
            # key val in each dico
            # if dico["trace"]["uid"] == "cnv_scatter_level_6":
            #    pprint(dico["file"]["dataframe"]["data"])
            #    exit()
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
                # print(dico.keys())
                # print(dico["trace"]["uid"])
                tmp["color"] = [
                    variants_color[var]
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
                # od["hovertext"] = list(np.repeat("", len(od["color"])))
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
                # DEBUG
                # if dico["trace"]["uid"] == "cnv_scatter_level_4":
                #    print(od)
                # exit()
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

        # check_data_plot(self.adapt_genes(dico["file"]["dataframe"]["data"]))
        return final

    def morbid_genes(self, genes: list) -> Generator:
        for g in genes:
            if g in self.df_morbid["genes"].to_list():
                yield "red"
            else:
                yield "lightgray"

    def merge_options(self):
        final = []
        data_list_list = self.adapt_data()

        # return self.scatter_cnv_level(
        #    data_list_list[0], data_list_list[1], data_list_list[2]
        # )
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

        # CHECK
        # for dico_data in final:
        #    check_data_plot(dico_data)
        # final[-1]["trace"]["marker"]["symbol"] = 0

        # pprint(final, sort_dicts=False)
        # exit
        for dico in final:
            check_data_plot(dico["file"]["dataframe"]["data"])
            # print(dico["trace"]["uid"])
            # print(dico["file"]["dataframe"]["data"].keys())
            # print("\n")
        return final
        # if len(final) == 1:
        #    return final[0]
        # else:
        #    return final

    def scatter_cnv_level(
        self, data, radius, level, hovertextformat=None, symbol=None, colorcolumn=None
    ):
        # if data.get("symbol") is not None:
        #    symbol = data.get("symbol")
        # else:
        #    symbol = data.get("val")
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
        # check_data_plot(data)
        return d
