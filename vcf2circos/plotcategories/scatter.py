from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import variants_color, check_data_plot
from collections import OrderedDict
from os.path import join as osj
import numpy as np
from pprint import pprint
from itertools import repeat


class Scatter_(Plotconfig):
    def __init__(self, plotconfig):
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
        gene_scatter["hovertext"] = dico["gene"] + dico["gene"]
        print(gene_scatter.keys())
        return gene_scatter

    def adapt_data(self, histo_data: list) -> list:
        final = []
        # list of dico from histo class
        for dico in histo_data:
            # key val in each dico
            # if dico["trace"]["uid"] == "cnv_scatter_level_6":
            #    pprint(dico["file"]["dataframe"]["data"])
            #    exit()
            if dico["trace"]["uid"].startswith("cnv_"):
                tmp = {}
                for key, val in dico["file"]["dataframe"]["data"].items():
                    if key not in ["start", "end", "color", "genes", "exons"]:
                        tmp[key] = list(np.repeat(val, 2))
                        # utils color match en fonciton du type
                        #    tmp_se.append(val)
                tmp["start"] = []
                for s, e in zip(
                    dico["file"]["dataframe"]["data"]["start"],
                    dico["file"]["dataframe"]["data"]["end"],
                ):
                    tmp["start"].append(s)
                    tmp["start"].append(e)
                tmp["color"] = [variants_color[var] for var in tmp["type"]]
                od = OrderedDict()
                od["chr_name"] = tmp["chr_name"]
                od["start"] = tmp["start"]
                od["val"] = tmp["val"]
                od["ref"] = tmp["ref"]
                od["alt"] = tmp["alt"]
                od["type"] = tmp["type"]
                od["color"] = tmp["color"]
                for key, val in tmp.items():
                    if key not in od.keys():
                        od[key] = val
                final.append([od, dico["radius"], dico["trace"]["uid"]])
            elif dico["trace"]["uid"] == "genes":
                dico[
                    "hovertextformat"
                ] = ' "<b>{}:{}<br>Gene: {}</b><br>{}".format(a[i,0], a[i,1], a[i,4], a[i,6])'
                dico["file"]["dataframe"]["data"] = self.adapt_genes(
                    dico["file"]["dataframe"]["data"]
                )
                final.append(
                    [
                        dico["file"]["dataframe"]["data"],
                        dico["radius"],
                        dico["trace"]["uid"],
                    ]
                )

                # check_data_plot(self.adapt_genes(dico["file"]["dataframe"]["data"]))
        return final

    def merge_options(self, histo_data):
        final = []
        data_list_list = self.adapt_data(histo_data)

        # return self.scatter_cnv_level(
        #    data_list_list[0], data_list_list[1], data_list_list[2]
        # )
        for data, radius, level in data_list_list:
            # if no more mutations in copy number level remove dict
            if data["chr_name"]:
                final.append(self.scatter_cnv_level(data, radius, level))
        # CHECK
        # for dico_data in final:
        #    check_data_plot(dico_data)
        final[-1]["trace"]["marker"]["symbol"] = 0
        return final

    def scatter_cnv_level(self, data, radius, level):
        if data.get("symbol") is not None:
            symbol = data.get("symbol")
        else:
            symbol = data.get("val")
        d = {}
        d["show"] = "True"
        d["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }

        d["sortbycolor"] = "False"
        d["colorcolumn"] = 6
        d["radius"] = radius
        d["hovertextformat"] = self.hovertextformat
        d["trace"] = {
            "hoverinfo": "text",
            "mode": "markers",
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
