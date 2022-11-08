from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import variants_color, check_data_plot
from collections import OrderedDict
from os.path import join as osj
import numpy as np
from pprint import pprint


class Scatter_(Plotconfig):
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
        self.hovertextformat = ' "<b>{}:{}</b> | {} > {}<br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,3], a[i,4], a[i,5], a[i,7])'
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

    def adapt_data(self, histo_data: list) -> list:
        final = []
        # list of dico from histo class
        for dico in histo_data:
            # key val in each dico
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
                final.append(
                    [
                        od,
                        dico["radius"],
                        dico["trace"]["uid"],
                    ]
                )
        return final

    def merge_options(self, histo_data):
        final = []
        data_list_list = self.adapt_data(histo_data)

        # return self.scatter_cnv_level(
        #    data_list_list[0], data_list_list[1], data_list_list[2]
        # )
        for data, radius, level in data_list_list:
            final.append(self.scatter_cnv_level(data, radius, level))
        # CHECK
        # for dico_data in final:
        #    check_data_plot(dico_data)
        return final

    def scatter_cnv_level(self, data, radius, level):
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
                "symbol": data["symbol"],
                "color": data["color"],
                "opacity": 1,
            },
            "uid": level,
        }
        check_data_plot(data)
        return d
