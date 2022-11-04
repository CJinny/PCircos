from pprint import pprint
from typing import Generator
from vcf2circos.plotcategories.plotconfig import Plotconfig
from os.path import join as osj
import pandas as pd
import os
import itertools
from collections import OrderedDict


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
        self.hovertextformat = ' "<b>{}:{}-{}</b><br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,6], a[i,8])'
        self.trace = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {"size": 5, "symbol": 0, "color": "gray", "opacity": 0.1},
        }
        self.layout = {
            "type": "path",
            "opacity": 1,
            "fillcolor": "gray",
            "line": {"color": "gray", "width": 5},
        }
        # TODO tile same as cytobandinfo in vfreader
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
        self.cytoband_data = {
            "show": "True",
            "file": {
                "path": "",
                "header": "infer",
                "sep": "\t",
                "dataframe": {"orient": "columns", "data": None},
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
                    "symbol": 0,
                    "color": None,
                    "opacity": 0,
                },  # 8
            },
            "layout": {
                "type": "path",
                "layer": "above",
                "opacity": 0,
                "line": {"color": None, "width": 0},
            },
        }

    def cytoband_histogram(self):
        pass

    def genes_histogram(self):
        pass

    def data_histogram_variants(self, cn) -> dict:
        # for each sv event regarding copy number
        # TODO get information in data as DATAFRAME
        dico = self.file.deepcopy()
        data = {
            "chr_name": [],
            "start": [],
            "end": [],
            "val": [],
            "ref": [],
            "alt": [],
            "type": [],
            "color": [],
            "hovertext": [],
            "symbol": [],
            "genes": [],
            "exons": [],
        }
        df_ = pd.DataFrame.from_dict(self.data).astype(
            {
                "Chromosomes": str,
                "Genes": str,
                "Exons": str,
                "Variants": object,
                "Variants_type": str,
                "CopyNumber": int,
                "Color": str,
            }
        )
        df_data = df_.loc[df_["CopyNumber"] == cn]
        start = []
        stop = []
        ref = []
        alt = []
        for items in list(
            self.extract_start_stop_ref_alt(
                df_data["Record"].to_list(),
                df_data["Variants"].to_list(),
                df_data["Variants_type"].to_list(),
            )
        ):
            # DEBUGG
            # print(*items)
            # for val in items:
            #    if isinstance(val, list):
            #        print(type(val[0]))
            #    else:
            #        print(type(val))
            # exit()
            start.append(items[0])
            stop.append(items[1])
            ref.append(items[2])
            alt.append(str(items[3][0]))
        data["chr_name"].extend(df_data["Chromosomes"].to_list())
        data["start"].extend(start)
        data["end"].extend(stop)
        data["val"].extend(list(itertools.repeat(2, len(df_data.index))))
        data["ref"].extend(ref)
        data["alt"].extend(alt)
        data["type"].extend(df_data["Variants_type"].to_list())
        data["color"].extend(list(itertools.repeat("grey", len(df_data.index))))
        # data["hovertext"].extend(list(itertools.repeat("", len(df_data.index))))
        data["hovertext"].extend(
            [
                "Genes ("
                + str(len(record.split(",")))
                + "): "
                + ",".join(record.split(",")[:5])
                for record in df_data["Genes"].to_list()
            ]
        )
        data["symbol"].extend(list(itertools.repeat(0, len(df_data.index))))
        data["genes"].extend(df_data["Genes"].to_list())
        data["exons"].extend(list(itertools.repeat("", len(df_data.index))))
        # data["info"].extend(list(self.dict_to_str(df_data["Variants"].to_list())))
        dico["dataframe"]["data"] = data
        return dico

    def dict_to_str(self, info_field: list) -> Generator:
        for info_dict in info_field:
            yield ";".join(
                [str(key) + "=" + str(value) for key, value in info_dict.items()]
            )

    def histo_cnv_level(self):
        for cn in list(set(self.data["CopyNumber"])):
            global_d = self.data_histogram_variants(cn)
            dico = {}
            dico["show"] = "True"
            dico["customfillcolor"] = "False"
            dico["file"] = global_d
            dico["sortbycolor"] = "False"
            dico["colorcolumn"] = 7
            radius = (
                self.rangescale[cn]
                + self.rangescale[cn]
                + self.options["Variants"]["rings"]["height"]
            ) / 2
            dico["radius"] = {
                "R0": radius,
                "R1": radius,
            }
            dico["hovertextformat"] = self.hovertextformat
            dico["trace"] = {
                "hoverinfo": "text",
                "mode": "markers",
                "marker": {
                    "size": 5,
                    "symbol": global_d["dataframe"]["data"["symbol"],
                    "color": "gray",
                    "opacity": 0.1,
                },
                "uid": "cnv_level_" + str(cn)
            }
            dico["layout"] = {
                "type": "path",
                "layer": "above",
                "opacity": 0.1,
                "fillcolor": "red",
                "line": {"color": "lightgray", "width": 5},
            }
            yield dico

    def merge_options(self, cytoband_data: dict) -> list:
        """
        func handle math and geometry need to take data in specific order
        chr start stop val OTHERWIS TROUBLE
        """
        histo_data = []

        cyto = {}
        cyto["chr_name"] = cytoband_data["chr_name"]
        cyto["start"] = cytoband_data["start"]
        cyto["end"] = cytoband_data["end"]
        # Remember to have val column in data otherwise it leads to crash]
        cyto["val"] = list(itertools.repeat(1, len(cytoband_data["chr_name"])))
        cyto["band_color"] = list(
            itertools.repeat("lightgray", len(cytoband_data["chr_name"]))
        )
        cyto["band"] = cytoband_data["band"]
        # Cytoband tiles 3  need fill data
        self.cytoband_data["file"]["dataframe"]["data"] = cyto

        self.cytoband_data["layout"]["line"]["color"] = cyto["band_color"]
        self.cytoband_data["trace"]["marker"]["color"] = cyto["band_color"]

        return list(self.histo_cnv_level())

        # def __call__(self):
        #    return pd.DataFrame.from_dict(self.data)

    def extract_start_stop_ref_alt(
        self, record: list, info_field: list, variant_type: list
    ) -> Generator:
        # infer type of var could be done before
        for i, info_dict in enumerate(info_field):
            if variant_type[i] != "OTHER":
                if "SV_start" in record[i].INFO and "SV_end" in record[i].INFO:
                    yield (int(info_dict.get("SV_start")), int(info_dict.get("SV_end")))
                elif "END" in record[i].INFO:
                    yield (
                        int(record[i].POS),
                        int(record[i].INFO["END"]),
                        record[i].REF,
                        record[i].ALT,
                    )
                elif "SVLEN" in record[i].INFO:
                    yield (
                        int(record[i].POS),
                        int(abs(record[i].INFO["SVLEN"][0])) + int(record[i].POS),
                        record[i].REF,
                        record[i].ALT,
                    )
                else:
                    print("Can't establish SV length, annotations missing EXIT")
                    exit()
            # SNVINDEL
            else:
                alternate = int(str(max([len(alt) for alt in record[i].ALT])))
                yield (
                    int(str(record[i].POS)),
                    int(str(record[i].POS)) + alternate,
                    record[i].REF,
                    record[i].ALT,
                )
