from pprint import pprint
from typing import Generator
from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import timeit, generate_hovertext_var, chr_valid

from os.path import join as osj
import pandas as pd
from itertools import chain, repeat
from collections import OrderedDict, Counter


class Histogram_(Plotconfig):
    """
    It need to create one histogram for each SV event FOR EACH SV height (from 0 copy number to 5), which will create the grey band between color dor
    """

    def __init__(self, plotconfig):
        self.plotconfig = plotconfig
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
        # print("#Range", self.rangescale)
        self.hovertextformat = ' "<b>{}:{}-{}</b><br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,6], a[i,8])'

        # self.hovertextformat = ""
        self.trace = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 5,
                "symbol": 0,
                "color": self.colors["INTERMEDIATE"],
                "opacity": 0.1,
            },
        }
        self.layout = {
            "type": "path",
            "opacity": 1,
            "fillcolor": self.colors["INTERMEDIATE"],
            "line": {"color": self.colors["INTERMEDIATE"], "width": 5},
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
            "hovertextformat": ' "<b>{}:{}-{}<br>{}</b>".format(a[i,0], a[i,1], a[i,2], a[i,5])',
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
        self.only_overlapping_ = self.options["Genes"]["only_snv_in_sv_genes"]

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def dict_to_str(self, info_field: list) -> Generator:
        for info_dict in info_field:
            yield ";".join(
                [str(key) + "=" + str(value) for key, value in info_dict.items()]
            )

    def only_snv_indels_in_sv(self, data):
        df_data = pd.DataFrame.from_dict(data).astype(
            {
                "chr_name": str,
                "start": int,
                "end": int,
                "val": int,
                "ref": str,
                "alt": str,
                "type": str,
                "color": str,
                "hovertext": str,
                "symbol": int,
                "genes": str,
                "exons": str,
            }
        )
        snv_indels_df = df_data[df_data["type"].isin(["SNV", "INDEL", "OTHER"])]
        sv_df = df_data[~df_data["type"].isin(["SNV", "INDEL", "OTHER"])]
        # pd.set_option("display.max_columns", None)
        # pd.set_option("display.width", None)
        # pd.set_option("display.max_colwidth", -1)
        # SNV / INDEL
        snv_indel_not_overlapp = []
        snv_indel_overlapp = []
        for i, si in snv_indels_df.iterrows():
            # WHOLE
            for j, v in sv_df.iterrows():
                if (
                    si["chr_name"] == v["chr_name"]
                    and si["start"] > v["start"]
                    and si["start"] < v["end"]
                ):
                    snv_indel_overlapp.append(
                        [si["chr_name"], si["start"], si["ref"], si["alt"]]
                    )
                # else:
                #    if (
                #        not [si["chr_name"], si["start"], si["end"]]
                #        in snv_indel_not_overlapp
                #    ):
                #        snv_indel_not_overlapp.append(
                #            [si["chr_name"], si["start"], si["end"]]
                #        )
        # print(snv_indel_not_overlapp)
        return snv_indel_overlapp

    def adapt_data(self, cn: int) -> dict:
        d_file = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": None},
        }
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

        df_data = self.df_data.loc[self.df_data["CopyNumber"] == cn]
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
            try:
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
            except IndexError:
                print("ERROR ", items)
                exit()

        data["chr_name"].extend(df_data["Chromosomes"].to_list())
        data["start"].extend(start)
        data["end"].extend(stop)
        data["val"].extend(list(repeat(2, len(df_data.index))))
        data["ref"].extend(ref)
        data["alt"].extend(alt)
        data["type"].extend(df_data["Variants_type"].to_list())
        data["color"].extend(
            list(repeat(self.colors["INTERMEDIATE"], len(df_data.index)))
        )
        # data["hovertext"].extend(list(itertools.repeat("", len(df_data.index))))
        data["hovertext"].extend(
            list(
                generate_hovertext_var(
                    df_data["Variants"],
                    full_annot=20,
                    true_annot=self.options["Variants"]["annotations"]["fields"],
                )
            )
        )
        # print(data["hovertext"])
        # data["hovertext"].extend(
        #    [
        #        "Genes ("
        #        + str(len(record.split(",")))
        #        + "): "
        #        + ",".join(record.split(",")[:5])
        #        for record in df_data["Genes"].to_list()
        #    ]
        # )
        data["symbol"].extend(list(repeat(0, len(df_data.index))))
        data["genes"].extend(df_data["Genes"].to_list())
        data["exons"].extend(list(repeat("", len(df_data.index))))
        # data["info"].extend(list(self.dict_to_str(df_data["Variants"].to_list())))

        d_file["dataframe"]["data"] = data
        return d_file

    def histo_cnv_level(self, cn: int) -> dict:
        d = {}
        d_file = self.adapt_data(cn)
        d["show"] = "True"
        d["customfillcolor"] = "False"
        d["file"] = d_file
        d["sortbycolor"] = "False"
        d["colorcolumn"] = 7
        try:
            radius = (
                self.rangescale[cn]
                + self.rangescale[cn]
                + self.options["Variants"]["rings"]["height"]
            ) / 2
        except TypeError:
            print(cn)
            print(d_file)
            exit()
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
                "color": self.colors["INTERMEDIATE"],
                "opacity": 0.1,
            },
            "uid": "cnv_scatter_level_" + str(cn),
        }
        d["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.1,
            "fillcolor": "red",
            "line": {"color": self.colors["INTERMEDIATE"], "width": 5},
        }
        return d

    def cytoband_tile(self, cytoband_data):
        dico_cyto = self.cytoband_data.copy()
        cyto = {}
        cyto["chr_name"] = cytoband_data["chr_name"]
        cyto["start"] = cytoband_data["start"]
        cyto["end"] = cytoband_data["end"]
        # Remember to have val column in data otherwise it leads to crash]
        cyto["val"] = list(repeat(1, len(cytoband_data["chr_name"])))
        cyto["band_color"] = list(
            repeat(self.colors["CYTOBAND"], len(cytoband_data["chr_name"]))
        )
        cyto["band"] = cytoband_data["band"]
        # Cytoband tiles 3  need fill data
        dico_cyto["file"]["dataframe"]["data"] = cyto

        dico_cyto["layout"]["line"]["color"] = cyto["band_color"]
        dico_cyto["trace"]["marker"]["color"] = cyto["band_color"]
        return dico_cyto

    def merge_options(self, cytoband_data: dict) -> list:
        """
        func handle math and geometry need to take data in specific order
        chr start stop val OTHERWIS TROUBLE
        """
        # exit()

        whole_cn = []
        # Histo_cnv_level
        for cn in list(set(self.data["CopyNumber"])):
            res = self.histo_cnv_level(cn)
            whole_cn.append(res)

        if self.only_overlapping_:
            print("#[INFO] SNV / indels Overlapping SV only")
            whole_var = {
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
            for dic in whole_cn:
                for key, val in dic["file"]["dataframe"]["data"].items():
                    whole_var[key].extend(val)
            snv_indel_overlapp = self.only_snv_indels_in_sv(whole_var)
            # print(snv_indel_overlapp)
            for dico in whole_cn:
                # print(dico["trace"]["uid"])
                if dico["trace"]["uid"] == "cnv_scatter_level_6":
                    df_ = pd.DataFrame.from_dict(dico["file"]["dataframe"]["data"])
                    # if at least one snv indel overlap a sv
                    if snv_indel_overlapp:
                        for wr in snv_indel_overlapp:
                            df_ = df_.loc[
                                (df_["chr_name"] == wr[0])
                                & (df_["start"] == wr[1])
                                & (df_["ref"] == wr[2])
                                & (df_["alt"] == wr[3])
                            ]
                            dico["file"]["dataframe"]["data"] = df_.to_dict("list")
                    else:
                        dico["file"]["dataframe"]["data"] = {
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
                # remaining_var.append(self.remove_snv_(snv_indel_overlapp))
        # Extra
        whole_cn.extend(self.generate_extra_plots_from_df())

        # Genes plots
        whole_cn.append(self.histo_genes())

        # cytoband tiles
        whole_cn.append(self.cytoband_tile(cytoband_data))
        return whole_cn

    def remove_snv_(self, list_to_remove):
        index_keep = []
        index_to_rm = []
        for vars in list_to_remove:
            # print(vars)
            for i, record in enumerate(self.data["Record"]):
                if (
                    record.CHROM == vars[0]
                    and record.POS == vars[1]
                    and record.REF == vars[2]
                    and str(record.ALT[0]) == vars[3]
                ):
                    index_keep.append(i)
        for j, rd in enumerate(self.data["Record"]):
            # print(rd.var_type)
            # print(rd)
            # print("\n")
            if rd.var_type == "snp" or rd.var_type == "indel":
                # print(rd)
                # print(j)
                if j not in index_keep:
                    # print("index to remove: " + str(j))
                    # print(record)
                    index_to_rm.append(j)

        self.df_data.drop(index=index_to_rm, inplace=True)
        self.data = self.df_data.to_dict("list")
        return self.data

    def process_gene_list(self, genes_list: list) -> Generator:
        for record in genes_list:
            if record:
                yield record.split(",")

    def histo_genes(self) -> dict:
        data = {}
        dico = {}
        # remove empty gene, df_data attribute of class basic data from plot config Parents class
        # gene_list = list(filter(lambda x: x != "", self.df_data["Genes"]))
        gene_list = list(
            set(
                list(
                    map(
                        str,
                        chain.from_iterable(
                            list(self.process_gene_list(self.df_data["Genes"]))
                        ),
                    )
                )
            )
        )
        # print(*self.df_genes.columns)
        # print(self.df_genes.head())
        ## select genes in or batch of variations (from refeseq assembly)
        df_filter = self.df_genes.loc[self.df_genes["gene"].isin(gene_list)]
        # print(*gene_list)
        # print(self.df_data["Genes"].head())

        # Set color
        for fields in df_filter.columns:
            if fields != "transcript":
                if fields == "color":
                    # data[fields] = list(self.morbid_genes(df_filter["gene"]))
                    data[fields] = list(
                        repeat(self.colors["GENES"], len(df_filter.index))
                    )
                else:
                    data[fields] = df_filter[fields].to_list()
        # pprint(data, sort_dicts=False)
        dico["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }
        dico["show"] = self.show
        dico["colorcolumn"] = 4
        dico["radius"] = {"R0": 0.96, "R1": 0.96}
        dico[
            "hovertextformat"
        ] = ' "<b>{}:{}-{}<br>Gene: {}</b><br>".format(a[i,0], a[i,1], a[i,2], a[i,5])'
        dico["trace"] = {
            "uid": "genes",
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 3,
                "symbol": 0,
                "color": data["color"],
                "opacity": 1,
            },
        }
        dico["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.2,
            "line": {"color": data["color"], "width": 3},
        }
        return dico

    def genes_omim_morbid(self):
        """ "
        If it's a morbid gene it will be colored in red in circos gene level
        done in static file in genes.<assembly>
        """
        pass

    def extract_start_stop_ref_alt(
        self, record: list, info_field: list, variant_type: list
    ) -> Generator:
        # infer type of var could be done before
        for i, info_dict in enumerate(info_field):
            if variant_type[i] not in ["OTHER", "SNV", "INDEL"]:
                if "SV_start" in record[i].INFO and "SV_end" in record[i].INFO:
                    yield (
                        int(info_dict.get("SV_start").split("|")[0]),
                        int(info_dict.get("SV_end").split("|")[0]),
                        record[i].REF,
                        record[i].ALT,
                    )
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

    def generate_extra_plots_from_df(self):
        extras = []
        if "gc" in self.options["Extra"]:
            # self.gcplus = pd.DataFrame(osj(self.options["static"], #"histogram_pos_chr"))
            for gc_ in ["histogram_pos_chr.txt", "histogram_neg_chr.txt"]:
                gc_dict = {
                    "show": "True",
                    "customfillcolor": "False",
                    "file": {
                        "path": osj(self.options["Static"], gc_),
                        "header": "infer",
                        "sep": "\t",
                    },
                    "sortbycolor": "False",
                    "colorcolumn": "None",
                    "radius": {"R0": 0.90, "R1": 0.94},
                    "hovertextformat": ' "Chromosome: {}<br>Start: {}<br>End: {}<br>LogFC:{}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ',
                    "trace": {
                        "hoverinfo": "text",
                        "mode": "markers",
                        "marker": {"size": 0, "opacity": 0},
                        "uid": "extra_gc",
                    },
                    "layout": {
                        "type": "path",
                        "opacity": 1,
                        "fillcolor": "blue",
                        "line": {"color": "blue", "width": 0},
                    },
                }
                extras.append(gc_dict)

        if "mappability" in self.options["Extra"]:
            data = pd.read_csv(
                osj(self.options["Static"], "dukeExcludeRegions.csv"),
                header=0,
                sep="\t",
            )
            data = data.loc[data["chr_name"].isin(chr_valid())]
            data["val"] = 2
            data["color"] = "red"
            data["ref"] = ""
            data["alt"] = ""
            # data["infos"] = ""
            # data["hovertext"] = ""
            # data["infos_dict"] = ""
            mappa_dict = {
                "show": "True",
                "customfillcolor": "False",
                "file": {
                    "path": "",
                    "header": "infer",
                    "sep": "\t",
                    "dataframe": {"orient": "columns", "data": data.to_dict("list")},
                },
                "sortbycolor": "False",
                "colorcolumn": 7,
                "radius": {"R0": 0.80, "R1": 0.84},
                "hovertextformat": ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Type:{}".format(a[i,0], a[i,1], a[i,2], a[i,6]) ',
                "trace": {
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {"size": 0, "opacity": 0},
                    "uid": "extra_mappability",
                },
                "layout": {
                    "type": "path",
                    "opacity": 1,
                    "fillcolor": "black",
                    "line": {"color": "black", "width": 0},
                },
            }
            extras.append(mappa_dict)
            print(mappa_dict)
            # exit()
        return extras
