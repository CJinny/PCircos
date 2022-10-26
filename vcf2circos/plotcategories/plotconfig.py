#!/usr/local/bin/python3

"""
aim: Epurate module and func, handle parameters of circos plot
"""

import sys
import re
import subprocess
import json
import os
import pandas as pd

from vcf2circos.vcfreader import VcfReader
from os.path import join as osj
from tqdm import tqdm
import gzip


# TODO commons file with utils function maybe one more for globals


class Plotconfig(VcfReader):
    """
    Options regroup options passed in args in json file otherwise
    it will be a empty dict,
    All func based on vcf input
    """

    def __init__(
        self,
        filename: str,
        options: dict,
        show: bool,
        file: dict,
        radius: dict,
        sortbycolor: bool,
        colorcolumn: int,
        hovertextformat: dict,
        trace_car: dict,
        data: list,
        layout: dict,
    ):
        super().__init__(filename, options)
        self.default_options = json.load(
            open("../demo_data/options.general.json", "r",)
        )
        if not self.options.get("General", {}).get("title", None):
            self.options["General"]["title"] = os.path.basename(
                self.get_metadatas().get("filename", "myCircos")
            )
        self.show = self.cast_bool(show)
        self.file = file
        self.radius = radius
        self.sortbycolor = self.cast_bool(sortbycolor)
        self.colorcolumn = colorcolumn
        self.hovertextformat = hovertextformat
        self.trace_car = trace_car
        self.data = self.process_vcf()
        self.layout = layout
        self.refgene_genes = osj(
            self.options["Static"],
            "Assembly",
            self.options["Assembly"],
            "genes." + self.options["Assembly"] + ".txt.gz",
        )
        self.refgene_exons = osj(
            self.options["Static"],
            "Assembly",
            self.options["Assembly"],
            "exons." + self.options["Assembly"] + ".txt.gz",
        )

    @staticmethod
    def cast_bool(value: bool) -> str:
        if value:
            return "True"
        else:
            return "False"

    def get_file_features():
        pass

    def get_snvindels_overlapping_sv(self):
        pass

    def vcf_options_default(self):
        pass

    def get_json(self) -> dict:
        """
        last func to be called, passed to Figure class to generate circos plot (main)
        """

        variants_position = (
            self.options.get("Variants", {}).get("rings", {}).get("position", 0.5)
        )
        variants_ring_height = (
            self.options.get("Variants", {}).get("rings", {}).get("height", 0.04)
        )
        variants_ring_space = (
            self.options.get("Variants", {}).get("rings", {}).get("space", 0.01)
        )

        pass

    def process_vcf(self) -> dict:
        """
        From vcfreader antony explode_category_file_dict_into_dataframe
        """
        # assert isinstance(
        #    self.options["Chromosomes"]["cytoband"], str
        # ) and os.path.exists(
        #        osj(
        #            self.options["Static"],
        #            "Assembly",
        #            self.options["Assembly"],
        #            "cytoband_hg19_chr_infos.txt.gz",
        #        )
        # )
        data = {
            "Chromosomes": [],
            "Genes": [],
            "Exons": [],
            "Variants": [],
            "CopyNumber": [],
        }
        # VCF parsed file from PyVCF3
        for val in self.vcf_reader:
            # data["Chromosomes"].append(val.CHROM)
            # print(val.INFO)

            # print(val.INFO["SV"])
            data["Chromosomes"].append("chr" + val.CHROM)
            data["Genes"].extend(self.get_genes_var(val))
            # TODO exons time consumming
            data["Variants"].append(val.INFO)
            data["CopyNumber"].append(self.get_copynumber_var())
        return data

    def get_copynumber_var(self, record: object) -> int:
        """
        take VCF variant object and return copy number for this variant as an integer
        """
        pass

    def find_record_gene(self, coord: list) -> list:
        """
        Greedy, for now need good info in vcf annotations
        """
        gene_list = []
        # Keep only first transcript per gene should be main
        refgene_genes = pd.read_csv(
            self.refgene_genes, sep="\t", header=0, compression="infer"
        ).drop_duplicates(subset=["gene"], keep="first")
        refgene_genes.loc[refgene_genes["chr_name"] == coord[0]]
        return refgene_genes

    def get_genes_var(self, record: object) -> list:
        gene_name = record.INFO.get("Gene_name")
        if isinstance(gene_name, str):
            gene_name = [gene_name]
        print(record.INFO)
        if gene_name is None:
            assert "SVLEN" in record.INFO
            gene_name = self.find_record_gene(
                [record.CHROM, record.POS, record.POS + record.INFO["SVLEN"]]
            )
        return gene_name

    def formatted_refgene(self, refgene: str, assembly: str) -> str:
        """
        Took refgene raw file from ucsc curated and create proper exon refgene, WITHOUT UTR(default choice)
        """
        df = pd.read_csv(refgene, sep="\t", header=None, compression="infer")
        output_genes = osj(os.path.dirname(refgene), "genes." + assembly + ".txt.gz")
        output_exons = osj(os.path.dirname(refgene), "exons." + assembly + ".txt.gz")
        df.columns = [
            "bin",
            "name",
            "chrom",
            "strand",
            "txStart",
            "txEnd",
            "cdsStart",
            "cdsEnd",
            "exonCount",
            "exonStarts",
            "exonEnds",
            "score",
            "name2",
            "cdsStartStat",
            "cdsEndStat",
            "exonFrames",
        ]
        with gzip.open(output_genes, "wb+") as out_g:
            with gzip.open(output_exons, "wb+") as out_e:
                out_g.write(
                    bytes(
                        "\t".join(
                            [
                                "chr_name",
                                "start",
                                "end",
                                "val",
                                "color",
                                "gene",
                                "transcript",
                            ]
                        )
                        + "\n",
                        "UTF-8",
                    )
                )
                out_e.write(
                    bytes(
                        "\t".join(
                            [
                                "chr_name",
                                "start",
                                "end",
                                "val",
                                "color",
                                "gene",
                                "exons",
                                "transcript",
                            ]
                        )
                        + "\n",
                        "UTF-8",
                    )
                )
                for i, row in tqdm(
                    df.iterrows(),
                    total=len(df.index),
                    desc="Formatting refgene file UCSC",
                    leave=False,
                ):
                    if row["name"].startswith("NM_"):
                        out_g.write(
                            bytes(
                                "\t".join(
                                    [
                                        row["chrom"],
                                        str(row["txStart"]),
                                        str(row["txEnd"]),
                                        "1",
                                        "lightgray",
                                        row["name2"],
                                        row["name"],
                                    ]
                                )
                                + "\n",
                                "UTF-8",
                            )
                        )
                        for i in range(len(row["exonStarts"].split(",")[:-1])):
                            exons_start = row["exonStarts"].split(",")
                            exons_end = row["exonEnds"].split(",")
                            out_e.write(
                                bytes(
                                    "\t".join(
                                        [
                                            row["chrom"],
                                            str(exons_start[i]),
                                            str(exons_end[i]),
                                            "1",
                                            "lightgray",
                                            row["name2"],
                                            "exon" + str(i + 1),
                                            row["name"],
                                        ]
                                    )
                                    + "\n",
                                    "UTF-8",
                                )
                            )
        return df


##############
# Commons func #TODO put in commons.py
def json_to_dict(jsonpath):
    with open(jsonpath) as json_file:
        return json.load(json_file)


def systemcall(command, log=None):
    """
    https://github.com/JbaptisteLam/DPNI/blob/main/src/utils/utils.py
    """
    print("#[SYS] " + command)
    p = subprocess.Popen(
        [command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    out, err = p.communicate()
    if not err:
        return out.decode("utf8").strip().split("\n")
    else:
        issues = err.decode("utf8").strip()
        try:
            re.search(r"(Warning|WARNING)", issues).group()
            print("--WARNING Systemcall--\n", err.decode("utf8").strip())
            return out.decode("utf8").strip().split("\n")
        except AttributeError:
            print("--ERROR Systemcall--\n", err.decode("utf8").strip())
            exit()
