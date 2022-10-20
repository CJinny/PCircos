#!/usr/local/bin/python3

"""
aim: Epurate module and func, handle parameters of circos plot
"""

import sys
import re
import subprocess
import json
import os

from vcf2circos.vcfreader import VcfReader
from os.path import join as osj


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
            open(
                "../demo_data/options.general.json",
                "r",
            )
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
        self.data = data
        self.layout = layout

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

        # VCF parsed file from PyVCF3
        chroms_list = []
        for val in self.vcf_reader:
            chroms_list.append(val.CHROM)
        return list(set(chroms_list))


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
