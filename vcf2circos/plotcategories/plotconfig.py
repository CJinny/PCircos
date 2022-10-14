#!/usr/local/bin/python3

"""
aim: Epurate module and func, handle parameters of circos plot
"""

import sys
import json
import os
from vcf2circos.vcfreader import VcfReader

print(__file__)
# sys.path.append(os.path.abspath(os.path.join("../", "demo_data")))
# sys.path.append(os.path.join(os.path.dirname(__file__), "demo_data"))#
print(sys.path)

# TODO commons file with utils function maybe one more for globals


class Plotconfig(VcfReader):
    """
    Options regroup options passed in args in json file otherwise
    it will be a empty dict,
    All func based on vcf input
    """

    def __init__(self, filename, options: dict = ...):
        if os.path.exists(
            "/maison/lamouche/dev_vcf2circos/demo_data/options.general.json"
        ):
            print(True)
        else:
            print(False)
        super().__init__(filename, options)
        self.default_options = json.load(
            open("/maison/lamouche/dev_vcf2circos/demo_data/options.general.json", "r")
        )
        if not self.options.get("General", {}).get("title", None):
            self.options["General"]["title"] = os.path.basename(
                self.get_metadatas().get("filename", "myCircos")
            )

    def get_snvindels_overlapping_sv(self):
        pass

    def vcf_options_default(self):
        pass

    def get_json(self):
        """
        last func to be called, passed to Figure class to generate circos plot (main)
        """
        pass
