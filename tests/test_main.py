import sys
from os.path import join as osj
import os


sys.path.append(osj(os.path.dirname(__file__), "..", "vcf2circos", "plotcategories"))
sys.path.append(osj(os.path.dirname(__file__), "..", "vcf2circos"))
print(__file__)
import typing
import pytest
import vcf

import typing
import os
import json

from plotconfig import Plotconfig
from datafactory import Datafactory


input_files = "test.vcf"
with open("test_config.json", "r") as f:
    options = json.load(f)
rangescale = []
val = options["Variants"]["rings"]["position"] + options["Variants"]["rings"]["space"]
rangescale.append(val)
for i in range(options["Variants"]["rings"]["nrings"]):
    val += (
        options["Variants"]["rings"]["height"] + options["Variants"]["rings"]["space"]
    )


def generate_testing_data():
    return Plotconfig(
        input_files,
        options,
        show=True,
        file=None,
        radius=None,
        sortbycolor=None,
        colorcolumn=6,
        hovertextformat=None,
        trace_car=None,
        data=None,
        layout=None,
        rangescale=rangescale,
        config_ring=options["Variants"]["rings"],
    )


def test_get_copynumber_type():
    pc = generate_testing_data()
    for variants in pc.data["Record"]:
        assert (
            pc.get_copynumber_type(variants)[0] in options["Color"].keys()
        ), "unknown SVtype"


# print(test_get_copynumber_type())


def test_get_genes_var():
    pc = generate_testing_data()
    for variants in pc.data["Record"]:
        print(variants)
        print(pc.get_genes_var(variants))


print(test_get_genes_var())
