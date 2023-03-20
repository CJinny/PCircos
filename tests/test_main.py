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

#Testfiles
input_files_copy = "test_copy_number.vcf"
input_files_genes = "test_gene_matcher.vcf"
input_files_transloc = "test_transloc.vcf"
input_files_empty = "test_empty.vcf"

with open("test_config.json", "r") as f:
    options = json.load(f)
rangescale = []
val = options["Variants"]["rings"]["position"] + options["Variants"]["rings"]["space"]
rangescale.append(val)
for i in range(options["Variants"]["rings"]["nrings"]):
    val += (
        options["Variants"]["rings"]["height"] + options["Variants"]["rings"]["space"]
    )

genes_check = ['ABHD12B', 'ACTR10', 'AP5M1', 'ARF6', 'ARID4A', 'ARMH4', 'ATG14', 'ATL1', 'BMP4', 'CCDC175', 'CCDC198', 'CDKL1', 'CDKN3', 'CGRRF1', 'CNIH1', 'DAAM1', 'DACT1', 'DDHD1', 'DLGAP5', 'DMAC2L', 'DNAAF2', 'ERO1A', 'EXOC5', 'FBXO34', 'FERMT2', 'FRMD6', 'FRMD6-AS1', 'FRMD6-AS2', 'GCH1', 'GMFB', 'GNG2', 'GNPNAT1', 'GPR135', 'GPR137C', 'JKAMP', 'KIAA0586', 'KLHDC1', 'KLHDC2', 'KTN1', 'KTN1-AS1', 'L2HGDH', 'L3HYPDH', 'LGALS3', 'LINC00216', 'LINC00519', 'LINC00520', 'LINC00640', 'LINC01500', 'LINC01588', 'LINC01599', 'LINC02310', 'LOC101927620', 'LOC101927690', 'LOC102723604', 'LOC105370489', 'LRR1', 'MAP4K5', 'MAPK1IP1L', 'MGAT2', 'MIR4308', 'MIR4504', 'MIR5580', 'MIR6076', 'NAA30', 'NEMF', 'NID2', 'NIN', 'OTX2', 'OTX2-AS1', 'PELI2', 'POLE2', 'PSMA3', 'PSMA3-AS1', 'PSMC6', 'PTGDR', 'PTGER2', 'PYGL', 'RN7SL1', 'RN7SL2', 'RN7SL3', 'RPL13AP3', 'RPL36AL', 'RPS29', 'RTRAF', 'SAMD4A', 'SAV1', 'SLC35F4', 'SOCS4', 'SOS2', 'STYX', 'TBPL2', 'TIMM9', 'TMEM260', 'TMX1', 'TOMM20L', 'TRIM9', 'TXNDC16', 'VCPKMT', 'DLGAP5', 'KTN1', 'LOC101927690', 'ARF6', 'ATG14', 'L2HGDH', 'SLC35F4', 'PSMA3-AS1', 'LINC01588', 'DDHD1', 'MIR4308', 'DAAM1', 'PELI2', 'PTGDR', 'LOC105370489', 'TXNDC16', 'TBPL2', 'ARID4A', 'TOMM20L', 'DNAAF2', 'ABHD12B', 'FERMT2', 'LINC02310', 'DACT1', 'ACTR10', 'KLHDC2', 'FRMD6-AS1', 'ERO1A', 'TMEM260', 'GPR137C', 'RPL13AP3', 'RPL36AL', 'LRR1', 'RTRAF', 'PTGER2', 'LGALS3', 'EXOC5', 'LINC01500', 'NID2', 'PSMC6', 'TIMM9', 'RN7SL3', 'AP5M1', 'GPR135', 'FRMD6-AS2', 'RN7SL2', 'NEMF', 'KTN1-AS1', 'FRMD6', 'ARMH4', 'CDKL1', 'KIAA0586', 'POLE2', 'NAA30', 'WDHD1', 'SOS2', 'BMP4', 'CDKN3', 'MIR4504', 'KLHDC1', 'MAPK1IP1L', 'LINC01599', 'OTX2-AS1', 'ATL1', 'FBXO34', 'CGRRF1', 'CNIH1', 'LINC00519', 'MIR5580', 'L3HYPDH', 'GNPNAT1', 'RN7SL1', 'LINC00216', 'LOC101927620', 'TRIM9', 'GMFB', 'CCDC198', 'MGAT2', 'RPS29', 'JKAMP', 'PYGL', 'SOCS4', 'STYX', 'LINC00520', 'DMAC2L', 'GCH1', 'TMX1', 'NIN', 'CCDC175', 'VCPKMT', 'PSMA3', 'LINC00640', 'GNG2', 'SAMD4A', 'MIR6076', 'OTX2', 'LOC102723604', 'MAP4K5', 'SAV1']

def generate_testing_data(input: str) -> object:
    return Plotconfig(
        input,
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
    pc = generate_testing_data(input_files_copy)
    for variants in pc.data["Record"]:
        assert (
            pc.get_copynumber_type(variants)[0] in options["Color"].keys()
        ), "unknown SVtype"


def test_get_genes_var():
    pc = generate_testing_data(input_files_genes)
    for variants in pc.data["Record"]:
        if variants.CHROM == "chr14" and variants.POS == 50000000:
            genes_list = pc.get_genes_var(variants).split(",")
            assert all(genes for genes in genes_list if genes in genes_check) and len(genes_list) == len(genes_check), "ERROR Genes finding have a problem \nGenes needed: "+str(len(genes_check))+" genes \nGenes: "+' ,'.join(genes_check)+"\nGenes provided: "+str(len(genes_list))+" genes \n Genes: "+' ,'.join(genes_list)

def test_datafactory_type():
    #pc = generate_testing_data(input_files_transloc)
    plots = []
    plots_expected = ["ideogram", "cytoband", "histogram", "link"]
    for k in Datafactory(input_files_transloc, options).plot_dict()["Category"]:
        plots.append(k)
    assert plots == plots_expected, "ERROR wrong category chart in list \nNeeded: "+' ,'.join(plots_expected)+"\nProvided: "+' ,'.join(plots)

def test_empty_vcf():
    pc = generate_testing_data(input_files_empty)
    assert len(pc.data["Variants"]) == 0, "No variant should be present, variants at least one "+str(pc.data["Variants"][0])
