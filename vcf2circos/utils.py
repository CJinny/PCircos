import json
import subprocess
import re
import pandas as pd
from os.path import join as osj
import os
import gzip
from tqdm import tqdm
import vcf2circos

# from natsort import natsort_keygen
from typing import Generator
from functools import wraps
import time
import numpy as np
import pyfiglet

# Globals
variants_color = {
    "INS": "red",
    "INV": "purple",
    "DEL": "orange",
    "DUP": "blue",
    "CNV": "brown",
    "BND": "blue",
    "SNV": "gray",
    "INDEL": "gray",
    "OTHER": "gray",
}


def launch():
    """
    https://ascii.co.uk/art/dna
    """
    print(
        """
-._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _   
    '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '.` : '.   
  '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.:   '.  '.
  : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.  : '.  `.
  '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.'   `.
         `-..,..-'       `-..,..-'       `-..,..-'       `         `
    """
    )
    print(pyfiglet.figlet_format("         vcf2circos", font="big"))
    print("Author: " + vcf2circos.__author__)
    print("Version: " + vcf2circos.__version__)
    print("Last update: " + vcf2circos.__date__)
    print("\n")


class Colorpal:
    def __init__(self, n=None):
        self.color = self.random_rgb()
        if n is not None:
            self.colorpal = [self.random_rgb() for c in range(0, n)]

    # def __call__(self):
    #    if len(self.colorpal) == 0:
    #        print("Specify length of color pal EXIT")
    #        exit()
    #    else:
    #        return self.colorpal

    def __iter__(self):
        try:
            for v in self.colorpal:
                yield v
        except AttributeError:
            print("Specify length of color pal EXIT")
            exit()

    def random_rgb(self):
        return tuple(map(str, list(np.random.choice(range(255), size=3))))


# Utils func
def json_to_dict(jsonpath: str) -> dict:
    """
    Load json file and create a dict
    """
    with open(jsonpath) as json_file:
        return json.load(json_file)


def check_data_plot(dico, list_keys=None):
    try:
        var_numb = len(dico["chr_name"])
    except KeyError:
        print(dico.keys())
        exit()
    for d_fields in dico:
        assert len(dico[d_fields]) == var_numb, (
            "Missing values in "
            + d_fields
            + " fields\n\t...chrom: "
            + str(var_numb)
            + " "
            + d_fields
            + ": "
            + str(len(dico[d_fields]))
        )
    if "end" in dico:
        pass
    # else:
    #    if list_keys is None:
    #        list_keys = ["chr_name", "start", "val", "ref", "alt", "type", "color"]
    #    assert list_keys == list(dico.keys())[:7], (
    #        "ISSUES wrong list data order, leads to crash in mathematical operations \n\t..."
    #        + ", ".join(list(dico.keys())[:7])
    #        + " ,..."
    #    )


def delete_multiple_element(list_object, indices):
    """
    from https://thispointer.com/python-remove-elements-from-list-by-index/
    """
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)


def map_annotations(field_annot):
    if field_annot is not None:
        uniq = list(set(str(field_annot).split("|")))
        if len(uniq) == 1:
            return uniq[0]
        else:
            return field_annot
    else:
        return "."


def generate_hovertext_var(
    variants_list, full_annot=None, true_annot=None
) -> Generator:
    # print(self.data["Variants"])
    # print(len(self.data["Variants"]))
    # print(len(self.data["Chromosomes"]))
    # exit()
    # dict containing INFO field for each var
    # print(variants_list)
    # for var in variants_list:
    #    yield "<br>".join(
    #        [
    #            ": ".join(
    #                [
    #                    str(value) if not isinstance(value, list) else str(value[0])
    #                    for value in pairs
    #                ]
    #            )
    #            for pairs in list(zip(var.keys(), var.values()))
    #        ]
    #    )
    # 30 longueur char
    # 15 hauteur annot
    for var in variants_list:
        tmp = []
        for i, pairs in enumerate(list(zip(var.keys(), var.values()))):
            # if pairs[0] == "OMIM_phenotype":
            # print(pairs[1])
            if true_annot:
                # If user want this annotations
                if pairs[0] not in true_annot:
                    # print(pairs[0] + " not in")
                    continue

            if full_annot is not None:
                if i == full_annot:
                    break
            else:
                if i == 15:
                    break
            if not isinstance(pairs[1], list):
                tmp.append(
                    ":".join([pairs[0], list(map(map_annotations, [pairs[1]]))[0]])
                )
            else:
                tmp.append(
                    ": ".join(
                        [pairs[0], ",".join(list(map(map_annotations, pairs[1])))]
                    )
                )
        to_add = []
        for items in tmp:
            if len(items) > 40:
                items = "".join(items[:40]) + "..."
                to_add.append(items)
            else:
                to_add.append(items)
        yield "<br>".join(to_add)
        # exit()
        # "SV_chrom",


# "SV_start",
# "SV_end",
# "FORMAT",
# "SpliceAI",
# "SPiP",
# "ACMG",
# "varankVarScore",
# "gene",
# "zygisity",
# "rsClinicalSignificance",
# "OMIM_ID",
# "OMIM_inheritance",
# "OMIM_phenotype"


# def generate_hovertext_var(variants_list) -> Generator:
#    # print(self.data["Variants"])
#    # print(len(self.data["Variants"]))
#    # print(len(self.data["Chromosomes"]))
#    # exit()
#    # dict containing INFO field for each var
#    for var in variants_list:
#        yield "<br>".join(
#            [
#                ": ".join(
#                    [
#                        str(value) if not isinstance(value, list) else str(value[0])
#                        for value in pairs
#                    ]
#                )
#                for pairs in list(zip(var.keys(), var.values()))
#            ]
#        )


def cast_svtype(svtype):
    """
    In case of pip in svtype
    """
    return svtype.split("|")[0]


def systemcall(command: str) -> list:
    """
    Call bash command
    https://github.com/JbaptisteLam/DPNI/blob/main/src/utils/utils.py

    In case of crash exit code 1 stop script
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


def formatted_refgene(refgene: str, assembly: str, ts=None) -> str:
    """
    In case of new version of refgene or new assembly\n
    transcripts: list in string format either 'NM_' or 'NR_ 'or both 'NM_,NR_'
    Took refgene raw file from ucsc curated and create proper exon refgene, WITHOUT UTR(default choice)
    https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=refSeqComposite&db=hg38
    ex: python -c 'from foo import hello; print hello()'
     zcat transcripts.hg38.txt.gz.tmp | grep chr_name > transcripts.hg38.sorted.txt && zcat transcripts.hg38.txt.gz.tmp | grep -v "chr_name" | sort -k1,1V -k2,2n >> transcripts.hg38.sorted.txt
     same for both 3 files except exons sort after chr and pos by exons
     zcat exons.hg38.txt.gz.tmp | grep chr_name > exons.hg38.sorted.txt && zcat exons.hg38.txt.gz.tmp | grep -v "chr_name" | sort -k1,1V -k2,2n >> exons.hg38.sorted.txt
    """
    if ts is not None:
        transcripts = ts.split(",")
    else:
        transcripts = "NM_"
    df = pd.read_csv(refgene, sep="\t", header=None, compression="infer")
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
    output_genes = os.path.dirname(refgene) + "/genes." + assembly + ".txt.gz.tmp"
    output_transcripts = (
        os.path.dirname(refgene) + "/transcripts." + assembly + ".txt.gz.tmp"
    )
    output_exons = os.path.dirname(refgene) + "/exons." + assembly + ".txt.gz.tmp"

    og = os.path.dirname(refgene) + "/genes." + assembly + ".txt.gz"
    ot = os.path.dirname(refgene) + "/transcripts." + assembly + ".txt.gz"
    oe = os.path.dirname(refgene) + "/exons." + assembly + ".txt.gz"

    with gzip.open(output_genes, "wb+") as out_g:
        with gzip.open(output_exons, "wb+") as out_e:
            with gzip.open(output_transcripts, "wb+") as out_t:
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
                            ]
                        )
                        + "\n",
                        "UTF-8",
                    )
                )
                out_t.write(
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
                    desc="Formatting " + os.path.basename(refgene) + " file UCSC",
                    leave=False,
                ):
                    out_g.write(
                        bytes(
                            "\t".join(
                                [
                                    row["chrom"],
                                    str(row["cdsStart"]),
                                    str(row["cdsEnd"]),
                                    "1",
                                    str(row["name2"]),
                                    str(row["name2"]),
                                ]
                            )
                            + "\n",
                            "UTF-8",
                        )
                    )
                    if [row["name"].startswith(ts) for ts in transcripts]:
                        out_t.write(
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
    # systemcall(
    #    "zcat "
    #    + output_genes
    #    + " | grep 'chr_name' > "
    #    + og
    #    + " && zcat "
    #    + output_genes
    #    + " | grep -v 'chr_name' "
    #    + output_genes
    #    + " | sort -k1,1V -k2,2n > "
    #    + og
    # )
    # systemcall(
    #    "grep chr_name "
    #    + output_genes
    #    + " > "
    #    + og
    #    + " && grep -v chr_name "
    #    + output_genes
    #    + " | sort -k1,1V -k2,2n > "
    #    + og
    # )
    #
    # systemcall(
    #    "grep chr_name "
    #    + output_exons
    #    + " > "
    #    + oe
    #    + " && grep -v chr_name "
    #    + output_exons
    #    + " | sort -k1,1V -k2,2n > "
    #    + oe
    # )
    return df


def timeit(func):
    """
    https://dev.to/kcdchennai/python-decorator-to-measure-execution-time-54hk
    """

    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f"Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds")
        return result

    return timeit_wrapper
