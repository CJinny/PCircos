import json
import subprocess
import re
import pandas as pd
from os.path import join as osj
import os
import gzip
from tqdm import tqdm
from natsort import natsort_keygen
from functools import wraps
import time

# Globals
variants_color = {
    "INS": "red",
    "INV": "purple",
    "DEL": "orange",
    "DUP": "blue",
    "CNV": "brown",
    "SNV": "gray",
    "INDEL": "gray",
    "OTHER": "gray",
}

# Utils func
def json_to_dict(jsonpath: str) -> dict:
    """
    Load json file and create a dict
    """
    with open(jsonpath) as json_file:
        return json.load(json_file)


def check_data_plot(dico):
    var_numb = len(dico["chr_name"])
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
    else:
        assert ["chr_name", "start", "val", "ref", "alt", "type", "color",] == list(
            dico.keys()
        )[:7], (
            "ISSUES wrong list data order, leads to crash in mathematical operations \n\t..."
            + ", ".join(list(dico.keys())[:7])
            + " ,..."
        )


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
                                    omim_morbid(row["name2"]),
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
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f"Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds")
        return result

    return timeit_wrapper
