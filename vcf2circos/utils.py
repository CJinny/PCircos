import json
import subprocess
import re
import pandas as pd
from os.path import join as osj
import os
import gzip

# Globals
variants_color = {
    "INS": "rouge",
    "INV": "purple",
    "DEL": "orange",
    "DUP": "blue",
    "CNV": "brown",
    "OTHER": "gray",
}

# Utils func
def json_to_dict(jsonpath: str) -> dict:
    """
    Load json file and create a dict
    """
    with open(jsonpath) as json_file:
        return json.load(json_file)


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


def formatted_refgene(self, refgene: str, assembly: str) -> str:
    """
    In case of new version of refgene or new assembly\n
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
