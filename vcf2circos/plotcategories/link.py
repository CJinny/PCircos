from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import (
    generate_hovertext_var,
    delete_multiple_element,
    timeit,
)
import re
from pprint import pprint
from itertools import repeat


class Link(Plotconfig):
    def __init__(self, plotconfig) -> None:
        self.plotconfig = plotconfig
        # self.chr_values_breakend
        # print(self.breakend_record)
        # print(self.breakend_genes)
        # self.merge_options()

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def correct_chosen_var(self, dico):
        values = {"chr1_name": [], "chr2_name": []}
        for key, val in dico.items():
            if key.endswith("_name"):
                for i, chrs in enumerate(val):
                    if chrs not in self.data["Chromosomes"]:
                        values[key].append(i)
        # for v in [*values["chr1_name"], [*values["chr2_name"]]]:
        #    print(v)
        for key, val in dico.items():
            delete_multiple_element(val, [*values["chr1_name"], *values["chr2_name"]])
        return dico

    def adapt_data(self):
        data = []
        # get total number of BND by searching in alt for all var
        for record in self.breakend_record:
            tmp = {}
            tmp[(record.CHROM, record.POS)] = {"record_info": record.INFO, "values": []}
            for alts in record.ALT:
                chrom_pos_tmp = self.extract_chr_pos_hover_from_bnd(
                    alts.__str__(), record
                )
                if not chrom_pos_tmp[0].startswith("chr"):
                    chrom_pos = ("chr" + chrom_pos_tmp[0], chrom_pos_tmp[1])
                else:
                    chrom_pos = (chrom_pos_tmp[0], chrom_pos_tmp[1])
                tmp[(record.CHROM, record.POS)]["values"].append(chrom_pos)
            data.append(tmp)
        return data

    def extract_chr_pos_hover_from_bnd(self, string: str, record: object) -> int:
        try:
            return (
                re.search(r"(?<=\[|\])[^:]+", string).group(),
                int(re.search(r"(?<=:)[0-9]+", string).group()),
                # generate_hovertext_var(record),
            )
        except AttributeError:
            print("ERROR record alt " + string + " malformatted")
            exit()
        # if "[" in rem.group():
        #    pass

        # chr_tmp = string.split(":")[0]
        # pos_tmp = string.split(":")[1]

    @timeit
    def merge_options(self) -> list:
        plot = {}
        data = {
            "chr1_name": [],
            "chr1_start": [],
            "chr1_end": [],
            "chr2_name": [],
            "chr2_start": [],
            "chr2_end": [],
            "color": [],
            "hovertext": [],
            "symbol": [],
            "genes": [],
        }
        tmp = self.adapt_data()
        for link in tmp:
            for key, values in link.items():
                for v in values["values"]:
                    if not key[0].startswith("chr"):
                        chr = "chr" + key[0]
                    else:
                        chr = key[0]
                    data["chr1_name"].append(chr)
                    data["chr1_start"].append(key[1])
                    data["chr1_end"].append(key[1])
                    data["chr2_name"].append(v[0])
                    data["chr2_start"].append(v[1])
                    data["chr2_end"].append(v[1])
                    tmp_gene = []
                    tmp_gene.extend(self.find_record_gene([key[0], key[1], key[1]]))
                    tmp_gene.extend(self.find_record_gene([v[0], v[1], v[1]]))
                    data["genes"].append(",".join(list(set(tmp_gene))))
                    data["hovertext"].extend(
                        list(generate_hovertext_var([values["record_info"]]))
                    )
                    data["color"].append(self.options["Color"]["BND"])
                    data["symbol"].append(0)
        # pprint(data, sort_dicts=False)

        plot["show"] = "True"
        plot["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }
        # R0 stay to zero and R1 follow rings position options
        plot["radius"] = {"R0": 0, "R1": self.options["Variants"]["rings"]["position"]}
        plot["sortbycolor"] = "False"
        plot["colorcolumn"] = 6
        plot["hovertextformat"] = [
            ' "<b>{}:{}-{}<br>{}:{}-{}</b><br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,3], a[i,4], a[i,5], a[i,7])',
            ' "<b>{}:{}-{}<br>{}:{}-{}</b><br><br>{}".format(a[i,3], a[i,4], a[i,5], a[i,0], a[i,1], a[i,2], a[i,7])',
        ]
        plot["trace"] = {
            "uid": "transloc",
            "hoverinfo": "text",
            "marker": {
                "size": 0,
                "symbol": data["symbol"],
                "opacity": 0,
                "color": data["color"],
            },
        }
        plot["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.8,
            "line": {"color": data["color"], "width": 2.5},
        }
        # print(data["chr1_name"])
        # print(data["chr2_name"])
        # exit()
        plot["file"]["dataframe"]["data"] = self.correct_chosen_var(data)
        return plot
