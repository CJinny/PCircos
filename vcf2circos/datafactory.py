from vcf2circos.plotcategories.histogram import Histogram_
from vcf2circos.plotcategories.ideogram import Ideogram
from vcf2circos.plotcategories.scatter import Scatter_
from vcf2circos.plotcategories.ring import Ring
from vcf2circos.plotcategories.cytoband import Cytoband
from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.plotcategories.link import Link
from pprint import pprint
import json


class Datafactory:
    def __init__(self, input_file, options) -> None:
        self.input_file = input_file
        self.options = options
        self.rangescale = []
        val = options["Variants"]["rings"]["position"] + options["Variants"]["rings"]["space"]
        self.rangescale.append(val)
        for i in range(options["Variants"]["rings"]["nrings"]):
            val += options["Variants"]["rings"]["height"] + options["Variants"]["rings"]["space"]
            self.rangescale.append(val)

    # Read vcf and process raw data to feed child class
    def plot_dict(self):
        pc = Plotconfig(
            filename=self.input_file,
            options=self.options.copy(),
            show=True,
            file=None,
            radius=None,
            sortbycolor=None,
            colorcolumn=6,
            hovertextformat=None,
            trace_car=None,
            data=None,
            layout=None,
            rangescale=self.rangescale,
            config_ring=self.options["Variants"]["rings"],
        )

        # Ugly as hell, if we wanna take only snv indel overlapping SV
        # Create plot object
        link = Link(pc)
        histogram = Histogram_(pc)
        cytoband = Cytoband(pc)
        link_me, chr_values = link.merge_options()
        data_histo = histogram.merge_options(cytoband.data_cytoband(chr_values))
        ideogram = Ideogram(pc)
        ring = Ring(pc, ["genes"])
        scatter = Scatter_(pc, data_histo.copy())

        # Final dict containing all plots informations
        js = {}
        js["General"] = ideogram.options["General"]
        # print("HISTO\n")
        # for ite in data_histo:
        #    print(ite["file"]["dataframe"]["data"].keys())

        js["Category"] = {
            "ideogram": ideogram.merge_options(chr_values),
            "ring": ring.create_ring(),
            "cytoband": cytoband.merge_options(chr_values),
            "histogram": data_histo,
            "scatter": scatter.merge_options(),
            "link": link_me,
        }

        # DEBUg
        # with open("without.json", "w+") as o:
        #    data = json.dumps(js, indent=4)
        #    o.write(data)
        # exit()
        # Adjustement in case of no data for example when use overlapping snv only
        remove_under = []
        plot_check = ["histogram", "scatter"]
        for plot_type in js["Category"]:
            if plot_type in plot_check:
                # Only for list
                if isinstance(js["Category"][plot_type], list):
                    for i, val in enumerate(js["Category"][plot_type]):
                        # print(val["file"]["dataframe"]["data"]["chr_name"])
                        if val["trace"]["uid"].startswith("cnv_") or val["trace"]["uid"].startswith(
                            "genes"
                        ):  # genes in cases of only translocation
                            if not val["file"]["dataframe"]["data"]["chr_name"]:
                                # print(val["file"]["dataframe"]["data"]["chr_name"])
                                remove_under.append((plot_type, i))
        ## Could remove only one ore need to build a copy
        if remove_under:
            for fields in remove_under:
                print(
                    "#[INFO] index of category to remove: "
                    + ", ".join([items[0] for items in remove_under])
                )
                # print("DELETE empty")
                del js["Category"][fields[0]][fields[1]]
                # print(js["Category"][remove_under[0][0]])
        remove = []
        for plot_type in js["Category"]:
            if plot_type in plot_check:
                if not js["Category"][plot_type]:
                    remove.append(plot_type)
            elif plot_type == "link":
                # exit()
                if not js["Category"][plot_type]["file"]["dataframe"]["data"]["chr1_name"]:
                    remove.append(plot_type)
        if remove:
            print("#[INFO] Whole category to remove: " + ", ".join(remove))
            for item in remove:
                del js["Category"][item]
                # mean only transloc or no variants
            if "scatter" in remove and not "link" in remove:
                del js["Category"]["ring"]
                js["Category"]["link"]["radius"] = {"R0": 0, "R1": 0.98}

        # DEBUG dico values
        # test_ = []
        # print("\n")
        # for item in js["Category"]["histogram"]:
        #    if "dataframe" in item["file"] and item["trace"]["uid"].endswith("level_5"):
        #        print(item["file"]["dataframe"]["data"])
        #        print(item["trace"]["uid"])
        #        print(item["hovertextformat"])
        #        print("\n")
        #        test_.append(item)
        #    elif "dataframe" in item["file"] and item["trace"]["uid"].startswith("cnv"):
        #        continue
        #    else:
        #        # print(item)
        #        test_.append(item)
        return js
