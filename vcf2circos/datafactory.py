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
        val = (
            options["Variants"]["rings"]["position"]
            + options["Variants"]["rings"]["space"]
        )
        self.rangescale.append(val)
        for i in range(options["Variants"]["rings"]["nrings"]):
            val += (
                options["Variants"]["rings"]["height"]
                + options["Variants"]["rings"]["space"]
            )
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

        histogram = Histogram_(pc)
        cytoband = Cytoband(pc)
        data_histo = histogram.merge_options(cytoband.data_cytoband())
        # If snv overlapping only
        # pc.data = histogram.data
        # pc.df_data = histogram.df_data

        # pc.data = histogram.data
        ideogram = Ideogram(pc)
        ring = Ring(pc, ["genes"])

        scatter = Scatter_(pc)
        link = Link(pc)
        js = {}
        js["General"] = ideogram.options["General"]

        js["Category"] = {
            "ideogram": ideogram.merge_options(),
            "ring": ring.create_ring(),
            "cytoband": cytoband.merge_options(),
            "histogram": data_histo,
            # "scatter": scatter.merge_options(data_histo),
        }
        with open("without.json", "w+") as o:
            data = json.dumps(js, indent=4)
            o.write(data)
        exit()
        remove_under = []
        for plot_type in js["Category"]:
            if plot_type == "histogram" or plot_type == "scatter":
                for i, val in enumerate(js["Category"][plot_type]):
                    # print(val["file"]["dataframe"]["data"]["chr_name"])
                    if not val["file"]["dataframe"]["data"]["chr_name"]:
                        # print(val["file"]["dataframe"]["data"]["chr_name"])
                        remove_under.append((plot_type, i))
        # Could remove only one ore need to build a copy
        if remove_under:
            # print("DELETE empty")
            del js["Category"][remove_under[0][0]][remove_under[0][1]]
            # print(js["Category"][remove_under[0][0]])

        remove = []
        for plot_type in js["Category"]:
            if plot_type == "histogram" or plot_type == "scatter":
                if not js["Category"][plot_type]:
                    remove.append(plot_type)
        for item in remove:
            del js["Category"][item]
        print(js["Category"].keys())

        # lm = link.merge_options()
        # sm = scatter.merge_options(data_histo)
        # if lm["file"]["dataframe"]["data"]["chr1_name"]:
        #    js["Category"]["link"] = lm
        #
        # histo_list = []
        # scatter_list = []
        # for hist in data_histo:
        #    if hist["file"]["dataframe"]["data"]["chr_name"]:
        #        histo_list.append(hist)
        #    else:
        #        print(hist["file"]["dataframe"]["data"])
        #        print(hist)
        #
        # if histo_list:
        #    js["histogram"] = histo_list
        #
        # for scatt in sm:
        #    if scatt["file"]["dataframe"]["data"]["chr_name"]:
        #        scatter_list.append(scatt)
        # if scatter_list:
        #    js["scatter"] = scatter_list
        ## pprint(js["Category"]["histogram"])
        # pprint(data_histo)

        # for key, val in js["Category"].items():
        #    # print(key)
        #    if isinstance(val, list):
        #        for dico_field in val:

        #            try:
        #                if not dico_field["file"]["dataframe"]["data"]:
        #                    print(key)
        #            except KeyError:
        #                print("errorlist", key)
        #    else:
        #        if key == "link":
        #            print(val["file"]["dataframe"]["data"])
        #        try:
        #            if not val["file"]["dataframe"]["data"]:
        #                print(key)
        #        except KeyError:
        #            print("error", key)
        # exit()
        return js
