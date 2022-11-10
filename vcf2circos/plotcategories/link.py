from vcf2circos.plotcategories.plotconfig import Plotconfig
import re
from pprint import pprint


class Link(Plotconfig):
    def __init__(self, plotconfig) -> None:
        self.plotconfig = plotconfig
        self.adapt_data()

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def adapt_data(self):
        data = []
        # get total number of BND by searching in alt for all var
        for record in self.breakend_record:
            tmp = {}
            tmp[(record.CHROM, record.POS)] = []
            for alts in record.ALT:
                chrom_pos = self.extract_chr_pos_from_bnd(alts.__str__())
                tmp[(record.CHROM, record.POS)].append(chrom_pos)
            data.append(tmp)

        return

    def extract_chr_pos_from_bnd(self, string: str) -> int:
        try:
            return (
                "chr" + re.search(r"(?<=\[|\])[^:]+", string).group(),
                int(re.search(r"(?<=:)[0-9]+", string).group()),
            )
        except AttributeError:
            print("ERROR record alt " + string + " malformatted")
            exit()
        # if "[" in rem.group():
        #    pass

        # chr_tmp = string.split(":")[0]
        # pos_tmp = string.split(":")[1]

    def merge_options(self) -> list:
        pass
