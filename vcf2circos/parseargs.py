import argparse
from argparse import RawTextHelpFormatter


class Parseargs:
    def __init__(self):
        pass

    # Define args
    def parseargs(self):
        parser = argparse.ArgumentParser(
            prog="python vcf2circos.py", formatter_class=RawTextHelpFormatter
        )
        parser.add_argument(
            "-i",
            "--input",
            type=str,
            required=True,
            help="Input vcf File\nVCF SHOULD be multiallelic split to avoid trouble in vcf2circos\nexample: bcftools -m -any <vcf>\nFormat will be autodetected from file path.\nSupported format:\n   'vcf.gz', 'vcf'",
        )
        parser.add_argument(
            "-o",
            "--output",
            type=str,
            required=True,
            help="""Output file.\nFormat will be autodetected from file path.\nSupported format:\n   'png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json'""",
        )
        parser.add_argument(
            "-e",
            "--export",
            type=str,
            required=False,
            help="""Export file.\nFormat is 'json'.\nGenerate json file from VCF input file""",
        )
        parser.add_argument(
            "-p",
            "--options",
            type=str,
            required=False,
            help="""Options file or string.\nFormat is 'json', either in a file or as a string.""",
        )
        parser.add_argument(
            "-n",
            "--notebook_mode",
            type=bool,
            required=False,
            help="""Notebook mode.\nDefault False""",
        )
        parser.add_argument(
            "-a",
            "--assembly",
            type=str,
            required=False,
            help="""Genome assembly to use for now values available (hg19, hg38)"""
        )
        # Parse args
        return parser.parse_args()
