# Standard imports
import copy
from genericpath import isfile
import re
from sre_compile import isstring
import textwrap
import vcf
import json
import os
from colour import Color
from os.path import join as osj
import sys


import random
import pandas as pd

TEXTWRAP_WIDTH = 60
TEXTWRAP_HOLDER = "..."

# Fixing PyVCF bug
# https://github.com/jamescasbon/PyVCF/pull/320
def _map(self, func, iterable, bad=[".", "", "NA", "-"]):
    """``map``, but make bad values None."""

    def _convert(x):
        if x in bad:
            return None
        try:
            return func(x)
        except Exception as e:
            # LOGGER.exception(e)
            return None

    return [_convert(x) for x in iterable]


vcf.Reader._map = _map

# End fixing


# Read file to dict
def file_to_dict(file: str = None, path: str = ""):
    # Result is dict in input
    result = file
    # Result is dict from the file
    if file and isstring(file):
        result = {}
        file_path = None
        if os.path.isfile(file):
            file_path = file
        elif os.path.isfile(os.path.dirname(path) + "/" + file):
            file_path = os.path.dirname(path) + "/" + file
        if file_path:
            with open(file_path, "r") as f:
                result = json.loads(f.read())
            f.close()
        else:
            print(f"[WARN] file {file} does NOT exist")
    return result


# explode file into dataframe
def explode_category_file_dict_into_dataframe(file_dict: dict = {}, path: str = ""):

    result = file_dict

    # data in tsv file
    category_path_file = None
    if os.path.isfile(result.get("path", None)):
        category_path_file = result.get("path")
    elif os.path.isfile(os.path.dirname(path) + "/" + result.get("path", None)):
        category_path_file = os.path.dirname(path) + "/" + result.get("path", "")

    # data in json file
    category_dataframe_file = None
    if os.path.isfile(str(result.get("dataframe", ""))):
        category_dataframe_file = result.get("dataframe", "")
    elif os.path.isfile(os.path.dirname(path) + "/" + str(result.get("dataframe", ""))):
        category_dataframe_file = (
            os.path.dirname(path) + "/" + result.get("dataframe", "")
        )

    # explode into dataframe
    if category_path_file:
        chr_info_pd = pd.read_csv(
            category_path_file,
            sep=result.get("sep", "\t"),
            header=result.get("head", "infer"),
        )
        category_path_file_dataframe = {}
        for i in chr_info_pd.keys():
            category_path_file_dataframe[i] = list(chr_info_pd[i])
        result["dataframe"] = {
            "orient": "columns",
            "data": category_path_file_dataframe,
        }
        result["dataframe"]["data"] = category_path_file_dataframe
        result["path"] = ""
    elif category_dataframe_file:
        with open(category_dataframe_file, "r") as f:
            category_dataframe_dict = json.loads(f.read())
        f.close()
        result["dataframe"] = category_dataframe_dict

    return result


def find_genes_for_record(record_data: dict = {}, genes_data: dict = {}):

    record_positions = []

    if "chr1_name" in record_data:
        record_positions.append(
            {
                "chr": record_data.get("chr1_name", "ERROR"),
                "start": record_data.get("chr1_start", 0),
                "end": record_data.get("chr1_end", record_data.get("chr1_start", 0)),
            }
        )
        record_positions.append(
            {
                "chr": record_data.get("chr_name", "ERROR"),
                "start": record_data.get("chr2_start", 0),
                "end": record_data.get("chr2_end", record_data.get("chr2_start", 0)),
            }
        )
    else:
        record_positions.append(
            {
                "chr": record_data.get("chr_name", "ERROR"),
                "start": record_data.get("start", 0),
                "end": record_data.get("end", record_data.get("start", 0)),
            }
        )

    record_genes = {}

    if len(record_data) and len(genes_data):
        for i in range(len(genes_data["chr_name"])):
            gene_chr = genes_data["chr_name"][i]
            gene_start = genes_data["start"][i]
            gene_end = genes_data["end"][i]
            gene_gene = genes_data["gene"][i]
            if "exon" in genes_data:
                gene_gene += "-" + genes_data["exon"][i]
            for record_position in record_positions:
                if (
                    record_position.get("chr") == gene_chr
                    and record_position.get("start") <= gene_end
                    and record_position.get("end") >= gene_start
                ):
                    record_genes[gene_gene] = {
                        "chr": gene_chr,
                        "start": gene_start,
                        "end": gene_end,
                    }

    return record_genes


VCF_TYPE_MAPPING = {
    "Float": "float",
    "Integer": "int",
    "Flag": "bool",
    "String": "str",
    "Character": "str",
}


class VcfReader:
    """VCF parser to extract data from vcf file

    Attributes:

    """

    def __init__(self, filename, options: dict = {}):
        """Construct a VCF Reader

        .. note::

        -filename: File filename handler returned by open.
        -options (dict): list of options
            This argument list options mainly to generate json for plot.
        """

        # params

        self.filename = filename
        self.options = options

        self.vcf_reader = vcf.VCFReader(
            filename=filename, strict_whitespace=True, encoding="utf-8"
        )
        # Fields descriptions
        self.fields = None
        self.get_fields()
        self.samples = None
        self.get_samples()
        self.metadatas = None
        self.get_metadatas()
        self.infos = None
        self.get_infos()
        self.formats = None
        self.get_formats()
        self.contigs = None
        self.get_contigs()

        ### options

        # General

        # general_default = {"width": 1500, "height": 1500, "title": None}

        # General title
        # add filename as title if not in options

        # Annotations

        variants_default = {
            "annotations": {"fields": ["chr", "pos", "ref", "alt"], "show_none": False},
            "rings": {"position": 0.50, "height": 0.04, "space": 0.01},
        }

        # force Variants section
        if not self.options.get("Variants", {}):
            self.options["Variants"] = variants_default
        else:
            # check Annotations section
            if not self.options.get("Variants", {}).get("annotations", {}):
                self.options["Variants"]["annotations"] = variants_default.get(
                    "annotations", {}
                )
            else:
                # annotations fields
                if (
                    not self.options.get("Variants", {})
                    .get("annotations", {})
                    .get("fields", None)
                ):
                    self.options["Variants"]["annotations"][
                        "fields"
                    ] = variants_default.get("annotations", {}).get("fields", [])
                # annotations none value
                if (
                    self.options.get("Variants", {})
                    .get("annotations", {})
                    .get("show_none", "None")
                    == "None"
                ):
                    self.options["Variants"]["annotations"][
                        "show_none"
                    ] = variants_default.get("annotations", {}).get("show_none", False)
            # check rings positions
            if not self.options.get("Variants", {}).get("rings", {}):
                self.options["Variants"]["rings"] = variants_default.get("rings", {})
            else:
                # rings first position
                if (
                    not self.options.get("Variants", {})
                    .get("rings", {})
                    .get("position", None)
                ):
                    self.options["Variants"]["position"] = variants_default.get(
                        "rings", {}
                    ).get("position", 0.5)
                # rings height
                if (
                    not self.options.get("Variants", {})
                    .get("rings", {})
                    .get("height", None)
                ):
                    self.options["Variants"]["height"] = variants_default.get(
                        "rings", {}
                    ).get("height", 0.04)
                # ring space
                if (
                    not self.options.get("Variants", {})
                    .get("rings", {})
                    .get("space", None)
                ):
                    self.options["Variants"]["space"] = variants_default.get(
                        "rings", {}
                    ).get("space", 0.01)

        self.progress_every = 100
        self.total_bytes = self.vcf_reader.total_bytes()
        self.read_bytes = 0

    def get_fields(self):
        """Get full fields descriptions

        This function is called a first time before variants insertion.

        Annotations fields are added here if they exist in the file.

        .. seealso:: :meth:`parse_fields` for basic default fields.

        :return: Tuple of fields.
            Each field is a dict with the following keys:
            `name, category, description, type`.
            Some fields have an additional constraint key when they are destined
            to be a primary key in the database.
            Annotations fields are added here if they exist in the file.
        :rtype: <tuple <dict>>
        """
        if self.fields is None:
            self.fields = tuple(self.parse_fields())

        return self.fields

    def get_variants(self):
        """Get variants as an iterable of dictionnaries

        "annotations" key is added here with the list of annotations if
        they exist in the file.

        .. seealso:: parse_variants()

        :return: Generator of full variants with "annotations" key.
        :rtype: <generator <dict>>
        """
        if self.fields is None:
            # This is a bad caching code ....
            self.get_fields()

        return self.parse_variants()

    def parse_variants(self, add_sample: bool = False):
        """Read file and parse variants

        1 variant is created for each alternative allele detected for each record.
        For each variant we add the corresponding genotype of each sample under
        the key "samples".

        Examples:
            For a record: `"REF": "A", "ALT": ["T", "C"]` with 2 samples with the
            following genotypes: 0/0 and 0/1 (ref/ref and ref/alt).
            We create 2 variants (because of the 2 alternative alleles),
            each with 2 samples with the following genotypes: 0 and 1.

            Where homozygous_ref = 0, heterozygous = 1, homozygous_alt = 2.

            We don't track which alternative allele is in the genotype of the
            sample.

        See Also:
            https://pyvcf.readthedocs.io/en/v0.4.6/INTRO.html
            https://pyvcf.readthedocs.io/en/latest/API.html#vcf.model._Call.gt_type

        See Also:
            :meth:`cutevariant.core.reader.abstractreader.AbstractReader.get_extra_variants`

        :return: Generator of variants.
        :rtype: <generator <dict>>
        """
        # loop over record
        vcf_reader = self.vcf_reader
        # Genotype format fields
        format_fields = set(map(str.lower, vcf_reader.formats))
        # Remove gt field (added manually later)
        format_fields.discard("gt")

        for i, record in enumerate(vcf_reader):

            self.read_bytes = vcf_reader.read_bytes()

            # split row with multiple alt
            for index, alt in enumerate(record.ALT):
                # Remap some columns
                variant = {
                    "chr": record.CHROM,
                    "pos": record.POS,
                    "ref": record.REF,
                    "alt": str(alt),
                    "rsid": record.ID,  # Avoid id column duplication in DB
                    "qual": record.QUAL,
                    "filter": "" if record.FILTER is None else ",".join(record.FILTER),
                }

                forbidden_field = ("chr", "pos", "ref", "alt", "rsid", "qual", "filter")

                # Parse info
                for name in record.INFO:
                    if name.lower() not in forbidden_field:
                        if isinstance(record.INFO[name], list):
                            variant[name.lower()] = ",".join(
                                [str(i) for i in record.INFO[name]]
                            )
                        else:
                            variant[name.lower()] = record.INFO[name]

                # Parse sample(s)
                if record.samples and add_sample:
                    variant["samples"] = []
                    for sample in record.samples:
                        # New sample data
                        sample_data = {
                            "name": sample.sample,
                            "gt": -1 if sample.gt_type is None else sample.gt_type,
                        }

                        # Load sample fields
                        # 1 genotype field per format
                        # In theory: All same fields for each sample
                        # print("FORMAT FIELD",format_fields, sample["GQ"])
                        for gt_field in format_fields:
                            try:
                                value = sample[gt_field.upper()]
                                if isinstance(value, list):
                                    value = ",".join(str(i) for i in value)
                                sample_data[gt_field] = value
                            except AttributeError:
                                # Some fields defined in VCF header by FORMAT data
                                # are not in genotype fields of records...
                                # LOGGER.debug(
                                #     "VCFReader::parse: alt index %s; %s not defined in genotype ", index, gt_field
                                # )
                                pass
                        variant["samples"].append(sample_data)
                yield variant

        self.read_bytes = self.total_bytes

    def parse_fields(self):
        """Extract fields informations from VCF fields

        .. note:: Fields used in PRIMARY KEYS have the constraint NOT NULL.
            By default, all other fields can have NULL values.

        :return: Generator of fields.
            Each field is a dict with the following keys:
            `name, category, description, type`.
            Some fields have an additional constraint key when they are destined
            to be a primary key in the database.
        :rtype: <generator <dict>>
        """

        yield {
            "name": "chr",
            "category": "variants",
            "description": "Chromosome",
            "type": "str",
            "constraint": "NOT NULL",
        }
        yield {
            "name": "pos",
            "category": "variants",
            "description": "Reference position, with the 1st base having position 1",
            "type": "int",
            "constraint": "NOT NULL",
        }
        yield {
            "name": "ref",
            "category": "variants",
            "description": "Reference base",
            "type": "str",
            "constraint": "NOT NULL",
        }
        yield {
            "name": "alt",
            "category": "variants",
            "description": "Alternative base",
            "type": "str",
            "constraint": "NOT NULL",
        }
        yield {
            "name": "rsid",
            "category": "variants",
            "description": "rsid unique identifier (see dbSNP)",
            "type": "str",
        }
        yield {
            "name": "qual",
            "category": "variants",
            "description": "Phred-scaled quality score for the assertion made in ALT: -10log10 prob(call in ALT is wrong).",
            "type": "int",
        }
        yield {
            "name": "filter",
            "category": "variants",
            "description": "Filter status: PASS if this position has passed all filters.",
            "type": "str",
        }

        # Read VCF
        vcf_reader = vcf.VCFReader(
            filename=self.filename, strict_whitespace=True, encoding="utf-8"
        )

        # Read VCF INFO fields
        for field_name, info in vcf_reader.infos.items():
            yield {
                "name": field_name.lower(),
                "category": "variants",
                "description": info.desc,
                "type": VCF_TYPE_MAPPING[info.type],
            }

        # Read VCF FORMAT fields
        for field_name, info in vcf_reader.formats.items():
            description = info.desc
            field_type = VCF_TYPE_MAPPING[info.type]

            if field_name == "GT":
                # Edit description of Genotype field
                description += (
                    " (0: homozygous_ref, 1: heterozygous, 2: homozygous_alt)"
                )
                field_type = VCF_TYPE_MAPPING["Integer"]

            yield {
                "name": field_name.lower(),
                "category": "samples",
                "description": description,
                "type": field_type,
            }

    def get_samples(self):
        """Return list of samples (individual ids)."""
        return self.samples

    def get_metadatas(self):
        """override from VCF"""
        output = {"filename": self.filename}

        for key, value in self.vcf_reader.metadata.items():
            output[key] = str(value)

        self.metadatas = output

        return output

    def get_infos(self):
        """override from VCF"""

        self.infos = self.vcf_reader.infos

        return self.infos

    def get_formats(self):
        """override from VCF"""

        self.formats = self.vcf_reader.formats

        return self.formats

    def get_contigs(self):
        """override from VCF"""

        self.contigs = self.vcf_reader.contigs

        return self.contigs

    def get_json(self):

        variants_position = (
            self.options.get("Variants", {}).get("rings", {}).get("position", 0.5)
        )
        variants_ring_height = (
            self.options.get("Variants", {}).get("rings", {}).get("height", 0.04)
        )
        variants_ring_space = (
            self.options.get("Variants", {}).get("rings", {}).get("space", 0.01)
        )

        # default json
        default_json = {
            "General": {"width": 1500, "height": 1500, "title": None},
            "Category": {
                "ideogram": {
                    "patch": {
                        "file": {
                            "path": "",
                            "header": "infer",
                            "sep": "\t",
                            "dataframe": {
                                "orient": "columns",
                                "data": {
                                    "chr_name": ["chr1", "chr2", "chr3"],
                                    "chr_size": [249250621, 243199373, 198022430],
                                    "chr_label": ["chr1", "chr2", "chr3"],
                                    "chr_color": ["pink", "rosybrown", "firebrick"],
                                },
                            },
                        },
                        "show": "True",
                        "degreerange": [0, 360],
                        "showfillcolor": "True",
                        "chrannotation": {
                            "show": "True",
                            "radius": {"R": 1.25},
                            "fonttype": "bold",
                            "textangle": {"angleoffset": 0, "anglelimit": 360},
                            "layout": {
                                "xref": "x",
                                "yref": "y",
                                "showarrow": False,
                                "font": {"size": 10, "color": "black"},
                            },
                        },
                        "customoptions": {
                            "customlabel": "True",
                            "customspacing": "False",
                            "customcolor": 3,
                        },
                        "npoints": 1000,
                        "radius": {"R0": 1.0, "R1": 1.1},
                        "layout": {
                            "type": "path",
                            "opacity": 0.9,
                            "layer": "above",
                            "line": {"color": "gray", "width": 2},
                        },
                    },
                    "majortick": {
                        "show": "True",
                        "spacing": 30000000,
                        "radius": {"R0": 1.1, "R1": 1.125},
                        "layout": {
                            "type": "path",
                            "opacity": 0.9,
                            "layer": "above",
                            "line": {"color": "black", "width": 1},
                        },
                    },
                    "minortick": {
                        "show": "True",
                        "spacing": 5000000,
                        "radius": {"R0": 1.1, "R1": 1.118},
                        "layout": {
                            "type": "path",
                            "opacity": 0.9,
                            "line": {"color": "black", "width": 0.5},
                        },
                    },
                    "ticklabel": {
                        "show": "True",
                        "spacing": 30000000,
                        "radius": {"R": 1.16},
                        "textformat": "Mb",
                        "textangle": {"angleoffset": -90, "anglelimit": 360},
                        "layout": {
                            "xref": "x",
                            "yref": "y",
                            "showarrow": False,
                            "font": {
                                "family": "Times New Roman",
                                "size": 8,
                                "color": "black",
                            },
                        },
                    },
                },
                "ring": [
                    # {
                    #     # genes
                    #     "radius": {
                    #         "R0": 0.90,
                    #         "R1": 0.99
                    #     },
                    #     "layout": {
                    #         "opacity": 0.0,
                    #         "fillcolor": "white",
                    #         "layer": "below",
                    #         "line": {
                    #             "color": "gray",
                    #             "width": 0
                    #         }
                    #     }
                    # },
                    # {
                    #     # exons
                    #     "radius": {
                    #         "R0": 0.90,
                    #         "R1": 0.99
                    #     },
                    #     "layout": {
                    #         "opacity": 0.0,
                    #         "fillcolor": "white",
                    #         "layer": "below",
                    #         "line": {
                    #             "color": "gray",
                    #             "width": 0
                    #         }
                    #     }
                    # },
                    {
                        # SNV
                        "radius": {
                            "R0": variants_position
                            + (7 * variants_ring_space)
                            + (6 * variants_ring_height),
                            "R1": variants_position
                            + (7 * variants_ring_space)
                            + (7 * variants_ring_height),
                        },
                        "layout": {
                            "opacity": 0.1,
                            "fillcolor": "gray",
                            "layer": "below",
                            "line": {"color": "gray", "width": 1},
                        },
                    },
                    {
                        # level 5
                        "radius": {
                            "R0": variants_position
                            + (6 * variants_ring_space)
                            + (5 * variants_ring_height),
                            "R1": variants_position
                            + (6 * variants_ring_space)
                            + (6 * variants_ring_height),
                        },
                        "layout": {
                            "opacity": 0.1,
                            "fillcolor": "lightgrey",
                            "layer": "below",
                            "line": {"color": "lightgrey", "width": 1},
                        },
                    },
                    {
                        # level 4
                        "radius": {
                            "R0": variants_position
                            + (5 * variants_ring_space)
                            + (4 * variants_ring_height),
                            "R1": variants_position
                            + (5 * variants_ring_space)
                            + (5 * variants_ring_height),
                        },
                        "layout": {
                            "opacity": 0.1,
                            "fillcolor": "lightgrey",
                            "layer": "below",
                            "line": {"color": "lightgrey", "width": 1},
                        },
                    },
                    {
                        # level 3
                        "radius": {
                            "R0": variants_position
                            + (4 * variants_ring_space)
                            + (3 * variants_ring_height),
                            "R1": variants_position
                            + (4 * variants_ring_space)
                            + (4 * variants_ring_height),
                        },
                        "layout": {
                            "opacity": 0.1,
                            "fillcolor": "lightgrey",
                            "layer": "below",
                            "line": {"color": "lightgrey", "width": 1},
                        },
                    },
                    {
                        # level 2
                        "radius": {
                            "R0": variants_position
                            + (3 * variants_ring_space)
                            + (2 * variants_ring_height),
                            "R1": variants_position
                            + (3 * variants_ring_space)
                            + (3 * variants_ring_height),
                        },
                        "layout": {
                            "opacity": 0.1,
                            "fillcolor": "white",
                            "layer": "below",
                            "line": {"color": "white", "width": 1},
                        },
                    },
                    {
                        # level 1
                        "radius": {
                            "R0": variants_position
                            + (2 * variants_ring_space)
                            + (1 * variants_ring_height),
                            "R1": variants_position
                            + (2 * variants_ring_space)
                            + (2 * variants_ring_height),
                        },
                        "layout": {
                            "opacity": 0.1,
                            "fillcolor": "lightgrey",
                            "layer": "below",
                            "line": {"color": "lightgrey", "width": 1},
                        },
                    },
                    {
                        # level 0
                        "radius": {
                            "R0": variants_position
                            + (1 * variants_ring_space)
                            + (0 * variants_ring_height),
                            "R1": variants_position
                            + (1 * variants_ring_space)
                            + (1 * variants_ring_height),
                        },
                        "layout": {
                            "opacity": 0.1,
                            "fillcolor": "lightgrey",
                            "layer": "below",
                            "line": {"color": "lightgrey", "width": 1},
                        },
                    },
                ],
            },
        }

        ### params from default
        params = copy.deepcopy(default_json)

        # General
        if self.options.get("General", None):
            params["General"] = self.options.get("General", None)

        # Ideogram

        # Cytoband
        if self.options.get("Chromosomes", {}).get("cytoband", None):

            # Find cytoband data from options json dataframe json
            self.options["Cytoband"] = file_to_dict(
                self.options.get("Chromosomes", {}).get("cytoband", None),
                self.options.get("File", ""),
            )

            # Explode data in dataframe
            self.options["Cytoband"] = explode_category_file_dict_into_dataframe(
                self.options["Cytoband"], self.options.get("File", "")
            )

            # Cytoband params
            cytoband = {
                "show": "True",
                "file": self.options["Cytoband"],
                "sortbycolor": "True",
                "colorcolumn": 3,
                "hovertextformat": ' "<b>{}</b>".format(a[i,0])',
                "trace": {
                    "uid": "cytoband",
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {
                        "size": 1,
                        "symbol": 0,  # 8
                        "color": "black",
                        "opacity": 1,
                    },
                },
                "layout": {
                    "type": "path",
                    "layer": "below",
                    "opacity": 1.0,
                    "line": {"color": "black", "width": 0},
                },
            }

            params["Category"]["cytoband"] = cytoband

            # Cytoband infos
            self.options["Cytoband_infos"] = copy.deepcopy(self.options["Cytoband"])

            cytoband_data = {
                "chr_name": self.options["Cytoband"]
                .get("dataframe", {})
                .get("data", {})
                .get("chr_name", []),
                "start": self.options["Cytoband"]
                .get("dataframe", {})
                .get("data", {})
                .get("start", []),
                "end": self.options["Cytoband"]
                .get("dataframe", {})
                .get("data", {})
                .get("end", []),
                "val": [1]
                * len(
                    self.options["Cytoband"]
                    .get("dataframe", {})
                    .get("data", {})
                    .get("chr_name", [])
                ),
                "band_color": ["lightgray"]
                * len(
                    self.options["Cytoband"]
                    .get("dataframe", {})
                    .get("data", {})
                    .get("chr_name", [])
                ),
                "band": self.options["Cytoband"]
                .get("dataframe", {})
                .get("data", {})
                .get(
                    "band",
                    [""]
                    * len(
                        self.options["Cytoband"]
                        .get("dataframe", {})
                        .get("data", {})
                        .get("chr_name", [])
                    ),
                ),
            }
            self.options["Cytoband_infos"]["dataframe"]["data"] = cytoband_data

            cytoband_infos = {
                "show": "True",
                "file": self.options["Cytoband_infos"],
                "colorcolumn": 4,
                "radius": {"R0": 1, "R1": 1.1},
                "hovertextformat": " \"<b>{}:{}-{}<br>{}{}</b>\".format(a[i,0], a[i,1], a[i,2], a[i,0].replace('chr', ''), ''.join(a[i,5:]))",
                # "hovertextformat": " \"<b>{}</b>\".format(a[i,0])",
                "trace": {
                    "uid": "cytoband_tile",
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {
                        "size": 0,
                        "symbol": 0,  # 8
                        "color": self.options["Cytoband_infos"]["dataframe"]["data"][
                            "band_color"
                        ],
                        "opacity": 0,
                    },
                    "hovertextformat": " \"<b>{}:{}-{}<br>{}{}</b>\".format(a[i,0], a[i,1], a[i,2], a[i,0].replace('chr', ''), ''.join(a[i,5:]))",
                    "trace": {
                        "uid": "cytoband_tile",
                        "hoverinfo": "text",
                        "mode": "markers",
                        "marker": {
                            "size": 0,
                            "symbol": 0,  # 8
                            "color": self.options["Cytoband_infos"]["dataframe"][
                                "data"
                            ]["band_color"],
                            "opacity": 0,
                        },
                    },
                },
            }

            cytoband_category_type = "histogram"
            if cytoband_category_type not in params["Category"]:
                params["Category"][cytoband_category_type] = []
            params["Category"][cytoband_category_type].append(cytoband_infos)

        # Genes
        if self.options.get("Genes", {}).get("data", None):

            # Find genes data from options json dataframe json
            self.options["Genes"]["data"] = file_to_dict(
                self.options.get("Genes", {}).get("data", None),
                self.options.get("File", ""),
            )

            # Explode data in dataframe
            self.options["Genes"]["data"] = explode_category_file_dict_into_dataframe(
                self.options.get("Genes", {}).get("data", None),
                self.options.get("File", ""),
            )

            # Store full refGene
            self.options["Genes_full"] = copy.deepcopy(self.options["Genes"]["data"])

        # Exons

        if self.options.get("Exons", None):

            # Find exons data from options json dataframe json

            self.options["Exons"]["data"] = file_to_dict(
                self.options.get("Exons", {}).get("data", None),
                self.options.get("File", ""),
            )

            # Explode data in dataframe
            self.options["Exons"]["data"] = explode_category_file_dict_into_dataframe(
                self.options.get("Exons", {}).get("data", None),
                self.options.get("File", ""),
            )

            # Store full refGene exon
            # self.options["Exons_full"] = copy.deepcopy(self.options["Exons"]["data"])

        # Additonal_annotations

        if self.options.get("Additonal_annotations", None):

            categories = []
            if isinstance(self.options.get("Additonal_annotations", None), str):
                self.options["Additonal_annotations"] = [
                    self.options.get("Additonal_annotations", None)
                ]
            if isinstance(self.options.get("Additonal_annotations", None), list):
                for category in self.options.get("Additonal_annotations", []):
                    categories.append(
                        file_to_dict(category, self.options.get("File", ""))
                    )
            else:
                print(
                    "[WARN] Additonal_annotations options not well formed. Use JSON file or JSON dict formed as a 'category'  in plotly"
                )

            for category in categories:
                for category_type in category:
                    if category_type not in params["Category"]:
                        params["Category"][category_type] = []
                    for category_dict in category[category_type]:
                        category_file_dict = explode_category_file_dict_into_dataframe(
                            category_dict.get("file", {}), self.options.get("File", "")
                        )
                        category_dict["file"] = category_file_dict
                        params["Category"][category_type].append(
                            copy.deepcopy(category_dict)
                        )

        ### Calculations

        # Font size
        # change font size depending on width and height of general params

        font_size = params["General"].get("height", 1000) / 100 or 10
        if font_size < 3:
            font_size = 3
        params["Category"]["ideogram"]["patch"]["chrannotation"]["layout"]["font"][
            "size"
        ] = font_size
        params["Category"]["ideogram"]["ticklabel"]["layout"]["font"]["size"] = (
            font_size - 2
        )

        ### create data for each variant type

        # Category structure and pattern
        categories = {
            "link": {
                "pattern": {
                    "show": "True",
                    "file": {
                        "path": "",
                        "header": "infer",
                        "sep": "\t",
                        "dataframe": {"orient": "columns", "data": {}},
                    },
                    "radius": {"R0": 0, "R1": 0.49},
                    "sortbycolor": "False",
                    "colorcolumn": 6,
                    # "hovertextformat": [
                    #     " \"<b>{}:{}-{}<br>{}:{}-{}</b><br><br>{}\".format(a[i,0], a[i,1], a[i,2], a[i,3], a[i,4], a[i,5], a[i,7].replace(', ', '<br> ').translate(str.maketrans({',': '<br>    ', '{': '{<br> ', '}': '<br>}'})))",
                    #     " \"<b>{}:{}-{}<br>{}:{}-{}</b><br><br>{}\".format(a[i,3], a[i,4], a[i,5], a[i,0], a[i,1], a[i,2], a[i,7].replace(', ', '<br> ').translate(str.maketrans({',': '<br>    ', '{': '{<br> ', '}': '<br>}'})))"
                    # ],
                    "hovertextformat": [
                        ' "<b>{}:{}-{}<br>{}:{}-{}</b><br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,3], a[i,4], a[i,5], a[i,7])',
                        ' "<b>{}:{}-{}<br>{}:{}-{}</b><br><br>{}".format(a[i,3], a[i,4], a[i,5], a[i,0], a[i,1], a[i,2], a[i,7])',
                    ],
                    "trace": {
                        "uid": "unknown",
                        "hoverinfo": "text",
                        "mode": "markers",
                        "marker": {"size": 0, "symbol": 0, "opacity": 0},
                    },
                    "layout": {
                        "type": "path",
                        "layer": "above",
                        "opacity": 0.8,
                        "line": {"color": "gray", "width": 5},
                    },
                },
                "data": [],
            },
            "scatter": {
                "pattern": {
                    "show": "True",
                    "file": {
                        "path": "",
                        "header": "infer",
                        "sep": "\t",
                        "dataframe": {"orient": "columns", "data": {}},
                    },
                    "radius": {"R0": 0.69, "R1": 0.78},
                    "sortbycolor": "False",
                    "colorcolumn": 6,
                    # "hovertextformat": " \"<b>{}:{}</b> | {} > {}<br>{}<br><br>{}\".format(a[i,0], a[i,1], a[i,3], a[i,4], a[i,5], a[i,7].replace(', ', '<br> ').translate(str.maketrans({',': '<br>    ', '{': '{<br> ', '}': '<br>}'})))",
                    "hovertextformat": ' "<b>{}:{}</b> | {} > {}<br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,3], a[i,4], a[i,5], a[i,7])',
                    "trace": {
                        "hoverinfo": "text",
                        "mode": "markers",
                        "opacity": 1,
                        "marker": {"size": 5, "symbol": 1, "opacity": 1},
                    },
                },
                "data": [],
            },
            "tile": {
                "pattern": {
                    "show": "True",
                    "file": {
                        "path": "",
                        "header": "infer",
                        "sep": "\t",
                        "dataframe": {"orient": "columns", "data": {}},
                    },
                    "colorcolumn": "ideogram",
                    "radius": {"R0": 0.6, "R1": 0.68, "min": 1, "max": 5},
                    # "hovertextformat": " \"<b>{}:{}-{}</b><br>{}<br><br>{}\".format(a[i,0], a[i,1], a[i,2], a[i,6], a[i,8].replace(', ', '<br> ').translate(str.maketrans({',': '<br>    ', '{': '{<br> ', '}': '<br>}'})))",
                    "hovertextformat": ' "<b>{}:{}-{}</b><br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,6], a[i,8])',
                    "trace": {
                        "hoverinfo": "text",
                        "mode": "markers",
                        "marker": {
                            "size": 5,
                            "symbol": 1,
                            "color": "red",
                            "opacity": 1,
                        },
                    },
                    "layout": {
                        "type": "path",  # ['circle', 'rect', 'path', 'line']
                        "layer": "above",  # above below
                        "opacity": 1,
                        "line": {"color": "gray", "width": 5},
                    },
                },
                "data": [],
            },
            "histogram": {
                "pattern": {
                    "show": "True",
                    "customfillcolor": "True",
                    "file": {
                        "path": "",
                        "header": "infer",
                        "sep": "\t",
                        "dataframe": {"orient": "columns", "data": {}},
                    },
                    # "sortbycolor": "True",
                    "colorcolumn": 7,
                    "radius": {"R0": 0.945, "R1": 0.99},
                    "hovertextformat": ' "<b>{}:{}-{}</b><br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,6], a[i,8])',
                    "trace": {
                        "hoverinfo": "text",
                        "mode": "markers",
                        "marker": {"size": 5, "color": "lightgray", "opacity": 0.1},
                    },
                    "layout": {
                        "type": "path",
                        "layer": "above",  # above below
                        "opacity": 0.1,
                        "fillcolor": "red",
                        "line": {"color": "lightgray", "width": 5},
                    },
                },
                "data": [],
            },
            "line": {
                "pattern": {
                    "show": "True",
                    "file": {
                        "path": "",
                        "header": "infer",
                        "sep": "\t",
                        "dataframe": {"orient": "columns", "data": {}},
                    },
                    "radius": {"R0": 0.85, "R1": 0.88},
                    "sortbycolor": "False",
                    "colorcolumn": "ideogram",
                    "hovertextformat": ' "Chromosome: {}<br>Position: {}<br>Value: {:.2f}".format(a[i,0], a[i,1], float(a[i,2]))',
                    "trace": {
                        "hoverinfo": "text",
                        "mode": "lines+markers",
                        "opacity": 0.9,
                        "marker": {"symbol": 0, "size": 3, "color": "black"},
                        "line": {
                            "color": "black",
                            "width": 2,
                            "shape": "linear",
                            "smoothing": 0,
                        },
                    },
                },
                "data": [],
            },
        }

        # Data structure variant type
        variants_data = {
            "breakpoint": {
                "category": "link",
                "nb": 0,
                "data": {
                    "chr1_name": [],
                    "chr1_start": [],
                    "chr1_end": [],
                    "chr_name": [],
                    "chr2_start": [],
                    "chr2_end": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {"R0": 0, "R1": variants_position},
            },
            "snv": {
                "category": "scatter",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (7 * variants_ring_space)
                    + (6.5 * variants_ring_height),
                    "R1": variants_position
                    + (7 * variants_ring_space)
                    + (6.5 * variants_ring_height),
                },
            },
            "cnv_level_5": {
                "category": "histogram",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (6 * variants_ring_space)
                    + (5.5 * variants_ring_height),
                    "R1": variants_position
                    + (6 * variants_ring_space)
                    + (5.5 * variants_ring_height),
                },
            },
            "cnv_scatter_level_5": {
                "category": "scatter",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (6 * variants_ring_space)
                    + (5.5 * variants_ring_height),
                    "R1": variants_position
                    + (6 * variants_ring_space)
                    + (5.5 * variants_ring_height),
                },
            },
            "cnv_level_4": {
                "category": "histogram",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (5 * variants_ring_space)
                    + (4.5 * variants_ring_height),
                    "R1": variants_position
                    + (5 * variants_ring_space)
                    + (4.5 * variants_ring_height),
                },
            },
            "cnv_scatter_level_4": {
                "category": "scatter",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (5 * variants_ring_space)
                    + (4.5 * variants_ring_height),
                    "R1": variants_position
                    + (5 * variants_ring_space)
                    + (4.5 * variants_ring_height),
                },
            },
            "cnv_level_3": {
                "category": "histogram",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (4 * variants_ring_space)
                    + (3.5 * variants_ring_height),
                    "R1": variants_position
                    + (4 * variants_ring_space)
                    + (3.5 * variants_ring_height),
                },
            },
            "cnv_scatter_level_3": {
                "category": "scatter",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (4 * variants_ring_space)
                    + (3.5 * variants_ring_height),
                    "R1": variants_position
                    + (4 * variants_ring_space)
                    + (3.5 * variants_ring_height),
                },
            },
            "cnv_level_2": {
                "category": "histogram",  # tile histogram
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (3 * variants_ring_space)
                    + (2.5 * variants_ring_height),
                    "R1": variants_position
                    + (3 * variants_ring_space)
                    + (2.5 * variants_ring_height),
                },
            },
            "cnv_scatter_level_2": {
                "category": "scatter",  # tile histogram
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (3 * variants_ring_space)
                    + (2.5 * variants_ring_height),
                    "R1": variants_position
                    + (3 * variants_ring_space)
                    + (2.5 * variants_ring_height),
                },
            },
            "cnv_level_1": {
                "category": "histogram",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (2 * variants_ring_space)
                    + (1.5 * variants_ring_height),
                    "R1": variants_position
                    + (2 * variants_ring_space)
                    + (1.5 * variants_ring_height),
                },
            },
            "cnv_scatter_level_1": {
                "category": "scatter",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (2 * variants_ring_space)
                    + (1.5 * variants_ring_height),
                    "R1": variants_position
                    + (2 * variants_ring_space)
                    + (1.5 * variants_ring_height),
                },
            },
            "cnv_level_0": {
                "category": "histogram",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (1 * variants_ring_space)
                    + (0.5 * variants_ring_height),
                    "R1": variants_position
                    + (1 * variants_ring_space)
                    + (0.5 * variants_ring_height),
                },
            },
            "cnv_scatter_level_0": {
                "category": "scatter",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position
                    + (1 * variants_ring_space)
                    + (0.5 * variants_ring_height),
                    "R1": variants_position
                    + (1 * variants_ring_space)
                    + (0.5 * variants_ring_height),
                },
            },
            "cnv": {
                "category": "histogram",
                "nb": 0,
                "data": {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "ref": [],
                    "alt": [],
                    "type": [],
                    "color": [],
                    "hovertext": [],
                    "symbol": [],
                    "genes": [],
                    "exons": [],
                },
                "radius": {
                    "R0": variants_position + 0.30,
                    "R1": variants_position + 0.00,
                },
            },
        }
        # variants_position + (1 * variants_ring_space) + (1 * variants_ring_height)
        genes_target = {}
        genes_target_type = {"BND": {}, "CNV": {}, "SNV": {}}
        exons_target = {}
        exons_target_type = {"BND": {}, "CNV": {}, "SNV": {}}

        ### fetch variant
        if True:
            for i, record in enumerate(self.parse_variants(add_sample=True)):

                record_data = {}

                if record.get("chr", "") in self.options.get("Chromosomes", {}).get(
                    "list", []
                ) or not self.options.get("Chromosomes", {}).get("list", []):

                    annotations_fields_to_show = {}
                    for annotation in self.options["Variants"]["annotations"]["fields"]:
                        if (
                            str(annotation).lower() != ""
                            and (
                                record.get(str(annotation).lower(), None)
                                or self.options["Variants"]["annotations"]["show_none"]
                            )
                            and str(annotation).lower() != "*"
                        ):
                            annotations_fields_to_show[annotation] = record.get(
                                str(annotation).lower(), None
                            )
                        elif annotation == "*":
                            for annotation_hovertext in record:
                                if (
                                    annotation_hovertext != ""
                                    and record.get(annotation_hovertext, None)
                                    and annotation_hovertext
                                    not in self.options["Variants"]["annotations"][
                                        "fields"
                                    ]
                                    and (
                                        record.get(annotation_hovertext, None)
                                        or self.options["Variants"]["annotations"][
                                            "show_none"
                                        ]
                                    )
                                ):
                                    annotations_fields_to_show[
                                        annotation_hovertext
                                    ] = record.get(annotation_hovertext, None)

                    annotations_fields_to_show_text = "<br>".join(
                        "{}: {}".format(key, value)
                        for key, value in annotations_fields_to_show.items()
                        if not isinstance(value, list)
                    )
                    record_data["hovertext"] = ""

                    record_types = []

                    if "svtype" in record:

                        # for svtype = BND, try to find if it's another CNV with ALT
                        if record.get("svtype", "") in ["CNV", "BND"]:
                            CNV_type = re.sub(
                                r"[^a-zA-Z]", "", record.get("alt", "").split(":")[0]
                            )
                            if CNV_type in ["INS", "INV", "DEL", "DUP", "TRA"]:
                                record["svtype"] = CNV_type
                        if record.get("svtype", "") in ["BND", "TRA"]:
                            record_data["chr1_name"] = record.get("chr", "")
                            record_data["chr1_start"] = record.get("pos", 0)
                            record_data["chr1_end"] = record.get("pos", 0)
                            find_point2 = "".join(
                                [
                                    s
                                    for s in record.get("alt", "")
                                    if s.isdigit() or s == ":" or s == "X" or s == "Y"
                                ]
                            ).split(":")
                            if len(find_point2) < 2:
                                if record.get(
                                    "chr2", record.get("chr", None)
                                ) and record.get("end", None):
                                    find_point2 = [
                                        re.sub(
                                            "chr",
                                            "",
                                            record.get("chr2", record.get("chr", None)),
                                        ),
                                        record.get("end", None),
                                    ]
                            record_data["chr_name"] = "chr" + find_point2[0]
                            record_data["chr2_start"] = int(find_point2[1])
                            record_data["chr2_end"] = int(find_point2[1])
                            record_data["color"] = "blue"
                            record_data["symbol"] = 0

                            record_data["genes"] = find_genes_for_record(
                                record_data,
                                self.options.get("Genes_full", {})
                                .get("dataframe", {})
                                .get("data", {}),
                            )
                            genes_target.update(record_data["genes"])
                            genes_target_type["BND"].update(record_data["genes"])

                            record_data["exons"] = find_genes_for_record(
                                record_data,
                                self.options.get("Exons_full", {})
                                .get("dataframe", {})
                                .get("data", {}),
                            )
                            exons_target.update(record_data["exons"])
                            exons_target_type["BND"].update(record_data["exons"])

                            # Add genes and exons info on hovertext
                            if (
                                self.options.get("Genes_full", {})
                                .get("dataframe", {})
                                .get("data", {})
                            ):
                                record_data["hovertext"] += (
                                    "Genes ("
                                    + str(len(record_data["genes"].keys()))
                                    + "): "
                                    + textwrap.shorten(
                                        ", ".join(record_data["genes"].keys()),
                                        width=TEXTWRAP_WIDTH,
                                        placeholder=TEXTWRAP_HOLDER,
                                    )
                                    + "<br>"
                                )

                            if (
                                self.options.get("Exons_full", {})
                                .get("dataframe", {})
                                .get("data", {})
                            ):
                                record_data["hovertext"] += (
                                    "Exons ("
                                    + str(len(record_data["exons"].keys()))
                                    + "): "
                                    + textwrap.shorten(
                                        ", ".join(record_data["exons"].keys()),
                                        width=TEXTWRAP_WIDTH,
                                        placeholder=TEXTWRAP_HOLDER,
                                    )
                                    + "<br>"
                                )

                            if record_data["hovertext"] != "":
                                record_data["hovertext"] += "<br>"

                            record_data["hovertext"] += annotations_fields_to_show_text

                            record_types.append({"breakpoint": record_data})

                        # CNV: INS INV DEL DUP
                        elif record.get("svtype", "") in [
                            "CNV",
                            "INS",
                            "INV",
                            "DEL",
                            "DUP",
                        ]:
                            # for svtype = CNV, try to find correct CNV with ALT
                            if record.get("svtype", "") in ["CNV"]:
                                CNV_type = re.sub(
                                    r"[^a-zA-Z]",
                                    "",
                                    record.get("alt", "").split(":")[0],
                                )
                                if CNV_type in ["CNV", "INS", "INV", "DEL", "DUP"]:
                                    record["svtype"] = CNV_type

                            record_data["chr_name"] = record.get("chr", "")
                            record_data["start"] = record.get("pos", 0)
                            if record.get("end", 0):
                                record_data["end"] = int(record.get("end", 0))
                            elif record.get("svlen", 0):
                                record_data["end"] = record_data["start"] + record.get(
                                    "svlen", 0
                                )
                            record_data["val"] = 2
                            level = 2
                            record_data["ref"] = record.get("ref", "")
                            record_data["alt"] = record.get("alt", "")
                            record_data["color"] = "gray"
                            record_data["type"] = "unknown"
                            if record.get("svtype", "") in ["INS"]:
                                record_data["type"] = "INS"
                                record_data["color"] = "purple"
                                record_data["val"] = 2
                            elif record.get("svtype", "") in ["INV"]:
                                record_data["type"] = "INV"
                                record_data["color"] = "orange"
                                record_data["val"] = 2
                            elif record.get("svtype", "") in ["DEL"]:
                                record_data["type"] = "DEL"
                                record_data["color"] = "red"
                                record_data["val"] = 2
                                level = 1
                            elif record.get("svtype", "") in ["DUP"]:
                                record_data["type"] = "DUP"
                                record_data["color"] = "blue"
                                record_data["val"] = 2
                            elif record.get("svtype", "") in ["CNV"]:
                                record_data["type"] = "CNV"
                                record_data["color"] = "blue"
                                record_data["val"] = 2

                            # Find copy number
                            copy_number = None
                            if re.sub("[^0-9]", "", record["alt"]) != "":
                                copy_number = int(
                                    re.sub("[^0-9]", "", record["alt"]) or "2"
                                )
                            elif (
                                re.sub("[^0-9]", "", str(record.get("cn", None))) != ""
                            ):
                                copy_number = int(
                                    re.sub("[^0-9]", "", str(record.get("cn", None)))
                                    or "2"
                                )
                            else:
                                for s in record.get("samples", []):
                                    cn = s.get("cn", None)
                                    if cn is not None:
                                        copy_number = cn
                            if copy_number is not None and copy_number >= 0:
                                level = copy_number
                            if level > 5:
                                level = 5

                            record_data["val"] = 2

                            record_data["symbol"] = 0

                            record_data["genes"] = find_genes_for_record(
                                record_data,
                                self.options.get("Genes_full", {})
                                .get("dataframe", {})
                                .get("data", {}),
                            )
                            genes_target.update(record_data["genes"])
                            genes_target_type["CNV"].update(record_data["genes"])

                            record_data["exons"] = find_genes_for_record(
                                record_data,
                                self.options.get("Exons_full", {})
                                .get("dataframe", {})
                                .get("data", {}),
                            )
                            exons_target.update(record_data["exons"])
                            exons_target_type["CNV"].update(record_data["exons"])

                            # Add genes and exons info on hovertext
                            if (
                                self.options.get("Genes_full", {})
                                .get("dataframe", {})
                                .get("data", {})
                            ):
                                record_data["hovertext"] += (
                                    "Genes ("
                                    + str(len(record_data["genes"].keys()))
                                    + "): "
                                    + textwrap.shorten(
                                        ", ".join(record_data["genes"].keys()),
                                        width=TEXTWRAP_WIDTH,
                                        placeholder=TEXTWRAP_HOLDER,
                                    )
                                    + "<br>"
                                )

                            if (
                                self.options.get("Exons_full", {})
                                .get("dataframe", {})
                                .get("data", {})
                            ):
                                record_data["hovertext"] += (
                                    "Exons ("
                                    + str(len(record_data["exons"].keys()))
                                    + "): "
                                    + textwrap.shorten(
                                        ", ".join(record_data["exons"].keys()),
                                        width=TEXTWRAP_WIDTH,
                                        placeholder=TEXTWRAP_HOLDER,
                                    )
                                    + "<br>"
                                )

                            if record_data["hovertext"] != "":
                                record_data["hovertext"] += "<br>"

                            record_data["hovertext"] += annotations_fields_to_show_text

                            # scatter
                            record_data_scatter_start = record_data.copy()
                            record_data_scatter_end = record_data.copy()
                            record_data_scatter_start.pop("end", None)
                            record_data_scatter_end.pop("end", None)
                            record_data_scatter_end["start"] = record_data["end"]

                            record_data["color"] = "gray"

                            if level >= 0:
                                cnv_level = "cnv_level_" + str(level)
                                cnv_scatter_level = "cnv_scatter_level_" + str(level)
                                record_types.append({cnv_level: record_data})
                                record_types.append(
                                    {cnv_scatter_level: record_data_scatter_start}
                                )
                                record_types.append(
                                    {cnv_scatter_level: record_data_scatter_end}
                                )

                    else:

                        # SNV InDel Multiallele
                        if True:
                            # record_data = {}
                            record_data["chr_name"] = record.get("chr", "")
                            record_data["start"] = record.get("pos", 0)
                            record_data["val"] = 1
                            record_data["ref"] = record.get("ref", "")
                            record_data["alt"] = record.get("alt", "")
                            type = "unknown"
                            color = "gray"
                            symbol = 0
                            if (
                                len(record.get("ref", "")) == 1
                                and len(record.get("alt", "")) == 1
                            ):
                                type = "SNV"
                                color = "blue"
                                symbol = 0
                            elif "," in record.get("ref", "") or "," in record.get(
                                "alt", ""
                            ):
                                type = "MultiAllele"
                                color = "gray"
                                symbol = 0
                            elif len(record.get("ref", "")) == len(
                                record.get("alt", "")
                            ):
                                type = "MNV"
                                color = "blue"
                                symbol = 0
                            elif len(record.get("ref", "")) > len(
                                record.get("alt", "")
                            ):
                                type = "DEL"
                                color = "red"
                                symbol = 0
                            elif len(record.get("ref", "")) < len(
                                record.get("alt", "")
                            ):
                                type = "INS"
                                color = "purple"
                                symbol = 0
                            record_data["type"] = type
                            record_data["color"] = color

                            record_data["symbol"] = symbol

                            record_data["genes"] = find_genes_for_record(
                                record_data,
                                self.options.get("Genes_full", {})
                                .get("dataframe", {})
                                .get("data", {}),
                            )
                            genes_target.update(record_data["genes"])
                            genes_target_type["SNV"].update(record_data["genes"])

                            record_data["exons"] = find_genes_for_record(
                                record_data,
                                self.options.get("Exons_full", {})
                                .get("dataframe", {})
                                .get("data", {}),
                            )
                            exons_target.update(record_data["exons"])
                            exons_target_type["SNV"].update(record_data["exons"])

                            # Add genes and exons info on hovertext
                            if (
                                self.options.get("Genes_full", {})
                                .get("dataframe", {})
                                .get("data", {})
                            ):
                                record_data["hovertext"] += (
                                    "Genes ("
                                    + str(len(record_data["genes"].keys()))
                                    + "): "
                                    + textwrap.shorten(
                                        ", ".join(record_data["genes"].keys()),
                                        width=TEXTWRAP_WIDTH,
                                        placeholder=TEXTWRAP_HOLDER,
                                    )
                                    + "<br>"
                                )

                            if (
                                self.options.get("Exons_full", {})
                                .get("dataframe", {})
                                .get("data", {})
                            ):
                                record_data["hovertext"] += (
                                    "Exons ("
                                    + str(len(record_data["exons"].keys()))
                                    + "): "
                                    + textwrap.shorten(
                                        ", ".join(record_data["exons"].keys()),
                                        width=TEXTWRAP_WIDTH,
                                        placeholder=TEXTWRAP_HOLDER,
                                    )
                                    + "<br>"
                                )

                            if record_data["hovertext"] != "":
                                record_data["hovertext"] += "<br>"

                            record_data["hovertext"] += annotations_fields_to_show_text

                            record_types.append({"snv": record_data})

                    for record_i in record_types:
                        for record_type in record_i:
                            for i in record_i[record_type]:
                                variants_data[record_type]["data"][i].append(
                                    record_i[record_type][i]
                                )
                            variants_data[record_type]["nb"] += 1

        # if self.options.get("Only_snv_in_sv_genes",False):
        if self.options.get("Genes", {}).get("only_snv_in_sv_genes", False):

            # find genes in BND + CNV and in SNV
            genes_to_keep_in_snv = []
            for gene in list(set(genes_target_type["BND"])) + list(
                set(genes_target_type["CNV"])
            ):
                if gene in list(set(genes_target_type["SNV"])):
                    genes_to_keep_in_snv.append(gene)

            # reconstruct snv data
            variant_data_snv_filtered = {}
            i = 0
            for snv_genes in variants_data["snv"]["data"]["genes"]:
                if not set(snv_genes.keys()).isdisjoint(genes_to_keep_in_snv):
                    for field in variants_data["snv"]["data"].keys():
                        if field in variant_data_snv_filtered:
                            variant_data_snv_filtered[field].append(
                                variants_data["snv"]["data"][field][i]
                            )
                        else:
                            variant_data_snv_filtered[field] = [
                                variants_data["snv"]["data"][field][i]
                            ]
                i += 1

            variants_data["snv"]["data"] = variant_data_snv_filtered

        ### Contigs

        # generate colors
        colors = list(
            Color("gray").range_to(Color("lightgrey"), len(self.get_contigs()))
        )
        # colors = list(Color("#333333",).range_to(Color("#ffffff"),len(self.get_contigs())))

        # create empty dataframe
        contig_dataframe = {}
        contig_dataframe["orient"] = "columns"
        contig_dataframe["data"] = {}
        contig_dataframe["data"]["chr_name"] = []
        contig_dataframe["data"]["chr_size"] = []
        contig_dataframe["data"]["chr_label"] = []
        contig_dataframe["data"]["chr_color"] = []

        # construct dataframe
        i = 0
        for contig in self.get_contigs():
            if contig in self.options.get("Chromosomes", {}).get(
                "list", []
            ) or not self.options.get("Chromosomes", {}).get("list", []):
                contig_dataframe["data"]["chr_name"].append(contig)
                contig_dataframe["data"]["chr_size"].append(
                    self.get_contigs()[contig].length
                )
                contig_dataframe["data"]["chr_label"].append(contig)
                contig_dataframe["data"]["chr_color"].append(str(colors[i]))
                i += 1

        # Add in params
        params["Category"]["ideogram"]["patch"]["file"]["dataframe"] = contig_dataframe

        # Genes

        if self.options.get("Genes", {}).get("data", None):

            gene_list_to_show = []
            if self.options.get("Genes", {}).get("list", None):
                gene_list_to_show = self.options.get("Genes", {}).get("list", None)
            if "targets" in gene_list_to_show:
                gene_list_to_show.extend(list(set(genes_target.keys())))

            # Gene List
            # filter with a list of genes

            if gene_list_to_show:
                gene_list = gene_list_to_show
                genes_data = {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "color": [],
                    "gene": [],
                }
                gene_i = 0
                # for gene in copy.deepcopy(self.options.get("Genes",{}).get("dataframe",{}).get("data",{})).get("gene",[]):
                for gene in copy.deepcopy(
                    self.options.get("Genes", {})
                    .get("data", {})
                    .get("dataframe", {})
                    .get("data", {})
                ).get("gene", []):
                    if gene in gene_list:
                        genes_data["chr_name"].append(
                            self.options["Genes"]["data"]["dataframe"]["data"][
                                "chr_name"
                            ][gene_i]
                        )
                        genes_data["start"].append(
                            self.options["Genes"]["data"]["dataframe"]["data"]["start"][
                                gene_i
                            ]
                        )
                        genes_data["end"].append(
                            self.options["Genes"]["data"]["dataframe"]["data"]["end"][
                                gene_i
                            ]
                        )
                        genes_data["val"].append(
                            self.options["Genes"]["data"]["dataframe"]["data"]["val"][
                                gene_i
                            ]
                        )
                        genes_data["color"].append(
                            self.options["Genes"]["data"]["dataframe"]["data"]["color"][
                                gene_i
                            ]
                        )
                        genes_data["gene"].append(
                            self.options["Genes"]["data"]["dataframe"]["data"]["gene"][
                                gene_i
                            ]
                        )
                    gene_i += 1
                self.options["Genes"]["data"]["dataframe"]["data"] = genes_data

            # Genes params

            genes = {
                "show": "True",
                "file": self.options["Genes"]["data"],
                "colorcolumn": 4,
                "radius": {"R0": 0.98, "R1": 0.98},
                "hovertextformat": ' "<b>{}:{}-{}<br>Gene: {}</b>".format(a[i,0], a[i,1], a[i,2], a[i,5])',
                "trace": {
                    "uid": "genes",
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {"size": 3, "symbol": 0, "color": "gray", "opacity": 1},
                },
                "layout": {
                    "type": "path",
                    "layer": "above",  # above below
                    "opacity": 1,
                    "line": {
                        "color": self.options["Genes"]["data"]["dataframe"]["data"][
                            "color"
                        ],
                        "width": 3,
                    },
                },
            }

            genes_category_type = "histogram"
            if genes_category_type not in params["Category"]:
                params["Category"][genes_category_type] = []
            params["Category"][genes_category_type].append(copy.deepcopy(genes))

            # Genes in scatter

            genes_scatter_start = copy.deepcopy(
                self.options["Genes"]["data"]["dataframe"]["data"]
            )
            genes_scatter_end = copy.deepcopy(
                self.options["Genes"]["data"]["dataframe"]["data"]
            )
            genes_scatter_start.pop("end", None)
            genes_scatter_end.pop("end", None)
            genes_scatter_end["start"] = self.options["Genes"]["data"]["dataframe"][
                "data"
            ]["end"]

            for i in genes_scatter_end:
                genes_scatter_start[i].extend(genes_scatter_end[i])

            genes_scatter = {
                "show": "True",
                "file": {
                    "path": "",
                    "header": "infer",
                    "sep": "\t",
                    "dataframe": {"orient": "columns", "data": genes_scatter_start},
                },
                "radius": {"R0": 0.98, "R1": 0.98},
                "colorcolumn": 3,
                "hovertextformat": ' "<b>{}:{}<br>Gene: {}</b>".format(a[i,0], a[i,1], a[i,4])',
                "trace": {
                    "hoverinfo": "text",
                    "mode": "markers",
                    "opacity": 1,
                    "marker": {"size": 5, "symbol": 0, "opacity": 1},
                },
            }

            genes_scatter_category_type = "scatter"
            if genes_scatter_category_type not in params["Category"]:
                params["Category"][genes_scatter_category_type] = []
            params["Category"][genes_scatter_category_type].append(
                copy.deepcopy(genes_scatter)
            )

        # Exons

        if self.options.get("Exons", None) and self.options.get("Exons", {}).get(
            "show", False
        ):
            gene_list_to_show = []
            if self.options.get("Genes", {}).get("list", None):
                gene_list_to_show = self.options.get("Genes", {}).get("list", None)
            if "targets" in gene_list_to_show:
                gene_list_to_show.extend(list(set(genes_target.keys())))

            # Gene List
            # filter with a list of exons

            if gene_list_to_show:
                gene_list = gene_list_to_show
                exons_data = {
                    "chr_name": [],
                    "start": [],
                    "end": [],
                    "val": [],
                    "color": [],
                    "gene": [],
                    "exon": [],
                }
                gene_i = 0
                # for gene in copy.deepcopy(self.options.get("Exons",{}).get("dataframe",{}).get("data",{})).get("gene"):
                for gene in copy.deepcopy(
                    self.options.get("Exons", {})
                    .get("data", {})
                    .get("dataframe", {})
                    .get("data", {})
                ).get("gene"):
                    if gene in gene_list:
                        exons_data["chr_name"].append(
                            self.options["Exons"]["data"]["dataframe"]["data"][
                                "chr_name"
                            ][gene_i]
                        )
                        exons_data["start"].append(
                            self.options["Exons"]["data"]["dataframe"]["data"]["start"][
                                gene_i
                            ]
                        )
                        exons_data["end"].append(
                            self.options["Exons"]["data"]["dataframe"]["data"]["end"][
                                gene_i
                            ]
                        )
                        exons_data["val"].append(
                            self.options["Exons"]["data"]["dataframe"]["data"]["val"][
                                gene_i
                            ]
                        )
                        exons_data["color"].append(
                            self.options["Exons"]["data"]["dataframe"]["data"]["color"][
                                gene_i
                            ]
                        )
                        exons_data["gene"].append(
                            self.options["Exons"]["data"]["dataframe"]["data"]["gene"][
                                gene_i
                            ]
                        )
                        exons_data["exon"].append(
                            self.options["Exons"]["data"]["dataframe"]["data"]["exon"][
                                gene_i
                            ]
                        )
                    gene_i += 1
                self.options["Exons"]["data"]["dataframe"]["data"] = exons_data

            # Exons params

            exons = {
                "show": "True",
                "file": self.options["Exons"]["data"],
                "colorcolumn": 4,
                "radius": {"R0": 0.96, "R1": 0.96},
                "hovertextformat": ' "<b>{}:{}-{}<br>Gene: {}<br>Exon: {}</b>".format(a[i,0], a[i,1], a[i,2], a[i,5], a[i,6])',
                "trace": {
                    "uid": "exons",
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {
                        "size": 3,
                        "symbol": 1,  # 8
                        "color": "gray",  # self.options["Exons"]["data"]["dataframe"]["data"]["color"],
                        "opacity": 1,
                    },
                },
                "layout": {
                    "type": "path",
                    "layer": "below",
                    "opacity": 1,
                    "line": {
                        "color": self.options["Exons"]["data"]["dataframe"]["data"][
                            "color"
                        ],
                        "width": 3,
                    },
                },
            }

            exons_category_type = "histogram"
            if exons_category_type not in params["Category"]:
                params["Category"][exons_category_type] = []
            params["Category"][exons_category_type].append(copy.deepcopy(exons))

        # Create param

        params_link = copy.deepcopy(params)
        variant_data_link = copy.deepcopy(variants_data)

        # for type in variants_data.keys():
        for type in variant_data_link.keys():

            # print(type)

            if variant_data_link[type]["nb"]:
                category_data = {}

                variant_data2 = copy.deepcopy(variant_data_link[type])
                category_data = copy.deepcopy(
                    categories[variant_data2["category"]]["pattern"]
                )

                if variant_data2["data"]:
                    category_data["trace"]["uid"] = type

                    category_data["file"]["dataframe"]["data"] = copy.deepcopy(
                        variant_data2["data"]
                    )
                    category_data["trace"]["marker"]["symbol"] = copy.deepcopy(
                        variant_data2["data"]["symbol"]
                    )

                    if variant_data_link[type]["category"] in ["tile", "histogram"]:
                        category_data["trace"]["marker"]["color"] = copy.deepcopy(
                            variant_data2["data"]["color"]
                        )[0]
                    else:
                        category_data["trace"]["marker"]["color"] = copy.deepcopy(
                            variant_data2["data"]["color"]
                        )

                    category_data["radius"] = copy.deepcopy(variant_data2["radius"])

                    if variant_data2["category"] not in params_link["Category"]:
                        params_link["Category"][variant_data2["category"]] = []
                    params_link["Category"][variant_data2["category"]].append(
                        copy.deepcopy(category_data)
                    )

        return params_link
