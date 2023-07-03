from genericpath import isfile
import json
import sys
import plotly.graph_objs as go
import plotly
import numpy as np
import pandas as pd
import maths
import colors
import copy
from Complex import Complex
from config import json_config, coord_config
import io
import os
import base64
from html import escape, unescape
from vcf2circos.utils import timeit, get_swap_dict
from ast import literal_eval
from scipy.spatial import KDTree
import inspect
from webcolors import CSS3_HEX_TO_NAMES, hex_to_rgb

pd.options.mode.chained_assignment = None  # default='warn'


# from dash_dict import *
def convert_rgb_to_names(rgb_tuple):
    # a dictionary of all the hex and their respective names in css3
    css3_db = CSS3_HEX_TO_NAMES
    names = []
    rgb_values = []
    for color_hex, color_name in css3_db.items():
        names.append(color_name)
        rgb_values.append(hex_to_rgb(color_hex))

    kdt_db = KDTree(rgb_values)
    distance, index = kdt_db.query(rgb_tuple)
    return names[index]


def process_legend(tr, trace_dict):
    if isinstance(tr["marker"]["color"], str):
        tr["marker"]["color"] = [tr["marker"]["color"]]
    for j, color in enumerate(tr["marker"]["color"]):
        # trace_dict.update(process_legend(color, tr, trace_dict))
        color = color_to_name(color)
        # For each item in dot
        for key, val in vars(tr)["_orphan_props"].items():
            # print(key, val)
            if isinstance(val, str):
                if key not in trace_dict[color]:
                    trace_dict[color][key] = val
            elif isinstance(val, np.ndarray):
                if key not in trace_dict[color]:
                    trace_dict[color][key] = []
                trace_dict[color][key].append(val[j])
            elif isinstance(val, list):
                if key not in trace_dict[color]:
                    trace_dict[color][key] = []
                trace_dict[color][key].append(val[j])
            elif isinstance(val, dict):
                if key not in trace_dict[color]:
                    trace_dict[color][key] = {}
                for k, v in val.items():
                    if isinstance(v, str) or isinstance(v, int) or isinstance(v, float):
                        if k not in trace_dict[color][key]:
                            trace_dict[color][key][k] = v
                    elif isinstance(v, np.ndarray):
                        if k not in trace_dict[color][key]:
                            trace_dict[color][key][k] = []
                        trace_dict[color][key][k].append(v[j])
                    elif isinstance(v, list):
                        if k not in trace_dict[color][key]:
                            trace_dict[color][key][k] = []
                        trace_dict[color][key][k].append(v[j])
        return trace_dict


# input_json_path = sys.argv[1]
def color_to_name(input):
    if isinstance(input, list) or isinstance(input, np.ndarray) or isinstance(input, tuple):
        l = []
        for colors in input:
            if colors.startswith("rgb"):
                # print(colors.replace("rgb", ""))
                col = literal_eval(colors.replace("rgb", ""))
                # print(colors)
                # print(col)
                l.append(convert_rgb_to_names(col))
            else:
                l.append(colors)
        # print(l)
        return l
    elif isinstance(input, str):
        if input.startswith("rgb"):
            col = literal_eval(input.replace("rgb", ""))
            return convert_rgb_to_names(col)
        else:
            return input
    else:
        return input


def merge_dict(basedict, *extradict):
    """
    this function updates basedict with extradict values, will apply to any nested values
    e.g:
    extradict = {'font': {'color': 'yellow'}}
    basedict = {'font': {'color': 'black', size: 16}}
    merge_dict(extradict, basedict):
        {'font': {'color': 'yellow', size: 16}}
    """
    for i in range(len(extradict)):
        for key in extradict[i].keys():
            if key not in basedict.keys():
                basedict[key] = extradict[i][key]
            else:
                if not isinstance(basedict[key], dict):
                    basedict[key] = extradict[i][key]

                else:
                    merge_dict(basedict[key], extradict[i][key])

    return basedict


class Figure(Complex):
    def __init__(self, *args, **kwargs):
        # not able to read dash_dict twice for some reason
        if "options" in kwargs:
            self.options = kwargs["options"].copy()
        if "input_json_path" in kwargs:
            self.input_json_path = kwargs["input_json_path"]

        if "dash_dict" in kwargs:
            self.dash_dict = kwargs["dash_dict"].copy()
            ## a switch
            self.config_dict = json_config.json2dict(self.dash_dict)
        else:
            self.config_dict = json_config.json2dict(kwargs["input_json_path"])

        ########################################################################
        # try replacing all self.get_chr_info() with self.chr_info

        # input data
        input_data = None
        if self.config_dict["Category"]["ideogram"]["patch"]["file"].get("dataframe", {}):
            input_data = self.config_dict["Category"]["ideogram"]["patch"]["file"].get(
                "dataframe", None
            )
        else:
            if isfile(self.config_dict["Category"]["ideogram"]["patch"]["file"]["path"]):
                input_data = self.config_dict["Category"]["ideogram"]["patch"]["file"]["path"]
            elif isfile(
                os.path.dirname(self.input_json_path)
                + "/"
                + self.config_dict["Category"]["ideogram"]["patch"]["file"]["path"]
            ):
                input_data = (
                    os.path.dirname(self.input_json_path)
                    + "/"
                    + self.config_dict["Category"]["ideogram"]["patch"]["file"]["path"]
                )

        if input_data is not None:
            self.chr_info = coord_config.chr_info(
                input_data,
                sep=self.config_dict["Category"]["ideogram"]["patch"]["file"]["sep"],
                custom_label=self.config_dict["Category"]["ideogram"]["patch"]["customoptions"][
                    "customlabel"
                ],
                custom_spacing=self.config_dict["Category"]["ideogram"]["patch"]["customoptions"][
                    "customspacing"
                ],
                custom_color=self.config_dict["Category"]["ideogram"]["patch"]["customoptions"][
                    "customcolor"
                ],
            )
        else:
            print('["ERROR"] Error in ideogram')

        ########################################################################

        self.ideogram_coord_config = self.get_ideogram_coord_config()

        self.layout_general = self.config_dict["General"]

        categories = self.config_dict["Category"]
        keyList = []
        for key in categories:
            if key in ["cytoband", "highlight", "annotation"]:
                if "file" not in categories[key]:
                    keyList.append(key)

        for key in keyList:
            categories.pop(key)

        # if cytoband is used, ideogram patch fill is automatically disabled!
        if "cytoband" in categories:
            categories["ideogram"]["patch"]["showfillcolor"] = False

        ###### for some reason annotation file path is not converted to io.StringIO format

        self.categories = categories

        self.ideogram = self.categories["ideogram"]

        self.ideogram_patch = self.ideogram["patch"]

        self.ideogram_majortick = self.ideogram["majortick"]
        self.ideogram_minortick = self.ideogram["minortick"]
        self.ideogram_ticklabel = self.ideogram["ticklabel"]

        self.ideogram_radius_dict = self.ideogram_patch["radius"]
        self.show_chr_annotation = self.ideogram_patch["chrannotation"]["show"]

        self.major_tick_radius_dict = self.ideogram_majortick["radius"]
        self.minor_tick_radius_dict = self.ideogram_minortick["radius"]
        self.tick_label_radius_dict = self.ideogram_ticklabel["radius"]
        self.SUM = self.get_ideogram_coord_config()["SUM"]

        self.degreerange = self.ideogram_patch["degreerange"]

        self.chr_color_dict = self.chr_info["chr_color_dict"]
        self.chr_label_dict = self.chr_info["chr_label_dict"]

        ############################################################
        # replace self.get_data_array_dict() & self.read_data() with ONE self.get_data attribute!!
        get_data_array_dict_ = {}
        get_read_data_ = {}

        for key in self.categories:
            if key == "ideogram":
                pass
            else:
                try:
                    get_data_array_dict_[key] = self.get_data_array_dict(key)
                except Exception:
                    pass

        self.get_data = get_data_array_dict_

    def generate_dash_dict(self):
        return self.config_dict

    def np_list_concat(self, x):
        if isinstance(x, list):
            try:
                return np.concatenate(x)
            except ValueError:
                return np.array(x)

        elif isinstance(x, np.ndarray):
            return x
        else:
            raise ValueError("input must be an ndarray or a list")

    def get_read_data(self, key):
        items = self.categories[key]
        sortindices = self.get_data_array_sortindex(key)

        def get_single_data(item, key, sortindex=None):
            assert isinstance(item, dict)

            if item["file"]["header"] in ["None", None, "none"]:
                unsorted_data = coord_config.read_data(
                    item["file"]["path"],
                    key,
                    # self.get_chr_info(),
                    self.chr_info,
                    sep=item["file"]["sep"],
                    header=None,
                )
            else:
                unsorted_data = coord_config.read_data(
                    item["file"]["path"],
                    key,
                    # self.get_chr_info(),
                    self.chr_info,
                    sep=item["file"]["sep"],
                )
            if item["sortbycolor"]:
                ## ONGOING
                assert sortindex is not None
                return unsorted_data[sortindex]
            else:
                return unsorted_data

        if isinstance(items, dict):
            return get_single_data(items, key, sortindex=sortindices)
        elif isinstance(items, list):
            return [
                *map(
                    lambda x, y: get_single_data(x, key, sortindex=y),
                    items,
                    sortindices,
                )
            ]

    def get_chr_info(self):
        # chr_info_file = self.config_dict['Category']['ideogram']['patch']['file']
        # custom_options = self.config_dict['Category']['ideogram']['patch']['customoptions']

        chr_info_dict = coord_config.chr_info(
            self.config_dict["Category"]["ideogram"]["patch"]["file"]["path"],
            sep=self.config_dict["Category"]["ideogram"]["patch"]["file"]["sep"],
            custom_label=self.config_dict["Category"]["ideogram"]["patch"]["customoptions"][
                "customlabel"
            ],
            custom_spacing=self.config_dict["Category"]["ideogram"]["patch"]["customoptions"][
                "customspacing"
            ],
            custom_color=self.config_dict["Category"]["ideogram"]["patch"]["customoptions"][
                "customcolor"
            ],
        )

        return chr_info_dict

    def get_ideogram_coord_config(self):
        ideogram_coord_config = coord_config.ideogram_coord_config(
            # self.get_chr_info(),
            self.chr_info,
            npoints=self.config_dict["Category"]["ideogram"]["patch"]["npoints"],
            show_major_tick=self.config_dict["Category"]["ideogram"]["majortick"]["show"],
            major_tick_spacing=self.config_dict["Category"]["ideogram"]["majortick"]["spacing"],
            show_minor_tick=self.config_dict["Category"]["ideogram"]["minortick"]["show"],
            minor_tick_spacing=self.config_dict["Category"]["ideogram"]["minortick"]["spacing"],
            show_tick_label=self.config_dict["Category"]["ideogram"]["ticklabel"]["show"],
            tick_label_spacing=self.config_dict["Category"]["ideogram"]["ticklabel"]["spacing"],
        )

        return ideogram_coord_config

    ## self.ideogram_complex() and self.ring_complex() is inherited

    def get_ideogram_theta(self):
        return self.ideogram_theta_list(
            self.ideogram_coord_config, self.SUM, degreerange=self.degreerange
        )

    def get_ideogram_complex(self):
        return self.ideogram_complex(
            self.ideogram_coord_config,
            self.SUM,
            self.degreerange,
            ideogram_radius_dict=self.ideogram_radius_dict,
        )

    # def get_ideogram_shapes(self):
    # return self.ideogram_path(self.get_ideogram_complex())

    def get_ring_complex(self):
        if isinstance(self.categories["ring"], dict):
            return self.ideogram_complex(
                self.ideogram_coord_config,
                self.SUM,
                degreerange=self.degreerange,
                ideogram_radius_dict=self.categories["ring"]["radius"],
            )
        elif isinstance(self.categories["ring"], list):
            return [
                *map(
                    lambda x: self.ideogram_complex(
                        self.ideogram_coord_config,
                        self.SUM,
                        degreerange=self.degreerange,
                        ideogram_radius_dict=x["radius"],
                    ),
                    self.categories["ring"],
                )
            ]

    def get_ring_paths_dict(self):
        def single_ring_dict(path, ring_dict):
            # ring data can only have one background color! therefore I join them
            # item = self.categories['ring']
            if isinstance(path, list):
                path = self.pathjoin(path)
            path_dict = dict(path=path)
            path_dict.update(ring_dict["layout"])
            return path_dict

        if isinstance(self.categories["ring"], dict):
            path = self.ideogram_path(self.get_ring_complex())
            path_dict = [single_ring_dict(path, self.categories["ring"])]

        elif isinstance(self.categories["ring"], list):
            path = [*map(lambda x: self.ideogram_path(x), self.get_ring_complex())]
            path_dict = [*map(lambda x, y: single_ring_dict(x, y), path, self.categories["ring"])]

        return path_dict

    def get_major_tick_path(self):
        major_tick_accum_coord_list = self.ideogram_coord_config["major_tick_accum_coord_list"]
        major_tick_theta = self.ideogram_tick_theta_list(
            self.ideogram_coord_config,
            major_tick_accum_coord_list,
            SUM=self.SUM,
            degreerange=self.degreerange,
        )
        major_tick_complex = self.tick_complex(
            major_tick_theta, tick_radius_dict=self.major_tick_radius_dict
        )

        return self.tick_path(major_tick_complex)

    def get_minor_tick_path(self):
        minor_tick_accum_coord_list = self.ideogram_coord_config["minor_tick_accum_coord_list"]
        minor_tick_theta = self.ideogram_tick_theta_list(
            self.ideogram_coord_config,
            minor_tick_accum_coord_list,
            SUM=self.SUM,
            degreerange=self.degreerange,
        )
        minor_tick_complex = self.tick_complex(
            minor_tick_theta, tick_radius_dict=self.minor_tick_radius_dict
        )

        return self.tick_path(minor_tick_complex)

    def get_ideogram_chrannot_theta(self):
        return self.ideogram_chrannot_theta(
            self.ideogram_coord_config, self.SUM, degreerange=self.degreerange
        )

    def get_ideogram_chrannot_complex(self):
        return self.ideogram_chrannot_complex(
            self.SUM,
            degreerange=self.degreerange,
            chr_annotation_radius_dict=self.ideogram_patch["chrannotation"]["radius"],
        )

    def get_tick_label_complex(self):
        tick_label_theta = self.ideogram_tick_label_theta_list(
            self.ideogram_coord_config, self.SUM, degreerange=self.degreerange
        )
        return self.tick_label_complex(
            tick_label_theta, tick_label_radius_dict=self.tick_label_radius_dict
        )

    def pathjoin(self, path_list):
        """this function will join any path_string_list into a single path_string, useful when fillcolor is the same"""
        return " ".join(path_list)

    def get_data_array_dict(self, key):
        assert key in self.categories
        items = self.categories[key]
        # print('key is {}'.format(key))
        # print(items)
        if key == "ideogram":
            raise ValueError("ideogram information should not be parsed in get_data_array()")

        else:

            def single_data_array(key, item):
                if "colorcolumn" not in item.keys():
                    item["colorcolumn"] = None
                if "sortbycolor" not in item.keys():
                    item["sortbycolor"] = False

                input_data = None

                # file_path=item['file']['path']
                if item["file"].get("dataframe", None) is not None:
                    if item["file"]["path"] == "":
                        input_data = item["file"].get("dataframe", {})
                else:
                    if isfile(item["file"]["path"]):
                        input_data = item["file"]["path"]
                    elif isfile(os.path.dirname(self.input_json_path) + "/" + item["file"]["path"]):
                        input_data = (
                            os.path.dirname(self.input_json_path) + "/" + item["file"]["path"]
                        )

                # item['file']['path']=file_path

                if item["file"]["header"] in ["None", None, "none"]:
                    return coord_config.data_array(
                        input_data,
                        key,
                        # self.get_chr_info(),
                        self.chr_info,
                        sep=item["file"]["sep"],
                        header=None,
                        colorcolumn=item["colorcolumn"],
                        sortbycolor=item["sortbycolor"],
                    )
                else:
                    return coord_config.data_array(
                        input_data,
                        key,
                        # self.get_chr_info(),
                        self.chr_info,
                        sep=item["file"]["sep"],
                        colorcolumn=item["colorcolumn"],
                        sortbycolor=item["sortbycolor"],
                    )

            if isinstance(items, dict):
                return single_data_array(key, items)
            elif isinstance(items, list):
                try:
                    return [*map(lambda x: single_data_array(key, x), items)]
                except Exception:
                    print(f'[WARN] warning on category "{key}"')
                    # print(key)

    def get_data_array(self, key):
        # if self.get_data_array_dict(key) is None:
        if self.get_data[key] is None:
            pass
        else:
            # if not isinstance(self.get_data_array_dict(key), list):
            if not isinstance(self.get_data[key], list):
                # return self.get_data_array_dict(key)['data_array']
                return self.get_data[key]["data_array"]
            else:
                # return [*map(lambda x: x['data_array'], self.get_data_array_dict(key))]
                return [*map(lambda x: x["data_array"], self.get_data[key])]

    def get_data_array_sortindex(self, key):
        # only used when sortbycolor is True! in other cases its just range(len(data_array))

        # if not isinstance(self.get_data_array_dict(key), list):
        if not isinstance(self.get_data[key], list):
            # return self.get_data_array_dict(key)['sortindex']
            return self.get_data[key]["sortindex"]
        else:
            # return [*map(lambda x: x['sortindex'], self.get_data_array_dict(key))]
            return [*map(lambda x: x["sortindex"], self.get_data[key])]

    def get_data_complexes(self, key, return_path=True):
        assert key in self.categories
        data_array = self.get_data_array(key)
        items = self.categories[key]

        # print('key is {}'.format(key))
        # print(key)
        # print(data_array)

        def single_data_complex(data_array, key, item):
            if key != "highlight":
                if key == "cytoband":
                    item["radius"] = self.ideogram_radius_dict
                elif key == "annotation" and item["customradius"]:
                    try:
                        assert isinstance(item["radiuscolumn"], int)
                    except AssertionError:
                        print(f"[WARN] Please enter a valid radiuscolumn under annotation")
                    try:
                        radiuscolumn = item["radiuscolumn"]
                    except IndexError:
                        print(
                            f"[WARN] the column you entered is out of bound, notice that column starts with 0, not 1 "
                        )
                    try:
                        data_array[:, radiuscolumn].astype("float")

                    ## DEBUG, the below should be ValueError, changing to Exception temporarily
                    except Exception:
                        print(radiuscolumn)
                        print(data_array[:, radiuscolumn])
                        print(f"[WARN] Please make sure to enter numeric value for radius column")

                    item["radius"] = {"R": data_array[:, radiuscolumn]}

                radius_dict = item["radius"]

            else:
                radius_dict = {
                    "R0": data_array[:, item["R0column"]],
                    "R1": data_array[:, item["R1column"]],
                }

            if key == "annotation":
                data_complex = self.data_complex(
                    self.ideogram_coord_config,
                    data_array,
                    key,
                    radius_dict,
                    self.SUM,
                    degreerange=self.degreerange,
                    custom_offset_degree=item["customoffsetdegree"],
                )
            else:
                data_complex = self.data_complex(
                    self.ideogram_coord_config,
                    data_array,
                    key,
                    radius_dict,
                    self.SUM,
                    degreerange=self.degreerange,
                    return_path=return_path,
                )
            return data_complex

        if isinstance(items, dict):
            return single_data_complex(data_array, key, items)
        elif isinstance(items, list):
            # print("key: "+key)
            # print('single_data_complex data array')
            # print(data_array)
            # print('single_data_complex items')
            # print(items)
            try:
                return [*map(lambda x, y: single_data_complex(x, key, y), data_array, items)]
            except Exception:
                print(f"[WARN] debugging data_array")
                print(data_array)  # got empty
                print(f"[WARN] debugging items")
                print(items)

    def get_hovertext(self, key):
        assert key in [
            "histogram",
            "line",
            "area",
            "scatter",
            "tile",
            "heatmap",
            "link",
            "ribbon",
            "twistedribbon",
        ]

        if not isinstance(self.get_data[key], list):
            hovertextformat = self.categories[key]["hovertextformat"]

            # a = self.get_read_data(key)
            a = self.get_data[key]["read_data"]
            hvtext = []
            if key not in ["link", "ribbon", "twistedribbon"]:
                assert a.shape[1] >= 3

                if key in ["histogram", "tile", "heatmap"]:
                    list_count = [*map(lambda x: len(x), self.get_data_complexes(key))]

                    for i in range(len(a)):
                        k = 0
                        while k < list_count[i]:
                            k += 1
                            hvtext.append(eval(hovertextformat))
                else:
                    for i in range(len(a)):
                        hvtext.append(eval(hovertextformat))

            else:
                # for link, ribbon and twistedribbon hovertext is a list of two elememnt
                assert a.shape[1] >= 6
                for i in range(len(a)):
                    hvtext.append(eval(hovertextformat[0]))
                    hvtext.append(eval(hovertextformat[0]))
                    hvtext.append(eval(hovertextformat[1]))
                    hvtext.append(eval(hovertextformat[1]))
        else:
            hvtext = []
            try:
                data_array_list = [*map(lambda t: t["read_data"], self.get_data[key])]

            except Exception:
                print(f"[WARN] debugging...")
                print(self.get_data[key])
                print(f"[WARN] end of debugging...")

            for j in range(len(data_array_list)):
                hovertextformat = self.categories[key][j]["hovertextformat"]
                hvtext.append([])
                a = data_array_list[j]

                if key not in ["link", "ribbon", "twistedribbon"]:
                    assert a.shape[1] >= 3

                    if key in ["histogram", "tile", "heatmap"]:
                        # print("DEV")
                        # print(f"j={j}")
                        # print(f"data complexes {self.get_data_complexes(key)}")
                        list_count = [*map(lambda x: len(x), self.get_data_complexes(key)[j])]

                        for i in range(len(a)):
                            k = 0
                            while k < list_count[i]:
                                k += 1
                                hvtext[j].append(eval(hovertextformat))

                    else:
                        for i in range(len(a)):
                            hvtext[j].append(eval(hovertextformat))

                else:
                    assert a.shape[1] >= 6
                    for i in range(len(a)):
                        hvtext[j].append(eval(hovertextformat[0]))
                        hvtext[j].append(eval(hovertextformat[0]))
                        hvtext[j].append(eval(hovertextformat[1]))
                        hvtext[j].append(eval(hovertextformat[1]))
        return hvtext

    def get_traces(self, key):
        ## always return a list, this makes concatenation easier
        # for line plot, there shouldn't be sortbycolor

        items = self.categories[key]
        complexes = self.get_data_complexes(key, return_path=False)
        data_arrays = self.get_data_array(key)
        hovertexts = self.get_hovertext(key)

        def single_trace(key, Complex, item, data_array, hovertext):
            assert key not in [
                "cytoband",
                "ideogram",
                "ring",
                "annotation",
                "highlight",
                "connector",
            ]
            # For scatter plot the complex is an ndarray, for lines it would be a list of ndarray separated by chromosomes
            # for other nonvisible plots, they can be concatenated into one ndarray
            assert isinstance(item, dict)

            # DEBUGGING
            # print(f"key={key}")
            # if key == "tile":
            #     print(item["trace"])

            if key != "line":
                if isinstance(Complex, list):
                    Complex = np.concatenate(Complex)

            if key == "line":
                # the only time when Complex is a list of ndarray
                trace = []
                index = np.cumsum([0] + [*map(lambda x: len(x), Complex)])

                def divide(l, index):
                    # this function creates a generator object for hovertext so I know how many hovertext element to take for each chromosome
                    for n in range(len(index) - 1):
                        yield l[index[n] : index[n + 1]]

                hovertext_generator = divide(hovertext, index)

                if item["colorcolumn"] in [None, "None"]:
                    for i in range(len(Complex)):
                        trace.append(
                            go.Scatter(
                                x=Complex[i].real,
                                y=Complex[i].imag,
                                text=next(hovertext_generator),
                                name=item["name"],
                            )
                        )
                        trace[-1].update(item["trace"])

                elif item["colorcolumn"] == "ideogram":
                    chr_label = data_array[:, 0]

                    # color = self.get_chr_info()['chr_fillcolor']
                    color = self.chr_info["chr_fillcolor"]
                    tmp_trace = item["trace"].copy()
                    for i in range(len(Complex)):
                        trace.append(
                            go.Scatter(
                                x=Complex[i].real,
                                y=Complex[i].imag,
                                text=next(hovertext_generator),
                            )
                        )

                        tmp_trace["line"]["color"] = color[i]
                        tmp_trace["marker"]["color"] = color[i]

                        trace[-1].update(tmp_trace)

            else:
                if Complex.ndim != 1:
                    Complex = Complex.ravel()

                if key in ["link", "ribbon", "twistedribbon"]:
                    index = []
                    for i in range(len(Complex)):
                        if i % 6 in [0, 1, 4, 5]:
                            index.append(i)
                    Complex = Complex[index]
                # if not "name" in item:
                #    item["name"] = "nolegend"
                #    item["showlegend"] = "False"
                trace = go.Scatter(x=Complex.real, y=Complex.imag, text=hovertext)
                # print(item['trace'])
                trace.update(item["trace"])

                if key == "scatter":
                    # if 'color' not in item['trace']['marker']:

                    chr_label = data_array[:, 0]

                    if item["colorcolumn"] == "ideogram":
                        # follows ideogram color
                        color = [*map(lambda x: self.chr_color_dict[x], chr_label)]

                    elif isinstance(item["colorcolumn"], int):
                        n = item["colorcolumn"]
                        assert isinstance(n, int)
                        color = colors.to_rgb(data_array[:, n])
                    else:
                        ### ONGOING
                        # assert item['colorcolumn'] in [None, 'None']
                        color = item["trace"]["marker"]["color"]

                    trace["marker"].update(color=color)
            # try:
            #    if isinstance(trace["marker"]["color"], str):
            #        if trace["marker"]["color"] != "gray":
            #            trace["name"] = trace["uid"]
            #    else:
            #        if not "gray" in trace["marker"]["color"]:
            #            trace["name"] = trace["uid"]
            # except ValueError:
            #    print(item)
            #    print(trace)
            #    exit()
            return trace

        if isinstance(items, dict):
            trace = single_trace(key, complexes, items, data_arrays, hovertexts)
            if not isinstance(trace, list):
                trace = [trace]
            return trace

        elif isinstance(items, list):
            return [
                *map(
                    lambda w, x, y, z: single_trace(key, w, x, y, z),
                    complexes,
                    items,
                    data_arrays,
                    hovertexts,
                )
            ]

    def trace(self):
        # aggregate all get_traces element into one trace, list variable
        trace = []
        for key in self.categories.keys():
            if key not in [
                "cytoband",
                "ideogram",
                "ring",
                "annotation",
                "highlight",
                "connector",
            ]:
                trace += self.get_traces(key)
        if hasattr(self, "options") and self.options:
            number_trace = []
            for tr in trace:
                if tr["uid"] in ["cytoband_tile", "transloc", None] or tr["uid"].startswith(
                    "extra_"
                ):
                    continue
                if isinstance(tr["marker"]["color"], list) or isinstance(
                    tr["marker"]["color"], tuple
                ):
                    number_trace.extend(tr["marker"]["color"])
                elif isinstance(tr["marker"]["color"], np.ndarray):
                    number_trace.extend(tr["marker"]["color"].tolist())
                else:
                    number_trace.append(tr["marker"]["color"])
            number_trace = color_to_name(list(set(number_trace)))
            # print(number_trace)
            # exit()
            trace_dict = {key: {} for key in number_trace}
            graph_obj = []
            # Dont know why but some dot are in double needed to fix that look above
            # For all CNV level dot
            for tr in trace:
                if tr["uid"] is not None:
                    if tr["uid"] == "genes" and tr["marker"]["size"] != 5:
                        continue
                    if tr["uid"].startswith("cnv_scatter") and tr["marker"]["opacity"] != 1:
                        continue
                    if tr["uid"] == "transloc":
                        tr["name"] = "BND"
                        graph_obj.append(tr)
                        continue
                    if tr["uid"] == "cytoband_tile":
                        tr["name"] = "CYTOBAND"
                        graph_obj.append(tr)
                        continue
                    if tr["uid"].startswith("extra_"):
                        graph_obj.append(tr)
                        continue
                    # if tr["uid"] not in ["genes, cytoband_tiles, transloc"]:
                    # Becarefull a value alone for SNV indels
                    if isinstance(tr["marker"]["color"], str):
                        tr["marker"]["color"] = [tr["marker"]["color"]]
                    for j, color in enumerate(tr["marker"]["color"]):
                        # Issues with cytoband only chr1 and 11 displayed
                        # trace_dict.update(process_legend(color, tr, trace_dict))
                        color = color_to_name(color)
                        # For each item in dot
                        for key, val in vars(tr)["_orphan_props"].items():
                            if isinstance(val, str):
                                if key not in trace_dict[color]:
                                    trace_dict[color][key] = val
                            elif isinstance(val, np.ndarray):
                                if key not in trace_dict[color]:
                                    trace_dict[color][key] = []
                                trace_dict[color][key].append(val[j])
                            elif isinstance(val, list):
                                if key not in trace_dict[color]:
                                    trace_dict[color][key] = []
                                trace_dict[color][key].append(val[j])
                            elif isinstance(val, dict):
                                if key not in trace_dict[color]:
                                    trace_dict[color][key] = {}
                                for k, v in val.items():
                                    if (
                                        isinstance(v, str)
                                        or isinstance(v, int)
                                        or isinstance(v, float)
                                    ):
                                        if k not in trace_dict[color][key]:
                                            trace_dict[color][key][k] = v
                                    elif isinstance(v, np.ndarray):
                                        if k not in trace_dict[color][key]:
                                            trace_dict[color][key][k] = []
                                        trace_dict[color][key][k].append(v[j])
                                    elif isinstance(v, list):
                                        if k not in trace_dict[color][key]:
                                            trace_dict[color][key][k] = []
                                        trace_dict[color][key][k].append(v[j])
            cast_color = get_swap_dict(self.options["Color"])
            trace_ = []
            for clrs, values in trace_dict.items():
                values["name"] = cast_color[clrs]
                trace_.append(go.Scatter(values))
            # exit()
            for val in graph_obj:
                if not "name" in val:
                    val["name"] = val["uid"]
                if "name" not in ["BND", "CYTOBAND"]:
                    val["showlegend"] = False
                trace_.append(val)
            return trace_
        else:
            return trace

    def get_paths_dict(self, key):
        # will join path_list into a path string if sortbycolor
        ### deal with:
        # histogram, ribbon, twistedribbon, if fillcolor is true, then linecolor==fillcolor!
        # cytoband, heatmap (sortbycolor=True)
        # tile, link (no fillcolor), color indicates to linecolor

        assert key not in ["scatter", "annotation", "line"]
        items = self.categories[key]
        data_arrays = self.get_data_array(key)
        data_complexes = self.get_data_complexes(key)

        def single_path(key, data_array, data_complex, item):
            if key in ["ribbon", "twistedribbon"]:
                interval_theta_array_0 = maths.to_theta(
                    data_array[:, 1:3], self.SUM, degreerange=self.degreerange
                )
                interval_theta_array_1 = maths.to_theta(
                    data_array[:, 4:6], self.SUM, degreerange=self.degreerange
                )
            else:
                interval_theta_array_0, interval_theta_array_1 = None, None

            if key == "highlight":
                path_list = self.data_path(
                    self.ideogram_coord_config,
                    key,
                    data_complex,
                    self.SUM,
                    degreerange=self.degreerange,
                    radius_dict={"R0": item["R0column"], "R1": item["R1column"]},
                    interval_theta_array_0=interval_theta_array_0,
                    interval_theta_array_1=interval_theta_array_1,
                )
            else:
                if key == "cytoband":
                    item["radius"] = self.ideogram_radius_dict
                path_list = self.data_path(
                    self.ideogram_coord_config,
                    key,
                    data_complex,
                    self.SUM,
                    degreerange=self.degreerange,
                    radius_dict=item["radius"],
                    interval_theta_array_0=interval_theta_array_0,
                    interval_theta_array_1=interval_theta_array_1,
                )

            if item["colorcolumn"] in [None, "None"]:
                path = " ".join(path_list)
                paths_dict = dict(path=path)
                paths_dict.update(item["layout"])

            # elif isinstance(item['colorcolumn'], int):
            else:
                # custom line color only, and no fill color: link, tile
                # custom fill color only: all others

                paths_dict = []
                if item["colorcolumn"] == "ideogram":
                    if not key == "area":
                        color = [*map(lambda x: self.chr_color_dict[x], data_array[:, 0])]

                    else:
                        # color = self.get_chr_info()['chr_fillcolor']
                        color = self.chr_info["chr_fillcolor"]

                elif isinstance(item["colorcolumn"], int):
                    n = item["colorcolumn"]
                    color = data_array[:, n]

                if key == "heatmap":
                    color = maths.val2heatmap(color, palatte_dict=item["palatte"])

                if key == "highlight":
                    o = item["opacitycolumn"]
                    opacity = data_array[:, o]

                    for i in range(len(color)):
                        paths_dict.append(
                            dict(
                                path=path_list[i],
                                fillcolor=color[i],
                                opacity=opacity[i],
                            )
                        )

                        paths_dict[i].update(item["layout"])

                elif key in [
                    "heatmap",
                    "cytoband",
                    "ribbon",
                    "twistedribbon",
                    "tile",
                    "link",
                    "area",
                    "histogram",
                ]:
                    for i in range(len(path_list)):
                        paths_dict.append(dict(path=path_list[i]))

                        paths_dict[i].update(copy.deepcopy(item["layout"]))

                        if key in ["tile", "link"]:
                            paths_dict[i]["line"]["color"] = color[i]
                        else:
                            paths_dict[i]["line"]["color"] = color[i]
                            paths_dict[i]["fillcolor"] = color[i]

            return paths_dict

        if isinstance(items, dict):
            if isinstance(single_path(key, data_arrays, data_complexes, items), dict):
                return [single_path(key, data_arrays, data_complexes, items)]
            else:
                return single_path(key, data_arrays, data_complexes, items)

        elif isinstance(items, list):
            ## ONGOING
            rtrn = [
                *map(
                    lambda x, y, z: single_path(key, x, y, z),
                    data_arrays,
                    data_complexes,
                    items,
                )
            ]
            return [*map(lambda a: a if isinstance(a, list) else [a], rtrn)]

            # return [*map(lambda x, y, z: single_path(key, x, y, z), data_arrays, data_complexes, items)]

        else:
            raise KeyError("{} file does not exist!".format(key))

    def get_annotations_dict(self):
        """layout['annotations'] can only append one text annotation at a time!! be careful"""

        # please make sure you only have one annotation file!

        assert "annotation" in self.categories.keys()
        try:
            assert isinstance(self.categories["annotation"], dict)
        except AssertionError:
            print("Please use one and only one annotation file")

        text_complex = self.get_data_complexes("annotation")
        text_array = self.get_data_array("annotation")[:, 2]
        text_theta = maths.to_theta(
            self.get_data_array("annotation")[:, 1],
            self.SUM,
            degreerange=self.degreerange,
        )

        textangle = self.angleconvert(
            text_theta,
            angleoffset=self.categories["annotation"]["textangle"]["angleoffset"],
            anglelimit=self.categories["annotation"]["textangle"]["anglelimit"],
        )

        if self.categories["annotation"]["customoffsetdegree"]:
            textangle += self.get_data_array("annotation")[:, 3]

        if self.categories["annotation"]["fonttype"] == "bold":
            if len(text_array) == 1:
                text = "<b>{}</b>".format(text_array)
            else:
                text = [*map(lambda x: "<b>{}</b>".format(x), text_array)]
        elif self.categories["annotation"]["fonttype"] == "italic":
            if len(text_array) == 1:
                text = "<i>{}</i>".format(text_array)
            else:
                text = [*map(lambda x: "<i>{}</i>".format(x), text_array)]
        elif self.categories["annotation"]["fonttype"] in [
            "bold+italic",
            "italic+bold",
        ]:
            if len(text_array) == 1:
                text = "<b><i>{}</i></b>".format(text_array)
            else:
                text = [*map(lambda x: "<b><i>{}</i></b>".format(x), text_array)]
        else:
            if len(text_array) == 1:
                text = "{}".format(text_array)
            else:
                text = [*map(lambda x: "{}".format(x), text_array)]

        if len(text_array) == 1:
            annotation_dict = dict(
                x=text_complex.real, y=text_complex.imag, text=text, textangle=textangle
            )
            annotation_dict.update(copy.deepcopy(self.categories["annotation"]["layout"]))
            return [annotation_dict]
        else:
            if not self.categories["annotation"]["customcolor"]:
                annotation_dict_list = [
                    *map(
                        lambda a, b, c: merge_dict(
                            dict(x=a.real, y=a.imag, text=b, textangle=c),
                            copy.deepcopy(self.categories["annotation"]["layout"]),
                        ),
                        text_complex,
                        text,
                        textangle,
                    )
                ]
            else:
                assert isinstance(self.categories["annotation"]["colorcolumn"], int)
                n = self.categories["annotation"]["colorcolumn"]
                assert n <= self.get_data_array("annotation").shape[1] - 1
                colors = self.get_data_array("annotation")[:, n].tolist()
                annotation_dict_list = [
                    *map(
                        lambda a, b, c, d: merge_dict(
                            dict(x=a.real, y=a.imag, text=b, textangle=c),
                            copy.deepcopy(self.categories["annotation"]["layout"]),
                            dict(font=dict(color=d)),
                        ),
                        text_complex,
                        text,
                        textangle,
                        colors,
                    )
                ]
            return annotation_dict_list

    def angleconvert(self, theta, angleoffset=-90, anglelimit=360, custom_offset_degree=0):
        # theta could be an ndarray or a list of ndarray
        # returns a processed degree data based on anglelimit
        # anglelimit is the degree limit before applying angleoffset and custom_offset!
        # custom_offset_degree is a 1D ndarray of the same shape as theta

        constant = 180 / np.pi
        assert isinstance(angleoffset, (int, float))
        assert isinstance(anglelimit, (int, float))

        if not isinstance(theta, list):
            degree = theta * constant

            for i in range(len(degree)):
                while degree[i] >= anglelimit:
                    degree[i] -= anglelimit

            degree += angleoffset + custom_offset_degree
            return degree
        else:
            assert isinstance(custom_offset_degree, list)
            assert len(theta) == len(custom_offset_degree)

            degree = [*map(lambda x: x * constant), theta]
            for i in range(len(degree)):
                for j in range(len(degree[i])):
                    while degree[i][j] >= anglelimit:
                        degree[i][j] -= anglelimit
                degree[i] += angleoffset + custom_offset_degree[i]
            return degree

    def layout(self):
        ## new version of plotly says layout['shapes'], ['annotations'] is a tuple, can't be used with list append
        ## I'm trying to add a temporary layout_shapes list variable collecting everything before assigning to layout['shapes']

        layout = go.Layout(self.layout_general)

        layout_shapes = []
        layout_annotations = []

        if self.ideogram_patch["show"]:
            if not self.ideogram_patch["showfillcolor"]:
                # if no show fillcolor, there is no need to separate pathstring into list, hense the join
                # seems like for dictionary update, the dict variable needs to be defined first

                ideogram_pathstring = self.pathjoin(self.ideogram_path(self.get_ideogram_complex()))

                layout_dict = self.ideogram_patch["layout"]
                layout_dict.update(dict(path=ideogram_pathstring))
                layout_shapes.append(layout_dict)
            else:
                ideogram_path_list = self.ideogram_path(self.get_ideogram_complex())

                for i in range(len(ideogram_path_list)):
                    layout_shapes.append(
                        dict(
                            path=ideogram_path_list[i],
                            # fillcolor=self.get_chr_info()['chr_fillcolor'][i]
                            fillcolor=self.chr_info["chr_fillcolor"][i],
                        )
                    )
                    layout_shapes[i].update(self.ideogram_patch["layout"])

        if self.ideogram_patch["chrannotation"]["show"]:
            chrannot_theta = self.get_ideogram_chrannot_theta()

            chrannot_theta = self.np_list_concat(chrannot_theta)

            chrannot_complex = maths.to_complex(
                chrannot_theta, self.ideogram_patch["chrannotation"]["radius"]["R"]
            )

            chrannot_complex = self.np_list_concat(chrannot_complex)

            chrannot_angleoffset = self.ideogram_patch["chrannotation"]["textangle"]["angleoffset"]
            chrannot_anglelimit = self.ideogram_patch["chrannotation"]["textangle"]["anglelimit"]

            if self.ideogram_patch["chrannotation"]["fonttype"] == "bold":
                chrannot_text = [
                    *map(
                        lambda x: "<b>{}</b>".format(x),
                        # self.get_chr_info()['chr_label']
                        self.chr_info["chr_label"],
                    )
                ]

            elif self.ideogram_patch["chrannotation"]["fonttype"] == "italic":
                chrannot_text = [
                    *map(
                        lambda x: "<i>{}</i>".format(x),
                        # self.get_chr_info()['chr_label']
                        self.chr_info["chr_label"],
                    )
                ]

            elif self.ideogram_patch["chrannotation"]["fonttype"] in [
                "bold+italic",
                "italic+bold",
            ]:
                chrannot_text = [
                    *map(
                        lambda x: "<b><i>{}</i></b>".format(x),
                        # self.get_chr_info()['chr_label']
                        self.chr_info["chr_label"],
                    )
                ]

            else:
                chrannot_text = [
                    *map(
                        lambda x: "{}".format(x),
                        # self.get_chr_info()['chr_label']
                        self.chr_info["chr_label"],
                    )
                ]

            textangle = self.angleconvert(
                chrannot_theta,
                angleoffset=chrannot_angleoffset,
                anglelimit=chrannot_anglelimit,
            )

            for i in range(len(chrannot_complex)):
                layout_annotations.append(
                    dict(
                        x=chrannot_complex[i].real,
                        y=chrannot_complex[i].imag,
                        text=chrannot_text[i],
                        textangle=textangle[i],
                    )
                )
                layout_annotations[i].update(self.ideogram_patch["chrannotation"]["layout"])

        if self.ideogram_majortick["show"]:
            layout_shapes.append(dict(path=self.get_major_tick_path()))
            layout_shapes[-1].update(self.ideogram_majortick["layout"])

        if self.ideogram_minortick["show"]:
            layout_shapes.append(dict(path=self.get_minor_tick_path()))
            layout_shapes[-1].update(self.ideogram_minortick["layout"])

        if self.ideogram_ticklabel["show"]:
            ticklabel_text = self.get_ideogram_coord_config()["tick_label_non_accum_list"]
            if isinstance(ticklabel_text, list):
                ticklabel_text = np.concatenate(ticklabel_text)

            ticklabel_coord = self.np_list_concat(
                self.get_ideogram_coord_config()["tick_label_accum_coord_list"]
            )
            ticklabel_theta = maths.to_theta(
                ticklabel_coord, self.SUM, degreerange=self.degreerange
            )
            ticklabel_angle = self.angleconvert(
                ticklabel_theta,
                angleoffset=self.ideogram_ticklabel["textangle"]["angleoffset"],
                anglelimit=self.ideogram_ticklabel["textangle"]["anglelimit"],
            )

            ticklabel_complex = np.concatenate(self.get_tick_label_complex())

            if self.ideogram_ticklabel["textformat"] == "Kb":
                ticklabel_text = [*map(lambda x: "{} Kb".format(x // 1000), ticklabel_text)]

            elif self.ideogram_ticklabel["textformat"] == "Mb":
                ticklabel_text = [*map(lambda x: "{} Mb".format(x // 1000000), ticklabel_text)]

            elif self.ideogram_ticklabel["textformat"] == "Gb":
                ticklabel_text = [*map(lambda x: "{} Gb".format(x // 1000000000), ticklabel_text)]

            else:
                raise ValueError("acceptable ideogram ticklabel textformats are: Kb & Mb & Gb")

            for i in range(len(ticklabel_complex)):
                layout_annotations.append(
                    dict(
                        x=ticklabel_complex[i].real,
                        y=ticklabel_complex[i].imag,
                        text=ticklabel_text[i],
                        textangle=ticklabel_angle[i],
                    )
                )
                layout_annotations[-1].update(self.ideogram_ticklabel["layout"])

            # due to the uncertain number of plot each type, we'll use extend [] instead of append!
        if "ring" in self.categories.keys():
            # the idea is to always to draw the ring background first, whereas highlight and custom annotation will be drawn last!
            # ONGOING, maybe to deprecate show=True for ring?

            layout_shapes.extend(self.get_ring_paths_dict())

        for key in self.categories.keys():
            if key in [
                "histogram",
                "cytoband",
                "area",
                "tile",
                "heatmap",
                "link",
                "ribbon",
                "twistedribbon",
                "connector",
            ]:
                # ONGOING, maybe to deprecate self.categories[key]['show']
                # if self.categories[key]['show']:

                if isinstance(self.get_paths_dict(key)[0], dict):
                    try:
                        layout_shapes.extend(self.get_paths_dict(key))
                    except Exception:
                        print("[WARN] Error trying to plot single {}".format(key))
                        break
                elif isinstance(self.get_paths_dict(key)[0], list):
                    try:
                        layout_shapes.extend(sum(self.get_paths_dict(key), []))
                    except Exception:
                        print("[WARN] Error trying to plot multiple {}".format(key))
                        print(self.get_paths_dict(key))
                        break

        if "highlight" in self.categories.keys():
            if self.categories["highlight"]["show"]:
                layout_shapes.extend(self.get_paths_dict("highlight"))

        if "annotation" in self.categories.keys():
            # print(self.categories['annotation'])
            # a = pd.read_csv(self.categories['annotation']['file']['path'])

            try:
                if self.categories["annotation"]["show"]:
                    layout_annotations.extend(self.get_annotations_dict())
            except Exception:
                layout_annotations.extend(self.get_annotations_dict())

        layout["shapes"] = layout_shapes
        layout["annotations"] = layout_annotations

        return layout

    @timeit
    def fig(self):
        try:
            trace_obj = self.trace()
            for datas in trace_obj:
                if hasattr(datas, "marker") and hasattr(datas, "text"):
                    if datas["marker"]["color"] is not None:
                        if len(datas["marker"]["color"]) != len(datas["text"]):
                            datas["marker"]["color"] = np.repeat(
                                self.options["Color"][datas["name"]], len(datas["text"])
                            )
            #with open("trace_layout_new.json", "w+") as trace:
            #    js = plotly.io.to_json(go.Figure(data=trace_obj, layout=self.layout()), pretty=True)
            #    trace.write(js)

            return go.Figure(data=trace_obj, layout=self.layout())
        except IndexError:
            # to deal with an issue like this when using dash (specifically line plot):
            ## self.trace() => [Scatter(x0), [Scatter(x1), Scatter(x2), Scatter(x3),]]
            return go.Figure(
                data=sum(
                    [
                        *map(
                            lambda t: [t] if not isinstance(t, list) else t,
                            self.trace(),
                        )
                    ],
                    [],
                ),
                layout=self.layout(),
            )
