from vcf2circos.plotcategories.plotconfig import Plotconfig


class Cytoband(Plotconfig):
    """
    "scatter":{
        "pattern":{
        ...
    },  "data":{
        ...
    }}
    """

    def __init__(
        self,
        filename,
        options,
        show,
        file,
        radius,
        sortbycolor,
        colorcolumn,
        hovertextformat,
        trace_car,
        data,
        layout,
    ):
        super().__init__(
            filename,
            options,
            show,
            file,
            radius,
            sortbycolor,
            colorcolumn,
            hovertextformat,
            trace_car,
            data,
            layout,
        )

        # Cytoband params
        # Need creation of dict in options attribute, regarding data input
        self.file = file
        self.colorcolumn = (3,)
        self.hovertextformat = (' "<b>{}</b>".format(a[i,0])',)
        self.trace = {
            "uid": "cytoband",
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 1,
                "symbol": 0,  # 8
                "color": "black",
                "opacity": 1,
            },
        }
        self.layout = {
            "type": "path",
            "layer": "below",
            "opacity": 1.0,
            "line": {"color": "black", "width": 0},
        }

    def merge_options(self):
        dico = {}
        dico["show"] = self.show
        dico["file"] = self.file
        dico["sortbycolor"] = self.sortbycolor
        dico["colorcolumn"] = self.colorcolumn
        dico["hovertextformat"] = self.hovertextformat
        dico["trace"] = self.trace
        dico["layout"] = self.layout
        return dico
