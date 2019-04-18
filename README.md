# PCircos
PCircos is a python package based on plotly which helps generating Circos plot.

## Dependencies

The following are needed for PCircos to work:

- [plotly](https://plot.ly/python/)
- [colorlover](https://github.com/jackparmer/colorlover)
- [dash](https://dash.plot.ly/)
- [dash colorscales](https://github.com/plotly/dash-colorscales)

To install dependencies
```
$ pip install plotly 
$ pip install colorlover
$ pip install dash==0.39.0
$ pip install dash-daq==0.1.0  
$ pip install dash_colorscales
$ pip install -i https://test.pypi.org/simple/ colorpicker-box
```
or if you don't have root privilege
```
$ sudo pip install plotly
$ sudo pip install colorlover
```

## Installation
```git clone https://github.com/CJinnny/PCircos.git```

## Usage
This package allows you to generate Circos plot in html format with customizable hover text. All you need to do is to write your custom json configuration file.
To run an example, try:

```python PCircos.py demo_data/demo_params.json```

## How to write your own json configuration file

- ### General
  informtion about figure width, height, title etc. You dont' have to include them in your json file and default information will be filled.
  
- ### Category

  information about ideogram (chromosome list and size) and cytoband (karyotype), ideogram annotation, ticks.
  
  - #### ideogram
    contains information about chromosome size, majortick, minortick, ticklabel. The only information you need to write is the file dictionary of chromosome size data. If you use a txt file, you should do `"sep": "\t"`, if you use csv file, you should do `"sep": ","`, it is advised that you keep the first row as header. 
    If you want to customize chromosome labels, do ```"customlabel": "True"``` and be sure to include custom labels are in the 3rd column in your chromosome size file. If you want to customize chromosome colors, do ```"customcolor": "True"``` and be sure to include custom color are in the 4th column in your chromosome size file. If you want to make ideogram patches thinner or thicker, change the radius parameters. If you want to make ideogram patches more transparent or opaque, change the opacity parameters in layout.
    
  - #### cytoband
    contains information about karyotype, please be sure to turn off ideogram fillcolor by doing ```"showfillcolor": "False"``` if you want to show the karyotype color instead
  
 - #### histogram
  contains information about histogram, you can plot one, or multiple histogram (use square bracket). If you want histogram colors to follow that of the ideogram, do ```"colorcolumn": "ideogram"```, or if you want single color , do ```"colorcolumn": "None"``` and specify the fillcolor in layout, or if you want custom color, do ```"colorcolumn": 5``` or whichever colum contains your custom color information. If you want to make your histogram pointing outward, make sure to set "R0" < "R1" and vice versa. 


## After creating your config.json file, type the following in the command line

```python PCircos.py /path_to_your_config.json/config.json```

Note that path_to_your_custom.json should be changed to where your config.json file is.
An html page would pop up on your browser

