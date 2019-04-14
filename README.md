# PCircos
PCircos is a python package based on plotly which helps generating Circos plot.

## Dependencies

The following are needed for PCircos to work:


python >= 3.6.8
[plotly][1]
[colorlover][2]

To install them 


## Installation
```git clone https://github.com/CJinnny/PCircos.git```

## Usage
This package allows you to generate Circos plot in html format with customizable hover text. All you need to do is to write your custom json configuration file.
To run an example, try:

```python PCircos.py demo_data/demo_params.json```


## After creating your config.json file, type the following in the command line:

```python PCircos.py /path_to_your_config.json/config.json```

Note that path_to_your_custom.json should be changed to where your config.json file is.
An html page would pop up on your browser



[1]: [https://plot.ly/python/]
[2]: [https://github.com/jackparmer/colorlover]
