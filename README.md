# Introduction

Package vcf2circos is a python package based on Plotly which helps generating Circos plot, from a VCF file or a JSON configuration file.

See documentation and code in [GitHub vcf2circos](https://github.com/bioinfo-chru-strasbourg/vcf2circos).

This package is based on [PCircos](https://github.com/CJinnny/PCircos) code


![Doc Circos](docs/dna_circos.png)
<br/>

# Installation

## Git clone

Download package source files.

```
$ git clone https://github.com/bioinfo-chru-strasbourg/vcf2circos.git .
```

## Pip

Compile source using Pip to generate binary "vcf2circos"

```
$ python -m pip install -e .
```

## Docker

Build docker image "vcf2circos:latest"

```
$ docker-compose build
```


<br/>

# Usage

## Binary 

```
$ vcf2circos --input config/Static/example.vcf.gz --options config/Static/options.json --output <outputpath>.html

```

## Docker

```
$ docker run -v $(pwd):/data vcf2circos:latest --input=demo_data/example.vcf.gz --output=/data/example.vcf.gz.html --options=demo_data/options.example.json

```

## Python 

```
$ python vcf2circos/__main__.py --input config/Static/example.vcf.gz --options config/Static/options.json --output <outputpath>.html

```

<br/>

# Help

```
usage: python vcf2circos.py [-h] -i INPUT -o OUTPUT [-e EXPORT] [-p OPTIONS]
                         [-n NOTEBOOK_MODE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file.
                        Format will be autodetected from file path.
                        Supported format:
                           'json', 'vcf'
  -o OUTPUT, --output OUTPUT
                        Output file.
                        Format will be autodetected from file path.
                        Supported format:
                           - html (dynamic): 'html'
                           - images: 'png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'eps'
                           - json (plotly): 'json'
  -e EXPORT, --export EXPORT
                        Export file.
                        Format is 'json'.
                        Generate json file from VCF input file
  -p OPTIONS, --options OPTIONS
                        Options file or string.
                        Format is 'json', either in a file or as a string.
```


<br/>

# Input

This package allows multiple input formats:
- VCF including SNV/InDel/SV (see [VCF specifications](https://samtools.github.io/hts-specs/VCFv4.2.pdf)). Header needs to contain contigs (in order of appearance). See [VCF example](demo_data/example.vcf.gz).
- JSON configuration file (see [PCircos](https://github.com/CJinnny/PCircos) documentation). See [JSON configuration example](demo_data/demo_params.json).

Format will be autodetected from file path.


<br/>

# Output

This package generates Circos plot in multiple formats (html, png, jpg, jpeg, webp, svg, pdf, eps, json):
- HTML file (format with customizable hover text). See [HTML example](docs/refontcircos.html)
- Image files (i.e. png, jpg, jpeg, webp, svg, pdf, eps). See [PNG example](docs/refontcircos.png) and [PDF example](docs/refontcircos.pdf)
- JSON Plotly file (see [Plotly documentation](https://plotly.com/)). 

Format will be autodetected from file path.

Output Circos plot sections from a VCF file:
![Doc Circos](docs/docs.circos.png)


<br/>

# Export

Circos plot generated from VCF file can be exported as JSON configuration file, for further use. See [JSON configuration export example](demo_data/example.vcf.gz.export.json)


<br/>

# Options 

Circos plot generated from a VCF file can be configured using a JSON options file. See [JSON options example](demo_data/options.example.json).

Here is an example of a JSON options file:
```
{
	"General": {
		"title": "",
		"width": 1200,
		"height": 1200,
		"plot_bgcolor": "white"
	},
	"Static": "config/Static",
	"Assembly": "hg19",
	"Chromosomes": {
		"cytoband": "True",
		"list": ["chr1", "chrX"]
	},
	"Genes": {
		"only_snv_in_sv_genes": false
	},
	"Variants": {
		"annotations": {
			"fields": ["SVTYPE", "SVLEN"]
		},
		"rings": {
			"position": 0.4,
			"height": 0.04,
			"space": 0.01,
			"nrings": 6
		}
	},
	"Extra": [
		"gc"
	]
}
```

<br/>

## Options format

Exemple of a data tab-delimited file (<b>STILL IN DEV</b>):
Overview of cytoband file, at terms it will be possible to add this kind of data above copy number level rings
```
chr_name  start     end       band_color  band
chr1      0         2300000   gneg        p36.33
chr1      2300000   5400000   gpos25      p36.32
chr1      5400000   7200000   gneg        p36.31
chr1      7200000   9200000   gpos25      p36.23
chr1      9200000   12700000  gneg        p36.22
chr1      12700000  16200000  gpos50      p36.21
chr1      16200000  20400000  gneg        p36.13
chr1      20400000  23900000  gpos25      p36.12
chr1      23900000  28000000  gneg        p36.11
```

<br/>

## General section

The "General" section is a Plotly General section, which configure main options of the Circos plot (e.g. title, size, back-ground color).

Example:
```
"General": {
    "title": "",
    "width": 1400,
    "height": 1400,
    "plot_bgcolor": "white"
}
```

<br/>

## Chromosomes section

The "Chromosomes" section defines information about chromosomes (e.g. contig, list of chromosomes).

Example:
```
"Chromosomes": {
    "list": ["chr7", "chr13", "chr12", "chr14", "chr15", "chrX", "chr1", "chr17"]
}
```


### List of chromosomes

The "list" option define the list of chromosomes to show in the Circos plot. Order of chromosome is still defined in the VCF header (in "contigs" section). If no chromosomes are listed, all chromosomes in the VCF header will be shown.

<br/>

## Genes section

The "Genes" section defines information about Genes (e.g. refGene data, list of genes to show). These information are used to annnotate variants (SNV and SV), and are used with algorithms highlight interesting information (e.g. only SNV on CNV genes). They also can be shown in the Circos plot (below Chromosomes ring)

```
    "only_snv_in_sv_genes": true
}
```
### List of genes

The "list" option defines the list of genes to show in the Circos plot, below Chromosomes/Cytoband ring. This list refers to the "gene" column in the data.

<br>
### Filter SNV on CNV genes

The "only_snv_in_sv_genes" option will select (and show) only SNV that are located on genes mutated with at least 1 SV.


<br/>

## Variants section

The "Variants" section defines varaints annotations to show in each variant hover text, and positions of the varaints rings.

Example:
```
"Variants": {
    "annotations": {
        "fields": ["SVTYPE", "SVLEN"],
    },
    "rings": {
        "position": 0.50,
        "height": 0.04,
        "space": 0.01,
        "nrings": 6
    }
}
```

### Annotations

The "annotations" option defines the annotations of variants to be shown.

The "fields" option configures the list of annotations in the hover text. If empty list if provided getting 15 first annotations in order of appearance in vcf info field. Moreover size of hover annotations is limited to 40 chars.

<br>

### Rings

The "rings" option defines the "position" and "height" of SNV and SV rings, "space" between rings and the number of ring in lightgray to display.



<br/>


# Contacts

 Medical Bioinformatics Applied to Diagnosis - Strasbourg University Hospital - France

[Website](https://www.chru-strasbourg.fr/service/bioinformatique-medicale-appliquee-au-diagnostic-unite-de/)

[GitHub](https://github.com/bioinfo-chru-strasbourg)

[bioinfo@chru-strasbourg.fr](mailto:bioinfo@chru-strasbourg.fr)


