#!/bin/bash
# modified from from https://gist.github.com/matsen/4263955
cd ../figures
# figure 1
inkscape --without-gui --export-png="useiggs.png" --export-dpi 900 useiggs.svg
convert -compress LZW -alpha remove useiggs.png useiggs.tiff
mogrify -alpha off useiggs.tiff
mv useiggs.tiff numbered/figure1.tiff
rm useiggs.png

# figure 2
convert -compress LZW -alpha remove wunifrac_Agp41_Bgp120_pcoa.png numbered/figure2.tiff

# figure 3

convert -compress LZW -alpha remove stacked_bars.png numbered/figure3.tiff

# figure 4

convert -compress LZW -alpha remove phi_heatmap_withlegend.png numbered/figure4.tiff

cd ../scripts

