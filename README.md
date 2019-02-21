[![Build Status](https://travis-ci.com/milescsmith/ReducedHeat.svg?branch=master)](https://travis-ci.com/milescsmith/ReducedHeat)
# ReducedHeat
Display gene expression along a given reduced dimension on a heatmap

Inspired by (and using code from) [t-SNE-Heatmaps](https://github.com/KlugerLab/t-SNE-Heatmaps) 
from the Kluger Lab, this creates a heatmap visualizing gene expression along a given dimensional 
reductions coordinates.  Unlike the example code given for t-SNE-Heatmaps, this works with any 
dimensional reduction (i.e. UMAP)

As an example, these are the marker genes (as determined by Seurat's FindAllMarkers() function, using MAST) 
for the cell types identified in the [Villani et al. Science 2017 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5775029/)
displayed along the first UMAP dimension:
![Villani et.al. plot](https://github.com/milescsmith/ReducedHeat/blob/master/newplot.png)

Currently works with Seurat 3 objects or a gene expression matrix.
