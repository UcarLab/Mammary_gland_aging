#!/bin/sh

cd ~/path/to/directory

##For Analysis. Sample Counts and Metadata file provided. 
cellphonedb method statistical_analysis counts.txt meta.txt --counts-data=gene_name --project-name=Project
#For Dotplot. Provide the means and pvalues text file from the previous step.

## Edit means and pvalues (Check Select Interactions). Interactions were picked manually on occasion so there is no specific way of automating this process.
cellphonedb plot dot_plot --means-path ~/path/to/means.txt --pvalues-path ~/path/to/pvalues.txt --output-path ~/path/to/output --output-name CellphoneDB_Dotplot.pdf
