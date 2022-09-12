#!/bin/sh

cd ~/path/to/directory

##For Analysis. Sample Counts and Metadata file provided. 
cellphonedb method statistical_analysis counts.txt meta.txt --counts-data=gene_name --project-name=Project
#For Dotplot. Provide the means and pvalues text file from the previous step.
cellphonedb plot dot_plot --means-path ~/path/to/means.txt --pvalues-path ~/path/to/means.txt --output-path ~/path/to/output --output-name CellphoneDB_Dotplot.pdf
