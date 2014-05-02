hicluster
=========

hierarchical clustering, PCA, RNA-seq distribution and heatmap rendering

This script is designed primarily for RNA-seq clustering analysis, and as such the
normalisation and filtering defaults are optimised for these kind of data. However, as
long as you are aware of the data manipulations required, the script should still be
capable of being run on other data.

The matrix being analysed can either be input as a tab- or space-delimited table, or
created from a list of cufflinks files. The format for a data table is as follows:

Gene:	alpha	beta	gamma	delta	epsilon	zeta	omega
Annie   1000    367     0       1312    0       400     2100
Bob     0       12      700     1234    700     12      12
Charlie	0       800     333     1378    333     760     14
Diego	0       200     666     1345    670     273     15
Eddie	0       153     700     1384    733     154     25
Faith	999     15      32      1389    33      12      2050
Gill	1000    330     0       1356    0       304     2012


Be sure to include the gene header row and sample header column. Values are fpkm


By providing a list of cufflinks files, a table can be constructed using the gene names
from column 1 and the gene values from column 9 of the cufflinks files. Sample names will
be the filenames.



usage: hicluster.py [-D <table.file>or -B <comma-delimited,list,of,cufflinks,files>] [optional arguments]

required (one of the following options is necessary):
  -D , --data_file      The data file for analysing (formatted as above)

  -B , --build_table    Provide a comma-delimited list of cufflinks files with
                        which to build the fpkm table for analysis.

optional arguments:
# output #
  -o , --output_file
                        output file name for results

# analysis #
  -z , --row_header
                        Not sure why this was put here
  -R , --row_method
                        The clustering method for rows (single, average,
                        complete, etc)
  -C , --column_method
                        The clustering method for columns (single, average,
                        complete, weighted, ward, centroid, etc)
  -r , --row_metric
                        The distance metric for rows (euclidean, correlation,
                        cosine, manhattan, etc)
  -c , --column_metric
                        The distance metric for columns (euclidean,
                        correlation, manhattan, etc)
  -P, --pca             Performs principal component analysis.
  -a , --anova
                        Perform ANOVA on 10 groups, and report genes with
                        P-value less than value specified
  -A , --filter_anova
                        Perform ANOVA and filter genes to keep only those with
                        P-value less than value given
  -s , --t_test
                        Perform student's t-test on two groups, and report
                        genes with P-value less than value specified
  -S , --filter_t
                        Perform student's t-test and filter genes to keep only
                        those with P-value less than value given

# viewing #
  -g , --color_gradient
                        The colour scheme (red_white_blue, red_black_sky,
                        red_black_blue, red_black_green, yellow_black_blue,
                        seismic, green_white_purple, coolwarm)
  -d, --distribution    Shows FPKM distribution of each sample before and
                        after normalisation
# filtering #
  -m , --magnitude
                        Filters out genes with magnitude of range less than
                        value given. Default = 2.5
  -L , --gene_list
                        Allows provision of a file containing a list of genes
                        for inclusion in the clustering (ie, will filter out
                        all genes NOT in the list provided). Otherwise, all
                        genes will be included.
  -p , --fpkm_max
                        Filters out genes with maximum fpkm greater than value
                        given. Default = 1 000 000
  -q , --fpkm_min
                        Filters out genes with maximum fpkm less than value
                        given. Default = 10
  -f, --filter_off      Turns off filtering.
  -y , --pretty_filter
                        Filters (non-normalised) matrix to remove genes whose
                        mean value between treatments is less than value
                        given. Try 2.5
  -k, --kill_PC1        removes first principal component

# data transformation #
  -t, --transform_off   Turns off log2(FPKM + 1) transformation (prior to
                        normalisation if selected).
  -n, --normalise_off   Turns off normalisation. Normalises by dividing by the
                        standard deviation.
  -u, --centering_off   Turns off gene centering. Centering subtracts the mean
                        from all values for a gene, giving mean = 0.
  -T, --transpose       Transpose the matrix. Columns should represent genes,
                        Rows samples
  -X, --sample_norm     Normalises samples instead of genes


clustering and visualisation code based on that by Nathan Salomonis - nsalomonis@gmail.com
(downloaded from http://code.activestate.com/recipes/578175-hierarchical-clustering-heatmap-python/
on April 10th, 2014)

The first commit of the python script is the original code by N. Salomonis.