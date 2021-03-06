#!/usr/bin/env python

""" Imports a whitespace-delimited expression matrix and produces a hierarchically
clustered heatmap. Can also build the expression matrix from cufflinks data file.

Multiple filtering options and matrix normalisation methods are now implemented.

T-test and ANOVA testing are also available to test individual rows/columns.

PCA and read distribution plots can be displayed to check for normalcy, clustering etc.



"""
import string
import time
import sys, os, re
import getopt
import argparse
from multiprocessing import Queue, Process
from operator import itemgetter
from itertools import chain, repeat



import numpy
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy import stats
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from scipy import interpolate
from statsmodels.stats.multicomp import pairwise_tukeyhsd


from genomepy import genematch

####################################################################################
# ############################### CLASSES ######################################## #
# #                                                                              # #
# #                                                                              # #

class Permutable(object):
    def __init__(self, cluster):
        if cluster._genes_as_rows:
            cluster.invert_matrix()
        self._genes_as_rows = False
        self.data_matrix = numpy.copy(cluster.data_matrix)
        self.sample_header = cluster.sample_header[:]
        self.gene_header = cluster.gene_header[:]
        self.column_header = cluster.column_header[:]
        self.row_header = cluster.row_header[:]
        self.column_metric = cluster.row_metric
        self.genenumber = cluster.genenumber
        self.samplesize = cluster.samplesize

    def randomise_samples(self):
        neworder = numpy.random.shuffle(self.sample_header)

    def reorder_matrix(self, groups=["_FL","_SP"], verbose=False):
        "reorders rows such that they are sorted accordign to specified groups"

        if verbose:
            print "Sorting data into %d groups" % (len(groups))

        # split matrix based on the specified groups;
        namelist = {}
        poslist = {}

        removedlist = []

        for s in self.sample_header[:]:
            for pattern in groups:
                if pattern not in namelist:
                    namelist[pattern] = []
                    poslist[pattern] = []
                # assign each sample to a group:
                if re.search(pattern,s) is not None:
                    namelist[pattern].append(s)
                    poslist[pattern].append(self.sample_header.index(s))
                    break
            else:       # group not found in any sample list
                removedlist.append(s)
                self.remove_sample(s)
        if len(removedlist) > 0:
            print("The following %d samples could not be matched to any group and were removed ==> \n%s\nGroups: %s"
                     % (len(removedlist), " ".join(removedlist), " ".join(groups)))
        grouporder = []
        limits = []
        boundarystone = 0
        for pattern in groups:
            try:
                if verbose:
                    print pattern, namelist[pattern]
                    print poslist[pattern]
                grouporder += poslist[pattern]
                boundarystone += len(poslist[pattern])
            except KeyError:
                print pattern, "None found!"
            limits.append(boundarystone)

        matrix_reord = self.data_matrix[grouporder,:]
        row_header_new = [self.sample_header[i] for i in grouporder]

        if verbose:
            print self.sample_header
            print "grouporder       :", grouporder
            print row_header_new
            print "limits (boundary):", limits


        self.data_matrix = matrix_reord
        self.row_header = row_header_new
        self.sample_header = row_header_new

        return limits

class Cluster(object):
    def __init__(self, datafile, exportPath=os.getcwd(), firstrow=True, genes_as_rows=False, \
                gene_metric='correlation', sample_metric='correlation', \
                gene_method='complete', sample_method='complete', \
                color_gradient='red_white_blue'):
        """Creates a matrix from the designated file, and assigns variables as appropriate
        """
        self.importData(datafile, firstrow)

        self.gene_metric    = gene_metric
        self.sample_metric  = sample_metric
        if gene_method == 'none':
            self.gene_method = None
        else:
            self.gene_method    = gene_method
        if sample_method == 'none':
            self.sample_method = None
        else:
            self.sample_method  = sample_method

        self._genes_as_rows = genes_as_rows
        self.refresh_headers()
        self.clean_header()

        if genes_as_rows:
            self.invert_matrix()

        self.samplesize, self.genenumber = self.check_size()

        self.vmax, self.vmin= self.getColorRange()

        if color_gradient == 'red_white_blue':
            self.cmap = plt.cm.bwr
        if color_gradient == 'red_black_sky':
            self.cmap = self._RedBlackSkyBlue()
        if color_gradient == 'red_black_blue':
            self.cmap = self._RedBlackBlue()
        if color_gradient == 'red_black_green':
            self.cmap = self._RedBlackGreen()
        if color_gradient == 'yellow_black_blue':
            self.cmap = self._YellowBlackBlue()
        if color_gradient == 'seismic':
            self.cmap = plt.cm.seismic
        if color_gradient == 'green_white_purple':
            self.cmap = plt.cm.PiYG_r
        if color_gradient == 'coolwarm':
            self.cmap = plt.cm.coolwarm

        self.averaged = False
        self.exportPath = exportPath

    def __str__(self):
        return "Data matrix of %d Samples and %d Genes" % (self.samplesize, self.genenumber)

    __rep__ = __str__

    def importData(self, filename, first_row=True):
        start_time = time.time()
        matrix=[]
        row_header=[]

        dataset_name = os.path.basename(filename)

        filename_h = open(filename, 'rb')
        for line in filename_h:
            t = line[:-1].split() ### remove end-of-line character - file is tab- or space-delimited
            if first_row:
                column_header = t[:]    # check later to see if first element needs to be removed
                first_row=False
            else:
                if 'X' not in t and '-' not in t: ### Occurs for rows with missing data
                    try:
                        s = map(float,t[1:])
                    except: # most common error for new analysis - creating a column containing cufflinks header names.
                        for element in t[1:]:
                            try:
                                float(element)
                            except ValueError:
                                print "Error importing table. Expected value, got %r instead!" % element
                    if (abs(max(s)-min(s)))>0: # only add to matrix if not all zeros
                        matrix.append(s)
                        row_header.append(t[0])

        # check appropriate number of items added to column header:
        if len(matrix[0]) == len(column_header) - 1:
            column_header = column_header[1:]

        time_diff = str(round(time.time()-start_time,1))
        try:
            verbalise('\n%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff))
        except Exception:
            verbalise("R", 'No data in input file.'); force_error

        self.data_matrix    = numpy.array(matrix)
        self.column_header  = column_header
        self.row_header     = row_header

    def clean_header(self):
        "removes long path name from row headers"

        new_headers = []
        for name in self.sample_header:
            try:
                newname = re.findall('/([^/]*)', name)[-2]
            except:
                newname = name
            new_headers.append(newname)

        if self._genes_as_rows:
            self.column_header = new_headers
        else:
            self.row_header = new_headers
        self.refresh_headers()

    def refresh_headers(self):
        if self._genes_as_rows:
            self.sample_header = self.column_header
            self.gene_header   = self.row_header
            self.row_metric    = self.gene_metric
            self.column_metric = self.sample_metric
            self.row_method    = self.gene_method
            self.column_method = self.sample_method
        else:
            self.sample_header = self.row_header
            self.gene_header   = self.column_header
            self.column_metric = self.gene_metric
            self.row_metric    = self.sample_metric
            self.column_method = self.gene_method
            self.row_method    = self.sample_method

    def invert_matrix(self):
        self.data_matrix = numpy.transpose(self.data_matrix)
        tempcol = self.column_header
        temprow = self.row_header
        tempcolmeth = self.column_method
        temprowmeth = self.row_method
        tempcolmet = self.column_metric
        temprowmet = self.row_metric

        self.column_header  = temprow
        self.row_header     = tempcol
        self.column_method  = temprowmeth
        self.row_method     = tempcolmeth
        self.column_metric  = temprowmet
        self.row_metric     = tempcolmet

        if self._genes_as_rows:
            self._genes_as_rows = False
        else:
            self._genes_as_rows = True

        self.refresh_headers()

    def set_genes_to_rows(self):
        if self._genes_as_rows is False:
            self.invert_matrix()

    def set_genes_to_columns(self):
        if self._genes_as_rows:
            self.invert_matrix()

    ## Create Custom Color Gradients
    #http://matplotlib.sourceforge.net/examples/plt_examples/custom_cmap.html
    def _RedBlackSkyBlue(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'green': ((0.0, 0.0, 0.9),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0)),

                 'blue':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }

        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    def _RedBlackBlue(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'blue':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }

        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    def _RedBlackGreen(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'blue': ((0.0, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }

        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    def _YellowBlackBlue(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'green': ((0.0, 0.0, 0.8),
                           (0.5, 0.1, 0.0),
                           (1.0, 1.0, 1.0)),

                 'blue':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }
        ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
        ### modulate between blue and cyan using the last y var in the first green tuple
        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    def getColorRange(self):
        """ Determines the range of colors, centered at zero, for normalizing cmap """
        vmax=self.data_matrix.max()
        vmin=self.data_matrix.min()

        if vmax * vmin < 0:     # ie number range spans +ve and -ve
            vmax = max([vmax, abs(vmin)])
            vmin = -1*vmax

        return vmax,vmin

    def check_size(self):
        if self._genes_as_rows:
            self.samplesize     = len(self.data_matrix[0])
            self.genenumber     = len(self.data_matrix)
        else:
            self.samplesize     = len(self.data_matrix)
            self.genenumber     = len(self.data_matrix[0])
        assert len(self.gene_header) == self.genenumber
        assert len(self.sample_header) == self.samplesize
        return self.samplesize, self.genenumber

    ################# Matrix manipulation methods ######################################

    def normaliseData(self, center=True, norm_var=True, log_t=True, sample_norm=False):
        "center, normalize and/or log transform the array x"

        verbalise("Y", "\n\nNormalising data according to input parameters.")
        verbalise("G", "Initial value range in matrix:       %s - %-4.3f" % ('{0: 3.3f}'.format(self.data_matrix.min()), self.data_matrix.max()))
        if log_t:
            count = 0
            k = self.data_matrix.min()
            for g in range(self.genenumber):
                for i in range(self.samplesize):
                    #if self.data_matrix[:,g][i] <= -1:
                    #    print self.data_matrix[:,g][i]
                    #if 0 < g < 3 and 0 < i < 10:
                    #    print "gene %d sample %d: %r (%s) %r" % (g, i, self.data_matrix[:,g][i], type(self.data_matrix[:,g][i]), numpy.log2(self.data_matrix[:,g][i] + 1))
                    self.data_matrix[:,g][i] = numpy.log2(self.data_matrix[:,g][i] - k + 1)

            verbalise("G", "log2(FPKM + k) transformed. New range: %-4.3f - %-4.3f" % (self.data_matrix.min(), self.data_matrix.max()))

        meanlist = []   # to store for later re-adjustment of matrix
        if center and sample_norm:
            for s in range(self.samplesize):
                ave_g = numpy.mean(self.data_matrix[s])
                for i in range(self.genenumber):
                    self.data_matrix[s,i] = self.data_matrix[s,i] - ave_g
                    meanlist.append(ave_g)
        elif center:
            # Centering data for each gene to have mean 0:
            for g in range(self.genenumber):
                # center the data so it has mean = 0 for each gene:
                ave_g = numpy.mean(self.data_matrix[:,g])
                for i in range(self.samplesize):
                    self.data_matrix[:,g][i] = self.data_matrix[:,g][i] - ave_g
                    meanlist.append(ave_g)
            verbalise("G", "Centered data. New value range:      %-4.3f - %-4.3f" % (self.data_matrix.min(), self.data_matrix.max()))

        if norm_var and sample_norm:
            # Normalising data so each *sample* has standard deviation of 1:
            for s in range(self.samplesize):
                #sum_sq = sum([ a ** 2 for a in (x[m])])
                stdev = numpy.std(self.data_matrix[s])
                for i in range(self.genenumber):
                    self.data_matrix[s,i] = self.data_matrix[s,i] / stdev
        elif norm_var:
            # Normalising data so each *gene* has standard deviation of 1"
            for g in range(self.genenumber):
                stdev = numpy.std(self.data_matrix[:,g])
                for i in range(self.samplesize):
                    self.data_matrix[:,g][i] = self.data_matrix[:,g][i] / stdev # dividing by variance gives identical result
            verbalise("G", "Normalised St Dev. New value range:  %-4.3f - %-4.3f" % (self.data_matrix.min(), self.data_matrix.max()))
        verbalise("G", "Final value range in matrix:         %s - %-4.3f" % ('{0: 3.3f}'.format(self.data_matrix.min()), self.data_matrix.max()))

        return meanlist

    def filter_genes(self, genelist):
        """takes a file and extracts the genes in the col_num specified column,
        then uses it to filter the matrix"""

        # create position-based list for filtering from gene-names:
        keeplist = []  # will be a position based list [ 0, 4, 23, 34, ... ]
        for gene in genelist:
            try:
                keeplist.append(self.gene_header.index(gene))
            except:     # some columns will be headers that would cause errors.
                verbalise("R", gene, "not found. Not keeping.")

        hitlist = range(len(self.gene_header))
        err_ct = 0
        for posn in keeplist: # will be the inverse of keeplist
            try:
                hitlist.remove(posn)
            except Exception as inst:
                err_ct += 1
        if err_ct > 0:
            verbalise("R", "There were %d errors encountered. Last error:\n" % (err_ct, inst))
        ## filter matrix with genelist:
        self.filter_matrix(hitlist)

    def remove_sample(self, sample):
        posn = self.sample_header.index(sample)

        if self._genes_as_rows:
            self.data_matrix = numpy.delete(self.data_matrix, posn, 1)
            print self.sample_header
            print self.column_header
            s1 = self.sample_header.pop(posn) # sample_header and column_header point to
            # the same object. Removing the name from one removes it from both!


        else:
            self.data_matrix = numpy.delete(self.data_matrix, posn, 0)
            s1 = self.sample_header.pop(posn)

        self.check_size()

    def filterData(self, mag=1, min_thresh=-1000000, max_thresh=1000000):
        """filters out any gene for which the magnitude of expression is less than mag,
        or for which the maximum value of all samples is less than min_thresh,
        or for which the maximum value is greater than max_thresh"""

        # will only work if genes are columns in matrix
        revert = False
        if self._genes_as_rows:
            self.invert_matrix()
            revert = True

        verbalise("G", "Filtering %d genes and %d samples:\nMin fold change: %.1f Min expression level (at least one sample): %d Max expression level: %d" % (self.genenumber, self.samplesize, mag, min_thresh, max_thresh))

        hitlist = []
        for g in range(self.genenumber):
            fpkm_max = max(self.data_matrix[:,g])
            fpkm_min = min(self.data_matrix[:,g])
            size = numpy.absolute(fpkm_max/(fpkm_min + 0.00001))
            #rms  = numpy.sqrt( sum(a**2 for a in x[:,g])/m )
            if size < mag or fpkm_max < min_thresh or fpkm_max > max_thresh :
                hitlist.append(g)

        self.filter_matrix(hitlist)

        # if matrix was inverted for gene removal, restore to its previous orientation:
        if revert:
            self.invert_matrix()

    def filter_matrix(self, hitlist):
        """given the positional list hitlist, will edit matrix and column headers
        by removing those positions"""
        # will only work if genes are columns in matrix
        revert = False
        if self._genes_as_rows:
            self.invert_matrix()
            revert = True

        # create new matrix and column_header without the columns in the hitlist:
        y = numpy.delete(self.data_matrix, hitlist, 1)
        verbalise("%d genes removed:" % len(hitlist))
        for gene in hitlist[::-1]:
            self.gene_header.pop(gene)

        self.column_header = self.gene_header
        self.data_matrix = y
        #verbalise("G", "The shape of the new matrix is", numpy.shape(y))
        self.samplesize, self.genenumber = self.check_size()
        verbalise("G", "there are now %d genes and %d samples" % (self.genenumber, self.samplesize))



        # if matrix was inverted for gene removal, restore to its previous orientation:
        if revert:
            self.invert_matrix()

    def reorder_matrix(self, groups=["_FL","_SP"], verbose=False):
        "reorders rows such that they are sorted accordign to specified groups"

        # will only work if genes are columns in matrix
        revert = False
        if self._genes_as_rows:
            self.invert_matrix()
            revert = True


        if verbose:
            verbalise("Y", "Sorting data into %d groups" % (len(groups)))

        # split matrix based on the specified groups;
        namelist = {}
        poslist = {}

        removedlist = []

        for s in self.sample_header[:]:
            for pattern in groups:
                if pattern not in namelist:
                    namelist[pattern] = []
                    poslist[pattern] = []
                # assign each sample to a group:
                if re.search(pattern,s) is not None:
                    namelist[pattern].append(s)
                    poslist[pattern].append(self.sample_header.index(s))
                    break
            else:       # group not found in any sample list
                removedlist.append(s)
                self.remove_sample(s)
        if len(removedlist) > 0:
            verbalise("R", "The following %d samples could not be matched to any group and were removed ==> \n%s\nGroups: %s"
                     % (len(removedlist), " ".join(removedlist), " ".join(groups)))
        grouporder = []
        limits = []
        boundarystone = 0
        for pattern in groups:
            try:
                if verbose:
                    print pattern, namelist[pattern]
                    print poslist[pattern]
                grouporder += poslist[pattern]
                boundarystone += len(poslist[pattern])
            except KeyError:
                verbalise("R", pattern, "None found!")
            limits.append(boundarystone)


        matrix_reord = self.data_matrix[grouporder,:]
        row_header_new = [self.sample_header[i] for i in grouporder]

        if verbose:
            print self.sample_header
            print "grouporder       :", grouporder
            print row_header_new
            print "limits (boundary):", limits


        self.data_matrix = matrix_reord
        self.row_header = row_header_new
        self.sample_header = row_header_new

        # if matrix was inverted for gene removal, restore to its previous orientation:
        if revert:
            self.invert_matrix()

        return limits

    def get_stats(self, gene, groups=["SP", "SL06", "SL12", "SL24","SL48", "SL96","FP06", "FP12", "FP24","FP48", "FP96", "FL"]):
        """
        For the given gene, returns lists of mean and standard deviation for each group.
        """

        try:
            pos = self.gene_header.index(gene)
        except:
            return [],[]

        limits = [0] + self.reorder_matrix(groups=groups)
        #print limits
        intervals = zip(limits[:-1],limits[1:])
        print "Debugging: ave and stdev:-"
        #print intervals
        #print self.data_matrix
        #print self.data_matrix[:,pos]
        #print [ (self.data_matrix[ i:j,pos ]) for i,j in intervals]
        ave = [ numpy.average(self.data_matrix[ i:j,pos ]) for i,j in intervals]
        print ave
        stdev = [ numpy.std(self.data_matrix[ i:j,pos ]) for i,j in intervals]
        print stdev
        return ave, stdev


    def average_matrix(self, groups=["SP", "SL06", "SL12", "SL24","SL48", "SL96","FP06", "FP12", "FP24","FP48", "FP96", "FL"]):
        """for each group in variable groups, calculates the average value in the matrix, and
        creates a new matrix showing average values"""

        # will only work if genes are columns in matrix
        revert = False
        if self._genes_as_rows:
            self.invert_matrix()
            revert = True

        limits = self.reorder_matrix(groups=groups)

        v = numpy.ones((len(groups),self.genenumber)) # create new array of size (y groups) and (n genes)

        for g in range(self.genenumber):
            v[0,g] = numpy.average(self.data_matrix[:limits[0],g])
            for i in range(len(groups)-1):
                v[i + 1,g] = numpy.average(self.data_matrix[limits[i]:limits[i + 1],g])

        self.row_header = groups
        self.column_header = self.gene_header
        self.data_matrix = v
        self.refresh_headers()

        # if matrix was inverted for gene removal, restore to its previous orientation:
        if revert:
            self.invert_matrix()
        self.averaged = True

    def export(self):
        "Saves datatable and columns to a text file"
        nfilename = self.exportPath[:-4] + '.data.tbl'

        export_text = open(nfilename,'w')
        column_header = '\t'.join(['gene']+self.column_header)+'\n' ### format column-names for export
        export_text.write(column_header)

        ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
        new_row_header = self.row_header[::-1]
        xt = self.data_matrix[::-1]

        ### Export each row in the clustered data matrix xt
        i = 0
        for row in xt:
            export_text.write('\t'.join([new_row_header[i]]+map(str, row))+'\n')
            i += 1
        export_text.close()

        ### Transpose text file for easier reading!
        oldfile_h = open(nfilename, 'rb')

        elements = [ line.split() for line in oldfile_h ]
        oldfile_h.close()

        biglist = []
        for splitline in elements:
            biglist.append(splitline)
        newarray = numpy.array(biglist)
        t_array = newarray.transpose()

        newfile_h = open(self.exportPath[:-4] + '.data.trans.tbl' , 'w')
        for row in t_array:
            newfile_h.write("\t".join(row) + "\n")
        newfile_h.close()

    def get_values(self, geneid):
        gindex = self.gene_header.index(geneid)
        if self._genes_as_rows:
            values = self.data_matrix[gindex]
        else:
            values = self.data_matrix[:,gindex]
        print geneid
        for sample, value in zip(self.sample_header, values):
            print "%-25s %.2f" % (sample, value)

    def randomise_samples(self):
        self.oldref = self.sample_header[:]
        numpy.random.seed(int(time.time()))
        neworder = numpy.random.shuffle(self.sample_header)

# #                                                                              # #
# #                                                                              # #
# ############################## END CLASSES ##################################### #
####################################################################################


####################################################################################
# ############################### FUNCTIONS ###################################### #
# #                                                                              # #
# #                                                                              # #

def define_arguments():
    parser = argparse.ArgumentParser(description="Performs heirarchical clustering")

    ### input options ###
    parser.add_argument("-D", "--data_file", type=str, help="The data file for analysing")
    parser.add_argument("-B", "--build_table", type=str, dest="build_list",
                        default=None, help="Provide a comma-delimited list of cufflinks files with which to build the fpkm table for analysis.")
    parser.add_argument("-T", "--genes_as_rows", action='store_true',
                        help="Select if top row of data file is sample names")
    ### output options ###
    parser.add_argument("-o", "--output_file", dest='filename', type=str, default=None,
                        help="output file path and root name for results")
    parser.add_argument("-e", "--export_table", action='store_true',
                        help="export transformed expression matrix")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="print descriptive output while processing")
    ### analysis options ###
    parser.add_argument("--groups", type=str,
                        default='SP,SL06,SL12,SL24,SL48,SL96,FL,FP06,FP12,FP24,FP48,FP96',
                        help="comma delimited list of (regex) groups for categorising samples.")
    parser.add_argument("--no_clustering", action='store_true',
                        help="Turn of clustering. Performs data transformation, filtering and analysis, then exits")
    parser.add_argument("-R", "--sample_method", type=str, default='complete',
            help="The clustering method for samples \n(single, average, complete, etc)")
    parser.add_argument("-G", "--gene_method", type=str, default='complete',
                        help="The clustering method for genes \n(single, average, complete, weighted, ward, centroid, etc)")
    parser.add_argument("-r", "--sample_metric", type=str, default='correlation',
        help="The distance metric for samples \n(euclidean, correlation, cosine, manhattan, etc)")
    parser.add_argument("-g", "--gene_metric", type=str, default='correlation',
        help="The distance metric for genes \n(euclidean, correlation, manhattan, etc)")
    parser.add_argument("-P", "--pca", action='store_true',
                        help="Performs principal component analysis.")
    parser.add_argument("-a", "--anova", type=float,
        help="Perform ANOVA on 10 groups, and report genes with q-value less than value specified")
    parser.add_argument("-A", "--filter_anova", type=float,
        help="Perform ANOVA and filter genes to keep only those with q-value less than value given")
    parser.add_argument("-s", "--t_test", type=float,
                        help="Perform student's t-test on two groups, and report genes with q-value less than value specified")
    parser.add_argument("-S", "--filter_t", type=float, dest="ttest_thresh",
                        help="Perform student's t-test and filter genes to keep only those with q-value less than value given")
    parser.add_argument("-K", "--kegg", type=float,
                        help="Perform KEGG pathway enrichment analysis on gene clusters. Outputs KEGG terms less than specified q-value")
    parser.add_argument("-E", "--go_enrichment", type=float,
                        help="Perform GO term enrichment analysis on gene clusters. Outputs GO terms less than specified q-value")
    parser.add_argument("--neighbours", type=str,
                        help="returns a list of the N closest neighbours to specified genes.")
    parser.add_argument("-N", "--num_neighbours", type=int, default=10,
                        help="specify the number of nearest neighbours to return. Default is 10.")
    parser.add_argument("--GSEA", type=int, help='performs gene set enrichment analysis')
    parser.add_argument("--processes", type=int, default=1,
            help="specify the number of processes to split the permutations over. Default is 1")
    parser.add_argument("--irizarry", action='store_true',
                        help="performs gene set enrichment analysis according to Irizarry's scale sensitive chi square statistic")
    parser.add_argument("--irizarry_z", action='store_true',
                help="performs gene set enrichment analysis according to Irizarry's z statistic")
    ### viewing options ###
    parser.add_argument("-c", "--color_gradient", type=str, dest="color_gradient",
                        default='red_white_blue', help="The colour scheme \n(red_white_blue, red_black_sky, red_black_blue, \nred_black_green, yellow_black_blue, seismic, \ngreen_white_purple, coolwarm)")
    parser.add_argument("-d", "--distribution", action='store_true', default=False,
                    help="Shows FPKM distribution of each sample before and after normalisation")
    parser.add_argument("--display_off", action='store_true',
                        help="Turn of displaying of clustering")
    parser.add_argument("--display_reverse", action='store_true',
                        help="Switch gene and sample axes on graphs")
    parser.add_argument("--bar_charts", type=str,
    help="List of genes for which you wish to create column graphs of average gene expression for.")
    ### filtering options ###
    parser.add_argument("-m", "--magnitude", type=float, dest="filter", default=2.5,
        help="Filters out genes with magnitude of range less than value given. Default = 1.13")
    parser.add_argument("-L", "--gene_list", type=str, dest="gene_list", default=None, help="Allows provision of a file containing a list of genes for inclusion in the clustering (ie, will filter out all genes NOT in the list provided). Otherwise, all genes will be included.")
    parser.add_argument("-p", "--fpkm_max", type=float, dest="fpkm_max", default=1000000,
        help="Filters out genes with maximum fpkm greater than value given. Default = 1 000 000")
    parser.add_argument("-q", "--fpkm_min", type=int, dest="fpkm_min", default=10,
        help="Filters out genes with maximum fpkm less than value given. Default = 10")
    parser.add_argument("-f", "--filter_off", action='store_true', help="Turns off filtering based on expression value. ")
    parser.add_argument("--kill_PC1", action='store_true', help="removes first principal component")
    parser.add_argument("--expression_peaks", type=float,
    help="finds genes with singular timepoint peaks of magnitude n greater than control timepoint")
    ### data transform options ###
    parser.add_argument("--randomise_samples", action='store_true', default=False,
                        help="Randomises samples immediately after cluster creation")
    parser.add_argument("-t", "--transform_off", action='store_true', default=False,
            help="Turns off log2(FPKM + 1) transformation (prior to normalisation if selected).")
    parser.add_argument("-n", "--normalise_off", action='store_true', default=False,
            help="Turns off normalisation. Normalises by dividing by the standard deviation.")
    parser.add_argument("-u", "--centering_off", action='store_true', default=False,
                        help="Turns off gene centering. Centering subtracts the mean from all values for a gene, giving mean = 0.")
    parser.add_argument("-X", "--sample_norm", action='store_true',
                        help='Normalises samples instead of genes')
    parser.add_argument("--show_averages", action='store_true',
                        help="Calculates the average value for each group and clusters based on this new matrix")

    return parser

def verbalise(arg1, *args):
    # define escape code: '\x1b[31m  %s  \x1b[0m'
    colordict = {'R':'\x1b[31m', 'G':'\x1b[32m',
         'Y':'\x1b[33m' ,'B':'\x1b[34m', 'M':'\x1b[35m' , 'C':'\x1b[36m' }
    if arg1 in colordict:
        argstring = " ".join([str(arg) for arg in args])
        if sys.stdout.isatty():
            color_code = colordict[arg1]
            end_color = '\x1b[0m'
        else:
            color_code = ""
            end_color = ""
    else:
        argstring = " ".join([arg1] + [arg for arg in args])
        color_code = ""
        end_color = ""

    print "%s%s%s" % (color_code, argstring, end_color)

def check_verbose(v=True):
    "allow optional printing with color conversion capability!"
    global verbalise
    if v:
        verbalise = verbalise
    else:
        verbalise = lambda *a: None

    return verbalise

def run_arguments(args):

    vfunc = check_verbose(True)

    ## create data table from cufflinks files:
    if args.build_list:
        data_table = create_table(args.build_list)
    else:
        data_table = os.path.realpath(args.data_file)

    ## create output folder and log file of arguments:
    timestamp = time.strftime("%b%d_%H.%M")
    if not args.filename:
        root_dir = os.path.dirname(data_table)
        newfolder = root_dir + "/hicluster." + timestamp
        os.mkdir(newfolder)  # don't have to check if it exists, as timestamp is unique
        filename = newfolder + "/hicluster." + timestamp + ".log"
    else:
        newfolder = os.path.realpath(args.filename)
        if os.path.exists(newfolder) is False:  # check to see if folder already exists...
            os.mkdir(newfolder)
        filename = newfolder + '/' + os.path.basename(args.filename) + '.' + timestamp + ".log"

    log_h = open(filename, 'w')
    log_h.write( "File created on %s\n" % (timestamp) )
    log_h.write( "Program called from %s\n\n" % (os.getcwd()) )
    for arg in str(args)[10:-1].split():
        log_h.write( "%s\n" % (arg) )
    log_h.close()

    ## create matrix:
    cluster = Cluster(data_table, exportPath=filename, firstrow=True, genes_as_rows=args.genes_as_rows, \
                gene_metric=args.gene_metric, sample_metric=args.sample_metric, \
                gene_method=args.gene_method, sample_method=args.sample_method, \
                color_gradient=args.color_gradient)

    if args.randomise_samples:
        cluster.randomise_samples()

    ####### FILTERING OF MATRIX ###########


    ## filter by provided gene list
    if args.gene_list:
        genelist = make_a_list(args.gene_list)
        print "Filtering to keep %d genes" % len(genelist)
        cluster.filter_genes(genelist)

    ## filter by magnitude, minimum values and maximum values
    if not args.filter_off:
        cluster.filterData(mag=args.filter, min_thresh=args.fpkm_min, max_thresh=args.fpkm_max)


    ####### NORMALISATION OPTIONS ##############

    ## show gene distribution
    if args.distribution:
        expression_dist(cluster)

    ## normalise matrix
    meanlist = cluster.normaliseData(center=not(args.centering_off), norm_var=not(args.normalise_off), log_t=not(args.transform_off), sample_norm=args.sample_norm)

    ## show distribution after normalisation
    if args.distribution:
        expression_dist(cluster, min_x=-2, max_x=2)

    ## export transformed data table:
    if args.export_table:
        cluster.export()

    if args.expression_peaks:
        print "\nChecking for expression peaks..."
        peaklist = expression_peaks(cluster, args.expression_peaks)
        cluster.filter_genes(peaklist)
        print "%d genes with expression peaks found" % (cluster.genenumber)
        out_h = open(filename[:-4] + ".expression_peaks.list", 'w')
        for gene in peaklist:
            out_h.write( "%-20s %s\n" % (gene, peaklist[gene]) )
        out_h.close()


    ####### ANALYSIS OPTIONS ##############

    ## Define groups for t-test/ANOVA etc
    groups = args.groups.split(',')

    ## Gene Set Enrichment Analysis (cf Subramanian et al. (2005) PNAS 102(43)
    if args.GSEA:
        nes, pvalues, qvalues = gene_set_enrichment(cluster, permutations=args.GSEA, processes=args.processes, display_on=not(args.display_off))
        out_h = open(filename[:-4] + ".GSEA.list", 'w')
        # write all positive scores with q values less than 0.25 to file:
        for pathway, qvalue in sorted(qvalues.iteritems(), key=itemgetter(1)):
            if qvalue <= 0.25 and nes[pathway] >= 0:
                out_h.write( "%-15s %-50s%s ES=%+6.3f P=%8.3e q=%8.4f\n" % (
                    pathway[0], pathway[1][:50], '...' if len(pathway[1])>50 else '   ',
                    nes[pathway], pvalues[pathway], qvalue))

        # write all negative scores with q values less than 0.25 to file:
        for pathway, qvalue in sorted(qvalues.iteritems(), key=itemgetter(1), reverse=True):
            if qvalue <= 0.25 and nes[pathway] < 0:
                out_h.write( "%-15s %-50s%s ES=%+6.3f P=%8.3e q=%8.4f\n" % (
                    pathway[0], pathway[1][:50], '...' if len(pathway[1])>50 else '   ',
                    nes[pathway], pvalues[pathway], qvalue))

        out_h.close()
        print "GSEA results (FDR < 0.25) have been saved to ", filename[:-4] + ".GSEA.list"

    if args.irizarry:
        print "Calculating Irizarry gene enrichment..."
        t_dict = find_degs(cluster, groups=groups)
        pvalues = {}
        paths = {}
        print "Calculating pathway enrichment scores..."
        paths, smallpaths = genematch.collect_kegg_pathways(minsize=10)
        paths.update(smallpaths)
        paths.update(genematch.collect_ipr_pathways(minsize=10))
        paths.update(genematch.collect_go_pathways(minsize=10) )
        for pathway in paths:
            pvalues[pathway] = irizarry_enrichment(t_dict, paths[pathway])
        out_h = open(filename[:-4] + ".GSEA_Irizarry.list", 'w')
        for pathway,pval in sorted(pvalues.iteritems(), key=itemgetter(1)):
            if pvalues[pathway]:
                out_h.write( "%-15s %-50s%s P=%8.4f\n" % (
                    pathway[0], pathway[1][:50], '...' if len(pathway[1])>50 else '   ',
                    pval))
        out_h.close()

    if args.irizarry_z:
        print "Calculating Irizarry gene enrichment..."
        t_dict = find_degs(cluster, groups=groups)
        #s_dict = signal_to_noise(cluster)
        pvalues = {}
        paths = {}
        print "Calculating pathway enrichment scores..."
        # extract all pathway genelists:
        paths, smallpaths = genematch.collect_kegg_pathways(minsize=10)
        paths.update(smallpaths)
        paths.update(genematch.collect_ipr_pathways(minsize=10))
        paths.update(genematch.collect_go_pathways(minsize=10) )
        # calculate enrichment score:
        for pathway in paths:
            pvalues[pathway] = irizarry_enrichment_z(t_dict, paths[pathway])

        # calculate q-values for the list of p-values:
        qvalues = p_to_q(pvalues.values(), display_on=not(args.display_off))
        for p,q in sorted(qvalues.iteritems()):
            if  0.95 < p <= 0.05:
                print "%-8.4f %-.4f" % (p, qvalues[p])

        out_h = open(filename[:-4] + ".GSEA_Irizarry_z.list", 'w')
        for pathway,pval in sorted(pvalues.iteritems(), key=itemgetter(1)):
            if pvalues[pathway] <= 0.05:
                out_h.write( "%-15s %-50s%s P=%7.4f q=%7.4f\n" % (
                    pathway[0], pathway[1][:50], '...' if len(pathway[1])>50 else '   ',
                    pval, qvalues[pval]))
        out_h.close()

    ## Principal Component Analysis:
    if args.pca:
        matrix_red = analyse_pca(cluster, data_table)

    if args.kill_PC1:
        cluster.data_matrix = matrix_red

    ## ANOVA analysis
    if args.filter_anova:
        a_dict = degs_anova(cluster, groups=groups)
        # collect q-values:
        q_dict = p_to_q(a_dict.values(), display_on=not(args.display_off))

        # output results to file:
        a_list = []
        out_h = open(filename[:-4] + ".ANOVA.list", 'w')
        out_h.write('Gene                p-value q-value\n')
        for gene in a_dict:
            if q_dict[a_dict[gene]] <= args.filter_anova:
                out_h.write( "%-12s %12.9f %.5f\n" % (gene, a_dict[gene], q_dict[a_dict[gene]] ))
                a_list.append(gene)
        out_h.close()
        print "Filtering matrix to %d genes with ANOVA P-value less than %.4f" % (len(a_list), args.filter_anova)
        cluster.filter_genes(a_list)
    elif args.anova:
        a_dict = degs_anova(cluster, groups=groups)
        q_dict = p_to_q(a_dict.values(), display_on=not(args.display_off))
        out_h = open(filename[:-4] + ".ANOVA.list", 'w')
        out_h.write('Gene                p-value q-value\n')
        for gene in a_dict:
            if q_dict[a_dict[gene]] <= args.anova:
                out_h.write( "%-12s %12.9f %.5f\n" % (gene, a_dict[gene], q_dict[a_dict[gene]]) )
        out_h.close()

    ## t-test analysis
    if args.ttest_thresh:
        t_dict = find_degs(cluster, groups=groups)
        q_dict = p_to_q([v[1] for v in t_dict.values()], display_on=not(args.display_off))
        # report output to file:
        t_list = []
        out_h = open(filename[:-4] + ".t_test.list", 'w')
        out_h.write('%-12s %-7s %-7s %s\n' % ('Gene','p-value', 'q-value', 'difference'))
        for gene in t_dict:
            if q_dict[t_dict[gene][1]] <= args.ttest_thresh:
                out_h.write( "%-12s %7.4f %.4f %.1f\n" % (gene, t_dict[gene][1],  q_dict[t_dict[gene][1]], t_dict[gene][2]))
                t_list.append(gene)
        out_h.close()
        print "Filtering matrix to %d genes with t-test P-value less than %.2f" % (len(t_list),args.ttest_thresh)
        cluster.filter_genes(t_list)
    elif args.t_test:
        t_dict = find_degs(cluster, groups=groups)
        q_dict = p_to_q([v[1] for v in t_dict.values()], display_on=not(args.display_off))
        out_h = open(filename[:-4] + ".t_test.list", 'w')
        out_h.write('%-12s %-7s %-7s %s\n' % ('Gene','p-value', 'q-value', 'difference'))
        for gene in t_dict:
            if q_dict[t_dict[gene][1]] <= args.t_test:
                out_h.write( "%-12s %7.4f %.4f %.2f\n" % (gene, t_dict[gene][1],  q_dict[t_dict[gene][1]], t_dict[gene][2]))
        out_h.close()

    ## report nearest neighbours:

    if args.neighbours:
        print "\nCalculating nearest neighbours..."
        # +1 to argument, since the first neighbour returned is the gene itself.
        neighbour_dict = find_nearest_neighbours(cluster, args.neighbours,
                                        numberofneighbours=args.num_neighbours + 1)

        out_h = open(filename[:-4] + ".nearest_neighbours.list", 'w')
        for neigh in neighbour_dict:
            t_neighbours, stat_str = neighbour_dict[neigh]
            out_h.write("## %s\n%s\n%s\n\n" % (neigh, stat_str, ",".join(t_neighbours[1:])) )
        out_h.close()

    ## bar-chart construction:
    if args.bar_charts:
        genelist = make_a_list(args.bar_charts, col_num=0)
        print "Constructing bar charts for %d genes: %s" % ( len(genelist), " ".join(genelist))
        bar_charts(cluster, genelist, groups=groups)

    ### re-format for publication purposes:
    ## transpose so genes are on y-axis:
    if args.display_reverse:
        cluster.invert_matrix()

    ## display average values for each group instead of individual values
    if args.show_averages is True:
        cluster.average_matrix(groups=groups)


    ################ perform hierarchical clustering ####################
    if cluster.genenumber > 1 and args.no_clustering is False:
        try:
            new_column_header, groups = heatmap(cluster, display=not(args.display_off), kegg=args.kegg, go=args.go_enrichment)
        except Exception as inst:
            print 'Error using %s ... trying euclidean instead\n%s' % (args.gene_metric, inst)

            cluster.gene_metric = 'euclidean'
            cluster.sample_metric = 'euclidean'
            cluster.refresh_headers()
            try:
                new_column_header, groups = heatmap(cluster, display=not(args.display_off), kegg=args.kegg, go=args.go_enrichment)
            except IOError:
                print 'Error with clustering encountered'
                new_column_header = ['']
                groups = ['']

    return cluster

def heatmap(cluster, display=True,kegg=False, go=False):
    ################# Perform the hierarchical clustering #################
    ###  heatmap function updated from original code hierarchical_clustering.py
    #    Copyright 2005-2012 J. David Gladstone Institutes, San Francisco California
    #    Author Nathan Salomonis - nsalomonis@gmail.com
    #
    #    Original message:
    #Permission is hereby granted, free of charge, to any person obtaining a copy
    #of this software and associated documentation files (the "Software"), to deal
    #in the Software without restriction, including without limitation the rights
    #to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    #copies of the Software, and to permit persons to whom the Software is furnished
    #to do so, subject to the following conditions:

    #THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    #PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    #HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
    #OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    #SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    verbalise("Y", "\nPerforming heirarchical clustering using %s for columns and %s for rows" % \
            (cluster.column_metric, cluster.row_metric))

    """
    This below code is based in large part on the protype methods:
    http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
    http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre

    cluster.data_matrix is an m by n ndarray, m observations, n genes
    """

    verbalise("G", "clustering %d genes and %d samples." % \
            (cluster.genenumber, cluster.samplesize ))

    ### Define the color gradient to use based on the provided name
    n = cluster.genenumber; m = cluster.samplesize

    ### Scale the max and min colors so that 0 is white/black
    vmax, vmin = cluster.getColorRange()
    verbalise("Vmax: %.2f\tVmin: %.2f" % (vmax, vmin))
    norm = mpl.colors.Normalize(vmin/2, vmax/2) ### adjust the max and min to scale these colors

    ### Scale the Matplotlib window size
    default_window_hight = 14.5 # 8.5
    default_window_width = 16
    fig = plt.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    color_bar_w = 0.015 ### Sufficient size to show

    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    #if cluster.row_method != None: w1 =
    [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.42,0.2,0.4] # [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix

    # axr, placement of row side colorbar
    [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y; axr_h = ax1_h
    width_between_axr_axm = 0.004

    # axc, placement of column side colorbar
    [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004

    # axm, placement of heatmap for the data matrix
    [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y; axm_h = ax1_h
    axm_w = axc_w

    # ax2, placement of dendrogram 2, on the top of the heatmap
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.09]

    # Compute and plot top dendrogram
    if cluster.column_method != None:
        start_time = time.time()
        d2 = dist.pdist(cluster.data_matrix.T)
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
        Y2 = sch.linkage(D2, method=cluster.column_method, metric=cluster.column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        # distances in linkage array can occasionally be negative due to floating point
        # calculation error. Therefore, clip third column inarray to remove negative
        # values (by making them zero):
        numpy.clip(Y2[:,2], 0, 100000, Y2[:,2])
        Z2 = sch.dendrogram(Y2)
        # assign to clusters based on distances within and between individuals in clusters:
        ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
        ax2.set_xticks([]) ### Hides ticks
        ax2.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        verbalise('Column clustering completed in %s seconds' % time_diff)
    else:
        ind2 = ['NA']*len(cluster.column_header) ### Used for exporting the flat cluster data

    # Compute and plot left dendrogram.
    if cluster.row_method != None:
        start_time = time.time()
        d1 = dist.pdist(cluster.data_matrix)
        D1 = dist.squareform(d1)  # full matrix
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
        Y1 = sch.linkage(D1, method=cluster.row_method, metric=cluster.row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z1 = sch.dendrogram(Y1, orientation='right')
        ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        ax1.set_xticks([]) ### Hides ticks
        ax1.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        verbalise('Row clustering completed in %s seconds' % time_diff)
    else:
        ind1 = ['NA']*len(cluster.row_header) ### Used for exporting the flat cluster data

    # Plot distance matrix.
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
    xt = cluster.data_matrix
    if cluster.column_method != None:

        ## create a list containing the order of the rearranged matrix:
        idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
        xt = xt[:,idx2]
        ind2 = [ ind2[i] for i in idx2] # similarly organises the cluster assignment list
        #ind2 = ind2[:,idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
    if cluster.row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        xt = xt[idx1,:]   # xt is transformed x
        ind1 = [ ind1[i] for i in idx1 ]
        #ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cluster.cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
    axm.set_xticks([]) ### Hides x-ticks
    axm.set_yticks([])

    # Add text for row data, rearrange column headers according to clustering:
    new_row_header=[]
    new_column_header=[]
    for i in range(cluster.data_matrix.shape[0]):
        if cluster.row_method != None:
            if len(cluster.row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(cluster.data_matrix.shape[1]-0.5, i, '  ' + cluster.row_header[idx1[i]])
            new_row_header.append(cluster.row_header[idx1[i]])
        else:
            if len(cluster.row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(cluster.data_matrix.shape[1]-0.5, i, '  ' + cluster.row_header[i]) ### When not clustering rows
            new_row_header.append(cluster.row_header[i])
    for i in range(cluster.data_matrix.shape[1]):
        if cluster.column_method != None:
            #axm.text(i, -0.9, ' ' + cluster.column_header[idx2[i]], rotation=270, verticalalignment="top") #  rotation could also be degrees
            new_column_header.append(cluster.column_header[idx2[i]])
        else: ### When not clustering columns
            #axm.text(i, -0.9, ' ' + cluster.column_header[i], rotation=270, verticalalignment="top")
            new_column_header.append(cluster.column_header[i])




    ## perform enrichment analyses, and modify gene header to show relevent go terms:
    if go:
        # create list of genes in each group:
        leafpairs = zip(ind2, new_column_header)
        genelistd = {}
        for group,geneid in leafpairs:
            try:
                genelistd[group].append(geneid)
            except KeyError:
                genelistd[group] = [geneid]

        verbalise("Y", "\nPerforming GO enrichment analysis for %d groups" % (len(genelistd)))
        out_h = open(cluster.exportPath[:-4] + ".GO_enrichment.list", 'w')

        out_h.write("%-4s %-11s %-7s %-7s %s\n" % ("Grp", "GO-term", "p-value", "q-value", "Definition"))
        go_monster = genematch.GO_maker()
        # perform hypergeometric test (one-sided fishers exact test) on each group:
        for group in genelistd:
            out_h.write('%s\n' % ("#" * 15))
            gops = genematch.go_enrichment(genelistd[group])
            qvalues = p_to_q([val[0] for val in gops.values()], display_on=False)
            for goterm in gops:
                if qvalues[gops[goterm][0]] <= go:
                    out_h.write( "%-4s %-11s %.5f %.5f %s\n" % (group, goterm, gops[goterm][0], qvalues[gops[goterm][0]], str(go_monster.define_go(goterm))) )
                    # modify gene header by appending go term to gene name
                    even_newer_header = new_column_header
                    new_column_header[:] = [ appendgo(geneid, goterm, go_monster) for geneid in new_column_header ]
                    """
                    for geneid in new_column_header:
                        if goterm in go_monster.findem(geneid.split()[0]):
                            geneid = " ".join([geneid, goterm, go_monster.define_go(goterm)[0]])
                            print geneid
                    """
        out_h.close()

    if kegg:
        leafpairs = zip(ind2, new_column_header)
        genelistd = {}
        for group,geneid in leafpairs:
            try:
                genelistd[group].append(geneid)
            except KeyError:
                genelistd[group] = [geneid]

        verbalise("Y", "\nPerforming KEGG pathway enrichment analysis for %d groups" % (len(genelistd)))
        out_h = open(cluster.exportPath[:-4] + ".KEGG_enrichment.list", 'w')

        out_h.write("Group KEGG pathway P-value\n")

        for group in genelistd:
            pathway_ps, gene_kos = genematch.kegg_pathway_enrichment(genelistd[group], pthresh=kegg)
            #gops = genematch.go_enrichment(genelistd[group])
            qvalues = p_to_q(pathway_ps.values(), display_on=False)
            for ko in pathway_ps:
                if qvalues[pathway_ps[ko]] <= kegg:
                    out_h.write( "# %-4s %-7s %.5f %s\n" % (group, goterm, pathway_ps[ko]) )

        even_newer_header = new_column_header
        new_column_header[:] = [ appendkegg(geneid, gene_kos) for geneid in new_column_header ]

        out_h.close()



    # Add text for column data:
    for i in range(cluster.data_matrix.shape[1]):
        if len(cluster.column_header)<200:
            axm.text(i - 0.3, -0.9, ' ' + new_column_header[i], rotation=270, verticalalignment="top")


    # Plot colside colors
    # axc --> axes for column side colorbar
    if cluster.column_method != None:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        cmap_c = mpl.colors.ListedColormap(['g', 'r', 'c', 'y', 'w', 'k', 'm'][:len(ind2)])
        dc = numpy.array(ind2, dtype=int)
        dc.shape = (1,len(ind2))
        im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
        axc.set_xticks([]) ### Hides ticks
        axc.set_yticks([])

    # Plot rowside colors
    # axr --> axes for row side colorbar
    if cluster.row_method != None:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        dr = numpy.array(ind1, dtype=int)
        dr.shape = (len(ind1),1)
        #print ind1, len(ind1)
        cmap_r = mpl.colors.ListedColormap(['g', 'r', 'c', 'y', 'w', 'k', 'm'][:len(ind2)])
        im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
        axr.set_xticks([]) ### Hides ticks
        axr.set_yticks([])

    # Plot color legend
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cluster.cmap, norm=norm, orientation='horizontal')
    axcb.set_title("colorkey")


    pdfname = cluster.exportPath[:-4] + '.pdf'
    cb.set_label("Differential Expression (log2 fold)")
    ind1 = ind1[::-1] # reverse order of flat cluster leaves to match samples in txt file.
    exportFlatClusterData(pdfname, new_row_header,new_column_header,xt,ind1,ind2)

    ### Render the graphic
    if len(cluster.row_header)>50 or len(cluster.column_header)>50 or max(len(term) for term in new_column_header) > 80:
        plt.rcParams['font.size'] = 7 #5
    else:
        plt.rcParams['font.size'] = 10 #8

    plt.savefig(pdfname)
    verbalise("B",'Exporting:',pdfname)
    filename = cluster.exportPath[:-4]+'.png'
    plt.savefig(pdfname, dpi=200) #,dpi=100
    if display:
        plt.show()

    # returns values to allow group-specific KEGG enrichment analysis
    return new_column_header, ind2

def find_nearest_neighbours(cluster, geneobj, numberofneighbours=10):
    "For a given gene, return the closest neighbours by distance"

    # because the distance calculation takes the longest time,
    # lists are processed within the function to make the pdist calculation only once
    genelist = make_a_list(geneobj)
    results_dict ={}

    # calculate distances:
    if cluster._genes_as_rows:
        dists    = dist.pdist(cluster.data_matrix)
    else:
        dists    = dist.pdist(cluster.data_matrix.T)
    sq_dists = dist.squareform(dists)

    for geneid in genelist:
        if geneid in cluster.gene_header:
            query    = cluster.gene_header.index(geneid)
            top10    = sorted(zip(sq_dists[query], cluster.gene_header),
                              key=itemgetter(0)
                              )[:numberofneighbours]

            # calculate the max distance for nearest neighbours cf w all distances
            neighdists = [ x[0] for x in top10 ]
            z =  (max(neighdists) - numpy.mean(sorted(sq_dists[query])[1:]))/float(numpy.std(sorted(sq_dists[query])[1:]))
            p_value = stats.norm.cdf(z)
            stat_str =  "Max(neigh): %.3f\tMean(data): %.3f\tstdev(data): %.4f\tz: %.3f\tP: %.5f" % (max(neighdists), numpy.mean(sq_dists[query]), float(numpy.std(sq_dists[query])), z, p_value)
            results_dict[geneid] = [ x[1] for x in top10 ], stat_str

    return results_dict

################# Export the flat cluster data #####################################

def exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """

    filename = string.replace(filename,'.pdf','.txt')
    export_text = open(filename,'w')
    column_header = string.join(['UID','row_clusters-flat']+new_column_header,'\t')+'\n' ### format column-names for export
    export_text.write(column_header)
    column_clusters = string.join(['column_clusters-flat','-']+ map(str, ind2),'\t')+'\n' ### format column-flat-clusters for export
    export_text.write(column_clusters)

    ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
    new_row_header = new_row_header[::-1]
    xt = xt[::-1]

    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_text.write(string.join([new_row_header[i],str(ind1[i])]+map(str, row),'\t')+'\n')
        i+=1
    export_text.close()

    ### Transpose text file for easier reading!
    oldfile_h = open(filename, 'rb')

    elements = [ line.split() for line in oldfile_h ]
    oldfile_h.close()

    biglist = []
    for splitline in elements:
        #print len(splitline)
        #print splitline
        biglist.append(splitline)
    newarray = numpy.array(biglist)
    #print numpy.shape(newarray)
    t_array = newarray.transpose()
    #print numpy.shape(t_array)
    #print newarray[:,0]

    newfile_h = open(filename[:-4] + "_transposed.txt" , 'w')
    for row in t_array:
        #print "The row is currently: %r" % row
        newfile_h.write("\t".join(row) + "\n")
    newfile_h.close()


    ### Export as CDT file
    filename = string.replace(filename,'.txt','.cdt')
    export_cdt = open(filename,'w')
    column_header = string.join(['UNIQID','NAME','GWEIGHT']+new_column_header,'\t')+'\n' ### format column-names for export
    export_cdt.write(column_header)
    eweight = string.join(['EWEIGHT','','']+ ['1']*len(new_column_header),'\t')+'\n' ### format column-flat-clusters for export
    export_cdt.write(eweight)

    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_cdt.write(string.join([new_row_header[i]]*2+['1']+map(str, row),'\t')+'\n')
        i+=1
    export_cdt.close()

################# Data construction or import methods ##############################

def create_table(build_list):
    """Takes a comma-delimited list of files and uses them to build the fpkm table for
    analysis"""
    # extract list of files:
    all_files = build_list.split(',')

    output_file = os.getcwd() + "/fpkm_table.tbl"
    verbalise("B", "saving built table to file ", output_file)

    awk_cmd = """awk '!($10~/FPKM/){\
    gene_sample[$1,FILENAME]=$9;\
    samples[FILENAME]=1;genes[$1]=1}\
    END{printf "%s\t", "genes";\
    for(g in genes){printf "%s\t", g};print "";\
    for(s in samples){printf "%s\t",s;\
    for(g in genes){printf "%s\t", gene_sample[g,s]};\
    print ""}}' """

    full_cmd = awk_cmd + " ".join(all_files) + ">" + str(output_file)

    # build table:
    os.system(full_cmd)

    # shorten filenames within table for easier reading!
    #output

    return output_file

################# Data analyis methods #############################################

def analyse_pca(cluster, three_dim=True):
    """
    performs principal component analysis on the matrix and saves output as png.
    No additional normalisation occurs at this step - all normalisation has been done
    on the matrix prior to PCA.
    """
    # create data array and name array:
    A = cluster.data_matrix
    names = cluster.row_header

    # assign colours to samples:
    colorconvert = {'F':'go', 'S':'co', 1:'ro', 2:'go', 3:'ko', 4:'bo', 5:'co', 6:'mo', 7:'yo', 8:'r<', 9:'g<', 10:'k<', 11:'b<', 12:'c<', 13:'m<', 14:'y<', 15:'rs', 16:'gs', 17:'ks', 18:'bs', 19:'cs', 20:'ms', 21:'ys' }
    colourlist = []
    for name in names:
        phase = re.search("(F|S)", name)
        if phase is not None:
            #print phase.groups()[0]
            colourlist.append(colorconvert[phase.groups()[0]])
        else:
            colourlist.append('ko')
    #print names, "\n", colourlist

    ############# PCA using numpy SVD decomposition ##################################
    print "#" * 30
    print "SVA analysis"
    U, s, Vt = numpy.linalg.svd(A, full_matrices=True)
    V = Vt.T

    # sort the PCs by descending order of the singular values (i.e. by the
    # proportion of total variance they explain)
    ind = numpy.argsort(s)[::-1]
    U = U[:, ind]
    s = s[ind]
    V = V[:, ind]
    S = numpy.diag(s)

    sumval = sum([ i ** 2 for i in s ])

    # if we use all of the PCs we can reconstruct the noisy signal perfectly

    # Mhat = numpy.dot(U, numpy.dot(S, V.T))
    # if we use only the first 2 PCs the reconstruction is less accurate
    # Mhat2 = numpy.dot(U[:, :2], numpy.dot(S[:2, :2], V[:,:2].T))

    # To remove the variance of the 1st PC, which is primarily associated with experimenter:
    matrix_reduced = numpy.dot(U[:,1:], numpy.dot(S[1:,1:], V[:,1:].T))
    #for checking decomposition is occurring properly:
    #print numpy.shape(U)
    #print numpy.shape(S)
    #print numpy.shape(Vt)
    #print numpy.shape(matrix_reduced)

    #print "#" * 30
    #print "SVD eigenvectors/loadings:"
    #print header[:var_num] , "\n"
    #print U        # need to work out appropriate way to calculate loadings!
    #print "#" * 30
    #print "checking distance of loadings (eigen vectors)"
    #for col in loadings[:,:]:
    #    print col
    #    print numpy.sqrt(sum([ a ** 2 for a in col ]))

    print "PCA explained variance:"
    print [ (z ** 2 / sumval) for z in s ]

    # * if M is considered to be an (observations, features) matrix, the PCs
    #   themselves would correspond to the rows of S^(1/2)*V.T. if M is
    #   (features, observations) then the PCs would be the columns of
    #   U*S^(1/2).

    #q_scores = numpy.dot(numpy.sqrt(S), V.T)
    q_scores = numpy.dot(U, numpy.sqrt(S))

    pp = PdfPages(cluster.exportPath[0:-4] + '.PCA.pdf')
    if three_dim:   # plot a three dimensional graph:
        fig = plt.figure(1)
        ax = fig.add_subplot(111, projection='3d')
        for idx in range(len(colourlist)):
            xs = q_scores[idx,0]
            ys = q_scores[idx,1]
            zs = q_scores[idx,2]
            name = re.search('[FS][LP][0-9]+',names[idx]).group(0)
            ax.scatter(xs, ys, zs, c=colourlist[idx][0], marker='o')
            ax.text(xs, ys, zs, name)

        ax.set_xlabel("PC1 (%.2f%%)" % (100.0 * (s[0]**2)/sumval))
        ax.set_ylabel("PC2 (%.2f%%)" % (100.0 * (s[1]**2)/sumval))
        ax.set_zlabel("PC3 (%.2f%%)" % (100.0 * (s[2]**2)/sumval))

        plt.savefig(pp, format='pdf')
        plt.show()
    else:   # plot two 2D graphs instead:
        for idx in range(len(colourlist)):
            fig = plt.figure(1)

            sub1 = fig.add_subplot(2,1,1)
            sub1.plot(q_scores[idx,0], q_scores[idx,1], colourlist[idx])
            plt.xlabel( "PC1 (%.2f%%)" % (100.0 * (s[0]**2)/sumval) )
            plt.ylabel( "PC2 (%.2f%%)" % (100.0 * (s[1]**2)/sumval) )
            sub1.annotate( names[idx], xy=(q_scores[idx,0], q_scores[idx,1]),xytext=(-15,10), xycoords='data', textcoords='offset points' )

            sub2 = fig.add_subplot(2,1,2)
            sub2.plot(q_scores[idx,0], q_scores[idx,2], colourlist[idx])
            plt.xlabel( "PC1 (%.2f%%)" % (100.0 * (s[0]**2)/sumval) )
            plt.ylabel( "PC3 (%.2f%%)" % (100.0 * (s[2]**2)/sumval) )
            sub2.annotate( names[idx], xy=(q_scores[idx,0],q_scores[idx,2]),xytext=(-15,10), xycoords='data', textcoords='offset points' )

        plt.savefig(pp, format='pdf')
        plt.show()

        plt.close()
    return matrix_reduced

def expression_dist(cluster, max_x=500, min_x=0, hist=False):
    "Creates histogram of gene expression"

    count = 0

    # create box and whisker plot of each sample:
    fig = plt.figure()
    subp = fig.add_subplot(111)
    something = subp.boxplot(cluster.data_matrix.transpose())

    # create labels for each sample:
    rng = [ i + 1 for i in range(len(cluster.row_header)) ]
    plt.xticks(rng , cluster.row_header)

    plt.show()

    if hist:
        for row in cluster.data_matrix:
            fig = plt.figure()
            subf = fig.add_subplot(111, title=cluster.row_header[count])
            n, bins, patches = subf.hist(row, 50, range=(min_x,max_x), histtype='stepfilled')
            count += 1
        plt.setp(labels, rotation=90)
        plt.show()

def find_degs(cluster, groups=["_FL","_SP"]):
    "finds DEGs using t-test and returns dictionary { Gene:t-value, P-value }"

    limits = [0] + cluster.reorder_matrix(groups=groups)
    t_dict = {}
    if len(limits) > 3:
        verbalise("R", "%d groups were specified. Only the first two will be compared in this t-test. Consider using ANOVA instead" % (len(limits)-1))
    verbalise("Y", "Performing t-test")

    for g in range(cluster.genenumber):
        exvalues = [cluster.data_matrix[limits[x]:limits[x+1],g] for x in range(len(limits) - 1)]
        difference = numpy.mean(exvalues[0]) - numpy.mean(exvalues[1])
        t_val, p_val = stats.ttest_ind(exvalues[0],exvalues[1])
        t_dict[cluster.gene_header[g]] = float(t_val), p_val, difference

    return t_dict

def degs_anova(cluster, groups=["SP", "SL06", "SL12", "SL24","SL48", "SL96", "FL", "FP06", "FP12", "FP24","FP48", "FP96" ], onegene=False, verbose=True):
    "finds DEGs using ANOVA and returns dictionary { Gene:P-value }"
    ### change to return {Gene: F-stat, Pval }
    # will only work if genes are columns in matrix
    revert = False
    if cluster._genes_as_rows:
        cluster.invert_matrix()
        revert = True

    A_dict = {}

    # set the positions at which group membership changes:
    limits = [0] + cluster.reorder_matrix(groups=groups)
    intervals = zip(limits[:-1],limits[1:])


    if onegene:
        try:
            g = cluster.gene_header.index(onegene)
            expn_data = [ (cluster.data_matrix[ i:j,g ]) for i,j in intervals]
            f_val, p_val = stats.f_oneway(*expn_data)
            A_dict[cluster.gene_header[g]] = p_val
        except:
            if verbose:
                print "Gene %s not found!" % (onegene)
    else:
        verbalise("Y", "Performing ANOVA for %d genes" % (cluster.genenumber))
        #print "limits for anova:",limits
        for g in range(cluster.genenumber):
            expn_data = [ (cluster.data_matrix[ i:j,g ]) for i,j in intervals]
            f_val, p_val = stats.f_oneway(*expn_data)
            A_dict[cluster.gene_header[g]] = p_val


    # if matrix was inverted for gene removal, restore to its previous orientation:
    if revert:
        cluster.invert_matrix()

    return A_dict

def signal_to_noise(cluster, groups=["_SP", "_FL"], verbose=False):
    """
    SNR = (u0 - u1) / (s0 + s1)

    snr_dict = {gene: signal to noise ratio}
    """

    # will only work if genes are columns in matrix
    revert = False
    if cluster._genes_as_rows:
        cluster.invert_matrix()
        revert = True

    snr_dict = {}
    limits = [0] + cluster.reorder_matrix(groups=groups, verbose=verbose)

    for g in range(cluster.genenumber):
        exvalues = [cluster.data_matrix[limits[x]:limits[x+1],g] for x in range(len(limits) - 1)]
        exmeans  = [numpy.mean(cluster.data_matrix[limits[x]:limits[x+1],g]) for x in range(len(limits) - 1)]
        exstdev  = [numpy.std(cluster.data_matrix[limits[x]:limits[x+1],g]) for x in range(len(limits) - 1)]

        if exstdev[0] < 0.2 * exmeans[0]:
            exstdev[0] = 0.2 * exmeans[0]
        if exstdev[1] < 0.2 * exmeans[1]:
            exstdev[1] = 0.2 * exmeans[1]

        snr_dict[cluster.gene_header[g]] =  (exmeans[0] - exmeans[1]) / sum(exstdev)

    # if matrix was inverted for gene removal, restore to its previous orientation:
    if revert:
        cluster.invert_matrix()

    return snr_dict

def enrichment_score( snr_dict, pathway_list, rho=1, pathway='unknown', display_on=False):
    """
    calculate the enrichment score for a pathway. Genes in the pathway are specified in
    pathway_list. The correlation of genotype to phenotype is the signal to noise ratio,
    but any other correlation can be used also. Rho is the scaling factor. According to
    Subramanian et al, if rho is 0, then the enrichment score scales to a classic
    Komogorov-Smirnov statistic.

    pathway variable is not used in the function, merely returned so that downstream
    functions know which pathway the pathway list came from
    """

    sorted_genes = ( (gene, value) for gene, value in sorted(snr_dict.iteritems(), key=itemgetter(1), reverse=True) )
    ES = [0]
    genespot = []
    correlation = []
    count = 0
    Nr = sum( abs(snr_dict[g])**rho for g in pathway_list if g in snr_dict)
    pathway_size = sum([1 for gene in pathway_list if gene in snr_dict])
    for gene, value in sorted_genes:
        if gene in pathway_list:
            ES.append(ES[-1] + (abs(snr_dict[gene])**rho)/Nr )
            genespot.append(count)
        else:
            ES.append( ES[-1] - 1. / (len(snr_dict) - pathway_size) )
        correlation.append(value)
        count += 1

    if max(ES[1:]) > abs(min(ES[1:])):
        maxdev = max(ES[1:])
    else:
        maxdev = min(ES[1:])


    if  display_on:
        plt.figure()
        ax1 = plt.subplot2grid((3,3),(0,0), rowspan=2, colspan=3)
        ax1.plot(range(len(ES) - 1), ES[1:], 'r', genespot, [0 for x in genespot], 'g^')
        ax1.set_ylabel('Enrichment Score (ES)')
        ax2 = plt.subplot2grid((3,3),(2,0), colspan=3)
        ax2.fill_between(range(len(correlation)), 0, correlation, facecolor='cyan')
        ax2.set_xlabel('Gene (ordered by signal-to-noise ratio)')
        ax2.set_ylabel('Correlation')
        ax1.axes.get_xaxis().set_ticks([])
        ax2.axes.get_xaxis().set_ticks([])
        plt.suptitle("GSEA for %-15s %-50s%s" %
            (pathway[0], pathway[1][:50], '...' if len(pathway[1])>50 else '   '))
        plt.show()
    return  maxdev

def nominal_p(score, distribution):
    "Given the empirical score, and permuted distribution, calculate the P value of score"

    if score >= 0:
        beatenby = len([x for x in distribution if x >= score])
    else:
        beatenby = len([x for x in distribution if x <= score])

    return 1.0 * beatenby / len(distribution)

def normalise_es(score, distribution):
    """
    "The global statistics Sk depend on size of gene set and therefore are not
    identically distributed. GSEA addresses this issue by normalizing Sk values to
    factor out the intrinsic gene set size dependence. The relevant normalization is a
    change or scale using the expected value of the positive (negative) null distribution
    statistic induced by sample permutation."
                                    - Tamayo et al. (2012) Stat Methods Med Res.

    This function normalises the observed enrichment score and the permuted es
    distribution by dividing by the mean of the distribution (conditional on the
    direction of the score, +ve or -ve.

    """

    # find average of all permutations from +ve distribution:
    exp_es_plus  = numpy.mean([x for x in distribution if  x >= 0])
    # find average of all permutations from -ve distribution:
    exp_es_minus = numpy.mean([x for x in distribution if x < 0])

    # normalise enrichment score:
    if score >= 0:
        nes = 1.0 * score / exp_es_plus
    else:
        nes = -1.0 * score / exp_es_minus

    # normalise permutation distribution:
    for res, posn in zip(distribution[:],range(len(distribution))):
        if res >= 0:
            distribution[posn] = res / exp_es_plus
        else:
            distribution[posn] = -1.0 * res / exp_es_minus
    return nes, distribution

def gene_set_enrichment(cluster, permutations=1000, processes=3, display_on=True, showbest=5):
    verbalise("C", "\n## Performing gene set enrichment analysis ##")
    snr_dict = signal_to_noise(cluster)

    verbalise("Y", "Calculating pathway enrichment scores...")
    opaths, smallpaths = genematch.collect_kegg_pathways(minsize=10)
    opaths.update(smallpaths)
    opaths.update(genematch.collect_ipr_pathways(minsize=10))
    opaths.update(genematch.collect_go_pathways(minsize=10) )
    # some pathways will have less genes than reported, as the gene expression for some
    # will mean they were filtered out. So the final pathway set must be checked for
    # genes in which there are insufficient number from this filtering:
    paths = {}
    es = {}
    for pathway in opaths:
        if sum([1 for gene in opaths[pathway] if gene in snr_dict]) >= 10:
            paths[pathway] = opaths[pathway]
            es[pathway] = enrichment_score(snr_dict, paths[pathway], pathway=pathway, display_on=False)
    #print paths.keys()[:5]
    # computing significance:
    verbalise("Y", "Permuting pathways %d times" % permutations)

    t = time.time()

    p = {}
    q = {}
    proc = {}
    rand_es_dict = {}
    for n in range(processes):
        p[n] = Permutable(cluster)
        q[n] = Queue()
        proc[n] = Process(target=permute_data, args=(p[n], q[n], (permutations + processes - 1) / processes, paths))
    # change last process to run with progress bar:
    proc[n] = Process(target=permute_data_progress, args=(p[n], q[n], (permutations + processes - 1) / processes, paths))

    for n in range(processes):
        proc[n].start()

    for n in range(processes):
        rand_es_dict[n] = q[n].get()
        proc[n].join()

        if (n+1) in rand_es_dict:
            if rand_es_dict[n] == rand_es_dict[n + 1]:
                verbalise("R", "*** WARNING: PERMUTATIONS ARE NOT RANDOM! ***")
    calc_time = time.time() - t
    m, s = divmod(calc_time, 60)
    h, m = divmod(m, 60)
    verbalise("Total calculation time: %d:%02d:%02d" % (h, m, s))

    rand_es = {}
    for path in rand_es_dict[0]:
        rand_es[path] = rand_es_dict[0][path]
        for n in range(processes)[1:]:
            rand_es[path] += rand_es_dict[n][path]

    # normalise scores for multiple testing comparison:
    nes      = {}
    rand_nes = {}
    for pathway in paths:
        nes[pathway], rand_nes[pathway] = normalise_es(es[pathway], rand_es[pathway])

    # calculate FDR q value for each test:
    qvalues = {}
    all_rand_nes = list(chain.from_iterable(rand_nes.values()))


    pp = PdfPages(cluster.exportPath[0:-4] + '.GSEA_permutations.pdf')
    bins = numpy.linspace(-2, 2, 100)
    cleaned_es = [x for x in nes.values() if x != 'nan']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    n, bins, patches = ax1.hist(all_rand_nes,
            weights=numpy.zeros_like(all_rand_nes) + 100. / len(all_rand_nes),
            bins=bins, facecolor='green', alpha=0.5)
    n, bins, patches = ax1.hist(cleaned_es,
            weights=numpy.zeros_like(cleaned_es) + 100. / len(cleaned_es),
            bins=bins, facecolor='red', alpha=0.5)
    plt.title('distribution of normalised enrichment scores')
    ax1.set_xlabel('ES')
    ax1.set_ylabel('#')
    plt.savefig(pp, format='pdf')
    if display_on:
        plt.show()

    # display ES for best candidates:
    if display_on and showbest > 0:
        for pathway, escore in [(pathway, escore) for pathway, escore in sorted(
                nes.iteritems(), key=itemgetter(1), reverse=True)][:showbest]:
            es = enrichment_score(snr_dict, paths[pathway], pathway=pathway, display_on=True)


    for pathway in paths:
        qvalues[pathway] = gsea_q(all_rand_nes, nes[pathway], nes)

    # calculate nominal normalised p values:
    pvalues = {}
    for pathway in paths:
        pvalues[pathway] = nominal_p(nes[pathway], rand_nes[pathway])

    return nes, pvalues, qvalues

def permute_data(cluster, queue, permutations, paths ):
    rand_es = {}
    for pathway in paths:
        rand_es[pathway] = []
    numpy.random.seed(os.getpid())
    for i in range(permutations):
        cluster.randomise_samples()
        snr_dict = signal_to_noise(cluster)
        for path in paths:
            rand_es[path].append(enrichment_score(snr_dict, paths[path]))
    queue.put(rand_es)

def randomise_and_analyse(cluster, processes=3):
    verbalise("Y", "Permuting pathways on %d cores" % processes)

    t = time.time()

    p = {}
    q = {}
    proc = {}
    rand_es_dict = {}
    for n in range(processes):
        p[n] = Permutable(cluster)
        q[n] = Queue()
        proc[n] = Process(target=permute_data, args=(p[n], q[n], (permutations + processes - 1) / processes, paths))
    # change last process to run with progress bar:
    proc[n] = Process(target=permute_data_progress, args=(p[n], q[n], (permutations + processes - 1) / processes, paths))

    for n in range(processes):
        proc[n].start()

    for n in range(processes):
        rand_es_dict[n] = q[n].get()
        proc[n].join()

        if (n+1) in rand_es_dict:
            if rand_es_dict[n] == rand_es_dict[n + 1]:
                verbalise("R", "*** WARNING: PERMUTATIONS ARE NOT RANDOM! ***")
    calc_time = time.time() - t
    m, s = divmod(calc_time, 60)
    h, m = divmod(m, 60)
    verbalise("Total calculation time: %d:%02d:%02d" % (h, m, s))

    rand_es = {}
    for path in rand_es_dict[0]:
        rand_es[path] = rand_es_dict[0][path]
        for n in range(processes)[1:]:
            rand_es[path] += rand_es_dict[n][path]

def permute_data_progress(cluster, queue, permutations, paths):

    # perform permutations:
    rand_es = {}
    for pathway in paths:
        rand_es[pathway] = []
    numpy.random.seed(os.getpid())
    tstart = time.time()
    pbcount = 0
    for i in range(permutations):
        if pbcount > 0:
            time_elapsed = time.time() - tstart
            time_per_cycle = 1. *  time_elapsed / pbcount
            time_remaining = (permutations - pbcount) * time_per_cycle
            m, s = divmod(time_remaining, 60)
            h, m = divmod(m, 60)
            if pbcount > 0:
                print '\r>> %d of %d cycles complete. ETA to completion: %d:%02d:%02d                      ' % (pbcount, permutations, h, m, s),
            sys.stdout.flush()
        else:
            time_remaining = permutations * 40
            m, s = divmod(time_remaining, 60)
            h, m = divmod(m, 60)
            print '>> 0 of %d cycles complete. Completion in approx: %d:%02d:%02d                      ' % (permutations, h, m, s),
            sys.stdout.flush()
        cluster.randomise_samples()
        snr_dict = signal_to_noise(cluster)
        for path in paths:
            rand_es[path].append(enrichment_score(snr_dict, paths[path]))
        pbcount += 1




    queue.put(rand_es)

def irizarry_enrichment(t_dict, pathway_list):
    """ a scale-sensitive test using the chi-square distribution to test for enrichment of
    a gene set. This method is (ostensibly) more sensitive than the GSEA K-S test
    (see Irizarry et al (2009) Stat Methods Med Res. 18(6) 565-575

    see also http://www.biostat.jhsph.edu/~ririzarr/688/enrichment.pdf for a helpful
    elaboration on the method.
    """


    t_hat = sum(t_dict[g][0] for g in pathway_list if g in t_dict) / len(pathway_list)
    # this score consistently gives only -ve values (chi square values should be +ve)
    # avoid using until fixed and/or confirmed.
    es_score = sum( [(t_dict[g][0] - t_hat)**2 - len(pathway_list) + 1 for g in pathway_list if g in t_dict] )     / numpy.sqrt(2 * (len(pathway_list) - 1))
    pvalue = 1 - stats.chi2.cdf(es_score, df=len(pathway_list) - 1)
    return es_score

def irizarry_enrichment_z(t_dict, pathway_list):
    """ a scale-sensitive test using the chi-square distribution to test for enrichment of
    a gene set. This method is (ostensibly) more sensitive than the GSEA K-S test
    (see Irizarry et al (2009) Stat Methods Med Res. 18(6) 565-575

    see also http://www.biostat.jhsph.edu/~ririzarr/688/enrichment.pdf for a helpful
    elaboration on the method.
    """


    t_hat = sum(t_dict[g][0] for g in pathway_list if g in t_dict) / len(pathway_list)
    es_score = numpy.sqrt(len(pathway_list)) * t_hat
    pvalue = 1 - stats.t.cdf(es_score, len(pathway_list) - 1)
    return pvalue

def bar_charts(cluster, genelist, groups=["SP", "SL06", "SL12", "SL24","SL48", "SL96", "FL", "FP06", "FP12", "FP24","FP48", "FP96" ], postfix=''):
    """creates bar plots of all genes in genelist, showing mean and error for each
    timepoint category specified in groups"""

    limits = cluster.reorder_matrix(groups)
    pp = PdfPages(cluster.exportPath[0:-4] + postfix + '.bar_plots.pdf')

    # get kegg pathways and NCBI values for each gene:
    ko_dict = genematch.cbir_to_pathway(genelist.keys())   # ko_dict = {gene:str(pathway)}
    go_monster = genematch.GO_maker()
    ncbi_terms = genematch.cbir_ncbi(genelist)

    for gene in genelist:
        # get gene details for later use:
        ignore, kotermdic = genematch.cbir_to_kegg([gene],reversedic=True)

        anova = degs_anova(cluster, onegene=gene, groups=groups)

        try:
            koterm = kotermdic[gene]
        except KeyError:
            koterm = 'no KO'

        genegos = go_monster.findem(gene)
        godesc = "".join([ "%s %s %s\n" % (g, genegos[g][1], genegos[g][0]) for g in genegos ])

        # calculate mean/SEM...
        if gene in cluster.column_header:
            pos = cluster.column_header.index(gene)
        else:
            continue
        gm = [groups[0]] * (limits[0])    # matrix of group names for Tukey's post hoc
        v = [numpy.average(cluster.data_matrix[:limits[0],pos])]   # averages
        se = [numpy.std(cluster.data_matrix[:limits[0],pos])/numpy.sqrt(limits[0]+1)]   #SEM
        for i in range(len(groups)-1):
            gm += [groups[i+1]] * (limits[i+1]-limits[i])
            v.append(numpy.average(cluster.data_matrix[limits[i]:limits[i + 1],pos]))
            se.append(numpy.std(cluster.data_matrix[limits[i]:limits[i + 1],pos])/numpy.sqrt(limits[i+1]-limits[i]+1))

        # calculate tukey's post-hoc values and plot:
        tfig, taxes = plt.subplots()

        try:
            posthoc = pairwise_tukeyhsd(cluster.data_matrix[:,pos],gm)
        except Exception as inst:
            verbalise("R", "Tukey calculation error - check that you have >1 value for each category.")
            print inst
            continue
        phimg = posthoc.plot_simultaneous(comparison_name='SP', \
            ax=taxes, ylabel='Groups', xlabel='Normalised Expression', \
            labelorder = ["SP", "SL06", "SL12", "SL24","SL48", "SL96", \
                          "FL", "FP06", "FP12", "FP24","FP48", "FP96" ])

        # plot_simultaneous does not correctly report the y-axis labels. So to fix:
        taxes.set_xticks(numpy.arange(13.0)*1)  # increase to gain all labels
        plt.tight_layout()                      # resets axes
        xlabels = taxes.get_xticklabels()       # gets values I need

        labelist = [xtick.get_text() for xtick in xlabels]  # creates an ordered list of labels
        labelist.pop(0)                         # removes first element (blank label)
        taxes.set_xticks(numpy.arange(12.0)*1)  # now create the right number of ticks
        taxes.set_xticklabels(labelist)         # reset with new names
        title_string = "%s %s(ANOVA P-value %.8f)\n%s\n KEGG ortholog %s:\n%s\n%s"
        taxes.set_title(title_string % (os.path.basename(cluster.exportPath[:-4]), gene, anova[gene], ncbi_terms[gene], koterm, ko_dict[gene], godesc), fontsize=12 )

        plt.tight_layout()
        plt.savefig(pp, format='pdf')
        #plt.show(phimg)
        plt.close()
        # print summary to file:
        tukeys_h = open(cluster.exportPath[:-4] + '.tukeys.txt','a')
        tukeys_h.write('Gene '  + str(gene) + ':\n')
        tukeys_h.write(str(posthoc) + '\n\n')
        tukeys_h.close()

        """
        # create box plot of expression values:
        ind = numpy.arange(len(groups))    # x-coords for bars
        width = 0.35                        # box width

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, v, width, color='r', yerr=se)

        # add details:
        ax.set_ylabel('Normalised Expression')
        ax.set_title('Gene Expression for %s (%s):\n %s\n%s' % (str(gene), koterm, ko_dict[gene], godesc), fontsize=12 )
        ax.set_xticks(ind+width)
        ax.set_xticklabels(groups)

        plt.tight_layout()
        plt.savefig(pp, format='pdf')
        plt.show()
        """
    pp.close()

def expression_peaks(cluster, magnitude, group1 = [ "SP", "SL06", "SL12", "SL24","SL48", "SL96" ], group2 = [ "FL", "FP06", "FP12", "FP24","FP48", "FP96" ]):
    """ Finds timepoint with max expression, then Checks adjacent timepoints to see if
    values lie with 2-50% of max expression value. Returns all genes that fit this
    pattern.
    Assumes first timepoint in each group is a control

    Instead of looking for fold change, this function now looks at the difference in
    peak size - this is a more appropriate metric for log-transformed data.
    """
    if cluster.averaged == False:
        cluster.average_matrix(group1 + group2)
    verbalise("G", cluster.sample_header)
    peaklist = {}

    for gene in range(cluster.genenumber):
        # for group 1:
        datalist = list(cluster.data_matrix[:,gene])
        maxexpression = max(datalist[:len(group1)])
        maxposn = datalist.index(maxexpression)

        # check fold change is sufficient:
        if maxexpression >= magnitude + datalist[0]:
            # check adjacent peaks are not too big:
            # difference of 5.64 corresponds to 2% of the untransformed fpkm value
            # difference of 1.00 corresponds to 50% of the untransformed fpkm value
            if maxposn == len(group1) - 1:
                if (maxexpression - 5.64 < datalist[maxposn - 1] < maxexpression - 1):
                    peaklist[cluster.gene_header[gene]] = group1[maxposn]

            elif (maxexpression * 0.02 < datalist[maxposn - 1] < maxexpression * 0.5) and \
                    (maxexpression * 0.02 < datalist[maxposn + 1] < maxexpression * 0.5):

                peaklist[cluster.gene_header[gene]] = group1[maxposn]

        # for group 2:
        maxexpression = max(datalist[len(group1):])
        maxposn = datalist.index(maxexpression)

        # check fold change is sufficient for reciprocal swap:
        if maxexpression >= magnitude * datalist[len(group1)]:
            # check adjacent peaks are not too big:
            try:
                if maxposn == len(group1+group2) - 1:
                    if (maxexpression * 0.02 < datalist[maxposn - 1] < maxexpression * 0.5):
                        peaklist[cluster.gene_header[gene]] = (group1 + group2)[maxposn]

                elif (maxexpression * 0.02 < datalist[maxposn - 1] < maxexpression * 0.5) and \
                           (maxexpression * 0.02 < datalist[maxposn + 1] < maxexpression * 0.5):

                    peaklist[cluster.gene_header[gene]] = (group1 + group2)[maxposn]
            except IndexError as inst:
                verbalise("R", inst)
                verbalise("R", datalist)
                verbalise("R", "Max is %.3f at position %d" % (maxexpression, maxposn))

    verbalise("G", len(peaklist),  "significant peaks found.")
    return peaklist

def p_to_q(pvalues, display_on=False, cut1s=False, set_pi_hat=False):
    """
    Given the list of pvalues, convert to pFDR q-values.
    According to Storey and Tibshirani (2003) PNAS 100(16) : 9440

    returns: { p-value:q-value }

    """
    # because fisher's exact test gives highly skewed pvalue dists (with P of 1)
    # it may be necessary to remove the 1s before analysing
    if cut1s:
        pvalues = [ps for ps in pvalues if ps < 1]

    # order p-values:
    pvalues.sort()

    # estimate pi0:
    # evaluate pi0 across the range of lambda:
    lamrange = numpy.arange(0,0.95,0.01)
    #pbeaters = [ sum( p > lam for p in pvalues) for lam in lamrange ]
    #denominator = [ (len(pvalues) * (1 - lam)) for lam in lamrange ]
    pi0_lam = [ (sum( p > lam for p in pvalues) / (len(pvalues) * (1 - lam))) for lam in lamrange ]
    #pi0_hardway = []

    #for i in range(len(pbeaters)):
    #    pi0_hardway += [ pbeaters[i] / denominator[i] ]
    #if pi0_lam != pi0_hardway:
    #    print "\n\n\npi0_lam is not the same as pi0_hardway!\n\n"
    #print "pi0_hardway length:", len(pi0_hardway)
    #print "p_values size:", len(pvalues)
    # fit cubic spline to data, then calculate value of pi0 for lambda = 1:
    tck = interpolate.splrep(lamrange, pi0_lam, s=3)
    splinecurve = interpolate.splev(numpy.arange(0,1.0,0.01), tck, der=0)
    pi0_hat = interpolate.splev(1, tck, der=0)
    tck_half = 0
    if pi0_hat > 1:
        tck_half = interpolate.splrep(lamrange[:85], pi0_lam[:85], s=3)
        spline_half = interpolate.splev(numpy.arange(0,1.0,0.01), tck_half, der=0)
        pi0_hat_half = interpolate.splev(1, tck_half, der=0)
        pi0_hat = pi0_hat_half
        verbalise("R", "pi0_hat > 1! Likely skewed P-value distribution. Converting to ", pi0_hat_half)
    if set_pi_hat:
        pi0_hat = set_pi_hat
    if display_on:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        try:
            n, bins, patches = ax1.hist(pvalues, bins=20, facecolor='green', label="P-values")
        except IndexError:
            ax1.plot(pvalues)
        plt.title('distribution of P-values')
        ax1.set_xlabel('lambda / P-value')
        ax1.set_ylabel('distribution #')
        plt.legend(loc=4)
        ax2 = ax1.twinx()
        ax2.plot(lamrange, pi0_lam, 'ro', numpy.arange(0,1.0,0.01), splinecurve, 'r', label='pi0_hat, s=3' )
        if tck_half != 0:
            ax2.plot(lamrange[:95], spline_half[:95], 'b', label='lambda < 0.85')
        ax2.set_ylabel('pi0_hat(lambda)')
        #ax1.plot(t, s1, 'b-')
        plt.legend(loc=1)
        plt.show()


    q_pm =  pi0_hat * pvalues[-1]    #  q(pm)
    # creates an ordered list of q(p(i)) values.
    q_pi_list = [q_pm] + [ (pi0_hat * len(pvalues)*pvalues[i])/i for i in range(len(pvalues)-1,1,-1)]
    # "The estimated q value for the ith most significant feature is q(p(i))"
    q_val = {}
    for i in range(len(pvalues)):
        q_val[pvalues[-1 * (i+1)]] = min(q_pi_list[:i+1])

    return q_val

def gsea_q(all_rand_nes, nes, nes_dict):
    """Calculate the q-value for a given normalised enrichment score, a given
    normalised permuted distribution and the normalised enrichment score distribution"""

    if nes >= 0:
        rand_gt_nes = float(len([ x for x in all_rand_nes if x >= nes if x >= 0]))
        num_rand    = float(len([ x for x in all_rand_nes if x >= 0 ]))
        alles_gt_nes= float(len([ x for x in nes_dict.values() if x >= nes if x >= 0 ]))
        num_alles   = float(len([ x for x in nes_dict.values() if x >= 0 ]))

        qvalue = (rand_gt_nes / num_rand) * (num_alles / alles_gt_nes)


    else:
        rand_lt_nes = float(len([ x for x in all_rand_nes if x <= nes if x <= 0]))
        num_rand    = float(len([ x for x in all_rand_nes if x <= 0 ]))
        alles_lt_nes= float(len([ x for x in nes_dict.values() if x <= nes if x <= 0 ]))
        num_alles   = float(len([ x for x in nes_dict.values() if x <= 0 ]))

        qvalue = (rand_lt_nes / num_rand) * (num_alles / alles_lt_nes)

    return qvalue

################# Miscellaneous methods #############################################

def make_a_list(geneobj, col_num=0):
    """given a path, list, dictionary or string, convert into a list of genes.
    col_num specifies the column from which to extract the gene list from."""

    # genefile can be given as a list of genes (say, from find_degs()... ), or as
    # a path to a file containing a list of genes.
    # The following builds a dictionary of genes from either input:
    if type(geneobj) is list:   # allows filtering from hicluster generated list of results.
        genelist = {}.fromkeys(geneobj,1)
    elif type(geneobj) is dict:
        genelist = geneobj # this of course will only work if the keys are the genes!
    elif type(geneobj) is str:   # assuming a filepath...
        if re.search("/",geneobj) is None:
            genelist = {}.fromkeys(geneobj.split(','),1)
        else:   # is a file path
            genefile_h = open(geneobj, 'rb')
            genelist = {}   # will be a dict of names { 'Cbir01255':1, 'CbirVgq':1, ... }
                            # used a dictionary to automatically remove any name duplications
            filegen = [ line.split() for line in genefile_h if len(line) > 0]

            genefile_h.close()

            for colset in filegen:
                try:
                    genelist[colset[col_num]]=1
                except IndexError:
                    verbalise("R", "Column %d not found in %s" % (col_num, str(colset)))
    else:
        genelist = {}

    return genelist

def appendgo(geneid, goterm, go_obj):
    if goterm in go_obj.findem(geneid.split()[0]):
        geneid = " ".join([geneid, goterm, go_obj.define_go(goterm)[0]])
    return  geneid

def appendkegg(geneid, ko_dictionary):
    if geneid in ko_dictionary:
        geneid = " ".join([geneid, ko_name[geneid]])
    return geneid

# #                                                                              # #
# #                                                                              # #
# ############################## END FUNCTIONS ################################### #
####################################################################################



if __name__ == '__main__':

    parser = define_arguments()

    args = parser.parse_args()

    run_arguments(args)

