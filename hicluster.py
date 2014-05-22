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

import matplotlib.pyplot as plt
import pylab
import matplotlib as mpl
import scipy
from scipy import stats
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy

import genematch

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

def heatmap(x, row_header, column_header, row_method,
            column_method, row_metric, column_metric,
            color_gradient, filename, display=True,
            kegg=False, go=False):

    print "\nPerforming hiearchical clustering using %s for columns and %s for rows" % (column_metric,row_metric)

    """
    This below code is based in large part on the protype methods:
    http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
    http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre

    x is an m by n ndarray, m observations, n genes
    """
    n = len(x[0]); m = len(x) # m  samples, n  genes
    print "Size of array: %d genes and %d samples.\nSize of header: %d genes\nSize of sample header: %d " % (n, m, len(column_header), len(row_header) )

    assert len(column_header) == n

    ### Define the color gradient to use based on the provided name
    n = len(x[0]); m = len(x)
    if color_gradient == 'red_white_blue':
        cmap=plt.cm.bwr
    if color_gradient == 'red_black_sky':
        cmap=RedBlackSkyBlue()
    if color_gradient == 'red_black_blue':
        cmap=RedBlackBlue()
    if color_gradient == 'red_black_green':
        cmap=RedBlackGreen()
    if color_gradient == 'yellow_black_blue':
        cmap=YellowBlackBlue()
    if color_gradient == 'seismic':
        cmap=plt.cm.seismic
    if color_gradient == 'green_white_purple':
        cmap=plt.cm.PiYG_r
    if color_gradient == 'coolwarm':
        cmap=plt.cm.coolwarm

    ### Scale the max and min colors so that 0 is white/black
    vmin=x.min()
    #print "line 77: vmin = %.5f\nm = %d\nn = %d" % (vmin, m, n)
    vmax=x.max()
    vmax = max([vmax,abs(vmin)])
    vmin = vmax*-1
    norm = mpl.colors.Normalize(vmin/2, vmax/2) ### adjust the max and min to scale these colors

    ### Scale the Matplotlib window size
    default_window_hight = 14.5 # 8.5
    default_window_width = 16
    fig = plt.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    color_bar_w = 0.015 ### Sufficient size to show

    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    #if row_method != None: w1 =
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
    if column_method != None:
        start_time = time.time()
        d2 = dist.pdist(x.T)
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
        Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
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
        print 'Column clustering completed in %s seconds' % time_diff
    else:
        ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data

    # Compute and plot left dendrogram.
    if row_method != None:
        start_time = time.time()
        d1 = dist.pdist(x)
        D1 = dist.squareform(d1)  # full matrix
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
        Y1 = sch.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z1 = sch.dendrogram(Y1, orientation='right')
        ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        ax1.set_xticks([]) ### Hides ticks
        ax1.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print 'Row clustering completed in %s seconds' % time_diff
    else:
        ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data

    # Plot distance matrix.
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
    xt = x
    if column_method != None:

        ## create a list containing the order of the rearranged matrix:
        idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
        xt = xt[:,idx2]
        ind2 = [ ind2[i] for i in idx2] # similarly organises the cluster assignment list
        #ind2 = ind2[:,idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
    if row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        xt = xt[idx1,:]   # xt is transformed x
        ind1 = [ ind1[i] for i in idx1 ]
        #ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
    axm.set_xticks([]) ### Hides x-ticks
    axm.set_yticks([])

    # Add text for row data, rearrange column headers according to clustering:
    new_row_header=[]
    new_column_header=[]
    for i in range(x.shape[0]):
        if row_method != None:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i, '  ' + row_header[idx1[i]])
            new_row_header.append(row_header[idx1[i]])
        else:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i, '  ' + row_header[i]) ### When not clustering rows
            new_row_header.append(row_header[i])
    for i in range(x.shape[1]):
        if column_method != None:
            #axm.text(i, -0.9, ' ' + column_header[idx2[i]], rotation=270, verticalalignment="top") #  rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            #axm.text(i, -0.9, ' ' + column_header[i], rotation=270, verticalalignment="top")
            new_column_header.append(column_header[i])




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

        print "Performing GO enrichment analysis for %d groups" % (len(genelistd))
        out_h = open(filename[:-4] + ".GO_enrichment.list", 'w')

        out_h.write("Group GOterm P-value\n")
        go_monster = genematch.GO_maker()
        # perform hypergeometric test (one-sided fishers exact test) on each group:
        for group in genelistd:
            gops = genematch.go_enrichment(genelistd[group])
            for goterm in gops:
                if gops[goterm] < 0.05: # this is a very weak threshold. Will need to do multiple test adjustment!
                    out_h.write( "%-4s %-7s %.5f %s\n" % (group, goterm, gops[goterm], str(go_monster.define_go(goterm))) )
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

        print "Performing KEGG pathway enrichment analysis for %d groups" % (len(genelistd))
        out_h = open(filename[:-4] + ".KEGG_enrichment.list", 'w')

        out_h.write("Group KEGG pathway P-value\n")

        for group in genelistd:
            pathway_ps = genematch.kegg_pathway_enrichment(genelistd[group])
            #gops = genematch.go_enrichment(genelistd[group])
            for ko in pathway_ps:
                if pathway_ps[ko] < 0.05: # this is a very weak threshold. Will need to do multiple test adjustment!
                    out_h.write( "%-4s %-7s %.5f %s\n" % (group, goterm, pathway_ps[ko]) )
        out_h.close()



    # Add text for column data:
    for i in range(x.shape[1]):
        if len(column_header)<200:
            axm.text(i - 0.3, -0.9, ' ' + new_column_header[i], rotation=270, verticalalignment="top")


    # Plot colside colors
    # axc --> axes for column side colorbar
    if column_method != None:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        dc = numpy.array(ind2, dtype=int)
        dc.shape = (1,len(ind2))
        im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
        axc.set_xticks([]) ### Hides ticks
        axc.set_yticks([])

    # Plot rowside colors
    # axr --> axes for row side colorbar
    if row_method != None:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        dr = numpy.array(ind1, dtype=int)
        dr.shape = (len(ind1),1)
        #print ind1, len(ind1)
        cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
        axr.set_xticks([]) ### Hides ticks
        axr.set_yticks([])

    # Plot color legend
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
    axcb.set_title("colorkey")


    pdfname = filename[:-4] + '.pdf'
    cb.set_label("Differential Expression (log2 fold)")
    ind1 = ind1[::-1] # reverse order of flat cluster leaves to match samples in txt file.
    exportFlatClusterData(pdfname, new_row_header,new_column_header,xt,ind1,ind2)

    ### Render the graphic
    if len(row_header)>50 or len(column_header)>50 or max(len(term) for term in new_column_header) > 80:
        plt.rcParams['font.size'] = 7 #5
    else:
        plt.rcParams['font.size'] = 10 #8

    plt.savefig(pdfname)
    print 'Exporting:',pdfname
    filename = filename[:-4]+'.png'
    plt.savefig(pdfname, dpi=200) #,dpi=100
    if display:
        plt.show()


    return new_column_header, ind2  # this is to allow group-specific KEGG enrichment analysis

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=x.max()
    vmin=x.min()
    if vmax<0 and vmin<0: direction = 'negative'
    elif vmax>0 and vmin>0: direction = 'positive'
    else: direction = 'both'
    if direction == 'both':
        vmax = max([vmax,abs(vmin)])
        vmin = -1*vmax
        return vmax,vmin
    else:
        return vmax,vmin

################# Export the flat cluster data #####################################

def export_expression_table(filename, new_column_header, new_row_header, xt):
    nfilename = filename[:-4] + '.data.tbl'
    "for when you need to create the data table but it's too big to analyse..."
    export_text = open(nfilename,'w')
    column_header = string.join(['gene']+new_column_header,'\t')+'\n' ### format column-names for export
    export_text.write(column_header)

    ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
    new_row_header = new_row_header[::-1]
    xt = xt[::-1]

    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_text.write(string.join([new_row_header[i]]+map(str, row),'\t')+'\n')
        i+=1
    export_text.close()

    ### Transpose text file for easier reading!
    oldfile_h = open(nfilename, 'rb')

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

    newfile_h = open(nfilename + ".transposed.txt" , 'w')
    for row in t_array:
        #print "The row is currently: %r" % row
        newfile_h.write("\t".join(row) + "\n")
    newfile_h.close()

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

################# Create Custom Color Gradients ####################################
#http://matplotlib.sourceforge.net/examples/plt_examples/custom_cmap.html

def RedBlackSkyBlue():
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

    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),

             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackGreen():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),

             'blue': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),

             'green':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def YellowBlackBlue():
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
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

################# Matrix manipulation methods ######################################

def normaliseData(x, center=True, norm_var=True, log_t=True, sample_norm=False):
    "center, normalize and/or log transform the array x"
    n = len(x[0]); m = len(x) # m = 6 samples, n = 6000 genes

    print "\n\nInitial minimum value in matrix: ", x.min()
    if log_t:
        for g in range(n):
            for i in range(m):
                x[:,g][i] = numpy.log2(x[:,g][i] + 1)

    print "log2(FPKM + 1) transformed data. New minimum value in matrix: ", x.min()


    meanlist = []   # to store for later re-adjustment of matrix
    if center and sample_norm:
        print "centering data for each sample to have mean 0"
        for s in range(m):
            ave_g = numpy.mean(x[s])
            for i in range(n):
                x[s,i] = x[s,i] - ave_g
                meanlist.append(ave_g)
    elif center:
        print "Centering data for each gene to have mean 0"
        for g in range(n):
            # center the data so it has mean = 0 for each gene:
            ave_g = numpy.mean(x[:,g])
            for i in range(m):
                x[:,g][i] = x[:,g][i] - ave_g
                meanlist.append(ave_g)
    print "new minimum value in matrix: ", x.min()

    if norm_var and sample_norm:
        print "Normalising data so each sample has standard deviation of 1"
        for s in range(m):
            #sum_sq = sum([ a ** 2 for a in (x[m])])
            stdev = numpy.std(x[s])
            for i in range(n):
                x[s,i] = x[s,i] / stdev
    elif norm_var:
        print "Normalising data so each gene has standard deviation of 1"
        for g in range(n):
            #sum_sq = sum([ a ** 2 for a in (x[:,g])])
            stdev = numpy.std(x[:,g])
            #print "1) mean: %.4f st.dev: %.4f" % ( numpy.mean(x[:,g]), stdev )
            for i in range(m):
                x[:,g][i] = x[:,g][i] / stdev # dividing by variance gives identical result
            #print "2) mean: %-7.4f st.dev: %-7.4f min: %-7.4f max: %-7.4f" % ( numpy.mean(x[:,g]), numpy.std(x[:,g]), numpy.min(x[:,g]), numpy.max(x[:,g]) )
    print "new minimum value in matrix: ", x.min()

    return meanlist # it is not nessary to return x, since it is modifying x

def filter_genes(x, column_header, row_header, genefile, col_num=0):
    """takes a file and extracts the genes in the col_num specified column,
    then uses it to filter the matrix"""

    # genefile can be given as a list of genes (say, from find_degs()... ), or as
    # a path to a file containing a list of genes.
    # The following builds a dictionary of genes from either input:
    if type(genefile) is list:   # allows filtering from hicluster generated list of results.
        genelist = {}.fromkeys(genefile,1)
    elif type(genefile) is dict:
        genelist = genefile # this of course will only work if the keys are the genes!
    else:   # assuming a filepath...
        genefile_h = open(genefile, 'rb')
        genelist = {}   # will be a dict of names { 'Cbir01255':1, 'CbirVgq':1, ... }
                        # used a dictionary to automatically remove any name duplications
        filegen = [ line.split() for line in genefile_h ]

        genefile_h.close()

        for colset in filegen:
            genelist[colset[col_num]]=1

    # create position-based list for filtering from gene-names:
    keeplist = []  # will be a position based list [ 0, 4, 23, 34, ... ]
    for gene in genelist:
        try:
            keeplist.append(column_header.index(gene))
        except:     # some columns will be headers that would cause errors.
            pass

    hitlist = range(len(column_header))
    for posn in keeplist: # will be the inverse of keeplist
        hitlist.remove(posn)

    ## filter matrix with genelist:
    y, column_header, row_header = filter_matrix(x, column_header, row_header, hitlist)

    return y, column_header, row_header

def filterData(x, column_header, row_header, mag=1.13, min_thresh=10, max_thresh=1000000):
    """filters out any gene for which the magnitude of expression is less than mag,
    or for which the maximum value of all samples is less than min_thresh,
    or for which the maximum value is greater than max_thresh"""

    n = len(x[0]); m = len(x) # m  samples, n  genes
    print "Filtering %d genes and %d samples:\nMin fold change: %.1f Min expression level (at least one sample): %d Max expression level: %d" % (n, m, mag, min_thresh, max_thresh)

    hitlist = []
    for g in range(n):
        fpkm_max = max(x[:,g])
        fpkm_min = min(x[:,g])
        size = numpy.absolute(fpkm_max/(fpkm_min + 0.00001))
        #rms  = numpy.sqrt( sum(a**2 for a in x[:,g])/m )
        if size < mag or fpkm_max < min_thresh or fpkm_max > max_thresh :
            hitlist.append(g)

    y, column_header, row_header = filter_matrix(x, column_header, row_header, hitlist)

    return y, column_header, row_header

def filter_matrix(x, column_header, row_header, hitlist):
    """given the positional list hitlist, will edit matrix and column headers
    by removing those positions"""

    # create new matrix and column_header without the columns in the hitlist:
    y = numpy.delete(x, hitlist, 1)
    print "Number of genes removed:", len(hitlist)
    for gene in hitlist[::-1]:
        column_header.pop(gene)

    n = len(y[0]); m = len(y) # m = 6 samples, n = 6000 genes
    print "there are now %d genes and %d samples" % (n, m)
    #print column_header

    return y, column_header, row_header

def filter_pretty(x, column_header, row_header, mag=2.5, group1="_F", group2="_S"):
    """searches for genes with magnitude mag different between two groups, as defined by
    group1 and group2 labels"""

    n = len(x[0]); m = len(x) # m  samples, n  genes
    print "Filtering %d genes and %d samples:\nMin fold change between groups: %.1f " % (n, m, mag)

    # split matrix into two based on the specified groups;
    group1list = []
    group1posns = []
    group2list = []
    group2posns = []
    for s in row_header:
        if re.search(group1,s) is not None:
            group1list.append(s)
            group1posns.append(row_header.index(s))
        elif re.search(group2, s) is not None:
            group2list.append(s)
            group2posns.append(row_header.index(s))
        else:       # group not found in sample list
            print "#" * 30
            print "%s could not be matched to group1 (%s) or group2 (%s)" % (s, group1, group2)
            print "#" * 30

    print group1list
    print group2list
    #print group1posns
    #print group2posns

    grouporder = group1posns + group2posns
    limit = len(group1posns)

    print "grouporder", grouporder
    print "limit", limit

    matrix_reord = x[grouporder,:]

    n = len(matrix_reord[0]); m = len(matrix_reord) # m  samples, n  genes
    print "Filtering %d genes and %d samples:\nMin fold change between groups: %.1f " % (n, m, mag)


    hitlist = []
    for g in range(n):
        mean1 = numpy.mean(matrix_reord[:limit,g])
        mean2 = numpy.mean(matrix_reord[limit:,g])
        ratio = mean1 / mean2
        if numpy.absolute( numpy.log2(ratio) / numpy.log2(mag) ) < 1 :
            hitlist.append(g)

    y, column_header, row_header = filter_matrix(matrix_reord, column_header, row_header, hitlist)

    return y, column_header, row_header

def reorder_matrix(x, column_header, row_header, groups=["_F","_S"]):
    "reorders rows such that they are sorted accordign to specified groups"

    n = len(x[0]); m = len(x) # m  samples, n  genes
    print "Sorting genes into %d groups" % (len(groups))

    # split matrix into two based on the specified groups;
    namelist = {}
    poslist = {}


    for s in row_header:
        for pattern in groups:
            if pattern not in namelist:
                namelist[pattern] = []
                poslist[pattern] = []
            # assign each sample to a group:
            if re.search(pattern,s) is not None:
                namelist[pattern].append(s)
                poslist[pattern].append(row_header.index(s))
                break
        else:       # group not found in any sample list
            print "#" * 30
            print "%s could not be matched to any group ==> %r" % (s, groups)
            print "#" * 30

    grouporder = []
    limits = []
    boundarystone = 0
    for pattern in groups:
        try:
            print pattern, namelist[pattern]
            grouporder += poslist[pattern]
            boundarystone += len(poslist[pattern])
        except KeyError:
            print pattern, "None found"
        limits.append(boundarystone)

    print "grouporder       :", grouporder
    print "limits (boundary):", limits

    matrix_reord = x[grouporder,:]
    row_header_new = [x for (y,x) in sorted(zip(grouporder,row_header))]

    return matrix_reord, column_header, row_header, limits

################# Data construction or import methods ##############################

def create_table(build_list):
    """Takes a comma-delimited list of files and uses them to build the fpkm table for
    analysis"""
    # extract list of files:
    all_files = build_list.split(',')

    output_file = os.getcwd() + "/fpkm_table.tbl"
    print "saving built table to file ", output_file

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

def importData(filename):
    start_time = time.time()
    matrix=[]
    row_header=[]
    first_row=True

    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]

    filename_h = open(filename, 'rb')
    for line in filename_h:
        t = line[:-1].split() ### remove end-of-line character - file is tab- or space-delimited
        if first_row:
            column_header = t[1:]
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
                if (abs(max(s)-min(s)))>0:
                    matrix.append(s)
                    row_header.append(t[0])


    # check appropriate number of items added to column header:
    if len(matrix[0]) != len(column_header):
        print "Fixing column header"
        filename_h = open(filename, 'rb')
        column_header = filename_h.next().split()[:]
        filename_h.close()

    time_diff = str(round(time.time()-start_time,1))
    try:
        print '\n%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    except Exception:
        print 'No data in input file.'; force_error

    return numpy.array(matrix), column_header, row_header

def clean_header(row_header):
    "removes long path name from row headers"

    new_headers = []
    for name in row_header:
        try:
            newname = re.findall('/([^/]*)', name)[-2]
        except:
            newname = name
        new_headers.append(newname)

    return new_headers

################# Data analyis methods #############################################

def analyse_pca(matrix, column_header, row_header, filename):
    """
    performs principal component analysis on the matrix and saves output as png.
    No additional normalisation occurs at this step - all normalisation has been done
    on the matrix prior to PCA.
    """
    # create data array and name array:
    A = matrix
    names = row_header

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
    U, s, Vt = numpy.linalg.svd(matrix, full_matrices=True)
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

    plt.show()

    return matrix_reduced

def expression_dist(matrix, column_header, row_header, filename, max_x=500, min_x=0, histo=False):
    "Creates histogram of gene expression"

    count = 0

    # create box and whisker plot of each sample:
    fig = plt.figure()
    subp = fig.add_subplot(111)
    something = subp.boxplot(matrix.transpose())

    # create labels for each sample:
    rng = [ i + 1 for i in range(len(row_header)) ]
    plt.xticks(rng , row_header)

    plt.show()

    if histo:
        for row in matrix:
            fig = plt.figure()
            subf = fig.add_subplot(111, title=row_header[count])
            n, bins, patches = subf.hist(row, 50, range=(min_x,max_x), histtype='stepfilled')
            count += 1
        plt.setp(labels, rotation=90)
        plt.show()

def find_degs(x, column_header, row_header, group1="_F", group2="_S"):
    "finds DEGs using t-test and returns dictionary { Gene:P-value }"

    """
    t_dict = {}
    n = len(x[0]); m = len(x) # m  samples, n  genes
    print "Sorting genes into group: %s and group: %s" % (group1, group2)

    # split matrix into two based on the specified groups;
    group1list = []
    group1posns = []
    group2list = []
    group2posns = []
    for s in row_header:
        if re.search(group1,s) is not None:
            group1list.append(s)
            group1posns.append(row_header.index(s))
        elif re.search(group2, s) is not None:
            group2list.append(s)
            group2posns.append(row_header.index(s))
        else:       # group not found in sample list
            print "#" * 30
            print "%s could not be matched to group1 (%s) or group2 (%s)" % (s, group1, group2)
            print "#" * 30

    print group1list
    print group2list
    #print group1posns
    #print group2posns

    grouporder = group1posns + group2posns
    limit = len(group1posns)

    print "grouporder", grouporder
    print "limit", limit

    matrix_reord = x[grouporder,:]
    """


    matrix_reord, column_header, reord_row_header, limits = reorder_matrix(x, column_header, row_header, groups=["_F","_S"])
    n = len(matrix_reord[0]); m = len(matrix_reord) # m  samples, n  genes

    print "Performing t-test"
    limit = limits[0]
    t_dict = {}
    print "limit for t-test",limit
    for g in range(n):
        t_val, p_val = stats.ttest_ind(matrix_reord[:limit,g], matrix_reord[limit:,g])
        t_dict[column_header[g]] = p_val

    return t_dict

def degs_anova(x, column_header, row_header, group1="_F", group2="_S"):
    "finds DEGs using ANOVA and returns dictionary { Gene:P-value }"

    matrix_reord, column_header, reord_row_header, limits = reorder_matrix(x, column_header, row_header, groups=["SL12", "SL24","SL48", "SL96","FP12", "FP24","FP48", "FP96", "FL", "SP"])     # groups=["SL12","SL24","SL48","SL96","FP12","FP24","FP48","FP96","SP12","FL12"])
    n = len(matrix_reord[0]); m = len(matrix_reord) # m  samples, n  genes

    # boundaries in matrix for each time point:
    #SL12,SL24,SL48,SL96,FP12,FP24,FP48,FP96,SP12,FL12 = limits[0], limits[1], limits[2], limits[3], limits[4], limits[5], limits[6], limits[7], limits[8], limits[9]
    #for i in range(len(limits)):


    print "Performing ANOVA"

    A_dict = {}
    print "limits for anova:",limits
    for g in range(n):
        f_val, p_val = stats.f_oneway(matrix_reord[:limits[0],g],\
            matrix_reord[limits[0]:limits[1],g], matrix_reord[limits[1]:limits[2],g],\
            matrix_reord[limits[2]:limits[3],g], matrix_reord[limits[3]:limits[4],g],\
            matrix_reord[limits[4]:limits[5],g], matrix_reord[limits[5]:limits[6],g],\
            matrix_reord[limits[6]:limits[7],g], matrix_reord[limits[7]:limits[8],g],\
            matrix_reord[limits[8]:limits[9],g])
        A_dict[column_header[g]] = p_val

    return A_dict

################# Miscellaneous methods #############################################

def appendgo(geneid, goterm, go_obj):
    if goterm in go_obj.findem(geneid.split()[0]):
        geneid = " ".join([geneid, goterm, go_obj.define_go(goterm)[0]])
    return  geneid

####################################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Performs heirarchical clustering")

    # input options
    parser.add_argument("-D", "--data_file", type=str, help="The data file for analysing")
    parser.add_argument("-B", "--build_table", type=str, dest="build_list", default=None, help="Provide a comma-delimited list of cufflinks files with which to build the fpkm table for analysis.")
    # output options
    parser.add_argument("-o", "--output_file", dest='filename', type=str, default=None, help="output file name for results")
    parser.add_argument("-e", "--export_table", action='store_true', help="export transformed expression matrix")
    # analysis options
    parser.add_argument("-z", "--row_header", type=str, dest="row_header", default=False, help="Not sure why this was put here")
    parser.add_argument("-R", "--row_method", type=str, dest="row_method", default='complete', help="The clustering method for rows \n(single, average, complete, etc)")
    parser.add_argument("-C", "--column_method", type=str, dest="column_method", default='complete', help="The clustering method for columns \n(single, average, complete, weighted, ward, centroid, etc)")
    parser.add_argument("-r", "--row_metric", type=str, dest="row_metric", default='correlation', help="The distance metric for rows \n(euclidean, correlation, cosine, manhattan, etc)")
    parser.add_argument("-c", "--column_metric", type=str, dest="column_metric", default='correlation', help="The distance metric for columns \n(euclidean, correlation, manhattan, etc)")
    parser.add_argument("-P", "--pca", action='store_true',  help="Performs principal component analysis.")
    parser.add_argument("-a", "--anova", type=float,  help="Perform ANOVA on 10 groups, and report genes with P-value less than value specified")
    parser.add_argument("-A", "--filter_anova", type=float, help="Perform ANOVA and filter genes to keep only those with P-value less than value given")
    parser.add_argument("-s", "--t_test", type=float,  help="Perform student's t-test on two groups, and report genes with P-value less than value specified")
    parser.add_argument("-S", "--filter_t", type=float, dest="ttest_thresh", help="Perform student's t-test and filter genes to keep only those with P-value less than value given")
    parser.add_argument("-K", "--kegg", action='store_true', help="Perform KEGG pathway enrichment analysis on gene clusters")
    parser.add_argument("-E", "--go_enrichment", action='store_true', help="Perform GO term enrichment analysis on gene clusters")
    # viewing options
    parser.add_argument("-g", "--color_gradient", type=str, dest="color_gradient", default='red_white_blue', help="The colour scheme \n(red_white_blue, red_black_sky, red_black_blue, \nred_black_green, yellow_black_blue, seismic, \ngreen_white_purple, coolwarm)")
    parser.add_argument("-d", "--distribution", action='store_true', default=False, help="Shows FPKM distribution of each sample before and after normalisation")
    parser.add_argument("--display_off", action='store_true', help="Turn of displaying of clustering")
    # filtering options
    parser.add_argument("-m", "--magnitude", type=float, dest="filter", default=2.5, help="Filters out genes with magnitude of range less than value given. Default = 1.13")
    parser.add_argument("-L", "--gene_list", type=str, dest="gene_list", default=None, help="Allows provision of a file containing a list of genes for inclusion in the clustering (ie, will filter out all genes NOT in the list provided). Otherwise, all genes will be included.")
    parser.add_argument("-p", "--fpkm_max", type=float, dest="fpkm_max", default=1000000, help="Filters out genes with maximum fpkm greater than value given. Default = 1 000 000")
    parser.add_argument("-q", "--fpkm_min", type=int, dest="fpkm_min", default=10, help="Filters out genes with maximum fpkm less than value given. Default = 10")
    parser.add_argument("-f", "--filter_off", action='store_true', help="Turns off filtering. ")
    parser.add_argument("-y", "--pretty_filter", type=float, default=None, help="Filters (non-normalised) matrix to remove genes whose mean value between treatments is less than value given. Try 2.5")
    parser.add_argument("--kill_PC1", action='store_true', help="removes first principal component")
    # data transform options
    parser.add_argument("-t", "--transform_off", action='store_true', default=False, help="Turns off log2(FPKM + 1) transformation (prior to normalisation if selected).")
    parser.add_argument("-n", "--normalise_off", action='store_true', default=False, help="Turns off normalisation. Normalises by dividing by the standard deviation.")
    parser.add_argument("-u", "--centering_off", action='store_true', default=False, help="Turns off gene centering. Centering subtracts the mean from all values for a gene, giving mean = 0.")
    parser.add_argument("-T", "--transpose", action='store_true',  help="Transpose the matrix. Columns should represent genes, Rows samples")
    parser.add_argument("-X", "--sample_norm", action='store_true', help='Normalises samples instead of genes')

    args = parser.parse_args()




    """ Running with cosine or other distance metrics can often produce negative Z scores
        during clustering, so adjustments to the clustering may be required.

    see: http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    see: http://docs.scipy.org/doc/scipy/reference/spatial.distance.htm
    color_gradient = red_white_blue|red_black_sky|red_black_blue|red_black_green|yellow_black_blue|green_white_purple'
    """

    ## create data table from cufflinks files:
    if args.build_list:
        data_table = create_table(args.build_list)
    else:
        data_table = args.data_file

    ## create log file of arguments:
    timestamp = time.strftime("%b%d_%H.%M")
    if not args.filename:
        root_dir = os.path.dirname(data_table)
        filename = root_dir + "/hicluster." + timestamp + ".log"
    else:
        filename = args.filename

    log_h = open(filename, 'w')
    log_h.write( "File created on %s\n" % (timestamp) )
    for arg in str(args).split():
        log_h.write( "%s\n" % (arg) )
    log_h.close()



    ## create matrix:
    matrix, column_header, row_header = importData(data_table)

    row_header = clean_header(row_header)

    ## transpose the matrix and swap the column and row headers.
    if args.transpose:
        matrix = numpy.transpose(matrix)
        tempcol = column_header
        column_header = row_header
        row_header = tempcol


    ####### FILTERING OF MATRIX ###########

    ## filter by differences between groups
    if args.pretty_filter:
        matrix, column_header, row_header = filter_pretty(matrix, column_header, row_header, mag=args.pretty_filter)

    ## filter by provided gene list
    if args.gene_list:
        print "Only keeping genes provided in", args.gene_list
        matrix, column_header, row_header = filter_genes(matrix, column_header, row_header, args.gene_list, col_num=0)

    ## filter by magnitude, minimum values and maximum values
    if not args.filter_off:
        matrix, column_header, row_header = filterData(matrix, column_header, row_header,  mag=args.filter, min_thresh=args.fpkm_min, max_thresh=args.fpkm_max)


    ####### NORMALISATION OPTIONS ##############

    ## show gene distribution
    if args.distribution:
        expression_dist(matrix, column_header, row_header, data_table)

    ## normalise matrix
    meanlist = normaliseData(matrix, center=not(args.centering_off), norm_var=not(args.normalise_off), log_t=not(args.transform_off), sample_norm=args.sample_norm)

    ## show distribution after normalisation
    if args.distribution:
        expression_dist(matrix, column_header, row_header, args.data_file, min_x=-2, max_x=2)

    ## export transformed data table:
    if args.export_table:
        export_expression_table(filename, column_header, row_header, matrix)

    ####### ANALYSIS OPTIONS ##############

    ## Principal Component Analysis:
    if args.pca:
        matrix_red = analyse_pca(matrix, column_header, row_header, data_table)

    if args.kill_PC1:
        matrix = matrix_red

    ## ANOVA analysis
    if args.filter_anova:
        a_list = []
        a_dict = degs_anova(matrix, column_header, row_header)
        out_h = open(filename[:-4] + ".ANOVA.list", 'w')
        for gene in a_dict:
            if a_dict[gene] <= args.filter_anova:
                out_h.write( "%-12s P-value: %.4f\n" % (gene, a_dict[gene]) )
                a_list.append(gene)
        out_h.close()
        print "Filtering matrix to %d genes with ANOVA P-value less than %.2f" % (len(a_list), args.filter_anova)
        matrix, column_header, row_header = filter_genes(matrix, column_header, row_header, a_list)
    elif args.anova:
        a_dict = degs_anova(matrix, column_header, row_header)
        out_h = open(filename[:-4] + ".ANOVA.list", 'w')
        for gene in a_dict:
            if a_dict[gene] <= args.anova:
                out_h.write( "%-12s P-value: %.4f\n" % (gene, a_dict[gene]) )
        out_h.close()

    ## t-test analysis
    if args.ttest_thresh:
        t_list = []
        t_dict = find_degs(matrix, column_header, row_header)
        out_h = open(filename[:-4] + ".t_test.list", 'w')
        for gene in t_dict:
            if t_dict[gene] <= args.ttest_thresh:
                out_h.write( "%-12s P-value: %.4f\n" % (gene, t_dict[gene]) )
                t_list.append(gene)
        out_h.close()
        print "Filtering matrix to %d genes with t-test P-value less than %.2f" % (len(t_list),args.ttest_thresh)
        matrix, column_header, row_header = filter_genes(matrix, column_header, row_header, t_list)
    elif args.t_test:
        t_dict = find_degs(matrix, column_header, row_header)
        out_h = open(filename[:-4] + ".t_test.list", 'w')
        for gene in t_dict:
            if t_dict[gene] <= args.t_test:
                out_h.write( "%-12s P-value: %.4f\n" % (gene, t_dict[gene]) )
        out_h.close()

    ## re-orient for publication purposes
    you_want = False    # change to True to swap the columns and rows (primarily visual).
    if you_want:
        matrix = numpy.transpose(matrix)
        tempcol = column_header
        column_header = row_header
        row_header = tempcol

    ## perform hierarchical clustering
    if len(matrix)>0:
        try:
            new_column_header, groups = heatmap(matrix, row_header, column_header, args.row_method, args.column_method, args.row_metric, args.column_metric, args.color_gradient, filename, display=not(args.display_off), kegg=args.kegg, go=args.go_enrichment)
        except Exception:
            print 'Error using %s ... trying euclidean instead' % args.row_metric

            try:
                new_column_header, groups = heatmap(matrix, row_header, column_header, args.row_method, args.column_method, 'euclidean', 'euclidean', args.color_gradient, filename, display=not(args.display_off), kegg=args.kegg, go=args.go_enrichment)
            except IOError:
                print 'Error with clustering encountered'
                new_column_header = ['']
                groups = ['']
"""
    if args.go_enrichment:

        leafpairs = zip(groups, new_column_header)
        genelistd = {}
        for group,geneid in leafpairs:
            try:
                genelistd[group].append(geneid)
            except KeyError:
                genelistd[group] = [geneid]

        print "Performing GO enrichment analysis for %d groups" % (len(genelistd))
        out_h = open(filename[:-4] + ".GO_enrichment.list", 'w')

        out_h.write("Group GOterm P-value\n")
        go_monster = genematch.GO_maker()

        for group in genelistd:
            #pathway_ps = genematch.kegg_pathway_enrichment(genelistd[group])
            gops = genematch.go_enrichment(genelistd[group])
            for goterm in gops:
                if gops[goterm] < 0.05:
                    out_h.write( "%-4s %-7s %.5f %s\n" % (group, goterm, gops[goterm], str(go_monster.define_go(goterm))) )
        out_h.close()

    if args.kegg:
        leafpairs = zip(groups, new_column_header)
        genelistd = {}
        for group,geneid in leafpairs:
            try:
                genelistd[group].append(geneid)
            except KeyError:
                genelistd[group] = [geneid]

        print "Performing KEGG pathway enrichment analysis for %d groups" % (len(genelistd))
        out_h = open(filename[:-4] + ".KEGG_enrichment.list", 'w')

        out_h.write("Group KEGG pathway P-value\n")

        for group in genelistd:
            pathway_ps = genematch.kegg_pathway_enrichment(genelistd[group])
            #gops = genematch.go_enrichment(genelistd[group])
            for ko in pathway_ps:
                if pathway_ps[ko] < 0.05:
                    out_h.write( "%-4s %-7s %.5f %s\n" % (group, goterm, pathway_ps[ko]) )
        out_h.close()
"""
