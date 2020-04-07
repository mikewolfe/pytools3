import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.cluster import hierarchy as sch
import argparse
#mpl.rc('ytick', labelsize=20)
#mpl.rc('xtick', labelsize=20)

def value_parser(value):
    """
    Parse a single value of the pvalue matrix in ipage
    Example:
    -1.046/-0.008
    Where the left value is significant enrichment and right value is
    significant de-enrichment

    """
    vals = value.split('/')
    vals = [float(val) for val in vals]
    min_index, min_value = min(enumerate(vals), key=lambda p : p[1])
    if min_index == 0:
        # make enrichments positive
        return -1*min_value
    else:
        return min_value

def pvmatrix_parser(pvmatrix_file, drop_col = None):
    """
    Parse the full ipage pvalue matrix
    """
    with open(pvmatrix_file, mode="r") as f:
        # remove the header
        f.readline()
        value_rows = []
        goterms = []
        for line in f:
            linearr = line.rstrip().split("\t")
            # parse all the values after the go terms
            vals = [value_parser(val) for val in linearr[1:]]
            if drop_col is not None:
                vals.pop(drop_col)
            value_rows.append(vals)
            # add the go term for that row
            goterms.append(linearr[0])
    return goterms, np.array(value_rows)

def add_missing_go(go_terms, go_term_set, values):
    """
    Add a value of zero for any missing go terms
    """
    # new go terms to add to current go terms
    new_go_terms = []
    this_go_terms = set(go_terms)
    # go through all of the total go terms
    for go_term in go_term_set:
        if go_term not in this_go_terms:
            # add go term if it is missing to end of new go term list
            new_go_terms.append(go_term)
    # add all the new go terms to the original list
    for go_term in new_go_terms:
        go_terms.append(go_term)
    # add zeros to the values for all of the new go terms
    values = np.vstack((values, np.zeros((len(new_go_terms), values.shape[1]))))
    return go_terms, values

def sort_go_terms(go_terms, values):
    # sort the values and go terms by the go terms
    values = values[np.argsort(go_terms),:]
    go_terms = np.sort(go_terms)
    return go_terms, values

def determine_clustering(data_matrix):
    # do heirarchical clustering of the go terms
    z = sch.linkage(data_matrix, method='average', metric='cityblock')
    d=sch.dendrogram(z,no_plot=True)
    return d['leaves']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a consolidate iPAGE heatmap")
    parser.add_argument("outfile", type=str, help="name of outputfile")
    parser.add_argument('--infiles', type=str, nargs="+", help="input pvalue matrices")
    parser.add_argument('--labels', type=str, nargs="+", default=None,
                        help="labels for each pvalue matrix")
    parser.add_argument('--category_labels', type=str, nargs="+", default=None,
                        help="category labels for each pvalue matrix")
    parser.add_argument('--pvalue_box', type=float, default=0.05, help="pvalue cutoff to put a box around")
    parser.add_argument('--max_pvalue', type=float, default=10, help="max log10 pvalue to show color for")
    parser.add_argument('--drop_category', type=int, default=None, help="drop a single category by number")
    parser.add_argument('--filter_gos', type=str, default=None, help="only plot these go terms")
    parser.add_argument('--cat_rotation', type=int, default=None, help="rotate category labels (int)")
    parser.add_argument('--fig_ratio', nargs=2, type=float,default=[2.5, 5], help="control size of figure")
    args = parser.parse_args()

    outfile = args.outfile
    infiles = args.infiles

    if args.labels is None:
        labels = np.arange(0, len(args.infiles))
    else:
        labels = args.labels
    go_terms, values = zip(*[pvmatrix_parser(this_file, drop_col = args.drop_category) for this_file in infiles ])
    if args.category_labels is None:
        cat_labels = np.arange(0, values[0].shape[1])
    else:
        cat_labels = args.category_labels
    if args.filter_gos is not None:
        good_gos = []
        with(open(args.filter_gos)) as inf:
                for line in inf:
                    good_gos.append(line.rstrip())
    else:
        good_gos = None
    go_terms = list(go_terms)
    values = list(values)
    go_term_set = set()
    for this_go in go_terms:
        go_term_set = go_term_set.union(set(this_go))
    
    # add missing go terms and sort all the go terms to be the same order:
    for i in range(len(go_terms)):
        this_go, this_val = add_missing_go(go_terms[i], go_term_set, values[i])
        values[i] = this_val
        go_terms[i] = this_go
        this_go, this_val = sort_go_terms(go_terms[i], values[i])
        values[i] = this_val
        go_terms[i] = this_go
    
    cluster_mat = np.column_stack(values)
    cluster = determine_clustering(cluster_mat)
    if good_gos is not None:
        new_cluster = []
        for val in cluster:
            if go_terms[1][val] in good_gos:
                new_cluster.append(val)
        cluster = np.array(new_cluster)

    
    # control the size of the figure
    # # make len(go_terms[0][cluster])/5 to go back to old way
    fig = plt.figure(figsize=(args.fig_ratio[0]*len(values),len(go_terms[0][cluster])/args.fig_ratio[1]))
    width_ratios = [5]*len(infiles)
    width_ratios.append(1)
    gs = gridspec.GridSpec(1, len(infiles)+1, width_ratios=width_ratios)
    ## set the maximum vals displayed on the heatmap, uncomment to see full range
    #max_val = np.max([np.max(this_val) for this_val in values])
    #min_val = np.min([np.min(this_val) for this_val in values])
    #val = np.max([abs(min_val), abs(max_val)])
    
    # for now maximum val displayed is -10 to 10
    val = args.max_pvalue
    
    # Plot each condition as a seperate heatmap
    for i in range(len(infiles)):
        ax = plt.subplot(gs[0,i])
        mat = ax.matshow(values[i][cluster,:], cmap=plt.cm.RdBu_r, vmin=(-1*val), vmax=val, aspect="auto", interpolation='nearest')
        # set up all the boxes around all the ticks
        plt.yticks(np.arange(0,values[i][cluster,:].shape[0]))
        ax.set_xticks(np.arange(-0.5, values[i][cluster,:].shape[1],1), minor=True)
        ax.set_yticks(np.arange(-0.5, values[i][cluster,:].shape[0],1), minor=True)
        if args.cat_rotation > 0:
            ha = "left"
        elif args.cat_rotation < 0:
            ha = "right"
        else:
            ha = "center"
        ax.set_xticklabels(cat_labels, rotation = args.cat_rotation, ha = ha)
        ax.set_title(labels[i])
        #ax.grid(True, which="minor", axis="both", ls="-", color="white", lw=2)
        plt.xticks(np.arange(0,values[i].shape[1]))
        print(values[i].shape[0])
        sigboxes = np.argwhere(np.abs(values[i][cluster,:]) > np.abs(np.log10(args.pvalue_box)))
        print(len(sigboxes))
        for box in sigboxes:
            box_x = box[1] - 0.5
            box_y = box[0] - 0.5
            ax.add_patch(patches.Rectangle((box_x, box_y), 1, 1, linewidth=2, edgecolor='black', facecolor='none'))
        if i == 0:
            ax.set_yticklabels(go_terms[i][cluster])
        else:
            ax.tick_params(labelleft='off')
        plt.yticks(rotation=0)
    plt.colorbar(mat, cax = plt.subplot(gs[0,i+1]))
    plt.savefig(outfile, bbox_inches='tight')
