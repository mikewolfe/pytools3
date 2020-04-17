import sys
import argparse
import fasta as fa
import numpy as np
import bed_utils as bed
import intervaltools as it
import logging
import random

class FIREfile(object):

    def __init__(self):
        self.data = {}
        self.names = []

    def __iter__(self):
        for name in self.names:
            yield (name, self.data[name])

    def __len__(self):
        return len(self.names)

    def __add__(self, other):
        newfile = FIREfile()
        for name, score in self:
            newfile.add_entry(name, score)
        for name, score in other:
            newfile.add_entry(name, score)
        return newfile

    def add_entry(self, name, score):
        self.data[name] = score
        self.names.append(name)

    def pull_value(self, name):
        return self.data[name]

    def discretize_quant(self, nbins=10):
        # first pull all the values
        all_vals = [val for name, val in self]
        all_vals = np.array(all_vals)
        quants = np.arange(0,100, 100.0/nbins)
        bins = []
        for quant in quants:
            bins.append(np.percentile(all_vals, quant))
        all_vals = np.digitize(all_vals, bins)
        for new_val, (name, old_val) in zip(all_vals, self):
            self.data[name] = new_val

    def shuffle(self):
        size = len(self)
        subset = np.random.permutation(size)
        shuffled_names = [self.names[val] for val in subset]
        self.names = shuffled_names


    def write(self, fname):
        with open(fname, mode="w") as outf:
            outf.write("name\tscore\n")
            for name, score in self:
                outf.write("%s\t%s\n"%(name, score))


def make_kfold_datasets(k, fastafile, firefile, outpre):
    logging.warning("Making %s fold datasets"%k)
    # shuffle all the fasta entries in place 
    np.random.shuffle(fastafile.names)
    # create k folds out of the names as a list of lists
    size = len(fastafile.names)
    folds = [firefile.names[i::k] for i in range(k)]
    # loop through each fold
    for test_idx in range(len(folds)):
        # make a seperate fasta file and fire file for the test fold
        this_test_fasta = fa.FastaFile()
        this_test_fire = FIREfile()
        for name in folds[test_idx]:
            this_test_fasta.add_entry(fastafile.pull_entry(name))
            this_test_fire.add_entry(name, firefile.pull_value(name))
        with open(outpre+"_test_%i.fa"%test_idx, mode="w") as outf:
            this_test_fasta.write(outf)
        this_test_fire.write(outpre+"_test_%i.txt"%test_idx)
        

        # make a sperate fasta and fire file for the train folds
        this_train_fasta = fa.FastaFile()
        this_train_fire = FIREfile()
        
        for train_idx in range(len(folds)):
            if train_idx != test_idx:
                for name in folds[train_idx]:
                    this_train_fasta.add_entry(fastafile.pull_entry(name))
                    this_train_fire.add_entry(name, firefile.pull_value(name))

        logging.warning("Writing fold %i"%test_idx)
        with open(outpre+"_train_%i.fa"%test_idx, mode="w") as outf:
            this_train_fasta.write(outf)
        this_train_fire.write(outpre+"_train_%i.txt"%test_idx)


def determine_start_end(feature, chrm, args):
    if args.centered:
        center = round((feature["start"] + feature["end"])/2)
        start = center
        end = center
    else:
        start = feature["start"]
        end = feature["end"]

    if feature["strand"] == "-":
        final_start = max(start - args.downstream, 0)
        final_end = min(end + args.upstream, len(chrm))
        final_rc = True
    else:
        final_start = max(start - args.upstream, 0)
        final_end = min(end + args.downstream, len(chrm))
        final_rc = False
    return(final_start, final_end, final_rc)

class ExcludedSearch(object):
    def __init__(self):
        self.exclusion_dict = {}

    def add_interval(self, chrm, start, end):
        intervals = self.exclusion_dict.get(chrm, it.Intervals())
        intervals.add_interval(it.Interval(start, end)) 
        self.exclusion_dict[chrm] = intervals

    def rand_loc(self, chrm, start, end, padding, max_trys=1000000):
        exclude = self.exclusion_dict.get(chrm, it.Intervals())
        val = np.random.randint(start, end)
        new_int = it.Interval(val-padding, val+padding)
        trys = 0
        while exclude.check_overlap(new_int) and trys <= max_trys:
            val = np.random.randint(start, end)
            new_int = it.Interval(val-padding, val+padding)
            trys +=1
        if trys >= max_trys:
            raise RuntimeError("Couldn't find a non-excluded location after %s tries"%max_trys)
        return val


def discretize_genome(wsize, stepsize, genome):
    out_bed = bed.BedFile()
    for chrm in genome:
        chrm_name = chrm.chrm_name()
        for i in range(0, len(chrm)-wsize, stepsize):
            start = i
            end = i + wsize
            name = "%s_%i_%i"%(chrm_name, start, end) 
            out_bed.add_entry(bed.BedEntry([chrm_name, start, end, name, 0]))
    return out_bed



def dense_sampling_main(args):

    # parse arguments
    args = parser.parse_args()
    # figure out random seed
    np.random.seed(args.seed)
    # read in genome
    genome = fa.FastaFile()
    logging.warning("reading in full genome")
    with open(args.fasta) as inf:
        genome.read_whole_file(inf)

    # read in bed file
    inbed = bed.BedFile()
    logging.warning("reading in bed")
    inbed.from_bed_file(args.bedfile)

    # discretize the genome by size of window and step of window
    outbed = discretize_genome(args.wsize, args.stepsize, genome)

    # convert input bed to an interval file by chromosome
    intervals = {chrm.chrm_name() : it.Intervals() for chrm in genome}

    for feature in inbed:
        this_interval = intervals[feature["chrm"]]
        this_interval.add_interval(it.Interval(feature["start"], feature["end"]))
        intervals[feature["chrm"]] = this_interval

    # figure out which intervals overlap and which don't

    logging.warning("determining which intervals overlap")
    positive_bed = bed.BedFile()
    negative_bed = bed.BedFile()
    for i, window in enumerate(outbed):
        if i % 10000 == 0:
            logging.warning("Checking interval %s"%i)
        this_chrm = window["chrm"]
        this_intervals = intervals[this_chrm]
        window_interval = it.Interval(window["start"], window["end"])
        perc_overlap = this_intervals.check_percent_overlap(window_interval)
        if perc_overlap >= args.perc_overlap:
            positive_bed.add_entry(window)
        else:
            negative_bed.add_entry(window)
    # make fire file
    fire = FIREfile()
    out_fasta = fa.FastaFile()
    for feature in positive_bed:
        this_chrm = feature["chrm"]
        this_name = feature["name"]
        this_start = feature["start"]
        this_end = feature["end"]
        fire.add_entry(this_name, args.default_score)
        out_fasta.add_entry(
                fa.FastaEntry(">" + this_name, genome.pull_entry(this_chrm).pull_seq(this_start, this_end)))

    for feature in negative_bed:
        this_chrm = feature["chrm"]
        this_name = feature["name"]
        this_start = feature["start"]
        this_end = feature["end"]
        fire.add_entry(this_name, args.rand_score)
        out_fasta.add_entry(
                fa.FastaEntry(">" + this_name, genome.pull_entry(this_chrm).pull_seq(this_start, this_end)))

    # write files
    
    if args.true_bed:
        positive_bed.write_bed_file(args.outpre + "_true.bed")
    if args.rand_bed:
        negative_bed.write_bed_file(args.outpre + "_rand.bed")

    if not args.no_fasta:
        with open(args.outpre + ".fa") as outf:
            out_fasta.write(outf)
        
    fire.write(args.outpre + "_fire.txt")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Take a bed format file and pull seqs for FIRE format")
    parser.add_argument('bedfile', type=str, help="input bed file")
    parser.add_argument('fasta', type=str, help="input fasta file")
    parser.add_argument('outpre', type=str, help="output file prefix")
    parser.add_argument('--centered', action="store_true", 
            help="make window centered on feature center")
    parser.add_argument('--upstream', type = int, default = 0,
            help = "size of window upstream of specified start location or center if centered is specified")
    parser.add_argument('--downstream', type = int, default = 0,
            help = "size of window downstream of specified end location or center if centered is specified")
    parser.add_argument('--kfold', type=int, default=None, 
                        help="create k fold datasets for CV")
    parser.add_argument('--nrand', type=int, default=0, 
            help="multiplier for number of random seqs to include")
    parser.add_argument('--seed', type=int, default=1234, 
            help="random seed for reproducibility")
    parser.add_argument('--rmchr', action="store_true", default=False, 
            help="rm chr string from peak chromosomes")
    parser.add_argument('--discretize', type=int, default=None, 
            help="discretize by value field?")
    parser.add_argument('--rand_score', type=int, default=0, 
            help="score to give for random seqs")
    parser.add_argument('--default_score', type=int, default=1, 
            help="score to give if the score is missing in the bed file")
    parser.add_argument('--rand_bed', action = "store_true",
            help="write a bed file for the random peaks?")
    parser.add_argument('--true_bed', action = "store_true",
            help="write a bed file for the true peaks?")
    parser.add_argument('--allow_rand_overlap', action = "store_true",
            help="allow random peaks to overlap? Default no overlap")
    parser.add_argument('--allow_true_overlap', action = "store_true",
            help="allow random peaks to overlap true peaks? Default no overlap")
    parser.add_argument('--no_fasta', action = "store_true",
            help="don't write fasta file")
    parser.add_argument('--dense', action = "store_true",
            help="Run program in dense mode")
    parser.add_argument('--wsize', type = int, default = 100,
            help="Size of window for dense mode, default = 100")
    parser.add_argument('--stepsize', type = int, default = 25,
            help="Step size for dense mode, default = 25")
    parser.add_argument('--perc_overlap', type = float, default = 0.5,
            help="Overlap cutoff for dense mode, default = 0.5")


    # parse arguments
    args = parser.parse_args()
    if args.dense:
        dense_sampling_main(args)
        sys.exit()
    # figure out random seed
    np.random.seed(args.seed)
    # read in genome
    genome = fa.FastaFile()
    logging.warning("reading in full genome")
    with open(args.fasta) as inf:
        genome.read_whole_file(inf)

    # read in bed file
    inbed = bed.BedFile()
    logging.warning("reading in bed")
    inbed.from_bed_file(args.bedfile)

    # check how much of the genome the regions cover
    genome_length = {}
    total_length = 0
    for chrm in genome:
        genome_length[chrm] = len(chrm)
        total_length += len(chrm)
    feature_length = 0
    for feature in inbed:
        this_chrm = genome.pull_entry(feature["chrm"])
        this_start, this_end, this_rc = determine_start_end(feature, this_chrm, args)
        feature_length += this_end - this_start

    fraction = feature_length / total_length
    if fraction >= 0.05:
        logging.warning("Input regions cover %0.2f %% of the genome. Consider using --dense"\
                %(fraction*100))
        

    # remove chr from chromosome names
    if args.rmchr:
        for feature in inbed:
            feature["chrm"] = feature["chrm"].replace("chr", "")

    # make output file
    outfasta = fa.FastaFile()
    realfire = FIREfile()
    fakefire = FIREfile()
    if args.rand_bed:
        rand_bed = bed.BedFile()

    # iterate over the features to create a database of intervals that are
    # already spanned (thus cannot be overlapped with random intervals)
    if args.nrand > 0:
        excluded_search = ExcludedSearch()
        if not args.allow_true_overlap:
            for feature in inbed:
                this_chrm = genome.pull_entry(feature["chrm"])
                this_start, this_end, this_rc = determine_start_end(feature, this_chrm, args)
                excluded_search.add_interval(feature["chrm"], 
                        this_start, this_end)

    for i, feature in enumerate(inbed):
        this_entry = fa.FastaEntry()
        this_chrm = genome.pull_entry(feature["chrm"])
        this_start, this_end, this_rc = determine_start_end(feature, this_chrm, args)
        this_entry.set_seq(this_chrm.pull_seq(this_start, this_end, rc = this_rc))
        this_name = "%s_%i_%i_%s_%s"%(feature["chrm"], 
                feature["start"], feature["end"], feature["strand"], feature["name"])

        this_entry.set_header(">"+this_name)
        outfasta.add_entry(this_entry)

        # have the ability to process files with no score
        this_score = feature["score"]
        if this_score == ".":
            this_score = args.default_score
        realfire.add_entry(this_entry.chrm_name(), this_score)

        # write out the true score instead of the discretized
        if args.true_bed:
            feature["name"] = this_name
            feature["strand"] = feature["strand"]
            feature["score"] = this_score
        
        # add random sequences if needed
        for j in range(args.nrand):
            this_wsize = round((this_end - this_start)/2)
            this_rand = fa.FastaEntry()
            this_seq = "N"*(this_end - this_start)
            while this_seq.count("N") > this_wsize:
                this_loc = excluded_search.rand_loc(feature["chrm"], 
                        0+this_wsize, len(this_chrm)-this_wsize, this_wsize)

                this_seq = this_chrm.pull_seq(this_loc-this_wsize, 
                        this_loc+this_wsize, rc = this_rc)
            this_rand.set_seq(this_seq)
            this_rand_name = this_name + "RAND%i"%j
            this_rand.set_header(">"+this_rand_name)
            outfasta.add_entry(this_rand)
            fakefire.add_entry(this_rand.chrm_name(), args.rand_score)
            if not args.allow_rand_overlap:
                excluded_search.add_interval(feature["chrm"], 
                        this_loc-this_wsize, this_loc+this_wsize)

            if args.rand_bed:
                rand_bed.add_entry(
                        bed.BedEntry(
                            [feature["chrm"], this_loc-this_wsize, this_loc+this_wsize,
                               this_rand_name, args.rand_score, feature["strand"]]))

    if args.discretize:
        realfire.discretize_quant(args.discretize)

    # if fake fire is empty it won't hurt anything
    finalfire = realfire + fakefire
    if args.kfold:
        finalfire.shuffle()
        make_kfold_datasets(args.kfold, outfasta, finalfire, args.outpre)
    else:
        with open(args.outpre+".fa", mode="w") as outf:
            outfasta.write(outf)
        finalfire.write(args.outpre+"_fire.txt")
    if args.rand_bed:
        rand_bed.write_bed_file(args.outpre + "_rand.bed")
    if args.true_bed:
        inbed.write_bed_file(args.outpre + "_true.bed")

