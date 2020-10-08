import fastq
import sys
import numpy as np
import argparse
import gzip

# If I want sampling with replacement, I will need  to read the whole file in
# and use this function
def sample_fastq_names(fqfile, fraction, replacement = False):
    size = len(fqfile)
    outsize = int(size * fraction)
    all_names = fqfile.names()
    sample = np.random.choice(np.arange(0,size),outsize,replace = replacement)
    return(all_names[sample])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Sample a (pair) of fastq file(s)")
    parser.add_argument("infile", type=str, help="path to input file")
    parser.add_argument("outfile", type = str, help="output prefix, can include path")
    parser.add_argument("--infilepair", type=str, help="path to input file pair")
    parser.add_argument("--outfilepair", type=str, help="path to input file pair")
    parser.add_argument('--fraction', type=float, default = 0.05, help= "fraction of reads to sample, default = 0.05")
    parser.add_argument('--seed', type=int, help= "seed for random number generator")
#    parser.add_argument('--withreplacement', action = "store_true",
#            help = "sample with replacement?")
    args = parser.parse_args()

    read1 = fastq.FastqFile(args.infile)
    if args.infilepair:
        if not args.outfilepair:
            raise ValueError("Must specify output file for R2 pair")
        read2 = fastq.FastqFile(args.infilepair)
        names = []

    if args.seed:
        np.random.seed(args.seed)
    if type(args.outfile) is str and args.outfile.endswith(".gz"):
        read1out = gzip.open(args.outfile, mode = "wb")
        prep_out = lambda x: str(x).encode()
    else:
        read1out = open(args.outfile, mode = "w")
        prep_out = lambda x: str(x)
    
    # paired processing
    if args.infilepair:     
        if type(args.outfilepair) is str and args.outfilepair.endswith(".gz"):
            read2out = gzip.open(args.outfilepair, mode = "wb")
            prep2_out = lambda x: str(x).encode()
        else:
            read2out = open(args.outfilepair, mode = "w")
            prep2_out = lambda x: str(x)

        for entry1, entry2 in zip(read1.stream_file(), read2.stream_file()):
            # flip a coin
            if np.random.binomial(1,args.fraction):
                read1out.write(prep_out(entry1))
                read2out.write(prep2_out(entry2))
        read1out.close()
        read2out.close()
    # single-end processing
    else:
        for entry in read1.stream_file():
            # flip a coin
            if np.random.binomial(1,args.fraction):
                read1out.write(prep_out(entry))
        read1out.close()
