import fasta as fa
from Bio.SeqUtils import MeltingTemp as mt
import numpy as np

def oligo_to_Tm(seq, c_seq = None, nn_table = mt.R_DNA_NN1, Na = 150, Tris = 0, Mg = 0, dNTPs = 0, dnac1  = 25, dnac2 = 25):
    """
    Using `biopythons` Tm nearest neighbor calculation for determining melting
    temperatures between DNA oligo in c_seq and sequence

    """
    try:
        return mt.Tm_NN(seq, c_seq = c_seq, Na = Na, Tris = Tris, Mg = Mg, dNTPs = dNTPs, nn_table = nn_table, dnac1 = dnac1, dnac2 = dnac2)
    except:
        ValueError
        return "NaN"

def hamming_distance(seq, oligo):
    return sum(b1 != b2 for b1, b2 in zip(seq, oligo))

def linear_weights(size, max_weight, slope):
    weights = [max_weight]
    this_weight = max_weight
    for val in range(1, size):
        weights.append(this_weight/slope)
        this_weight /= slope

    weights = np.array(weights)
    # make weights sum to size
    weights = (weights/sum(weights))*size
    print(weights)
    return weights

def triangle_weights(size, max_weight, slope):
    weights = [max_weight]
    this_weight = max_weight
    # descending arm
    for val in range(1, size//2):
        weights.append(this_weight/slope)
        this_weight /= slope
    # Ascending arm
    weights2 = weights[::-1]
    # add middle zero if odd
    if size % 2 > 0:
        weights = weights.append(0)
    weights.extend(weights2)
    weights = np.array(weights)

    # make weights sum to size
    weights = (weights/sum(weights))*size
    print(weights)
    return weights

def weighted_hamming_distance(seq, oligo, weights):
    return sum([weight*(b1 != b2) for b1,b2,weight in zip(seq, oligo, weights)])


def seq_to_Tm(seq, nn_table = mt.DNA_NN3, Na = 150, Tris = 0, Mg = 0, dNTPs = 0, dnac1 = 25, dnac2 = 25):
    """
    Using `biopythons` Tm nearest neighbor calculation for determining melting
    temperatures

    Following along with what NEB uses for their Tm calculations this means
    using the nn_table for DNA from SantaLucia 1998 (DNA_NN3) but also using
    the Owczarzy et al 2004 salt correction formula

    For RNA:DNA hybrids use table R_DNA_NN1. Sequence must be the RNA.

    For RNA:RNA hybrids use table RNA_NN3.
    """
    return mt.Tm_NN(seq, Na = Na, Tris = Tris, Mg = Mg, dNTPs = dNTPs, nn_table = nn_table, saltcorr = 6, dnac1 = dnac1, dnac2 = dnac2)

def generate_windowed_tm(ffile, wsize, nn_table = mt.DNA_NN3, Na = 150, Tris = 0, Mg = 0, dNTPs = 0, strand = "+", dnac1 = 25, dnac2 = 25):
    out = {}
    for seq in ffile:
        i = 0
        max_size = len(seq)
        start = []
        end = []
        subseqs = []
        Tm = []
        for window in seq.window_iter(wsize):
            start.append(i)
            end.append(min(i + wsize, max_size))
            if strand == "-":
                this_window = fa.complement(window)[::-1]
            else:
                this_window = window
            subseqs.append(this_window)

            Tm.append(seq_to_Tm(this_window, nn_table = nn_table, Na = Na, Mg = Mg, Tris = Tris, dNTPs = dNTPs, dnac1 = dnac1, dnac2 = dnac2))
            i += 1
        out[seq.chrm_name()] = [start, end, subseqs, Tm]
    return out

def generate_windowed_oligo(ffile, oligo, strand = "+", tm_args = None):
    out = {}
    wsize = len(oligo)
    l_weights = linear_weights(wsize, 1, 1.5)
    t_weights = triangle_weights(wsize, 1, 1.5)
    for seq in ffile:
        i = 0
        max_size = len(seq)
        start = []
        end = []
        subseqs = []

        if tm_args is not None:
            Tm = []
        hamming = []
        lw_hamming = []
        tw_hamming = []

        for window in seq.window_iter(wsize):
            start.append(i)
            end.append(min(i + wsize, max_size))
            if strand == "-":
                this_window = fa.complement(window)[::-1]
            else:
                this_window = window
            subseqs.append(this_window)

            if tm_args is not None:
                Tm.append(oligo_to_Tm(this_window, c_seq = oligo, *tm_args))
            hamming.append(hamming_distance(this_window, oligo))
            lw_hamming.append(weighted_hamming_distance(this_window, oligo, l_weights))
            tw_hamming.append(weighted_hamming_distance(this_window, oligo, t_weights))
            i += 1

        if tm_args is not None:
            out[seq.chrm_name()] = [start, end, subseqs, hamming, lw_hamming, tw_hamming, Tm]
        else:
            out[seq.chrm_name()] = [start, end, subseqs, hamming, lw_hamming, tw_hamming]
    return out

def to_tsv_oligo(windowed_tm, outf, strand = "+"):
    with open(outf, mode = "w") as fhandle:
        for contig in windowed_tm.keys():
            for start, end, subseq, Tm, hamming, w_hamming in zip(*windowed_tm[contig]):
                fhandle.write("\t".join([contig, str(start), str(end), subseq, str(Tm), strand, str(hamming), str(w_hamming)]) + "\n")

def to_tsv(windowed_tm, outf, strand = "+"):
    with open(outf, mode = "w") as fhandle:
        for contig in windowed_tm.keys():
            for start, end, subseq, Tm, in zip(*windowed_tm[contig]):
                fhandle.write("\t".join([contig, str(start), str(end), subseq, str(Tm), strand]) + "\n")

def to_bw(windowed_tm, wsize, ffile, outfile):
    import bwtools as bw
    import numpy as np
    chrm_dict = {}
    Tm_dict = {}
    for seq in ffile:
        chrm_dict[seq.chrm_name()] = len(seq) - wsize + 1
        Tm_dict[seq.chrm_name()] = np.array(windowed_tm[seq.chrm_name()][3])
    bw.write_arrays_to_bigwig(outfile, Tm_dict, chrm_dict, res = 1)



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser("Tools to calculate Tms of sequences") 
    parser.add_argument('infile', type=str, 
            help="file to convert from")
    parser.add_argument('outfile', type=str, help="file to convert to")
    parser.add_argument('--oligo', type = str, help = "oligo to calculate against")
    parser.add_argument('--window', type = int, default = 6, help = "size of window to slide over genome. Default = 6")
    parser.add_argument('--Na', type = int, default = 150, help = "Concentration of Sodium in mM. Default = 150")
    parser.add_argument('--Tris', type = int, default = 0, help = "Concentration of Tris in mM. Default = 0")
    parser.add_argument('--Mg', type = int, default = 0, help = "Concentration of Mg in mM. Default = 0")
    parser.add_argument('--dNTPs', type = int, default = 0, help = "Concentration of dNTPs in mM. Default = 0")
    parser.add_argument('--nuc_type', type = str, default = "DNA:DNA", help = "Type of nucleic acid interaction. Default = DNA:DNA; Possible: RNA:RNA RNA:DNA")
    parser.add_argument('--conc1', type = int, default = 25, help = "Concentration of short oligo. Default = 25 nM")
    parser.add_argument('--conc2', type = int, default = 25, help = "Concentration of binding site. Default = 25 nM")
    parser.add_argument('--strand', type = str, default = "+", help = "Which strand of the nucleic acid? Default = '+'")
    parser.add_argument('--out_type', type = str, default = ".tsv", help = "Which output type? .tsv or .bw. Default = .tsv")

    args = parser.parse_args() 

    ffile = fa.FastaFile()
    with open(args.infile, mode = "r") as inf:
        ffile.read_whole_file(inf)

    # function factory for type of interaction
    ab_nuc = {"DNA:DNA" : mt.DNA_NN3, "RNA:DNA" : mt.R_DNA_NN1, "RNA:RNA": mt.RNA_NN3}

    if args.oligo is not None:
        tm_args = {"nn_table": ab_nuc[args.nuc_type], "Na": args.Na, "Tris": args.Tris, \
                "Mg": args.Mg, "dNTPs": args.dNTPs, \
                "dnac1": args.conc1, "dnac2": args.conc2}

        to_tsv_oligo(generate_windowed_oligo(ffile, args.oligo,args.strand, None),args.outfile, args.strand)
    else:
        if args.out_type == ".tsv":
            to_tsv(generate_windowed_tm(ffile, args.window, \
                    nn_table = ab_nuc[args.nuc_type], \
                    Na = args.Na, Tris = args.Tris, \
                    Mg = args.Mg, strand = args.strand, dNTPs = args.dNTPs, \
                    dnac1 = args.conc1, dnac2 = args.conc2), args.outfile, args.strand)
        else:
            to_bw(generate_windowed_tm(ffile, args.window, \
                    nn_table = ab_nuc[args.nuc_type], \
                    Na = args.Na, Tris = args.Tris, \
                    Mg = args.Mg, strand = args.strand, dNTPs = args.dNTPs, \
                    dnac1 = args.conc1, dnac2 = args.conc2), args.window, ffile, args.outfile)

    


