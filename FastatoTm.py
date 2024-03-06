import fasta as fa
from Bio.SeqUtils import MeltingTemp as mt

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
    parser = argparse.ArgumentParser("Tools to work with BigWig files") 
    parser.add_argument('infile', type=str, 
            help="file to convert from")
    parser.add_argument('outfile', type=str, help="file to convert to")
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

    


