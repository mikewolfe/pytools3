from Bio import SeqIO
import sys
import argparse

def convert_strandval(strandval):
    if strandval == 1:
        out = "+"
    elif strandval == -1:
        out = "-"
    else:
        out = "NA"
    return out

def out_tsv(gb):
    sys.stdout.write("chr\tstart\tend\tstrand\tlocus_tag\tprotein_id\tproduct\tseq\ttranslation\n")
    for feature in gb.features:
        outstr = ""
        if feature.type == "CDS":
            outstr += str(gb.name) + "\t"
            outstr += str(feature.location.start) + "\t"
            outstr += str(feature.location.end) + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\t"
            outstr += str(feature.qualifiers.get("locus_tag", ["NA"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("protein_id", ["NA"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("product", ["NA"])[0]) + "\t"
            # pull the nucleotide sequence
            nuc_seq = feature.extract(gb.seq)
            this_translation = nuc_seq.translate(table = "Bacterial", cds = True)
            translation = str(feature.qualifiers['translation'][0])
            # check that the translation matches the expected product. Have to drop the stop codon in
            # the nucleotide translation
            if not (this_translation == translation):
                raise ValueError("Nucleotide sequence doesn't match translation\n%s\n%s"%(this_translation, translation))
            outstr += str(nuc_seq) + "\t"
            outstr += translation + "\n"
            sys.stdout.write(outstr)


def out_bed(gb, features = ["CDS"]):
    for feature in gb.features:
        if feature.type in "CDS":
            outstr = ""
            outstr += str(gb.name) + "\t"
            outstr += str(feature.location.start) + "\t"
            outstr += str(feature.location.end) + "\t"
            outstr += str(feature.qualifiers.get("locus_tag", ["NA"])[0]) + "\t"
            outstr += "." + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\n"
            sys.stdout.write(outstr)

def out_fasta(gb, features = ["CDS"]):
    import fasta as fa
    out_fasta = fa.FastaFile()
    for feature in gb.features:
        if feature.type in features:
            header = ">" + str(feature.qualifiers.get("locus_tag", ["NA"])[0]) + \
            " " +  \
            str(feature.qualifiers.get("product", ["NA"])[0])
            seq = str(feature.extract(gb.seq))
            entry = fa.FastaEntry(header = header, seq = seq)
            out_fasta.add_entry(entry)
    out_fasta.write(sys.stdout)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a genbank file for features")
    parser.add_argument("infile", type=str, help="path to input file")
    parser.add_argument("--outfmt", type=str, help='''output format. default = bed. 
                                                    Supported: [fa, fasta, bed, tsv]
                                                    Note that tsv is a special format
                                                    that only pulls CDS's and their translation''',
                        default = "bed")
    parser.add_argument('--features', type=str, nargs="+", help='''
    feature types to parse for. Can specify multiple. default = "CDS"
    ''',
            default = "CDS")
    args = parser.parse_args()

    fname = args.infile
    ftype = args.outfmt
    features = args.features
    with open(fname, mode = "r") as inf:
        gb = SeqIO.read(inf, "genbank")
    if ftype == "fasta" or ftype == "fa":
        out_fasta(gb, features)
    elif ftype == "bed":
        out_bed(gb, features)
    else:
        out_tsv(gb)

