#!/usr/bin/python
def U00096_2_to_ATCC_47076(loc):
    if loc > 547832:
        loc = loc +  1
    return loc

def ATCC_47076_to_U00096_2(loc):
    if loc > 547832:
        loc = loc -  1
    return loc

def U00096_2_to_U00096_3(loc):
    if loc > 257899 and loc <= 547832:
        # adding the 777 basepair IS1 insertion at U00096.2
        # coordinate 257899
        loc = loc + 776
    elif loc > 547832 and loc <= 1298719:
        # adding the single basepair insertion at U00096.2
        # coordinate 547832
        loc = loc + 776 + 1
    elif loc > 1298719 and loc <= 2171386:
        # adding the 1200 bp IS5 insertion at U00096.2
        # coordinate 1298719
        loc = loc + 776 + 1 + 1199
    elif loc > 2171386 and loc < 3558478:
        # adding the 2 bp insertion at U00096.2
        # coordinate 2171386
        loc = loc + 776 + 1 + 1199 + 2
    elif loc >=  3558478:
        # adding the 1 bp deletion at U00096.2
        # coordinate 3558478
        loc = loc + 776 + 1 + 1199 + 2 - 1
    return loc

def U00096_3_to_U00096_2(loc):
    if loc > 257899 and loc <= 258675:
        # if within the IS1 insertion, convert to the basepair before
        loc = 257899
    elif loc > 258675 and loc < 548609:
        # after the IS1 insertion, just subtract the insertion up to the
        # next insertion
        loc = loc - 776
    elif loc >= 548609 and loc <= 1299496:
        # if at the single bp insertion or after, subtract the insertion and the
        # IS1 insertion
        loc = loc - 776 - 1
    elif loc > 1299496 and loc <= 1300696:
        # if within the IS5 insertion convert to the bp before
        loc = 1298719
    elif loc > 1300696 and loc <= 2173362:
        # if after the IS5 insertion, subtract the IS1, IS5 and 1 bp insertion
        loc = loc - 776 - 1 - 1199
    elif loc > 2173362 and loc <= 2173364:
        # if within the 2 bp insertion, act like you are the bp before
        loc = 2171386
    elif loc > 2173364 and loc < 3560456:
        loc = loc - 776 - 1 - 1199 - 2
    elif loc >=  3560456:
        loc = loc - 776 - 1 - 1199 - 2 + 1
    return loc

def ATCC_47076_to_U00096_3(loc):
    return(U00096_2_to_U00096_3(ATCC_47076_to_U00096_2(loc)))

def U00096_3_to_ATCC_47076(loc):
    return(U00096_2_to_ATCC_47076(U00096_3_to_U00096_2(loc)))

if __name__ == "__main__":
    import argparse
    import sys
    import gfftools as gff

    parser = argparse.ArgumentParser(description="Convert genomic coordinates")
    parser.add_argument('number', type=str, help="coordinate to convert, if - uses stdin")
    parser.add_argument('in_genome', type=str, help="input genome [U00096_2, U00096_3, ATCC_47076]")
    parser.add_argument('out_genome', type=str, help="output genome [U00096_2, U00096_3, ATCC_47076]")
    parser.add_argument('--gff_out', type=str, help="treat number as gff file and output to gff specified")
    parser.add_argument('--bed_out', type=str, help="treat number as bed file and output to bed specified")

    func_dict = {"U00096_2_to_ATCC_47076":U00096_2_to_ATCC_47076,
                 "U00096_3_to_ATCC_47076":U00096_3_to_ATCC_47076,
                 "U00096_3_to_U00096_2":U00096_3_to_U00096_2,
                 "ATCC_47076_to_U00096_2":ATCC_47076_to_U00096_2,
                 "ATCC_47076_to_U00096_3":ATCC_47076_to_U00096_3,
                 "U00096_2_to_U00096_3":U00096_2_to_U00096_3}
    out_chrom_name = {"U00096_3": "U00096.3",
                      "ATCC_47076": "ATCC.47076",
                      "U00096_2":"U00096.2"}
    args = parser.parse_args()
    conversion_func = func_dict[args.in_genome + "_to_" +args.out_genome]

    if args.number == "-":
        for value in sys.stdin:
            value = int(float(value.rstrip()))
            sys.stdout.write(str(conversion_func(value))+"\n")
    elif args.gff_out:
        input_gff = gff.GffData()
        input_gff.parse_gff_file(args.number)
        for line in input_gff.data:
            line.genome_name = out_chrom_name[args.out_genome]
            line.start = conversion_func(line.start)
            line.end = conversion_func(line.end)
        input_gff.write_gff_file(args.gff_out)
    elif args.bed_out:
        with open(args.number, mode="r") as in_f:
            with open(args.bed_out, mode="w") as out_f:
                for line in in_f:
                    linearr = line.rstrip().split("\t")
                    linearr[0] = out_chrom_name[args.out_genome]
                    linearr[1] = str(conversion_func(int(linearr[1])))
                    linearr[2] = str(conversion_func(int(linearr[2])))
                    out_f.write("\t".join(linearr)+"\n")
    else:
        sys.stdout.write(str(conversion_func(int(args.number))) + "\n")



    
    
