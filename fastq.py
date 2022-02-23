"""
fastq module to manipulate fastq files
Written by Michael Wolfe
"""

import gzip
from itertools import islice

def check_unique(fastq_file):
    name_dict = {}
    count_dict = {}
    this_file = FastqFile(fastq_file)
    for entry in this_file.stream_file():
        name_count = name_dict.get(entry.name, 0)
        name_dict[entry.name] = name_count + 1 
    for name, name_count in name_dict.items():
        histogram_count = count_dict.get(name_count, 0)
        count_dict[name_count] = histogram_count + 1
    print("Name_count\tHistogram_count")
    for key in sorted(count_dict):
        print("%s\t%s"%(key, count_dict[key]))

        




class FastqEntry(object):
    """
    Stores information for a single fastq entry

    An example of a fastq entry is below

    @SEQ_ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    +
    !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
    """
    def __init__(self, name = "", seq = "", opt = "", qual = ""):
        self.name = name
        self.seq = seq
        self.opt = opt
        self.qual = qual

    def __str__(self):
        return "%s\n%s\n%s\n%s\n"%(self.name, self.seq, self.opt, self.qual)

class FastqFile(object):
    def __init__(self, infile = ""):
        self.infile = infile
        self.names = []
        self.data = {}

    def read_whole_file(self):
        if type(self.infile) is str and self.infile.endswith(".gz"):
            fhandle = gzip.open(self.infile, mode = "rt")
        else:
            fhandle = open(self.infile, mode = "r")
        # assumes well-formed fastq
        while True:
            entry = [line.rstrip() for line in islice(fhandle,4)]
            if entry:
                name, seq, opt, qual = entry
                self.data[name] = FastqEntry(name, seq, opt, qual)
                self.names.append(name)
            else:
                break
        fhandle.close()
        if len(self.names()) != len(self.data):
            raise ValueError("Duplicate read names found in file")

    def stream_file(self):
        if type(self.infile) is str and self.infile.endswith(".gz"):
            fhandle = gzip.open(self.infile, mode = "rt")
        else:
            fhandle = open(self.infile, mode = "r")
        # this is a terrible way to do this, but I can't think of a better way
        while True:
            entry = [line.rstrip() for line in islice(fhandle,4)]
            if entry:
                name, seq, opt, qual = entry
                yield FastqEntry(name, seq, opt, qual)
            else:
                break

    def write(self, outfile):
        if type(outfile) is str and outfile.endswith(".gz"):
            fhandle = gzip.open(outfile, mode = "wb")
            for entry in self:
                fhandle.write(str(entry).encode())
        else:
            fhandle = open(outfile, mode = "w")
            for entry in self:
                fhandle.write(str(entry))

    def __iter__(self):
        for name in self.names:
            yield self.data[name] 

    def __len__(self):
        return len(self.names)

    def names(self):
        return self.names

    def pull_entry(self, name):
        return self.data[name]

    def add_entry(self, entry):
        self.data[entry.name] = entry
        self.names.append(entry.name)

    def subset_by_names(self, names):
        out = FastqFile()
        for name in names:
            out.add_entry(self.data[entry.name])
        return out

if __name__ == "__main__":
    import sys
    infile = sys.argv[1]
    if infile == "0":
        infile = int(infile)

    outfile = sys.argv[2]
    if outfile == "1":
        outfile = int(outfile)


    this_file = FastqFile(infile)
    if len(sys.argv) > 3:
        check_unique(infile)
    else:
        for entry in this_file.stream_file():
            print(entry)
