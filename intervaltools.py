import numpy as np
import bed_utils
import gfftools

def overlap(int1, int2):
    return int1[0] <= int2[1] and int2[0] <= int1[1]

class Interval(object):

    def __init__(self, start=0, end=0, value=None):

        self.start = start
        self.end = end
        self.value = value

    def __repr__(self):
        return "<Interval object %s; %s-%s; %s>"%(hex(id(self)), self.start, self.end, self.value)

    def overlap(self, other):
        return overlap((self.start, self.end), (other.start, other.end))


class Intervals(object):

    def __init__(self):
        self.data = []

    def __iter__(self):
        for entry in self.data:
            yield entry

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return "<Intervals object %s; %s ints; %s - %s>"%(id(self), len(self), self.data[0], self.data[-1])

    def add_interval(self, interval):
        self.data.append(interval)

    def clean_up(self):
        self.data = list(set(self.data))
        self.data.sort(cmp=lambda a,b: cmp(a.start, b.start))

    def check_overlap(self, ointerval):
        for sinterval in self:
            if sinterval.overlap(ointerval):
                return True
        return False

    def to_array(self, array=None,  array_length=None):
        if array is None and array_length is None:
            raise ValueError("Need array length to convert to array")
        if array is None:
            array = np.zeros(array_length, dtype=bool)
        for entry in self:
            array[entry.start:entry.end] = True
        return array

    def random(self, array_length, no_overlap=False, num=1, bad_locs=None):
        new_intervals = Intervals()
        already_used = np.zeros(array_length, dtype=bool)
        if bad_locs:
            already_used = np.logical_or(already_used, bad_locs)
        for entry in self:
            entry_length = entry.end-entry.start
            found = 0
            while found < num:
                newstart = np.random.randint(0, high=array_length-entry_length)
                if no_overlap and np.sum(already_used[newstart:newstart+entry_length]) > 0:
                    continue
                else:
                    new_intervals.add_interval(Interval(newstart, newstart+entry_length))
                    if no_overlap:
                        already_used[newstart:newstart+entry_length] = True
                    found += 1
        return new_intervals

    def probes_overlap(self, other, genome_length):
        self_all_locs = self.to_array(genome_length)
        other_all_locs = self.to_array(genome_length)
        return np.sum(np.logical_and(self_all_locs, other_all_locs))

    def overlapping_intervals(self, other):
        self.clean_up()
        other.clean_up()
        new_self = Intervals()
        new_other = Intervals()
        for sint in self:
            for oint in other:
                if sint.overlap(oint):
                    new_self.add_interval(sint)
                    new_other.add_interval(oint)
        new_self.clean_up()
        new_other.clean_up()
        return new_self, new_other

    def from_other_ftype(self, otherobject, func):
        for entry in otherobject:
            start, end, value = func(entry)
            self.add_interval(Interval(start, end, value))

def array_overlap(array1, array2):
    return np.sum(np.logical_and(array1, array2))

def convert_file_to_intervals(fname, function_factory):
    new_interval = Intervals()
    if fname.endswith(".bed") or fname.endswith(".narrowPeak"):
        this_file =bed_utils.BedFile()
        this_file.from_bed_file(fname)
        new_interval.from_other_ftype(this_file, function_factory[".bed"])

    elif fname.endswith(".gff"):
        this_file =gfftools.GffData()
        this_file.parse_gff_file(fname)
        new_interval.from_other_ftype(this_file, function_factory[".gff"])
    else:
        raise ValueError("%s filetype not supported yet"%fname)
    return new_interval

def sampled(intervals1, intervals2, length = 4639675, n=2000):

    samples_as_extreme = 0
    array1 = intervals1.to_array(array_length=length)
    array2 = intervals2.to_array(array_length=length)
    array1size = (np.sum(array1) + 0.0)
    array2size = (np.sum(array2)+0.0)
    actual = (array_overlap(array1, array2)/array2size)/(array1size/length)
    for i in range(n):
        this_random = intervals1.random(array_length=length)
        random_array = this_random.to_array(array_length=length)
        r_arraysize = (np.sum(random_array)+0.0)
        sampled_val = (array_overlap(random_array, array2)/array2size)/(r_arraysize/length)
        if sampled_val >= actual:
            samples_as_extreme += 1
        if i % 100 == 0:
            print("On sample %i, actual: %s, random: %s %s"%(i, actual, sampled_val, np.sum(random_array)))
    return samples_as_extreme /(n+0.0)

if __name__ == "__main__":
    import sys

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    n = int(sys.argv[3])
    genome_length = 4639675

    possible_ftypes = {".bed": lambda entry: (entry["start"], entry["end"], entry["score"]), 
            ".gff": lambda entry: (entry.start-1, entry.end, ".")}

    intervals1 = convert_file_to_intervals(file1, possible_ftypes)
    intervals2 = convert_file_to_intervals(file2, possible_ftypes)
    overlaps1, overlaps2 = intervals1.overlapping_intervals(intervals2)
    print("File 1: %s"%file1)
    print("File 2: %s"%file2)
    print("For each below: number in interval; number overlapping; fraction overlapping")
    print("Overlap 1: %s %s %s"%(len(intervals1), len(overlaps1), (len(overlaps1)+0.0)/len(intervals1)))
    print("Overlap 2: %s %s %s"%(len(intervals2), len(overlaps2), (len(overlaps2) + 0.0)/len(intervals2)))
    array1 = intervals1.to_array(array_length=genome_length)
    array1size = (np.sum(array1) + 0.0)
    array2 = intervals2.to_array(array_length=genome_length)
    array2size = (np.sum(array2) + 0.0)
    overlap = array_overlap(array1, array2)
    print("Probes Overlap 1: %s %s %s"%(int(np.sum(array1)), int(overlap), (overlap+0.0)/np.sum(array1)))
    print("Probes Overlap 2: %s %s %s"%(int(np.sum(array2)), int(overlap), (overlap+0.0)/np.sum(array2)))
    print("Summary Stat: %s"%((overlap/array2size)/(array1size/genome_length)))
    print("Pvalue %s for %s random samples of File 1"%(sampled(intervals1, intervals2, n=n), n))




