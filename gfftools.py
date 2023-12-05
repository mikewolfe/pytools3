#!usr/bin/python

# Tools for managing gff files
# These are really just files that contain a set of 9 tab delimited fields:
#  Genome_name Data_origin Site_type start end . +/-(direction) . comments

# additional functions for splitting the last field of a gff entry
def newSplit(value):
    import shlex
    lex = shlex.shlex(value)
    lex.quotes = '"'
    lex.whitespace_split = True
    lex.commenters = ''
    return list(lex)

def make_comment_dict(gff_entry):
    gff_entry.comment_dict = {}
    keyvalue = gff_entry.comments.split(";")
    for pair in keyvalue:
        key_value= newSplit(pair)
        if len(key_value) < 2:
            value = ""
        else:
            key = key_value[0]
            value = "".join(key_value[1:])
            value = value.replace('"', "")
        gff_entry.comment_dict[key] = value

def overlap(x1, x2, y1, y2):
    return(x1 <= y2 and y1 <= x2)

class GffEntry:
  """
  A simple container class for the equivalent of a gff file line
  """

  FORMAT_STRING = "%s\t%s\t%s\t%i\t%i\t.\t%s\t.\t%s"

  def __init__(self, line=None):

    if (line):
      self.parse_gff_line(line)
    else:
      self.genome_name = ""
      self.data_origin = ""
      self.site_type = ""
      self.start = 0
      self.end = 0
      self.direction = "+"
      self.comments = ""


  def parse_gff_line(self, line):
    """
    Set this entry's values to those of a line from a gff file
    """

    datarray = line.split("\t")
    self.genome_name = datarray[0]
    self.data_origin = datarray[1]
    self.site_type = datarray[2]
    self.start = int(datarray[3])
    self.end = int(datarray[4])
    self.direction = datarray[6]
    self.comments = " ".join(datarray[8:])
    self.comments = self.comments.replace("\t", " ")
    self.comments = self.comments.replace("\n", "")

  def __repr__(self):
    """
    Return a formatted gff line, which can be used to reconstitute the object or be written directly to a gff file
    """

    return GffEntry.FORMAT_STRING % (self.genome_name, self.data_origin, self.site_type, self.start, self.end, self.direction, self.comments)



class GffData:
  """
  Class for storing and manipulating gff data
  """

  def __init__(self):
    self.data = []
    self.index = 0

  def __iter__(self):
    return self

  def __next__(self):
    if self.index == len(self.data):
      self.index = 0
      raise StopIteration
    self.index = self.index + 1
    return self.data[self.index-1]



  def clear_db(self):
    self.data = []

  def parse_gff_file(self,filename, clear=True):
    """
    Parse a gff file and store the lines in self.data

    The current contents of this object are overwritten iff clear 
    """

    if (clear):
      self.clear_db()

    instr = open(filename, "r")
    for line in instr:
      newline = GffEntry(line)
      self.data.append(newline)

  def cleanup(self):
    """
    Remove all duplicate entries and sort based on starting position
    """

    self.data = list(set(self.data))
    self.data.sort(cmp=lambda a,b: cmp(a.start, b.start))

  def write_gff_file(self, filename):
    """
    Write the current contents of my data to a file
    """

    ostr = open(filename, "w")

    for line in self:
      ostr.write("%s\n" % line)

    ostr.close()

  def addline(self, genome_name, data_origin, site_type, start, end, direction, comments):
    """
    Add a line with the given data
    """

    newobj = GffEntry()
    newobj.genome_name = genome_name
    newobj.data_origin = data_origin
    newobj.site_type = site_type
    newobj.start = start
    newobj.end = end
    newobj.direction = direction
    newobj.comments = comments

    self.data.append(newobj)

  def to_bedfile(self, filename, ID_field = "ID"):

    ostr = open(filename, "w")
    #chrom,start,end,name,score,strand
    outfmt="%s\t%i\t%i\t%s\t%s\t%s\n"
    for line in self:
      comment_dict = make_comment_dict(line)
      this_line=outfmt%(line.genome_name, line.start-1, line.end, line.comment_dict[ID_field],".",line.direction)
      ostr.write(this_line)

    ostr.close()




  def find_entry(self, findfunc, findall=False):
    """
    Return the line or lines for which findfunc is true when given a gff item

    If findall is false, only the first such entry is returned

    Findfunc should take a gffline object as its only argument
    """

    matches = filter(findfunc, self.data)

    if len(matches) == 0:
      return []

    if (findall):
      return matches
    else:
      return matches[0]

  def add_entry(self, new_entry):
    # Add an externally constructed GffEntry object

    self.data.append(new_entry)

if __name__ == "__main__":
    # convert gff to bed file
    import sys
    infile = sys.argv[1]
    outfile = sys.argv[2]
    ID_field = sys.argv[3]
    ingff = GffData()
    ingff.parse_gff_file(infile)
    ingff.to_bedfile(outfile, ID_field)
