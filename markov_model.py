import numpy as np
from itertools import product

def initialize_matrix(alphabet, order):
    """
    Helper function to intialize an empty transition matrix based on an alphabet
    and order of the Markov Model

    Arguments:
        alphabet (list of char) - possible sequence values
        order (int) - order of the markov model (i.e 3 indicates triplets)
        
    Returns:
        matrix (np.array, float) - zeros of size needed for model
        rownames (np.array, unicode) - names for the rows
        colnames (np.array, unicode) - names for the columns
    """
    rownames = []
    colnames = []
    # use itertools cartesian product to get all possible combinations of a
    # given order
    for i, val in enumerate(product(alphabet, repeat=order)):
        rownames.append("".join(val))
    for j, name in enumerate(alphabet):
        colnames.append(name)
    rownames = np.array(rownames)
    colnames = np.array(colnames)
    matrix = np.zeros((len(rownames), len(colnames)))
    return (matrix, rownames, colnames)

class TransitionMatrix(object):
    """
    Object to store a Transition Matrix for a markov model with appopriate
    methods to train and sample from the matrix
    """

    def __init__(self, alphabet=["A","G","T","C"], order=3, prior=1):
        """
        Class attributes:
        alphabet (list) - possible sequence values
        order (int) - order of the transition matrix
        prior (int or matrix of appropriate size) - prior to add to matrix
                                                    counts. No error checking is 
                                                    done on matrix size
        matrix (np.array, float) - numpy matrix containing transition matrix
        rownames (np.array, unicode) - rownames which equal all possible perm
                                       of order self.order
        colnames (np.array, unicode) - values same as self.alphabet
        colname_map (dict) - mapper from name to index
        rowname_map (dict) - mapper from name to index
        """
        # alphabet
        self.alpha = alphabet
        # order
        self.order = order
        self.prior = prior
        matrix, rownames, colnames = initialize_matrix(self.alpha, self.order)
        # transition matrix
        self.matrix = matrix
        # rownames
        self.rownames = rownames
        # colnames
        self.colnames = colnames
        self.rowname_map = {}
        # make a dictionary to map rownames to corresponding index in matrix
        for i, val in enumerate(self.rownames):
            self.rowname_map[val] = i
        # ditto for colnames
        self.colname_map = {}
        for i, val in enumerate(self.colnames):
            self.colname_map[val] = i

    def index_by_row(self, key):
        """
        Accessor for mapping rowname to index
        """
        return self.rowname_map[key]

    def index_by_col(self, key):
        """
        Accessor for mapping colname to index
        """
        return self.colname_map[key]

    def increment_val(self, rowkey, colkey, val = 1):
        """
        Increment a value in the transition matrix by val
        """
        row = self.index_by_row(rowkey)
        col = self.index_by_col(colkey)
        self.matrix[row,col] += val

    def normalize(self):
        """
        Normalize the transtion matrix from a count matrix to a frequency
        matrix. Prior is added from self.prior before normalization.
        This operation is irreversible and loses information
        """
        self.matrix += self.prior
        for i, row in enumerate(self.matrix):
            self.matrix[i,] = row/np.sum(row)

    def set_seed(self, seed):
        """
        TO DO use this function in parallelization. 
        Set the numpy random seed and store what it was set to in the class
        attributes. Might be needed if used in parallel. 
        """
        self.seed = seed
        np.random.seed(seed)

    def emit(self, instr):
        """
        Emit a single alphabet value using the transition matrix. Transition matrix
        MUST be normalized before using emit, otherwise nonsense answers may result.
        Error checking for normalization should be done outside this function
        """
        outstr = np.random.choice(self.colnames, p=self.matrix[self.index_by_row(instr),:])
        return(outstr)


class MarkovModel(object):
    """
    Object to store a Markov Model for both training and generation of new
    sequences
    """

    def __init__(self, alphabet=["A", "G", "T", "C"], order=3, prior=1):
        """
        Args:
        alphabet - possible sequence values passed to TransitionMatrix
        order - order of model, passed to TransitionMatrix
        prior - prior for model, passed to Transition Matrix
        Attributes:
        self.tm (TransitionMatrix) - holds transition matrix needed for the model
        self.normalize (bool) - has the matrix been normalized?
        """
        self.tm = TransitionMatrix(alphabet, order, prior)
        self.normalized = False

    def train(self, seqs, normalize=True):
        """
        Function to train a MarkovModel based off a set of sequences

        Arguments:
            seqs (list of char) - list of sequences to train on
            normalize (bool) - normalize the transition matrix after training?
        """

        wind_size = self.tm.order
        for seq in seqs:
            seq_len = len(seq)
            if seq_len < wind_size + 1:
                continue
            for i in range(0, seq_len-wind_size-1):
                self.tm.increment_val(seq[i:i+wind_size], seq[i+wind_size])
        if normalize:
            self.normalize()

    def normalize(self):
        if not self.normalized:
            self.tm.normalize()
            self.normalized = True

    def generate(self, size, seed):
        """
        Function to generate a new sequence based off a markov model

        Arguments:
            size (int) - size of sequence to create
            seed (char) - a character seed of size self.tm.order to start 
                          sequence generation
        Returns:
            generated sequence of length size
        """
        if size < self.tm.order + 1:
            raise ValueError("Size of sequence to generate must be larger than order + 1 size: %s order: %s"%(size, self.tm.order))
        new_seq = [val for val in seed]
        for i in range(0,size-self.tm.order):
            new_seq.append(self.tm.emit(seed))
            seed = "".join(new_seq[-self.tm.order:])
        return "".join(new_seq)

    def generate_using_seqs(self, seqs):
        """
        Function to generate a new sequences based off a markov model and 
        a set of sequences. Returns an empty sequence for any sequences not
        long enough.

        Arguments:
            seqs (list of char) - list of sequences to generate new seqs
                                  from
        Returns:
            list of generated sequences each of the same length as the original
            sequences. Uses as a seed the first self.tm.order characters of 
            each original sequence
        """
        final_seqs = []
        for seq in seqs:
            seq_len = len(seq)
            if seq_len > self.tm.order:
                seq_seed = seq[0:self.tm.order]
                yield self.generate(seq_len, seq_seed)
            else:
                yield ""
