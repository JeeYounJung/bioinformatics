def pairwise_edit_distances(seqs, sub_matrix):
    """Computes a pairwise distance matrix between a set of sequences
       by computing a global alignment (with linear gap penalty)
       between each pair and then computing the edit distances implied 
       by the resuling alignments.
    
    Args:
        seqs: a list of sequences (strings)
        sub_matrix: a substitution matrix, represented as a dictionary.
                    the matrix should include space scorees.
    Returns:
       A symmetric distance matrix represented as a list of lists.
    """    
    n = len(seqs)
    d = [[0] * n for i in range(n)]
    for j in range(n):
        for i in range(j):
            d[i][j] = d[j][i] = num_edits(seqs[i], seqs[j], sub_matrix)
    return d

def num_edits(seq1, seq2, sub_matrix):
    """Computes the edit distance between two sequences after aligning
       them globally with a linear gap penalty and given substitution matrix.
    Args:
        seq1: the first sequence (a string)
        seq2: the second sequence (a string)
        sub_matrix: a substitution matrix, represented as a dictionary.
                    the matrix should include space scorees.
    Returns:
       The edit distance implied by the alignment of the two sequences.
    """
    alignment = align_global([seq1], [seq2], sub_matrix)
    return sum(map(num_edits_in_column, transpose_alignment(alignment)))

def num_edits_in_column(column):
    """Computes the total number of pairwise edits (substitutions or insertions/deletions)
       within a column of a pairwise or multiple alignment.
    Args:
        column: a column of an alignment represented as a string
    """
    return sum(column[i] != column[j] for j in range(len(column)) for i in range(j))

def align_global(x, y, substitution_matrix):
    """Computes a global pairwise alignment of sequences/alignments x and y
    
    Uses a linear gap scoring function with the space score encoded in the given substitution matrix.
    In the case of multiple optimal alignments, this function returns a "middle-road" 
    alignment (prefers match/mismatch first).
    
    Args:
        x: an alignment/sequence represented as a list of strings. A single 
           unaligned sequence should be represented by a list with a single element (the sequence).
        y: an alignment/sequence represented as a list of strings. A single 
           unaligned sequence should be represented by a list with a single element (the sequence).
        substitution_matrix: the substitution matrix, represented as a dictionary 
    Returns:
        An alignment, represented as a list of aligned strings of the same length.
    """
    # create lists of the columns of the two alignments
    xt, yt = transpose_alignment(x), transpose_alignment(y)
    num_rows, num_cols = len(xt) + 1, len(yt) + 1
    x_all_spaces, y_all_spaces = '-' * len(x), '-' * len(y)
    m = matrix(num_rows, num_cols) # score matrix
    t = matrix(num_rows, num_cols) # traceback pointer matrix

    # constants representing the traceback pointers, with higher
    # values representing higher priority for traceback in case of ties
    M, IX, IY = 3, 2, 1
    
    # A function that scores the alignment of two columns
    def score_column_pair(x_column, y_column):
        return sum(substitution_matrix[a, b] for a in x_column for b in y_column)

    # initialization
    m[0][0] = 0
    for i in range(1, num_rows):
        m[i][0] = m[i - 1][0] + score_column_pair(xt[i - 1], y_all_spaces)
        t[i][0] = IX
    for j in range(1, num_cols):
        m[0][j] = m[0][j - 1] + score_column_pair(x_all_spaces, yt[j - 1])
        t[0][j] = IY
  
    # main fill
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            m[i][j], t[i][j] = max((m[i - 1][j - 1] + score_column_pair(xt[i - 1], yt[j - 1]),     M),
                                   (m[i - 1][j]     + score_column_pair(xt[i - 1], y_all_spaces), IX),
                                   (m[i][j - 1]     + score_column_pair(x_all_spaces, yt[j - 1]), IY))

    # traceback
    columns = []
    i, j = num_rows - 1, num_cols - 1
    while i > 0 or j > 0:
        traceback = t[i][j]
        if traceback == M:
            i -= 1
            j -= 1            
            columns.append(xt[i] + yt[j])
        elif traceback == IX:
            i -= 1
            columns.append(xt[i] + y_all_spaces)
        elif traceback == IY:
            j -= 1
            columns.append(x_all_spaces + yt[j])
        else:
            assert False
    
    return transpose_alignment(reversed(columns)) if columns else [""] * (len(x) + len(y))

def matrix(num_rows, num_cols, initial_value=None):
    """Returns a matrix (a list of rows, each of which is a list) 
    with num_rows and num_cols and with initial_value in each entry"""
    return [[initial_value] * num_cols for i in range(num_rows)]

def transpose_alignment(alignment):
    """Returns a column-based alignment from a row-based alignment or vice versa"""
    return list(map(''.join, zip(*alignment)))


