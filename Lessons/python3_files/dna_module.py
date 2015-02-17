"""
    This module contains functions for computing various DNA sequence statistics.
"""

def compute_pairwise_similarity(seq1, seq2, skip_gaps = True):
    """
        Compute the pairwise similarity between two sequences.
        Positional arguments:
            1. seq1 and seq2, the two sequence strings to compare. These sequences must be the same length!!
        Keyword arguments:
            1. skip_gaps: whether positions which are a gap in one sequence but not the other should be skipped (default: True)
    """
    # assert(len(seq1) == len(seq2) )
    same = 0.
    total = 0.
    for i in range(len(seq1)):
        if skip_gaps and (seq1[i] == "-" or seq2[i] == "-"):
            continue
        else:
            total += 1.
            if seq1[i] == seq2[i]:
                same += 1.
    
    return same / total
            
                
        
    


def compute_at_content(seq, digits = 5, percent = False):
    """ 
        Compute AT content from a provided sequence. By default, result is rounded to 5 significant digits.
        Positional arguments: 
            1. seq:  DNA sequence from which GC-content should be calculated
        Keyword arguments:
            1. digits: number of significant digits for results (default: 5)
            2. percent: return value as a percent rather than decimal (default: False)
    """
    a = seq.lower().count("a")
    t = seq.lower().count("t")
    at = round( float(a + t) / length(seq), digits)
    if percent:
        return at*100
    else:
        return at
        
        
        
def compute_gc_content(seq, digits = 5, percent = False):
    """ 
        Compute GC content from a provided sequence. By default, result is rounded to 5 significant digits.
        Positional arguments: 
            1. seq:  DNA sequence from which GC-content should be calculated
        Keyword arguments:
            1. digits: number of significant digits for results (default: 5)
            2. percent: return value as a percent rather than decimal (default: False)
    """
    g = seq.upper().count("G")
    c = seq.upper().count("C")
    gc = float(g + c) / length(seq)
    if percent:
        return round(gc * 100, digits)
    else:
        return round(gc, digits)
        

    
    
    
    
    