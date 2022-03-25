# NW.py
# HW1, Computational Genomics, Spring 2022
# andrewid: nmehtani

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np
import math
def ReadFASTA(filename):
    fp=open(filename, 'r')
    Sequences={}
    tmpname=""
    tmpseq=""
    for line in fp:
        if line[0]==">":
            if len(tmpseq)!=0:
                Sequences[tmpname]=tmpseq
            tmpname=line.strip().split()[0][1:]
            tmpseq=""
        else:
            tmpseq+=line.strip()
    Sequences[tmpname]=tmpseq
    fp.close()
    return Sequences

def create_matrices(numRows, numCols, gap):
    
    ## Create Array
    alignMtx = np.ndarray(shape=(numRows, numCols))
    alignMtx.fill(-math.inf)
    
    tracebackMtx = np.ndarray(shape=(numRows, numCols), dtype=object)

    ## Initialize
    alignMtx[0, 0] = 0
    
    for i in range(1, numRows):
        alignMtx[i, 0] = gap*i
        tracebackMtx[i, 0] = 2
    
    for j in range(1, numCols):
        alignMtx[0, j] = gap*j
        tracebackMtx[0, j] = 1
        
#     print(alignMtx)
    return alignMtx, tracebackMtx

def traceback_alignments(seq1, seq2, tracebackMtx):

    n1 = len(seq1)
    n2 = len(seq2)
    
    i, j = n1, n2
    currentCell = tracebackMtx[i, j]
    
    alignment1, alignment2 = '', ''
    while i != 0 or j != 0:
        if currentCell == 0:
            
            alignment1 += seq1[i - 1]
            alignment2 += seq2[j - 1]
            
            i -= 1
            j -= 1
            currentCell = tracebackMtx[i, j]

        elif currentCell == 1:
            
            alignment1 += '-'
            alignment2 += seq2[j - 1]
            
            j -= 1

#             alignment1 += '-'
#             alignment2 += seq1[j-1]
        
#             j -= 1
            
            currentCell = tracebackMtx[i, j]

        elif currentCell == 2:
            alignment1 += seq1[i - 1]
            alignment2 += '-'
            
            i -= 1

#             alignment1 += seq2[i-1] 
#             alignment2 += '-'

#             i -= 1 
            
            currentCell = tracebackMtx[i, j] 
        print(i,j)
    
    return alignment1[::-1], alignment2[::-1]

def needleman_wunsch(seq1, seq2):
    """Find the global alignment for seq1 and seq2
    Returns: 3 items as so:
    the alignment score, alignment in seq1 (str), alignment in seq2 (str)
    """

    n1, n2 = len(seq1), len(seq2)

    numRows = n1 + 1
    numCols = n2 + 1
    
    match, mismatch, gap = 1, -2, -1
    
    alignMtx, tracebackMtx = create_matrices(numRows, numCols, gap)
    
    for i in range(1, numRows):
        for j in range(1, numCols):
            if seq1[i - 1] == seq2[j - 1]:
                aligned = alignMtx[i - 1, j - 1] + match
            else:
                aligned = alignMtx[i - 1, j - 1] + mismatch
            
            topScore = alignMtx[i - 1, j] + gap # Seq1[i] w/ '-'; top
            leftScore = alignMtx[i, j - 1] + gap # Seq2[j] w/ '-'; left
            
            
            scores = np.array([aligned, leftScore, topScore])
            maxScore = max(scores)
            
            alignMtx[i, j] = maxScore
            tracebackMtx[i, j] = np.argmax(scores)
            
    alignScore = alignMtx[n1, n2]
    alignments = traceback_alignments(seq1, seq2, tracebackMtx)

    return int(alignScore), alignments[0], alignments[1]


if __name__=="__main__":
    Sequences=ReadFASTA(sys.argv[1])
    assert len(Sequences.keys())==2, "fasta file contains more than 2 sequences."
    seq1=Sequences[list(Sequences.keys())[0]]
    seq2=Sequences[list(Sequences.keys())[1]]

    score, align1, align2 = needleman_wunsch(seq1, seq2)

    print('Score: ', score)
    print('Seq1: ', align1)
    print('Seq2: ', align2)
