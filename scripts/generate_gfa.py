#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 23:06:54 2020

@author: Jordan
"""

import random
import math
import re

comp = {"A":"T", "C":"G", "G":"C", "T":"A"}
def rev_comp(seq):
    return "".join(comp[n] for n in reversed(seq))

def dirichletvariate(alphas):
    gams = [random.gammavariate(alpha, 1) for alpha in alphas]
    denom = sum(gams)
    return tuple(gam/denom for gam in gams)

def geomvariate(p):
    return 1 + int(random.expovariate(-math.log(1 - p)))

alphabet = "ACGT"

def random_seq(length):
    return "".join(random.choice(alphabet) for i in range(length))

def mutate(seq, mut_rate, indel_rate):
    nts = list(seq)
    for i in range(len(nts)):
        if random.random() < indel_rate * 0.5:
            nts[i] = ""
            continue
        if random.random() < mut_rate:
            nts[i] = random.choice(alphabet)
        if random.random() < indel_rate * 0.5:
            nts[i] += random.choice(alphabet)
    return "".join(nts)

def random_seq_with_repeats(length, identity = 0.9, expected_lib_size = 10.0):
    library = []
    library_length = 0
    while library_length < length:
        chunk_len = geomvariate(expected_lib_size / length)
        library_length += chunk_len
        library.append(random_seq(chunk_len))
    
    # give them mixture weights
    weights = dirichletvariate([1 for i in range(len(library))])
    
    # factor in mutation to the same base (still valid for high identities)
    # and evolution distance of 2T from library ancestor
    mut_rate = 0.8 * (1.0 - identity) * 2.0 / 3.0
    indel_rate = 0.2 * (1.0 - identity) * 2.0 / 3.0
    
    seqs = []
    seq_len = 0
    while seq_len < length:
        lib_chunk = random.choices(library, weights)[0]
        chunk = mutate(lib_chunk, mut_rate, indel_rate)
        seq_len += len(chunk)
        seqs.append(chunk)
    
    return "".join(seqs)[:length]

def sample_reads(seq, cov, min_read_len, exp_read_len, error_rate):
    assert(exp_read_len > min_read_len)
    reads = []
    total_cov = 0
    while total_cov < len(seq) * cov:
        i = random.randrange(len(seq))
        read_len = min_read_len + int(random.gammavariate(5.0, (exp_read_len - min_read_len) / 5.0))
        if read_len > len(seq):
            continue
        if i + read_len > len(seq):
            read = seq[i:] + seq[:read_len - (len(seq) - i)]
        else:
            read = seq[i:i+read_len]
        read = mutate(read, error_rate * .8, error_rate * .2)
        total_cov += len(read)
        if random.random() < .5:
            reads.append(read)
        else:
            reads.append(rev_comp(read))
    return reads

def overlap(left, right):
    
    mat = [[-10000000000 for j in range(len(right) + 1)] for i in range(len(left) + 1)]
    for i in range(len(left) + 1):
        mat[i][0] = 0
    
    for i in range(1, len(left) + 1):
        for j in range(1, len(right) + 1):
            if left[i - 1] == right[j - 1]:
                match = 1
            else:
                match = -1
            mat[i][j] = max(
                    mat[i - 1][j - 1] + match,
                    mat[i - 1][j] - 4,
                    mat[i][j - 1] - 4
               )
    
    
    max_j = 0
    for j in range(1, len(right) + 1):
        if mat[len(left)][j] > mat[len(left)][max_j]:
            max_j = j
    
    score = mat[len(left)][max_j]
    
    
#    print("\t" + "\t".join(right))
#    for i in range(len(mat)):
#        if i > 0:
#            print(left[i - 1], end = "")
#        print("\t" + "\t".join(str(v) for v in mat[i]))
    
    # start with a dummy op
    cigar_ops = [['N', 0]]
    i = len(left)
    j = max_j
    while j > 0:
        if left[i - 1] == right[j - 1]:
            match = 1
        else:
            match = -1
            
        if mat[i][j] == mat[i - 1][j - 1] + match:
            op = 'M'
            i -= 1
            j -= 1
        elif mat[i][j] == mat[i - 1][j] - 4:
            op = 'D'
            i -= 1
        else:
            op = 'I'
            j -= 1
        
        if op == cigar_ops[-1][0]:
            cigar_ops[-1][1] += 1
        else:
            cigar_ops.append([op, 1])
    
    cigar = ""
    for op, length in cigar_ops[len(cigar_ops) - 1:0:-1]:
        cigar += str(length) + op
    
    return score, cigar

def flip_cigar(cigar):
    flipped = ""
    for match in re.finditer("(\d+)([MDI])", cigar):
        length = int(match.group(1))
        op = match.group(2)
        if op == "D":
            op = "I"
        elif op == "I":
            op = "D"
        flipped += str(length) + op
    return flipped

if __name__ == "__main__":
    
    #random.seed(1)
    
    ref_length = 200
    ref_chunk_identity = .9
    ref_expected_num_chunks = 3
    coverage = 3
    min_read_length = 40
    mean_read_length = 60
    read_error_rate = 0.05
    min_score = 8
    
    seq = random_seq_with_repeats(ref_length, ref_chunk_identity, ref_expected_num_chunks)
    
    reads =  sample_reads(seq, coverage, min_read_length, mean_read_length, read_error_rate)
    
    print("H\tVN:Z:1.0")
    for i in range(len(reads)):
        print("S\t{}\t{}".format(i, reads[i]))
    
    
    for i in range(len(reads)):
        for j in range(i + 1, len(reads)):
            for i_left_side in [True, False]:
                for j_left_side in [True, False]:
                    if i_left_side:
                        seq_i = rev_comp(reads[i])
                    else:
                        seq_i = reads[i]
                    if j_left_side:
                        seq_j = reads[j]
                    else:
                        seq_j = rev_comp(reads[j])
                    
                    #print("{} {}: {} {}".format(i, j , i_left_side, j_left_side))
                    
                    score, cigar = overlap(seq_i, seq_j)
                    if score >= min_score:
                        
                        if random.random() < .5:
                            cigar = flip_cigar(cigar)
                            label1 = j
                            label2 = i
                            if j_left_side:
                                strand2 = '-'
                            else:
                                strand2 = '+' 
                            if i_left_side:
                                strand1 = '+'
                            else:
                                strand1 = '-'
                        else:
                            label1 = i
                            label2 = j
                            if i_left_side:
                                strand1 = '-'
                            else:
                                strand1 = '+'
                            if j_left_side:
                                strand2 = '+'
                            else:
                                strand2 = '-'
                        
                        print("L\t{}\t{}\t{}\t{}\t{}".format(label1, strand1, label2, strand2, cigar))
                    
                    
            
            
