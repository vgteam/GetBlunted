#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 18:16:04 2020

@author: Jordan
"""
import re
import sys

comp = {"A":"T", "C":"G", "G":"C", "T":"A"}


def rev_comp(seq):
    return "".join(comp[n] for n in reversed(seq))


class Edge:
    def __init__(self, target, rev, length, backward_length):
        self.target = target
        self.rev = rev
        self.length = length
        self.backward_length = backward_length
    
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        if self.rev:
            r = "-"
        else:
            r = "+"
        return "({}, {}, {}, {})".format(self.target.name, r, self.length, self.backward_length)


class SeqNode:
    def __init__(self, sequence, name):
        self.sequence = sequence
        self.name = name
        self.left = []
        self.right = []
        
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        return "{} {}\n\tleft: {}\n\tright: {}".format(self.name, self.sequence, \
                " ".join(str(e) for e in self.left), " ".join(str(e) for e in self.right))
            

class GFA:
    def __init__(self, filename):
        self.nodes = {}
        with open(filename) as f:
            for l,line in enumerate(f):
                if line.startswith("S"):
                    try:
                        linetype, name, seq = line.strip().split()[0:3]
                        self.add_node(seq, name)
                    except Exception as e:
                        sys.stderr.write(str(e))
                        exit("Error parsing line %d\n\t%s" % (l,line))
        with open(filename) as f:
            for line in f:
                if line.startswith("L"):
                    tokens = line.strip().split()
                    linetype, name1, orientation1, name2, orientation2 = tokens[:5]
                    cigar = ""
                    if len(tokens) > 5:
                        cigar = tokens[5]
                    rev1 = (orientation1 == "-")
                    rev2 = (orientation2 == "-")
                    
                    self.add_edge(name1, rev1, name2, rev2, cigar)
    
    def add_node(self, seq, name):
        self.nodes[name] = SeqNode(seq, name)
    
    def add_edge(self, name1, rev1, name2, rev2, cigar):
        reversing = rev1 != rev2
        node1 = self.nodes[name1]
        node2 = self.nodes[name2]
        len1 = self.cigar_length(cigar, True)
        len2 = self.cigar_length(cigar, False)
        adj1 = None
        if rev1:
            adj1 = node1.left
        else:
            adj1 = node1.right
        adj1.append(Edge(node2, reversing, len1, len2))
        adj2 = None
        if rev2:
            adj2 = node2.right
        else:
            adj2 = node2.left
        adj2.append(Edge(node1, reversing, len2, len1))
    
    def cigar_length(self, cigar, is_ref = False):
        total_len = 0
        i = 0
        j = 0
        while j < len(cigar):
            while not cigar[j].isalpha():
                j += 1
            op_len = int(cigar[i:j])
            op = cigar[j]
            if op == 'M' or op == 'X' or (op == 'D' and is_ref) or (op == 'I' and not is_ref):
                total_len += op_len
            j += 1
            i = j
        return total_len
    
    def __str__(self):
        return "\n".join(str(self.nodes[name]) for name in self.nodes)
    
    
class Subsequence:
    def __init__(self, node, name, begin, end, rev):
        self.node = node
        self.name = name
        self.begin = begin
        self.end = end
        self.rev = rev
    
    def __str__(self):
        if self.rev:
            r = "-"
        else:
            r = "+"
        return "{}[{}:{}]{}".format(self.name, self.begin, self.end, r)

    def __repr__(self):
        return self.__str__()


def parse_translation_table(filename, out_gfa):
    table = {}
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip()
            if len(line) == 0:
                continue
            tokens = line.split()
            assert(len(tokens) == 2)
            subseqs = []
            for subseq in tokens[1].split(","):
                match = re.match("([^\[]+)\[([0-9]+):([0-9]+)\]([\+\-])", subseq)
                name, begin, end, rev = match.group(1), int(match.group(2)), int(match.group(3)), match.group(4) == '-'
                subseqs.append(Subsequence(out_gfa.nodes[tokens[0]], name, begin, end, rev))
            table[tokens[0]] = subseqs
    return table


# every sequence from the input GFA is present at least once
def sequences_are_exhaustive(in_gfa, table):
    success = True

    intervals = {}
    for nodename in table:
        for subseq in table[nodename]:
            if subseq.name not in intervals:
                intervals[subseq.name] = []
            intervals[subseq.name].append((subseq.begin, subseq.end))
    
    for nodename in in_gfa.nodes:
        if nodename not in intervals:
            print("input GFA node {} is missing from translation table".format(nodename), file = sys.stderr)
            success = False
        
        node_intervals = intervals[nodename]
        node_intervals.sort()
        furthest = 0
        for begin, end in node_intervals:
            if begin > furthest:
                print("sequence {} of length {} has no coverage in interval [{}:{})".format(nodename, len(in_gfa.nodes[nodename].sequence), furthest, begin), file = sys.stderr)
                success = False

            furthest = max(furthest, end)
        seq_len = len(in_gfa.nodes[nodename].sequence)
        if furthest != seq_len:
            print("sequence {} of length {} has no coverage in interval [{}:{})".format(nodename, len(in_gfa.nodes[nodename].sequence), furthest, seq_len), file = sys.stderr)
            success = False

    return success

# the sequences that the table indicates are matched between the two GFAs indeed match
def table_sequences_are_consistent(in_gfa, out_gfa, table):
    success = True

    for nodename in table:
        if nodename not in out_gfa.nodes:
            print("translation table node {} is not in output GFA".format(nodename), file = sys.stderr)
            success = False
        if len(out_gfa.nodes[nodename].sequence) == 0:
            print("output node {} does not have a sequence".format(nodename), file = sys.stderr)
            success = False
        
    for nodename in out_gfa.nodes:
        if nodename not in table:
            print("output GFA node {} is missing from translation table".format(nodename), file = sys.stderr)
            success = False
        
        out_node = out_gfa.nodes[nodename]
        for subseq in table[nodename]:
            in_node = in_gfa.nodes[subseq.name]
            if subseq.rev:
                out_seq = rev_comp(out_node.sequence)
                strand = "reverse"
            else:
                out_seq = out_node.sequence
                strand = "forward"
            if in_node.sequence[subseq.begin:subseq.end] != out_seq:
                if len(out_node.sequence) < 400 and len(in_node.sequence) < 400:
                    print("sequence {} of output node {} doesn't match interval [{}:{}) of {} strand of input node {} sequence {}: {}".format(out_node.sequence, out_node.name, subseq.begin, subseq.end, strand, in_node.name, in_node.sequence, in_node.sequence[subseq.begin:subseq.end]), file = sys.stderr)
                elif len(out_node.sequence) < 400 and len(in_node.sequence) > 400:
                    print("sequence {} of output node {} doesn't match interval [{}:{}) of {} strand of input node {} sequence [too long to print]: {}".format(out_node.sequence, out_node.name, subseq.begin, subseq.end, strand, in_node.name, in_node.sequence[subseq.begin:subseq.end]), file = sys.stderr)
                elif len(out_node.sequence) > 400 and len(in_node.sequence) < 400:
                    print("sequence [too long to print] of output node {} doesn't match interval [{}:{}) of {} strand of input node {} sequence {}: {}".format(out_node.name, subseq.begin, subseq.end, strand, in_node.name, in_node.sequence, in_node.sequence[subseq.begin:subseq.end]), file = sys.stderr)
                else:
                    print("sequence [too long to print] of output node {} doesn't match interval [{}:{}) of {} strand of input node {} sequence [too long to print]: {}".format(out_node.name, subseq.begin, subseq.end, strand, in_node.name, in_node.sequence[subseq.begin:subseq.end]), file = sys.stderr)

                success = False

    return success


def is_blunt(out_gfa):
    for nodename in out_gfa.nodes:
        node = out_gfa.nodes[nodename]
        for edge in node.left + node.right:
            if edge.length != 0 or edge.backward_length != 0:
                print("node {} has edge {}, which is not blunt".format(nodename, edge), file = sys.stderr)
                return False
    return True

def adjacencies_are_exhaustive(in_gfa, out_gfa, table):
    
    rev_table = {}
    for nodename in table:
        for subseq in table[nodename]:
            if subseq.name not in rev_table:
                rev_table[subseq.name] = []
            rev_table[subseq.name].append((subseq, nodename))
    
    for nodename in in_gfa.nodes:
        in_node = in_gfa.nodes[nodename]
        
        full_length_paths = []
        
        stack = []
        for subseq, out_nodename in rev_table[nodename]:
            if subseq.begin == 0:
                stack.append((subseq, subseq.rev, False))
        
        
        while len(stack) != 0:
            subseq, rev, walked = stack[-1]
            if walked:
                stack.pop()
                continue
            stack[-1] = (subseq, rev, True)
            
            if subseq.end == len(in_node.sequence):
                full_length_paths.append([(s.node, r) for s, r, w in stack if w])
                continue
            if rev:
                adj = subseq.node.left
            else:
                adj = subseq.node.right
            
            for edge in adj:
                for target_subseq in table[edge.target.name]:
                    on_rev = (edge.rev != rev)
                    if target_subseq.name == nodename and target_subseq.begin == subseq.end and on_rev == target_subseq.rev:
                        stack.append((target_subseq, on_rev, False))
                        
       
        if len(full_length_paths) == 0:
            print("node {} doesn't exist as a walk of paths in the output GFA".format(nodename), file = sys.stderr)
            return False
        
        # because of the POA realignment, the only thing that we can really guarantee is that
        # some prefix/suffix pair are aligned to each other
        for edge in in_node.left + in_node.right:
            found_aligned_base = False
            for full_length_path in full_length_paths:
                 for node, node_rev in full_length_path:
                     expect_rev = node_rev != edge.rev
                     for subseq in table[node.name]:
                         if subseq.name == edge.target.name and expect_rev == subseq.rev:
                             found_aligned_base = True
            if not found_aligned_base:
                print("didn't find any aligned bases in expected orientation between {} and {} despite having overlap edge".format(nodename, edge.target.name), file = sys.stderr)
                print("this can happen with unusual POA graphs, but it probably indicates an error", file = sys.stderr)
                
    
    return True

if __name__ == "__main__":
    
    if len(sys.argv) != 4:
        print("usage: validate_bluntification.py input.gfa output.gfa translation_table.txt", file = sys.stderr)
        sys.exit()
    
    in_gfa = GFA(sys.argv[1])
    out_gfa = GFA(sys.argv[2])
    table = parse_translation_table(sys.argv[3], out_gfa)
    
#    print(in_gfa)
#    print(out_gfa)
#    print(table)
    
    assert(is_blunt(out_gfa))
    assert(table_sequences_are_consistent(in_gfa, out_gfa, table))
    assert(sequences_are_exhaustive(in_gfa, table))
    assert(adjacencies_are_exhaustive(in_gfa, out_gfa, table))
    
    print("No invalid output detected")
    
        
    