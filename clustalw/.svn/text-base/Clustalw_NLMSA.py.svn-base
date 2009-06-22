#Author Eden Elos
#May28,08

#convert clustalw2 output format to multiple alignment format(maf)
#create NLMSA from Clustalw alignment file

#! /usr/bin/env python2.5

import os
import glob
from pygr import cnestedlist, seqdb

class ClustalwResidues:
    
    """
    A single clustalw 'residue'/alignment block containing
    upto 60 residues from each sequence
    """
    
    def __init__(self, no_seq, seq_names, seqs, start_indices):
        self.no_seq = no_seq
        self.seq_names = seq_names
        self.seqs = seqs
        self.start_indices =  []
        #self.start_indices=start_indices
        for i in range(0, len(start_indices)):
            self.start_indices.append(start_indices[i])

    def get_no_seqs(self):
        """
         returns the number of sequences in the alignment block
        """
        return self.no_seq

    def get_names(self):
        """
          returns list of sequence names
        """
        return self.seq_names

    def get_seqs(self):
        """
         returns list of sequence
        """
        return self.seqs

    def get_start_indices(self):
        """
          returns list of  start indices for each
          sequence in the alignment block 
        """
        return self.start_indices

    def reset_start_indices(self):
        """
          reset the start indices if the sequence in the
          alignment block is all gaps this is done by
          subtracting 1 from corresponding start_indices...
        """
        ungap_ct = self.ungapped_count()
        for i in range(0, self.no_seq):
            if ungap_ct[i] == 0:
                self.start_indices[i] -= 1

    def get_end_indices(self):
        """
         returns list of end indices of each sequence in the alignment block.
        """
        ungapped_lens = self.ungapped_count()
        end_indices = []
        for i in range(0,self.no_seq):
            #should we have -1 here... convention is 0-based indexing
            if ungapped_lens[i] != 0:
                end_indices.append(self.start_indices[i] + ungapped_lens[i]-1)
            else:
                end_indices.append(self.start_indices[i] + ungapped_lens[i])
        return end_indices
        
    def gap_count(self):
        """
        returns the number of gaps in each sequence
        """
        gaps = []
        for s in self.seqs:
            gaps.append(s.count("-"))

        return gaps

    def ungapped_count(self):
        """
        returns the ungapped length
        """

        ungapped_len = []
        gaps = self.gap_count()
        
        for i in range(0,len(self.seqs)):
            ungapped_len.append(len(self.seqs[i])-gaps[i])

        return ungapped_len

         
def read_clustalw(f):
    """
    Read aligned sequences from a CLUSTALW alignment
    """
    fp = open(f, "r")
    lines = fp.readlines()
    fp.close()
    
    assert lines[0].startswith('CLUSTAL '), lines[0]
    lines = lines[3:]
    seq_counter = 0 #identify the number of sequences
    
    while True:
        if lines[seq_counter][:16].strip():
            seq_counter += 1
        else:
            break
    clustal_res_list = [] #holds alignment blocks as a list of clustalResidues
    start_indices = [] #holds the start_indices of each sequence in the alignment blocks

    for i in range(0, seq_counter):
        start_indices.append(0)

    for i in range(0, len(lines), seq_counter+2):
        seq_names = []
        seqs = []
        for j in range(0, seq_counter):
            line = lines[i+j].strip()
            ls = line.split()
            #sometimes, especially at the last block, the length value might not be there
            #if the whole string of characters is composed on gaps only
            #not required by the standard ?
            if len(ls) == 3:
                name, seq, length = ls
            else:
                name, seq = ls
                length = 0
                
            seq_names.append(name)
            seqs.append(seq)
        cl = ClustalwResidues(seq_counter, seq_names, seqs, start_indices)
        ungapped_len = cl.ungapped_count()

        for j in range(0, len(start_indices)):
            start_indices[j] += ungapped_len[j]
                     
        #check the current cl to reset the start_indices cl...
        #in case a sequence in the current cl is completely a gap,
        #the start_index has to be decreased by one
        cl.reset_start_indices()
        clustal_res_list.append(cl)
        
    return clustal_res_list

def calc_total_length(clustal_res_list):
    total_lengths = clustal_res_list[0].ungapped_count()
    
    for cl_res in clustal_res_list[1:]:
        ls = cl_res.ungapped_count()
        total_lengths = [total_lengths[i] + ls[i] for i in range(0, len(ls))]

    return total_lengths
        
        
def write_MAF(fp, clustal_res_list):
    """ takes a list of ClustalwResidues/ clustalw alignment blocks and constructs the corresponding
        and writes MAF file format
        refer to the MAF file format at http://genome.ucsc.edu/FAQ/FAQformat#format5
        
        N.B. Each block may have a subset of the input sequences.
            All sequence columns and rows must contain at least one nucleotide
            (no columns or rows that contain only insertions)
        
    """
    f = open(fp, "w")

    #to be changed ....
    #there are some information that are not found in Clustalw output
    
    #f.write("track name=sample visibility=pack\n")
    f.write("##maf version=1 scoring=tba.v8\n")
    f.write("# ********\n\n")
    
    total_lengths = calc_total_length(clustal_res_list)
    

    
    for cl_res in clustal_res_list:
        seq_count = cl_res.get_no_seqs()
        names = cl_res.get_names()
        seqs = cl_res.get_seqs()
        start_indices = cl_res.get_start_indices()
        #print start_indices," ****"
        ungapped_len = cl_res.ungapped_count()

        #make sure at least a pair of sequences aligned in the block
        ct = 0
        for i in range(0, seq_count):
            if ungapped_len[i] > 0:
                ct += 1
        if ct > 1:
            f.write("a   score=00000.0\n")
            for i in range(0,seq_count):
                if ungapped_len[i] > 0:
                    line_st="s "+names[i]+" "+str(start_indices[i])+" "+\
                             str(ungapped_len[i])+" + "+str(total_lengths[i])+" "
                    line_st+=seqs[i]
                    f.write(line_st)
                    f.write("\n")
            f.write("\n")
    f.close()

def build_interval_list(a, b):
    """
    Hacky code to extract all ungapped aligned subintervals from a
    pair of aligned sequences.
    """
    interval_list = []

    a_start = None
    b_start = None

    a_count = b_count = 0
    for i in range(0, len(a)):
        if a[i] == '-' or b[i] == '-':
            if a_start is not None:           # want to end at i-1
                interval_list.append((a_start, a_count, b_start, b_count))

                a_start = b_start = None
        else:
            if a_start is None:
                a_start = a_count
                b_start = b_count

        if a[i] != '-':
            a_count += 1
        if b[i] != '-':
            b_count += 1

    if a_start is not None:
        interval_list.append((a_start, a_count, b_start, b_count))

    assert a_count == len(a.replace('-', ''))
    assert b_count == len(b.replace('-', ''))
        
    return interval_list


def create_NLMSA_clustalw(clu_f, seqDb, al):
    """
        takes clustalw alignment file as input and creates and returns NLMSA
    """
    clustal_res_list = read_clustalw(clu_f)
    genome_names = clustal_res_list[0].get_names() 
    genomes_dict = {}
    
    for genome in genome_names:
        genomes_dict[genome] = seqDb[genome]

    al += genomes_dict[genome_names[0]]

    for clu_res in clustal_res_list:    
        # build list of aligned sub-intervals

        seq = clu_res.get_seqs()
        start_indices = clu_res.get_start_indices()
        end_indices = clu_res.get_end_indices()

                
        for i in range(0, len(seq)):
            start1 = start_indices[i]
            stop1 = end_indices[i]
            if start1 != stop1:
                genome1_ival = genomes_dict[genome_names[i]][start1:stop1+1]
                genome1_ival_str = seq[i]
        
                for j in range(i+1, len(seq)):
                    start2 = start_indices[j]
                    stop2 = end_indices[j]
                    if start2 != stop2:
                        genome2_ival = genomes_dict[genome_names[j]][start2:stop2+1]
                        genome2_ival_str = seq[j]
                
                        interval_list = build_interval_list(genome1_ival_str,
                                                            genome2_ival_str)

                        for (a, b, x, y) in interval_list:
                            ival1 = genome1_ival[a:b]
                            ival2 = genome2_ival[x:y]
                            al[ival1] += ival2

    al.build()
    return al

    
