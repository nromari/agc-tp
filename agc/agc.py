#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""
# Modules

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
#import nwalign3 as nw

__author__ = "Noura Romari"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Noura Romari"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Noura Romari"
__email__ = "n_a_e@hotmail.fr"
__status__ = "Developpement"

#==============================================================
# Definition des fonctions
#==============================================================

def isfile(path):
    """verifie que l'existence du fichier
      :Parameters:
          path: chemin du fichier
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Recupère les arguments du programme
      Returns: un objet contenant les arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile,
                        required=True, help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):

    """take a fastq file as argument and return a sequence generator,
    with a minimale sequence length = minseqlen.
    """

    with gzip.open(amplicon_file, "rt") as fasta_file :
        for line in fasta_file :
            if not line.startswith(">") :
                seq = line.strip
                if len(seq)>=minseqlen :
                    yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):

    """take a fastq file, a minimale sequence length,
    and a sequence minimale count as arguments.
    return a dictionnary with sequences and the count of each one.
    """

    seq_dic = {}
    for seq in read_fasta(amplicon_file, minseqlen) :
        if seq in seq_dic :
            seq_dic[seq] += 1
        else :
            seq_dic[seq] = 1
    sec_dic_tr = sorted(seq_dic.items(), key=lambda item: item[1], reverse=True)
    for seq, count in seq_dic_tr.items():
        if count >= mincount :
            yield [seq, count]


def get_chunks(sequence, chunk_size) :

    """take a sequence and a chunk_size as arguments,
       return a chunk-size length sub-sequences list
    """

    chunk_seq = []
    i=1
    while i*chunk_size < len(sequence) :
        sub_seq = sequence[i*chunk_size-chunk_size:i*chunk_size]
        chunk_seq.append(sub_seq)
    i += 1
    return chunk_seq


def cut_kmer(sequence, kmer_size) :

    """take a sequence and a kmer_size as arguments,
       return a kmer genrator
    """

    for i in range(len(sequence) - kmer_size + 1) :
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size) :

    """Arguments : kmer index dictionnary, a sequence, an id sequence,
	   and a kmer size.
	   return a kmer dictionnary
	"""

    kmer_dict = {}
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer]=[id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size) :

    """Return the 8 most common sequences of a kmer dictionnary
    """

    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size)
           if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list) :

    """Return the idendity percent af an alignement
    """

    nb_base = 0
    for i in range(0, len(alignment_list[0])) :
        if alignment_list[0][i] == alignment_list[1][i] :
            nb_base += 1
    return nb_base/len(alignment_list[0])*100


def get_unique(ids):

    """
    used only for test
    """

    return {}.fromkeys(ids).keys()


def detect_chimera(perc_identity_matrix) :

    """Take an identity matrice between 2 sequences (ref and new) as argument
       return "True" if the new sequence is a chimera, and "False" if not
    """

    std_dev = []
    booleen_1 = False
    booleen_2 = False

    for i in range(len(perc_identity_matrix)) :
        std_dev.append(std(perc_identity_matrix[i]))
        perc_0 = perc_identity_matrix[i][0]
        perc_1 = perc_identity_matrix[i][1]
        if perc_0 > perc_1 :
            booleen_1 = True
        if perc_0 < perc_1 :
            booleen_2 = True
   
    if statistics.mean(std_dev) > 5 and booleen_1 and booleen_2 :
        return True
    else:
        return False


def common(lst1, lst2):

    """Return the common items between two lists
    """

    return list(set(lst1) & set(lst2))



def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size) :
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size) :
    pass

def write_OTU(OTU_list, output_file) :
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()
