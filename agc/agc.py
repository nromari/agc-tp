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

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Noura Romari"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Noura Romari"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Noura Romari"
__email__ = "n_a_e@hotmail.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
	parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
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
	   
	   Parametre
       ----------
       
       amplicon_file : fastq file with amplicon sequences
       minseqlen : the minimale length of the selected amplicon sequences
       
       Returns
       --------
       
       seq_gen
    """

	with open(amplicon_file, "r") as fasta_file :
		for line in fasta_file :
			if not line.startswith(">") :
				seq = line.strip
			if len(seq)>=minseqlen :
					yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
	
	"""take a fastq file, a minimale sequence length, 
	   and a sequence minimale count as arguments.
	   return a dictionnary with sequences and the count of each one.
	   
	   Parametre
           ----------
       
       	   amplicon_file : fastq file with amplicon sequences
           minseqlen : the minimale length of the selected amplicon sequences
           mincount : the minimale occurence for an informative sequence
       
           Returns
           --------
       
           seq_dic
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
