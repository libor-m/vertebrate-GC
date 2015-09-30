"""input is 
- faidx indexed fasta sequence
- window size

output is
- stats table: seq_name win_start win_end A C G T 

Use pysam.Fastafile to access the indexed data.

TODO: output in wig file, GC content
"""
import sys
import collections
import itertools
import argparse

import pysam

def fai(file):
    """read and parse fai into a list of tuples
    """
    with open(file) as f:
        lines = f.readlines()
        return [tuple(l.strip().split('\t')) for l in lines]
    
def fasta(file):
    """load index, 
    use pysam to access all sequences mentioned in the index
    """
    idx = fai(file + '.fai')
    fa = pysam.Fastafile(file)
    
    for seq, _, _, _, _ in idx:
        yield seq, seq, fa.fetch(seq)

def windows(data, size):
    """walk a window along a sequecne
    """
    for i in xrange(0, len(data) - size + 1, size):
        yield i, data[i:(i + size)]

def argparser():
    args = argparse.ArgumentParser(description="Count base counts in windows along chromosomes.")
    args.add_argument("--window", type=int, default=100000,
        help="Window size (default 100k).")
    args.add_argument("--bases", default="ACGTacgt",
        help="Bases to count (default: ACGTacgt")
    args.add_argument("-q", "--quiet", action="store_true",
        help="Do not report progress.")
    args.add_argument("--ucase", action="store_true",
        help="Uppercase input (to remove soft masking).")
    args.add_argument("input", 
        help="Multifasta with .fai index (can be gzipped).")
    return args

def main():
    parser = argparser()
    args = parser.parse_args()

    # TODO: overlapping windows
    for seq, _, data in fasta(args.input):
        
        # output current sequence name
        if not args.quiet:
            sys.stderr.write(seq)
        
        # walk by windows of given size
        for start, win in windows(data, args.window):
            # uppercase each window if required
            if args.ucase:
                win = win.upper()

            # count different bases in window
            c = collections.Counter(win)
            
            # print the results for each window
            print "\t".join([seq, str(start)] + [str(c[base]) for base in args.bases])
            
            # progress indicator
            if not args.quiet:
                sys.stderr.write('.') 

if __name__ == '__main__':
    main()
