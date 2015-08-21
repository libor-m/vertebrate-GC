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

def main():
    if len(sys.argv) < 2:
        print >> sys.stderr, "use: %s faidxed_fasta [window_size=100K]" % sys.argv[0]
        sys.exit(1)

    if len(sys.argv) > 2:
        win_size = int(sys.argv[2])
    else:
        win_size = 100000

    # TODO: overlapping windows
    for seq, _, data in fasta(sys.argv[1]):
        
        # output current sequence name
        sys.stderr.write(seq)
        
        # walk by windows of given size
        for start, win in windows(data, win_size):
            # count different bases in window
            # TODO: add repeat masking obedient flag
            try:
            	c = collections.Counter(c.upper() for c in win)
            except:
                print >> sys.stderr, win
            
            # print the results for each window
            print "\t".join([seq, str(start)] + [str(c[base]) for base in ['A', 'C', 'G', 'T']])
            
            # progress indicator
            sys.stderr.write('.') 

if __name__ == '__main__':
    main()
