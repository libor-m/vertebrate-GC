"""input is 
- gzipped fasta sequence
- window size

output is
- stats table: seq_name win_start win_end A C G T 

TODO: output in wig file, GC content
"""
import sys
import gzip
import collections
import itertools

def fasta(file):
    """multifasta stream parser
    for each sequence in file yields 
    (name, details, seq_generator)
    seq_generator returns the sequence base by base
    """
    def catseq(file):
        """return concatenated lines base by base
        until empty line is found
        """
        for line in file:
            line1 = line.strip()
            if line1 == "":
                return
            for base in line1:
                yield base

    for line in file:
        if line[0] == '>':
            seq_name = line.split()[0][1:]
            yield (seq_name, line.strip()[1:], catseq(file))

def windows(it, size):
    """yield separate iterator for each size 
    elements from it

    TODO: use itertools tee to implemnt overlapping windows
    combined with counter and subtractions
    FIXME: this loses one base a time!
    """
    for item in it:
        yield itertools.islice(it, size)

def main():
    if len(sys.argv) < 2:
        print >> sys.stderr, "use: %s fasta.gz [window_size=100K] [skip=window_size]" % sys.argv[0]
        sys.exit(1)

    if len(sys.argv) > 2:
        win_size = int(sys.argv[2])
    else:
        win_size = 100000
    # TODO: overlapping windows

    file = gzip.open(sys.argv[1])
    for seq, _, data in fasta(file):
        # walk by windows of given size
        for i, win in enumerate(windows(data, win_size)):
            # count different bases in window
            # TODO: add repeat masking obedient flag
            c = collections.Counter(c.upper() for c in win)
            # print the results for each window
            print "\t".join([seq, str(i * win_size)] + [str(c[base]) for base in ['A', 'C', 'G', 'T']])
            # progress indicator
            sys.stderr.write('.') 

if __name__ == '__main__':
    main()
