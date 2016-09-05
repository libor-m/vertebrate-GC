"""Calculate histogram of bases in given
regions.

Input:
  genome file (gzipped, indexed)
  region file (gff, bed)
Output:
  tab delimited table, long format
  seq name, region index, base, count
"""

import sys
import collections

import pysam
# import fai
import pybedtools

def main():
    if len(sys.argv) < 3:
        sys.exit("use: %s genome regions" % sys.argv[0])
        
    genome = pysam.Fastafile(sys.argv[1])
    # genome = fai.fasta(sys.argv[1])
    regions = pybedtools.BedTool(sys.argv[2])
    
    for idx, r in enumerate(regions):
        seq = genome.fetch(str(r.chrom), r.start, r.end)
        # seq = genome.fetch_sequence(r.chrom, r.start, r.end)
        
        if len(seq) == 0:
            print >> sys.stderr, r.chrom, r.start, r.end, "len == 0"
            break
        
        cnt = collections.Counter(seq)
        
        for k, v in cnt.iteritems():
            items = [r.chrom, r.start, r.end, idx, k, v]
            print "\t".join(map(str, items))
    
if __name__ == "__main__":
    main()
  
  

