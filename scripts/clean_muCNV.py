"""
Cleans a muCNV result to make it parsable
$1 - what contig is being parsed
$2 - file to be parsed

outputs to stdout, at which point you'd want to pipe to vcf-sort | bgzip > final.vcf.gz
"""
import re
import sys
import gzip

# Need these lengths
lengths = "/users/u233287/scratch/ccdg_analysis/ref_annotations/grch38.vcf.contigs.txt"
# Fixing the contig header information
my_ctg = ""
with open(lengths, 'r') as fh:
    for line in fh:
        if line[13:].startswith(sys.argv[1] + ','):
            my_ctg = line

if my_ctg == "":
    sys.stderr.write("COULDN'T FIND CONTIG %s\n" % (sys.argv[1]))
    exit(1)
   

def header_fix(line):
    bad1 = "##INFO=<ID=DPOverlap"
    if line.startswith(bad1):
        return "##INFO=<ID=DPoverlap" + line[len(bad1):]
    
    bad2 = "##INFO=<ID=DP2Overlap"
    if line.startswith(bad2):
        return "##INFO=<ID=DP2overlap" + line[len(bad2):]
    
    #if line.startswith(bad2): and line.endswith(
    bad3 = "##INFO=<ID=SPLIT,Number=1"
    if line.startswith(bad3):
        return '##INFO=<ID=SPLIT,Number=.,Type=String,Description="Breakpoints estimated by split reads (start/end pairs in N,N+1 indices)">\n'

    bad7 = "##INFO=<ID=CLIP,Number=1"
    if line.startswith(bad7):
        return '##INFO=<ID=CLIP,Number=.,Type=String,Description="Breakpoints estimated by soft clips (start/end pairs in N,N+1 indices)">\n'
        
    bad4 = "##INFO=<ID=PRE, Number=1"
    if line.startswith(bad4): 
        return "##INFO=<ID=PRE,Number=." + line[len(bad4):]

    bad5 = "##INFO=<ID=POST, Number=1"
    if line.startswith(bad5):
        return "##INFO=<ID=POST,Number=." + line[len(bad5):]
    
    bad6 = "##FORMAT=<ID=DD,Number=A"
    if line.startswith(bad6):
        return "##FORMAT=<ID=DD,Number=." + line[len(bad6):]
    
    bad7 = '##INFO=<ID=DPCNT,Number=1,Type=String,Description="2-D Depth and Read Count clustering">'
    if line.startswith(bad7):
        line = '##INFO=<ID=DPCNT,Number=.,Type=String,Description="2-D Depth and Read Count clustering">\n'
        line += '##INFO=<ID=DPOverlap,Number=1,Type=Float,Description="Overlap between 1-D depth clusters">\n'
        line += '##INFO=<ID=DP2Overlap,Number=.,Type=Float,Description="Overlap between 2-D depth clusters">'
    bad8 = '##INFO=<ID=DP2,Number='
    if line.startswith(bad8):
        line = '##INFO=<ID=DP2,Number=.,Type=String,Description="2-D Depth clustering">'
        
    if line.startswith("#CHROM"):
        extra = '##INFO=<ID=DPCNToverlap,Number=1,Type=String,Description="2-D Depth and Read Count clustering">\n'
        line = extra + my_ctg + line
    return line

def entry_fix(line):
    data = line.strip().split('\t')
    data[7] = re.sub("[\d];[\d]", ",", data[7].replace('(','').replace(')',''))
    data[7] = data[7].replace(";;", ";")
    data[7] = data[7].replace(':DPCNToverlap', ';DPCNToverlap')
    #data[0] = "chr" + data[0]
    return "\t".join(data)

def get_upto_fmt(line):
    """
    Return up to the FORMAT column (faster processing?)
    """
    col = 0
    pos = 0
    while col <= 7:
        pos = line.index('\t', pos) + 1
        col += 1
    return pos
        
def main():
    """
    """
    debug = 10
    with gzip.GzipFile(sys.argv[2], 'r')  as fh:
        for line in fh:
            line = line.decode()
            line = re.sub(":,:", ":.:", line)
            if line.startswith("#"):
                line = header_fix(line)
                sys.stdout.write(line)
            else:
                debug -= 1
                pos = get_upto_fmt(line)
                edited = entry_fix(line[:pos])
                sys.stdout.write(edited + "\t" + line[pos:])
            #if debug <= 0: break

def test():
    for info in sys.stdin:
        info = entry_fix(info)
        print(info)

if __name__ == '__main__':
    #test()
    main()
