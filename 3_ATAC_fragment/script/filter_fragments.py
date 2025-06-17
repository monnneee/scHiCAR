#!/usr/bin/python3
import getopt, sys

def usage():
    print(sys.argv[0] + ' -f fragment -b barcode -o output_file')
    print(sys.argv[0] + ' -h # get help info')

try:
    opts, args = getopt.getopt(sys.argv[1:], "f:b:o:h")
except getopt.GetoptError:
    usage()
    sys.exit(2)

fragment = barcode = output_file = None

for a, o in opts:
    if a == '-f':
        fragment = o
    elif a == '-b':
        barcode = o
    elif a == '-o':
        output_file = o
    elif a == '-h':
        usage()
        sys.exit()

if not all([fragment, barcode, output_file]):
    usage()
    sys.exit(1)

def file2set(filehandle):
    s = set()
    for line in filehandle:
        l = line.strip().split()
        if len(l) >= 2:
            s.add(l[1])
    return s

with open(barcode, 'rt') as bfile:
    barcode_set = file2set(bfile)

with open(fragment, 'rt') as fragfile, open(output_file, 'wt') as outfile:
    for line in fragfile:
        parts = line.strip().split('\t')
        if len(parts) >= 4 and parts[3] in barcode_set:
            outfile.write(line)