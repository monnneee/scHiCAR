#!/usr/bin/python3
import getopt, sys

def usage():
    print(sys.argv[0] + ' -p pair -b barcode -o output_file')
    print(sys.argv[0] + ' -h # get help info')

try:
    opts, args = getopt.getopt(sys.argv[1:], "p:b:o:h")
except getopt.GetoptError:
    usage()
    sys.exit(2)

pair = barcode = output_file = None

for a, o in opts:
    if a == '-p':
        pair = o
    elif a == '-b':
        barcode = o
    elif a == '-o':
        output_file = o
    elif a == '-h':
        usage()
        sys.exit()

if not all([pair, barcode, output_file]):
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

with open(pair, 'rt') as pairfile, open(output_file, 'wt') as outfile:
    for line in pairfile:
        if line.startswith('#'):
            continue
        bar = line[:18]
        if bar in barcode_set:
            outfile.write(line)