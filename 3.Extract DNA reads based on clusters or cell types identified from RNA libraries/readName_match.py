#!/usr/bin/python3
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], "l:s:o:h")

def usage():
        print(sys.argv[0] + ' -l list_file -s set_file -o output_file')
        print(sys.argv[0] + ' -h #get help info')


for a,o in opts:
        if a in ('-l'):
                list_file = o
        elif a in ('-s'):
                set_file = o
        elif a in ('-o'):
                output_file = o
        elif a in ('-h'):
                usage()
                sys.exit()
def file2set(filename):
    s = set()
    for line in filename:
        l = line.strip('\n').split('\t')
        s.add(l[0])
    return s

with open(list_file, 'rt') as file1, open(set_file, 'rt') as file2, open(output_file, 'wt') as file3:
    set1 = file2set(file2)

    output_lines = (line[1:] for line in file1 if line[1:29] in set1)

    file3.writelines(output_lines)
