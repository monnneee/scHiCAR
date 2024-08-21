import getopt,sys
import time
import gzip

def from_file_to_barcode_list(filename):
    '''
    Load the barcode whitelist into the memory
    '''

    # f = open(filename,"r")
    f = gzip.open(filename,'rt')
    barcodes = []
    while(True):
        line = f.readline().rstrip('\n')
        if not line:
            break
        barcodes.append(line)
    f.close()
    print()
    return barcodes

def extract_from_line(line):
    '''
    This function is used to extract barcode count info from a line.
    '''
    curr = line.strip().split(" ")
    return (int(curr[0]),curr[1])



def if_one_mismatch(barcode_set,barcode,f_out_barcode_map):
    bp = ['A','G','T','C']
    b_found = False

    for i in range(len(barcode)):
        for each in bp:
            if each == barcode[i].upper():
                continue
            newbarcode = barcode[:i]+each+barcode[i+1:] ## create 1 miss-mismatch barcode
            if newbarcode in barcode_set:
                if b_found: # return in the loop 
                    return 0,"" ## more than one hit in the whitelist
                b_found = True
                corresponding_whitelist_barcode = newbarcode ## first one
    if b_found: # walk all combination 
        return 1,corresponding_whitelist_barcode ## only one mismatch
    return 0,"" ## no one mismatch found in whilelist


def catagorize_barcode(line,barcode_set,f_out_barcode_map,barcodes_slices,slicer):
    '''
    The following code is used to handle each barcode in the fragment file.
    It will insert an info entry into the barcode_info list.
    An info entry has a form of (barcode, # of fragments, 0/1 mismatch).
    '''
    (number_of_fragments, barcode) = extract_from_line(line.rstrip('\n'))

    corresponding_whitelist_barcode = ''
    corresponding_whitelist_barcode_slice = ''
    one_mismatch = 0

    if barcode in barcode_set: ## perfect match
        match_type = 0
    else:
        for i in range(len(slicer)):
            if i == 0:
                barcode_piece = barcode[:slicer[i]]
            else:
                barcode_piece = barcode[sum(slicer[:i]):sum(slicer[:(i+1)])]


            # Check if the current barcode slice has a perfect mismatch in the whitelist
            if barcode_piece in barcodes_slices[i]:
                corresponding_whitelist_barcode += barcode_piece
                continue

            # Check if the current barcode slice has unique corresponding one mismatch in the whitelist
            one_mismatch,corresponding_whitelist_barcode_slice = if_one_mismatch(barcodes_slices[i],barcode_piece,f_out_barcode_map)
            if not one_mismatch:
                match_type = 2 # Cannot find barcode correction
                return (barcode, number_of_fragments, match_type)
            else:
                corresponding_whitelist_barcode +=corresponding_whitelist_barcode_slice

        f_out_barcode_map.write(" ".join([barcode,corresponding_whitelist_barcode,"\n"]))
        match_type = 1 # Not a perfect match but can find a correction with 1-mismatch tolerence for each slice

    return (barcode, number_of_fragments, match_type)

def find_barcode_info(fragmentsfilename, barcodes, barcode_map_file, barcode_log_file, OUTPUT_FILENAME,slicer):

    # Slicer Check

    if(sum(slicer)!=len(barcodes[0])):
        print('Invalid slicer input values! Barcode length is %d.' % len(barcodes[0]))
        return

    # Build barcode index

    barcodes_slices = [[] for _ in range(len(slicer))]
    for each in barcodes:
        for i in range(len(slicer)):
            if i == 0:
                barcodes_slices[i].append(each[:slicer[i]])
            else:
                barcodes_slices[i].append(each[ sum(slicer[:i]) : sum(slicer[:(i+1)])])


        #barcode_set.append(set(barcodes)) # A set structre makes the lookup faster.

    barcode_set = set(barcodes)
    if len(barcodes) == len(barcode_set):
        print("%d barcodes are provided. All of them are unique." % len(barcodes))
    else:
        print("%d barcodes are provided. %d of them are unique." % (len(barcodes),len(barcode_set)))

    product = 1
    for i in range(len(slicer)):
        barcodes_slices[i] = set(barcodes_slices[i])
        product = product*len(barcodes_slices[i])
        print("For slice part %d, %d of them are unique" % (i,len(barcodes_slices[i])))

    print('There are %d possible slice combinations. There are %d provided barcodes.' % (product,len(barcode_set)))
    if product == len(barcode_set):
        print('All possible barcodes are provided!')
    else:
        print('Warning: not all possible barcodes are provided!')



    # f = gzip.open(fragmentsfilename,"r")
    f = open(fragmentsfilename,"r")
    f_out_info = open(OUTPUT_FILENAME,"w+")
    f_out_barcode_map = open(barcode_map_file,"w+")
    f_out_barcode_log = open(barcode_log_file,"w+")
    # f = open(outputfile,"w")
    print('Catagorizing barcodes...')

    frag_info = [0,0,0]
    barcode_info = [0,0,0]
    match_type = 0
    # 0 - perfect match
    # 1 - unique 1-mismatched
    # 2 - 2+ misamtch or 1 mismatch to 2+ barcodes.
    num_of_fragments = 0 # Number of fragments in this line.

    total_fragements = 0

    while(True):
        line = f.readline().rstrip('\n')
        if not line:
            break
        total_fragements+=1
        (barcode, num_of_fragments,match_type) = catagorize_barcode(line,barcode_set,f_out_barcode_map,barcodes_slices,slicer)
        f_out_info.write(" ".join([barcode,str(num_of_fragments),str(match_type),"\n"]))
        frag_info[match_type]+=num_of_fragments
        barcode_info[match_type]+=1
    f.close()
    f_out_info.close()
    f_out_barcode_map.close()

    # Display results

    # print('\n')
    # print('Number of Lines in Total: %d' % (total_fragements))
    # print('Number of Barcodes Provided: %d' % len(barcodes))

    # print("0   mismatch: %d fragments from %d barcodes" % (frag_info[0],barcode_info[0]))
    # print("1   mismatch: %d fragments from %d barcodes" % (frag_info[1],barcode_info[1]))
    # print("2+  mismatch: %d fragments from %d barcodes" % (frag_info[2],barcode_info[2]))
    # print("Match with 2: %d fragments from %d barcodes" % (frag_info[3],barcode_info[3]))
    f_out_barcode_log.write('\n')
    f_out_barcode_log.write('Total sequenced barcode: %d \n' % (total_fragements))
    f_out_barcode_log.write('whitelist bacrode Provided: %d \n' % len(barcodes))
    f_out_barcode_log.write("perfect match: %d fragments from %d barcodes \n" % (frag_info[0],barcode_info[0]))
    f_out_barcode_log.write("1_mismatch: %d fragments from %d barcodes \n" % (frag_info[1],barcode_info[1]))
    f_out_barcode_log.write("not found : %d fragments from %d barcodes \n" % (frag_info[2],barcode_info[2]))
    f_out_barcode_log.write("%.2f reads  from %.2f barcode  are valided \n" % ((frag_info[0] + frag_info[1])/sum(frag_info), (barcode_info[1]+barcode_info[0])/ total_fragements))
    f_out_barcode_log.close()

    return



def write_output(barcode,fragments):
    '''
    This function writes line info into a file.
    Var fragments refers to the number of fragments matched with the given barcode.
    Var type: 0 refers to perfect match & 1 refers to 1 mismatch.
    '''
    f_out = open(OUTPUT_FILENAME,"w+")

    f_out.close()

"""
# the single file version
# read commandline arguments, first
fullCmdArguments = sys.argv
# - further arguments
argumentList = fullCmdArguments[1:]
unixOptions = "b:f:h0:"
gnuOptions = ["barcodes=", "fragments=","help","output="]
try:
    arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
except getopt.error as err:
    # output error, and return with an error code
    print (str(err))
    sys.exit(2)
OUTPUT_FILENAME = "barcode_info.txt"
if(len(arguments)==0):
    # barcodes = from_file_to_barcode_list('all-barcodes-test.txt')
    barcodes = from_file_to_barcode_list('barcode-output.txt')
    fragement_filename = 'sorted-fragments-test.txt'
for currentArgument, currentValue in arguments:
    if currentArgument in ("-b", "--barcodes"):
        # print(currentValue)
        barcodes = from_file_to_barcode_list(currentValue)
    elif currentArgument in ("-f", "--fragments"):
        # print(currentValue)
        fragement_filename = currentValue
    elif currentArgument in ("-o", "--output"):
        OUTPUT_FILENAME = currentValue
    elif currentArgument in ("-h", "--help"):
        print("-b --barcodes  : barcodes.txt")
        print("-f --fragments : fragments.txt")
        print("-o --output    : output-filename (default output.txt)")
        exit()
"""

## update the argument using through snakemake

'''
fragement_filename = snakemake.input[0]
barcodes = from_file_to_barcode_list(snakemake.config['white_list']) ## change the position ?
barcode_map_file = snakemake.output['map']
barcode_log_file = snakemake.output['log']
OUTPUT_FILENAME =  snakemake.output['sum']
'''

# Using testing files here
# fragement_filename = './Sci-HiC-ATAC-2_raw_barcode_count.txt'
# barcodes = from_file_to_barcode_list('./sciHCAR_whitelist_20Oct25_test.txt.gz') # change the position ?
# barcode_map_file = 'map_2.txt'
# barcode_log_file = 'log_2.txt'
# OUTPUT_FILENAME =  'output_2.txt'

# Initialize slicer
slicer = [6,6,6]




fragement_filename = snakemake.input[0]
barcodes = from_file_to_barcode_list('linear_sciRNA_18bp_barcode.txt.gz') ## change the position ## use the fulllength
barcode_map_file = snakemake.output['map']
barcode_log_file = snakemake.output['log']
OUTPUT_FILENAME =  snakemake.output['sum']
# find_barcode_info(fragement_filename,barcodes,barcode_map_file, barcode_log_file,OUTPUT_FILENAME)
find_barcode_info(fragement_filename,barcodes,barcode_map_file, barcode_log_file,OUTPUT_FILENAME,slicer)
