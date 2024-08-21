import gzip
import itertools

# def from_file_to_barcode_list(filename):
#     '''
#     Load the barcode whitelist into the memory
#     '''
#    
#     # f = open(filename,"r")
#     f = gzip.open(filename,'rt')
#     barcodes = []
#     while(True):
#         line = f.readline().rstrip('\n')
#         if not line:
#             break
#         barcodes.append(line)
#     f.close()
#     print()
#     return set(barcodes)

def from_file_to_barcode_list(filename):  ## add the whitelist check
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
    return set(barcodes)


def newbarcode_white_list_dic(barcode_map):
    barcode_dic = {}
    with open(barcode_map,'r') as f:
        for line in f:
            (mismatch_barcode,white_list_barcode) = line.strip().split(" ")
            barcode_dic[mismatch_barcode] = white_list_barcode
    return barcode_dic

# def update_fastq(r1,r2,out_r1,out_r2, barcode_dic ): ## process two files
def update_fastq(r1,out_r1, barcode_dic ,white_list_barcode,log_file): ## process one files
    """"
    modify the barcode
    """
    f_r1 = open(r1, 'r')
    f_out_r1 = open(out_r1, 'w') 
    num_of_fragments = 0 
    keeped_reads = 0
    while True:
        num_of_fragments+=1
        cur_r1_name = f_r1.readline().strip()[1:] # remove @
        cur_r1_read = f_r1.readline().strip()
        cur_r1_plus = f_r1.readline().strip()
        cur_r1_qual = f_r1.readline().strip()
        # 
        # cur_r2_name = f_r2.readline().strip()[1:]
        # cur_r2_read = f_r2.readline().strip()
        # cur_r2_plus = f_r2.readline().strip()
        # cur_r2_qual = f_r2.readline().strip()
    
        if cur_r1_name == "" : break
                
        cur_barcode = (cur_r1_name[:28].upper())  ## makesure the barcode length

        if   cur_barcode in white_list_barcode : ## only the whitelist or corrected barcode are remained
        #if  cur_barcode in barcode_dic or cur_barcode in white_list_barcode : ## only the whitelist or corrected barcode are remained
         #   if cur_barcode in barcode_dic:
          #      cur_r1_name = (barcode_dic[cur_barcode]+ cur_r1_name[28:]) ## 28bp cell barcode
        
            f_out_r1.write('@' + cur_r1_name +"\n")
            f_out_r1.write(cur_r1_read+"\n")
            f_out_r1.write(cur_r1_plus+"\n")
            f_out_r1.write(cur_r1_qual+"\n")     
            keeped_reads+=1
    
        # f_out_r2.write('@' + cur_r1_name +"\n")
        # f_out_r2.write(cur_r1_read+"\n")
        # f_out_r2.write(cur_r1_plus+"\n")
        # f_out_r2.write(cur_r1_qual+"\n")    


    f_r1.close()
    f_out_r1.close()
    
    f_out_barcode_log = open(log_file,"w+")
    f_out_barcode_log.write(" %d reads, %d keeped" % (num_of_fragments, keeped_reads))
    f_out_barcode_log.close()

r1 = str(snakemake.input[0])
barcode_dic = newbarcode_white_list_dic(snakemake.input[1])
# r2 = str(snakemake.input['r2'])

# white_list_barcode = from_file_to_barcode_list('sciHCAR_whitelist_20Oct25_test.txt.gz')
white_list_barcode = from_file_to_barcode_list('sciHiCAR_18bp_barcode_440k.txt.gz')

out_r1 = snakemake.output[0]
# out_r2 = snakemake.output['r2']
log_file = snakemake.log[0]

# update_fastq(r1,r2,out_r1,out_r2,barcode_dic )
update_fastq(r1,out_r1,barcode_dic ,white_list_barcode  ,log_file )

