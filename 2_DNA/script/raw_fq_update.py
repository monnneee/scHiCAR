import gzip


def update_fastq(r1,r2, out_r1,out_r2 ): ## process two files

    f_r1 = gzip.open(r1, 'rt')
    f_r2 = gzip.open(r2, 'rt')
    f_out_r1 = open(out_r1, 'w')
    f_out_r2 = open(out_r2, 'w')

    while True:
        cur_r1_name = f_r1.readline().strip()[1:]
        cur_r1_read = f_r1.readline().strip()
        cur_r1_plus = f_r1.readline().strip()
        cur_r1_qual = f_r1.readline().strip()
        
        cur_r2_name = f_r2.readline().strip()[1:]
        cur_r2_read = f_r2.readline().strip()
        cur_r2_plus = f_r2.readline().strip()
        cur_r2_qual = f_r2.readline().strip()
    
        if cur_r1_name == "" : break
        
        barcode1 = cur_r1_read[36:42] ## 6bp tn5 index sequence
        barcode2 = cur_r1_read[0:6] ## 6bp tn5 index sequence
        barcode3 = cur_r2_read[8:14] ## 6bp tn5 index sequence
        UMI = cur_r2_read[0:8] ## 6bp tn5 index sequence
        # i5_i7 = (cur_r2_read[154:170])  ## current_150bp_formate, extract 154-170 total 16bp barcode
        #i5_i7 = cur_r2_read[100:106] + cur_r2_read[-10:] ## first 6 bp + last 10 bp
        cell_barcode = barcode1 + barcode2 + barcode3
        
        cur_r1_read = cur_r1_read[42:] ## trim first 31 bp adaptor sequence
        cur_r2_read = cur_r2_read[31:] #+ tn5_index_i5 + tn5_index_i7  ## the i5 and i7 index is used for plate demutiplex.
        
        cur_r1_qual = cur_r1_qual[42:] # + cur_r1_qual[0:6]+ cur_r2_qual[0:6] ## update quality
        cur_r2_qual = cur_r2_qual[31:]
        
        cur_r2_name = cur_r2_name.replace("/2","/1")
        
        cur_r1_name = ('@'+cell_barcode + ':' + UMI + ':' + cur_r1_name + ' '  + ' 1:' ) # + tn5_index_i5+tn5_index_i7)
        cur_r2_name = ('@'+cell_barcode + ':' + UMI + ':' + cur_r2_name + ' '  + ' 2:' ) #+ tn5_index_i5+tn5_index_i7)
        
        
        f_out_r1.write(cur_r1_name+"\n")
        f_out_r1.write(cur_r1_read+"\n")
        f_out_r1.write(cur_r1_plus+"\n")
        f_out_r1.write(cur_r1_qual+"\n")     
    
        f_out_r2.write(cur_r2_name+"\n")
        f_out_r2.write(cur_r2_read+"\n")
        f_out_r2.write(cur_r2_plus+"\n")
        f_out_r2.write(cur_r2_qual+"\n")
        
    f_r1.close()
    f_r2.close()
    f_out_r1.close()
    f_out_r2.close()
    
r1 = str(snakemake.input['r1'])
r2 = str(snakemake.input['r2'])
out_r1 = str(snakemake.output['r1'])
out_r2 = str(snakemake.output['r2'])

update_fastq(r1,r2, out_r1,out_r2 )
