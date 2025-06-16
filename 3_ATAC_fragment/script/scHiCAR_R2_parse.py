import sys
import gzip
import pysam
import os
import collections 
import time

"""
This script is adapted from https://github.com/r3fang/snATAC/blob/master/bin/snATAC_pre
"""


def is_sorted_queryname(header):
    """
    Check if bam fiel is sorted by read name.
    """
    if("HD" in header):
        if("SO" in header["HD"]):
            if(header["HD"]["SO"] == "queryname"):
                return True
    return False

class qc(object):
    """A quality control object that has the following attributes:
    Attributes:
        total: total number of sequenced fragments.
        mapped: number of mappable fragments.
        chrM: number of fragments mapped to chrM.
        paired: number of fragments are paired.
        single: number of fragments are from single read.
        proper_paired: number of paired reads are properly paired.
        usable: number of usable fragments.
        uniq: number of unique fragments.
        isize: average insert size distribution.
    """
    def __init__(self):
        """Return a qc object"""
        self.id = 0
        self.total = 0
        self.mapped = 0
        self.paired = 0
        self.proper_paired = 0
        self.singlemap = 0
        self.proper_flen = 0
        self.usable = 0
        self.uniq = 0
        self.final = 0

class fragment(object):
    """A fragment object that has the following attributes:
    Attributes:
        chrom: chromsome name
        start: start position
        end: end position
        mapq: mapping quality
        is_proper_pair: whether properly paired
        is_single: whether it is a single read
        is_secondary: whether it is a secondary alignment read
        flen: fragment length
    """
    def __init__(self, qname, chrom, pos, flen, mapq, is_proper_pair, strand):
        """Return a qc object"""
        self.qname = qname
        self.chrom = chrom
        self.pos = pos
        self.flen = flen
        self.mapq = mapq
        self.is_proper_pair = is_proper_pair
        self.strand = strand


def group_reads_by_barcode_bam(input_bam ,coverage):
    """ Group reads based on the barcodes
    
    Args:
        input_bam: a bam file
    Returns:
        Generator that contains reads sharing the same barcode
    """
    if not os.path.exists(input_bam): 
        print(("Error @group_reads_by_barcode_bam: " + input_bam + " does not exist!"));
 
    read_group_list = []; 
    pre_barcode = "";
    samfile = pysam.AlignmentFile(input_bam, "rb");
    if not is_sorted_queryname(samfile.header):
    	raise ValueError('need sort bam by reads')
    for cur_read in samfile:
#      if cur_read.is_read1:  
        cur_barcode = cur_read.qname.split(":")[0];
        if cur_barcode == pre_barcode:
            read_group_list.append(cur_read)
        else:
            if pre_barcode != "":
                if len(read_group_list) > coverage: ## two reads per fragment
                	yield (x for x in read_group_list) # 
            read_group_list = [cur_read] # add the first read
            pre_barcode = cur_barcode
    yield (x for x in read_group_list)  # reads from the last barcode

def read1ByName(read_list):
    """ Pair reads based on read names
    
    Args:
        read_list: a list of reads that share the same barcode
    Returns:
        Generator contains read pairs from the same fragment
        and a boolen variable indicates whether it is supplementary alignment
    """
    # modify to only report the paired reads
    for read1 in read_list:
        # read until read1 is not a secondary alignment
        while read1.qstart>0 or read1.is_secondary:
            # yield (read1, None, False, True)
            # print( 'read1 pass')
            try:
                #print "Warning: skipping ", read1.qname;
                read1 = next(read_list);
            except:
            	break
        yield (read1)  # is pair and both not supplementary


def read1ToFragment(read1):
    """ convert read pairs to fragments
    
    Args:
        read1: R1 read
        read2: R2 read
    Returns:
        Generator that contains a fragment object
    """
    try:
        read1.qname;
    except ValueError as e:
        sys.exit('read_pair_to_fragment: can not get read1 or read2 name!');   
    barcode = read1.qname.split(":")[0];
    mapq = read1.mapq;
    try:
        # strand1 = "-" if read1.is_reverse else "+";    
        chrom2 = read1.reference_name;
        start2 = read1.reference_start;
        strand2 = "-" if read1.is_reverse else "+";    
        # it is possible that flen1 is None  

        flen2 = read1.reference_length if read1.reference_length != None else 0
        end2   = read1.reference_end;
        return fragment(read1.qname, chrom2, start2, flen2, mapq, False, strand2);
    except ValueError as e:
        return fragment(read1.qname, None,   None,   None,   mapq, False);


def main():

    from argparse import ArgumentParser
    # parameters
    
    parser = ArgumentParser(description='snATAC-seq preprocessing')
    parser.add_argument('-i', '--input', help='input bam file', required=True)
    parser.add_argument('-o', '--output', help='output bed/bed.gz file', required=True)
    parser.add_argument('-m', '--mapq', help='min mappability score [10]', required=True)
    parser.add_argument('-f', '--flen', help='maximum fragment length [2000]', required=True)
    parser.add_argument('-l', '--log', help='logfile', required=True)
    parser.add_argument('-c', '--cov', help='minimal_coverage', required=True)
    options = parser.parse_args()		

 
    # input parsing
    input_bam = options.input
    output_bed = options.output
    min_mapq = int(options.mapq)
    max_flen = int(options.flen)
    coverage = int(options.cov)
    fout = open(output_bed, "w")
    log_fle = options.log

    min_mapq = 20
    max_flen = 1200
    
    check_point = 0;
    start_time = time.time();
    qc_dict = collections.defaultdict(qc);   
    # num_barcode = 0;
    print('start...')
    # print (time.time())
    
    for read_group in group_reads_by_barcode_bam(input_bam,coverage):
        frag_list = [];    
        for read1 in read1ByName(read_group):   
            frag = read1ToFragment(read1);       
            # extract the barcode
            barcode = frag.qname.split(":")[0].upper();
            umi = frag.qname.split(":")[1].upper();
            umi = umi[0:1] + umi[6:7];
            # only for printing the progress
            check_point += 1;
            if check_point%1000000 == 0:
                print(("%d\t tags, %s seconds " % (check_point, time.time() - start_time)));  
            # total number of sequencing fragments (exclude supplementary alignments)
            qc_dict[barcode].total += 1;
            ## 1. Filter non-uniquely mapped fragments
            if read1.mapq < min_mapq:
                continue; # filter out low mapq          
            qc_dict[barcode].mapped += 1;
            # 2. if it is single read, keep it only if keep_unpaired is true
            #    if it is paired, keep it only if it is approperly paired             
            # 3. check fragment size
            if frag.flen < max_flen:
                qc_dict[barcode].proper_flen += 1;
                qc_dict[barcode].singlemap += 1;
                #continue
            # 4. combine single and paired as fragments
            frag_list.append((frag.chrom, str(frag.pos), str(frag.pos+frag.flen), barcode, umi, frag.strand));              
        
        # 5. remove duplicate fragments
        # num_barcode += 1
        qc_dict[barcode].usable = len(frag_list);
        frag_list_uniq = set(frag_list); # remove duplicated fragments
        qc_dict[barcode].uniq = len(frag_list_uniq);

        for item in frag_list_uniq: ## write the current barcode files
                fout.write("\t".join(list(item))+"\n")
        del frag_list, frag_list_uniq;

    fout.close()
                
    with open(log_fle, "w") as fout:
        # fout.write("Total number of unique barcodes:             %d\n"%num_barcode)
        fout.write("Total number of barcodes writed:             %d\n"%len(qc_dict))
        fout.write("TN - Total number of fragments:              %d\n"%sum([qc_dict[key].total for key in qc_dict]))
        fout.write("UM - Total number of mapped:                 %d\n"%sum([qc_dict[key].mapped for key in qc_dict]))
        fout.write("PP - Total number of singlemap:              %d\n"%sum([qc_dict[key].singlemap for key in qc_dict]))
        fout.write("PL - Total number of proper frag len:        %d\n"%sum([qc_dict[key].proper_flen for key in qc_dict]))
        fout.write("US - Total number of usable fragments:       %d\n"%sum([qc_dict[key].usable for key in qc_dict]))
        fout.write("UQ - Total number of unique fragments:       %d\n"%sum([qc_dict[key].uniq for key in qc_dict]))
    return 0

if __name__ == '__main__':
    main()

    
