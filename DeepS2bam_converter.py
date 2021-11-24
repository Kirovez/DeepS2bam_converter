from Bio import SeqIO
import pandas as pd
import pysam
import numpy as np
import os
import re
import sys

class methylPoint():
    def __init__(self, chrom, pos, strand, pos_in_strand, readname, read_strand, prob_0, prob_1, called_label, k_mer):
        self.chrom, self.pos, self.strand, self.pos_in_strand, self.readname, self.read_strand, self.prob_0, self.prob_1, self.called_label, self.k_mer = chrom, pos, strand, pos_in_strand, readname, read_strand, prob_0, prob_1, called_label, k_mer


def deepsignalTab2dic(deepsignal_tab):
    dsp_dic = {} # read_id:reference:[sorted(methylated_signals)]
    dsp_dic_length = {} # read_id:[methylation positions]
    minus_strand = 0
    plus_strand = 0
    with open(deepsignal_tab) as inDsp:
        for lines in inDsp:
            chrom, pos, strand, pos_in_strand, readname, read_strand, prob_0, prob_1, called_label, k_mer = lines.rstrip().split('\t')
            if called_label == '1':# and strand == "+": ##methylated only
                if readname not in dsp_dic:
                    dsp_dic[readname] = {}
                    dsp_dic_length[readname] = []
                dsp_dic_length[readname].append(pos)
                if chrom not in dsp_dic[readname]:
                    dsp_dic[readname][chrom] = []
                dsp_dic[readname][chrom].append(methylPoint(chrom, pos, strand, pos_in_strand, readname, read_strand, prob_0, prob_1, called_label, k_mer))
                if strand == "+":
                    plus_strand +=1
                else:
                    minus_strand +=1
    print("Total number of Cme detected on plus reads: {0} and minus reads: {1}".format(plus_strand, minus_strand))
    print("Total number of reads possesing Cme signals: {}".format(len(dsp_dic)))

    return [dsp_dic,dsp_dic_length]


def getPerChromosomeAlgns(bam_file):
    ref_seq_per_algns = {}  # chromosome: [algnsegment1,.....]
    cnt = 0
    for num_al, algns in enumerate(pysam.AlignmentFile(bam_file, 'rb')):
        if not algns.is_supplementary and not algns.is_secondary:
            chrom = algns.reference_name
            if chrom not in ref_seq_per_algns:
                ref_seq_per_algns[chrom] = []
            ref_seq_per_algns[chrom].append(algns)
            cnt += 1

    print("Number of alignment selected", cnt)
    print("\n".join([str(len(ref_seq_per_algns[i])) for i in ref_seq_per_algns]))
    return (ref_seq_per_algns)


## 4. collect reference sequence
def replace_reference_sequence(genome_fasta, ref_seq_per_algns, dsp_dic):
    new_ref_seq_per_algns = {}
    cnt_reads_not_in_depstable = 0
    for seq in SeqIO.parse(genome_fasta, 'fasta'):
        if seq.id in ref_seq_per_algns:
            if seq.id not in new_ref_seq_per_algns:
                new_ref_seq_per_algns[seq.id] = []
            seq_str = str(seq.seq)
            print(seq.id)
            print(len(ref_seq_per_algns[seq.id]))
            for algns in ref_seq_per_algns[seq.id]:
                if algns.query_name in dsp_dic:
                    if seq.id in dsp_dic[algns.query_name]:
                        poss = [int(p.pos) for p in dsp_dic[algns.query_name][seq.id]]
                        ref_star = min(poss)  # algns.reference_start
                        ref_end = max(poss)  # algns.reference_end
                        # if ref_star < algns.reference_start:
                        algns.reference_start = ref_star  # take the most left position
                        # if ref_end < algns.reference_end:
                        #     ref_end = algns.reference_end

                        seq_new_str = seq_str[int(ref_star):int(ref_end) + 1]
                        #print(int(ref_star),int(ref_end))
                        # if len(seq_new_str) > len(algns.query_sequence):
                        algns.query_sequence = seq_new_str
                        algns.cigar = ((0, len(algns.query_sequence)),)

                        new_ref_seq_per_algns[seq.id].append(algns)
                        # print(algns.reference_end, ref_end)
                        # algns.reference_end = int(ref_end)
                else:
                    cnt_reads_not_in_depstable += 1
    print('Number of reads that were not found in deepsignal-plant table', cnt_reads_not_in_depstable)
    return new_ref_seq_per_algns


def _getMmTag(algn_old, dsp_dic):
    algn = algn_old

    cnt_not_found_cme = 0
    cnt_passed_cme = 0
    ### find all c IN REFERENCE SEQUENCE
    ref = algn.query_sequence
    target_read = algn.query_name
    target_chrom = algn.reference_name
    if dsp_dic[target_read][target_chrom][0].strand == '-':
        all_C_in_reference = [c.start() + algn.reference_start for c in re.finditer('G', str(ref).upper())]
        # print(all_C_in_reference)
        # print(target_chrom)
        # print(all_C_in_reference)
    else:
        all_C_in_reference = [c.start() + algn.reference_start for c in re.finditer('C', str(ref).upper())]

    # 4. determine number of methylated positions C
    cme_C_in_reference = []
    ml_C = []

    for cme in dsp_dic[target_read][target_chrom]:
        # print(cme.strand)
        # print(cme.k_mer)
        # print('reference', str(ref).upper()[int(cme.pos) - 5 - algn.reference_start:int(cme.pos) + 10 - algn.reference_start])
        try:
            cme_C_in_reference.append(all_C_in_reference.index(int(cme.pos)) + 1)
            cnt_passed_cme += 1
            ml_C.append(str(round(float(cme.prob_1) * 255)))
        except:
            if cme.chrom == 'NC_003076.8':
                print(cme.chrom, cme.pos, cme.strand, cme.k_mer)
            # print(int(cme.pos) > all_C_in_reference[-1])
            cnt_not_found_cme += 1


    ## write Mm string
    # if algn.is_reverse:
    #     mm_str = ['C-m']
    # else:
    mm_str = ['C+m']

    cme_C_in_reference = sorted(cme_C_in_reference)
    for i, cme in enumerate(cme_C_in_reference):
        if i == 0:
            # if algn.is_reverse:
            #     mm_str.append(str(cme + 1))
            # else:
                mm_str.append(str(cme - 1))
        else:
            mm_str.append(str(cme_C_in_reference[i] - cme_C_in_reference[i - 1] - 1))
            # if algn.is_reverse:
            #     mm_str.append(str(cme_C_in_reference[i] - cme_C_in_reference[i - 1]))

    len_mm = len(mm_str) - 1
    if algn.is_reverse:
        mm_str[-1] = str(int(mm_str[-1]) + 1)
        mm_str = ['C+m'] + mm_str[1:][::-1]
        mm_str_toret = ",".join(mm_str)
    else:
        mm_str_toret = ",".join(mm_str)
    # if algn.is_reverse:
    #     mm_str_toret += ';C-h;'
    # else:
    mm_str_toret += ';C+h;'


    Ml = ",".join(ml_C)
    algn.tags += (('Mm', mm_str_toret), ('Ml', Ml))
    #algn.flag = 0
    return [algn, cnt_not_found_cme, cnt_passed_cme]


def writeBam(modified_algns, header, bam_name = "tpm.bam"):
    with pysam.AlignmentFile(bam_name, "wb", header=header) as outf:
        for algn in modified_algns:
            outf.write(algn)

def main(bam_file, genome_fasta, deepsignal_tab):
    # main(genome_fasta, deepsignal_tab, bamfile)
    out_bam_name = deepsignal_tab + ".bam"

    ## 2.convert deep signal tab into dictionary
    dsp_dic = deepsignalTab2dic(deepsignal_tab)  # [dsp_dic(read_id:reference:[sorted(methylated_signals)] , []]
    print("DeepSignal-plant table was parsed and dictionary was obtained")

    ## 3. parse bam file
    ref_seq_per_algns = getPerChromosomeAlgns(bam_file)  # chromosome: [algnsegment1,.....]
    print("Alignments have been collected")

    # 4. replace SEQ by reference seq in each segment and change start and end coordinated to min and max (for each read in dsp table)
    ref_seq_per_algns_replaced = replace_reference_sequence(genome_fasta, ref_seq_per_algns, dsp_dic[0])
    print("New SEQ was added to the alignments from sam")

    # 5. add Mm and Ml
    reads_wo_meth = 0
    reads_wi_meth = 0
    cnt_not_found_cme, cnt_passed_cme = {True:0, False:0}, 0
    modified_algns = []
    print([len(ref_seq_per_algns_replaced[i]) for i in ref_seq_per_algns_replaced])
    for chrom in ref_seq_per_algns_replaced:
        print(chrom)
        for i, algns in enumerate(ref_seq_per_algns_replaced[chrom]):
            readname = algns.query_name
            if i%1000 == 0:
                print(i, 'alignment checked')
            if readname in dsp_dic[0]:
                if chrom in dsp_dic[0][readname]:
                    reads_wi_meth += 1
                    algns, cnf, cp = _getMmTag(algns, dsp_dic[0])
                    cnt_not_found_cme[algns.is_reverse] += cnf
                    cnt_passed_cme += cp
                    modified_algns.append(algns)
            else:
                reads_wo_meth += 1
    print("Number of reads with detected methylation", reads_wi_meth)
    print("Number of reads without detected methylation", reads_wo_meth)
    print('Total number of Cm positions mapped', cnt_passed_cme)
    print('Number of Cm positions not found in reference for - strand:', cnt_not_found_cme[True], 'for + strand:', cnt_not_found_cme[False])

    writeBam(modified_algns, pysam.AlignmentFile(bam_file, 'rb').header, bam_name = out_bam_name)
    os.system('samtools sort -@ 100 {0} > {1}'.format(out_bam_name, out_bam_name + "_sorted.bam"))
    os.system('samtools index -@ 100 {}'.format(out_bam_name + "_sorted.bam"))

if __name__ == '__main__':
    genome_fasta = sys.argv[1]
    deepsignal_tab = sys.argv[2]
    bam_file = sys.argv[3]
    main(bam_file, genome_fasta, deepsignal_tab)

