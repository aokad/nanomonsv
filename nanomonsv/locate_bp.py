#! /usr/bin/env python3

import sys, os, subprocess, shutil, statistics 
from .logger import get_logger

logger = get_logger(__name__)

def locate_bp(consensus_file, output_file, reference_fasta, debug):

    #fasta_file_ins = pysam.FastaFile(reference_fasta)

    with open(consensus_file, 'r') as hin, open(output_file, 'w') as hout, \
        open(output_file + ".locate_bp.log", 'w') as hout_log:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            temp_key, tconsensus = F[0], F[1]
            print(temp_key + '\n' + tconsensus, file = hout_log)
            chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, cid = temp_key.split(',')
            #start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2)       
            #bret = get_refined_bp(tconsensus, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, hout_log)
            nanomonsv_get_refined_bp_output = output_file + "." + temp_key + ".nanomonsv_get_refined_bp"
            subprocess.check_call([
                "nanomonsv_get_refined_bp", tconsensus, reference_fasta, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, nanomonsv_get_refined_bp_output
            ], stdout = hout_log, stderr = subprocess.DEVNULL)
            
            with open(nanomonsv_get_refined_bp_output) as f:
                bret = f.read().split("\t")
            os.remove(nanomonsv_get_refined_bp_output)
            if len(bret) == 3:
                bp_pos1, bp_pos2, inseq = bret 
                print(bp_pos1, bp_pos2, inseq, file = hout_log)
                print('', file = hout_log)
                print(f"{chr1}\t{bp_pos1}\t{dir1}\t{chr2}\t{bp_pos2}\t{dir2}\t{inseq}\t{mode}_{cid}", file = hout)
                # print('\t'.join([chr1, str(bp_pos1), dir1, chr2, str(bp_pos2), dir2, inseq]), file = hout)

    if not debug: os.remove(output_file + ".locate_bp.log")
