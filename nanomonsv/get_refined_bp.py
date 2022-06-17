#! /usr/bin/env python3

import sys, os, statistics 
import pysam

from . import smith_waterman
from .my_seq import reverse_complement
from .logger import get_logger

logger = get_logger(__name__)

#def get_refined_bp(contig, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, h_log, rd_margin = 20, i_margin = 500):
def get_refined_bp(contig, reference_fasta, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, output_file, rd_margin, i_margin):
    
    fasta_file_ins = pysam.FastaFile(reference_fasta)
    
    start1 = max(1, start1)
    start2 = max(1, start2)

    if mode != "i":

        ref_len1 = fasta_file_ins.get_reference_length(chr1)
        ref_len2 = fasta_file_ins.get_reference_length(chr2)

        bstart1, bend1 = max(int(start1) - rd_margin, 1), int(end1) + rd_margin
        bstart2, bend2 = max(int(start2) - rd_margin, 1), int(end2) + rd_margin

        if ref_len1 < bend1: bend1 = ref_len1
        if ref_len2 < bend2: bend2 = ref_len2    

        region1_seq = fasta_file_ins.fetch(chr1, bstart1 - 1, bend1)
        region2_seq = fasta_file_ins.fetch(chr2, bstart2 - 1, bend2)

        if dir1 == '-': region1_seq = reverse_complement(region1_seq)
        if dir2 == '+': region2_seq = reverse_complement(region2_seq)

        sret = smith_waterman.sw_jump(contig, region1_seq, region2_seq)
        if sret is None:
            with open(output_file, "w") as f:
                f.write("")
            return(None)
        score, contig_align, region1_align, region2_align, contig_seq, region_seq = sret

        bp_pos1 = bstart1 + region1_align[1] - 1 if dir1 == '+' else bend1 - region1_align[1] + 1 
        bp_pos2 = bstart2 + region2_align[0] - 1 if dir2 == '-' else bend2 - region2_align[0] + 1

        if contig_align[2] - contig_align[1] == 1:
            inseq = '---'
        elif contig_align[2] - contig_align[1] > 1:
            inseq = contig[(contig_align[1]):(contig_align[2] - 1)]
            if dir1 == '-': inseq = reverse_complement(inseq)
        else:
            logger.warning("Alignment inconsistent!!")

        print(score, contig_align, region1_align, region2_align)
        print(contig_seq)
        print(region_seq)

        with open(output_file, "w") as f:
            f.write("{pos1}\t{pos2}\t{inseq}".format(pos1 = bp_pos1, pos2 = bp_pos2, inseq = inseq))
    else:
    
        contig_start = contig[:min(i_margin, len(contig))]
        contig_end = contig[-min(i_margin, len(contig)):]

        region_seq = fasta_file_ins.fetch(chr1, max(0, start1 - 100), end2 + 100)
        sret = smith_waterman.sw_jump(region_seq, contig_start, contig_end)
        if sret is None:
            with open(output_file, "w") as f:
                f.write("")
            return(None)
        score, region_align, contig_start_align, contig_end_align, region_seq, contig_seq = sret

        bp_pos1 = max(0, start1 - 100) + region_align[1]
        bp_pos2 = max(0, start1 - 100) + region_align[2]

        inseq_start = contig_start_align[1]
        inseq_end = len(contig) - (len(contig_end) - contig_end_align[0] + 1)
        inseq = contig[inseq_start:(inseq_end + 1)]

        print(score, region_align, contig_start_align, contig_end_align)
        print(region_seq)
        print(contig_seq)

        with open(output_file, "w") as f:
            f.write("{pos1}\t{pos2}\t{inseq}".format(pos1 = bp_pos1, pos2 = bp_pos2, inseq = inseq))
        
def get_refined_bp_IF(args):
    get_refined_bp(args.contig, args.reference_fasta, args.chr1, args.start1, args.end1, args.dir1, args.chr2, args.start2, args.end2, args.dir2, args.mode, args.output_file, args.rd_margin, args.i_margin)

def main():
    import argparse
    from .version import __version__

    parser = argparse.ArgumentParser(prog = "nanomonsv_get_refined_bp")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    parser.add_argument("contig", default = None, type = str, help = "")
    parser.add_argument("reference_fasta", type = str, help = "")
    parser.add_argument("chr1", type = str, help = "")
    parser.add_argument("start1", type = int, help = "")
    parser.add_argument("end1", type = int, help = "")
    parser.add_argument("dir1", type = str, help = "")
    parser.add_argument("chr2", type = str, help = "")
    parser.add_argument("start2", type = int, help = "")
    parser.add_argument("end2", type = int, help = "")
    parser.add_argument("dir2", type = str, help = "")
    parser.add_argument("mode", type = str, help = "")
    parser.add_argument("output_file", type = str, help = "")
    parser.add_argument("--rd_margin", type = int, help = "", default = 20)
    parser.add_argument("--i_margin", type = int, help = "", default=500)

    parser.set_defaults(func = get_refined_bp_IF)

    args = parser.parse_args()
    args.func(args)
