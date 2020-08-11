#! /usr/bin/env python

import sys, os
from subprocess import Popen, PIPE
import pysam

from .logger import get_logger
logger = get_logger(__name__)

def is_exists(input_file):
    
    if not os.path.exists(input_file):
        logger.error("Input not exists: %s" % input_file)
        sys.exit(1)


def bam_format_check(bam_file):

    try:
        bam_t = pysam.AlignmentFile(bam_file)
    except:
        logger.error("BAM format error: %s" % bam_file)
        sys.exit(1)

def fasta_format_check(fasta_file):

    try:
        fasta_t = pysam.FastaFile(fasta_file)
    except:
        logger.error("FASTA format error: %s" % fasta_file)
        sys.exit(1)
       
 
def is_tool(executable):

    from shutil import which
    if which(executable) is None:
        logger.error("Executable does not exist: " + executable)
        sys.exit(1) 

    return True


def libssw_check():

    # check whether libssw.so is in LD_LIBRARY_PATH
    sLibPath = ""
    for ld_path in os.environ["LD_LIBRARY_PATH"].split(':'):
        # print ld_path
        if os.path.exists(ld_path + "/libssw.so"):
            sLibPath = ld_path # + "/libssw.so"
            break
    if sLibPath == "":
        logger.error("Cannot find libssw.so in LD_LIBRARY_PATH")
        sys.exit(1)

