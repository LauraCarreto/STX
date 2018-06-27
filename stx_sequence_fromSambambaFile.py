#!/usr/bin/env python
'''
This script gets the sequence of a particular stx gene present in a given ecoli sample.Only one confirmed type should be present.
'''

import sys
import os
import argparse
import json
import glob
import tempfile
import subprocess
import shutil
#import pysam
import fileinput


__version__ = '0.2'
__date__ = '07April2017'
__author__ = 'adapted from ulf.schaefer@phe.gov.uk by laura.dasilvacarreto@ouh.nhs.uk'

#COMPONENTNAME = "stx_subtyping_bowtie2"
ref_to_uniq_pos_dict = {'stx2a_gi|14892|emb|X07865.1|':{1174:'G', 1179:'C', 1191:'G'},
                        'stx2c_gi|15718404|dbj|AB071845.1|':{1054:'A', 1173:'A', 1191:'A'},
                        'stx2d_gi|30313370|gb|AY095209.1|':{1037:'C', 1173:'A', 1178:'T'},
                        'stx2b_gi|49089|emb|X65949.1|':{1060:'C', 1083:'A', 1158:'C'},
                        'stx2e_gi|8346567|emb|AJ249351.2|':{966:'C', 1108:'G', 1192:'G'},
                        'stx2f_gi|254939478|dbj|AB472687.1|':{1015:'C', 964:'G', 897:'G'},
                        'stx2g_gi|30909079|gb|AY286000.1|':{914:'C', 938:'A', 1079:'C'},
                        'stx1a_gi|147832|gb|L04539.1|':{468:'A', 528:'A', 712:'G'},
                        'stx1c_gi|535088|emb|Z36901.1|':{741:'T', 905:'G', 922:'T'},
                        'stx1d_gi|28192582|gb|AY170851.1|':{478:'C', 501:'G', 508:'T'}}

# -------------------------------------------------------------------------------------------------

def parse_args():
    """
    Parse arguments

    Parameters
    ----------
    no inputs

    Returns
    -------
    oArgs: obj
        arguments object
    """

    sDescription = 'version %s, date %s, author %s' %(__version__, __date__, __author__)

    parser = argparse.ArgumentParser(description=sDescription)

    parser.add_argument("--workflow",
                        "-w",
                        type=str,
                        metavar="STRING",
                        required=True,
                        dest="workflow",
                        help="The workflow for this sample.")

    parser.add_argument("--input",
                        "-i",
                        type=str,
                        metavar="STRING",
                        dest="input",
                        required=True,
                        help="The folder where the processed fastqs are found.")

    parser.add_argument("--ref",
                        "-r",
                        type=str,
                        metavar="STRING",
                        dest="ref",
                        required=True,
                        help="The folder where the reference file is found.")

    oArgs = parser.parse_args()
    return oArgs

# -------------------------------------------------------------------------------------------------

def main():
    '''
    Main funtion, creates all logs and result files

    Parameters
    ----------
    no inputs

    Returns
    -------
    0 on success else 1
    '''

    oArgs = parse_args()

    bDirExisted = False

    dResult = {}

    bamFileName = glob.glob(oArgs.input + "/bowtie2_out/*.sorted.bam")
    #print "bamFileName", bamFileName
    aTmp = [x.strip() for x in os.path.basename(bamFileName[0]).split(".")]
    #print aTmp

    sNGSSampleID = aTmp[0]
    dResult['samid'] = sNGSSampleID
    dResult['wf'] = 'escherichia_coli_typing'
    dResult['wf_version'] = 'ngsservice'
    #print sNGSSampleID

    sTmpFolder = oArgs.input
    bowtie2_output_dir = os.path.join(sTmpFolder, 'bowtie2_out')

    if not os.path.exists(bowtie2_output_dir):
        os.makedirs(bowtie2_output_dir)

    ref = oArgs.ref + "/stx_genes_with_flanking_regions.fa"

    check_reference(ref)
    make_ref_list(ref)
    length_reference(ref)

    ##NOTE:comment out next line if running sambamba_coverage to get sambambafile
    sambambafile = '%s/%s.sambambaOut.txt' % (bowtie2_output_dir, sNGSSampleID)

    sambamba_stats(oArgs.ref, sambambafile, bowtie2_output_dir, sNGSSampleID)

    sequence_from_sambamba(sambambafile, bowtie2_output_dir, sNGSSampleID)


    return 0

# end of main --------------------------------------------------------------------------------------

def sequence_from_sambamba(sambambafile, bowtie2_output_dir, sNGSSampleID):

    '''
    Obtain the base sequence from the Sambamba depth base file

    Parameters
    ----------
    coveragefile: str (global)
        coverage file output (txt) from sambamba_coverage

    Returns
    -------
    cov_dict: dictionary
        a dictionary of the references with respective coverage
    '''

    #open sambamba file and iterate through lines
    #make dictionary with pos:base

    for reference in references:
        dCoverage = {} #dictionary of pos= [covTotal, covA, covC, covG, covT] cov= depth

        lines = iter(fileinput.input([sambambafile]))
        countLines = 0
        for line in lines:
            if line.startswith(reference):
                countLines += 1
                details = line.split("\t")
                pos = int(details[1])
                covTotal = int(details[2])
                covA = int(details[3])
                covC = int(details[4])
                covG = int(details[5])
                covT = int(details[6])

                #add position, total coverage and base coverages to dCoverage
                dCoverage[pos]=[covTotal, covA, covC, covG, covT]

        #dictionary with reference length
        lenLines = iter(fileinput.input([lengthfile]))
        for line in lenLines:
            if line.startswith(reference):
                details = line.split ("\t")
                refLen = details[1]
            #print refLen

        uniq_pos = ref_to_uniq_pos_dict[reference]
        #print uniq_pos


        dNoCoverage = {} #dictionary of pos covered less than 10x
        base_seq='' #empty string for base sequence
        #get base present more than 90% and with at least 10 reads
        for pos in range(0,int(refLen)):
            if pos in dCoverage.keys()and dCoverage[pos][0]>= 5:
                covA = dCoverage[pos][1]
                covC = dCoverage[pos][2]
                covG = dCoverage[pos][3]
                covT = dCoverage[pos][4]
                dCovlist = {'A':covA, 'C':covC,'G':covG, 'T':covT}
                max_base = max(dCovlist, key=lambda k: dCovlist[k]) #key with max reads

                #if dCovlist[max_base] >= 0.9: #90% filter NOT WORKING; check FROM HERE
                base_seq += max_base

                #else:
                    #print 'rBase', dCovlist[max_base]
                    #base_seq += 'N'

            else:
                base_seq += 'N'
                if pos in dCoverage.keys()and dCoverage[pos][0] < 10:
                    dNoCoverage[pos] = dCoverage[pos]
                else:
                    dNoCoverage[pos] = [0,0,0,0,0]
                print pos, dNoCoverage[pos]

            if pos in uniq_pos:
                print pos, covA, covC, covG, covT
                print 'uniq_pos',pos, max_base, dCoverage[pos]

        #print sequence to fasta file
        stxRefName = reference.split('_')[0]
        sampleID = sNGSSampleID.split('_')[0]
        sequence_to_fasta(stxRefName, sampleID, base_seq, bowtie2_output_dir)

    return

#---------------------------------------------------------------------------------------------------
def sequence_to_fasta(stxRefName, sampleID, base_seq, bowtie2_output_dir):

    fastaID = stxRefName+'_'+sampleID
    print 'Hello '+ fastaID
    print 'seq len', len(base_seq)
    print base_seq
    print 'N', base_seq.count('N')

    fastaFile = os.path.join(bowtie2_output_dir, fastaID+'.fa') # output file
    print 'fasta file', fastaFile

    #open output file
    openFile = open(fastaFile, 'w')
    #write line with sample name
    print >> openFile, '>'+ fastaID
    #append a sequence
    print >> openFile, base_seq

    openFile.close()

    return

#---------------------------------------------------------------------------------------------------

def sambamba_coverage(bamfile, bowtie2_output_dir, s_name):

    '''
    Create a coverage file using Sambamba

    Parameters
    ----------
    bamfile: str (global)
        BAM file
    bowtie2_output_dir: str
        outout directory for bwa
    s_name: str
        sample name

    Returns
    -------
    result: Sambamba depth base file (.txt)

    '''

    global sambambafile
    sambambafile = '%s/%s.sambambaOut.txt' % (bowtie2_output_dir, s_name)

    sCmd= "sambamba depth base %s > %s" % (bamfile, sambambafile)
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()


    return
#---------------------------------------------------------------------------------------------------
def sambamba_stats(ref_dir, sambambafile, bowtie2_output_dir, s_name):

    '''
    Obtain coverage statistics from Sambamba depth base file

    Parameters
    ----------
    sambambafile: str (global)
        coverage file output (txt) from sambamba_coverage

    Returns
    -------
    cov_dict: dictionary
        a dictionary of the references with respective coverage
    cov_stats: dictionary
        a dictionary of the references with respective coverage statistsics (min, average, max, percent coverage)
    '''

    #open coverage file and make list with values in 3rd column (COV) from lines the start with reference
    #open length_reference file and create a dictionary of references and respective length
    global lengthfile
    lengthfile = ref_dir +  "/stx_genes_with_flanking_regions.fa.fai"

    dCoverage = {} #reference = list of position coverages
    dPosition = {} #reference = list of positions
    dRefLen = {} #reference = reference length
    dCovStats ={} #reference = list with min, average, max, percent coverage

    for reference in references:
        #print reference
        covLines = iter(fileinput.input([sambambafile]))
        countLines = 0
        for line in covLines:
            if line.startswith(reference):
                countLines += 1
                details = line.split("\t")
                pos = int(details [1])
                cov = int(details[2])
                #for each reference, add position coverages to dCoverage
                if reference in dCoverage:
                    dCoverage[reference].append(cov)
                    dPosition[reference].append(pos)
                else:
                    dCoverage[reference] = [cov]
                    dPosition[reference] = [pos]
        #print reference, countLines

        #dictionary with reference length
        lenLines = iter(fileinput.input([lengthfile]))
        for line in lenLines:
            if line.startswith(reference):
                details = line.split ("\t")
                leng = details[1]
                #add length to dRefLen
                if reference in dRefLen:
                    next
                else:
                    dRefLen[reference] = int(leng)


    for reference in dCoverage:
        min_cov = min(dCoverage[reference])
        average_cov = sum(dCoverage[reference]) / dRefLen[reference]
        max_cov = max(dCoverage[reference])
        percent_cov = (countLines * 100)/ dRefLen[reference]

        dCovStats[reference]= [min_cov, average_cov, max_cov, percent_cov]

    for reference in sorted(dCovStats.keys()):
        print reference.split('_')[0], dCovStats[reference]
        print 'Average coverage', dCovStats[reference][1]
        print 'Lines in file', countLines
        print 'Reference length', dRefLen[reference]

    return

#---------------------------------------------------------------------------------------------------

def length_reference(ref):
    '''
    Uses samtools faidx to create a file summarising the length of the references

    Parameters
    ----------
    ref: str
        reference file

    Returns
    -------
    result: boolean
        True if present and already indexed or if present and indexing successful
        False if absent or indexing fails
    '''

    sCmd= "samtools faidx %s" % (ref)
    #print sCmd
    #sys.exit()
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    return

#---------------------------------------------------------------------------------------------------

def make_ref_list(ref):

    sCmd = "cat %s | grep '^>' | sed 's/^>//g'" % (ref)
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    global references
    references = p_stdout.strip().split("\n")
    references.sort()
    #print references

    return

#---------------------------------------------------------------------------------------------------
def check_reference(ref_file):
    '''
    Checks whether the reference is present and if it has been indexed with bowtie2-build.
    Indexes if it has not been.

    Parameters
    ----------
    ref_file: str
        reference file
    logger: obj
        logger object

    Returns
    -------
    result: boolean
        True if present and already indexed or if present and indexing successful
        False if absent or indexing fails

    '''

    sCmd= "/home/laura/bowtie2-2.3.1/bowtie2-build %s %s" % (ref_file, ref_file)
    #print sCmd
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()


# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    sys.exit(main())
