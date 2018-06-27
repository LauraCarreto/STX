#!/usr/bin/env python
'''
This script runs the stx subtyping assay for a given ecoli sample.
'''

import sys
import os
import argparse
import json
import glob
from lxml import etree
import tempfile
import subprocess
import shutil
import pysam
import fileinput


__version__ = '0.1'
__date__ = '09Aug2016'
__author__ = 'ulf.schaefer@phe.gov.uk'

COMPONENTNAME = "stx_subtyping"
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
                        help="The folder where the procesed fastqs are found.")

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
    sResFolder = os.path.join(oArgs.input, COMPONENTNAME)
    try:
        os.mkdir(sResFolder)
    except OSError:
        bDirExisted = True

    dResult = {}

    aFastqFileName = glob.glob(oArgs.input + "/*.processed.R*.fastq*")
    print "aFastqFileName", aFastqFileName

    aFastqFileName.sort()
    sFileIn1 = ""
    sFileIn2 = ""
    sFileIn1 = aFastqFileName[0]
    sFileIn2 = aFastqFileName[1]
    aTmp = [x.strip() for x in os.path.basename(sFileIn1).split(".")]
    [sNGSSampleID, sWorkflow] = [".".join(aTmp[:-6]), ".".join(aTmp[-6:-4])]
    dResult['samid'] = sNGSSampleID
    dResult['wf'] = aTmp[-6]
    dResult['wf_version'] = aTmp[-5]


    sTmpFolder = oArgs.input
    bowtie2_output_dir = os.path.join(sTmpFolder, 'bowtie2_out')

    if not os.path.exists(bowtie2_output_dir):
        os.makedirs(bowtie2_output_dir)

    ref = oArgs.ref + "/stx_genes_with_flanking_regions.fa"

    #check_reference(ref)

    #map_to_stx(ref, sFileIn1, sFileIn2, bowtie2_output_dir, sNGSSampleID)
    sam_to_bam(bowtie2_output_dir, sNGSSampleID)

    results_dir = os.path.join(sTmpFolder, 'results')

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    dResult['matches'] = run_pileup_analysis(ref, bowtie2_output_dir, sNGSSampleID, results_dir)

    # write results xml file
    oXMLRoot = etree.Element("ngs_sample", id=dResult['samid'])
    etree.SubElement(oXMLRoot,
                     "workflow",
                     value=dResult['wf'],
                     version=dResult['wf_version'])
    oXMLResultsNode = etree.SubElement(oXMLRoot, "results")
    if len(dResult['matches']) <= 0:
            oXMLNode = etree.SubElement(oXMLResultsNode,
                                        COMPONENTNAME,
                                        success='false',
                                        reason='no matches found')
    else:
        oXMLNode = etree.SubElement(oXMLResultsNode, COMPONENTNAME, success='true')
        for match in dResult['matches']:
            _ = etree.SubElement(oXMLNode, "stx", value=match)
    with open(os.path.join(sResFolder, "%s.results.xml" % (sNGSSampleID)), 'w') as fOutFile:
        fOutFile.write(etree.tostring(oXMLRoot, pretty_print=True))

    return 0

# end of main --------------------------------------------------------------------------------------

def check_reference(ref_file):
    '''
    Checks whether the reference is present and if it has been indexed with bwa.
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

    sCmd= "bowtie2-build %s %s" % (ref_file, ref_file)
    print sCmd
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()


# --------------------------------------------------------------------------------------------------

def check_reference_prepared(ref):
    '''
    Checks if reference is prepared.

    Parameters
    ----------
    ref: str
        absolute path to reference fasta file

    Returns
    -------
    True if all files are present, else False
    '''

    return all([os.path.exists(ref + '.' + i) for i in ['amb', 'ann', 'bwt', 'pac', 'sa']])

# --------------------------------------------------------------------------------------------------

def run_pileup_analysis(ref, bowtie2_output_dir, s_name, results_dir):
    '''
    Runs the pipeup analysis and returns all matched stx genes.

    Parameters
    ----------
    ref: str
        reference genome file
    bowtie2_output_dir: str
        outout directory for bwa
    s_name: str
        sample name
    results_dir: str
        results folder
    logger: obj
        logger object

    Returns
    -------
    aMatchedStxs: list
        a list of all the stx genes that were matched
    '''

    bamfile = '%s/%s.unique.sorted.bam' % (bowtie2_output_dir, s_name)

    sCmd = "cat %s | grep '^>' | sed 's/^>//g'" % (ref)
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    references = p_stdout.strip().split("\n")


    """
    references = ['stx2a_gi|14892|emb|X07865.1|',
                  'stx2c_gi|15718404|dbj|AB071845.1|',
                  'stx2d_gi|30313370|gb|AY095209.1|',
                  'stx2b_gi|49089|emb|X65949.1|',
                  'stx2e_gi|8346567|emb|AJ249351.2|',
                  'stx2f_gi|254939478|dbj|AB472687.1|',
                  'stx2g_gi|30909079|gb|AY286000.1|',
                  'stx1a_gi|147832|gb|L04539.1|',
                  'stx1c_gi|535088|emb|Z36901.1|',
                  'stx1d_gi|28192582|gb|AY170851.1|']
                  """

    aMatchedStxs = []

    for reference in references:
        #print reference
        mapped_bases = make_pile_up(bamfile, reference)
        matched_stx = parse_pileup(mapped_bases, reference)
        if matched_stx is not None:
            aMatchedStxs.append(matched_stx.split('_')[0])

    return aMatchedStxs

# --------------------------------------------------------------------------------------------------

def make_pile_up(bamfile, reference):
    '''
    Makes a pileup dictionary from the bamfile and the reference.

    Parameters
    ----------
    bamfile: str
        path to bam file to analyse
    reference: str
        path to reference file

    Returns
    -------
    mapped_bases: dict
        maps a position to a list of bases mapped there
    '''

    bamfile = pysam.Samfile(bamfile, "rb")
    mapped_bases = {}
    for pileupcolumn in bamfile.pileup(reference):
        #print pileupcolumn
        for pileupread in pileupcolumn.pileups:
            #print pileupread
            #print '\tbase in base %s = %s' % (pileupcolumn.pos, pileupread.alignment.seq[pileupread.qpos])
            pos = pileupcolumn.pos
            print pileupread.query_position

            if pileupread.query_position:
                base = pileupread.alignment.seq[pileupread.query_position].split()
                if pos in mapped_bases:
                    mapped_bases[pos].append(base[0])
                else:
                    mapped_bases[pos] = base
    #for each in mapped_bases:
    #    print each, mapped_bases[each]=
    return mapped_bases

# --------------------------------------------------------------------------------------------------

def parse_pileup(mapped_bases, reference):
    '''
    Makes a pileup dictionary from the bamfile and the reference.

    Parameters
    ----------
    mapped_bases: dict
        as returned by make_pile_up()
    reference: str
        path to reference file

    Returns
    -------
    if found:
    reference: str
        the gene name
    else
    None
    '''

    uniq_pos = ref_to_uniq_pos_dict[reference]
    i = 0
    for each in mapped_bases:
        #print each, mapped_bases[each]
        if each in uniq_pos:

            depth = len(mapped_bases[each])
            uniq_mapping = mapped_bases[each].count(uniq_pos[each])
            if depth > 10:
                #print uniq_mapping / depth
                if uniq_mapping / float(depth) >= 0.9:
                    i += 1
                    #print each, mapped_bases[each]

    if i == 3:
        return reference
        #print '%s is present in sample' % reference

            #print each, mapped_bases[each]
    else:
        return None

# --------------------------------------------------------------------------------------------------

def sam_to_bam(bowtie2_output_dir, s_name):
    '''
    Converts a sam to a sorted indexed bam file.

    Parameters
    ----------
    bowtie2_output_dir:
        output dir for bwa
    s_name: str
        sample name
    logger: obj
        logger object

    Returns
    -------
    no returns
    '''

    #
    #sCmd = 'samtools view -h -F 4 -bS -o %s/%s.unique.bam %s/%s.sam' \
    #       % (bowtie2_output_dir, s_name, bowtie2_output_dir, s_name)
    #print "sCmd", sCmd
    #
    #
    #p = subprocess.Popen(sCmd, shell=True, stdin=None,
    #                     stdout=subprocess.PIPE,
    #                     stderr=subprocess.PIPE, close_fds=True)
    #(p_stdout, p_stderr) = p.communicate()
    #
    sCmd ='samtools sort %s/%s.unique.bam -o %s/%s.unique.sorted.bam' \
          % (bowtie2_output_dir, s_name, bowtie2_output_dir, s_name)

    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    sCmd = 'samtools index %s/%s.unique.sorted.bam' % (bowtie2_output_dir, s_name)
    print sCmd
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    if os.path.exists('%s/%s.unique.sorted.bam' % (bowtie2_output_dir, s_name)):
        os.remove('%s/%s.sam' % (bowtie2_output_dir, s_name))
        os.remove('%s/%s.unique.bam' % (bowtie2_output_dir, s_name))

    return

# ---------------------------------------------------------------------------------------------------

def map_to_stx(ref, fastq_read1, fastq_read2, bowtie2_output_dir, s_name):#, logger):
    '''
    Uses Bowtie2 mapper to map the fastqs to the reference in a subprocess.

    Parameters
    ----------
    ref: str
        reference genome fasta file
    fastq_read1: str
        r1 fastq file (abs path to)
    fastq_read2: str
        r2 fastq file (abs path to)
    bowtie2_output_dir:
        output dir for Bowtie2
    s_name: str
        sample name
    logger: obj
        logger object

    Returns
    -------
    no returns
    '''

    #sCmd = 'bowtie2 -x %s -1 %s -2 %s --fr --no-unal --minins 300 --maxins 1100 -k 99999 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 > %s/%s.tmp' % (ref, fastq_read1, fastq_read2, bowtie2_output_dir, s_name)
    sCmd = 'bowtie2 --fr --minins --no-unal 300 --maxins 1100 -k 99999 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x %s  -1 %s -2 %s -S %s/%s.tmp'\
           % (ref, fastq_read1, fastq_read2, bowtie2_output_dir, s_name)

    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, close_fds=True)

    (p_stdout, p_stderr) = p.communicate()

    tmp = os.path.join(bowtie2_output_dir, s_name +'.tmp') # temporary sam output
    #print "tmp", tmp
    sam = os.path.join(bowtie2_output_dir, s_name +  '.sam')
    #print "sam", sam

    i = open(tmp)
    o = open(sam, 'w')
    remove_secondary_mapping_bit(tmp, sam)
    i.close()
    o.close()

    return

# ---------------------------------------------------------------------------------------------------
def remove_secondary_mapping_bit(sam, sam_parsed):

    """
    Function
    Takes a SAM file and deducts 256 from the second column(FLAG) that unset the
    secondary alignment bit score
    NB: reads with bit(250) set are not reported when using Samtools pileup

    The option for method:
    sam[string]: SAM file
    sam_parsed[string]: parsed SAM file
    """
    lines = iter(fileinput.input([sam]))
    sam_parsed_file = open(sam_parsed, "w")
    headers = []
    body = []

    for line in lines:
        if line.startswith('@'):
            sam_parsed_file.write(line)
        else:
            # chomp line
            line = line.rstrip('\n')
            details = line.split("\t")
            flag = int(details[1])
            if flag > 256:
                details[1] = str(flag - 256)
            print >> sam_parsed_file, '\t'.join(details)
    sam_parsed_file.close()

#------------------------------------------------------------------------------------------------------------

def read_config():
    '''
    Read the config set by the modulefile from the environment into a dict

    Parameters
    ----------
    no inputs

    Returns
    -------
    d: dict
        config dictionary
    '''

    d = {}

    d['stx_ref'] = os.environ['STX_REFERENCE']

    return d

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
