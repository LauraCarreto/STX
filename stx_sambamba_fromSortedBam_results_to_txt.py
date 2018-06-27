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

import numpy as np
import matplotlib.pyplot as plt


__version__ = '0.2'
__date__ = '07April2017'
__author__ = 'adapted from ulf.schaefer@phe.gov.uk by laura.dasilvacarreto@ouh.nhs.uk'

COMPONENTNAME = "stx_subtyping_bowtie2"
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
    length_reference(ref)

    coverageThresholds = [0.9, 0.4, 0.3, 0.2]

    #print to file in dir one level up the input dir - appends subtyping results from multiple samples
    allSamplesDir = os.path.normpath(os.path.join(oArgs.input, ".."))
    resultsFile = open(os.path.join(allSamplesDir, "subtyping_Allresults.txt"), 'a')
    resultsFile.write(sNGSSampleID.rstrip('\n'))

    for covThresh in coverageThresholds:

        dResult['matches'] = run_pileup_analysis(ref, bowtie2_output_dir, sNGSSampleID, covThresh)
        #print 'Coverage threshold', covThresh, 'MatchedStxs', aMatchedStx

        ##write results xml file
        #oXMLRoot = etree.Element("ngs_sample", id=dResult['samid'])
        #etree.SubElement(oXMLRoot,
                        #"workflow",
                        #value=dResult['wf'],
                        #version=dResult['wf_version'])
        #oXMLResultsNode = etree.SubElement(oXMLRoot, "results")

        if len(dResult['matches']) <= 0:
                result = 'None'
                #oXMLNode = etree.SubElement(oXMLResultsNode,
                                        #COMPONENTNAME,
                                        #success='false',
                                        #reason='no matches found')
        else:
            result = dResult['matches']
            result = ''.join(result)#map(str, result))
            result = result.replace('stx','')

            #oXMLNode = etree.SubElement(oXMLResultsNode,
                                        #COMPONENTNAME,
                                        #success='true')

            #for match in dResult['matches']:
                #_ = etree.SubElement(oXMLNode, "stx", value=match)

        #with open(os.path.join(sResFolder, "%s.results_%scov.xml" % (sNGSSampleID, covThresh)), 'w') as fOutFile:
            #fOutFile.write(etree.tostring(oXMLRoot, pretty_print=True))

        #print results to txt file
        resultsFile.write('\t'+result.rstrip('\n')) #append without changing line

    resultsFile.write('\n') #change of line for next sample
    resultsFile.close()

    sambamba_coverage (bamfile, bowtie2_output_dir,sNGSSampleID)

    ##NOTE:comment out next line if running sambamba_coverage to get sambambafile
    sambambafile = '%s/%s.sambambaOut.txt' % (bowtie2_output_dir, sNGSSampleID)

    sambamba_stats(oArgs.ref, sambambafile, bowtie2_output_dir, sNGSSampleID)

    #save to file
    averageCoverageFile = open(os.path.join(allSamplesDir, "coverage_All.txt"), 'a')
    averageCoverageFile.write(sNGSSampleID.rstrip('\n'))
    for cov in covList:
        averageCoverageFile.write('\t'+str(cov).rstrip('\n'))
        #resultsFile.write('\n') #change of line for next sample
    averageCoverageFile.write('\n') #change of line for next sample
    averageCoverageFile.close()

    return 0

# end of main --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

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

    overall_cov = 0
    type1_cov = 0
    type2_cov = 0

    r_count = 0
    t1_count = 0
    t2_count = 0

    type1 = ['stx1a', 'stx1c', 'stx1d']

    for reference in references:
        if reference in dCoverage:
            r_count +=1
            average_cov = dCovStats[reference][1]
            #calculate overall average coverage as running average
            overall_cov = round((overall_cov * (r_count-1) + average_cov) / r_count, 0)
            #print reference.split('_')[0], 'average', average_cov
            if reference.split('_')[0] in type1:
                t1_count +=1
                type1_cov = round((type1_cov * (t1_count-1) + average_cov) / t1_count, 0)
            else:
                t2_count +=1
                type2_cov = round((type2_cov * (t2_count-1) + average_cov) / t2_count, 0)

        else:
            print reference.split('_')[0], 'has no aligned reads'

    global covList # to print to file in main
    covList = [type1_cov, type2_cov, overall_cov]

    print 'overall average coverage', overall_cov
    print 'type 1 coverage', type1_cov
    print 'type 2 coverage', type2_cov

    data_cov = []
    data_pos = []
    empty = np.zeros((1,5))
    label_list= ['stx1a', 'stx1c', 'stx1d', 'stx2a', 'stx2b', 'stx2c', 'stx2d', 'stx2e', 'stx2f', 'stx2g']
    dLabel = {'stx1a':[468, 528, 712],
            'stx1c':[741, 905, 922],
            'stx1d':[478, 501, 508],
            'stx2a':[1174, 1179, 1191],
            'stx2b':[1060, 1083, 1158],
            'stx2c':[1054, 1173, 1191],
            'stx2d':[1037, 1173, 1178],
            'stx2e':[966, 1108, 1192],
            'stx2f':[897, 964, 1015],
            'stx2g':[914, 938, 1079],
            }

    for reference in references:
        #append data to plot:
        if dCoverage.has_key(reference):
            data_cov.append(dCoverage[reference])
            data_pos.append(dPosition[reference])
        else:
            data_cov.append(empty)
            data_pos.append(empty)
    #print data


    #boxplots
    fig1, ax1 = plt.subplots()
    ax1.set_title('Comparison of alignment depth across references\n'+'(overall average coverage = %s)' % overall_cov)
    ax1.set_xlabel('Reference')
    ax1.set_ylabel('Depth')
    ax1.boxplot(data_cov)
    xtickNames = plt.setp(ax1, xticklabels= label_list)
    plt.setp(xtickNames, rotation=45, fontsize=8)
    plt.subplots_adjust(bottom = 0.15)
    #plt.text(2, 0.65, 'overall average coverage = %s' % overall_cov)

    from matplotlib.backends.backend_pdf import PdfPages
    bp = PdfPages('%s/%s.boxplots.pdf' % (bowtie2_output_dir, s_name))
    bp.savefig()
    bp.close()

    #line plots
    fig2 = plt.figure()
    fig2.add_subplot(111, frameon = False)
    plt.suptitle('Alignment with different stx gene sequences')
    plt.tick_params(labelcolor = 'none', top = 'off', bottom = 'off', left = 'off', right = 'off')
    plt.xlabel('Nucleotide position')
    plt.ylabel('Depth')

    for i in range (1,11):
        ax = fig2.add_subplot(4,3,i)
        ax.plot(data_pos[i-1], data_cov[i-1], 'k')
        ax.set_title(label_list[i-1], size=8)
        #ax.set_ylim([0,400])
        ax.set_xlim([0,1500])
        uniqLineList = dLabel[label_list[i-1]]
        for line in uniqLineList:
            ax.axvline(x=line, color='r', linestyle='--')

        if i in (2,3,5,6):
            ax.tick_params(axis ='x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off' )
            ax.tick_params(axis ='y', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off' )

        if i in (1,4,7):
            ax.tick_params(axis ='x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off' )

    fig2.subplots_adjust(hspace = 0.3, wspace = 0.3)

    lp = PdfPages('%s/%s.lineplots.pdf' % (bowtie2_output_dir, s_name))
    lp.savefig()
    lp.close()


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

    return
# --------------------------------------------------------------------------------------------------

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


# --------------------------------------------------------------------------------------------------

def run_pileup_analysis(ref, bowtie2_output_dir, s_name, covThresh):
    '''
    Runs the pipeup analysis and returns all matched stx genes.

    Parameters
    ----------
    ref: str
        reference genome file
    bowtie2_output_dir: str
        output directory
    s_name: str
        sample name
    covThresh: float
        coverage threshold for subtyping


    Returns
    -------
    aMatchedStxs: list
        a list of all the stx genes that were matched
    '''
    global bamfile
    bamfile = '%s/%s.unique.sorted.bam' % (bowtie2_output_dir, s_name)

    sCmd = "cat %s | grep '^>' | sed 's/^>//g'" % (ref)
    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    global references
    references = p_stdout.strip().split("\n")
    references.sort()

    """
    references = ['stx1a_gi|147832|gb|L04539.1|',
                  'stx1c_gi|535088|emb|Z36901.1|',
                  'stx1d_gi|28192582|gb|AY170851.1|',
                  'stx2a_gi|14892|emb|X07865.1|',
                  'stx2b_gi|49089|emb|X65949.1|',
                  'stx2c_gi|15718404|dbj|AB071845.1|',
                  'stx2d_gi|30313370|gb|AY095209.1|',
                  'stx2e_gi|8346567|emb|AJ249351.2|',
                  'stx2f_gi|254939478|dbj|AB472687.1|',
                  'stx2g_gi|30909079|gb|AY286000.1|']
                  """

    global aMatchedStxs
    aMatchedStxs = []

    for reference in references:
        #print reference
        mapped_bases = make_pile_up(bamfile, reference)
        matched_stx = parse_pileup(mapped_bases, reference, covThresh)

        if matched_stx is not None:
            aMatchedStxs.append(matched_stx.split('_')[0])

    print 'For coverage threshold', covThresh, ': MatchedStxs', aMatchedStxs

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
        for pileupread in pileupcolumn.pileups:
            pos = pileupcolumn.pos
            if pileupread.qpos:

                base = pileupread.alignment.seq[pileupread.qpos].split()


                if pos in mapped_bases:
                    mapped_bases[pos].append(base[0])
                else:
                    mapped_bases[pos] = base

    return mapped_bases

# --------------------------------------------------------------------------------------------------

def parse_pileup(mapped_bases, reference, covThresh):
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
            #print reference.split('_')[0], each, uniq_mapping, depth

            count_dic = {}
            nts = ['A', 'T', 'C', 'G']
            for nt in nts:
                count = mapped_bases[each].count(nt)
                count_dic[nt] = count

                #if count > 10:
                    #print nt, count_dic[nt], '=',  float(100*count_dic[nt]/depth), '%'

            if depth > 10:
                if uniq_mapping / float(depth) >= covThresh:
                    i += 1
                    #print each, mapped_bases[each]

    #print reference, i, '\n'

    if i == 3:
        #print '%s is present in sample' %reference.split('_')[0], '\n'
        return reference

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


    sCmd = 'samtools view -h -F 4 -bS -o %s/%s.unique.bam %s/%s.sam' \
           % (bowtie2_output_dir, s_name, bowtie2_output_dir, s_name)
    #print "sCmd", sCmd
    #sys.exit()

    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    sCmd ='samtools sort %s/%s.unique.bam -o %s/%s.unique.sorted.bam' \
          % (bowtie2_output_dir, s_name, bowtie2_output_dir, s_name)
    #print "sCmd", sCmd
    #sys.exit()

    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    (p_stdout, p_stderr) = p.communicate()

    sCmd = 'samtools index %s/%s.unique.sorted.bam' % (bowtie2_output_dir, s_name)
    #print "sCmd", sCmd
    #sys.exit()

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

    sCmd = '/home/laura/bowtie2-2.3.1/bowtie2 --fr --minins --no-unal 300 --maxins 1100 -k 99999 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x %s  -1 %s -2 %s -S %s/%s.tmp'\
          % (ref, fastq_read1, fastq_read2, bowtie2_output_dir, s_name)

    p = subprocess.Popen(sCmd, shell=True, stdin=None,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, close_fds=True)

    (p_stdout, p_stderr) = p.communicate()

    tmp = os.path.join(bowtie2_output_dir, s_name +'.tmp') # temporary sam output
    sam = os.path.join(bowtie2_output_dir, s_name +  '.sam')

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
    prefixes = ['@HWI', '@D0']

    for line in lines:
        if line.startswith('@'):
            if not line.startswith(tuple(prefixes)):
            #if not line.startswith('@D0'):
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
