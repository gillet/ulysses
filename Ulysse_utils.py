# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 13:45:41 2012

@author: ingridl & AG
"""
#!/usr/bin/python
__doc__       = """ texte doc """


#import linecache

import os, sys
import csv
import argparse
try:
    import numpy as np
except ImportError:
    print "\tError: NumPy is not installed. \n\
    Please see: http://www.numpy.org/\n\n"
    sys.exit()
import math
import itertools
import time
from bisect import bisect_left

import Ulysse_stats as Ualex

from sys import getsizeof
from itertools import chain
from collections import deque

#-----------------------------------------------------------------------------
def testFile(fil):
        try:
	   with open(fil) as f:
	       return 1
        except IOError as e:
	    return 0
    
#------------------------------------------------------------------------------
def write_VCF(vcffile, dicres, lib, params):
    """ Write a VCF 4.2 output file 

    :param vcffile: filename of output file.vcf
    :type vcffile: string
    :dicres: dictionnary key=sv type, value = tuple (name of inputfile, fdr_for_svtype)
    
    copyright__ = 
    Copyright (C) 2013 - Tim te Beek
    Copyright (C) 2013 - Wai Yi Leung
    Copyright (C) 2013 AllBio (see AUTHORS file)

    __desc__ = 'Convert Ulysse output to pseudo .vcf file format.'
    __created__ = "Mar 18, 2013"
    __author__ = "tbeek (adapted)"    
    
    """
    
    #We will append the vcf file, empty it before:
    if os.path.isfile(vcffile):
        os.remove(vcffile)
    
    
    Header = True
    for sv, tupsv in dicres.iteritems():
        #print "TUPSV", tupsv
        tsvfile = tupsv[0]
        fdr_threshlod = tupsv[1]        
        with open(tsvfile) as reader:
            # Parse file
            dictreader = _parse_tsvfile(reader)
            #print dictreader.fieldnames
    
            # Write out file
            _format_vcffile(dictreader, vcffile, Header, sv, fdr_threshlod, lib, params)
    
        # Quick output
#        with open(vcffile) as reader:
#            print reader.read(1000)
        Header = False



def _parse_tsvfile(readable):
    '''
    Read readable using csv.Sniffer and csv.DictReader
    :param readable: open file.tsv handle to read with csv.DictReader
    :type readable: file
    '''
    prev, curr = 0, 0
    while True:
        line = readable.readline()
        if not line.startswith('#'):
            # lets start from prev # line, without the hash sign
            readable.seek(prev + 1)
            break
        else:
            prev = curr
            curr = readable.tell()

    # Determine dialect
    curr = readable.tell()
    dialect = csv.Sniffer().sniff(readable.read(3000))
    #dialect = 'excel-tab'
    readable.seek(curr-1)

    # Read file
    dictreader = csv.DictReader(readable, dialect=dialect)
    return dictreader


_tsv_fields = ('Library', '(pair of) chromosome(s)', 
               'ID', 'nbRP', 'nbA', 'nbB', 'left_borderA', 
               'right_borderA', 'deltaA', 'cen_posA', 'left_borderB', 
               'right_borderB', 'deltaB', 'cen_posB', 'SV_size_min', 
               'SV_size_max', 'p-value', 'cov', 'AvrQual', 'pval')
               
#('Chr1', 'Pos1', 'Orientation1',
#               'Chr2', 'Pos2', 'Orientation2',
#               'Type', 'Size', 'Score',
#               'num_Reads', 'num_Reads_lib',
#               'ERR031544.sort.bam')

# 'Chr1': '1',
# 'Pos1': '269907',
# 'Orientation1': '39+39-',
# 'Chr2': '1',
# 'Pos2': '270919',
# 'Orientation2': '39+39-',
# 'Type': 'DEL',
# 'Size': '99',
# 'Score': '99',
# 'num_Reads': '38',
# 'num_Reads_lib': '/home/allbio/ERR031544.sort.bam|38',
# 'ERR031544.sort.bam': 'NA'

_vcf_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO\n')


def _format_vcffile(dictreader, vcffile, Header, svtype, fdr_threshlod, lib, params):
    '''
    Create a pseudo .vcf file based on values read from DictReader instance.
    :param dictreader: DictReader instance to read data from
    :type dictreader: csv.DictRedaer
    :param vcffile: output file.vcf filename
    :type vcffile: string
    '''
    stats, chrDicos = read_stats(params["stats"])
    median = str(stats["median"])
    #print "toto"
    with open(vcffile, mode='a') as writer:
        if Header:
            date = time.strftime("%d/%m/%Y")
            writer.write("""##fileformat=VCFv4.1\n\
##fileDate="""+date+"""\n\
##ALT=<ID=DEL,Description="Deletion">\n\
##ALT=<ID=DUP,Description="Tandem Duplication">\n\
##ALT=<ID=INV,Description="Inversion">\n\
##ALT=<ID=INS,Description="Insertion">\n\
##ALT=<ID=sINS,Description="Small Insertion">\n\
##ALT=<ID=RT,Description="Reciprocal Translocation">\n\
##ALT=<ID=NRT,Description="Non Reciprocal Translocation">\n\
##FILTER=<ID=fdr,Description="SV below FDR">\n\
##FILTER=<ID=q"""+str(params['mapq'])+""",Description="RP average mapping quality below """+str(params['mapq'])+"""">\n\
##FILTER=<ID=q"""+str(params['mapqx'])+""",Description="RP average mapping quality below """+str(params['mapqx'])+"""">\n\
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="SOMATIC, GERMLINEor UNKNOWN status">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of Structural variant">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n\
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">\n\
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\n\
##INFO=<ID=LIB,Number=1,Type=String,Description="Name of the input file">\n\
##INFO=<ID=CHR,Number=1,Type=String,Description="(pair of) chromosome(s)">\n\
##INFO=<ID=ID,Number=1,Type=String,Description="Internal ID">\n\
##INFO=<ID=nbRP,Number=1,Type=Integer,Description="Number of RP describing the SV">\n\
##INFO=<ID=nbA,Number=1,Type=Integer,Description="nb of RP describing the SV on chromosome A">\n\
##INFO=<ID=nbB,Number=1,Type=Integer,Description="nb of RP describing the SV on chromosome B">\n\
##INFO=<ID=LBA,Number=1,Type=Integer,Description="5' border of the SV on chromosome A">\n\
##INFO=<ID=RBA,Number=1,Type=Integer,Description="3' border of the SV on chromosome A">\n\
##INFO=<ID=DA,Number=1,Type=Integer,Description="SV range on chromosome A between right and left borders">\n\
##INFO=<ID=CPA,Number=1,Type=String,Description="Four digit code giving the position of reads with respect to centromere (L for left arm and R for right arm) and strand orientation (+/-). The first two digits for read1 and the last two digits for read2">\n\
##INFO=<ID=LBB,Number=1,Type=Integer,Description="5' border of the SV on chromosome B">\n\
##INFO=<ID=RBB,Number=1,Type=Integer,Description="3' border of the SV on chromosome B">\n\
##INFO=<ID=DB,Number=1,Type=Integer,Description="SV range on chromosome B between right and left borders">\n\
##INFO=<ID=CPB,Number=1,Type=String,Description="Four digit code giving the position of reads with respect to centromere (L for left arm and R for right arm) and strand orientation (+/-). The first two digits for read1 and the last two digits for read2">\n\
##INFO=<ID=SVMI,Number=1,Type=Integer,Description="Minimum estimated SV size">\n\
##INFO=<ID=SVM,Number=1,Type=Integer,Description="Maximum estimated SV size">\n\
##INFO=<ID=PBAL,Number=1,Type=String,Description="For SV with two junctions only. P-value of the binomial test (see methods)">\n\
##INFO=<ID=COV,Number=1,Type=String,Description="Coverage +-10kb">\n\
##INFO=<ID=AVRQUAL,Number=1,Type=Float,Description="Average read mapping quality">\n\
##INFO=<ID=PVAL,Number=1,Type=Float,Description="P-value of the SV statistical validation test">\n\
""")
##FILTER=<ID=LowQual,Description="INV with only 1PS in one of the 2 clusters">\n\


            writer.write('#{}'.format('\t'.join(_vcf_fields)))
        info_fields = {'LIB':'Library', 'CHR':'(pair of) chromosome(s)', 'ID':'ID', 
                       'nbPS':'nbRP', 'nbA':'nbA', 'nbB':'nbB', 'LBA':'left_borderA', 
                       'RBA':'right_borderA', 'DA':'deltaA', 'CPA':'cen_posA', 
                       'LBB':'left_borderB', 'RBB':'right_borderB', 'DB':'deltaB',
                       'CPB':'cen_posB', 'SVMI':'SV_size_min', 'SVM':'SV_size_max', 
                       'PVAL':'pval', 'AVRQUAL':'AvrQual', 'COV':'cov', 'PBAL':'p-balanced'}
        inv_info_fields = {v:k for k, v in info_fields.items()}
        #print "inv_info_fields", inv_info_fields
                        #'PBAL':'p-balanced', 'COV':'cov',
        order_fields = ['Library', '(pair of) chromosome(s)', 'ID',
               'nbRP', 'nbA', 'nbB', 'left_borderA', 
               'right_borderA', 'deltaA', 'cen_posA', 'left_borderB', 
               'right_borderB', 'deltaB', 'cen_posB', 'SV_size_min', 
               'SV_size_max', 'p-balanced', 'cov', 'AvrQual', 'pval']

        alreadyIn = []
        output_vcf = []
        for line in dictreader:
            #print "line", line
            #CHROM = line['Chr1']
            CHROM = line[_tsv_fields[1]].split("-")[0]
            POS = int(line[_tsv_fields[6]])
            ID = line[_tsv_fields[1]]+"."+line[_tsv_fields[2]]
            #ID='.'
            REF = "N"
            if svtype in ['NRT', 'sINS', 'RT']:
                svtypeVCF = 'INS'
            else:
                svtypeVCF = svtype
            ALT = "<"+svtypeVCF+">"

            try:
                QUAL = int(abs(-10*math.log(float(line[_tsv_fields[19]]))))
            except ValueError:
                QUAL = 0
 
            #print "LINE", line
            if float(line[_tsv_fields[19]]) > fdr_threshlod:
                FILTER = 'fdr'
            #elif (int(line['nbA'])<2 or int(line['nbB'])<2) and svtype =='INV':
            #    FILTER = 'LowQual'
            elif float(line[_tsv_fields[18]]) < int(params['mapqx']) and svtype in ['NRT', 'sINS', 'RT', 'INS']:
                FILTER = 'q'+str(params['mapqx'])
            elif float(line[_tsv_fields[18]]) < int(params['mapq']) and svtype in ['DUP', 'DEL', 'INV']:
                FILTER = 'q'+str(params['mapq'])
            else:
                FILTER = 'PASS'

  
            INFOListe = []
            INFOListe.append("UNKNOWN")
            if svtype in ['NRT', 'sINS', 'RT']:
                svtypeCorrec = 'INS'
            else:
                svtypeCorrec = svtype
            INFOListe.append("SVTYPE="+svtypeCorrec)
            #print line['left_borderB'], line['left_borderA']
            INFOListe.append("END="+str(max(int(line['left_borderB']), int(line['left_borderA']))))
            INFOListe.append("CIEND=-"+median+","+median)            
            INFOListe.append("CIPOS=-"+median+","+median)
            for key in order_fields:
                INFOListe.append(str(inv_info_fields[key])+"="+str(line[key]))
            #print "INFOListe", INFOListe
            INFOListe[-1] = "PVAL="+str("{:.5e}".format(float(line['pval'])))
            #print line
            INFOListe[7] = '.'.join(['ID='+line['(pair of) chromosome(s)'], line['ID']])

            INFO = ';'.join(INFOListe)
            nameSV = CHROM + "-" + line['left_borderA'] + "-" + line['right_borderA'] + "-" + FILTER
            
            # Create record
            
            if nameSV not in alreadyIn:    
                output_vcf.append([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO])
                alreadyIn.append(nameSV)

        # Sort all results
        output_vcf.sort()
        output = "\n".join(["\t".join(map(str,vcf_row)) for vcf_row in output_vcf])
        # Write record
        writer.write(output)
        if output:
            writer.write("\n")        


#-------------------------------------------------------------------------
def allsame(seq):
    """determine if all members of a list are the same"""
    return min([elem==seq[0] for elem in seq]+[True])
#-------------------------------------------------------------------------
def getcommonstart(seqlist):
    """retrieve the common start of the members of a list"""
    if len(seqlist)>1:
        letters = itertools.izip(*seqlist)                 # let's be lazy
        common = itertools.takewhile(allsame, letters)     # still lazy
        return ''.join([letters[0] for letters in common]) # merge back
        
    else:
        return ''
    
#-------------------------------------------------------------------------
def parser(sv_type):

    message = "\n######## ULYSSES V1.0 #############\n\n Detection of "+sv_type
    parser = argparse.ArgumentParser(description=message)
    parser.add_argument("-p", metavar='parameter_file_name',
                        help='if not given, set to ulysses_params',
                        default='ulysses_params')

    args = parser.parse_args()

    print "Parameter file :", args.p,"\n\n"

    return os.path.abspath(args.p)

#----------------------------------------------------------------------------
def getori(nb):
	ori = "-"
	if not nb :
		ori = "+"
	return ori

#--------------------------------------------------------------------------
def Sortread(name, ori1, chr1, pos1, ori2, chr2, pos2, length):
    """ Sort reads in PS from BAM file. Reads ordered according to the
    chromosome for inter-chromosomal PS. Reads ordered according to position
    for intra-chromosome PS."""

    if chr1 == chr2 :
        if int(pos2) < int(pos1):
            tempo= (name, ori2, chr2, pos2, ori1, chr1, pos1, length)
        else:
            tempo= (name, ori1, chr1, pos1, ori2, chr2, pos2, length)
    else:
        if chr1 > chr2 :
            tempo=(name,ori2,chr2,pos2,ori1,chr1,pos1,0)
        else:
            tempo=(name,ori1,chr1,pos1,ori2,chr2,pos2,0)

    return tempo

#--------------------------------------------------------------------------
def getCoord(read1, chr1, chr2):
    """ Get the mate coordinates with reads in correct orientation"""

    pos1=read1.pos + 1 #0-based leftmost coordinate
    pos2=read1.pnext + 1
    ori1 = getori(read1.is_reverse)
    ori2 = getori(read1.mate_is_reverse)

    tempo = Sortread("toto",ori1, chr1, pos1, ori2, chr2, pos2, 0 )
    return tempo[1], tempo[2], tempo[3], tempo[4], tempo[5], tempo[6]

#--------------------------------------------------------------------------
def ps_fictive(PS):
    """ Returns True if PS name has all elements of A PS fictive name"""

    Fic = "mock.TN.NApb.chr22-chr21.2PS-0PS.PS"
    Fic = "1209:8995:58751"
    if isinstance(PS[0], str):
        if Fic in PS[0]:
            return True
    else:
        return False
#    return False


#--------------------------------------------------------------------------

def getPSSetListe(dicos):
    """After dectection, and before Affichage, can take the SV candidats
    dict and report the set of PS"""
    PS = []
    for key, val in dicos.iteritems():
        for sv in val:
        #print val
            PS.extend(sv[0])
            PS.extend(sv[1])
    #return set(PS)
    return PS


#--------------------------------------------------------------------------
def tuples2list(vals):
    """Transform complex lists to complex tuples"""
    l=[]
    for val in vals:
        if type(val) == tuple:
            l.append(list(val))
        else:
            l.append(val)
    return l


#--------------------------------------------------------------------------
def list2tuples(vals):
    """Transform complex lists to complex tuples"""
    l=[]
    for val in vals:
        if type(val) == list:
            l.append(tuple(val))
        else:
            l.append(val)
    return tuple(l)


#--------------------------------------------------------------------------
def get_range(liste):
    """extract list of chromosomes from liste retrieved in run info file"""

    list_x = []
    liste = liste.replace(" ","")
    elem = liste.split(",")
    for xsome in elem:
        if "-" in xsome:
            list_x += range(int(xsome.split("-")[0]), int(xsome.split("-")[1])+1)
        else:
            list_x += [xsome]
#    print "chromosome range ", list_x
    return list_x
#--------------------------------------------------------------------------
def get_run_info(fichier):
    """retrieve parameters for detection in the parameter file given in
    argument"""
    pars = {}
    for line in open(fichier,"r"):
        if "=" in line:
            pars[line.split("=")[0].lower()] = line[:-1].split("=")[1]

    checkf = {
    "in":"name of input file. example human for human_X X name of chromosome",\
    "mapq":"minimum mapping quality of reads to select a PS",\
    "out":"name of output file",\
    "stats":"name of statistics file indicating mean=X, median=Y and MAD=Z",\
    "only_stats":"If true, only statistical tests are performed",\
    "NSv":"Maximum number of expected SV candidate",\
    "annotation":"name of file containing coordinates of genetic elements \
    (centromere and telomere are necessary)",\
    "field_chr":"column in \"annotation\" for chromosome number",\
    "field_start":"column in \"annotation\" for element start",\
    "field_end":"column in \"annotation\" for element end",\
    "field_type":"column in \"annotation\" for feature type",\
    "field_sep":"field separator",\
    "range":"enumerate chromosomes example : 1-22,X,Y",\
    "n":"factor for definition of d and ics parameters",\
    "fdr":"limit for deletion fdr", "mapq":"minimum mapping quality for \
     the 2 reads of a given PS to be considered"}
    
    for cle, valeur in pars.items():
        if cle == "range" and valeur != "all":
            pars[cle] = get_range(valeur)
        elif cle in ['mapq', 'mapqx', 'nsv', 'field_type', 'field_start', 'field_end']:
            
            pars[cle] = int(valeur)

    
    pars["fdr"] = float(pars["fdr"])

    if len(pars) < len(checkf):
        print "Error: missing infos in ", fichier
        print "Required :\n"
        for cle in checkf:
            if cle not in pars.keys():
                print cle, checkf[cle]
        return 0

    return pars

#--------------------------------------------------------------------------
def create_files(sv_type, fichier, deletion = False):
    """creates header of output files"""
    #deletion, duplication, inversions
    if deletion:
        out = open(fichier + "_" + sv_type + "_bySV.csv", "w")
        out.write("\"Library\";\"(pair of) chromosome(s)\";\"ID\";\"nbRP\";\"nbA\";\
\"nbB\";\"left_borderA\";\"right_borderA\";\"deltaA\";\"cen_posA\";\
\"left_borderB\";\"right_borderB\";\"deltaB\";\"cen_posB\";\"SV_size_min\";\
\"SV_size_max\";\"p-balanced\";\"\n")
        out.close()
        out = open(fichier +"_" +sv_type + "_byRP.csv", "w")
        out.write("\"Library\";\"(pair of) chromosome(s)\";\"ID\";\"nbRP\";\"nbA\";\
\"nbB\";\"left_borderA\";\"right_borderA\";\"deltaA\";\"cen_posA\";\
\"left_borderB\";\"right_borderB\";\"deltaB\";\"cen_posB\";\"SV_size_min\";\
\"SV_size_max\";\"p-balanced\";\"RP\";\"str1\";\"chr1\";\"pos1\";\"st2r\";\
\"chr2\";\"pos2\";\"L\"\n")
        out.close()
    else:
        out = open(fichier + "_" + sv_type + "_bySV.csv", "w")
        out.write("\"Library\";\"(pair of) chromosome(s)\";\"ID\";\"nbRP\";\"nbA\";\
\"nbB\";\"left_borderA\";\"right_borderA\";\"deltaA\";\"cen_posA\";\
\"left_borderB\";\"right_borderB\";\"deltaB\";\"cen_posB\";\"SV_size_min\";\
\"SV_size_max\";\"p-balanced\";\"cov\"\n")
        out.close()
        out = open(fichier +"_" +sv_type + "_byRP.csv", "w")
        out.write("\"Library\";\"(pair of) chromosome(s)\";\"ID\";\"nbRP\";\"nbA\";\
\"nbB\";\"left_borderA\";\"right_borderA\";\"deltaA\";\"cen_posA\";\
\"left_borderB\";\"right_borderB\";\"deltaB\";\"cen_posB\";\"SV_size_min\";\
\"SV_size_max\";\"p-balanced\";\"cov\";\"RP\";\"str1\";\"chr1\";\"pos1\";\"str2\";\
\"chr2\";\"pos2\";\"L\"\n")
        out.close()

#--------------------------------------------------------------------------
def read_stats(statfile):
    """ Open statfile containing mean, median and MAD for the data
    """
    chrDicos = {}
    stats = {}
    reels = ["mean", "median", "mad", "stdev"]
    entiers = ["ndup", "ninv", "nsins", "ndel", "ninter", "genome_length",
               "nsdel","rl"]
    with open(statfile,"r") as statfile:
        for l in statfile :
            elem = l.lower().split()[0]
            if elem.startswith("ps_type"):
                 stats[elem] = l.split()[2]
            if elem == "read_length":
                 stats["rl"] = l.split()[2]
            if elem.startswith("chromosome_prefix"):
                 try:
                     stats[elem] = l.split()[2]
                 except:
                     stats[elem] = ""
            if elem in reels :
                 stats[elem] = float(l.split()[2])
            elif elem in entiers:
                stats[elem] = int(l.split()[2])
            elif l.startswith("length"):
                chrx = l.split()[1]
                chrDicos[chrx] = int(l.split()[2])


    return stats, chrDicos

#--------------------------------------------------------------------------
def cleanRange(c):
    csplit = c.split(",")
    l = []
    for elt in csplit:
        css = elt.split("-")
        if len(css) == 1:
            if type(css[0])  == int:
                l.append(int(css[0]))
            else:
                print "\n\n\t************** Error : Range of chromosomes is not well defined\n"
                sys.exit()
        elif len(css) == 2:
            try:
                ll = range(int(css[0]), int(css[1])+1)
                l.extend(ll)
            except:
                print "\n\n\t************** Error : Range of chromosomes is not well defined\n"
                sys.exit()
        else:
            print "\n\n\t************** Error : Range of chromosomes is not well defined\n"
            sys.exit()
    return l
                
                
#--------------------------------------------------------------------------
def prepare_detection(params):
    """ Read stat files, read run info
    """

    stats, chrDicos = read_stats(params["stats"])

    if params["range"] == "all":
        liste = chrDicos.keys()
        liste.sort()
        params["range"] = liste
    else:
        l = cleanRange(params["range"])
        liste = []
        for cx in l:
            print "---------------------------------- ", stats["chromosome_prefix"]+str(cx)
            liste.append(stats["chromosome_prefix"]+str(cx))
        params["range"] = liste

    return stats, chrDicos

#--------------------------------------------------------------------------
def get_subtelo_limits(params):
    """ Retrieve the coordinates of subtelomeres from annotation file"""
    subtelo=[]
    if os.path.isfile(params["annotation"]):
        with open(params["annotation"], "r") as elems:
            temp = [pi[:-1] for pi in elems if "subtelomere" in pi]
        for s in temp:
            if params["field_sep"] == ";":
                s = s.split(";")
            else:
                s = s.split()
            if s[params["field_type"]] == "subtelomere":
                bras = "R"
                debut = int(params["field_start"]) - 1
                fin = int(params["field_end"]) - 1
                start = int(s[debut])
                end = int(s[fin])
                coord = start
                #distinguish between left and right subtelomeres
                #arbitrary, feets for yeasts and humans
                if end < 70000 :
                    bras = "L"
                    coord = end
                subtelo.append([s[int(params["field_chr"])-1], bras, coord])
    return subtelo

#--------------------------------------------------------------------------
def get_centro_coords(params):
    """ Retrieve the coordinates of centromeres from annotation file"""
    centrom=[]

    if params["annotation"] != 'None': #eval(params["annotation"]):
        
        with open(params["annotation"], "r") as elems:
            
            temp = [pi[:-1] for pi in elems if "centromere" in pi.lower()]
        for pi in temp:
            ii = pi.split()

            
#            if params["field_sep"] == ";":
#                ii = pi.split(";")
            if params["field_sep"] not in ['Tb', '\t', ' ', 'tb']:
                ii = pi.split(params["field_sep"])
                
                
            if len(ii)<2:
                sys.exit("\n *** Error: Annotation file is note correctly formated or separator not specified ***\n")
            
            if ii[int(params["field_type"])-1] == "centromere":
                centrom.append([ii[int(params["field_chr"])-1], int(ii[int(params["field_start"])-1]),
                            int(ii[int(params["field_end"])-1])]) #chr, start, end
                            
    
        print"Centromeres:" #, centrom
        Ualex.printplus(centrom)
    else:
        print "\nCentromeres positions not given. Distinction between INS and RT will not be optimal\n"
    return centrom

#--------------------------------------------------------------------------
def Subtelo(subtelo,  chrx,  coord):
    """ retourne 1 si la coordonnee est dans une region subtelomerique,
selon fichier subtelo_limit (dernier gene essentiel)"""
    sub = 0
    for i_s in subtelo:
        if chrx == i_s[0]:
            if (i_s[1] == "L") and (coord < i_s[2]):
                sub = 1
            if (i_s[1] == "R") and (coord > i_s[2]):
                sub = 1
        if sub == 1:
            break
    return sub

#--------------------------------------------------------------------------
def Centro_Inter(centrom,  chr1,  s1,  coord1,  chr2,  s2,  coord2):
    """If centromeres coordinates are known, returns the arm (L/R) and
    orientation (+/-) of the two reads of a PS.
    If centromeres coordinates are not known, returns N for arm and i
    for orientation"""

    cen = 0
    combi = ["Ni", "Ni"]
    i = 0

    while i < len(centrom) and cen < 2:
        start = centrom[i][1]
        end = centrom[i][2]
        if chr1 == centrom[i][0]:
            if (coord1 < start):
                combi[0] = "L"+s1
            elif (coord1 > end ):
                combi[0] = "R"+s1
            cen += 1
        if chr2 == centrom[i][0]:
            if (coord2 < start ):
                combi[1] = "L"+s2
            elif (coord2 > end ):
                combi[1] = "R"+s2
            cen += 1
        i = i+1

    return combi[0] + combi[1]

#--------------------------------------------------------------------------
def sortListsInDicos(dicos, pos):
    """ Sort list of lists in a dico """
    for cle, val in dicos.iteritems():
        l = dicos[cle]
        l.sort(key=lambda x: x[pos])
        dicos[cle] = l
    return dicos
#--------------------------------------------------------------------------
def parameter_ok(val1, val2, param):
    """Retourne True if difference between val1 and val2 < param"""
    return (abs(val1 - val2) < param)

#--------------------------------------------------------------------------
def binary_search(a, x, lo=0, hi=None):   # can't use a to specify default for hi
    """Very fast search of x in sorted array
    Return index if found, else returns -1"""
    hi = hi if hi is not None else len(a) # hi defaults to len(a)
    pos = bisect_left(a,x,lo,hi)          # find insertion position
    return (pos if pos != hi and a[pos] == x else -1) # don't walk off the end

#--------------------------------------------------------------------------
def overlap(xsome, candidats, icand, signe, d, opposite):
    """Definit si le PS cand a des copains dans une autre orientation
        range tous dans overreads (pour savoir si subtelo_transloc)
        et uniquement ceux dans orientation opposee dans overopps.
        Donc c'est 2d"""

    cand = xsome[icand]
    coord1 = cand[3]
    coord2 = cand[6]
    overreads = []
    overopps = []

    signeOpp = opposite[signe]

    oris = ["++", "--", "+-", "-+"]
    oris.remove(signe)

    for ori in [signeOpp]:
    #for ori in oris:
        for i in range(len(candidats[ori])):
            li = xsome[candidats[ori][i]]

#            if ps_fictive(cand) and ps_fictive(li):
#                print "overlap ulysse 0 ", d, cand, li, parameter_ok(coord1, li[3], 2*d), parameter_ok(coord2, li[6], 2*d)



            #if parameter_ok(coord1, li[3], d) and  parameter_ok(coord2, li[6], d):
            if parameter_ok(coord1, li[3], 2*d) or parameter_ok(coord2, li[6], 2*d):
                overreads.append(candidats[ori][i])
                #if cand[0] == "2:70:4864:8181#0":
                #    print "overlap ", cand, li
                if ori == opposite[signe]:
                    overopps.append(candidats[ori][i])

    #determine les orientations dans lesquelles il y a des

#    if ps_fictive(cand):
#        print "Ulysse overlap ", cand[0], overreads, overopps
    return overreads, overopps

#--------------------------------------------------------------------------
def overlap_intra(xsome, candidats, icand, signe, d, opposite):
    """Definit si le PS cand a des copains dans une autre orientation
        range tous dans overreads (pour savoir si subtelo_transloc)
        et uniquement ceux dans orientation opposee dans overopps."""

    cand = xsome[icand]
    coord1 = cand[3]
    coord2 = cand[6]
    overreads = []
    overopps = []

    oris = ["++", "--"]
    oris.remove(signe)
    for ori in oris:
        for i in range(len(candidats[ori])):
            li = xsome[candidats[ori][i]]
            if (abs(coord1-li[3]) < d) and  (abs(coord2-li[6]) < d):
                overreads.append(candidats[ori][i])
                #if cand[0] == "2:70:4864:8181#0":
                #    print "overlap ", cand, li
                if ori == opposite[signe]:
                    overopps.append(candidats[ori][i])

    #determine les orientations dans lesquelles il y a des

    return overreads, overopps


#--------------------------------------------------------------------------
def TestArmOri(xsome, centrom, allowed, ps1, ps2):
    """ If centromeres coodinates are known, returns list of PS compatibles
    with Translocation (compat) and list of PS compatibles with Insertion.
    If centromeres coordinates are unknown, compat will contain all PS because
    in this case, no distinction can be made"""
    comp = []
    notcomp = []
    bk1 = {}
    for ps in ps1:
#        print "$$$$ ",ps, ps1, type(ps)
        tr = xsome[ps]

        b = Centro_Inter(centrom, tr[2], tr[1], tr[3], tr[5], tr[4], tr[6])
        if b in bk1:
            bk1[b].append(ps)
        else:
            bk1[b] = [ps]

    bk2 = {}
    for ps in ps2:
#        print "$$$$ ",ps, ps2, type(ps)
        tr = xsome[ps]
        b = Centro_Inter(centrom, tr[2], tr[1], tr[3], tr[5], tr[4], tr[6])

        if b in bk2:
            bk2[b].append(ps)
        else:
            bk2[b] = [ps]


    for b1 in bk1:
        ok = False
	# POURQUOI AND NOT OK
        #while b < len(bk2cle) and not ok:
	# Si on laisse not ok on n'a pas tous les PS comp
        for b2 in bk2:

            for i in allowed:
                if (b1, b2) == i or (b2, b1) == i:
                    #comp.append([bk1[b1], bk2[b2], b1, b2])
                    ok = True
#        if not ok:
#            notcomp.append([bk1[b1], bk2[b2], b1, b2])
        if ok:
            comp=[bk1[b1], bk2[b2], b1, b2]
        else:
            notcomp=[bk1[b1], bk2[b2], b1, b2]
#        if Fictive:
#            print "TestArmOri ",tr, comp, notcomp, len(comp),len(notcomp)

    return comp, notcomp

#--------------------------------------------------------------------------
def Limites(xsome, s1, s2):

    """ Definit les extremites des SV"""

    deb1, fin1, deb2, fin2 = Extremites(xsome, s1)
    d1, f1, d2, f2 = Extremites(xsome, s2)

    if deb1 == -1 and d1 > 0:
        return d1, f1, d2, f2

    elif d1 == -1 and deb1 > 0:
        return deb1, fin1, deb2, fin2

    elif d1 > 0 and deb1 > 0:
        deb1 = min(d1, deb1)
        fin1 = max(fin1, f1)
        deb2 = min(d2, deb2)
        fin2 = max(fin2, f2)
        return deb1, fin1, deb2, fin2
    else:
        return -1, -1, -1, -1


#------------------------------------------------------------------------------
def ExtensionSV_old(candid, select, d, ids):
    """ For Duplication and Deletion SV
    Search for PS not in the first homogeneous group at right end"""
    #recupere tous les ids du groupe select
    pop = [select[i][7] for i in candid]
    #recherche le min et max ids du groupe
    idsv = (min(pop), max(pop))
    #recupere tous les coord left du groupe select
    pop = [select[i][3] for i in candid]
    #recherche le min des lpop
    leftv = (min(pop), max(pop))
    #recupere tous les coord right du groupe select
    pop = [select[i][6] for i in candid]
    #recherche le min et le max des rpop
    rightv = (min(pop), max(pop))

    temp = []
#    j = candid[-1] + 1
    k = candid[-1] + 1
    #reads compatibles a droite
    #VIRER CONDITION CHEVAUCHE STRICT
#    while j < len(select) and (select[j][3] < rightv[0]) :
#    while j < len(select) :
    ok = True
    while ok:
        if select[k:] == []:
            ok = False
        for j, val in enumerate(select[k:], start = k):
            #condition left coordinates
            a_gauche = parameter_ok(val[3], leftv[0], d)
            #condition rigth coordinates
    #        a_droite = parameter_ok(select[j][6], rightv[1], d)
            a_droite = parameter_ok(val[6], rightv[1], d)
            #condition insert size
            idsok = parameter_ok(val[7], idsv[1], ids)
            #OR dans v9 mais je ne comprends plus pourquoi ce n'est pas AND
            #ca ne change rien dans la detection des deletions dans manip LJD
            idsok = idsok or parameter_ok(val[7], idsv[0],  ids)
            #idsok = idsok and parameter_ok(select[j][7], idsv[0],  ids)
    
            if idsok and a_gauche and a_droite :
                temp.append(j)
            j = j + 1
            if abs(val[3] - leftv[0])> 2*d or j>=len(select):
                ok = False

    return temp

#------------------------------------------------------------------------------
def ExtensionSV(candid, select, d, ids, cladic):
    """ For Duplication and Deletion SV
    Search for PS not in the first homogeneous group at right end"""
    #print "in Extension SV dup"    
    #recupere tous les ids du groupe select
    pop = [cladic[i][7] for i in candid]
    #recherche le min et max ids du groupe
    idsv = (min(pop), max(pop))
    #recupere tous les coord left du groupe select
    pop = [cladic[i][3] for i in candid]
    #recherche le min des lpop
    leftv = (min(pop), max(pop))
    #recupere tous les coord right du groupe select
    pop = [cladic[i][6] for i in candid]
    #recherche le min et le max des rpop
    rightv = (min(pop), max(pop))

    temp = []
#    j = candid[-1] + 1
    k = candid[-1] + 1
    #print "--->", len(select), k
    #reads compatibles a droite
    #VIRER CONDITION CHEVAUCHE STRICT
#    while j < len(select) and (select[j][3] < rightv[0]) :
#    while j < len(select) :
    ok = True
    while ok:
        if select[k:] == []:
            ok = False
        for j, val in enumerate(select[k:], start = k):
            #print "***--->", j, val
            #condition left coordinates
            a_gauche = parameter_ok(val[3], leftv[0], d)
            #condition rigth coordinates
    #        a_droite = parameter_ok(select[j][6], rightv[1], d)
            a_droite = parameter_ok(val[6], rightv[1], d)
            #condition insert size
            idsok = parameter_ok(val[7], idsv[1], ids)
            #OR dans v9 mais je ne comprends plus pourquoi ce n'est pas AND
            #ca ne change rien dans la detection des deletions dans manip LJD
            idsok = idsok or parameter_ok(val[7], idsv[0],  ids)
            #idsok = idsok and parameter_ok(select[j][7], idsv[0],  ids)
    
            if idsok and a_gauche and a_droite :
                temp.append(j)
            j = j + 1
            #print "    --->", len(select), k, j
            if abs(val[3] - leftv[0])> 2*d or j>=len(select):
                ok = False
    #print "exiting Extension SV dup"
    return temp
#------------------------------------------------------------------------------
def conv(g):
    """split concatenarted liste to a liste"""
    return [int(x) for x in g.split(".")]

#------------------------------------------------------------------------------
def inde_groups(sel, tri):
    """create groups with indep PS????"""
    tri.sort()
    if tri and tri not in sel:
        sel.append(tri)
        
#------------------------------------------------------------------------------
def is_number(s):
   try:
       float(s)
       return True
   except:
       return False       
#------------------------------------------------------------------------------
def inde_groupsMinPS(sel, tri, minPS):
    """create groups with indep PS????"""
    tri.sort()
    if tri and tri not in sel and len(tri)>=minPS:
        sel.append(tri)
        
##--------------------------------------------------------------------------
#def testStat(n1, n2):
#    """Fais un test binomial pour repartition des PS dans les deux orientations
#        et retourne la p-value. Si p-value < 0.05 H0 rejetee: la distrib
#        n'est pas homogène"""
#    r=pyper.R()
#    pvalue = -1
#    if n1+n2 > 0:
#        #Faire un test binomial x = PS++ ou PS-- et n = PS tot p= 0.5
#        x = str(max(n1, n2))
#        n = str(n1+n2)
#        r("res=binom.test("+x+", "+n+", p=0.5)[]")
#        p=r('res$p.value')
#        pvalue=float(p.split()[2])
#    return pvalue
#--------------------------------------------------------------------------
def testStatSciPy(n1, n2):
    """Fais un test binomial pour repartition des PS dans les deux orientations
        et retourne la p-value. Si p-value < 0.05 H0 rejetee: la distrib
        n'est pas homogène"""
    x = max(n1, n2)
    #n = n1+n2
    y = min(n1, n2)
    #return stats.binom_test(x,n)
    #print "min/max", y, x
    return y/(x+0.0)
 




def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements} """
                    
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        #if verbose:
        #    print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)



  
#--------------------------------------------------------------------------
  
def getEltPosition(myListe, elt):
    """return the position of an element 'elt' in 'myListe' """
    #print "elt", elt
    #print "myliste", myListe
    for i,element in enumerate(myListe):
        if element == elt:
            #print "position", i
            return i
    return -1



   
#--------------------------------------------------------------------------

def filter_list(L):
    """ Remove sublists !! Carefull, is slow with large lists """
    return [x for x in L if not any(set(x)<=set(y) for y in L if x is not y)]
    

#--------------------------------------------------------------------------

def filter_listOfLists(L):
    """ Remove sublists that contain lists!! Carefull, is slow with large lists """    
    return [x for x in L if not any(set(map(tuple,x))<=set(map(tuple,y)) for y in L if x is not y)]  
 
#--------------------------------------------------------------------------

def addMeanSVQuality(bySVfileName, byPSfileName, dicQual):
    """ function that can add mean PS SV quality in any SV detection file pair (SV and PS) """
    #
    ###print "\nStart removing duplicates in files", bySVfileName, "and", byPSfileName
    #open byPS file
    bPStmpMem = csv.reader(open(byPSfileName), delimiter=';', quotechar='"')
    bPStmp = []
    for row in bPStmpMem:
        bPStmp.append(row)
    bPSraw = bPStmp[1:] #remove header
    bPSheader = bPStmp[0]
    #
    #
    #open bySV file
    bSVtmpMem = csv.reader(open(bySVfileName), delimiter=';', quotechar='"')
    bSVtmp = []
    for row in bSVtmpMem:
        bSVtmp.append(row)
    bSVraw = bSVtmp[1:] #remove header
    bSVheader = bSVtmp[0]
    #
    #creast short names
    bSV = bSVraw
    bPS = bPSraw
    #
    # clean headers
    bSVheader = [x.replace("    \"","") for x in bSVheader]
    bPSheader = [x.replace("    \"","") for x in bPSheader]
    bSVheader = [x.replace("\"","") for x in bSVheader]
    bPSheader = [x.replace("\"","") for x in bPSheader]
    #
    #
    #make liste of involved chrs
    chrList = list(set([x[1] for x in bSV]))
    #
    #
    SVfull = {} #dicos for by PS
    SVsimple = {} #dicos for by SV
    for cChr in chrList: #initialize the dicos
        SVfull[cChr] = []
        SVsimple[cChr] = []
    #
    #split PS by chr
    for PS in bPS:
       SVfull[PS[1]].append(PS) #fill each key (chrom) with PS in the dic
    #
    #split SV by chr
    for iSV in bSV:
        SVsimple[iSV[1]].append(iSV)
    #
    #split PS by SV ID
    for chromo,PSs in SVfull.iteritems():
        split = {}
        for PS in PSs:
            if PS[2] in split.keys():
                split[PS[2]].append(PS)
            else:
                split[PS[2]] = [PS]
        SVfull[chromo] = split #delete the list of PS and replace it by this clean dic "ID SV = PSs"
    #
    #
    #get position of PS ID in the header #PS nom called RP !!
    po = getEltPosition(bPSheader, 'RP')
    #print "po", po
    #
    #Add qualities to byPS files
    SVQualDic = {}
    for chromo,SVids in SVfull.iteritems():
        for ID, PSs in SVids.iteritems():
            quality = []
            for PS in PSs:
                #print "PS[18]", PS[18], PS, bPSheader
                qual = dicQual[PS[po]]
                quality.extend(qual)
            qual = round(np.mean(quality),1)
            SVQualDic[chromo + "-" + ID] = qual
            for PS in PSs:
                PS.append(str(qual))
    #
    #Add qualities to bySV files
    #
    for chromo, SVids in SVsimple	.iteritems():
        for sv in SVids:
            sv.append(str(SVQualDic[sv[1]+"-"+sv[2]]))
    #
    #append headers
    bPSheader.append("AvrQual")
    bSVheader.append("AvrQual")
    #
    #byPSfile = open("tumor.chr20.bam_Ulysse_duplications_byPS.2.csv", "w")
    byPSfile = open(byPSfileName, "w")
    byPSfile.write("%s\n" % ';'.join(bPSheader))
    #
    for chromo,SVids in SVfull.iteritems():
        for ID, PSs in SVids.iteritems():
            quality = []
            for PS in PSs:
                byPSfile.write("%s\n" % ';'.join(PS))
    #rien
    #bySVfile = open("tumor.chr20.bam_Ulysse_duplications_bySV.2.csv", "w")
    bySVfile = open(bySVfileName, "w")
    bySVfile.write("%s\n" % ';'.join(bSVheader))
    for chromo, SVids in SVsimple	.iteritems():
        for SV in SVids:
            bySVfile.write("%s\n" % ';'.join(SV))


#--------------------------------------------------------------------------


def cleanAnySVFilePair(bySVfileName, byPSfileName):
    """ function that can take any SV detection pair files (SV and PS), search for sub-groups, 
    filter them out and replace the original files by clean ones """
    #
    ###print "\nStart removing duplicates in files", bySVfileName, "and", byPSfileName
    #open byPS file
    bPStmpMem = csv.reader(open(byPSfileName), delimiter=';', quotechar='"')  
    #print byPSfileName
    bPStmp = []
    for row in bPStmpMem:
        bPStmp.append(row)
    bPSraw = bPStmp[1:]
    bPSheader = bPStmp[0]
    #print "header", bPSheader
    #
    #open bySV file
    bSVtmpMem = csv.reader(open(bySVfileName), delimiter=';', quotechar='"') 
    #bSVtmp = [x.rstrip() for x in bySVfile.readlines()]
    bSVtmp = []
    for row in bSVtmpMem:
        bSVtmp.append(row)
    #
    bSVraw = bSVtmp[1:]
    bSVheader = bSVtmp[0]
    #
    #
    bSV = bSVraw
    bPS = bPSraw
    #
    #
    #make liste of involved chrs
    chrList = list(set([x[1] for x in bSV])) 
    #print chrList
    #
    SVfull = {}
    SVsimple = {}
    for cChr in chrList: #initialize
        SVfull[cChr] = []
        SVsimple[cChr] = []
    #
    #split PS by chr
    for PS in bPS:
       #print PS
        SVfull[PS[1]].append(PS) #fill each key (chro) with PS in the dic
    #
    #split SV by chr
    for iSV in bSV:
        SVsimple[iSV[1]].append(iSV)
    #
    #split PS by SV ID
    for chromo,PSs in SVfull.iteritems():
        split = {}
        for PS in PSs:
            if PS[2] in split.keys():
                split[PS[2]].append(PS)
            else:
                split[PS[2]] = [PS]
        SVfull[chromo] = split #delete the list of PS and replace it by this clean dic "ID SV = PSs"
    #
    #   
    #   
    #
    SVtoKeep=[]
    SVsimpletoKeep=[]
    #
    #get position o PS ID in the header
    po = getEltPosition(bPSheader, 'PS')
    #
    #
    for chromo,SVs in SVfull.iteritems(): #pour chaque chromosome
        SVPSs = []
        SVPSsDic = {}
        for idSV, PSs in SVs.iteritems(): #pour chaque sv et ses PS
            #print PSs
            SVPSs.append([x[po] for x in PSs]) #dresse la liste de toutes les PS (leur ID) impliquee dans une SV
            SVPSsDic[idSV] = [x[po] for x in PSs] #fait dictionnaire: ID_SV=[liste des PS ID] #meme liste que dans SVfull mais avec PS ID seulement
        values = SVPSsDic.values()
        for x in values:
            x.sort()
        #
        #
        #ff = filter_list(map(list,list(set(map(tuple,values))))) #remove sub SV
        #
        #
        #rien    
        tobeRemoved = []
        appendedElements = []
        #
        for i, elt in enumerate(SVPSsDic.values()): #pour chaque SV "i" et ses PS "elt"
            elt.sort()
            #if (elt not in ff) and (elt not in tobeRemovedFull):
            if elt in appendedElements:
                tobeRemoved.append(i)
            else:
                appendedElements.append(elt)
        #rien
        #rien
        listeIDtoRemove = [SVPSsDic.keys()[x] for x in tobeRemoved]
        listeIDtoKeep = [x for x in SVPSsDic.keys() if x not in listeIDtoRemove]
        #rien
        for SVID in listeIDtoKeep:
            SVtoKeep.append(SVs[SVID])
        #rien
        for nSV in SVsimple[str(chromo)]:
            if nSV[2] in listeIDtoKeep:
                SVsimpletoKeep.append(nSV)
    #do a "uniq of byPS"
    byPSListe = []
    for SV in SVtoKeep:
        for i, PS in enumerate(SV):
            #print PS
            byPSListe.append(tuple(PS))
    byPSListUniq = list(set(byPSListe))
    #rien 
    #rien
    #rien
    #byPSfile = open("tumor.chr20.bam_Ulysse_duplications_byPS.stats2.csv", "w")
    byPSfile = open(byPSfileName, "w")
    byPSfile.write("%s\n" % ';'.join(bPSheader))
    #for SV in SVtoKeep:
    #    for PS in SV:
    for PS in byPSListUniq:
        byPSfile.write("%s\n" % ';'.join(PS))
    #rien
    #bySVfile = open("tumor.chr20.bam_Ulysse_duplications_bySV.stats2.csv", "w")
    bySVfile = open(bySVfileName, "w")
    bySVfile.write("%s\n" % ';'.join(bSVheader))
    for SV in SVsimpletoKeep:
        bySVfile.write("%s\n" % ';'.join(SV))
    #rien
    #print "Finished removing duplicates in files", bySVfileName, "and", byPSfileName, "\n"

     
#--------------------------------------------------------------------------

def cleanAnySVFilePairOLD(bySVfileName, byPSfileName):
	""" function that can take any SV detection pair files (SV and PS), search for sub-groups, 
	filter them out and replace the original files by clean ones """
	
	###print "\nStart removing duplicates in files", bySVfileName, "and", byPSfileName
	bPStmpMem = csv.reader(open(byPSfileName), delimiter=';', quotechar='"') 
	
	bPStmp = []
	for row in bPStmpMem:
		bPStmp.append(row)
 	bPSraw = bPStmp[1:]
 	#print bPSraw
	bPSheader = bPStmp[0]

	bSVtmpMem = csv.reader(open(bySVfileName), delimiter=';', quotechar='"') 	
	#bSVtmp = [x.rstrip() for x in bySVfile.readlines()]
  	bSVtmp = []
	for row in bSVtmpMem:
		bSVtmp.append(row)

 	bSVraw = bSVtmp[1:]
	bSVheader = bSVtmp[0]
	#

	bSV = bSVraw
	bPS = bPSraw
	#
	#
	chrList = list(set([x[1] for x in bSV]))
 	#print chrList
	#
	SVfull = {}
	SVsimple = {}
	for cChr in chrList:
		SVfull[cChr] = []
		SVsimple[cChr] = []
	#
	#split PS by chr
	for PS in bPS:
     		#print PS
		SVfull[PS[1]].append(PS)
	#
	#split SV by chr
	for iSV in bSV:
		SVsimple[iSV[1]].append(iSV)
	#
	#split PS by SV ID
	for chromo,PSs in SVfull.iteritems():
		split = {}
		for PS in PSs:
			if PS[2] in split.keys():
				split[PS[2]].append(PS)
			else:
				split[PS[2]] = [PS]
		SVfull[chromo] = split
	#
	#   
	#   
	#
	SVtoKeep=[]
	SVsimpletoKeep=[]
	#
	for chromo,SVs in SVfull.iteritems(): #pour chaque chromosome
		SVPSs = []
		SVPSsDic = {}
		for idSV, PSs in SVs.iteritems(): #pour chaque sv et ses PS
			#print PSs
			SVPSs.append([x[17] for x in PSs]) #dresse la liste de toutes les PS impliquee dans SV
			SVPSsDic[idSV] = [x[17] for x in PSs] #fait dictionnaire: ID_SV=[liste des PS ID]
		ff = filter_list(SVPSsDic.values()) #remove sub SV
		#print len(ff)
		#rien	
		tobeRemoved = []
		for i, elt in enumerate(SVPSsDic.values()): #pour chaque SV "i" et ses PS "elt"
			#print "ELT", elt #TODO: notworking here: elt is a liste of PS ID
#			for singleElt in elt:
#				if singleElt not in ff:
#					tobeRemoved.append(i)
			if elt not in ff:
					tobeRemoved.append(i)
     
		#rien
		listeIDtoRemove = [SVPSsDic.keys()[x] for x in tobeRemoved]
		listeIDtoKeep = [x for x in SVPSsDic.keys() if x not in listeIDtoRemove]
		#rien
		for SVID in listeIDtoKeep:
			SVtoKeep.append(SVs[SVID])
		#rien
		for nSV in SVsimple[str(chromo)]:
			if nSV[2] in listeIDtoKeep:
				SVsimpletoKeep.append(nSV)

	#do a "uniq of byPS"
	byPSListe = []
	for SV in SVtoKeep:
		for i, PS in enumerate(SV):
			#print PS
			byPSListe.append(tuple(PS[1]))
	byPSListUniq = list(set(byPSListe))
 
 
 
	#byPSfile = open("tumor.chr20.bam_Ulysse_duplications_byPS.stats2.csv", "w")
	byPSfile = open(byPSfileName, "w")
	byPSfile.write("%s\n" % ';'.join(bPSheader))
	#for SV in SVtoKeep:
	#	for PS in SV:
	for PS in byPSListUniq:
		byPSfile.write("%s\n" % ';'.join(PS))

	#bySVfile = open("tumor.chr20.bam_Ulysse_duplications_bySV.stats2.csv", "w")
	bySVfile = open(bySVfileName, "w")
	bySVfile.write("%s\n" % ';'.join(bSVheader))
	for SV in SVsimpletoKeep:
		bySVfile.write("%s\n" % ';'.join(SV))

	#print "Finished removing duplicates in files", bySVfileName, "and", byPSfileName, "\n"

#--------------------------------------------------------------------------

def reOrder(bySVfileName, byPSfileName, d):
    """ function that can add mean PS SV quality in any SV detection file pair (SV and PS) """
    #
    #
    #open byPS file
    bPStmpMem = csv.reader(open(byPSfileName), delimiter=';', quotechar='"')
    bPStmp = []
    for row in bPStmpMem:
        bPStmp.append(row)
    bPSraw = bPStmp[1:] #remove header
    bPSheader = bPStmp[0]
    #
    #
    #open bySV file
    bSVtmpMem = csv.reader(open(bySVfileName), delimiter=';', quotechar='"')
    bSVtmp = []
    for row in bSVtmpMem:
        bSVtmp.append(row)
    bSVraw = bSVtmp[1:] #remove header
    bSVheader = bSVtmp[0]
    #
    #creast short names
    bSV = bSVraw
    bPS = bPSraw
    #
    #
    #make liste of involved chrs
    chrList = list(set([x[1] for x in bSV]))
    #
    #
    SVfull = {} #dicos for by PS
    SVsimple = {} #dicos for by SV
    SVsimplebyID={}
    for cChr in chrList: #initialize the dicos
        SVfull[cChr] = []
        SVsimple[cChr] = []
        SVsimplebyID[cChr] = {}
    #
    #split PS by chr
    for PS in bPS:
       SVfull[PS[1]].append(PS) #fill each key (chrom) with PS in the dic
    #
    #split SV by chr
    for iSV in bSV:
        SVsimple[iSV[1]].append(iSV)
    #split by SV ID
    for chromo, SVids in SVsimple    .iteritems():
        for sv in SVids:
            SVsimplebyID[chromo][sv[2]] = sv
    #
    #
    #split PS by SV ID
    for chromo,PSs in SVfull.iteritems():
        split = {}
        for PS in PSs:
            if PS[2] in split.keys():
                split[PS[2]].append(PS)
            else:
                split[PS[2]] = [PS]
        SVfull[chromo] = split #delete the list of PS and replace it by this clean dic "ID SV = PSs"
    #
    #
    #get strand of reads from the header
    po1 = getEltPosition(bPSheader, 'str1')
    po2 = getEltPosition(bPSheader, 'str2')
    coordL = getEltPosition(bPSheader, 'pos1')
    coordR = getEltPosition(bPSheader, 'pos2')
    #
    #Split the SV in the byPS dicos (SVfull) by groups of orientations: c1 and c2
    for chromo,SVids in SVfull.iteritems():
        for ID, PSs in SVids.iteritems():
            c1, c2 = [], []
            flag = {}
            n = 0
            for PS in PSs:
                ori = PS[po1]+PS[po2]
    #
                if n == 0: #pour le premiere PS
                    c1.append(PS)
                    n+=1
                    flag[ori] = 1
    #
                elif ori not in flag.keys(): #si on a change d orientation
                    c2.append(PS)
                else: #c est c est tjr la meme que la premiere
                    c1.append(PS)
    #
            SVfull[chromo][ID] =  { 'c1':c1, 'c2':c2 } #remplacer PSs pour ce dicos c1: PS groupe 1 , c2:PS groupe 2
    #
    #check if orientation is correct, if not reverse
    #print "toto"
    for chromo,SVids in SVfull.iteritems():
    #
        for ID, groups in SVids.iteritems():
            #print groups
            c1 = groups['c1']
            c2 = groups['c2']
            #print c1
            c1_coordsL = int(np.mean([int(x[coordL]) for x in c1]))
            c1_coordsR = int(np.mean([int(x[coordR]) for x in c1]))
            c2_coordsL = int(np.mean([int(x[coordL]) for x in c2]))
            c2_coordsR = int(np.mean([int(x[coordR]) for x in c2]))
            
            
            if abs(c1_coordsR-c2_coordsR) < 2*d and abs(c1_coordsL-c2_coordsL) >= 2*d:
                
                #reverse the by SV file
                svraw = SVsimplebyID[chromo][ID]
                #print "svraw", svraw
                PS = svraw #.split(";")
                PS[1] = PS[1].split("-")[1] + "-" + PS[1].split("-")[0]
                PS[4], PS[5] = PS[5], PS[4]
                PS[6], PS[7], PS[8], PS[9], PS[10], PS[11], PS[12], PS[13] = \
                    PS[10], PS[11], PS[12],PS[13], PS[6], PS[7], PS[8], PS[9]
                #sv = ";".join(";")
                SVsimplebyID[chromo][ID] = PS
                PS= ""
                
                #reverse the byPS file
                g1, g2 = [], []
                for rawPS in c1: #reverse cluster with orientation one
                    PS = rawPS 
                    PS[1] = PS[1].split("-")[1] + "-" + PS[1].split("-")[0]
                    PS[4], PS[5] = PS[5], PS[4]
                    PS[6], PS[7], PS[8], PS[9], PS[10], PS[11], PS[12], PS[13] = \
                    PS[10], PS[11], PS[12],PS[13], PS[6], PS[7], PS[8], PS[9]
                    PS[21], PS[22], PS[23], PS[24], PS[25], PS[26] = \
                    PS[24], PS[25], PS[26], PS[21], PS[22], PS[23]
                    g1.append(PS)
                SVfull[chromo][ID]['c1']= g1
                for rawPS in c2: #reverse the other cluster
                    PS = rawPS 
                    PS[1] = PS[1].split("-")[1] + "-" + PS[1].split("-")[0]
                    PS[4], PS[5] = PS[5], PS[4]
                    PS[6], PS[7], PS[8], PS[9], PS[10], PS[11], PS[12], PS[13] = \
                    PS[10], PS[11], PS[12],PS[13], PS[6], PS[7], PS[8], PS[9]
                    PS[21], PS[22], PS[23], PS[24], PS[25], PS[26] = \
                      PS[24], PS[25], PS[26], PS[21], PS[22], PS[23]
                    g2.append(PS)
                SVfull[chromo][ID]['c2']= g2
    #
    #
    #byPSfile = open("tumor.chr20.bam_Ulysse_insertions_byPS.2.csv", "w")
    byPSfile = open(byPSfileName, "w")
    byPSfile.write("%s\n" % ';'.join(bPSheader))
    #
    for chromo,SVids in SVfull.iteritems():
        for ID, groups in SVids.iteritems():
            for cluster, PSs in groups.iteritems():
                for PS in PSs:
                    byPSfile.write("%s\n" % ';'.join(PS))
    #rien
    #bySVfile = open("tumor.chr20.bam_Ulysse_insertions_bySV.2.csv", "w")
    bySVfile = open(bySVfileName, "w")
    bySVfile.write("%s\n" % ';'.join(bSVheader))
    for chromo, SVids in SVsimplebyID.iteritems():
        for ID, SV in SVids.iteritems():
            bySVfile.write("%s\n" % ';'.join(SV))
            

#--------------------------------------------------------------------------
def getN(fileIn, tipe, list_chr_real, prefix):
    """ get the real ndisc """
    
    if tipe == "INV":
        num = 2
    else:
        num = 1
    if ''.join(list_chr_real)=='all':
        with open(fileIn, "r") as f:
            line = f.readline().rstrip().split(" ")
            dicos = {}
            while line:
                line = f.readline().rstrip().split(" ") 
                if line == ['']:
                       break
                dicos[line[0]] = line[num]
        
        return sum([int(dicos[x]) for x in dicos.keys()])
    else:
        #print num, list_chr_real, tipe
        if list_chr_real[0] == "noGood666":
            return 100
        else:
            with open(fileIn, "r") as f:
                lineT = f.readline().rstrip().split(" ")
                dicos = {}
                while lineT:
                    lineT = f.readline() #.rstrip().split(" ") 
                    #print "lineT:", lineT, repr(lineT), type(lineT)
                    if lineT.rstrip() == '':
                    	   break
                    if lineT:
                        line = lineT.rstrip().split(" ")
                        dicos[line[0]] = line[num]
            #print dicos
            #return sum([int(dicos[prefix+ str(x)]) for x in list_chr_real])
            return sum([int(dicos[str(x)]) for x in list_chr_real])

#--------------------------------------------------------------------------
def getMinPSvalue(minPS):
    """Get the minimum value in the dictionnary minPS containing min nb of PS
    for each interSV type"""
    liste = [minPS["ins"]*2, minPS["tn"], minPS["tr"]]
    return min(liste)

#--------------------------------------------------------------------------
def add2Coord(dico, cle, valeur):
    if cle not in dico:
        dico[cle] = [valeur]
    else:
        dico[cle].append(valeur)

#--------------------------------------------------------------------------
def Extremites(xsome, liste):
    """ Definit les extremites d'un groupe de PS"""
    # initialise une liste avec coord min et max pour chr1 et chr2
    if liste != []:
        l = xsome[liste[0]]
        deb1 = l[3]
        fin1 = l[3]
        deb2 = l[6]
        fin2 = l[6]

    #determine les min et max des coordonnees sur les deux chromosomes
        for i in range(1, len(liste)):
            l = xsome[liste[i]]
            p1 = l[3]
            p2 = l[6]
            deb1 = min(p1, deb1)
            fin1 = max(fin1, p1)
            deb2 = min(p2, deb2)
            fin2 = max(fin2, p2)

    else:
        deb1 = fin1 = deb2 = fin2 = -1

    return deb1, fin1, deb2, fin2

#------------------------------------------------------------------------------
def sizeSV(liste, select, typesv):
    """returns min, max of SV and median of Insert Size"""
    taille=[]
#    start = select[liste[0]][3]
#    end = select[liste[0]][6]
#    for i in liste:
#        taille.append(int(select[i][7]))
#        if select[i][6] > end :
#            end = select[i][6]
#        if select[i][3] < start :
#            start = select[i][3]
#
#    return np.median(taille), start, end
    #print liste
    
    #initialize
    minl = select[liste[0]][3]
    maxl = select[liste[0]][3]
    minr = select[liste[0]][6]
    maxr = select[liste[0]][6]    
    
    for i in liste:
        taille.append(int(select[i][7]))
        if select[i][6] > maxr :
           maxr  = select[i][6]
        elif select[i][6] < minr :
           minr  = select[i][6]
        
        if select[i][3] < minl :
            minl = select[i][3]
        elif select[i][3] > maxl :
            maxl  = select[i][3]

    if typesv == "del":
        return np.median(taille), maxl, minr
    else:
        return np.median(taille), minl, maxr


#------------------------------------------------------------------------------
def write_byPS(out, tup, deletion=False):
    """ Write a line of the _byPS file"""
    #print tup
    if deletion:
        out.write("\"%s\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%s\";\
\"%s\";\
\"%s\";\
\"%d\";\
\"%s\";\
\"%s\";\
\"%d\";\
\"%d\"\n" % tup)
    else:
         out.write("\"%s\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%s\";\
\"%s\";\
\"%s\";\
\"%s\";\
\"%d\";\
\"%s\";\
\"%s\";\
\"%d\";\
\"%d\"\n" % tup)

#   % (manip, xsome, ndup, nbps, nbps, -1, start, end, end-start, "NA", -1, -1, -1, "NA", sizemin, sizemax, -1, k[1], k[2], k[3], k[4], k[5], k[6], k[7]))
#   %s;\ manip
#   %s;\ pair
#   %d;\ ID
#   %d;\ nbPS
#   %d;\ nbA
#   %s;\ nbB
#   %d;\ left_borderA
#   %d;\ right_border_A
#   %d;\ deltaA
#   %s;\ cen_strA
#   %d;\ left_borderB
#   %d;\ right_borderB
#   %d;\ deltaB
#   %s;\ cen_strB
#   %d;\ SV_size_min
#   %d;\ SV_size_max
#   %s;\ p-value
#   %s;\ PS
#   %s;\ str1
#   %s;\ chr1
#   %d;\ pos1
#   %s;\ str2
#   %s;\ chr2
#   %d;\ pos2
#   %d\ Insert size

#------------------------------------------------------------------------------
def write_bySV(out, tup, deletion=False):
    """ Write a line of the SV file"""
    #print tup
    if deletion:
        out.write("\"%s\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%s\"\n" % tup)
    else:
        out.write("\"%s\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%d\";\
\"%d\";\
\"%s\";\
\"%s\"\n" % tup)

#   %s;\ manip
#   %s;\ pair
#   %d;\ ID
#   %d;\ nbPS
#   %d;\ nbA
#   %d;\ nbB
#   %d;\ left_borderA
#   %d;\ right_border_A
#   %d;\ deltaA
#   %s;\ cen_strA
#   %d;\ left_borderB
#   %d;\ right_borderB
#   %d;\ deltaB
#   %s;\ cen_strB
#   %d;\ SV_size_min
#   %d;\ SV_size_max
#   %s;\ p-value



    

 
