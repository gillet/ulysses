#!/usr/bin/env python
# -*- coding: utf-8 -*-




import os, sys
try:
    import pysam as pys
except ImportError:
    print "\tError: Pysam is not installed. \n\
    Please see: https://github.com/pysam-developers/pysam\n\n"
    sys.exit()
try:
    import numpy as np
    import numpy.ma as ma
except ImportError:
    print "\tError: NumPy is not installed. \n\
    Please see: http://www.numpy.org/\n\n"
    sys.exit()

import argparse


import Ulysse_utils as U
import Ulysse_stats as Ualex




#-----------------------------------------------------------------------------------

def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """

    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m

#-----------------------------------------------------------------------------------
def getori(nb):
    """Convert (True/False) result from read.is_reverse to Ulysse read 
    orientation format (-/+) """
    ori = "-"
    if not nb :
        ori = "+"
    return ori

#-----------------------------------------------------------------------------------

def getReadsInfos(fil, chromInfo, nbMax, RL, insertL, insertS, typeps, DUP):
    """Reads reads from the BAM file for the given chromosome. Number of reads
    read is proportionnal to the length of chromosome. Lecture begins at the 
    middle of the chromosome until the maxNbReads is reached
    """
    chromN = chromInfo[0]
    chromL = chromInfo[1]
    genomLen = chromInfo[2]
    
    maxNbReads = int(nbMax * chromL / genomLen)
    with pys.Samfile(fil,'rb') as bam:
        cp = 0
        for read1 in bam.fetch(chromN, chromL/2, chromL-1):
            cp +=1
            #print cp
            if cp >= maxNbReads:
                bam.close()
                return 1
                
            if U.is_number(read1.alen):
                RL.append(read1.alen)

            if read1.is_read1 and not read1.is_unmapped and not read1.mate_is_unmapped \
            and read1.tid!=-1 and read1.rnext!=-1:
                
                chr1 = bam.getrname(read1.tid)
                chr2 = bam.getrname(read1.rnext)
                
                #if not read1.is_duplicate and read1.is_proper_pair and\
                if not read1.is_duplicate and read1.is_paired and\
                read1.tid >= 0 and read1.rnext >= 0:
                    chr1 = bam.getrname(read1.tid)
                    chr2 = bam.getrname(read1.rnext)
                    ori1, chr1, pos1, ori2, chr2, pos2 = U.getCoord(read1, 
                                                                    chr1, chr2)
                    if chr1 == chr2 and pos1 != pos2 :
                        
                        insertL.append(abs(read1.tlen))
                        typeps[ori1+ori2] += 1
                        insertS.append(abs(read1.tlen))
                        orient = ''.join([ori1, ori2])
                        if orient == '-+' or orient == '+-':
                            DUP[orient].append(abs(pos1-pos2))
            if cp >= nbMax:
               return 1
        if cp==0:
            print "Warning: sampling entire chromosome", chromN, "reads from the start"
            for read1 in bam:
                    
                if U.is_number(read1.alen):
                    RL.append(read1.alen)
    
                if read1.is_read1 and not read1.is_unmapped and not read1.mate_is_unmapped \
                and read1.tid!=-1 and read1.rnext!=-1:
                    
                    chr1 = bam.getrname(read1.tid)
                    chr2 = bam.getrname(read1.rnext)
                    
                    #if not read1.is_duplicate and read1.is_proper_pair and\
                    if not read1.is_duplicate and read1.is_paired and\
                    read1.tid >= 0 and read1.rnext >= 0:
                        chr1 = bam.getrname(read1.tid)
                        chr2 = bam.getrname(read1.rnext)
                        ori1, chr1, pos1, ori2, chr2, pos2 = U.getCoord(read1, 
                                                                        chr1, chr2)
                        if chr1 == chr2 and pos1 != pos2 :
                            
                            insertL.append(abs(read1.tlen))
                            typeps[ori1+ori2] += 1
                            insertS.append(abs(read1.tlen))
                            orient = ''.join([ori1, ori2])
                            if orient == '-+' or orient == '+-':
                                DUP[orient].append(abs(pos1-pos2))
                            #print read1.qname, read1.tlen
                if cp >= 100000000:
                   return 1            
    #return RL, insertL, insertS, typeps, DUP  
    return 1

#-----------------------------------------------------------------------------------
def Read10MBAM(fil, chrDicos) :
    """ Read first 1000000 PS in BAM file to estimate mean, mad and median. 
    The input file is a BAM file already sorted and indexed by samtools sort
    and samtools index """
    
   
    insertL = []
    insertS = []
    typeps = {"++":0,"--":0,"+-":0,"-+":0}
    DUP = {'-+':[], '+-':[]}
    RL = []
    
    nbMax = 2000000
    genomLen = sum(chrDicos.values())
    for chrom in chrDicos:
        chromInfo = (chrom, chrDicos[chrom], genomLen)
        #print fil, chromInfo , nbMax, len(insertL)
        p =getReadsInfos(fil, chromInfo, nbMax, RL, insertL, insertS, typeps, DUP)

    if typeps["-+"] > typeps["+-"]:
        type_in = "MP"
        DUP = DUP['+-']
    else:
        type_in = "PE"
        DUP = DUP['-+']
    
    insertL = np.array(insertL)
    mean = np.mean(insertL)
    median = np.median(insertL)
    mad = MAD(insertL)
    sd = np.std(insertL)
    RL = int(np.mean(RL))
    
    DUP = [x for x in DUP if x < median-3*mad]
    if median > 1000:
        if DUP:        
            DUPcutoff = int(np.median(DUP) + 6*MAD(DUP))
        else:
            DUPcutoff = 0
    else:
        DUPcutoff = 0
    del insertL
    
    return mean, median, mad, sd, type_in, insertS, str(RL), DUPcutoff
    
#----------------------------------------------------------------------------
def isDup(ori1, ori2, inverse, length, median, chr1, pos1, chr2, pos2, subtel, DUPcutoff):
    """Return True if a given PS is candidate for a Duplication"""
    
    ps1sub = U.Subtelo(subtel, chr1, pos1)
    ps2sub = U.Subtelo(subtel, chr2, pos2)
    if (ori1+ori2 == inverse) and (length > median) and (ps1sub + ps2sub <= 1) \
        and length > DUPcutoff:
        return True
    else:
        return False
#----------------------------------------------------------------------------
def isDel(ori1, ori2, orientation, length, sizemindel):
    """Return True if a given PS is candidate for a Deletion"""

    if (ori1+ori2 == orientation) and (length > sizemindel) :
        return True
    else:
        return False
#----------------------------------------------------------------------------
def isSmallDup(ori1, ori2, orientation, length, sizemaxdup):
    """Return True if a given PS is candidate for a Small Duplication"""
    
    if (ori1+ori2 == orientation) and (length < sizemaxdup) :
        return True
    else:
        return False
#----------------------------------------------------------------------------
def isSVCandidate(rang, name, ori1, chr1, pos1, ori2, chr2, pos2, length, ps_type,
                  median, mad, mean, subtel, params, nbcand,
                  insert, DUPcutoff):
    """Return true if a given PS is candidate for any SV type"""

    selec = False
    if ps_type == "MP":
        orientation = "-+"
        inverse = "+-"
    else:
        orientation = "+-"
        inverse = "-+"
        
    if ori1+ori2 == inverse:
        selec = isDup(ori1, ori2, inverse, length, median,  chr1, pos1, chr2, 
                pos2, subtel, DUPcutoff)
        if selec and rang == 1:
            #nbcand["DUP"] += 1
            nbcand["DUP"][chr1] = nbcand["DUP"].get(chr1, 0) + 1
            
    elif ori1+ori2 == orientation:
        sizemindel = median + float(params["n"])*mad
        selec = isDel(ori1, ori2, orientation, length, sizemindel)
        if selec and rang == 1:
            nbcand["DEL"] += 1
            insert["del"].append(length)
        else:
            sizemaxdup = median - float(params["n"])*mad
            if sizemaxdup > 0:
                selec = isSmallDup(ori1, ori2, orientation, length, sizemaxdup)
                if selec and rang == 1:
                    nbcand["SINS"] += 1
                    insert["ins"].append(length)

    elif ((ori1+ori2 == "--") or (ori1+ori2 == "++")):
        selec = True
	if rang == 1:
            nbcand["INV"][chr1] = nbcand["INV"].get(chr1, 0) + 1
        #nbcand["INV"] += 1
        

    return selec
        
#----------------------------------------------------------------------------
def ReadBAM(fil, mapq, median, mad, mean, sd, subtel, params, nbcand, 
            ps_type, DUPcutoff) :
    """
	Read Files in BAM format and generate a format one line per PS to be
     compatible with Ulysse detection tools.  
	The input file is a BAM file already sorted and indexed by samtools sort and samtools index
    """
    
    
    insert = {"del":[], "ins":[]}
    with pys.Samfile(fil,'rb') as bam:
        
        
        diff = pys.Samfile(params["in"]+"_xdiff","wb", template = bam)
        filx = {}
        for cx in bam.references:
            filx[cx] = pys.Samfile(params["in"]+"_"+cx,"wb", template = bam)
            
        ## PREPROCESSING: get reformat.table and coverage estimation while reading BAM file ######################################
        chrDicos, genomelength = Ualex.getChrList(params["in"], False) #F can be an opened pysam file or just the path or the file
        
        #make dic for coverage and reformat.table list
        reformat_dist=[]
        compteur = 0
            
        for read1 in bam:
            compteur+=1
            if not read1.is_duplicate and \
               not read1.mate_is_unmapped and not read1.is_unmapped and \
               read1.tid >= 0 and read1.rnext >= 0:
                if compteur%1000000 == 0 and compteur != 0:
                    print "Processed", '{0:,}'.format(compteur), "reads"
                    
                if read1.is_read1:
                    rang = 1
                else:
                    rang = 2

                chr1=bam.getrname(read1.tid)
                chr2=bam.getrname(read1.rnext)
                ori1, chr1, pos1, ori2, chr2, pos2 = U.getCoord(read1, 
                                                                chr1, chr2)
                Ualex.updateReformat(read1, reformat_dist, bam)
                if chr1 == chr2 and pos1 != pos2:
                    length = abs(read1.tlen)

                    selected = isSVCandidate(rang, read1.qname, ori1, chr1, 
                                             pos1, ori2, chr2, pos2,
                                             length, ps_type, median, mad,
                                             mean, subtel, params,
                                             nbcand, insert, DUPcutoff)
                    if selected:
                        filx[chr1].write(read1)

                elif chr1 != chr2:
                    pair = [chr1, chr2]
                    pair.sort()
                    nbcand["INTER"]["-".join(pair)] = \
                    nbcand["INTER"].get("-".join(pair), 0) + 1
                    diff.write(read1)

        bam.close()        
        diff.close()
        for cx in filx:
                filx[cx].close()
    return chrDicos, genomelength, reformat_dist, insert
#-------------------------------------------------------------------------
def createParamFile(args):
    """create the parameter file if no one is specified in argument"""
    
    bamf = os.path.abspath(args.bamfile)
    paramf =  os.path.abspath(args.p)
    with open(paramf,"w") as file:
        file.write("in="+bamf+"\n")
        #file.write("mapq="+str(args.mapq)+"\n")
        file.write("mapq=35\n")
        file.write("mapqx=20\n")        
        out = bamf
        if args.out:
          out = args.out
        file.write("out="+out+"_Ulysses\n")
        file.write("vcf=False\n")
        file.write("stats="+bamf+"_stats.txt\n")
        file.write("NSV=10000\n")
        file.write("fdr=0.01\n")
        file.write("range=all\n")
        file.write("n="+str(args.n)+"\n")
        file.write("annotation=None\n")
        file.write("######### Informations for annotation file\n") #a virer car que GFF?
        file.write("field_chr=1\n")
        file.write("field_type=3\n")
        file.write("field_start=4\n")
        file.write("field_end=5\n")
        file.write("field_sep=tb\n")	#a virer
#-------------------------------------------------------------------------
def writeStatsFile(params, type_in, prefix, mean, median, mad, sd, nbcand, 
                   genomelength, chrDicos, RL):
    """create the statistics file"""
                       
    with open(params["stats"],"w") as stats:

        stats.write("#Estimation on the first 100000 PS\n")
        stats.write("PS_type : "+type_in+"\n")
        stats.write("Read_length : "+RL+"\n")
        stats.write("Chromosome_prefix : "+prefix+"\n")
        stats.write("Mean : %.2f\n" %(mean))
        stats.write("Median : %.2f\n" %(median))
        stats.write("MAD : %.2f\n" %(mad))
        stats.write("stdev : %.2f\n" %(sd))
        for tipe in nbcand:
            if tipe not in ["INTER", "DUP", "INV"] :
                stats.write("n%s : %d\n" % (tipe, nbcand[tipe]))

        stats.write("nDUP : %d\n" % (sum(nbcand["DUP"].values())))
        stats.write("nINV : %d\n" % (sum(nbcand["INV"].values())))
        stats.write("nInter : %d\n" % (sum(nbcand["INTER"].values())))
        stats.write("genome_length : %d\n" % (genomelength))
        for chrx in chrDicos:
            stats.write("length %s %d\n" % (chrx, chrDicos[chrx]))

# Write number of discordant ps (per chromosomes and per pairs)
    with open(params["in"]+"_PScandidate_intra.csv","w") as countPS:
        countPS.write("%s\n" % ("chr DUP INV"))
        for chromo in chrDicos.keys():
            countPS.write("%s %d %d\n" % (chromo, nbcand["DUP"].get(chromo,0),\
                                                nbcand["INV"].get(chromo,0))) 

           
# Write number of discordant ps (per chromosomes and per pairs)
    with open(params["in"]+"_PScandidate_per_Xsome_pairs.csv","w") as countPS:
        countPS.write("%s %s\n" % ("V1", "V2"))        
        for pair in nbcand["INTER"]:
            countPS.write("%s %d\n" % (pair, nbcand["INTER"][pair]))

    #print params["stats"], "created"

#-------------------------------------------------------------------------
def parser():
    """Parse command line argument"""
    message = "\tULYSSES V1.0 ReadBam, BAM/SAM file \
pre-processing before detection of Structural Variations"
    
    print "\n\n\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n"    
    parser = argparse.ArgumentParser(description=message)
    print "\t\t\t\t---  ULYSSES ReadBAM 20141023  ---     \n\n"
    print "  Ulysses: Accurate detection of rare structural variations from high coverage genome sequencing\n"
    print "\t\t Alexandre Gillet, Hugues Richard, Gilles Fischer and Ingrid Lafontaine\n"
    print "\t\t\t http://www.lcqb.upmc.fr/ulysses\n\n"
    print "\t\t\t\t Copyright UPMC - CNRS\n\n"
    



    parser.add_argument("bamfile", help='Name of the original BAM/SAM file')
    parser.add_argument("-p", metavar='parameter_file_name',
                        help='file containing all parameter values. default \
ulysses_params', default='ulysses_params')

    parser.add_argument("-n", metavar='mulitplicative factor for MAD and d(n)', type=int,
                        default=6, help='Interdistance and insert size consistency factor')
    
#    parser.add_argument("-mapq", metavar='mapping_quality_threshold',
#                        type=int, default=20, help='minimal mapping quality of reads')

    parser.add_argument("-out", metavar='prefix_of_output_file',
                        help='if not given, library_file_name_Ulysses')

    parser.add_argument("-stats", metavar='statistics_file_name',
                        help='name of statistics file created during library processing')

    args = parser.parse_args()
    
    print message+"\n\nBAM file :", args.bamfile
    print "Parameter file :", args.p,"\n\n"
    print "oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n"

    if args.p == "ulysses_params" or not os.path.exists(os.path.abspath(args.p)) :
        createParamFile(args)
        #createParamFile(os.path.abspath(args.bamfile), os.path.abspath(args.p),
        #                 args)
                                       
            
    return os.path.abspath(args.bamfile), os.path.abspath(args.p)
    

#-------------------------------------------------------------------------
if (__name__) == "__main__":
    """1. Estimate stats for PS library (mean, median, mad) and deduce lim 
    parameters for removing PE in a Mate Pairs library on the first 100000 PS.
    2. Write Candidate PS on text files    
    1st arg = BAM file 2nd arg = run_info file"""
    
    ok = False
    installed = False
        
    bamfile, paramfile = parser()
    if os.path.exists(bamfile):
        os.system("date")
        if not os.path.exists(bamfile+".bai"):
            print "Indexing file "+bamfile
            pys.index(bamfile)
            if not os.path.exists(bamfile+".bai"):
                print "\n*** Warning:", os.path.split(bamfile)[1], "is not sorted"
            else:
                print "\n",os.path.split(bamfile)[1], "indexed\n"
        
        #get list of all chromosomes
        chrDicos, genomelength = Ualex.getChrList(bamfile, True)
        
        mean, median, mad, sd, type_in, insertS, RL, DUPcutoff = Read10MBAM(bamfile, chrDicos)
        os.system("date")
        
        


        params = U.get_run_info(paramfile)

        chrnames = chrDicos.keys()
        chrnames.sort()
        
                
        params["ps_type"] = type_in
            
        subtelo=U.get_subtelo_limits(params)
        
        del insertS
            
       #Complete lecture of BAM file
   
        nbcand = {"INV":{}, "DEL":0, "DUP":{}, "SINS":0, "INTER":{}}
                  
        chrDicos, genomelength, reformat_dist, insert = ReadBAM(params["in"],
                                                        params["mapq"], median,
                                                        mad, mean, sd, subtelo,
                                                        params, nbcand, 
                                                        params["ps_type"], 
                                                        DUPcutoff)
        prefix = U.getcommonstart(chrnames)
        writeStatsFile(params, type_in, prefix, mean, median, mad, sd, nbcand, 
                       genomelength, chrDicos, RL)
     
        #Ecriture du fichier reformat pour calcul stats de deletions
        fName, msg = Ualex.reduceReformat(reformat_dist, params["in"], params["out"])
        print msg
        


        os.system("date") 
        
    else:
        print "Error : ", bamfile, "doesn't exist"
