# -*- coding: utf-8 -*-



import os
import sys
try:
    import pysam
except ImportError:
    print "\tError: Pysam is not installed. \n\
    Please see: https://github.com/pysam-developers/pysam\n\n"
    sys.exit()
import random
import Ulysse_utils as U
from subprocess import Popen, PIPE, call

########################################################################################################################################################
def runRForClusterThreshold(tipe, ISmean, ISmed, ISmad, ISsd, genomeLength, 
                            nDiscordant, numberOfInterChromosomalFile, 
                            seuil_cluster, script, functionsR, list_chr_length, 
                            list_chr_names, detectionFile, bam_input_name, 
                            dist_table, average_cov, sval, seuil_fdr, dirWork,
                            RL, n_param, numberOfIntraChromosomalPSFile):
    """ 2 functions:
	1- detectionFile="None": Computes the #PS threshold for different SV types (DUP, INV, TR, INS)
	2- detectionFile="aExistingDetectionFile": add a score to each SV based on expected number of SVs
    """
    
    
    
    
    #Compile C code for R (to speed up pIS and dn calculation)
    f = os.path.dirname(os.path.realpath(__file__))+"/probamaxdist.so"
    if not os.path.isfile(f):
        com = "R CMD SHLIB " + os.path.dirname(os.path.realpath(__file__))+"/probamaxdist.c"
        print "\nUlysses will now compile C code for R"    
        print com
        os.system(com)
        print "Compilation Done\n"

    param = [ISmean, ISmed, ISmad, ISsd, nDiscordant, genomeLength, numberOfInterChromosomalFile, \
    seuil_cluster, functionsR, tipe, list_chr_length, list_chr_names, detectionFile, dist_table, \
    sval, seuil_fdr, RL, n_param, numberOfIntraChromosomalPSFile, f]

    #Popen will fail if special characters are present in chromosome names
    param[11] = '"'+param[11]+'"'

    param_str = map(str, param)
    com = ["Rscript", str(script)]
    com.extend(param_str)
    commande = " ".join(com)
        
    proc= Popen(commande, stdout=PIPE, stderr=PIPE, shell=True)
    procOutput = proc.communicate()

    #print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
    #print procOutput
    #print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "

    output = procOutput[0]
    error = procOutput[1]    

#    print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
#    print "output:"
#    print output
#    print "error:"
#    print error
#    print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "


    if proc.returncode != 0:
        #if len(procOutput[-1])>40:
        #    messag1= "Execution of R script was stopped for the following reason:\n"
        #    messag=messag1+procOutput[-1]
        #else:
        
        messag = "Execution of R script stopped for unexpected reasons\n"
        messag = error
        return 1, messag,1


    if error!="":
         return 1, error+"\n"+output,1

    sSplitLines = filter(None, output.split("\n"))
    sSplitVal = [x.replace("[1] ", "") for x in sSplitLines]  
    #print "sSplitVal", sSplitVal     
       
    if detectionFile=="ClusterLimitSize":
        rawTh = [x for x in sSplitVal if "MINPS" in x]
        clusterTh = int(rawTh[0].split("=")[1])
        #clusterTh = int(rawTh[0].split(":")[1])
    else:
        clusterTh = "NA"
        pval_seuil = 1
        rawPVAL = [x for x in sSplitVal if "pval_seuil=" in x]
        pval_seuil = float(rawPVAL[0].split("=")[1])
        #print "sSplitVal", sSplitVal
    infoMsg = sSplitVal[-1].replace("\"", "")
    
    
    if not 'pval_seuil' in locals():
        pval_seuil = 1
    return clusterTh, infoMsg, pval_seuil

########################################################################################################################################################

def filterClusterSize(f, s):
	"""Removes from file f all SV with less or equal number of PS than threshold s """
	detection = open(f).readlines()
	fout = open(f, "w")
	for l in detection:
		if l.split(";")[0]=="manip":
			fout.write(l)
		elif l!='' and int(l.split(";")[3])>=int(s):
			fout.write(l)
	return 1

########################################################################################################################################################

def getChrList(bamF, silent=False):
	"""returns a dictionnary with chr length from a BAM file opened with pysam"""
	if silent:
		print "---------------------------------------------"
		print "Reading chromosome sizes from BAM header:"
 

	F = pysam.Samfile(bamF, "rb")
	chrStats=dict(zip(F.references, map(int, F.lengths)))
	genomelength = sum(chrStats.values())
	
	if silent:
		printplus(chrStats)
		print
		print "Genome length is:", '{0:,}'.format(genomelength)
		print
	return chrStats, genomelength

########################################################################################################################################################

def getAChr(t):
	"""return up to 10 chromosomal ranges to compute coverage"""
	chrCovDict = {k: t[k] for k in random.sample(t.keys(), min(len(t.keys()), 10))}

	for chromosome, startEnd in chrCovDict.iteritems():
		window = 10000
		if startEnd > window:
			st = random.randint(0, startEnd-window+1)
			ed = st + window
			chrCovDict[chromosome] = [st, ed]
		else:
			chrCovDict[chromosome] = [0, startEnd]
	print "Randomly choosen postions for coverage estimation are:"
	printplus(chrCovDict)
	print "---------------------------------------------"
	print

	return chrCovDict

########################################################################################################################################################
def printplus(obj):
    """
    Pretty-prints the object passed in.

    """
    # Dict
    if isinstance(obj, dict):
        for k, v in sorted(obj.items()):
            if isinstance(v, int):
            	v = '{0:,}'.format(v)
            print u'{0}: {1}'.format(k, v)

    # List or tuple
    elif isinstance(obj, list) or isinstance(obj, tuple):
        for x in obj:
            print x

    # Other
    else:
        print obj



########################################################################################################################################################

def updateCov(samFile, read, chrCovDict, covDir, compteur):
	"""Evaluate coverage"""
 
    
	chrom1 = samFile.getrname(read.tid)
	chrom2 = samFile.getrname(read.rnext)

	#print compteur
	if compteur%1000000 == 0 and compteur != 0:
		print "Processed", '{0:,}'.format(compteur), "reads"
	compteur+=1


	if chrom1 in chrCovDict.keys():

		startEv = chrCovDict[chrom1][0]
		endEv = chrCovDict[chrom1][1]

		if chrom1 == chrom2 and read.tlen > 0 and read.tlen < 100000:

			sPS, ePS = min([read.pos+1, read.pnext+1]), max([read.pos+1, read.pnext+1])
			#cov_part = float(getOverlapLen(sPS, ePS, startEv, endEv)) / (endEv - startEv)

			cov_part = getOverlapLen(sPS, ePS, startEv, endEv)
			#if cov_part >0:
			#	print sPS, ePS, startEv, endEv, cov_part
			covDir[chrom1] = covDir.get(chrom1, 0) + cov_part

	return compteur
     #return covDir, compteur

########################################################################################################################################################

def getOverlapLen(sPS, ePS, startEnv, endEnv):
	"""Get overlap between window (startEnv-endEnv) and a PS (sPS-ePS)"""
	if ePS>startEnv and ePS<endEnv and sPS<startEnv:
		return ePS-startEnv+1
	elif sPS<startEnv and ePS>endEnv:
		return endEnv-startEnv+1
	elif sPS<endEnv and sPS>startEnv and ePS>endEnv:
		return endEnv-sPS+1
	elif sPS>startEnv and sPS<endEnv and ePS>startEnv and ePS<endEnv:
		return ePS-sPS+1
	else:
		return 0

########################################################################################################################################################

def reduceReformat(reformat, FBam, paramOut):
    """ Makes a dict.table file for pval calculation (sort | uniq -c of the reformat.table.list)"""
#    fName = paramOut+".dist.table"
    fName = FBam+"_dist.table"
    sep = "_--_--_--_"
    try:
        with open(fName):
            pass
            msg = fName + " file already exists."
        return fName, msg
    except IOError:
    		# Build a dictionary of lines and their associated counts.
    		counts = {}
    		#input_file = open("/path/to/file", "r")
    		for PS in reformat:
    			counts[PS] = counts.get(PS, 0) + 1
    
    		# Build a list of [(lineA, countA), (lineB, countB), ... ]
    		sorted_counts = list(counts.items())
    
    		# Sort the list by (count, line) in reverse order.
    		sorted_counts.sort(lambda a,b: -cmp((a[1], a[0]), (b[1], b[0])))
    
    		# write the lines
    		with open(fName, "w") as f:
    			for line, count in sorted_counts:
    				#print line           
    				li = line.split(sep)
    				lw = [str(count), li[0], li[1], li[2], li[3]+"\n"]
    				f.write('\t'.join(lw))
    		msg = os.path.split(fName)[1] + " has been writen"
    return fName, msg


########################################################################################################################################################
def updateReformat(read, reformat, bam):
	""" Makes a reformat.table.list """
 	sep = "_--_--_--_"
	if len(reformat)<100000:
		if read.tlen > 0: #only keep intrachromosomal PS
			if read.is_reverse == 1:
				sign1 = "-"
			else:
				sign1 = "+"
			if read.mate_is_reverse == 1:
				sign2 = "-"
			else:
				sign2 = "+"
			chr1 = bam.getrname(read.tid)
			chr2 = bam.getrname(read.rnext)
			l = str(chr1)+sep+str(chr2)+sep+sign1+"...SIGN..."+sign2+sep+str(read.tlen) #tid not the chrName: pysam.Samfile.getrname(id)
			reformat.append(l)
	return reformat

########################################################################################################################################################

def makeCov(reads, startEv, endEv):
	"""Evaluate coverage"""
	val=0
	for read in reads:
		if read.is_read1 and read.is_paired and not read.mate_is_unmapped and \
		read.tlen < 100000 and read.tid==read.rnext:
			sPS, ePS = min([read.pos+1, read.pnext+1]), max([read.pos+1, read.pnext+1])
			cov_part = getOverlapLen(sPS, ePS, startEv, endEv)
			val+=cov_part
	cov = val / ((endEv-startEv)+1)

	return cov


########################################################################################################################################################

def cleanCov(covDir, chrCovDict):
	"""Estimate average coverage from raw coverage extracted from BAM """


#	chrVuList = covDir.keys()
	if min(covDir.values()) == 0:
		sys.exit("Error: Could not find any read on any of the windows from chromosomes choosen from BAM header")

	#print
	#print "Raw coverage:"
	#printplus(covDir)
	#print

	for chrom, raw_cov in covDir.iteritems():
		window_length = chrCovDict[chrom][1]-chrCovDict[chrom][0]
		covDir[chrom] = raw_cov / window_length

	print
	print "---------------------------------------------"
	print "Ended raw coverage aquisition:"
	printplus(covDir)
	print
	print "Cleaning raw coverage"

	cov_values = covDir.values()

	#print cov_values
	mean_cov  = sum(cov_values)/len(cov_values)
	std_dev = ((sum([x*x for x in cov_values]) / len(cov_values)) - (mean_cov * mean_cov))**.5

	for val in cov_values:
		if not abs(val - mean_cov) <= std_dev:
			cov_values.remove(val)

	print
	print "Filtered coverages are:", ", ".join([str(x) for x in cov_values])
	print

	average_coverage = sum(cov_values) / len(cov_values)

	print "The experimental coverage in fragment is:", int(average_coverage)
	print "---------------------------------------------"

	return int(average_coverage)

########################################################################################################################################################

def filterThisDetectionFileWithEstimatedClusterSize(detectionFile, ISmean, ISmed, ISmad, ISsd, l, n, seuil, scriptRClusterSize,functionsR, tipe, list_chr_length, list_chr_names):
    """Estimates the min cluster size and than filters the detection file"""

    s=runRForClusterThreshold(ISmean, ISmed, ISmad, ISsd, l, n, detectionFile, seuil, scriptRClusterSize,functionsR, tipe, list_chr_length, list_chr_names) #Get the PS/SV threshold
    if s != "Error":
        print "\nAt", round(float(seuil)*100,1), "% certitude, the estimated cluster minimal size is:", s, "PS\n"
    else:
        print "Error in threshold calculation"
        #filterClusterSize(detectionFile, s) #overwrites the detection file with the filtered one
    return 1

########################################################################################################################################################

def EstimatedClusterSize(detectionFile, ISmean, ISmed, ISmad, ISsd, l, n, seuil, scriptRClusterSize,functionsR, tipe, list_chr_length, list_chr_names):
	"""Estimates the min cluster size and than filters the detection file"""
	s=runRForCluster(ISmean, ISmed, ISmad, ISsd, l, n, detectionFile, seuil, scriptRClusterSize,functionsR, tipe, list_chr_length, list_chr_names) #Get the PS/SV threshold
	return s

########################################################################################################################################################

def filterDeletionsWithPQval(detectionFile, bam_input_name, reformat_dist, scriptRPval, average_cov, sval, functionsR, seuil, dirWork):
	"""Filter deletions detection file by FDR """

	dist_table, msg = reduceReformat(reformat_dist, bam_input_name, dirWork)
	print msg
	#dist_table = "/media/data/NA12878/ERR0012_3XX.noDUP.chr20.dist.table"
	#average_cov = 30
	print "\nFiltering deletions at", round(float(seuil)*100,1), "% certitude\n"


	call(["Rscript", "--slave --vanilla ",str(scriptRPval),"--args", str(average_cov), str(dist_table), str(detectionFile), str(sval), str(functionsR), str(seuil)])
	return 1


########################################################################################################################################################

def runRForCluster(ISmean, ISmed, ISmad, ISsd, l, n, cluster, seuil, script,functionsR, tipe, list_chr_length, list_chr_names):
	""" Fonction pour trouver le nombre de clusters minimum par SV. Utilise les parametres de la librairie (med, mad, sd), le nombre de PS discordantes, et le fichier de detection
	"""

	output=call("Rscript --silent "+str(script)+" --args "+str(ISmean)+" "+str(ISmed)+" "+str(ISmad)+" "+str(ISsd)+" "+str(n)+" "+str(l)+" "+str(cluster)+" "+str(seuil)+" "+str(functionsR)+" "+str(tipe)+" "+str(list_chr_length)+" "+str(list_chr_names), shell=True)

	if output=="":
		return "Error"
	else:
		return output


########################################################################################################################################################
def addCovToDelFile3(bam, delFile):
	""" adds local coverage to deletion detecton file """
	F = pysam.Samfile(bam, "rb")
	fDel = open(delFile)
	lines = fDel.readlines()


	if  "cov" in lines[0].rstrip():
			print "Coverage already present in "+delFile
			return 1
	header = lines[0].rstrip()+";\"cov\""
	lToWrite=[]
	lToWrite.append(header)

	deletionss = lines[1:]
	deletions = [x.rstrip() for x in deletionss]

	compteur = 0
	total = len(deletions)
	thForCov = 10000000

	for deletion in deletions:
		dele = deletion.rstrip().split(";")
		#TODO: REMOVE "chr" in chrom
		#chrom = "chr"+dele[1][1:-1]
		chrom = dele[1][1:-1]
		left = max(0,int(dele[6][1:-1])-500)
		right = max(0,int(dele[7][1:-1])+500)
		idDel = dele[2][1:-1]

		if abs(right-left)>thForCov:
			print "!!! Warning: Deletion number", idDel, "on chromosome", chrom, "is bigger \
			than", thForCov, "pb. Coverage will be estimated from the 10 first kb"
			avr=[]
			#print chrom, left, right
			for pileup in F.pileup(chrom, left, left+10000):
				avr.append(pileup.n)
				if len(avr) != 0:
					d_cov = (sum(avr)/len(avr))
					#print d_cov
				else:
					d_cov= 0

		else:

			avr=[]
			#print chrom, left, right
			for pileup in F.pileup(chrom, left, right):
				avr.append(pileup.n)
			if len(avr) != 0:
				d_cov = (sum(avr)/len(avr))
				#print d_cov
			else:
				d_cov= 0
		deletion = deletion + ";\""+ str(d_cov)+ "\""
		lToWrite.append(deletion)

		compteur +=1
		if compteur%200==0:
			print "Processed local coverage for deletion", compteur, "/", total

	fDel.close()
	delF = open(delFile, "w")
	for deletion in lToWrite:
		#print deletion   ########################################, "++++++++++"
		delF.write("%s\n" % deletion)

	print "---------------------------------------------"
	print "Added local coverage to deletion file"

	return 1
########################################################################################################################################################
def addCovToDelFile(bam, delFile):
	""" adds local coverage to deletion detecton file """
	print "\n---------------------------------------------"
	print "Starting coverage calculation"

	if not os.path.exists(bam+".bai"):
		pysam.index(bam)
		if not os.path.exists(bam+".bai"):
			print "\n\n**** Warning:", os.path.split(bam)[1], "is not sorted  ****"
			print "**** Sorting the file", os.path.split(bam)[1], "\n"
			pysam.sort(bam, bam+".sorted")
			pysam.index(bam+".sorted.bam")
			bam = bam+".sorted.bam"
	
 
 
 	F = pysam.Samfile(bam, "rb")
	fDel = open(delFile)
	lines = fDel.readlines()

	#print "addCov", delFile, lines[0].rstrip()
	if  "cov" in lines[0].rstrip():
			print "Coverage already present in "+delFile
			return 1
	header = lines[0].rstrip()+"\"cov\""
	lToWrite=[]
	lToWrite.append(header)

	deletionss = lines[1:]
	deletions = [x.rstrip() for x in deletionss]

	compteur = 0
	total = len(deletions)
	#thForCov = 10000000
	thForCov = 10000

	dicosCov = {}

	for deletion in deletions:
		dele = deletion.rstrip().split(";")
		#TODO: REMOVE "chr" in chrom
		#chrom = "chr"+dele[1][1:-1]
		chrom = dele[1][1:-1]
  		svID = dele[1][1:-1]+"-"+dele[2][1:-1]
		left = max(0,int(dele[6][1:-1]))
		right = max(0,int(dele[7][1:-1]))
		idDel = dele[2][1:-1]

		if abs(right-left)>thForCov:
			#print "!!! Warning: Deletion number", idDel, "on chromosome", chrom, "is bigger \
			#than", thForCov, "pb. Coverage will be estimated from the 10 first kb"
			avr=[]
			#print chrom, left, right
			for pileup in F.pileup(chrom, left, left+10000):
				avr.append(pileup.n)
				if len(avr) != 0:
					d_cov = (sum(avr)/len(avr))
					#print d_cov
				else:
					d_cov= 0

		else:

			avr=[]
			left, right = min(left, right), max(left, right)
			for pileup in F.pileup(chrom, left, right):
				avr.append(pileup.n)
			if len(avr) != 0:
				d_cov = (sum(avr)/len(avr))
				#print d_cov
			else:
				d_cov= 0
		dicosCov[svID] = d_cov
  
		deletion = deletion + ";\""+ str(d_cov)+ "\""
		lToWrite.append(deletion)

		compteur +=1
		if compteur%200==0:
			print "Processed local coverage for SV", compteur, "/", total

	fDel.close()
	delF = open(delFile, "w")
	for deletion in lToWrite:
		#print deletion   ########################################, "++++++++++"
		delF.write("%s\n" % deletion)

	#also add coverage to the byPS file 
	fDelbyPS = delFile.replace("bySV", "byRP")
	fDel = open(fDelbyPS)
	lines = fDel.readlines()
	header_tmp = lines[0].rstrip().split(";")
 	header_tmp.insert(17,"\"cov\"")
   	header = ";".join(header_tmp)
     	#print "header", header
 
 	PSToWrite=[]
	PSToWrite.append(header)
 
 	deletionss = lines[1:]
	deletions = [x.rstrip() for x in deletionss]
 
 	for PS in deletions:
		dele = PS.rstrip().split(";")
  		svID = dele[1][1:-1]+"-"+dele[2][1:-1]
		coco = "\"" + str(dicosCov[svID]) + "\""
		dele.insert(17, coco)
		#print "dele", dele
		PS = ";".join(dele)
		PSToWrite.append(PS)
	fDel.close()
	delF = open(fDelbyPS, "w")
	for PS in PSToWrite:
		
		delF.write("%s\n" % PS)
 
	print "Added local coverage to "+delFile
	print "---------------------------------------------"

	return 1
#----------------------------------------------------------------------------
def define_MinPS(tipe, stats, seuil_cluster, list_chr_length, list_chr_names, 
               detectionFile, paramsIn, reformat_dist, sval, seuil_fdr, 
               dirWork, RL, n_param, list_chr_real=['noGood666']):
    """ Calls runRForClusterThreshold to get the minimum number of PS required
    for SV"""

    #Set minimum nb of ps for detection of insertion
    #Attention nb_ps_min_ins = nb total de ps par jonction
    
    
        
    numberOfInterChromosomalPSFile = paramsIn+"_PScandidate_per_Xsome_pairs.csv"
    numberOfIntraChromosomalPSFile = paramsIn+"_PScandidate_intra.csv"
    dirprog = os.path.split(os.path.abspath(sys.argv[0]))[0]
    scriptRClusterSize = dirprog+"/scriptRClusterSize.R"
    functionsR = dirprog+"/Rfunctions.R"

    #average_cov is no longer needed. set to 0
    average_cov = 0
    seuil_fdr = float(seuil_fdr)
    
    if tipe == "INTER":
        ps_min_ins, msg, pval_seuil = runRForClusterThreshold("INS", stats["mean"], stats["median"], 
                                       stats["mad"], stats["stdev"], stats["genome_length"], stats["ninter"],
                                       numberOfInterChromosomalPSFile, seuil_cluster, 
                                       scriptRClusterSize, functionsR, 
                                       list_chr_length, list_chr_names, detectionFile,
                                       paramsIn, reformat_dist, average_cov, sval,
                                       seuil_fdr, dirWork, RL, n_param,
                                       numberOfIntraChromosomalPSFile)                                
        print msg
        
    
        ps_min_tr, msg, pval_seuil = runRForClusterThreshold("TR", stats["mean"], stats["median"], 
                                       stats["mad"], stats["stdev"], stats["genome_length"], stats["ninter"],
                                       numberOfInterChromosomalPSFile, seuil_cluster, 
                                       scriptRClusterSize, functionsR, 
                                       list_chr_length, list_chr_names, detectionFile,
                                       paramsIn, reformat_dist, average_cov, sval,
                                       seuil_fdr, dirWork, RL, n_param,
                                       numberOfIntraChromosomalPSFile)
    
        print msg
    
    
        ps_min_tn, msg, pval_seuil = runRForClusterThreshold("TN", stats["mean"], stats["median"], 
                                       stats["mad"], stats["stdev"], stats["genome_length"], stats["ninter"],
                                       numberOfInterChromosomalPSFile, seuil_cluster, 
                                       scriptRClusterSize, functionsR, 
                                       list_chr_length, list_chr_names, detectionFile,
                                       paramsIn, reformat_dist, average_cov, sval,
                                       seuil_fdr, dirWork, RL, n_param,
                                       numberOfIntraChromosomalPSFile)
                                       
        print msg
    
    
        return ps_min_ins, ps_min_tr, ps_min_tn, msg
    else:
        if tipe == "INV":
            n = U.getN(paramsIn+"_PScandidate_intra.csv", "INV", list_chr_real, stats["chromosome_prefix"])
        else:
            n = U.getN(paramsIn+"_PScandidate_intra.csv", "DUP", list_chr_real, stats["chromosome_prefix"])
       
	#print "Define MINPS for", tipe #, "-n"+tipe.lower(), stats["n"+tipe.lower()]
        ps_min, msg, pval_seuil = runRForClusterThreshold(tipe, stats["mean"], stats["median"], 
                               stats["mad"], stats["stdev"], stats["genome_length"], 
                               n, numberOfInterChromosomalPSFile, 
                            #stats["n"+tipe.lower()], numberOfInterChromosomalPSFile, 
                               seuil_cluster, scriptRClusterSize, functionsR, 
                               list_chr_length, list_chr_names, detectionFile,
                               paramsIn, reformat_dist, average_cov, sval,
                               seuil_fdr, dirWork, RL, n_param,
                               numberOfIntraChromosomalPSFile)
                               
	#print "$$$$$$$$$$$$$$$$$$$ ", ps_min, msg
        return ps_min, msg


#----------------------------------------------------------------------------
def functstats(tipe, stats, nDiscordant, seuil_cluster, list_chr_length, list_chr_names, 
               detectionFile, paramsIn, reformat_dist, sval, seuil_fdr, dirWork,
               RL, n_param):
    """ Calls runRForClusterThreshold to get the p-values for deletions and 
    scores for the other SV types"""

    #CALCUL DES P-VALEURS pour les deletions ou score pour les autres

    #Set minimum nb of ps for detection of insertion
    #Attention nb_ps_min_ins = nb total de ps par jonction

    numberOfInterChromosomalPSFile = paramsIn+"_PScandidate_per_Xsome_pairs.csv"
    numberOfIntraChromosomalPSFile = paramsIn+"_PScandidate_intra.csv"
    
    dirprog = os.path.split(os.path.abspath(sys.argv[0]))[0]
    scriptRClusterSize = dirprog+"/scriptRClusterSize.R"
    functionsR = dirprog+"/Rfunctions.R"    

    #average_cov is no longer needed. set to 0
    average_cov = 0
    
#    print tipe, stats["mean"], stats["median"], \
#                                   stats["mad"], stats["stdev"], stats["genome_length"], \
#                                   nDiscordant, numberOfInterChromosomalPSFile,\
#                                   seuil_cluster, scriptRClusterSize, functionsR, \
#                                   list_chr_length, list_chr_names, detectionFile,\
#                                   paramsIn, reformat_dist, average_cov, sval,\
#                                   seuil_fdr, dirWork, RL, n_param,\
#                                   numberOfIntraChromosomalPSFile
    
    ps_min, msg, pval_seuil = runRForClusterThreshold(tipe, stats["mean"], stats["median"], 
                                   stats["mad"], stats["stdev"], stats["genome_length"], 
                                   nDiscordant, numberOfInterChromosomalPSFile,
                                   seuil_cluster, scriptRClusterSize, functionsR, 
                                   list_chr_length, list_chr_names, detectionFile,
                                   paramsIn, reformat_dist, average_cov, sval,
                                   seuil_fdr, dirWork, RL, n_param,
                                   numberOfIntraChromosomalPSFile)
    print "Statistics done"
    return ps_min, msg, pval_seuil
