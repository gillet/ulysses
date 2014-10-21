#!/usr/bin/R
###
### Set of functions for pvalue processing on deletion events
###

options(warn=-1)
#invisible(suppressMessages(require(zoo)) || suppressMessages(install.packages("zoo")))
#suppressMessages(install.packages("BiocInstaller", repos="http://www.bioconductor.org/packages/2.11/bioc"))
suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
invisible(suppressWarnings(suppressMessages(require(IRanges) || if(!require(IRanges)) {source("http://www.bioconductor.org/biocLite.R"); try(biocLite("IRanges"), silent = T)})))
invisible(suppressWarnings(suppressMessages(require(qvalue) || if(!require(qvalue)) {source("http://www.bioconductor.org/biocLite.R"); try(biocLite("qvalue"), silent = T)})))
options(warn=0)




###
### Compute the empirical distribution of Insert size
###
GetISemp <- function(dist.tabulated, perc.clean = 0.999){
	#filt.table = subset(dist.tabulated, chr.F == chr.R & chr.F != "chrmt" & mate.orient == "+.-")
	#filt.table = subset(dist.tabulated, chr.F == chr.R, substr(mate.orient, 1,1) != substr(mate.orient, 3,3) & mate.dist > 400)
	filt.table = subset(dist.tabulated, chr.F == chr.R)
	dist.tab = as.data.frame(xtabs( Freq ~ mate.dist, data = filt.table))
	dist.tab$mate.dist = as.integer(as.character(dist.tab$mate.dist))
	dist.tab = dist.tab[order(dist.tab$mate.dist), ]
	dist.tab$cFreq = cumsum(dist.tab$Freq)
	d.bound = dist.tab$mate.dist[min(which(dist.tab$cFreq >= (perc.clean * tail(dist.tab$cFreq,n=1))))]
	IS.emp = subset(dist.tab, mate.dist <= d.bound)[,c("mate.dist", "Freq")]
	names(IS.emp) = c("IS", "Freq")
	##Distribution and cumulative.
	IS.emp$prob = IS.emp$Freq / sum(IS.emp$Freq)
	IS.emp$cprob = cumsum(IS.emp$prob)
	IS.emp$icprob = rev(cumsum(rev(IS.emp$prob)))
	return(IS.emp)
}

####################################################################################################################

###
### Get the mean coverage for each row in a table
### given a Rlelist of coverage
###
ComputeMeanCoverageOnTable <- function(table, covrle, c.chr = "chr", c.start = "start_min", c.end = "end_max", c.ID = "Del_id"){
	###Compute the mean coverage
  RD = RangedData(IRanges(start = table[,c.start], end = table[,c.end]), space = table[,c.chr])
  RD.list = ranges(RD)
  chr.loop = intersect(names(covrle), names(RD))
  mean.covs.list = lapply(chr.loop, function(x)  {
    data.frame(chr = x, start = start(RD.list[[x]]), end = end(RD.list[[x]] ),
               cov = viewMeans(Views(covrle[[x]], RD.list[[x]]))) } )
  mean.covs = do.call("rbind", mean.covs.list)
  t.out = merge(table, mean.covs, by.x = c(c.chr, c.start, c.end), by.y = c("chr", "start", "end"))
  return(t.out)
}

####################################################################################################################

### Compute the pvalue of getting at least n fragments with a insert size larger than
### s out of cov fragments
### the empirical distribution is given by IS.emp
Pvaldist = function(nobs, cov, s, IS.emp){
  cov.i = round(as.numeric(cov))
  probs = sapply( 1:length(nobs), function(x) {
    is = suppressWarnings(min(which(IS.emp$IS >= as.numeric(s[x]))))
    if(is== Inf){is = which(IS.emp$IS == max(IS.emp$IS))}
    qs = IS.emp$icprob[is]
    #pbinom(cov.i[x] -1 , cov.i[x], qs, lower.tail = FALSE)
    pbinom(nobs[x] -1 , cov.i[x], qs, lower.tail = FALSE)
    })
  #cat(probs)
  return(probs)
}



####################################################################################################################

###
### After pval computations, makes FDR,
### filters the files and write new files
###
filterAndSaveFiles <- function(detection, seuil, tipe, detectionFile) {
	## ADDS A BUG FOR R CRASH TO TEST R ERROR HANDLING IN PYTHON
  #add 2 extrem pvals (0 and 1) to the pval list (prevents qvalue to crash in extrem pval distribution cases)
  if(debug==TRUE) {print("filterAndSaveFiles1")}
  if(debug==TRUE) {print(Sys.time())}
  #print(head(detection))
  #print(detection$pval)
  
  #pvals <- as.numeric(c(detection$pval, 0,1)) 
  pp = as.numeric(as.character(detection$pval))
  pvals <- c(pp,1)
 
  #print(detection)
  #qvals <- qvalue(pvals) #FDR (qvalues)
  qvals <- qvalue(pvals, pi0.method="bootstrap") #FDR (qvalues)
  if(debug==TRUE) {print("filterAndSaveFiles1.1")}
  
  #check if FDR crashed or not
  if(is.numeric(qvals)==TRUE) {
    msg = paste("!!! Error No FDR performed for ", tipe, ". Either all ", tipe, " are significant or all are not significant !!!")
    message(msg)
    write.table(detection, file = detectionFile, dec=".", sep=";", row.names = FALSE, quote=FALSE) #write to remove " symbol from file
    detectionFile.byPS.name <- sub("_bySV.csv", "_byPS.csv", detectionFile) #byPS file
    detectionFile.byPS <- read.csv2(detectionFile.byPS.name, check.names = FALSE) #read byPS file
    write.table(detectionFile.byPS, file = detectionFile.byPS.name, dec=".", sep=";", row.names = FALSE, quote=FALSE) #write to remove " symbol from file
    stop()
  }
  if(debug==TRUE) {print("filterAndSaveFiles1.2")}
  ###Add all pvalues to original detection file
  #BY SV
  #sort nicely by chromosomes and SV ID
  detectionO <- changeName(detection, "chr") #column names with special characters are a problem for "order" => change tmp the names
  detectionO <- detectionO[with(detectionO, order(chr, ID)), ] #sort nicely
  detectionO <- restoreName(detectionO) #restores names
  write.table(detectionO, file = detectionFile, dec=".", sep=";", row.names = FALSE, quote=FALSE) #write in a new file
  
  if(debug==TRUE) {print("filterAndSaveFiles2")}
  if(debug==TRUE) {print(Sys.time())}
  
  #BY PS
  f.byPS <- sub("_bySV.csv", "_byPS.csv", detectionFile) #byPS file
  cur.byPS <- read.csv2(f.byPS, check.names = FALSE) #read byPS file
  cur.byPS[["pval"]] <- NULL #remove existing pval
  #the merge sometimes does not work because of column 2: "(pair of) chromosome(s)"
  #change temporarely its name
  cur.byPS <- changeName(cur.byPS, "pair")
  #add pval to byPS file
  #if(tipe=="DEL" || tipe=="sINS") { #DEL have the coverage column, not the others
    tmp <- data.frame(ID=detection[["ID"]], pair=detection[,2], pval=detection[["pval"]], cov=detection[["cov"]])
    cur.byPS.pval <- merge(cur.byPS, tmp, by=c("pair", "ID"))
    #if(debug==TRUE){print(head(cur.byPS.pval))}
    cur.byPS.pval <- cur.byPS.pval[,c(3,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,27,28,19,20,21,22,23,24,25,26)] #deletions have the coverage column
    nn = names(cur.byPS.pval)
    nn[18]="cov"
    names(cur.byPS.pval) = nn
#   } else {
#     
#     tmp <- data.frame(ID=detection[["ID"]], pair=detection[,2], pval=detection[["pval"]])
#     cur.byPS.pval <- merge(cur.byPS, tmp, by=c("pair", "ID"))
#     cur.byPS.pval <- cur.byPS.pval[,c(3,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,26,18,19,20,21,22,23,24,25)]
#     
#   }
 

  if(debug==TRUE) {print("filterAndSaveFiles2.1")}
  if(debug==TRUE) {print(Sys.time())}
  
  #put the correct name for column 2
  cur.byPS.pval <- restoreName(cur.byPS.pval)
  #write the filtered byPS with pval file
  write.table(cur.byPS.pval, file = f.byPS, dec=".", sep=";", row.names = FALSE, quote=FALSE)
  
  if(debug==TRUE) {print("filterAndSaveFiles3")}
  if(debug==TRUE) {print(Sys.time())}
  
  #get pvalue cut-off for FDR with threshold "seuil.del"
  tseuil <- suppressWarnings(max(qvals$pvalues[qvals$qvalues<seuil]))
  
  msg = paste("The p-value cut-off estimated to filter", tipe, "at a threshold of (", seuil*100, "%) is", tseuil)
  write(msg, stdout())  
  #message(msg)
  msg = paste("pval_seuil=",tseuil, sep="")
  write(msg, stdout())
  
 
  if(tseuil== -Inf){
    msg = paste("\n!!! Warning, No filtering performed for", tipe, ". All", tipe, "appear to be significant with FDR threshold", seuil, ". You can also try to increase th FDR threshold to increase specificity.\n")
    message(msg)
    write.table(detection, file = detectionFile, dec=".", sep=";", row.names = FALSE, quote=FALSE) #write to remove " symbol from file
    detectionFile.byPS.name <- sub("_bySV.csv", "_byPS.csv", detectionFile) #byPS file
    detectionFile.byPS <- read.csv2(detectionFile.byPS.name, check.names = FALSE) #read byPS file
    write.table(detectionFile.byPS, file = detectionFile.byPS.name, dec=".", sep=";", row.names = FALSE, quote=FALSE) #write to remove " symbol from file
    stop()
  }
  #print("toto2")
  #print(detection)
  #filter deletions
  #f.qval <- subset(detection, pval<=tseuil) #filter deletions
  #print(length(detection[,1]))
  #print(tseuil)
  #print(detection)
  f.qval <- subset(detection, pval<=tseuil)
  #f.qval <- subset(detection, pval<tseuil)
  #print(length(f.qval[,1]))
  #print("toto1")
  #Make some countings
  numBefore <- length(detection[,1])
  numAfter <- length(f.qval[,1])
  nfiltrees <- numBefore - numAfter #number of removed deletions
  
  if(debug==TRUE) {print("filterAndSaveFiles4")}
  if(debug==TRUE) {print(Sys.time())}
  
  #sort nicely by chromosomes and SV ID
  detectionB <- changeName(f.qval, "chr") #column names with special characters are a problem for "order" => change tmp the names
  detectionB <- detectionB[with(detectionB, order(chr, ID)), ] #sort nicely
  detectionB <- restoreName(detectionB) #restores names
  detectionFileOut <- sub("_bySV.csv", "_bySV.stats.csv", detectionFile) #build a bySV stat file name
  #write.csv2(f.qval, file = cur.f, row.names = FALSE, quote=FALSE) #used to arase existing file
  write.table(detectionB, file = detectionFileOut, dec=".", sep=";", row.names = FALSE, quote=FALSE) #write in a new file
  
  # also filter the byPS file and add pval and cov
  id_remaining <- paste(f.qval[,2], f.qval[["ID"]], sep=".") #ID of remaining deletions
  f.out.byPS <- sub("_bySV.csv", "_byPS.csv", detectionFile) #byPS file
  cur.byPS <- read.csv2(f.out.byPS, check.names = FALSE) #read byPS file
  cur.byPS[["pval"]] <- NULL #remove existing pval
  
  ##the merge sometimes does not work because of column 2: "(pair of) chromosome(s)"
  ##change temporarely its name
  cur.byPS <- changeName(cur.byPS, "pair")
  
  ##keep only PS that belong to non filtered deletions
  cur.byPS.filtered <- subset(cur.byPS, paste(pair, ID, sep=".") %in% id_remaining)
  
  if(debug==TRUE) {print("filterAndSaveFiles5")}
  if(debug==TRUE) {print(Sys.time())}
  
   ##add pval to byPS file
#  if(tipe=="DEL" || tipe=="sINS") {
    tmp <- data.frame(ID=f.qval[["ID"]], pair=f.qval[,2], pval=f.qval[["pval"]], cov=f.qval[["cov"]])
    suppressWarnings(cur.byPS.filtered.qval <- merge(cur.byPS.filtered, tmp, by=c("pair", "ID")))
    if(debug==TRUE){print(head(cur.byPS.filtered.qval))}
    cur.byPS.filtered.qval <- cur.byPS.filtered.qval[,c(3,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,29,19,28,20,21,22,23,24,25,26,27)] #deletions have the coverage column
    nn = names(cur.byPS.filtered.qval)
    nn[18]="cov"
    names(cur.byPS.filtered.qval) = nn
#   } else {
#     tmp <- data.frame(ID=f.qval[["ID"]], pair=f.qval[,2], pval=f.qval[["pval"]])
#     cur.byPS.filtered.qval <- merge(cur.byPS.filtered, tmp, by=c("pair", "ID"))
#     cur.byPS.filtered.qval <- cur.byPS.filtered.qval[,c(3,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,26,18,19,20,21,22,23,24,25)]
#   }
                                                      
  ##put the correct name for column 2
  cur.byPS.filtered.qval <- restoreName(cur.byPS.filtered.qval)
  
  ##write the filtered byPS with pval file
  detectionFileOutbbyPS <- sub("_byPS.csv", "_byPS.stats.csv", f.out.byPS) #bySV file
  write.table(cur.byPS.filtered.qval, file = detectionFileOutbbyPS, dec=".", sep=";", row.names = FALSE, quote=FALSE)
  
  #update detection file (both byPS and bySV) with pvalues
  #write.table(detection, file = detectionFile, dec=".", sep=";", row.names = FALSE, quote=FALSE)

  #write(paste("seuil: ", seuil, sep=""), stdout())
  msg <- paste("FDR analysis successfully performed on files ", detectionFile, " and ", f.out.byPS, ". ", numAfter, " SV remain (", nfiltrees, " SV removed)", sep="")
  write(msg, stdout())
}



####################################################################################################################

###
### This function maybe produces pvalues which are to conservative (do we really have the IS.emp
### distribution for the fragments observed on any bp ?)
###
ComputePvalDeletions <- function(tab.del, IS.emp, s.cols = 5:19, col.cov = "cov" ){
  tab.out = tab.del
  cols.comp = colnames(tab.del)[s.cols]
  #cat(cols.comp, "\n")
  for (ccol in cols.comp){
    #Get size from name
    cs = as.integer(sub("s", "", ccol)) * 100
    ccount = sub("-", "0", sapply(strsplit(tab.del[[ccol]], "/", fixed = TRUE), function(x) x[1]))
    ccount = as.integer(ccount)
    c.pval = Pvaldist( ccount, tab.del[[col.cov]], cs, IS.emp)
    cpname = paste("pval", ccol, sep = ".")
    tab.out[[cpname]] = c.pval
  }
  return(tab.out)
}

####################################################################################################################

###
### Parse the set of coverage files to a list of Rle (for an easy computation of coverage) 
###
GetCovAsRle <- function(d.covs){
	f.list = dir(path = paste("./",d.covs, sep = "") , pattern = ".csv", full.names = TRUE)
	chr.l = sub("histo_", "",sapply(strsplit(basename(f.list),".",fixed = TRUE), head, n = 1))
	Rle.cov.list = list()
	for(i in 1:length(f.list)){
		cchr = chr.l[i]
		Rle.cov.list[[cchr]] = Rle(scan(f.list[i], what = list(pos = 1, cov = 1), quiet = TRUE)$cov)
		}
	return(Rle.cov.list)
}

####################################################################################################################

###
### FDR from pval to get threshold for pval
###
getFalseDicoveryRate <- function(tab.del.pval, seuil){
	s.cols = grep("pval", colnames(tab.del.pval)) 	#Get the pval col names	and their positions
	cols.list = colnames(tab.del.pval)[s.cols]
	for(i in 1:length(cols.list)){
		#Pour chaque seuil,		
		id = cols.list[i]
		#calculer les qval de ce seuil
		pvals = as.numeric(c(tab.del.pval[["pval"]]),0,1)
		qvals = qvalue(pvals)
		seuil = max(qvals$pvalues[qvals$qvalues<=seuil])
	}
	return(seuil)
}




####################################################################################################################

###
### Adds nbObs inter fragment to ClusterAll as well as chr lengths
###
getClusterAllCleanForINTER <- function(xdiff, detection, d, l, t.n, detectionFile, tipe, cAll) {
    
  xdiff.n <- names(xdiff) #name the chr pair column "(pair of) chromosome(s)"
  xdiff.n[xdiff.n == "V1"] <- names(detection)[2]		
  names(xdiff) <- xdiff.n
    
  #in the detection file, pairs can be in whatever order, not in xdiff
  #build a xdiff2 df in order to be able to easely take inverted chr combinations into consideration
  xdiff2 <- data.frame("(pair of) chromosome(s)"=do.call(paste,c(xdiff[c("X2","X1")], sep="-")), PSdisc=xdiff[["V2"]], check.names=F) 
  xdiff1 <- data.frame("(pair of) chromosome(s)"=xdiff[["(pair of) chromosome(s)"]],PSdisc=xdiff[["V2"]], check.names=F)
  
  #add the number of discordant reads (in xdiff1 & 2) on the stats of detection (cAll)
  clusterAll1 <- merge(cAll, xdiff1, by = names(detection)[2])
  clusterAll2 <- merge(cAll, xdiff2, by = names(detection)[2])
  
  clusterAll <- rbind(clusterAll1,clusterAll2)
    
  #split pairs of chrom in 2 columns in the stats df in order to add chr lengths
  clusterAll <- cbind(data.frame(do.call(rbind,strsplit(as.character(clusterAll[["(pair of) chromosome(s)"]]), "-")), check.names=F), clusterAll) 
  
  #change the chr column names to X1 and X2
  neu <- names(clusterAll)
  neu[1:2] <- c("X1","X2")
  names(clusterAll) <- neu
  
  #add chr length
  clusterAll[["X1len"]] <- sapply(clusterAll[["X1"]], function(x) {getChrLen(x, t.n)}) #add chrA length
  clusterAll[["X2len"]] <- sapply(clusterAll[["X2"]], function(x) {getChrLen(x, t.n)}) #add chrB length
  
  clusterAll <- subset(clusterAll, X1len!=0 & X2len!=0) #remove wired things (depreciated)
  
  return(clusterAll)
}




####################################################################################################################

###
### Adds nbObs inter fragment to ClusterAll as well as chr lengths
###
getClusterAllDetailForINTRA <- function(xdiff, detection, d, l, t.n, detectionFile, tipe, cAll) {
    
  #rename the "(pair of) chromosome(s)" of cAll to be able to merge with xdiff
  names(cAll)=c("chr", "nbPS", "freq")
  
  #SV were not found in some chr, remove them from xdiff
  disp = cAll$chr
  xdiff.t = subset(xdiff, chr %in% disp, select=c("chr", tipe))
  names(xdiff.t)=c("chr", "PSdisc")
  
  #add the number of discordant reads (in xdiff1 & 2) on the stats of detection (cAll)
  clusterAll <- merge(cAll, xdiff.t, by = "chr")

  #add chr length
  clusterAll[["len"]] <- sapply(clusterAll[["chr"]], function(x) {getChrLen(x, t.n)}) #add chrA length

  #Sort nicely by chr
  clusterAll = clusterAll[order(clusterAll[,1]),]

  return(clusterAll)
}




####################################################################################################################

###
### Calculate a pseudo-Mates interactions
###
getISddist <- function(IS.mean, sd ){
  #options(digits=14)
  proba.tmp=round(rnorm(1e5, mean= IS.mean, sd), d=0)
  proba = proba.tmp[proba.tmp>=0]
  #ISdist = prop.table(table(proba.tmp)) #The distribution from the normal function
  ISdist = prop.table(table(proba)) #The distribution from the normal function
  ISddist = as.data.frame(ISdist)
  ISddist[,1] = as.integer(as.character(ISddist[,1]))
  names(ISddist) = c("dist", "proba")

  for(x in seq(min(ISddist$dist), max(ISddist$dist))){ #fill the table with missing values
    if(is.element(x, ISddist$dist)==FALSE){
      ISddist = rbind(ISddist, c(x, 0))
      #ISddist = insertRow(ISddist, c(x, 0))
    }
  }
  ISddist = ISddist[with(ISddist, order(dist)),]
  return(ISddist)
}



####################################################################################################################
###
### insert row in DF
###
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}


####################################################################################################################
###
### Calculate pIS from a distribution
###
getPIS <- function(ISddist, ISc, nPS){
  #compute pIS for DUP and INV
  #pIS <- sum(sapply(1:(nrow(ISddist)-ISc+1), function(i){ ISddist$proba[i]*(sum(ISddist$proba[i:(i+ISc-1)])^(nPS))} ), na.rm=T)
  pIS <- sum(sapply(1:(nrow(ISddist)-ISc+1), function(i){ ISddist$proba[i]*(sum(ISddist$proba[i:(i+ISc-1)])^(nPS))} ), na.rm=T)
  return(pIS)
}

### Calculate pIS from a distribution part 2
getPIS.2 <- function(ISddist, ISc, nPS){
  #compute pIS for DUP and INV
  nISdist = nrow(ISddist)
  pIS <- sum(sapply(1:nISdist, function(i){ ISddist$proba[i]*(sum(ISddist$proba[max(1,i-ISc+1):min(nISdist, i+ISc-1)])^(nPS-1))} ), na.rm=T)
  return(pIS)
}


####################################################################################################################
###
### Returns the total number os expected SV for interchromosomal SV.
### It is basicaly the sum of predicted SV for each pair of chromosomes
###
getC10kFor2ChrCombinations <- function(xdiff, tipe){
  #print(tipe)
  #
  if(tipe=="INS"){ 
    return(sum(apply(xdiff, 1, getPINS, tipe=tipe)))
  } else if(tipe=="TN") {
    return(sum(apply(xdiff, 1, getPTN, tipe=tipe)))
  } else {
    return(sum(apply(xdiff, 1, getPTR, tipe=tipe)))
  }
}


####################################################################################################################
###
### Get the minps threshold so that no more than "seuil" SV will be expected during detection
###
fixeMinps <- function(tipe, ISddist, ISc, l, d, n, seuil, xdiff){
  c10k <- seuil+1 #sets the initial number of expexted SV (c10k) to more than the threshold "seuil" 
  minps <- 2 #we first test the number of expected SV for clusters of size "minps"
  
  ###### DUP, DEL, INV, sINS (only 1 chr)
  if(tipe=="DUP" | tipe=="INV" | tipe=="DEL" | tipe=="sINS"){
    while(c10k > seuil & minps < 1000){  #calculate the expected number of clusters and stop when c10k < to seuil
      pIS <- getPIS.2(ISddist, ISc, minps) #get the proba of "minps" to have a compatible IS
      ###FIXME change pIS vector
      pcomp <- getPCOMP(l,d,tipe, pIS, minps,1) #proba for "minps" to be overlapping
      c10k <- choose(n,minps)*pcomp #number of expected overlapping clusters of size "minps"
      #print(paste(minps, c10k, n, sep=" ----" ))
      minps <- minps + 1  #number of expected overlapping clusters of size "minps"
    }
 
  } else {
    ###### INS, TN, TR (2 chr cobination)
    
    while(c10k > seuil & minps < 1000){ #calculate the expected number of clusters and stop when c10k < to seuil
      xdiff[["nbPS"]] <- minps #adds "minps" to the interchromosomal stats (xdiff)
      c10k <- getC10kFor2ChrCombinations(xdiff, tipe) #returns t
      minps <- minps + 1 #number of expected overlapping clusters of size "minps"
    }
  }
  if(minps>2){
    #write(minps-1, stdout())
    #write(c10k, stdout())
    #write(seuil, stdout())
    #msg <- paste("MINPS for ", tipe, " is : ", minps-1, sep="")
    msg <- paste("MINPS_", tipe, "=", minps-1, sep="")
    write(msg, stdout())
  } else {
    #write(2, stdout())
    #write(c10k, stdout())
    #write(seuil, stdout())
    #msg <- paste("MINPS for ", tipe, " is : ", 2, sep="")
    msg <- paste("MINPS_", tipe, "=",2, sep="")
    write(msg, stdout())
  }
  
}


fixeMinps.2 <- function(tipe, ISddist, ISc, l, d, n, seuil, xdiff){
  c10k <- seuil+1 #sets the initial number of expexted SV (c10k) to more than the threshold "seuil" 
  minps <- 2 #we first test the number of expected SV for clusters of size "minps"
  nPS = 10
  if(debug==TRUE) {print("INfixeMinPS.1")}
  if(debug==TRUE) {print(tipe)}
  
  ###### DUP, DEL, INV, sINS (only 1 chr)
  if(tipe=="DUP" | tipe=="INV" | tipe=="DEL" | tipe=="sINS"){
    if(debug==TRUE) {print("INfixeMinPS.1.1")}
    if(median(ISddist$dist)>1000) {
      if(debug==TRUE) {print("INfixeMinPS.1.1.1")}
      #print("start compute pIS fast")
      pISs <- ComputePIS.fast(ISddist, round(min(ISc, median(ISddist$dist))-20), nPSmax = nPS+20) #get all pIS, from 2PS to nPS (Warning, first pIS is for clusters of 2, not 1)
      #print("end start compute pIS fast")
    } else {
      if(debug==TRUE) {print("INfixeMinPS.1.2")}
      #print(median(ISddist$dist))
      #print(min(ISc,median(ISddist$dist))-10)
      pISs <- ComputePIS.fast(ISddist, round(min(ISc,median(ISddist$dist))-20) , nPSmax = nPS+20) #get all pIS, from 2PS to nPS (Warning, first pIS is for clusters of 2, not 1)
      if(debug==TRUE) {print("INfixeMinPS.1.2.1")}
    }
    i = 1 # make index
    #print(pISs)
    if(debug==TRUE) {print("INfixeMinPS.1.3")}
    while(c10k > seuil & minps < 1000){  #calculate the expected number of clusters and stop when c10k < to seuil
      #pIS <- getPIS.2(ISddist, ISc, minps) #get the proba of "minps" to have a compatible IS
      if(debug==TRUE) {print(c10k)}
      if(debug==TRUE) {print(seuil)}
      if(debug==TRUE) {print(minps)}
      pIS <- pISs[i]
#       if(debug==TRUE) {print("INfixeMinPS.1.5")}
      pcomp <- getPCOMP(l,d,tipe, pISs, minps, pISs[2]) #proba for "minps" to be overlapping
#       if(debug==TRUE) {print("INfixeMinPS.1.6")}
      c10k <- choose(n,minps)*pcomp #number of expected overlapping clusters of size "minps"
      #print(paste(minps, c10k, n, sep=" ----" ))
#       if(debug==TRUE) {print("INfixeMinPS.1.7")}
      minps <- minps + 1  #number of expected overlapping clusters of size "minps"
      i <- i + 1
    }
    if(debug==TRUE) {print("INfixeMinPS.1.8")} 
  } else {
    ###### INS, TN, TR (2 chr cobination)
    if(debug==TRUE) {print("INfixeMinPS.inter1")}
    
    while(c10k > seuil & minps < 1000){ #calculate the expected number of clusters and stop when c10k < to seuil
      xdiff[["nbPS"]] <- minps #adds "minps" to the interchromosomal stats (xdiff)
      #if(debug==TRUE) {print(head(xdiff, 3))}
      c10k <- getC10kFor2ChrCombinations(xdiff, tipe) #returns t
      if(debug==TRUE) {cat(paste(minps, "PS: ", c10k, " SV\n", sep=""))}
      
      minps <- minps + 1 #number of expected overlapping clusters of size "minps"
    }
    if(debug==TRUE) {print("INfixeMinPS.inter2")}
  }
  if(minps>2){
    #write(minps-1, stdout())
    #write(c10k, stdout())
    #write(seuil, stdout())
    #msg <- paste("MINPS for ", tipe, " is : ", minps-1, sep="")
    msg <- paste("MINPS_", tipe, "=", minps-1, sep="")
    write(msg, stdout())
    
    #write("msg tototototototo sofijseiofj ", stdout())
    
  } else {
    #write(2, stdout())
    #write(c10k, stdout())
    #write(seuil, stdout())
    #msg <- paste("MINPS for ", tipe, " is : ", 2, sep="")
    #write("msg t646546 6546sf 684ototototototo sofijseiofj ", stdout())
    msg <- paste("MINPS_", tipe, "=", 2, sep="")
    write(msg, stdout())
  }
}

####################################################################################################################
###
### Calculate all pcomp for INS, TR, TN 
###
getAllpCOMPINTER <- function(tipe, l, d, clusterAll, v=1) {
  if(v==1){ # for minps evaluation
    sapply(1:length(clusterAll[,1]), function(x) {
      minps <- clusterAll[x,4]
      l1 <-  clusterAll[x,9]
      l2 <- clusterAll[x,10]
      #pcompVal <- getPCOMPINTER.new(d,tipe, minps, l1, l2) #proba for "minps" to be overlapping
      pcompVal <- getPCOMPINTER.new.hugues(d,tipe, minps, l1, l2) #proba for "minps" to be overlapping
      #print(pcompVal)
    })
  } else { #for stats, clusterAll is already a bit modified
    sapply(1:length(clusterAll[,1]), function(x) {
      minps <- clusterAll[x,4]
      l1 <-  clusterAll[x,7]
      l2 <- clusterAll[x,8]
      #pcompVal <- getPCOMPINTER.new(d,tipe, minps, l1, l2) #proba for "minps" to be overlapping
      pcompVal <- getPCOMPINTER.new.hugues(d,tipe, minps, l1, l2) #proba for "minps" to be overlapping
      #print(pcompVal)
    })
  }
}


####################################################################################################################
###
### Calculate all pcomp for INS, TR, TN 
###
getAllpCOMPINTRA <- function(tipe, l, d,ISddist, ISc, IS.sd, IS.mad, clusterAll, v=1) { 
  if(v==1){ # for minps evaluation
    #nPS <- max(clusterAll[,2])
    #pISs <- ComputePIS.fast(ISddist, ISc, nPSmax = nPS)
    sapply(1:length(clusterAll[,1]), function(x) {
      minps <- clusterAll[x,2]
      l1 <-  clusterAll[x,5]
      pIS <- getPIS.2(ISddist, ISc, minps)
      pcompVal <- getPCOMP(l1,d,tipe, pIS, minps, 1) #proba for "minps" to be overlapping
      #print(pcompVal)
    })
  } else { #for stats, clusterAll is already a bit modified
    #nPS <- max(clusterAll[,2])
    #pISs <- ComputePIS.fast(ISddist, ISc, nPSmax = nPS)
    sapply(1:length(clusterAll[,1]), function(x) {
      minps <- clusterAll[x,2]
      pIS <- getPIS.2(ISddist, ISc, minps)
      l1 <-  clusterAll[x,5]
      #proba for "minps" to be overlapping
      pcompVal <- getPCOMP(l1,d,tipe, pIS, minps, 1)
      #print(pcompVal)
    })
  }
}

getAllpCOMPINTRA.2 <- function(tipe, l, d,ISddist, ISc, IS.sd, IS.mad, clusterAll, v=1) {
  #if(debug==TRUE) {print(paste(mean(ISddist$dist), ISc, sep = "]---["))}
  #if(debug==TRUE) {print("getAllpCOMPINTRA.2.1")}
  
  if(v==1){ # for minps evaluation
    nPS <- max(clusterAll[,2])
    
    if(median(ISddist$dist)>1000) {
      ###TODO: change scale of ISddist
      pISs <- ComputePIS.fast(ISddist, ISc/10, nPSmax = nPS)
      #pISs <- ComputePIS.fast(ISddist, ISc/10, nPSmax = nPS)
    } else {
      pISs <- ComputePIS.fast(ISddist, ISc, nPSmax = nPS)
    }

    sapply(1:length(clusterAll[,1]), function(x) {
      minps <- clusterAll[x,2]
      l1 <-  clusterAll[x,5]
      pIS <- pISs[minps]
      pcompVal <- getPCOMP(l1,d,tipe, pISs, minps, pISs[2]) #proba for "minps" to be overlapping
      #pcompVal <- getPCOMP.hugues(l1,d,tipe, pISs, minps, pISs[2])
      #print(pcompVal)
    })
  } else { #for stats, clusterAll is already a bit modified
    nPS <- max(clusterAll[,2])
    #if(debug==TRUE) {print("getAllpCOMPINTRA.2.2")}
    #computePIS
    if(median(ISddist$dist)>1000) {
      pISs <- ComputePIS.fast(ISddist, ISc/10, nPSmax = nPS)
    } else {
      pISs <- ComputePIS.fast(ISddist, ISc, nPSmax = nPS)
    }
    #if(debug==TRUE) {print("getAllpCOMPINTRA.2.3")}
    
    #pISs <- ComputePIS.fast(ISddist, ISc, nPSmax = nPS)
    sapply(1:length(clusterAll[,1]), function(x) {
      #if(debug==TRUE) {print("getAllpCOMPINTRA.2.4")}
      minps <- clusterAll[x,2]
      #if(debug==TRUE) {print("getAllpCOMPINTRA.2.5")}
      #pIS <- getPIS.2(ISddist, ISc, minps)
      #pIS <- pISs[minps-1]
      l1 <-  clusterAll[x,5]
      #if(debug==TRUE) {print("getAllpCOMPINTRA.2.6")}
      #proba for "minps" to be overlapping
      #pcompVal <- getPCOMP(l1,d,tipe, pISs, minps, pISs[2])
      pcompVal <- getPCOMP.hugues(l1,d,tipe, pISs, minps, pISs[2]) #function(l,d, tipe, pISvec, n = 2, pIDde2)
      #if(debug==TRUE) {print("getAllpCOMPINTRA.2.7")}
      #print(pcompVal)
    })
  }
}

####################################################################################################################
###
### Calculate pcomp for INS, TR, TN (2-chr only)
###
getPCOMPINTER.OLD <- function(d, tipe, minps, l1, l2){
  if(tipe=="INS"){
    pcomp = (((2*d/l1)^(minps-1))*((2*d/l2)^(minps-2)) + ((2*d/l1)^(minps-2))*((2*d/l2)^(minps-1)))/2
    #pcomp = ((2*d/l1)^(2*(minps-1)))*((2*d/l2)^(minps-1)) + ((2*d/l2)^(2*(minps-1)))*((2*d/l1)^(minps-1)) 
  } else if(tipe=="TR") {
    #pcomp = (d/l1)^(minps-1)*(d/l2)^(minps-1)
    pcomp = ((2*d/l1)^(2*(minps-1)))*(2*d/l2)^(2*(minps-1))
  } else if(tipe=="TN") {
    #pcomp = (d/l1)^(minps-1)*(d/l2)^(minps-1)
    pcomp = (d*d/(l1*l2))^(minps-1)
  } 
  return(pcomp) 
}

###New version with the value of minPS fixed
getPCOMPINTER.new <- function(d, tipe, minps, l1, l2){
  ###Compute the proba that minps-1 are within a range d
  ###of the first one
  if(tipe=="INS"){
    pcomp = ( (2*d/l1)^(2*(minps-1)) ) * ( (2*d/l2)^(minps-1) )
  } else if(tipe=="TR") {
    pcomp = ( (2*d/l1)^(2*(minps-1)) )*(2*d/l2)^(2*(minps-1))
  } else if(tipe=="TN") {
    pcomp = (d*d/(l1*l2))^(minps-1)
  } 
  return(pcomp) 
}

####################################################################################################################
###
### Calculate pcomp for DEL, DUP, INV (single chr only)
###
getPCOMP <- function(l,d, tipe, pISvec, n=2, pISde2){
  p=1/l
  if(tipe=="DUP"| tipe=="DEL" | tipe=="sINS"){
    #pcomp = (p * d * pIS)^(n-1)  # Proba of having a cluster of size f 
    pcomp = ((p * d )^(n-1)) * pISvec[n]  # Proba of having a cluster of size f 
    #} else if(tipe=="INV") {
  } else {
    pcomp = (p*2*d*pISde2)*(((p*d)^(n-2))*pISvec[n-1]) 
  } 
  return(pcomp)
}

####Hugues version for probas functions
#### on distances between PSs and proba of an SV pIS 

### Compute the probability to have x Ps within a distance dn in a genome of length L
### for x = nPS
ComputePdistPS <- function(dn, L, nPS){
  mu = dn / L
  if (nPS == 1){
    return(1)
  }
  else{
    pDN =  mu^(nPS - 1) * ( nPS - (nPS-1)*mu )
    return(pDN)
  }
}

###Simulation of groups of x PS according to IS dist. and empirical proba to get them at 
###
pIS_simu  <- function(IS, IScn, nPS, nsim = 1000){
  mean(sapply(1:nsim, function(x) {diff(range(sample(IS$d, nPS, prob=IS$p))) <= IScn} ))
}

###Simulation for having nPS within a distance dn in a genome of length L
pdisPS_simu <- function(dn, L, nPS , nsim = 10000 ){
  mu = dn / L
  mean(sapply(1:nsim, function(x) {diff(range(runif(nPS))) <= mu } ))
}

pdisPS_simu.fast <- function(dn, L, nPS, nsim = 10000){
  pest = 0
  .C("pdist_simu", 
     as.integer(dn),
     as.integer(L),
     as.integer(nPS),
     as.double(pest),
     as.integer(nsim)
  )[[4]]
}


###New version, with correction on the computation of consistent distances dn
###!!!!WARNING WE ASSUME a vector pIS here
getPCOMP.hugues <- function(l,d, tipe, pISvec, n = 2, pISde2){
  p=1/l
  if(tipe=="DUP"| tipe=="DEL" | tipe=="sINS"){
    pcomp = ComputePdistPS(d, l , n) * pISvec[n]
  } else {
    pcomp = ComputePdistPS(2*d,l, 2) * pISde2 * ComputePdistPS(d, l, n-1) * pISvec[n-1]
  } 
  return(pcomp)  
}

###New version with the value of minPS fixed
getPCOMPINTER.new.hugues <- function(d, tipe, minps, l1, l2){
  ###Compute the proba that minps-1 are within a range d
  ###of the first one
  #d=d/2
  if(tipe=="INS"){
    #pcomp = (ComputePdistPS(2*d,l1, 2*minps) * ComputePdistPS(d, l2, minps)) #+(ComputePdistPS(d,l1, 2*minps) * ComputePdistPS(2*d, l2, minps)))/2
    #pcomp = ((ComputePdistPS(2*d,l1, 2*minps) * ComputePdistPS(2*d, l2, minps)) + (ComputePdistPS(2*d,l1, 2*minps) * ComputePdistPS(2*d, l2, minps)))/2
    pcomp = ((ComputePdistPS(2*d,l1, 2*(minps-2)) * ComputePdistPS(2*d, l2, minps-2)) + (ComputePdistPS(2*d,l1, minps-2) * ComputePdistPS(2*d, l2, 2*(minps-2))))/2
    #pcomp = ComputePdistPS(2*d,l1, 2*minps) * ComputePdistPS(2*d, l2, minps)
  } else if(tipe=="TR") {
    pcomp = ComputePdistPS(2*d,l1, 2*minps) * ComputePdistPS(2*d, l2, 2*minps)
  } else if(tipe=="TN") {
    pcomp = ComputePdistPS(d,l1, minps) *  ComputePdistPS(d,l2, minps)
    
  } 
  return(pcomp) 
}



####################################################################################################################
###
### From the ClusterAll file (frequency by nbPS classes from detection file), get all pcomp values
###
getAllpCOMP <- function(tipe, l, d, ISddist, ISc, IS.sd, IS.mad, clusterAll) {
  dfpcomp=data.frame(nbPS=0, pcomp=0)
  #nPS <- max(clusterAll[,2])
  #pISs <- ComputePIS.fast(ISddist, ISc, nPSmax = nPS)
  
  for(i in clusterAll$nbPS){
    minps <- as.numeric(i)
    pIS <- getPIS.2(ISddist, ISc, minps) #get the proba of "minps" to have a compatible IS
    #print(pIS)
    pcompVal <- getPCOMP(l,d,tipe, pIS, minps) #proba for "minps" to be overlapping
    tmp = data.frame(nbPS=minps, pcomp=pcompVal)
    dfpcomp=rbind(dfpcomp, tmp)
  }
  #dfpcomp <- dfpcomp[-1,]
  #print("---")
  #print(dfpcomp)
  #print("---")
  return(dfpcomp)
}

getAllpCOMP.2 <- function(tipe, l, d, ISddist, ISc, IS.sd, IS.mad, clusterAll) {
  dfpcomp=data.frame(nbPS=0, pcomp=0)
  
  nPS <- max(clusterAll$nbPS)
  pISs <- ComputePIS.fast(ISddist, ISc, nPSmax = nPS)
  
  for(i in clusterAll$nbPS){
    minps <- as.numeric(i)
    #pIS <- getPIS.2(ISddist, ISc, minps) #get the proba of "minps" to have a compatible IS
    #print(pIS)
    pIS <- pISs[minps-1]
    pcompVal <- getPCOMP(l,d,tipe, pIS, minps) #proba for "minps" to be overlapping
    tmp = data.frame(nbPS=minps, pcomp=pcompVal)
    dfpcomp=rbind(dfpcomp, tmp)
  }
  #dfpcomp <- dfpcomp[-1,]
  #print("---")
  #print(dfpcomp)
  #print("---")
  return(dfpcomp)
}
####################################################################################################################

###
### Calculate the expected number of clusters
###
getNbClusters <- function(n,d,l, ISc,f, pis, tipe, unOuDeuxD=1){ 
  #f = cluster size of overlapping PS
  #pis = probability of picking 2PS with compatible IS
  p=1/l  # proba of each position in the genome
  #if(substring(tipe,1,1)=="D" || substring(tipe,1,1)=="I"){}
  
  #pj = sum(head(pi,length(rollapply(pi, ISc, sum, by=1)))*pis^(f-1)) #proba of having f compatible PS (considering only size not position)	
  pDUP = (p * d * pis)^(f-1)	# Proba of having a cluster of size f 
  pINV = (p*2*d*pis)*((p*d*pis)^(f-2)) #
  #pINV = (p*2*d*pj)^(f-1) #
  
  if(tipe=="DUP"| tipe=="DEL" | tipe=="sINS"){
  ##For the distribution on the number of duplicates, go for a poisson of parameter pdup (when n is not too high)
	nbclusters = choose(n,f) * pDUP # 
  	return(round(nbclusters))
  } else {
  	nbclusters = choose(n,f) * pINV # expected number of cluster of size f, as many -- as ++
  	return(round(nbclusters))
  }
}

####################################################################################################################

###
### Calculate cluster size stats in the detection file
###
getClusterStats <- function(cl, tipe){
  t = read.csv2(cl, h=T, check.names = FALSE)
  #print("tata")
#   if(tipe=="DUP" | tipe=="INV"){
# 	  u = as.data.frame(table(t$nbPS))
# 	  #print(u)
# 	  names(u) = c("nbPS", "freq")
# 	  #u$nbPS = as.numeric(u$nbPS)
#   } else {
  #print(head(t))
  #yyA <- data.frame( "(pair of) chromosome(s)" = t[[ "(pair of) chromosome(s)"]], nbPS=t$nbA, check.names=F )
  #print(yyA)
  #yyB <- data.frame( "(pair of) chromosome(s)" = t[[ "(pair of) chromosome(s)"]], nbPS=t$nbB, check.names=F )
  #yy <- rbind(yyA, yyB)
  #u = aggregate(rep(1, nrow(yy)), by = list(x = yy[[ "(pair of) chromosome(s)"]], nbPS = yy$nbPS), sum)
  
  yy <- data.frame( "(pair of) chromosome(s)" = t[[ "(pair of) chromosome(s)"]], nbPS=t$nbPS, check.names=F )
  u = aggregate(rep(1, nrow(yy)), by = list(x = yy[[ "(pair of) chromosome(s)"]], nbPS = yy$nbPS), sum)
  names(u) = c("(pair of) chromosome(s)", "nbPS", "freq")
#  }
  u$freq = as.numeric(u$freq)
  return(u)
}

####################################################################################################################

###
### Calculate cluster size stats in the detection file
###
getClusterStatsINTER <- function(cl, tipe){
  t = read.csv2(cl, h=T, check.names = FALSE)
  #print(head(t))
  t$nbAB = paste(t$nbA, t$nbB, sep="-")
  
  yy <- data.frame( "(pair of) chromosome(s)" = t[[ "(pair of) chromosome(s)"]], nbPS=t$nbPS, nAB=t$nbAB, check.names=F )
  u = aggregate(rep(1, nrow(yy)), by = list(x = yy[[ "(pair of) chromosome(s)"]], nAB = yy$nAB), sum)
  #print(head(u))
  names(u) = c("(pair of) chromosome(s)", "nAB", "freq")
  #  }
  u$freq = as.numeric(u$freq)
  return(u)
}
##########################################################################################
#acces à une taille:
getChrLen <- function(num, allChr){
	cu.chr= as.vector(num)
	#if(('chr' %in% cu.chr) == FALSE){
	#  cu.chr = paste('chr', cu.chr, sep="")
	#}
  return(as.numeric(as.vector(allChr$tailles[allChr$chr==cu.chr])))
  }

##########################################################################################
#calculate the expected number of clusters for a combination
getPINS <- function(x, tipe){
	#print(x)
	ndis=as.numeric(x['V2']) #number of discordant fragments
	li=as.numeric(x['X1len'])
	lj=as.numeric(x['X2len'])
	d=as.numeric(x['d'])
	l=as.numeric(x['l'])
	PS=as.numeric(x['nbPS'])
	
  #print(paste(ndis, li, lj, d, l, PS, sep=" -- "))
  
  #pc1=((2*d/li)^(PS-1))*(2*d/lj)^(PS-2)
  #pc2=((2*d/lj)^(PS-1))*(2*d/li)^(PS-2)
  #s=choose(ndis,PS)*(pc1+pc2)/2 #only half of the cases give INS, the other half is TR
  
  #pc1 = ((2*d/li)^((PS-1)))*(2*d/lj)^(PS-2)
  #pc2 = ((2*d/lj)^((PS-1)))*(2*d/li)^(PS-2)
  #s=choose(ndis,PS)*(pc1+pc2)/2
  
	pIS = getPCOMPINTER.new.hugues(d, tipe, PS, li, lj)
	s=choose(ndis,PS)*pIS
  
  
		
	return(round(s))
}

###########################################################################################
#calculate the expected number of clusters for a combination
getPTR <- function(x, tipe){
	ndis=as.numeric(x['V2'])
	li=as.numeric(x['X1len'])
	lj=as.numeric(x['X2len'])
	d=as.numeric(x['d'])
	l=as.numeric(x['l'])
	PS=as.numeric(x['nbPS'])
  
	#pc1=((d/li)^(PS-1))*(d/lj)^(PS-1)
	#pc2=((d/lj)^(PS-1))*(d/li)^(PS-1)
	#s=choose(ndis,PS)*(pc1+pc2)/4 #only half of the cases give TR, the other half is INS
  
	#pc1 = ((2*d/li)^(2*(PS-1)))*((2*d/lj)^(2*(PS-1)))
	#s=choose(ndis,PS)*pc1 #only half of the cases give TR, the other half is INS
  
	pIS = getPCOMPINTER.new.hugues(d, tipe, PS, li, lj)
	s=choose(ndis,PS)*pIS
  
	
  
	return(round(s))
}

##########################################################################################

#calculate the expected number of clusters for a combination
getPTN <- function(x, tipe){
	ndis=as.numeric(x['V2'])
	li=as.numeric(x['X1len'])
	lj=as.numeric(x['X2len'])
	d=as.numeric(x['d'])
	l=as.numeric(x['l'])
	PS=as.numeric(x['nbPS'])
  
# 	pc1=((d/li)^(PS-1))*(d/lj)^(PS-1)
# 	s=choose(ndis,PS)*pc1/2 #only half of genome concerned
#   
# 	pc1=((d*d)/(li*lj))^(PS-1)
# 	s=choose(ndis,PS)*pc1

  pIS = getPCOMPINTER.new.hugues(d, tipe, PS, li, lj)
  s=choose(ndis,PS)*pIS
	  
	return(round(s))
}

##########################################################################################

#change name of column 2 "(pair of) chromosome(s)" to "chr"
changeName <- function(x, name) {
  tmp <- names(x)
  tmp[2] <- name
  names(x) <- tmp
  return(x)
}


##########################################################################################
#restore name of column 2 to "(pair of) chromosome(s)" 
restoreName <- function(x) {
  tmp <- names(x)
  tmp[2] <- "(pair of) chromosome(s)" 
  names(x) <- tmp
  return(x)
}




## Probability that x PS are all within a range 
### Dynamic programming approach
###Initialise the probability of 2 PS at distance IScn
###Just an upper triangular matrix by construction
ConsistentPS.2 = function(ISdist, IScn){
  ndists = diff(range(ISdist$dist)) + 1
  pmat = matrix(0, nrow=ndists, ncol = ndists)
  for (i in 1:(ndists-1)){
    pmat[i,i] = ISdist$proba[i] * ISdist$proba[i]
    for (j in (i + 1):min(ndists, i+IScn-1)){
      pmat[i,j] = 2 *ISdist$proba[i] * ISdist$proba[j]
    }
  }
  return( pmat )
}


##From one number of PS to the next one
ConsistentPS.next = function(pmat, IS, IScn){
  ndists = diff(range(IS$dist)) 
  pmat_next =  matrix(0, nrow=ndists, ncol = ndists)
  for (i in 1:(ndists-1))
  {
    for (j in (i+1):min(ndists, i+IScn-1))
    {
      ##overcount some extension already in the inclusion case
      pmat_next[i,j] = IS$proba[i] * sum(pmat[i:j, j ]) +
        IS$proba[j] * sum(pmat[i, i:j ] ) +
        pmat[i,j] * sum(IS$proba[ (i+1):(j-1) ] )
    }
  }
  return(pmat_next)
}

ConsistentPS.next.fast = function(pmat, IS, IScn){
  ndists = diff(range(IS$dist)) + 1
  return( sapply(1:(ndists-1), function(i) {
    sapply((i+1):min(ndists, i+IScn-1), function(j){
      IS$proba[i] * sum(pmat[i:j, j ]) + 
        IS$proba[j] * sum(pmat[i, i:j ] ) +
        pmat[i,j] * sum(IS$proba[ (i+1):(j-1) ] ) })
  }))
}

ComputePIS.fast = function(IS, IScn, nPSmax = 2){
  ndist = diff(range(IS$dist)) 
  .C("nPsdist",
     as.double(as.vector(IS$proba)),
     as.double(vector("double", nPSmax)),
     as.integer(ndist),
     as.integer(nPSmax),
     as.integer(IScn))[[2]]
}


###Compute the probabilities of consistent PS up to nPS
ComputePIS = function(IS, IScn, nPSmax=2){
  pIS = rep(1, nPSmax)
  #cat("Dynamic Programming for pIS up to n =", nPSmax, "\n")
  pmat = ConsistentPS.2(IS, IScn)
  pIS[2] = sum(pmat)
  i = 3
  while (i <= nPSmax){
    #cat("#")
    pmat_next = ConsistentPS.next(pmat, IS, IScn)
    pIS[i] = sum(pmat_next)
    pmat = pmat_next
    i = i+1
  }
  #cat("\n")
  return(pIS)
}

###Simulation of groups of x PS according to IS dist. and empirical proba to get them at 
###
pIS_simu  <- function(IS, IScn, nPS, nsim = 1000){
  mean(sapply(1:nsim, function(x) {diff(range(sample(IS$d, nPS, prob=IS$p))) <= IScn} ))
}



### a function to compare results up to nPS PS 
### starting from ISdist, sigma the dispersion of the IS, 
### n the parameter of IScn, and a max of nPS PS
Compare_pIS_approx <- function(ISdist, sigma, n, plot.template = "Evaluate_PIS", nPS = 10, nreplicate = 30, nsim = 1000){
  IScn = n*sigma
  pIS.alex.corr = sapply(1:nPS, function(x){ getPIS.2(ISdist, IScn, nPS= x)})
  PIS.DP = ComputePIS.fast(ISdist, IScn, nPSmax = nPS)
  #cat("Simulations\n")
  #PIS.simu = sapply(1:nreplicate, function(x){sapply(1:nPS, function(x) pIS_simu(ISdist, IScn, x, nsim = nsim)) })
  #m.simu = melt(t(PIS.simu))
  #colnames(m.simu) = c("rep", "nPS", "pIS")
  #m.simu$type = "Simulation"
  m.simu = pIS.simurep(ISdist, IScn, nPS= nPS, nreplicate= nreplicate)
  #cat("Merging results\n")
  m.out = m.simu[, c("nPS", "pIS", "type")]
  m.out = rbind(m.out, data.frame(nPS = 1:nPS, pIS = pIS.alex.corr, type = "Approx"))
  m.out = rbind(m.out, data.frame(nPS = 1:nPS, pIS = PIS.DP, type = "Dyn.Prog."))
  plotout = paste(plot.template, ".sigma_", sigma, ".n_", n, ".pdf" , sep = "")
  pdf(plotout)
  boxplot(pIS ~ nPS, data = m.simu, ylim = c(0,1))
  points(1:nPS, pIS.alex.corr, col = "red", pch = 19)
  points(1:nPS, PIS.DP, col = "blue", pch = 19)
  legend("bottomleft", c("Original", "Approx", "Dyn. Prog."), col = c("darkred", "red","blue"), pch = 19)
  dev.off()
  return (m.out)
}

pIS.simurep = function(ISdist, IScn, nPS, nsim= 1000, nreplicate = 20){
  PIS.simu = sapply(1:nreplicate, function(x){sapply(1:nPS, function(x) pIS_simu(ISdist, IScn, x, nsim = nsim)) })
  m.simu = melt(t(PIS.simu))
  colnames(m.simu) = c("rep", "nPS", "pIS")
  m.simu$type = "Simulation"
  return(m.simu)
}

###Basic function for looking at the effect on
Compare_Approx_gaussian <- function(mu, sigma, n, nPS = 10){
  eps_keep = 1e-6
  ISmin = round(max(1, qnorm(eps_keep, mean = mu, sd = sigma, lower.tail = TRUE)))
  ISmax = round(qnorm(eps_keep, mean = mu, sd = sigma, lower.tail = FALSE))
  dvals = seq(ISmin, ISmax)
  ISdist = data.frame( dist = dvals, proba = dnorm(dvals, mean =  mu, sd = sigma) )
  ISdist$proba[1] = ISdist$proba[1] + 1 - sum(ISdist$proba)
  pl.t = paste("Gaussian.mu_", mu, sep = "")
  a = Compare_pIS_approx(ISdist, sigma, n, plot.template = pl.t, nPS = nPS)
  return(a)
}

