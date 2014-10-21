# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 13:45:41 2012

@author: ingridl
"""

import os, sys
from operator import itemgetter
try:
    import pysam
except ImportError:
    print "\tError: Pysam is not installed. \n\
    Please see: https://github.com/pysam-developers/pysam\n\n"
    sys.exit()
import Ulysse_utils as U
import Ulysse_stats as Ualex
import datetime
now = datetime.datetime.now()


# on a besoin de , meme si on a d, car meme si les conditions sont respectees
# a droite et a gauche pour le d, on peut se retrouver avec des IS qui ont des
# tailles tres differentes et qui ne peuvent pas "coller" avec la meme SV

#--------------------------------------------------------------------------
def ReadFilesBAM(params, cx, subtel, median, ori, dicQual):
    """
    Read BAM Files of discordant PS 
    """

    clasamex = []
    tempo = [] 
    
    
    with pysam.Samfile(params["in"]+"_"+cx, 'rb') as bam:
    #Retrieve PS Only read1 is written after BAM filtering
        for read1 in bam :

            chr1=bam.getrname(read1.tid)
            ori1, chr1, pos1, ori2, chr2, pos2 = U.getCoord(read1, chr1, chr1)
            if ori1+ori2 == ori:
                ps1sub = U.Subtelo(subtel, chr1, pos1)
                ps2sub = U.Subtelo(subtel, chr1, pos2)

                if (abs(read1.tlen) > median) and (ps1sub + ps2sub <= 1) :
                    tempo.append([read1.qname, ori1, chr1, pos1, ori2, chr2,\
                                   pos2, abs(read1.tlen)])
                    if read1.qname in dicQual:
                        dicQual[read1.qname].append(read1.mapq)
                    else:
                        dicQual[read1.qname] = [read1.mapq]           

    clasamex = sorted(tempo, key = itemgetter(3))
    print "Processing", cx, "- Nb of discordant PS:", \
    str(len(clasamex)), "-", now.strftime("%Y-%m-%d %H:%M")
    #print "Nb of discordant PS", len(clasamex), "("+cx+")"  #,  clasamex
    return clasamex
    
#------------------------------------------------------------------------------
def group_overlaping_old(pos, Lim, clasamex, g_overlap, last, d):
    """ Make group of overlaping PS"""

    end = clasamex[pos][6]
    j = pos + 1
    while j < Lim :
        #Adding +d may implies that the ExtensionSV in create_dup is no more
        #necessary (to be checked)
        if clasamex[j][3] <= end + d :
            g_overlap.append(j)
            last = j
        j = j + 1
    return last
    
#------------------------------------------------------------------------------
def group_overlaping(pos, Lim, clasamex, g_overlap, last, d):
    """ Make group of overlaping PS"""

    end = clasamex[pos][6]
    j = pos + 1
    ok = True
    while ok and j<Lim:
    #while j < Lim :
        #Adding +d may implies that the ExtensionSV in create_dup is no more
        #necessary (to be checked)
        #print "IN THE BOUUUUCLE"
        if clasamex[j][3] <= end + d :
            g_overlap.append(j)
            last = j
            j = j + 1
            if j >= Lim:
                #print "SORTIE DE BOUCLE"
                return last
        else:
            ok = False
            j = j + 1
            if j >= Lim:
                #print "SORTIE DE BOUCLE"
                return last
    #print "SORTIE DE BOUCLE"    
    return last
#------------------------------------------------------------------------------
def TestHomog(ref, candid, clasamex, d, ids):
    """ Make homogeneous PS groups with d and ids parameterss"""

    #Ps candid est-il compatible IS size par rapport a PS ref
    cond_ids = U.parameter_ok(clasamex[candid][7], clasamex[ref][7], ids)
    #Ps est-il compatible sur coordonnees left
    cond_d_left = U.parameter_ok(clasamex[candid][3], clasamex[ref][3], d)
    #Ps est-il compatible sur coordonnees right
    cond_d_right = U.parameter_ok(clasamex[candid][6], clasamex[ref][6], d)
    rajout = cond_ids and (cond_d_left and cond_d_right)

    return rajout

#------------------------------------------------------------------------------
def HomogAll(candid, temp, clasamex, d, ids):
    """ Define PS from temp compatible with PS candid"""

    tl = [temp[candid]]
    for i in range(len(temp)):
        ps_homogen = TestHomog(temp[candid], temp[i], clasamex, d, ids)
        if temp[candid] != temp[i] and ps_homogen:        
            tl.append(temp[i])
            
    tl.sort()
    return tl

#------------------------------------------------------------------------------
def UniqPS(liste, candid):
    """ Return True if PS from candid is not already in liste"""
    uniq = True
    j = 0
    while j < len(candid) and uniq :
        l = 0
        while l < len(liste) and uniq:
            if len(candid) < len(liste[l]):
                uniq = candid[j] not in liste[l]
            l = l + 1
        j = j + 1
    return uniq

#--------------------------------------------------------------------------
def already_SV(groupes, candid):
    """Check if candid is already classified in groupes."""
    already = [indiv for indiv in groupes if indiv == candid]
    return already

#------------------------------------------------------------------------------
def largest_group(clasamex, d, ids, temp, allPS, ps_min):
    """Return maxPS: largest sub-group of homogeneous PS within temp
       Return allPS : true if maxPS is the original group temp
       Return sous : all the sub-groups"""
    maxPS = []
    sous = []
    nj = 0
    while nj < len(temp) and not allPS :
        te = HomogAll(nj, temp, clasamex, d, ids) 
        if len(te) >= ps_min and not already_SV(sous,te) :
            sous.append(te)
            if len(te) == len(temp) :
                #Si tous les PS sont compatibles = 1 groupe complet
                allPS = True
                maxPS = temp
            if len(te) > len(maxPS):
                maxPS = te
        nj = nj + 1
    
    return maxPS, allPS, sous
    
#------------------------------------------------------------------------------
def write_DupbyPS(g, clasamex, manip, ndup, chrx, medianSV, start, end, d):
    """Write Additional duplication to output file _byPS
    and determines extremities and min/max PS length """
    
    size = int(medianSV + d)
#    n = 0
#    for toto in (manip, str(chrx), ndup, len(g), len(g), -1, start, end, end-start, "NA", -1, -1, -1, "NA", size, size,"NA"):
#        print n, toto, type(toto)
#        n = n + 1
    out = open(manip+"_duplications_byPS.csv", "a")
    for nj in g:
        k = clasamex[nj]
        manipw = os.path.split(manip)[1].split("_Ulysses")[0]
        U.write_byPS(out, (manipw, chrx, ndup, len(g), 
                           len(g), -1, start, end, end-start, "NA", -1, -1,
                           -1, "NA", size, size, -1, "NA",  k[0], k[1], k[2], 
                           k[3], k[4], k[5], k[6], k[7]))
    out.close()

#------------------------------------------------------------------------------
def create_dup(maxPS, clasamex, d, ids, manip, chrx, ndup, cladic):
    """ Define an SV after adding compatible PS not catched at first round
        because of the first PS"""
    
    ajout = U.ExtensionSV(maxPS, clasamex, d, ids, cladic)
    if ajout != []:
        maxPS = maxPS+ajout
    medianSV, start, end = U.sizeSV(maxPS, clasamex, "dup")
    write_DupbyPS(maxPS, clasamex, manip, ndup, chrx, medianSV, start, end, d)
    outg = open(manip+"_duplications_bySV.csv", "a")

    U.write_bySV(outg, (os.path.split(manip)[1], chrx, ndup, len(maxPS), 
                        len(maxPS), -1, start, end, end-start, "NA", -1, -1,
                        -1, "NA", max(1,medianSV-d), medianSV+d, -1, "NA"))
    outg.close()


#def cleanDup(maxPS, clasamex, d, ids, manip, chrx, ndup)

#------------------------------------------------------------------------------
def Duplication(clasamex, mad, median, manip, chrx, ndup, ps_min):
    """Detection of duplication >= MEAN PS SIZE"""
    d = median + 6*mad
    ids= 6*mad
    L = len(clasamex)
    #
    #clasamex = sorted(clasamex, key = lambda x : (x[3], x[6]))
    #Put PS in dico instead of list to spead up U.ExtensionSV function
    clasamexDic = {}
    for i,PSS in enumerate(clasamex):
        clasamexDic[i] = PSS
        
    i = 0
    while i < L-1 :
        g_overlap = [i]
        last = i
        #make group of overlaping PS
        last = group_overlaping(i, L, clasamex, g_overlap, last, d)
        if len(g_overlap) >= ps_min :
            #Cherche le plus grand groupe de PS tous compatibles entre eux
            allPS = False
            maxPS, allPS, sous = largest_group(clasamex, d, ids, g_overlap,
                                               allPS, ps_min)
          
            if len(maxPS) >= ps_min:
 
                # Prendre tous les groupes d'au moins deux PS compatibles et
                # non-chevauchants entre eux
                ndup = ndup + 1
                create_dup(maxPS, clasamex, d, ids, manip, chrx, ndup, clasamexDic)
                last = max(maxPS)
                already = maxPS
                #Si plus de 2 groupes
                if not allPS and len(sous) > 1:
                    for sg in sous:
                        sgnr = [ii for ii in sg if ii not in already]
                        if len(sgnr) >= ps_min :
                            #if sg != maxPS:
                            #on autorise le "reuse" des PS entre des
                            # sous-groupes homogenes au sein d'un groupe
                            #si et seulement si il y a au moins 2 PS non
                            # encore utilises dans les sous-groupes deja
                            # definis remain = [isg for isg in sg if isg
                            # not in already] if len(remain) > 1 :
                            ndup = ndup+1
                            create_dup(sgnr, clasamex, d, ids, manip, chrx,
                                       ndup, clasamexDic)
                            already += sg
                            if last < max(sg):
                                last = max(sg)
#                        else:
#                            print "Rejected Dup ", ndup, chrx, len(sgnr)
#            else:
#                print "Largest group is smaller than tolerated group ", len(maxPS)

        if last > i :
            i = last
        else:
            i = i + 1
            
#    #Remove Sub-groups
#    U.cleanAnySVFilePair(manip+"_duplications_bySV.csv", manip+"_duplications_byPS.csv")
    return ndup


#------------------------------------------------------------------------------
def runDetectionDup(params, stats, chrDicos, list_chr_real):
    
# START ALEX PARAMETERS #######################################################

    list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    list_chr_names = "SEP".join(chrDicos.keys())
    
    

    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + float(params["n"]) * stats["mad"]

    #Determine le nb minimum de ps pour chaque type de SV
    detectionFile = "ClusterLimitSize"
    
    ps_min, infoMsg = Ualex.define_MinPS("DUP", stats, params["nsv"], 
        list_chr_length, list_chr_names, detectionFile, params["in"], 
        params["in"]+".dist.table", sval, params["fdr"],params["out"], stats["rl"], 
                                 params["n"], list_chr_real)


# END ALEX PARAMETERS #######################################################

    subtelo=U.get_subtelo_limits(params)
    
    dicQual = {}
                                 
    U.create_files("duplications", params["out"])
    dup = 0
    
    if stats["ps_type"] == "MP":
        inverse = "+-"
    else:
        inverse = "-+"
        
        
    with open(params["out"]+".DUP.report.out","a") as fileo:
        fileo.write("\n----------- \tDetection of Duplication\n")
        fileo.write("\t\tInsert Size Consistency (isc) = %.2f\n" % (6*stats["mad"]))
        fileo.write("\t\tMinimal Insert Size (s) = %d\n" % (stats["median"]))
        fileo.write("\t\tMaximum distance between PS extremities (d) = %.2f\n"\
        % (stats["median"]+6*stats["mad"]))
        fileo.write("\t\tMinimum number of PS to define a Duplication = %d\n\n"\
        % (ps_min))

    print "MINPS for duplications:", ps_min
#    for cx in range(1, 17):        
    for chx in params["range"]:
        if not os.path.isfile(params["in"]+"_"+chx):
            continue
        #select reads in compatible orientation for duplication 
        # "+-" for MP reads and "-+" for PE reads
        #print "Processing ", chx
        classorix = ReadFilesBAM(params, chx, subtelo, stats["median"],
                                 inverse, dicQual)
                                 
#        import pickle
#        output = open('qualities.pkl', 'wb')
#        pickle.dump(dicQual, output)
#        output.close()
        
        
###################################################################################################
#### TEMPORAIRE A CAUSE DES CHIMERES (BAM file: flag pas ok, il peut y avoir 2 primary reads)
        classorix = list(set([tuple(i) for i in classorix]))
######################################################################################################                                 
                                 
        dup = Duplication(classorix, stats["mad"], stats["median"], 
                          params["out"], chx, dup, ps_min)
                          
    print "Detection done"
    with open(params["out"]+".DUP.report.out","a") as fileo:
        fileo.write("Detection processed\n")
        fileo.write("results written to "+params["out"]+"_duplications_by[PS/SV].csv\n")

    #Remove Sub-groups
    U.cleanAnySVFilePair(params["out"]+"_duplications_bySV.csv", params["out"]+"_duplications_byPS.csv")    
    
    #add qualities
    U.addMeanSVQuality(params["out"]+"_duplications_bySV.csv", params["out"]+"_duplications_byPS.csv", dicQual)
    
    toto, msg, pval_seuil_dup = Ualex.functstats("DUP", stats, stats["ndup"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_duplications_bySV.csv", 
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])
    print msg
                                 
    with open(params["out"]+".DUP.report.out","a") as fileo:
        fileo.write(msg+"\n")
    return pval_seuil_dup

#------------------------------------------------------------------------------

def runStatsDup(params, stats, chrDicos):
    n = int(params["n"])

# START ALEX PARAMETERS #######################################################
    list_chr_names = "SEP".join(params["range"])
    list_chr_length = "SEP".join([ str(chrDicos[x]) for x in params["range"] ])
    #list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    #list_chr_names = "SEP".join(chrDicos.keys())
    
    
    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + n * stats["mad"]
    

    print "Statistics on Duplications : Dps = ", stats["ndup"]
    print "False Discovery Rate (fdr) = ", params["fdr"]
    
    toto, msg, pval_seuil_dup = Ualex.functstats("DUP", stats, stats["ndup"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_duplications_bySV.csv", 
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])
                                 
    print msg
    return pval_seuil_dup
#------------------------------------------------------------------------------
def launch(paramfile, onlyStatPerform, list_chr_real):    
    
    if os.path.isfile(paramfile):
        print "parameter file is :", paramfile
        params, stats, chrDicos = U.prepare_detection("duplications", paramfile,
                                                      "NA")
        
        #if ''.join(list_chr_real)!='all':
        #    list_chr_real = [ stats["chromosome_prefix"] + str(x) for x in list_chr_real ]
        
        if onlyStatPerform:
            pval_seuil_dup = runStatsDup(params, stats, chrDicos)
        else:
            pval_seuil_dup = runDetectionDup(params, stats, chrDicos, list_chr_real)                                             
        return pval_seuil_dup
    else:
        print "Error :", paramfile, "doesn't exist"

#------------------------------------------------------------------------------
if (__name__)  ==  "__main__":
    paramfile = U.parser("duplications")
    launch(paramfile)
