# -*- coding: utf-8 -*-
__date__      = "2011/02/21"
#__doc__       = """ texte doc """

import os, sys
import gc
try:
    import pysam
except ImportError:
    print "\tError: Pysam is not installed. \n\
    Please see: https://github.com/pysam-developers/pysam\n\n"
    sys.exit()
try:
    import numpy as np
except ImportError:
    print "\tError: NumPy is not installed. \n\
    Please see: http://www.numpy.org/\n\n"
    sys.exit()
from operator import itemgetter
import operator
import Ulysse_utils as U
import Ulysse_stats as Ualex
import datetime



#--------------------------------------------------------------------------
def ReadFilesBAM(params, cx, orix, dicQual):
    """
    Read BAM Files of discordant PS 
    """
    tempo = []
    ntot = 0
    

    with pysam.Samfile(params["in"]+"_"+cx, 'rb') as bam:
    

    #Retrieve PS Only read1 is written after BAM filtering
        for read1 in bam :
            if read1.is_read1:
                if read1.qname in dicQual:
                    dicQual[read1.qname].append(read1.mapq)
                else:
                    dicQual[read1.qname] = [read1.mapq]
                    
                chr1=bam.getrname(read1.tid)
                ori1, chr1, pos1, ori2, chr2, pos2 = U.getCoord(read1, chr1, chr1)
                if ori1+ori2 == orix and read1.mapq :
                    tempo.append([read1.qname, ori1, chr1, pos1, ori2, chr2, \
                                   pos2, abs(read1.tlen)])

    clasamex = sorted(tempo, key = itemgetter(3))
    return clasamex, ntot
#----------------------------------------------------------------------------                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
def ReadFilesPaired(fil, orix) :
    """ Read Files in format obtained after filtering of data formatted like
    Paired.uniquelymapped_noDuplicates.SV through filterPS.py
    """
    tempo = []
    ntot = 0
    with open(fil, 'r') as gile:
        for line in gile:
            l = line[:-1].split(";")
            sor = l[1] + l[4]
            if sor == orix :
                tempo.append(l[0:3] + [int(l[3]), l[4], l[5], int(l[6]), int(l[7])])
                ntot = ntot + 1
    clasamex = sorted(tempo, key = itemgetter(3))
    return clasamex, ntot

#----------------------------------------------------------------------------
def all_conditions(ps1, ps2, d, ids, mimx) :
    """ Return True if conditions on right / left coordinates are ok
    and conditions on PS length are ok
    """

    left_ok = U.parameter_ok(ps1[3], ps2[3], d)
    right_ok = U.parameter_ok(ps1[6], ps2[6], d)
    ik = U.parameter_ok(ps2[7], mimx[0], ids)
    ik = ik and U.parameter_ok(ps2[7], mimx[1], ids)
    #if U.ps_fictive(ps1) and U.ps_fictive(ps2):
        #print "all_conditions ", ps1[0], ps2[0], left_ok, right_ok, ik
    #if "Deletion" in ps1[0] or "Deletion" in ps2[0]:
        #print left_ok, right_ok, ik

    #overlap_ok = ps2[3] <  ps1[6] inutile car on a deja des ps chevauchants
    #return left_ok and right_ok and ik and overlap_ok
    return left_ok and right_ok and ik
#-----------------------------------------------------------------------------
def min_max_ps_size(val, mimx) :
    """ Return list containing min and max insert size"""
    mimx[0] = min(mimx[0], val)
    mimx[1] = max(mimx[1], val)
#------------------------------------------------------------------------------
def original_group(liste):
    """ Sub-group equals original group"""
    selected = liste
    final = [liste]
    return np.array(final), selected
#------------------------------------------------------------------------------
def one_group(select, liste, sel, endmax, ngmax, d, ids, listeDic):
    """ Only one Sub-group"""
    if endmax in sel[ngmax]:
        listex = U.ExtensionSV(sel[ngmax], select, d, ids, listeDic)
        sel[ngmax] += listex
    selected = [liste[i] for i in sel[ngmax]]
    final = [[liste[i] for i in sel[ngmax]]]
    return np.array(final), np.array(selected)
#------------------------------------------------------------------------------
def n_groups(select, liste, sel, endmax, ngmax, nb, d, ids, listeDic):
    """ More than one Sub-group: give all independant ones"""
    selected = [liste[i] for i in sel[ngmax]]
    final = [[liste[i] for i in sel[ngmax]]]
    listex = []
    while nb :
        del nb[ngmax]
        for g in nb.keys():
            #if g != ngmax:
            t = [liste[i] for i in sel[g] if i in sel[ngmax]]
            #Si un PS se retrouve dans le groupe de taille max
            if t or nb[g] < 2 :
                del nb[g]

        if nb :
            ngmax = nb.keys()[0]
            for g in nb:
                if nb[g] > nb[ngmax]:
                    ngmax = g
            if nb[ngmax] > 1 :
                final.append([liste[i] for i in sel[ngmax]])
                selected += [liste[i] for i in sel[ngmax]]


    for fi in final:
        if endmax in fi:
            listex = U.ExtensionSV(fi, select, d, ids, listeDic)
            fi += listex
            selected += listex
#    if listex:
#        print "Il y a eu rajout. A la fin PS_Select_Del "

    return np.array(final), np.array(selected)

#------------------------------------------------------------------------------
def create_groups(select, liste, sel, nb, ngmax, endmax, d, ids, listeDic):
    """ Define all homogeneous group or subgroups. return the list of groups
        and list of PS included in the groups"""

    selected = []
    final = [[]]
    if nb[ngmax] > 1 :
        if len(sel[ngmax]) == len(liste):
            #No PS eliminated
            final, selected = original_group(liste)
        elif len(sel[ngmax]) == 1 :
            #Only one subgroup
            final, selected = one_group(select, liste, sel, endmax, ngmax, d, ids, listeDic)
        elif len(sel[ngmax]) >= 2 :
            #2 or more sub-groups
            final, selected = n_groups(select, liste, sel, endmax, ngmax, nb, d, ids, listeDic)

    return final, selected

#------------------------------------------------------------------------------
#comprehension list
def Select_PS_Del(select, liste, d, ids, listeDic):
    """Selectionne les PS qui remplissent conditions sur d
    start(-)2 - start(-)1 < d end(+)2 - end(+)1   < d"""

    L = len(liste)
    sel = []
    nb = {}
    ngmax = 0
    endmax = 0
    i = 0
    for i, li in enumerate(liste[:-1]):
        #Min and Max PS insert sizes
        if U.ps_fictive(select[li]):
            print "DEL_PS_fictive", select[li]
        mimax = [select[li][7], select[li][7]]
        endmax = i
        num = L - i
        tri = [i]
        #Create homogeneous group of PS
        #Keep in endmax the left far most PS
        for j, lj in enumerate(liste[i + 1:]):

            if all_conditions(select[li], select[lj], d, ids, mimax):
                tri.append(i+j+1)

                min_max_ps_size(select[lj][7], mimax)

            else:
                num = num - 1
            if select[lj][6] > select[liste[endmax]][6]:
                endmax = i+j+1
            if num < 1:
                break

        #Keep PS with at least one compatible PS
        sel.append(np.array(tri))
        nb[i] = num
        if num > len(sel[ngmax]):
            ngmax = i
        if sel and len(sel[ngmax]) == L :
            break


    return create_groups(select, liste, sel, nb, ngmax, endmax, d, ids, listeDic)

#------------------------------------------------------------------------------
def AlreadyDel(deletion, elem):
    """ Return True if elem already classified in deletion"""
    temp = []
    for d, dele in enumerate(deletion):
        temp = [i for i in dele if i == elem]
    if temp == []:
        return False
    else:
        return True

#------------------------------------------------------------------------------
def WriteSV(manip, outf, xsome, deletion, select, d, mediane, typesv):
    """ Write Deletion candidates in output files"""
    out = open(outf+"_"+typesv+"_byRP.csv", "a")
    out2 = open(outf+"_"+typesv+"_bySV.csv", "a")
    ndup = 0
#    for j in deletion:
#        l = deletion[j]
    for l in deletion:
        ndup = ndup + 1
        medianSV, start, end = U.sizeSV(l, select, "del")
        if typesv == "deletions":
            size = medianSV - mediane
        else:
            size = medianSV + d
        for kk in l :
            nbps = len(l)
            k = select[kk]
            U.write_byPS(out, (manip, xsome, ndup, nbps, nbps, -1, 
                               start, end, abs(end-start), "NA", -1, -1, -1,
                               "NA", size, size, -1, k[0], k[1], k[2], k[3],
                               k[4], k[5], k[6], k[7]), True)


                

        U.write_bySV(out2, (manip, xsome, ndup, nbps, nbps, -1, 
                               start, end, abs(end-start), "NA", -1, -1, -1,
                               "NA", size, size, -1), True)

    out.close()
    out2.close()

#------------------------------------------------------------------------------
def Deletion(params, xsome, liste, stats, chrDicos, ps_min, listeDic):

    """ Detection of deletion"""
    
    manip = params["in"]
    outf = params["out"]
    median = stats["median"]
    mad = stats["mad"]
    n = int(params["n"])
    #TO REMOVE
    #n = 10
    

    ids = n*mad # difference max pour taille des inserts entre PS formant une deletion (= ids)
    d = median + n*mad # difference max pour coordonnees des inserts d'un meme cote d'une deletion

#    deletion = {}
    deletion = []  #pour utiliser numpy
    ntot = 0
    #Recupere les reads >sizemin
    #Prend taille de l'insert (la colonne j[7])
    sizemin = median + n*mad # taille minimale des PS (=s)
#    if len(liste) !=0:
#        print "avant", len(liste)
    #print sizemin
    #print liste[0], liste[0][7], liste[0][7]>sizemin
    #print liste[1], liste[1][7], liste[1][7]>sizemin
    select2 = [j for j in liste if j[7]  > sizemin]
    select = sorted(select2, key=operator.itemgetter(3))
#    if len(select) !=0:
#        print "apres", len(select)
    #print select[0]
    #k = 0
    #L = len(select)
    #while k < L :
    bypass = -1
    E = enumerate(select)
    for k, selk in E:
        if AlreadyDel(deletion, k) or k < bypass :
            continue
            #k = k + 1
        else:
            bypass = -1
            temp = [k]
            #j = k + 1
            #reads chevauchants uniquement
            K = enumerate(select[k + 1:])
            temp0 = temp + [k+j+1 for j,selj in K if selj[3] <= selk[6]]
            temp = [j for j in temp0 if not AlreadyDel(deletion, j) ]
            #Selection des deletions
            #verifie que les coordonnees des differentes PS ne different pas de plus de d
            if len(temp) == 1 :
                continue
                #k = k + 1
            else:
#                deletion_groups , selected = Select_PS_Del(select, temp, d, ids)
                deletion_groups , selected = Select_PS_Del(select, 
                                                           np.array(temp), d,
                                                           ids, listeDic)
                #aucune deletion detectee on repart de k+1

                if deletion_groups == [[]] and selected == [] :
                    continue
                    #k = k + 1
                else:
                    if len(deletion_groups) == 1 and len(deletion_groups[0]) >= ps_min :
                        #une seule deletion formee par tout ou un sous-groupe de PS chevauchants
                        ntot = ntot + 1
#                        deletion[ntot] = deletion_groups[0]
                        deletion.append(deletion_groups[0])
                        
                        #on repart du premier PS apres le groupe de del
                    else:
                        #plusieurs deletions
                        for nd in deletion_groups:
                            if len(nd) >= ps_min:                            
                                ntot = ntot + 1
                                deletion.append(nd)

                    #On repart du premier PS apres le PS le plus "a droite" dans les groupes de del
                    #Ca merdait ici????
                    #k = max(selected) + 1
                    bypass = max(selected) + 1
                    del deletion_groups
    #print "++++", len(deletion)
#    for g in deletion:
#        if len(g)==4:
#            print g
    WriteSV(os.path.split(manip)[1], outf, xsome, deletion, select, d, median, "deletions")
    del deletion
    gc.collect()


def Deletion_v2(params, xsome, liste, stats, chrDicos, ps_min, listeDic):

    """ Detection of deletion"""
    
    manip = params["in"]
    outf = params["out"]
    median = stats["median"]
    mad = stats["mad"]
    n = int(params["n"])
    

    ids = n*mad # difference max pour taille des inserts entre PS formant une deletion (= ids)
    d = median + n*mad # difference max pour coordonnees des inserts d'un meme cote d'une deletion

#    deletion = {}
    deletion = []  #pour utiliser numpy

    #Recupere les reads >sizemin
    #Prend taille de l'insert (la colonne j[7])
    sizemin = median + n*mad # taille minimale des PS (=s)
    #print "avant", len(liste)
    select2 = [j for j in liste if j[7]  > sizemin]
    #print select2[0]
    #print "apres", len(select2)

    
     #sort the list by 1) coord left  coord
    select = sorted(select2, key=operator.itemgetter(3))

    deletion =[]
    id_ps = []
    
    #for each PS, search for the compatible ones in the same orientation
    #it makes homogeneous groups
    for i, psi in enumerate(select):
        #print "+++", i, psi
        #a ps is obviously compatible with at least itself...
        del_i = [i]
        id_ps_tmp=[]
        
        if i in set(id_ps):
            continue
            #break

        #PS are sorted primarily by left coord and secondarily by right coord
        #Therefore, for a given PS, compatible PS are the ones right before
        #(to the left) and right after (to the right). We search for compatible
        #PS to the left (and right) until we find an uncompatible PS on the
        #left coord.
        #initialize for left and right search
        i1 = i+0

        #left search
        while 1:
            if i1<len(select)-1:
                i1=i1+1
            else:
                break
            #print psi[3], select[i1][3], d
            coord_left_ok = (abs(psi[3]-select[i1][3]) <= d)
            if coord_left_ok:
                IS_ok = (abs(psi[7]-select[i1][7]) <= ids)
            else:
                break
            if coord_left_ok and IS_ok and psi != select[i1]:
                #print "comp"
                del_i.append(i1)
                #print del_i
                id_ps_tmp.append(i)
               
        if len(del_i)>=ps_min:
            deletion.append(del_i)
            for ps in id_ps_tmp:
                id_ps.append(ps)
    #Remove sub groups
    deletion = U.filter_list(deletion)    
    #print select[0], d, ids
    #print select[1], d, ids
    
    #print id_ps 
           
    #In groups of deletion, all PS are compatible with the first one
    #now we need to extract the largest or the largests subgroups(s)
    del_true = [] #Warning for now only keeps largest subgroup (not all combinations (but should be very very rare))
    for g in deletion:
        PSs = [select[x] for x in g] #get PS from their ID
        ISs = [x[7] for x in PSs] #list if the IS
        if abs(min(ISs)-max(ISs) <= ids):
            del_true.append(g) #case were all IS are compatible
        else:
            PSs= sorted(PSs, key=operator.itemgetter(7))
            sub_groups = []
            for i, ps in enumerate(PSs):
                i1 = i+0
                tmp = [i]
                while 1:
                    if i1<len(PSs)-1:
                        i1=i1+1
                    else:
                        break
                    #print psi[3], select[i1][3], d
                    IS_ok = (abs(ps[7]-PSs[i1][7]) <= ids)
                    if IS_ok:
                        tmp.append(i1)
                    else:
                        break  
                if len(tmp)>=ps_min:
                    sub_groups.append(tmp)
            if sub_groups:
                del_true.append(max(sub_groups, keu=len))
    
    #merge groups that should be together    
    del_fuse = [] 
    gtmp = []
    G = enumerate(del_true)    
    for i, g in G:
        if i < len(del_true)-2: #index value error check
            if not gtmp: #if no temporary group of merged groups
                if (g[-1] in del_true[i+1]) and (g[-2] in del_true[i+1]): #check if groups overlapps with next one
                    gtmp = list(set(g + del_true[i+1])) #merge them in new temporary group
                    gtmp.sort() # sort it (coords and IDs are both increasing)
                else:
                    del_fuse.append(g) #if they do not overlapp, merge them

            else: #if their is already a temporary merge group
                if (gtmp[-1] in g) and (gtmp[-2] in g): #check if group and temporary group overlapp
                    gtmp = list(set(g + gtmp)) #merge them in new temporary group
                    gtmp.sort()
                else:
                    del_fuse.append(gtmp) #if not, write the temporary merge group
                    gtmp = []
                    if (g[-1] in del_true[i+1]) and (g[-2] in del_true[i+1]): #and check if current group overlapp with the next one
                        gtmp = list(set(g + del_true[i+1])) #if they do, creat new temporary group
                        gtmp.sort()
                    else:
                        del_fuse.append(g) #if not, just write the current gruop
        else:
            if gtmp:
                del_fuse.append(gtmp)
            else:
                del_fuse.append(g)
    ###
    
    #print len(g)
    print "-->", len(deletion)
    print "--->", len(del_true)
    #print del_fuse    
    print "--->", len(del_fuse)
    #print deletion
    WriteSV(os.path.split(manip)[1], outf, xsome, del_fuse, select, d, median, "deletions")
    del deletion
    gc.collect()



#------------------------------------------------------------------------------
def SmallIns(params, xsome, liste,stats, chrDicos, ps_min, listeDic):
    """ Detection of small insertions (initially duplication)"""
    
    manip = params["in"]
    outf = params["out"]
    median = stats["median"]
    mad = stats["mad"]
    #n = stats["nsins"]
    n = int(params["n"])
    #print "median", "n", "mad", median, n, mad
    sizemax = median - n*mad # taille minimale des PS (=s)
    ids = n*mad # difference max pour taille des inserts entre PS formant une deletion (= ids)
    d = median + n*mad # difference max pour coordonnees des inserts d'un meme cote d'une deletion
    dup = []
    ntot = 0
    #Recupere les reads >sizemin
    #Prend taille de l'insert (la colonne j[7])
    select = [j for j in liste if j[7]  < sizemax]
    #print "sizemax", sizemax
    k = 0
    L = len(select)
    #print "LLLLL", L
    #while k < L :
    bypass = -1
    E = enumerate(select)
    for k, selk in E:
        if AlreadyDel(dup, k) or k < bypass :
            continue
            #k = k + 1
        else:
            bypass = -1
            temp = [k]
            #j = k + 1
            #reads chevauchants uniquement
            K = enumerate(select[k + 1:])
            temp0 = temp + [k+j+1 for j,selj in K if selj[3] <= selk[6]]
            temp = [j for j in temp0 if not AlreadyDel(dup, j) ]
            #Selection des deletions
            #verifie que les coordonnees des differentes PS ne different pas de plus de d
            if len(temp) == 1 :
                continue
                #k = k + 1
            else:
                dup_groups , selected = Select_PS_Del(select, np.array(temp), d, ids, listeDic)
                #aucune deletion detectee on repart de k+1
                if dup_groups == [[]] and selected == [] :
                    continue
                    #k = k + 1
                else:
                    if len(dup_groups) == 1 and len(dup_groups[0]) > ps_min :
                        #une seule deletion formee par tout ou un sous-groupe de PS chevauchants
                        ntot = ntot + 1
                        dup.append(dup_groups[0])
                        #on repart du premier PS apres le groupe de del
                    else:
                        #plusieurs deletions
                        for nd in dup_groups:
                            if len(nd) > ps_min:                            
                                ntot = ntot + 1
                                dup.append(nd)

                    #On repart du premier PS apres le PS le plus "a droite" dans les groupes de del
                    #Ca merdait ici????
                    #k = max(selected) + 1
                    bypass = max(selected) + 1
                del dup_groups

#    print "Small Dup Final", dup, xsome
    WriteSV(os.path.split(manip)[1], outf, xsome, dup, select, d, median, "small_insertions")
    del dup
    #Remove Sub-groups
    #print "OUTF", outf
    #U.cleanAnySVFilePair(outf+"_small_insertions_bySV.csv", outf+"_small_insertions_byPS.csv")
    
    gc.collect()

#------------------------------------------------------------------------------        
def runStatsDel(params, stats, chrDicos):    
    
    
    n = int(params["n"])

# START ALEX PARAMETERS #######################################################

    list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    list_chr_names = "SEP".join(chrDicos.keys())
    
    
    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + n * stats["mad"]
    
    
    print "Statistics on Deletion : Dps = ", stats["ndel"]
    print "False Discovery Rate (fdr) = ", params["fdr"]
    #Calcul stats only
    pval_seuil_sins = -1
    toto, msg, pval_seuil_del = Ualex.functstats("DEL", stats, stats["ndel"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_deletions_bySV.csv", 
                                 params["in"], params["in"]+"_dist.table", 
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])
    print msg

    if stats["median"] > 10 and (stats["median"] - int(params["n"])*stats["mad"] > 0):
        toto, msg, pval_seuil_sins = Ualex.functstats("sINS", stats, stats["ndel"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_small_insertions_bySV.csv", 
                                 params["in"], params["in"]+"_dist.table", 
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])   

                                 
    print msg
    #print "pvals", pval_seuil_del, pval_seuil_sins
    return pval_seuil_del, pval_seuil_sins
#------------------------------------------------------------------------------        
def runDetectionDel(params, stats, chrDicos):    
    
    
    if stats["ps_type"] == "MP":
        orientation = "-+"
    else:
        orientation = "+-"

    n = int(params["n"])
    #print "NNNN", n


    U.create_files("deletions", params["out"], True)
    #if stats["ps_type"] == "MP":
    U.create_files("small_insertions", params["out"], True)


# START ALEX PARAMETERS #######################################################
    list_chr_names = "SEP".join(params["range"])
    list_chr_length = "SEP".join([ str(chrDicos[x]) for x in params["range"] ])
    
    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + n * stats["mad"]
    

    #Determine le nb minimum de ps pour chaque type de SV
    detectionFile = "ClusterLimitSize"
    ps_min, infoMsg = Ualex.define_MinPS("DEL", stats, params["nsv"], 
        list_chr_length, list_chr_names, detectionFile, params["in"] ,
        params["in"]+"_dist.table", sval, params["fdr"], params["out"], stats["rl"], 
        params["n"], params["range"])

    with open(params["out"]+".DEL.report.out","a") as fileo:
        fileo.write("\n-----------\tDetection of Deletion\n")
        fileo.write("\t\tInsert Size Consistency (isc) = %.2f\n" % (n * stats["mad"]))
        fileo.write("\t\tMinimal Insert Size = %d\n" % (stats["median"] + n*stats["mad"]))
        fileo.write("\t\tMaximum distance between PS extremities (D) = %.2f\n"\
        % (stats["median"] + n*stats["mad"]))
        fileo.write("\t\tFDR threshold (d) = %.2f\n" % (params["fdr"]))
        fileo.write("\n Minimum number of PS to define a deletion : %d" % (ps_min))
           
    if ps_min == 1:
        ps_min = 2   
    #ps_min =3
    print "\nMINPS for a deletion:", ps_min


# END ALEX PARAMETERS #######################################################
    dicQualDEL = {}
    dicQualsINS = {}
    for cx in params["range"]:
#    for cx in range(1, 17):
        if not os.path.isfile(params["in"]+"_"+cx):
            continue
        #print "\n Scanning ", cx
        #os.system("date")
        classorix, ntot  = ReadFilesBAM(params, cx, orientation, dicQualDEL)
        #classorix = list(set([tuple(i) for i in classorix]))

        
        #print "\nScanning ", cx,"with", len(classorix), "PS"
        print "Processing", cx, "- Nb of discordant PS:", \
    str(len(classorix)), "-", datetime.datetime.now() 
        
        classorixDic = {}
        for i,PSS in enumerate(classorix):
            classorixDic[i] = PSS        
        
        Deletion(params, cx, classorix, stats, chrDicos, ps_min, classorixDic)
        #Deletion_v2(params, cx, classorix, stats, chrDicos, ps_min, classorixDic)
        del classorix
        gc.collect()

    if stats["median"] > 10000 :
        #if stats["median"] > 1500 :
               
        ps_min, infoMsg = Ualex.define_MinPS("sINS", stats, params["nsv"], 
        list_chr_length, list_chr_names, detectionFile, params["in"] ,
        params["in"]+"_dist.table", sval, params["fdr"], params["out"], stats["rl"], 
        params["n"])
        print "\nMINPS for a small insertion:", ps_min


        with open(params["out"]+".sINS.report.out","a") as fileo:
            fileo.write("\n-----------\tDetection of Small Insertion\n")
            fileo.write("\t\tInsert Size Consistency (ics) = %.2f\n" % (n * stats["mad"]))
            fileo.write("\t\tMaximal Insert Size  = %d\n" % (n *stats["median"]))
            fileo.write("\t\tMaximal number of SV candidate = %s\n" % (params["nsv"]))
            fileo.write("\t\tMaximum distance between PS extremities (D) = %.2f\n"\
                         % (stats["median"] + n*stats["mad"]))
            fileo.write("\n Minimum number of PS to define a Small Insertion : %d" % (ps_min))

        for cx in params["range"]:
                if not os.path.isfile(params["in"]+"_"+cx):
                    continue
                classorix, ntot = ReadFilesBAM(params, cx, orientation, dicQualsINS)
                print "Processing", cx, "- Nb of discordant PS:", \
    str(len(classorix)), "-", datetime.datetime.now() 
                classorixDicIns = {}
                for i,PSS in enumerate(classorix):
                    classorixDicIns[i] = PSS 
                SmallIns(params, cx, classorix, stats, chrDicos, ps_min, classorixDicIns)
                #os.system("date")
                del classorix
                gc.collect()


    print "Detection done\n\nReport file : "+params["out"]+".DEL.report.out"
    with open(params["out"]+".DEL.report.out","a") as fileo:
        fileo.write("Detection processed\n")
        fileo.write("results written to "+params["out"]+"_deletions_by[PS/SV].[stats.]csv\n")
        print "results written to "+params["out"]+"_deletions_by[PS/SV].csv\n"
    if stats["median"] > 1000 :
            #print "Detection done\n Report file : "+params["out"]+".DEL.report.out"
            with open(params["out"]+".sINS.report.out","a") as fileo:
                fileo.write("\and "+params["out"]+"_small_insertions_by[PS/SV].csv\n")
                fileo.write("results written to "+params["out"]+"_deletions_by[PS/SV].[stats.]csv")
                print "results written to ", params["out"]+"_deletions_by[PS/SV].csv"

    

    os.system("date")
    
    #BAM = bam d'origine, prefixe fichier res del
    pval_seuil_sins = -1
    Ualex.addCovToDelFile(params["in"], params["out"]+"_deletions_bySV.csv")
    
        
    #add qualities
    U.addMeanSVQuality(params["out"]+"_deletions_bySV.csv", params["out"]+"_deletions_byRP.csv", dicQualDEL)
    
        
    toto, msg, pval_seuil_del = Ualex.functstats("DEL", stats, stats["ndel"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_deletions_bySV.csv", 
                                 params["in"], params["in"]+"_dist.table", 
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])
    print msg
                                 
    if int(stats["median"]) > 10 and (int(stats["median"]) - int(params["n"])*int(stats["mad"]) > 0):
        Ualex.addCovToDelFile(params["in"], params["out"]+"_small_insertions_bySV.csv")
        #add qualities
        U.addMeanSVQuality(params["out"]+"_small_insertions_bySV.csv", params["out"]+"_small_insertions_byRP.csv", dicQualsINS)
        
        toto, msg, pval_seuil_sins = Ualex.functstats("sINS", stats, stats["ndel"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_small_insertions_bySV.csv", 
                                 params["in"], params["in"]+"_dist.table", 
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])
        #print "pval_seuil_sins", pval_seuil_sins
                                 

    print msg
    return pval_seuil_del, pval_seuil_sins

#------------------------------------------------------------------------------
def launch(params, stats, chrDicos):
    
    if params["only_stats"]:
        pval_seuil_del, pval_seuil_sins = runStatsDel(params, stats, chrDicos)
    else:
        pval_seuil_del, pval_seuil_sins = runDetectionDel(params, stats, chrDicos)
    return pval_seuil_del, pval_seuil_sins
