# -*- coding: utf-8 -*-
__date__      = "2011/02/21"
#__doc__       = """ texte doc """

import os, sys
import datetime
import itertools
#now = datetime.datetime.now()


try:
    import pysam
except ImportError:
    print "\tError: Pysam is not installed. \n\
    Please see: https://github.com/pysam-developers/pysam\n\n"
    sys.exit()
try:
    from numpy import median
    import numpy as np
except ImportError:
    print "\tError: NumPy is not installed. \n\
    Please see: http://www.numpy.org/\n\n"
    sys.exit()

try:
    import networkx as nx
except ImportError:
    print "\tError: networkx package in not installed. \n\
    Please see: http://networkx.lanl.gov/index.html"

import Ulysse_utils as U
import Ulysse_stats as Ualex

#--------------------------------------------------------------------------
def ReadFilesBAM(params, cx, clasame, dicQual):
    """
    Read BAM Files of discordant PS 
    """
    
    with pysam.Samfile(params["in"]+"_"+cx, 'rb') as bam:
    

    #Retrieve PS Only read1 is written after BAM filtering
        for read1 in bam :
            if read1.is_read1:

                chr1=bam.getrname(read1.tid)
                ori1, chr1, pos1, ori2, chr2, pos2 = U.getCoord(read1, chr1, chr1)
                                                                                    
                if ori1+ori2 == "--" or ori1+ori2 == "++":
                    clasame[ori1+ori2].append([read1.qname, ori1, chr1, pos1, \
                                              ori2, chr2, pos2, abs(read1.tlen)])
                    if read1.qname in dicQual:
                        dicQual[read1.qname].append(read1.mapq)
                    else:
                        dicQual[read1.qname] = [read1.mapq]



#------------------------------------------------------------------------------
def ReadFilesPaired(fil, clasame, pslist) :
    """ Read Files in format obtained after filtering of data formatted like
        Paired.uniquelymapped_noDuplicates.SV through filterPS.py
    """
    with open(fil, 'r') as gile:
        for line in gile:
            l = line.split(";")
            sor = l[1]+l[4]
            if sor !=  "+-" :
                clasame[sor].append( l[0:3] + [int(l[3]), l[4],
                                      l[5], int(l[6]), int(l[7])])
                if sor  ==  "--" or sor  ==  "++":
                    pslist.append(l[0:3] + [int(l[3]), l[4], l[5], int(l[6]),
                                  int(l[7])])

#------------------------------------------------------------------------------
def NoOverlap(cand, compo, d):
    """Definit si un PS dans l'ensemble compo chevauche le PS cand avec
    une autre orientation des reads"""
    coord1 = cand[3]
    coord2 = cand[6]
    patt = cand[1]+cand[4]
    no = 0
    overreads = []
    others = compo.keys()
    others.remove(patt)
    for pat in others:
        for li in compo[pat]:
            if (abs(coord1-li[3]) <= d) and  (abs(coord2-li[6]) <= d) :
                no += 1
                overreads.append(li)
    return no, overreads

#------------------------------------------------------------------------------
def extremites(clasamex, liste, signe):
    """determine le min et max des coordonnees"""
    
    deb = clasamex[signe][liste[0]][3]
    fin = clasamex[signe][liste[0]][6]
    maxleft = deb
    minright = fin
    if len(liste) == 1:
        taille = int(clasamex[signe][liste[0]][7])
    else:
        taille = []
        for i in range(1, len(liste)):
            #print i
            #print liste[i]
            #print clasamex
            #print clasamex[signe][liste[i]]
            taille.append(int(clasamex[signe][liste[i]][7]))

            if clasamex[signe][liste[i]][3] < deb :
                deb = clasamex[signe][liste[i]][3]
            if clasamex[signe][liste[i]][3] > maxleft:
                maxleft = clasamex[signe][liste[i]][3]
                
            if clasamex[signe][liste[i]][6] > fin :
                fin = clasamex[signe][liste[i]][6]
            if clasamex[signe][liste[i]][6] < minright :
                minright = clasamex[signe][liste[i]][6]
        taille = median(taille)
            
    return deb, fin, maxleft, minright, taille

#------------------------------------------------------------------------------
def candidatInvMoinsPlus_PE(m, p, ids, d):
    """Return True if pair of PS (+,+) et (-,-) ok
    to form the inversion junction"""
    
    coord_left_ok = (m[3] > p[3]) and U.parameter_ok(p[3], m[3], 2*d)
    coord_right_ok = (m[6] > p[6]) and U.parameter_ok(p[6], m[6], 2*d)
    ids_ok = U.parameter_ok(m[7], p[7], ids)
    
    if U.ps_fictive(m) and U.ps_fictive(p) :
        print "candidatInvMoinsPlus ", m, p, m[3], p[3], m[6], p[6], (m[6] < p[6]), U.parameter_ok(p[6], m[6], 2*d)
    if coord_left_ok and coord_right_ok and ids_ok :
        return True
    else:
        return False

#------------------------------------------------------------------------------
def candidatInvMoinsPlus_MP(m, p, ids, d):
    """Return True if pair of PS (+,+) et (-,-) ok
    to form the inversion junction"""
    
    coord_left_ok = (m[3] < p[3]) and U.parameter_ok(p[3], m[3], 2*d)
    coord_right_ok = (m[6] < p[6]) and U.parameter_ok(p[6], m[6], 2*d)
    ids_ok = U.parameter_ok(m[7], p[7], ids)
    
    if U.ps_fictive(m) and U.ps_fictive(p) :
        print "candidatInvMoinsPlus ", m, p, m[3], p[3], m[6], p[6], (m[6] < p[6]), U.parameter_ok(p[6], m[6], 2*d)
    if coord_left_ok and coord_right_ok and ids_ok :
        return True
    else:
        return False

#------------------------------------------------------------------------------
def candidatInvMemeSigne(c1, c2, d, ids):
    """Return True if pair of PS in the same orientation (+, +) or (-, -) ok
    to be part of the inversion"""
    maxL = max(c1[3], c2[3])
    minR = min(c1[6], c2[6])
    coord_left_ok = U.parameter_ok(c1[3], c2[3], d)
    coord_right_ok = U.parameter_ok(c1[6], c2[6], d)
    ids_ok = U.parameter_ok(c1[7], c2[7], ids)
    if maxL < minR and  coord_left_ok  and coord_right_ok and ids_ok :
        return True
    else:
        return False

#------------------------------------------------------------------------------
def afficheCand(l1, l2, fil):
    """Print PS candidates in the _byPS.csv file"""
    out = open(fil, "a")
    for i in (l1, l2):
        out.write("%s;%s;%s;%d;%s;%s;%d;%d;"% (i[0], i[1], i[2], i[3], i[4],
                                               i[5], i[6], i[7]))
    out.write("\n")
    out.close()

#------------------------------------------------------------------------------
def common(list1, list2):
    """If members of list1 are present in list2, create a the intersection list
    between list1 and list2"""
    add = True
    ad = ()
    if list1 not in list2: ##THIS DOES NOT WORK AS EXPECTED ???!!
        for k in list2:
            inter = [j for j in list1 if j in k]
            if inter:
                ad = (inter, k)
        if ad :
            list2.remove(ad[1])
            list2.append(ad[0])
            add = False
    else:
        add = False

    if add :
        list1.sort()
        list2.append(list1)

#------------------------------------------------------------------------------
def sumCom(summary, g1, g2):
    """Returns True if [g1,g2] are already defined in summary"""
    already = False
    for si in summary:
        if (si[0]  ==  g1) and (si[1]  ==  g2):
            already = True
    return already


#------------------------------------------------------------------------------
def rmKeyDic(dicos, minPS):
    """ Removes small groups from dicos"""
    nDicos = {}
    for key, val in dicos.iteritems():
        if len(val)>=minPS:
            nDicos[key]=val
    return nDicos

#------------------------------------------------------------------------------
def rmSubInv(liste):
    """remouve INV contained in each other """
    noDups =  map(U.tuples2list, set(map(U.list2tuples,liste)))
    cleanListe = []
    for i, elt in enumerate(noDups):
        li, ri = set(elt[0]), set(elt[1])
        ok = True
        for j,eltj in enumerate(noDups):
            lj, rj = set(eltj[0]), set(eltj[1])
            if li.issubset(lj) and ri.issubset(rj) and i!=j:
                ok = False
                break
        if ok:
            cleanListe.append(elt)
    
    return cleanListe

#----------------------------------------------------------------------------

def inversions_Etape1et2(clasamex, d, ids, ps_type, ps_min):
    """Step 1 and 2 of detection

    si les mate-pairs ne sont pas classes par coordonnees il faut parcourir
    tous les reads -- pour chaque ++
    recherche tous les -- qui vont bien avec les ++
    dictionnaire cle = "-, -" et valeurs : copains +, +
    dictionnaire cle = "+, +" et valeurs : copains -, -
    
    
    
    Cherche les PS compatibles en orientation opposees
    2 dicos (un pour les -- l autre pour les ++:
    ID PS : [liste des compatibles]"""

    candmoins = {}
    candplus = {}
    for i, canmi in enumerate(clasamex["--"]):
        for j, canpj in enumerate(clasamex["++"]):
            if ps_type == "PE":
                candidat = candidatInvMoinsPlus_PE(canmi, canpj, ids, d)
            else:
                candidat = candidatInvMoinsPlus_MP(canmi, canpj, ids, d)
                
            if candidat:
#                afficheCand(canmi, canpj, "candidats-moins-plus")
                if i not in candmoins:
                    candmoins[i] = [j]
                else:
                    candmoins[i].append(j)
                if j not in candplus:
                    candplus[j] = [i]
                else:
                    candplus[j].append(i)
    
    #Remove small groups    
    candmoins = rmKeyDic(candmoins, ps_min)
    candplus = rmKeyDic( candplus, ps_min)
    #print "etape1-2", len(candmoins), len(candplus)

    return candmoins, candplus


#------------------------------------------------------------------------------
def inversions_Etape4(clasamex, cand, signe, d, ids, minPS):
    """ Group PS of same signs ("++" or "--") according to their PS in opposite
    orientation that they have in common
    /!\ It is not because the PS are compatible with the same opposite PS that they
    are themselves compatible !!
    """
    homogene = {}
    oppose = "--"
    if signe  ==  "--":
        oppose = "++"
    for un in cand.keys():
        toto = inversions_Etape5(clasamex, cand[un], oppose, d, ids, minPS) #fait des sousgroups homogenes dans l orientation oppose a la PS 'un'
        #if un == 175:
        #    print "ca commence", toto
        if toto :
            toto2=[]
            for elt in toto:
                toto2.append(elt)
                #if len(elt)>=minPS:
                 #   toto2.append(elt)
            
            if toto2:
                #print "toto", toto2
                homogene[un] = toto2
#	if U.ps_fictive(clasamex[signe][un]):
#		    print un, "inversions Etape4 ", clasamex[signe][un][0], cand[un], toto
    #Renvoie tous les PS homogenes de signe oppose
    return homogene

#------------------------------------------------------------------------------
def inversions_Etape5(clasamex, liste, signe, d, ids, minPS):
    """ Make homogenous sub-groups with PS in liste
    liste = liste des PS compatibles en orientation opposee de la PS 'un' dans Etape4"""
    #i = 0
    L = len(liste)
    #stop = 0
    sel = []
    for i in liste:
        
        tri = [i]
        for j in liste:
            
            cand_ok = candidatInvMemeSigne(clasamex[signe][i],
                                           clasamex[signe][j], d, ids)
#            if (i == 16 and j==175) or (j == 16 and i==175):
#                print i, j, cand_ok, clasamex[signe][i], clasamex[signe][j]
#	    if U.ps_fictive(clasamex[signe][i]) and U.ps_fictive(clasamex[signe][j]):
#		    print "inversions Etape5 ", i, clasamex[signe][i][0], j, clasamex[signe][j][0], cand_ok
            if (i != j) and cand_ok :
                tri.append(j)
            #j = j + 1

        #ne garde que les PS qui ont au moins un autre PS compatible
        #ne garde que les groupes qui ne sont pas deja formes
        #U.inde_groupsMinPS(sel, tri, minPS)
        U.inde_groups(sel, tri)
        if len(tri)  ==  L :
            break
    sel = U.filter_list(sel)   
    
    #verify that everybody is compatible with each other inside each 'group  of compatible PS'
    #if not keep the largest subgroupintersection 
    nsel = []
    for g in sel:
        t=[]
        for i, PS in enumerate(g):
            for j, otherPS in enumerate(g):
                if j>i:
                    cand_ok = candidatInvMemeSigne(clasamex[signe][PS],
                                           clasamex[signe][otherPS], d, ids)
                    if cand_ok:                                           
                        t.append((PS, otherPS))
        G=nx.Graph()
        G.add_edges_from(t)
        fused = nx.find_cliques_recursive(G)
        if fused:
            #if len(g) != len(max(fused, key=len)):
            #    print 'avant - fused',len(g), len(max(fused, key=len))
#            if 16 in g and 175 in g:
#                print "ET LA ???", g, max(fused, key=len)
            if len(max(fused, key=len)) > minPS:
                nsel.append(list(max(fused, key=len)))
              

    return nsel

#------------------------------------------------------------------------------
def inversions_Etape6(classorix, homo, groupes, signe, summary, minPS):
    """Keep only links between PS (-, -) and (+, +) present in intermediate
    file"""
    if groupes:
	#pour chaque groupe de PS de meme signe homogene
        for g in groupes:
            if g:
                g.sort()
		#pour chaque PS dans un groupe homogene
            new = []
            for i in g:
		#Regarde si les PS copines signe opposes sont partagees par d'autres PS de meme signe
#		print "Etape6 homo.keys ", homo.keys()
                if i in homo.keys() :
                    #print "copains oppose :", homo_signe[i]
                    for l in homo[i]:
                        common(l, new)
                else:
                    print "no opposite partners"
                    #new = homo[i]
#		if U.ps_fictive(classorix[signe][i]):
#			print "Etape6 fictive ", i, classorix[signe][i], new, signe

            #Remove small groups
            new2=[]
            for s in new:
                if len(s)>=minPS:
                    new2.append(s)

            for s in new2:
                #print "-----> SSS", s
                if s:
                    if signe  ==  "++":
                        if not sumCom(summary, s, g):
                            summary.append([s, g])
                        else:
                            if not sumCom(summary, g, s):
                                summary.append([g, s])
#    else:
#        print "alone"



#------------------------------------------------------------------------------
def inversions_Etape6_v2(classorix, homo, groupes, signe, summary, minPS):
    """Keep only links between PS (-, -) and (+, +) present in intermediate
    file"""
    #print "GROUPS", groupes
    if groupes:
	#pour chaque groupe de PS de meme signe homogene
        for g in groupes:
            if g:
                g.sort()
		#pour chaque PS dans un groupe homogene
            newD = {}
            for ind, i in enumerate(g):
		#Regarde si les PS copines signe opposes sont partagees par d'autres PS de meme signe
                if i in homo.keys() :
                    li = map(list, map(set, homo[i]))
                    #print "LI", li
                    newD[ind] = li
                else:
                    print "no opposite partners"
            new2 = reduce(lambda x, y: U.intersectall(x, y, minPS), newD.values())  #this is like below but with a reduce to limite the number of combinations
            
            
            
#            #the first time i used this approach, i used LIST(comb) which is memory Horrible !!!!
#            comb = itertools.product(*newD.values()) #combinations of groups, it is a generator
#            #not so good to do the comb because if i have 100PS in the groupe g and each PS is compatible with 4 opposite groupes
#            #it makes 4^100 possible combinations !! It is impossible to compute !!!
#            #intersections = map(list, [set.intersection(*x) for x in comb])
#            grp=[]             #list of groups that have an intersection of len>=minPS for each PS
#            for elt in comb:
#                tmp = list(set(elt[0]).intersection(*elt))
#                if len(tmp)>=minPS:
#                    if tmp not in grp:
#                        grp.append(tmp)
#            new2 = list(set(tuple(i) for i in grp)) #remove duplicates
                    
            

                            
#this is version 1 with the bug                            
#            flattened  = [val for sublist in newD.values() for val in sublist] #flattened list: put all lists of lists in the same list
#            #the flattened list is a mistake because i lose the info of which PS is compatible with each group !!
#            skt = list(set(tuple(i) for i in flattened)) #remove duplicates
#            toRem = []
#            for ieme, i in enumerate(skt):
#            	for jeme, j in enumerate(skt):
#            		if ieme < jeme:
#            			l = [len(i), len(j)]
#            			if l[0] != l[1]:
#            				#set([i, j][np.argmin(l)]).issubset(set([i, j][np.argmax(l)])), i, j, l
#            				if set([i, j][np.argmin(l)]).issubset(set([i, j][np.argmax(l)])):
#            					toRem.append([i, j][np.argmax(l)])
#            rien = [skt.remove(x) for x in list(set(tuple(i) for i in toRem))] #remove from skt the superlists
#            #Remove small groups
#            #new2=[x for x in new if len(x) >= minPS]
#            new2=[x for x in skt if len(x) >= minPS]
#            #for s in new:
#            #    if len(s)>=minPS:
#            #        new2.append(s)

            for s in new2:
                if s:
                    if signe  ==  "++":
                        if not sumCom(summary, s, g):
                            summary.append([s, g])
                    else:
                        if not sumCom(summary, g, s):
                            summary.append([g, s])
#    else:
#        print "alone"

#------------------------------------------------------------------------------
def inversions_fusion(params, clasamex, liste, chrox, nb, d, minPS, n):
    """Fuses inversions that are split due to  2-junctions constraints """
    #1: get everybodies coordinates
    
    dikiExtr = {}
    for ieme, i in enumerate(liste):
        deb, fin, maxlm, minrm, medSVm = extremites(clasamex, i[0], "--")
        start, end, maxlp, minrp, medSVp = extremites(clasamex, i[1], "++")
        left_in = max(deb, fin)
        right_in = min(start, end)
        left_out = maxlm
        right_out = minrp
        dikiExtr[ieme] = {}
        dikiExtr[ieme]['left_in'] = left_in
        dikiExtr[ieme]['right_in'] = right_in
        dikiExtr[ieme]['left_out'] = left_out
        dikiExtr[ieme]['right_out'] = right_out
        
        
    #2: identify compatible invs, creat all 2 by 2 comparisons    
    toFuse = []    
    for ieme, i in enumerate(liste):
        left_in = dikiExtr[ieme]['left_in']
        right_in = dikiExtr[ieme]['right_in']
        left_out = dikiExtr[ieme]['left_out']
        right_out = dikiExtr[ieme]['right_out']
        for jeme, j in enumerate(liste):
            left_inj = dikiExtr[jeme]['left_in']
            right_inj = dikiExtr[jeme]['right_in']
            left_outj = dikiExtr[jeme]['left_out']
            right_outj = dikiExtr[jeme]['right_out']

            overlapLeft = max(left_out, left_outj)<min(left_in, left_inj) #check overlap 
            overlapRight = max(right_inj, right_in)<min(right_outj, right_out)

            c1 = abs(left_inj-left_in)<d #check proximity of coordinates
            c2 = abs(right_inj-right_in)<d
            c3 = abs(left_outj-left_out)<d
            c4 = abs(right_outj-right_out)<d
            
            p1 = max(right_inj, right_in) < min(right_outj, right_out)
            p2 = max(left_outj, left_out) < min(left_inj, left_in)
            #print "OVERLAP", overlapLeft, overlapRight, max(right_inj, right_in) < min(right_outj, right_out), max(left_outj, left_out) < min(left_inj, left_in)
            if ieme < jeme: #the compatibility test
                if p1 and p2 and \
                c1 and c2 and c3 and c4 and \
                overlapLeft and overlapRight:
                    toFuse.append([ieme, jeme])
    #3: fuse compatible invs as maximal cliques
    G=nx.Graph()
    G.add_edges_from(toFuse)
    fusedG = nx.find_cliques_recursive(G)
    if fusedG:
        f = [val for sublist in fusedG for val in sublist] #flatten list
        singletonsINV = list(set(range(0,len(liste))).difference(set(f))) #search singletons (INV compatible with no other one)
        fusedG = fusedG + [[x] for x in singletonsINV] #add singletons
#        if n==2:
#            print "liste", liste
#            print "toFuse", toFuse
#            print "fusedG", fusedG
        fusedListe = []
        for elt in fusedG:
            mm = []
            pp = []
            for SV in elt:
                mm.extend(liste[SV][0])
                pp.extend(liste[SV][1])
                mm = list(set(mm))
                pp = list(set(pp))
            fusedListe.append([mm, pp])
#        print "LONGUEURS", len(liste), len(fusedListe)
        return fusedListe, True
    else:
        return liste, False
            
  

#------------------------------------------------------------------------------
def inversions_Affiche(params, clasamex, liste, chrox, nb, d, minPS):
    """Write results to outpufiles """
    out = open(params["out"]+"_inversions_byRP.csv", "a")
    outp = open(params["out"]+"_inversions_bySV.csv", "a")
    manip = os.path.split(params["in"])[1]
    #print "Affiche Inversion ", liste
    #ff = 1
    for i in liste:
        
        n1 = len(i[0])
        n2 = len(i[1])
        nbPS = n1 + n2
        #print nbPS, n1, n2
        if  i[0] and i[1] and nbPS >= minPS :
            
            nb += 1
            deb, fin, maxlm, minrm, medSVm = extremites(clasamex, i[0], "--")
            start, end, maxlp, minrp, medSVp = extremites(clasamex, i[1], "++")
            left_in = max(deb, fin)
            right_in = min(start, end)
            left_out = maxlm
            right_out = minrp

            # fin = max right des ps --
            # start = min left des ps ++
            
            minsize = min(abs(minrp -maxlm), abs(fin - start))
            
            # minrp = min right ++
            # maxlm = max left --

            maxsize = max(abs(minrp -maxlm), abs(fin - start))
            
            pvalue = U.testStatSciPy(n1, n2)
            #print pvalue
#            field=("%s","%s","%d","%d","%d","%d","%d","%d","%d","%d","%d","%d","%d","%d","%d","%s","%f","%s","s;","%s","%d","%s","%s","%d","d")
#            n = 0
#            for toto in (out_f, str(chrox), nb, nbPS, n1, n2, left_out, left_in, fin-deb, "NA", right_out, right_in, end-start, "NA", minsize, maxsize, pvalue):
#                print n, field[n], toto, type(toto)
#                n = n + 1
            U.write_bySV(outp, (manip, chrox, nb, nbPS, len(i[0]), len(i[1]), 
                                left_out, left_in, abs(left_in - left_out), 
                                "NA", right_out, right_in, abs(right_out - right_in),
                                "NA", minsize, maxsize, pvalue, "NA"))
            for s in i[0]:
                si = clasamex["--"][s]
                #if ff ==1:
                #    print "ff", s, si[0]
                U.write_byPS(out, (manip, chrox, nb, nbPS, len(i[0]), len(i[1]), 
                                left_out, left_in, abs(left_in - left_out), 
                                "NA", right_out, right_in, abs(right_out - right_in),
                                "NA", minsize, maxsize, pvalue, "NA", si[0], si[1],
                                si[2], si[3], si[4], si[5], si[6], si[7]))
#		if U.ps_fictive(si):summary
#			print "Affiche Inversion ",si[0]
            for s in i[1]:
                si = clasamex["++"][s]
                U.write_byPS(out, (manip, chrox, nb, nbPS, len(i[0]), len(i[1]), 
                                left_out, left_in, abs(left_in - left_out), 
                                "NA", right_out, right_in, abs(right_out - right_in),
                                "NA", minsize, maxsize, pvalue, "NA", si[0], si[1],
                                si[2], si[3], si[4], si[5], si[6], si[7]))
#		if U.ps_fictive(si):
#			print "Affiche Inversion ",si[0]
            #ff +=1
    out.close()
    outp.close()
    return nb
 


#------------------------------------------------------------------------------

def runStatsInv(params, stats, chrDicos):
    n = int(params["n"])

# START ALEX PARAMETERS #######################################################
    list_chr_names = "SEP".join(params["range"])
    list_chr_length = "SEP".join([ str(chrDicos[x]) for x in params["range"] ])
    #list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    #list_chr_names = "SEP".join(chrDicos.keys())
    
    
    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + n * stats["mad"]
    

    print "Statistics on Inversions : Dps = ", stats["ninv"]
    print "False Discovery Rate (fdr) = ", params["fdr"]
    
    toto, msg, pval_seuil_inv = Ualex.functstats("INV", stats, stats["ninv"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_inversions_bySV.csv", 
                                 params["in"], params["in"]+".dist.table", 
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])
                                 
    print msg  
    return pval_seuil_inv

#-------------------------------------------------------------------------
def runDetectionInv(params, stats, chrDicos, list_chr_real):
# START ALEX PARAMETERS #######################################################

    list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    list_chr_names = "SEP".join(chrDicos.keys())
    
    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + float(params["n"]) * stats["mad"]

    #Determine le nb minimum de ps pour chaque type de SV
    detectionFile = "ClusterLimitSize"
    #print params["in"]+"_PScandidate_intra.csv", "INV", list_chr_real
     
    #print n
    
     
    ps_min, infoMsg = Ualex.define_MinPS("INV", stats, params["nsv"], 
        list_chr_length, list_chr_names, detectionFile, params["in"], 
        params["in"]+".dist.table", sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"], list_chr_real)

# END ALEX PARAMETERS #######################################################

    dicQual = {}
    D = int(stats["median"]) + int(params["n"])*float(stats["mad"])
    IDS = 2*D
    #print "IDS", IDS
    #IDS = 1000*D #TODO REMOVE
    U.create_files("inversions", params["out"])
#    prefix = stats["chromosome_prefix"]
#TO BE REMOVED FOR TEST RAD27 
    #TODO
    #ps_min=2
    print "MINPS for inversions:", ps_min
    with open(params["out"]+".INV.report.out","a") as fileo:
        fileo.write("\n----------- \tDetection of Inversion\n")
        fileo.write("\t\tInsert Size Consistency (isc) = %d\n" % (IDS))
        fileo.write("\t\tMaximum distance between PS extremities (d) = %.2f\n"\
        % (D))
        fileo.write("\t\tMinimum number of PS to define an Inversion = %d\n\n"\
        % (ps_min))
        
    for chrx in params["range"]:
        if not os.path.isfile(params["in"]+"_"+chrx):
            continue
        
        classorix = {"++":[], "--":[], "-+":[]}
        ReadFilesBAM(params, chrx, classorix, dicQual)
        #print(dicQual)
        
#        for pair, PSS in classorix.iteritems():
#        #print "pair", pair
#            clasx = list(set([tuple(i) for i in PSS]))
#            classorix[pair]=clasx
        
        #print "classorix", classorix['--']
        print "Processing", chrx, "- Nb of discordant PS:", str(len(classorix['--'])+len(classorix['++'])), \
        "-", datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        #TODO: TOREMOVE
        if int(stats["median"])<1000:
            ps_min = 1
        #ps_min = 1
#        print "D1", D
        #Pour chaque moins, moins recuperer copains plus, plus. et inversement
        candimoins, candiplus = inversions_Etape1et2(classorix, D, IDS, 
                                                     stats["ps_type"], ps_min)
        #
#        print "D2", D
        #Faire des sous-groupes homogenes avec les plus, plus copains d'un moins, moins
        homo_moins = inversions_Etape4(classorix, candimoins, "--", D, IDS, ps_min)
        #print homo_moins
#        print "D3", D
        #Faire des sous-groupes homogenes avec les moins, moins copains d'un plus, plus
        homo_plus = inversions_Etape4(classorix, candiplus, "++", D, IDS, ps_min)
        
#        import pickle
#                
#        ids = zip(*classorix['--'])[0]
#        PS1 = 'HWI-D00473:127:C5T3JACXX:2:2215:11082:100048'
#        PS2 = 'HWI-D00473:127:C5T3JACXX:3:2310:15747:5973'
#        
#        
#        print ids.index(PS1)
#        print ids.index(PS2)
        


#        pickle.dump(homo_moins, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/homo_moins.p", "wb"))
#        pickle.dump(homo_plus, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/homo_plus.p", "wb"))
#        print "D4", D
        #Avec les -,- (avec copains +, +) Faire des groupes homogenes de -,-
        #cherche les sousgroupes de -- compatbles entre eux
         
        groupe_moins = inversions_Etape5(classorix, homo_moins.keys(), "--", D,
                                         IDS, ps_min)
#        print "D5", D
#        for g in groupe_moins:
#            if 175 in g and 16 in g:
#                print "GGG pas bon", g
        #Avec les +,+ (avec copains -, -) Faire des groupes homogenes de +,+
        #cherche les sousgroupes de ++ compatbles entre eux
        groupe_plus = inversions_Etape5(classorix, homo_plus.keys(), "++", D,
                                        IDS, ps_min)
#        pickle.dump(groupe_moins, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/groupe_moins.p", "wb"))
#        pickle.dump(groupe_plus, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/groupe_plus.p", "wb"))        
#        print "D6", D
##########################################
#        out = open(params["out"]+"_inversions_byRP_GROUPESmoins.csv", "w")
#        nb=0
#        for i in groupe_moins:
#            n1 = len(i)
#            n2 = len(i)
#            nbPS = n1 + n2
#            nb+=1
#            if  i:
#                nb+=1
#                deb, fin, maxlm, minrm, medSVm = 0,0,0,0,0
#                start, end, maxlp, minrp, medSVp = 0,0,0,0,0
#                left_in = max(deb, fin)
#                right_in = min(start, end)
#                left_out = maxlm
#                right_out = minrp
#                minsize = min(abs(minrp -maxlm), abs(fin - start))
#                maxsize = max(abs(minrp -maxlm), abs(fin - start))
#                
#                pvalue = 0
#                for s in i:
#                    si = classorix["--"][s]
#                    #if ff ==1:
#                    #    print "ff", s, si[0]
#                    U.write_byPS(out, ("test", chrx, nb, nbPS, len(i), len(i), 
#                                    left_out, left_in, abs(left_in - left_out), 
#                                    "NA", right_out, right_in, abs(right_out - right_in),
#                                    "NA", minsize, maxsize, pvalue, "NA", si[0], si[1],
#                                    si[2], si[3], si[4], si[5], si[6], si[7]))
#        out.close()
#
#        out = open(params["out"]+"_inversions_byRP_GROUPESplus.csv", "w")
#        
#        nb=0
#        for i in groupe_plus:
#            n1 = len(i)
#            n2 = len(i)
#            nbPS = n1 + n2
#            nb+=1
#            if  i:
#                deb, fin, maxlm, minrm, medSVm = 0,0,0,0,0
#                start, end, maxlp, minrp, medSVp = 0,0,0,0,0
#                left_in = max(deb, fin)
#                right_in = min(start, end)
#                left_out = maxlm
#                right_out = minrp
#                minsize = min(abs(minrp -maxlm), abs(fin - start))
#                maxsize = max(abs(minrp -maxlm), abs(fin - start))
#                
#                pvalue = 0
#                for s in i:
#                    si = classorix["++"][s]
#                    #if ff ==1:
#                    #    print "ff", s, si[0]
#                    U.write_byPS(out, ("test", chrx, nb, nbPS, len(i), len(i), 
#                                    left_out, left_in, abs(left_in - left_out), 
#                                    "NA", right_out, right_in, abs(right_out - right_in),
#                                    "NA", minsize, maxsize, pvalue, "NA", si[0], si[1],
#                                    si[2], si[3], si[4], si[5], si[6], si[7]))
#        out.close()
#        print "DOOOONE"
    

##############################

        #Sortir tous les sous-groupes homogenes (-, -) et (+, +) avec les liens presents dans le fichier intermediaire
        resume = []
        inversions_Etape6_v2(classorix, homo_moins, groupe_moins, "--", resume, ps_min)
        inversions_Etape6_v2(classorix, homo_plus, groupe_plus, "++", resume, ps_min)

        
        
        resumeClean = rmSubInv(resume)
        
        
#        pickle.dump(params, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/params.p", "wb"))
#        pickle.dump(classorix, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/classorix.p", "wb"))        
#        pickle.dump(resumeClean, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/resumeClean.p", "wb"))
#        pickle.dump(chrx, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/chrx.p", "wb"))
#        pickle.dump(D, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/D.p", "wb"))
#        pickle.dump(ps_min, open("/home/alex/NGS_data/ulysses3.0/WT/inv/splitDetection/testChr1/ps_min.p", "wb"))
#        print "PICKLE DONE"
        
        fusionHappend = True
        n = 1
        
#        resumeClean, fusionHappend = inversions_fusion(params, classorix, resumeClean, chrx, 0, D, ps_min, n)
#        print "fusion2"
#        resumeClean, fusionHappend = inversions_fusion(params, classorix, resumeClean, chrx, 0, D, ps_min, n)

        while fusionHappend:
            resumeClean, fusionHappend = inversions_fusion(params, classorix, resumeClean, chrx, 0, D, ps_min, n)
            n+=1

        
        
        #print "resume", resumeClean
        #Affichage
        #print "$$$$$$   AFFICHAGE "
        
        nbinv = inversions_Affiche(params, classorix, resumeClean, chrx, 0, D, ps_min)
        
    #add qualities
    U.addMeanSVQuality(params["out"]+"_inversions_bySV.csv", params["out"]+"_inversions_byRP.csv", dicQual)    
    
    print "Detection done\n"
    with open(params["out"]+".INV.report.out","a") as fileo:
        fileo.write("Detection processed\n")
        fileo.write("results written to "+params["out"]+"_inversions_by[PS/SV].csv\n")    
    toto, msg, pval_seuil_inv = Ualex.functstats("INV", stats, stats["ninv"], params["nsv"], 
                                 list_chr_length, list_chr_names, 
                                 params["out"]+"_inversions_bySV.csv", 
                                 params["in"], params["in"]+".dist.table", 
                                 sval, params["fdr"], params["out"], stats["rl"], 
                                 params["n"])
    with open(params["out"]+".INV.report.out","a") as fileo:
        fileo.write(msg+"\n")                             
    print msg
    return pval_seuil_inv

#------------------------------------------------------------------------------
def launch(params, stats, chrDicos):
    """Detection of inversions sys.argv[1] = Name of parameter file
    """
    if params["only_stats"]:
        pval_seuil_inv = runStatsInv(params, stats, chrDicos)
    else:
        pval_seuil_inv = runDetectionInv(params, stats, chrDicos, params["range"])
    return pval_seuil_inv
