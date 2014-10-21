# -*- coding: utf-8 -*-
__date__      = "2011/02/21"
#__doc__       = """ texte doc """

import os, sys
import datetime
import itertools
now = datetime.datetime.now()


try:
    import pysam
except ImportError:
    print "\tError: Pysam is not installed. \n\
    Please see: https://github.com/pysam-developers/pysam\n\n"
    sys.exit()
try:
    from numpy import median
except ImportError:
    print "\tError: NumPy is not installed. \n\
    Please see: http://www.numpy.org/\n\n"
    sys.exit()

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
    2 dicos (un pour les -*- l autre pour les ++:
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
    
    Cherche les sous groupse de PS compatiblent entre elles en orientation opposee de chaque PS """
    homogene = {}
    oppose = "--"
    if signe  ==  "--":
        oppose = "++"
    for un in cand.keys():
        toto = inversions_Etape5(clasamex, cand[un], oppose, d, ids, minPS) #fait des sousgroups homogenes dans l orientation oppose a la PS 'un'
        if toto :
            toto2=[]
            for elt in toto:
                if len(elt)>=minPS:
                    toto2.append(elt)
            
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
#	    if U.ps_fictive(clasamex[signe][i]) and U.ps_fictive(clasamex[signe][j]):
#		    print "inversions Etape5 ", i, clasamex[signe][i][0], j, clasamex[signe][j][0], cand_ok
            if (i != j) and cand_ok :
                tri.append(j)
            #j = j + 1

        #ne garde que les PS qui ont au moins un autre PS compatible
        #ne garde que les groupes qui ne sont pas deja formes
	U.inde_groupsMinPS(sel, tri, minPS)
        if len(tri)  ==  L :
            break
        
                    
    return sel

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
    #print groupes
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
                    #for l in homo[i]:
                    #    common(l, new)
                    #print "homo[i]", homo[i]
                    li = map(set, homo[i])
                    new.append(li)
                else:
                    print "no opposite partners"
                    #new = homo[i]
#		if U.ps_fictive(classorix[signe][i]):
#			print "Etape6 fictive ", i, classorix[signe][i], new, signe
            comb=list(itertools.product(*new))
            #print comb
            intersections = map(list, [set.intersection(*x) for x in comb])
            #print intersections

            #Remove small groups
            new2=[x for x in intersections if len(x) >= minPS]
            #for s in new:
            #    if len(s)>=minPS:
            #        new2.append(s)

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
def inversions_Affiche(params, clasamex, liste, chrox, nb, d, minPS):
    """Write results to outpufiles """
    out = open(params["out"]+"_inversions_byPS.csv", "a")
    outp = open(params["out"]+"_inversions_bySV.csv", "a")
    manip = os.path.split(params["in"])[1]
    #print "Affiche Inversion ", liste

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
            
            minsize = abs(fin - start)
            
            # minrp = min right ++
            # maxlm = max left --

            maxsize = minrp -maxlm
            
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
        
###################################################################################################
#### TEMPORAIRE A CAUSE DES CHIMERES (BAM file: flag pas ok, il peut y avoir 2 primary reads)
        for pair, PSS in classorix.iteritems():
        #print "pair", pair
            clasx = list(set([tuple(i) for i in PSS]))
            classorix[pair]=clasx
######################################################################################################
        print "Processing", chrx, "- Nb of discordant PS:", str(len(classorix['--'])+len(classorix['++'])), \
        "-", now.strftime("%Y-%m-%d %H:%M")
        #TODO: TOREMOVE
        if int(stats["median"])<1000:
            ps_min = 1
        #ps_min = 1
        #Pour chaque moins, moins recuperer copains plus, plus. et inversement
        candimoins, candiplus = inversions_Etape1et2(classorix, D, IDS, 
                                                     stats["ps_type"], ps_min)
        
        #Faire des sous-groupes homogenes avec les plus, plus copains d'un moins, moins
        homo_moins = inversions_Etape4(classorix, candimoins, "--", D, IDS, ps_min)
        
        
        #Faire des sous-groupes homogenes avec les moins, moins copains d'un plus, plus
        homo_plus = inversions_Etape4(classorix, candiplus, "++", D, IDS, ps_min)
        
        #Avec les -,- (avec copains +, +) Faire des groupes homogenes de -,-
        #cherche les sousgroupes de -- compatbles entre eux
        groupe_moins = inversions_Etape5(classorix, homo_moins.keys(), "--", D,
                                         IDS, ps_min)
        #Avec les +,+ (avec copains -, -) Faire des groupes homogenes de +,+
        #cherche les sousgroupes de ++ compatbles entre eux
        groupe_plus = inversions_Etape5(classorix, homo_plus.keys(), "++", D,
                                        IDS, ps_min)
        
        #Sortir tous les sous-groupes homogenes (-, -) et (+, +) avec les liens presents dans le fichier intermediaire
        resume = []
        inversions_Etape6_v2(classorix, homo_moins, groupe_moins, "--", resume, ps_min)
        inversions_Etape6_v2(classorix, homo_plus, groupe_plus, "++", resume, ps_min)
        
        
        resumeClean = rmSubInv(resume)
        
        #print "resume", resumeClean
        #Affichage
        #print "$$$$$$   AFFICHAGE "
        
        nbinv = inversions_Affiche(params, classorix, resumeClean, chrx, 0, D, ps_min)
        
    #add qualities
    U.addMeanSVQuality(params["out"]+"_inversions_bySV.csv", params["out"]+"_inversions_byPS.csv", dicQual)    
    
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
def launch(paramfile, onlyStatPerform, list_chr_real):
    """Detection of inversions sys.argv[1] = Name of parameter file
    """

    if os.path.isfile(paramfile):
        params, stats, chrDicos = U.prepare_detection("inversions", 
                                                      paramfile, "NA")
        #print "ggg", list_chr_real, len(list_chr_real)
        #if ''.join(list_chr_real)!='all':
        
        #if len(list_chr_real)==3:
#            if list_chr_real[0]=='a' and list_chr_real[1]=='l' and list_chr_real[2]=='l':
 #               list_chr_real = [ stats["chromosome_prefix"] + str(x) for x in list_chr_real ]
        
        if onlyStatPerform:
            pval_seuil_inv = runStatsInv(params, stats, chrDicos)
        else:
            pval_seuil_inv = runDetectionInv(params, stats, chrDicos, list_chr_real)  
        return pval_seuil_inv
    else:
        print "Error :", paramfile, "doesn't exist"
        
#------------------------------------------------------------------------------
if (__name__)  ==  "__main__":
    paramfile = U.parser("inversions")
    launch(paramfile)
