# -*- coding: utf-8 -*-
""" texte doc """



import Ulysse_utils as U
import os, sys, random
try:
    import pysam
except ImportError:
    print "\tError: Pysam is not installed. \n\
    Please see: https://github.com/pysam-developers/pysam\n\n"
    sys.exit()
import operator
import datetime
try:
    import networkx as nx
except ImportError:
    print "\tError: networkx package in not installed. \n\
    Please see: http://networkx.lanl.gov/index.html"



#--------------------------------------------------------------------------
def ReadFilesBAM(params, clasx, ps_type):
    """
    Read BAM Files of discordant PS
    """
#    names = []
    dicQual = {}
    with pysam.Samfile(params["in"]+"_xdiff", 'rb') as bam:


    #Retrieve PS Only read1 is written after BAM filtering
        for read1 in bam :
            if read1.is_read1:
                if read1.qname in dicQual:
                    dicQual[read1.qname].append(read1.mapq)
                else:
                    dicQual[read1.qname] = [read1.mapq]
    
                chr1=bam.getrname(read1.tid)
                chr2=bam.getrname(read1.rnext)
                ori1, chr1, pos1, ori2, chr2, pos2 = U.getCoord(read1, chr1, chr2)
                clasx.append([read1.qname, ori1, chr1, pos1, ori2, chr2, pos2, abs(read1.tlen)])

    return dicQual

#------------------------------------------------------------------------------
#TEMPORAIRE A VIRER (lecture fichier ulysse)
def enumerate_X(xsome, names):
    """ Store original names of xsome in list names
    """
    if xsome.isdigit():
        xsome = int(xsome)
    if xsome not in names:
        names.append(xsome)
#------------------------------------------------------------------------------
#TEMPORAIRE A VIRER (lecture fichier ulysse)
def get_X_number(names):
    """ Attribute numbers to xsome in list names
        Return dict ame2numb
    """
    names.sort()
    print "Original chromosome names : ", names
    return dict(zip(map(str, names), map(str, range(1, len(names) + 1))))

#------------------------------------------------------------------------------
#TEMPORAIRE A VIRER (lecture fichier ulysse)
def get_X_original(pair, names):
    """ Retrieve the original xsome names from dict names
    """
    inv = dict(zip(names.values(), names.keys()))
    ch1 = pair.split("-")[0]
    ch2 = pair.split("-")[1]
    return inv[ch1] + "-" + inv[ch2]
#--------------------------------------------------------------------------
def orderPS(liste, names):
    """Pour que l'ordre d'apparition des chromosomes dans un PS donne soit
    toujours croissant. Uniquement pour PS impliquant deux chromosomes
    differents"""

    liste[2] = names[liste[2]]
    liste[5] = names[liste[5]]
    if int(liste[5]) < int(liste[2]):
        #inverse les chromosomes
        liste[2],  liste[5] = liste[5],  liste[2]
        #inverse l'orientation
        liste[1],  liste[4] = liste[4],  liste[1]
        #inverse les coordonnees
        liste[3],  liste[6] = liste[6],  liste[3]
#--------------------------------------------------------------------------
#TEMPORAIRE A VIRER (lecture fichier ulysse)
def ReadFilesPaired(fil,  clasx):
    """
    Read Files in format obtained after filtering of data formatted
    like Paired.uniquelymapped_noDuplicates.SV
    """
    names = []
    with open(fil, 'r') as gile:
    #Retrieve PS
        for line in gile :
            ps = line[:-1].split(";")
            clasx.append( ps[0:3] + [int(ps[3]),  ps[4],  ps[5],  int(ps[6])])
            enumerate_X(ps[2], names)
            enumerate_X(ps[5], names)

        #enumerate xsome for analysis
        number_X = get_X_number(names)
        #order reads in PS and give X number
        for ps in clasx :
            orderPS(ps, number_X)

    return number_X

#--------------------------------------------------------------------------
def Filter(clasx,  xsome, chrs):
    """Store PS according to chromosome pairs"""
    #print chrs
    for i in clasx:
        if i[2] in chrs and i[5] in chrs:
            p = i[2]+"-"+i[5]
            
            if p in xsome:
                xsome[p].append(i)
            else:
                xsome[p] = [i]

#--------------------------------------------------------------------------
def isSubtelo_Cap(xsome, centrom, subtelo, sv, interne):
    """Determine if coordinates are subtelomeric."""

    oksub = []
    okins = []
    bks = "-"
    bki = "-"
    for i in sv:
        lis = xsome[i]

        centro = U.Centro_Inter(centrom, lis[2], lis[1], lis[3],
                                           lis[5], lis[4], lis[6])
        #print "Avant Centro Inter ", lis[2], lis[1], lis[3], lis[5], lis[4], lis[6], centro
        n = False
        for cas in interne:
            if centro == cas[0]:
                n = True

########################################################################
#Uniquement pour les genomes dont on connait les coordonnees subetlo
                if subtelo:
                    nb = centro.index(cas[2])
                    if nb == 0:
                        sub = U.Subtelo(subtelo, lis[2], lis[3])
                    else:
                        sub = U.Subtelo(subtelo, lis[5], lis[6])
                    if sub == 1:
                        oksub.append(i)
                        bks = centro
########################################################################
                else:
                    oksub.append(i)
                    bks = centro

        if subtelo:
            is_subtelo1 = U.Subtelo(subtelo, lis[5], lis[6])
            is_subtelo2 = U.Subtelo(subtelo, lis[2], lis[3])
            #pour insertion
            if (not n) and (is_subtelo1 == 0) and (is_subtelo2 == 0):
                okins.append(i)
                bki = centro
        # ELIF BY ALEX #################################################
        # Sans ca, les NRT se classent dans INS lorsqu'il n y a pas
        # de fichier de subtelo ni de centro --> a faire valider par Ingrid
        elif not subtelo and centro == "NiNi":
            #print "ca a la air d etre ca"
            oksub.append(i)
            bks = centro
        ################################################################

        else:
            if not n:
#                okins.append(i)
#                bki = centro
                oksub.append(i)
                bks = centro

#        if U.ps_fictive(lis):
#            print "isSubtelo ",lis[0], n, centro

    return oksub, bks, okins, bki

#--------------------------------------------------------------------------
def homosameORI(xsome, liste, d, minPS, pairOppCoordListe):
    """Make homogeneous sub-groups with PS in a given orientation."""
    #print "====", liste
    i = 0
    L = len(liste)
#    print "homosame 0", [xsome[j][0] for j in liste]
    stop = 0
    sel = []
    while i < L and stop == 0:
        tri = [liste[i]]
        cai = xsome[liste[i]]
        Fictive = U.ps_fictive(cai)
        j = 0

        while j < L:
            caj = xsome[liste[j]]
            coord_left_ok = U.parameter_ok(cai[3], caj[3], d)
            coord_right_ok = U.parameter_ok(cai[6], caj[6], d)
            if (i != j) and  coord_left_ok and coord_right_ok  :
                tri.append(liste[j])
#            if U.ps_fictive(cai) and U.ps_fictive(caj):
#                print "homosame 1",liste[i], cai[0], liste[j], caj[0], tri, coord_left_ok, coord_right_ok
            j = j+1

        #if len(tri) > 1:
            #ne garde que les groupes qui ne sont pas deja formes
        left = 1
        for k in sel:
            if len( [l for l in tri if l not in k] ) == 0:
                left = 0
        if left > 0:
            sel.append(tri)

        if len(tri) == L:
            stop = 1

        i = i+1

    #print "homosame sel ", sel

    homo = []
#    print "Sel 1",sel
    if sel:
        for i in sel:
            if len(i)>= minPS:
                left = 1
                for j in sel:
                    if (i != j) and  len( [l for l in i if l not in j] ) == 0:
                        left = 0
                if left > 0:
                    i.sort()
                    homo.append(i)
#        if Fictive:
#            print "homosame 2 toto ", homo
#
#
    #print "homosame homo", homo
    return homo



#--------------------------------------------------------------------------
def homosame(xsome, liste, d, minPS, pairOppCoordListe):
    """Make homogeneous sub-groups with PS in a given orientation."""

    #get the PS IDs
    liste2 = map(U.conv, pairOppCoordListe)
   
    #sort the list by 1) coord left and 2) right coord
    liste2 = sorted(liste2, key=operator.itemgetter(1, 2))
    

    seli =[]
    sel2 =[]
    #for each PS, search for the compatible ones in the same orientation
    for i, psi in enumerate(liste2):
        #print "+++", i, psi
        #a ps is obviously compatible with at least itself...
        tri = [psi[0]]
        tri2 = [psi]
        #PS are sorted primarily by left coord and secondarily by right coord
        #Therefore, for a given PS, compatible PS are the ones right before
        #(to the left) and right after (to the right). We search for compatible
        #PS to the left (and right) until we find an uncompatible PS on the
        #left coord.
        # THIS SHOULD BE DONE WITH BINARY SEARCH !!!
        #initialize for left and right search
        i1 = i+0
        i2 = i+0

        #left search
        while 1:
            if i1<len(liste2)-1:
                i1=i1+1
            else:
                break
            coord_left_ok = (abs(psi[1]-liste2[i1][1]) <= d)
            if coord_left_ok:
                coord_right_ok = (abs(psi[2]-liste2[i1][2]) <= d)
            else:
                break
            if coord_left_ok and coord_right_ok and psi[0] != liste2[i1][0]:
               #tri.append(liste[i1])
               tri.append(liste2[i1][0])
               tri2.append(liste2[i1])

        #right search
        while 1:
            if i2>0:
                i2=i2-1
            else:
                break
            coord_left_ok = (abs(psi[1]-liste2[i2][1]) <= d)
            if coord_left_ok:
                coord_right_ok = (abs(psi[2]-liste2[i2][2]) <= d)
            else:
                break
            if coord_left_ok and coord_right_ok and psi[0] != liste2[i2][0]:
               #tri.append(liste[i2])
               tri.append(liste2[i2][0])
               tri2.append(liste2[i2])

        #Sort PS IDs nicely
        tri2.sort(key=lambda x: x[0])
        tri.sort()

        #only append non-already found groups that are big enough
        if ( tri not in seli ) and len(tri) >= minPS:
            seli.append(tri)

        #only append non-already found groups that are big enough
        if ( tri2 not in sel2 ) and len([x[0] for x in tri2]) >= minPS:
            sel2.append(tri2)

    return seli



def homoBinaryRight(psi, liste, d, listeAll, ieme):
    """   """
    orderAll = listeAll
    lenorder = len(liste) #longueur de la liste en cours
    if lenorder== 1:
        return ieme
    iddic = lenorder/2 #index de la valeur "centrale" de la liste
    psj = liste[iddic] #valeur de la valeure "centrale"
    if ieme=="NA": #cas ou c est une nouvelle liste (pas encore recurence)
        ieme = iddic
    #Si la valeur cherchee est plus petite que la valeur a la position actuelle
    #dans la liste, il faut prendre la sous partie gauche de la liste
    if psi+d <= psj:
        #calculer l indice de correction pour le calcule de la position de la
        #valeur mediane dans la future liste par rapport a la toute premiere
        #liste (orderAll)
        #Si on va a gauche, la correction ne depend que de la parite de la
        #longueur de la nouvelle liste
        if len(liste[:iddic])%2 == 0:
            correc = 0
        else:
            correc = 1
        return homoBinaryRight(psi, liste[:iddic], d, orderAll, ieme-iddic/2-correc)
    #Si la valeur cherchee est plus grande que la valeur a la position actuelle
    #dans la liste il faut prendre la sous partie doite de la liste
    elif psi+d >= psj:
        #calculer l indice de correction pour le calcule de la position de la
        #valeur mediane dans la future liste par rapport a la toute premiere
        #liste (orderAll)
        #Si on va a droite, la correction depend de la parite de la
        #longueur de la nouvelle liste mais aussi de la parite de la longueur
        # de la liste actuelle
        if lenorder%2 ==1 and len(liste[iddic:])%2 == 0:
            correc = 1
        else:
            correc = 0
        return homoBinaryRight(psi, liste[iddic:], d, orderAll, ieme+iddic/2+correc)


def homoBinaryLeft(psi, liste, d, listeAll, ieme):
    """   """
    orderAll = listeAll
    lenorder = len(liste) #longueur de la liste en cours
    #print liste, ieme, psi
    if lenorder < 3:
        return ieme
    iddic = lenorder/2 #index de la valeur "centrale" de la liste
    psj = liste[iddic] #valeur de la valeure "centrale"
    if ieme=="NA": #cas ou c est une nouvelle liste (pas encore recurence)
        ieme = iddic
    #Si la valeur cherchee est plus petite que la valeur a la position actuelle
    #dans la liste, il faut prendre la sous partie gauche de la liste
    if psj >  psi-d :
        #calculer l indice de correction pour le calcule de la position de la
        #valeur mediane dans la future liste par rapport a la toute premiere
        #liste (orderAll)
        #Si on va a gauche, la correction ne depend que de la parite de la
        #longueur de la nouvelle liste
        if len(liste[:iddic])%2 == 0:
            correc = 0
        else:
            correc = 1
        return homoBinaryLeft(psi, liste[:iddic], d, orderAll, ieme-iddic/2-correc)
    #Si la valeur cherchee est plus grande que la valeur a la position actuelle
    #dans la liste il faut prendre la sous partie doite de la liste
    elif psj <=  psi-d:
        #calculer l indice de correction pour le calcule de la position de la
        #valeur mediane dans la future liste par rapport a la toute premiere
        #liste (orderAll)
        #Si on va a droite, la correction depend de la parite de la
        #longueur de la nouvelle liste mais aussi de la parite de la longueur
        # de la liste actuelle
        if lenorder%2 ==1 and len(liste[iddic:])%2 == 0:
            correc = 1
        else:
            correc = 0
        return homoBinaryLeft(psi, liste[iddic:], d, orderAll, ieme+iddic/2+correc)

#--------------------------------------------------------------------------
def homosameBinary(xsome, liste, d, minPS, pairOppCoordListe):
    """Make homogeneous sub-groups with PS in a given orientation."""

    #print "1", datetime.datetime.now()
    #get the PS IDs
    liste2 = map(U.conv, pairOppCoordListe)
    #print "gggg", liste2
    #print minPS
    #liste_light = [x[0] for x in liste]
    
    sel2 = []

    #sort the list by 1) coord left and 2) right coord
    liste2 = sorted(liste2, key=operator.itemgetter(1, 2))
    
    idy = [x[0] for x in liste2]
    cleft = [x[1] for x in liste2]
    cright = [x[2] for x in liste2]
    
    #d = int(d)
    
    for i, psi in enumerate(idy):
        #a ps is obviously compatible with at least itself...
        tri2 = []
        #First find the position where to start search
        j = homoBinaryRight(cleft[i], cleft[i:], d, cleft, 'NA')
        #print i, j, cleft[i], cleft[i:], d
        if j == 'NA':
            j = i
        else:
            j = i + j
        if j>i:
            while j>i:
                tri2.append([idy[j], cleft[j], cright[j]])
                j = j - 1
        j = homoBinaryLeft(cleft[i], cleft[:i+1], d, cleft, 'NA')
        if j == 'NA':
            j = i
        if j<=i:
            while j<=i:
                #print i, j, idy, cleft, cright
                tri2.append([idy[j], cleft[j], cright[j]])
                j = j + 1   
        #print "tri2", tri2
        #tri2 = U.filter_list(tri2)
        tri2 = sorted(tri2, key=operator.itemgetter(2))
        #print "tri2", tri2
        if len(tri2)<minPS:
            continue
        
        diky = {}
        for y, val in enumerate(tri2):
            diky[val[0]]= y
        #print "diky", diky, i, psi
        posy = diky[psi]
        
        iDy = [x[0] for x in tri2]
        cLeft = [x[1] for x in tri2]
        cRight = [x[2] for x in tri2]
        
        tri3 = []
        j = homoBinaryRight(cRight[posy], cRight[posy:], d, cRight, 'NA')
        if j == 'NA':
            j = posy
        else:
            j = posy + j
        if j>i:
            while j>i:
                tri3.append([iDy[j], cLeft[j], cRight[j]])
                j = j - 1
        j = homoBinaryLeft(cLeft[posy], cLeft[:posy+1], d, cLeft, 'NA')
        if j == 'NA':
            j = posy
        if j <= posy:
            while j <= posy:
                tri3.append([iDy[j], cLeft[j], cRight[j]])
                j = j + 1 
        #tri3 = U.filter_list(tri3)
        tri3 = sorted(tri3, key=operator.itemgetter(1,2))
        
        #only append non-already found groups that are big enough
        if ( tri3 not in sel2 ) and len([x[0] for x in tri3]) >= minPS:
            ss = [x[0] for x in tri3]
            sel2.append(ss)
    
    #sel2 = U.filter_list(sel2)
    
    #print "2", datetime.datetime.now()
   
      
    return sel2        
    
    
    
#    seli =[]
#    sel2 =[]
#    #for each PS, search for the compatible ones in the same orientation
#    for i, psi in enumerate(liste2):
#        #print "+++", i, psi
#        #a ps is obviously compatible with at least itself...
#        tri = [psi[0]]
#        tri2 = [psi]
#        #PS are sorted primarily by left coord and secondarily by right coord
#        #Therefore, for a given PS, compatible PS are the ones right before
#        #(to the left) and right after (to the right). We search for compatible
#        #PS to the left (and right) until we find an uncompatible PS on the
#        #left coord.
#        # THIS SHOULD BE DONE WITH BINARY SEARCH !!!
#        #initialize for left and right search
#        i1 = i+0
#        i2 = i+0
#
#        #search to the left 
#        while 1:
#            if i1<len(liste2)-1:
#                i1=i1+1
#            else:
#                break
#            coord_left_ok = (abs(psi[1]-liste2[i1][1]) <= d)
#            if coord_left_ok:
#                coord_right_ok = (abs(psi[2]-liste2[i1][2]) <= d)
#            else:
#                break
#            if coord_left_ok and coord_right_ok and psi[0] != liste2[i1][0]:
#               #tri.append(liste[i1])
#               tri.append(liste2[i1][0])
#               tri2.append(liste2[i1])
#
#        #search to the right
#        while 1:
#            if i2>0:
#                i2=i2-1
#            else:
#                break
#            coord_left_ok = (abs(psi[1]-liste2[i2][1]) <= d)
#            if coord_left_ok:
#                coord_right_ok = (abs(psi[2]-liste2[i2][2]) <= d)
#            else:
#                break
#            if coord_left_ok and coord_right_ok and psi[0] != liste2[i2][0]:
#               #tri.append(liste[i2])
#               tri.append(liste2[i2][0])
#               tri2.append(liste2[i2])
#
#        #Sort PS IDs nicely
#        tri2.sort(key=lambda x: x[0])
#        tri.sort()
#
#        #only append non-already found groups that are big enough
#        if ( tri not in seli ) and len(tri) >= minPS:
#            seli.append(tri)
#
#        #only append non-already found groups that are big enough
#        if ( tri2 not in sel2 ) and len([x[0] for x in tri2]) >= minPS:
#            sel2.append(tri2)
#
#
#    seli = U.filter_list(seli)



##
##
#    l1 = set(map(tuple, homo))
#    l2 = set(map(tuple, U.filter_list(seli)))
#    if not l1 == l2:
#        print        
#        print homo
#        print "---"
#        print U.filter_list(seli)
#        print
#        
#        
#
#    #print "homosame homo", l1 == l2 
#    #print "homosame new",  seli
#    return seli

#--------------------------------------------------------------------------
def homoOpp(xsome, cand, pairopp, d, minPS, pairOppCoord):
    """Group PS in ori X with respect to their common PS in ori oppX.
    xsome: PS info dictionnary xsome[ps id] ==> PS reformated
    cand: groupe of compatible PS in a given orientation
    pairopp: dictionnary of compatible PS in opposite direction"""
    #print "homoOpp1"
    LPairOpp, LPairOppCoord = {}, {} #dict: key=list of compatible opposite pairs, value=list of PS with the compatble list key
    for un in cand: #for each PS in cand, make dico of the PS compatible in opposite orientation
        tmp = pairopp[un]
        tmp2 = pairOppCoord[un]
        
        if tuple(tmp) not in LPairOpp: #if new group of opposite ori compatible PS, creat new key
            LPairOpp[tuple(tmp)] = [un]
            LPairOppCoord[tuple(tmp2)] = [un]
        else: #else, we just save that this PS is also compatible with an existing key (gtoup of oposite PS)
            LPairOpp[tuple(tmp)].append(un)
            LPairOppCoord[tuple(tmp2)].append(un)
    #print "homoOpp2", len( LPairOpp.keys())
    oppoComp = {} #dict: key=list of compatible opposite pairs, value=list the lists of PS groups made from the key
    for oppG, val in LPairOpp.iteritems():
        tutu = homosame(xsome, list(oppG), d, minPS, pairOppCoord[val[0]]) #make subgroups out of the group oppG
        oppoComp[oppG] = tutu
    
    #print "homoOpp3"
    homo_opp = {}
    for un in cand:
        tutu = oppoComp[tuple(pairopp[un])] #the opposite groups lists are already calculated

        
        if tutu != []:
        #if len(tutu)>=minPS:
            #homo_opp[un]=toto[0]
            homo_opp[un] = tutu
#            if U.ps_fictive(xsome[un]):
#                print "homoOpp ", xsome[un][0], homo_opp[un]
    
    #print "homoOpp4"
#Renvoie tous les PS homogenes de signe oppose
    return homo_opp


#--------------------------------------------------------------------------
def sumCom_bkp(g_homo, new, k, j, inter_kj):
    """Check if k and j are already in a homogenous group"""

    already = False
#    print "SumCom old", g_homo, new
#    print "SumCom new", k, j, inter_kj

    for i, gval in enumerate(g_homo):
#        if 1435 in gval:
#            print "inSumCom new[i]  ",len(new[i]), new[i]
#            print "inSumCom inter_kj", len(inter_kj), inter_kj



#        print "inSumCom ", inter_kj == new[i]
        if (inter_kj == new[i]):
#            print "inSumCom ", i, gval, g_homo, k, j

            if k not in g_homo[i]:
                g_homo[i].append(k)
            if j not in g_homo[i]:
                g_homo[i].append(j)
            already = True

#        elif (inter == new[i]):
#            g_homo[i] += [j for j in l if j not in g_homo[i]]
#            already = True
    return already


#--------------------------------------------------------------------------
def sumCom(g_homo, new, k, j, inter_kj):
    """Check if k and j are already in a homogenous group"""

    already = False
    inter = True

#    print "SumCom old", g_homo, new
#    print "SumCom new", k, j, inter_kj

    for i, gval in enumerate(g_homo):
#        if 1435 in gval:
#            print "inSumCom new[i]  ",len(new[i]), new[i]
#            print "inSumCom inter_kj", len(inter_kj), inter_kj

        newS = set(map(tuple,new[i]))
        inter_kjS = set(map(tuple,inter_kj))
        #print "set:", i, j, inter_kjS, newS

        #if inter_kjS.issubset(newS) and newS != inter_kjS:
        if newS.issubset(inter_kjS) and newS != inter_kjS:
            #Mettre la PS dans l ancien group si c est subset
            #print "ici", inter_kj, new[i]
            if k not in g_homo[i]:
                g_homo[i].append(k)
            if j not in g_homo[i]:
                g_homo[i].append(j)
            #Le reste, ca fait une nouveau g_homo (retourner False)
            inter = False
            return already, list(inter_kjS.difference(newS)), inter


        elif (inter_kj == new[i]):
            #print "inSumCom ", inter_kj == new[i]
            #print "inSumCom ", i, gval, g_homo, k, j

            if k not in g_homo[i]:
                g_homo[i].append(k)
            if j not in g_homo[i]:
                g_homo[i].append(j)
            already = True

#        elif (inter == new[i]):
#            g_homo[i] += [j for j in l if j not in g_homo[i]]
#            already = True
    return already, [], inter
#--------------------------------------------------------------------------
def sumPSinSV(sv):
    """Return the number of PS defining SV."""
    if sv:
        return len(sv[0])+len(sv[1])
    else:
        return 0
#--------------------------------------------------------------------------
def oriPair(l1, l2, xsome):
    """Determine les orientations de 2 PS et les renvoie en tuple"""
    ori1 = xsome[l1][1]+xsome[l1][4]
    if l2:
        ori2 = xsome[l2][1]+xsome[l2][4]
    else:
        ori2 = "NA"
    return (ori1, ori2)

#--------------------------------------------------------------------------
def is_subset(ga, gb):
    """ Check if ga is included in gb. or gb is included in ga"""
    if ga <= gb:
        subset = set(ga).issubset(set(gb))
    else:
        subset = set(gb).issubset(set(ga))
    return subset

#--------------------------------------------------------------------------
def sumSV(SV, g1, g2, oripair):
    """Check if groups of PS are already classified together to form
    potential SV."""

    already = False
    #Check if group has not been classified in a previous scan in another orientation
    opair = oripair[1]
    if opair in SV:
        for si in SV[opair]:
            for s1 in si[1]:
                sens1 = is_subset(si[0], g1) and is_subset(s1, g2)
                sens2 = is_subset(si[0], g2) and is_subset(s1, g1)
                if sens1 or sens2 :
                    already = True
                    break

    return already

#--------------------------------------------------------------------------
def add2Dict_old(dico, cle, svs, deb1, fin1, deb2, fin2, pvalue, xsome, typesv, d):
    """Rajoute une SV a un dictionnaire si elle n'est pas deja presente"""

    newsv = [svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]

    define_SV_borders(newsv, svs[0], svs[1], xsome, typesv, d)
    if cle not in dico:
        dico[cle] = [newsv]
    elif newsv not in dico[cle]:
        dico[cle].append(newsv)


#--------------------------------------------------------------------------
def add2Dict(dico, cle, svs, deb1, fin1, deb2, fin2, pvalue, xsome, typesv, d):
    """Rajoute une SV a un dictionnaire si elle n'est pas deja presente
    ADDED ALEX: L'Etape 4 ne respecte plus l ordre ++ avant -- et
    +- avant -+ car c etait trop lent. C est bcp plus rapide de le faire ici
    """
    
    ###newsv = [svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]
    
    ###if cle == ('--', '++') or cle == ('-+', '+-'):
        #newsv = [svs[1], svs[0], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]
        #if len(svs[0]) == 0:
            #newsv = [svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]
        ###define_SV_borders(newsv, svs[0], svs[1], xsome, typesv, d)
        ###newsv = [svs[1], svs[0], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]
        #else:
            #newsv = [svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]
            #print "AAAA", newsv, svs
         #   define_SV_borders(newsv, svs[0], svs[1], xsome, typesv, d)
         #   newsv = [svs[1], svs[0], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]

    if len(svs[0]) == 0:
        newsv = [svs[1], svs[0], svs[2], svs[3], deb2, fin2, deb1, fin1, pvalue]            
        define_SV_borders(newsv, svs[1], svs[0], xsome, typesv, d)

    else:
        newsv = [svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2, fin2, pvalue]
        define_SV_borders(newsv, svs[0], svs[1], xsome, typesv, d)

    if cle not in dico:
        dico[cle] = [newsv]
    elif newsv not in dico[cle]:
        dico[cle].append(newsv)


#--------------------------------------------------------------------------
def relative_position(ps1, ps2, d, ps_type):
    """ determiner quel est le bras implique pour chaque chromosome,
        afin de savoir quelle coordonnee doit etre plus grande que l'autre
        (sur chaque chromosome) en fait,  selon la diapo 23 du
        1stRUN-clb5workfile.pptx,  il semble que quel que soit le bras
        implique,  c'est toujours coord read - < coord read + il faut juste
        determiner + et -!!!
    """
    ins = False
    transloc = False

    if ps_type == "MP":
        signe = "-"
    else:
        signe = "+"

    co1 = (ps1[1] == signe) and (ps1[3] < ps2[3])
    co2 = (ps1[4] == signe) and (ps1[6] < ps2[6])

    do1 = (ps2[1] == signe) and (ps2[3] < ps1[3])
    do2 = (ps2[4] == signe) and (ps2[6] < ps1[6])


    # si la difference entre les coordonnees sur un meme
    # chromosome est inferieure a 2*d et coordonnees non
    # chevauchante mets dans transloc

    # si la difference entre les coordonnees sur un meme
    # chromosome est
    # inferieure a 2*d et si coordonnees chevauchantes
    # d'un cote ou de l'autre met dans insertion

    d1 = U.parameter_ok(ps1[3], ps2[3], 2*d)
    d2 = U.parameter_ok(ps1[6], ps2[6], 2*d)

    if (d1 and d2):
        if (co1 or do1) and (co2 or do2):
            transloc = True
        elif (co1 or do1) or (co2 or do2):
            ins = True
    elif (d1 or d2) :
        ins = True

    if U.ps_fictive(ps1) and U.ps_fictive(ps2):
        print "relative position #", ps1, "#", ps2, transloc, ins, d1, d2, co1, co2, do1, do2

    return transloc, ins

#--------------------------------------------------------------------------
def candidatOppose(xsome, j, overopps, d, ps_type):

    """ Ne garde dans overopps (liste des PS dans orientation opposee au PS
        xsome[j]) que les PS qui repondent aux criteres sur d et sur orienta
        tion et position ok pour translocations et insertions. Pour les
        translocations il faut repasser dans candidatOpposeTransloc
        """
    cai = xsome[j]
#    if U.ps_fictive(cai):
#        print "candidatOppose 0", cai[0]

    for kind in overopps:

        k = xsome[kind]
        transloc, ins = relative_position(cai, k, d, ps_type)

        if not transloc and not ins:
            overopps.remove(kind)
        if U.ps_fictive(cai) and U.ps_fictive(k):
            print "candidatOppose", j, cai[0], kind, k[0], transloc, ins

#--------------------------------------------------------------------------
def candidatOpposeTransloc(xsome, compat, d, ps_type):
    """ Test supplementaire par rapport candidatOppose
        pour voir si ok transloc"""

#    ori1 = compat[0][3]
#    ori2 = compat[0][2]
    ori1 = compat[3]
    ori2 = compat[2]

    newcompat = []
    ins = []

    #print "candidatOppTransloc ", compat
    #fait des sous-groupes homogenes
    sous_groupes = {}

    for cai in compat[0]:    #pour chaque PS dans ori1
        ps1 = xsome[cai]
        tempo=[]
        for caj in compat[1]:   #test si ok transloc avec toutes les PS opp
            ps2 = xsome[caj]
            transloc, put_ins  = relative_position(ps1, ps2, d, ps_type)
            if transloc: #si PS opp ok, rajoute sous-groupe
                tempo.append(caj)
            elif put_ins :
                ins.append(caj)
            if U.ps_fictive(ps1) and U.ps_fictive(ps2):
                print "candidatOpposeTransloc 1 #",cai," : ", ps1[0], "#",caj, ps2[0], tempo


        # si le sous-groupe est plus petit que le groupe des PS opp et qu'il
        # n'est pas deja definit, on le rajoute.
        if (0 < len(tempo) < len(compat[1])):
            #cle-chaîne a partir de tempo
            tempo_key = ";".join([str(il) for il in tempo])
            sous_groupes.setdefault(tempo_key,[]).append(cai)
        else:
        #Si tout le groupe est homogene, rajoute tel quel
            newcompat=compat

#        if U.ps_fictive(ps1):
#            print "candidatOpposeTransloc 2 #",cai, ps1[0], sous_groupes, newcompat, "compat de depart ", compat

#    if not sous_groupes:
#        newcompat=compat
#    else:
    if sous_groupes:
        for sg_key in sous_groupes:
            #reconvertit la cle en liste de PS homogene
            sg = [int(il) for il in sg_key.split(";")]
            # s'il existe des sous-groupes, renvoyer tous ceux qui contiennent au
            # au moins 3 PS au total
            #if len(sg) + len(sous_groupes[sg_key]) > 2:
            # au moins 2 PS au total
            if len(sg) + len(sous_groupes[sg_key]) > 1:
                newcompat=[sous_groupes[sg_key], sg, ori1, ori2]

    #if 851 in list(set(ins)):
    #    print "inside candidatOpposeTransloc   transloc:", newcompat, "put_ins", list(set(ins))
    
#    print "candidatOpposeTransloc 3 #",cai, ps1[0], sous_groupes, newcompat
    return newcompat, list(set(ins))
#--------------------------------------------------------------------------
def type_SV(deb1, fin1, deb2, fin2, diff_max_geneconv, diff_min_other,
            diff_max_other):
    """ Selon la distance entre les coordonnees sur un meme chromosome,
    distingue les SNP, les GENE conversion"""

    #print "type_SV"
    #parametres de decision
    inf_gconv = 500  # en deca, gene conversion
    inf_snp = 53      # en deca, snp
    # determine la diff max entre coordonnees des deux orientations
    # (en comparant entre eux les min et max)
    diff1 = abs(fin1-deb1)
    diff2 = abs(fin2-deb2)
    type_sv = "Z"

    #si diff coordonnees sur 1 chromosome < inf_snp pb,  classe en snp
    snp = (diff1 < inf_snp) or  (diff2 < inf_snp)

    # si diff coordonnees sur 1 chromosome > inf_snp et sur l'autre < inf_gconv
    # classe en geneconv
    gconv1 = (diff1 < inf_gconv) and  (diff2 > inf_snp)
    gconv2 = (diff2 < inf_gconv) and  (diff1 > inf_snp)

    gconv3 = (diff1 <= diff_max_geneconv)
    gconv3 = gconv3 and (diff_min_other <=  diff2 <= diff_max_other)
    gconv4 = (diff2 <= diff_max_geneconv)
    gconv4 = gconv4 and (diff_min_other <=  diff1 <= diff_max_other)

    if snp:
        type_sv = "snp"
    elif gconv1 or gconv2 or gconv3 or gconv4:
        type_sv = "gconv"

    return type_sv

#--------------------------------------------------------------------------
def isTransloc(xsome, svs, diff_max_geneconv, diff_max_other, diff_min_other,
               geneconv, transloc, snp, subtelcap, d, centrom, subtelo, interne):
    """ Test conditions for Translocation."""
#    for svs in compat:
    cle = oriPair(svs[0][0], svs[1][0], xsome)

    
    n1 = len(svs[0])
    n2 = len(svs[1])

    #print "bbbbbbbbbb", svs
    deb1, fin1, deb2, fin2 = U.Limites(xsome, svs[0], svs[1])
    tsv = "Z"
    if diff_max_geneconv > 1:
        tsv = type_SV(deb1, fin1, deb2, fin2, diff_max_geneconv,
                 diff_max_other, diff_min_other)
    pbal = U.testStatSciPy(n1, n2)
    #print "isTransloc ", svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2, fin2, pbal
    if pbal >= 0.1:
    #if pbal >= 0:    
    #if pbal <= 0.1:
        liste = transloc
        tipe = "tr"
    else:
        liste = subtelcap
        #virer le cote seul pour NRT
        svs[0] = max(svs[0], svs[1], key=len)
        svs[1] = []
        sub1, bks, ins1, bki = isSubtelo_Cap(xsome, centrom, subtelo, svs[0], interne)
        deb1, fin1, deb2, fin2 = U.Extremites(xsome, sub1)
        tipe ="tn"
        
    if tsv == "Z":
        if tipe == "tr":
            add2Dict(liste, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome,
             "tr", d)
        else:
            if sub1:
                add2Dict(liste, cle, [sub1, [], bks, "NA"], deb1, fin1, deb2,
                     fin2, "NA", xsome, "tn", d)
        
        
        
        #add2Dict(liste, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome,
        #         "tr", d)
#        transloc.append((svs[0], svs[1], svs[2], svs[3], deb1, fin1,
#                      deb2, fin2, pvalue))
    elif tsv == "snp":
        add2Dict(snp, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome, "tr",
                 d)
#        snp.append((svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2,
#                 fin2, pvalue))
    elif tsv == "gconv":
        add2Dict(geneconv, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome,
                 "tr", d)
#        geneconv.append((svs[0], svs[1], svs[2], svs[3], deb1, fin1,
#                     deb2, fin2, pvalue))

    if U.ps_fictive(xsome[svs[0][0]]) and U.ps_fictive(xsome[svs[1][0]]):
        print "isTransloc "
        for sv in transloc[cle]:
            print "isTransloc ", cle, sv


#--------------------------------------------------------------------------
def isInsert(xsome, svs, diff_max_geneconv, diff_max_other, diff_min_other,
             geneconv, insert2, snp, subtelcap, d, centrom, subtelo, interne):
    """ Test conditions for Insertion."""

    n1 = len(svs[0])
    n2 = len(svs[1])
    cle = oriPair(svs[0][0], svs[1][0], xsome)

    deb1, fin1, deb2, fin2 = U.Limites(xsome, svs[0], svs[1])
    tsv = "Z"
    if diff_max_geneconv > 1:
        tsv = type_SV(deb1, fin1, deb2, fin2, diff_max_geneconv,
                  diff_max_other, diff_min_other)
    pbal = U.testStatSciPy(n1, n2)
    #print "isInsert2 ", svs[0], svs[1], svs[2], svs[3], deb1, fin1, deb2, fin2,pbal
    if pbal >= 0.1:
    #if pbal >= 0:
    #if pbal <= 0.1:
        liste = insert2
        tipe = "ins"
    else:
        liste = subtelcap
        svs[0] = max(svs[0], svs[1], key=len)
        svs[1] = []
        sub1, bks, ins1, bki = isSubtelo_Cap(xsome, centrom, subtelo, svs[0], interne)
        deb1, fin1, deb2, fin2 = U.Extremites(xsome, sub1)
        tipe ="tn"
    
    if tsv == "Z":
        if tipe == "ins":
            add2Dict(liste, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome,
             "ins", d)
        else:
            if sub1:
                add2Dict(liste, cle, [sub1, [], bks, "NA"], deb1, fin1, deb2,
                     fin2, "NA", xsome, "tn", d)

                     
            #add2Dict(liste, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome,
            #     "ins", d)
    elif tsv == "snp":
        add2Dict(snp, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome, "ins",
                 d)
    elif tsv == "gconv":
        add2Dict(geneconv, cle, svs, deb1, fin1, deb2, fin2, pbal, xsome,
                 "ins", d)
#    if U.ps_fictive(xsome[svs[0][0]]) and U.ps_fictive(xsome[svs[1][0]]):
#        print "isInsert ", xsome[svs[0][0]][0],xsome[svs[1][0]][0], tsv

#--------------------------------------------------------------------------
def isPutativIns(xsome, svs,  diff_max_geneconv, diff_max_other, diff_min_other,
             geneconv, insert1, snp, subtelcap, d, centrom, subtelo, interne):


    cle = oriPair(svs[0][0], svs[1][0], xsome)

    #Visiblement inutile ???
    if cle == ('--','++') or cle == ('-+','+-'):
        svs[0], svs[1] = svs[1], svs[0]


    n1 = len(svs[0])
    n2 = len(svs[1])

    deb1, fin1, deb2, fin2 = U.Limites(xsome, svs[0], svs[1])
    tsv = "Z"
    if diff_max_geneconv > 1:
        tsv = type_SV(deb1, fin1, deb2, fin2, diff_max_geneconv,
                  diff_max_other, diff_min_other)
        #print "isSubtel ", sub1, [], bks, "NA", deb1, fin1, deb2, fin2, "NA"

    pbal = U.testStatSciPy(n1, n2)
    if pbal >= 0.1:
    #if pbal >= 0:
    #if pbal <= 0.1:
        liste = insert1
        tipe = "ins"
    else:
        liste = subtelcap
        svs[0] = max(svs[0], svs[1], key=len)
        svs[1] = []
        sub1, bks, ins1, bki = isSubtelo_Cap(xsome, centrom, subtelo, svs[0], interne)
        deb1, fin1, deb2, fin2 = U.Extremites(xsome, sub1)
        tipe ="tn"

    if tsv == "Z":
        if tipe == "ins":
            add2Dict(liste, cle, svs, deb1, fin1, deb2, fin2, "NA", xsome,
             "ins", d)
        else:
            if sub1:
                add2Dict(liste, cle, [sub1, [], bks, "NA"], deb1, fin1, deb2,
                     fin2, "NA", xsome, "tn", d)
            
            
            #add2Dict(liste, cle, svs, deb1, fin1, deb2, fin2, "NA", xsome,
            # tipe, d)

    elif tsv == "snp":
        add2Dict(snp, cle, svs, deb1, fin1, deb2, fin2, "NA", xsome,
             "ins", d)

    elif tsv == "gconv":
        add2Dict(geneconv, cle, svs, deb1, fin1, deb2, fin2, "NA", xsome,
             "ins", d)

#--------------------------------------------------------------------------
def isSubtel(xsome, centrom, subtelo, interne, sv, diff_max_geneconv,
             diff_max_other, diff_min_other, subtelcap, geneconv, insert1,
             snp, d, minPS):
    """ Test conditions for Subtelomleric capture."""
    ok = False
    sub1, bks, ins1, bki = isSubtelo_Cap(xsome, centrom, subtelo, sv[0], interne)
    #print "isSubtel ", sub1, bks, ins1, bki
    n1 = len(sub1)
    if U.ps_fictive(xsome[sv[0][0]]):
        print "isSubtelo 0", xsome[sv[0][0]][0], sub1, bks, ins1, bki, "n1:", n1
    tsv = "NULL"
    if (n1 >= 2 ):
        ok = True
    if ok:
        #print "CHELOU"
        deb1, fin1, deb2, fin2 = U.Extremites(xsome, sub1)
        #print "SUB1", sub1
        cle = oriPair(sub1[0], [], xsome)
        tsv = "Z"
        if diff_max_geneconv > 1:
            tsv = type_SV(deb1, fin1, deb2, fin2, diff_max_geneconv,
                      diff_max_other, diff_min_other)
            #print "isSubtel ", sub1, [], bks, "NA", deb1, fin1, deb2, fin2, "NA"
        if tsv == "Z":
            if sub1:
                add2Dict(subtelcap, cle, [sub1, [], bks, "NA"], deb1, fin1, deb2,
                     fin2, "NA", xsome, "tn", d)
        elif tsv == "snp":
            if sub1:
                add2Dict(snp, cle, [sub1, [], bks, "NA"], deb1, fin1, deb2, fin2,
                     "NA", xsome, "tn", d)
        elif tsv == "gconv":
            if sub1:
                add2Dict(geneconv, cle, [sub1, [], bks, "NA"], deb1, fin1, deb2,
                     fin2, "NA", xsome, "tn", d)

    if not ok and ins1 != []:
        if U.ps_fictive(xsome[sv[0][0]]):
            print "isSubtelo 0 bis", xsome[sv[0][0]][0], sub1, bks, ins1, bki
#        if n1 >= minPS["ins"] and len(ins1) >= minPS["ins"]:
        if len(ins1) >= minPS["ins"]:
            deb1, fin1, deb2, fin2 = U.Extremites(xsome, ins1)
            tsv = "Z"
            cle = oriPair(ins1[0], [], xsome)
            if diff_max_geneconv > 1:
                tsv = type_SV(deb1, fin1, deb2, fin2, diff_max_geneconv,
                          diff_max_other, diff_min_other)
            #print "isInsert1 ", ins1, [], bks, "NA", deb1, fin1, deb2, fin2, "NA"
            if tsv == "Z":
                add2Dict(insert1, cle, [ins1, [], bki, "NA"], deb1, fin1, deb2,
                         fin2, "NA", xsome, "ins", d)
            elif tsv == "snp":
                add2Dict(snp, cle, [ins1, [], bki, "NA"], deb1, fin1, deb2, fin2,
                         "NA", xsome, "ins", d)
            elif tsv == "gconv":
                add2Dict(geneconv, cle, [ins1, [], bki, "NA"], deb1, fin1, deb2,
                         fin2, "NA", xsome, "ins", d)

#    if U.ps_fictive(xsome[sv[0][0]]):
#        print "isSubtelo ", xsome[sv[0][0]][0], ok, tsv


#--------------------------------------------------------------------------
def dichotoSearch(order, pos1, id1, d, orderAll, ieme):
    """Dicotomy search for compatible PS  """
    #initialisation de la fonction
    try:
        lenorder = len(order) #longueur de la liste en cours
        iddic = lenorder/2 #index de la valeur "centrale" de la liste
        cai1 = order[iddic][1] #valeur de la valeure "centrale"
        if ieme=="NA": #cas ou c est une nouvelle liste (pas encore recurence)
            ieme = iddic
    except: #si probleme, retourne une liste vide
        return []

    #cas ou on peut s arreter (liste de longueur 1 avec une coord pas compatible)
    if lenorder == 1 and (cai1-2*d > pos1 or cai1+2*d < pos1):
        return []

    #autres cas particuliers (liste de taille 1 et 2)
    elif lenorder <= 2:
        opp = scanOrder(orderAll, pos1, ieme, d, id1)
        if lenorder == 2:
            opp2 = scanOrder(orderAll, pos1, ieme-1, d, id1)
            return opp+opp2
        else:
            return opp

    #Si la valeur cherchee est plus petite que la valeur a la position actuelle
    #dans la liste, il faut prendre la sous partie gauche de la liste
    elif cai1-2*d > pos1:
        #calculer l indice de correction pour le calcule de la position de la
        #valeur mediane dans la future liste par rapport a la toute premiere
        #liste (orderAll)
        #Si on va a gauche, la correction ne depend que de la parite de la
        #longueur de la nouvelle liste
        if len(order[:iddic])%2 == 0:
            correc = 0
        else:
            correc = 1

        return dichotoSearch(order[:iddic], pos1, id1, d, orderAll, ieme-iddic/2-correc)

    #Si la valeur cherchee est plus grande que la valeur a la position actuelle
    #dans la liste il faut prendre la sous partie doite de la liste
    elif cai1+2*d < pos1:
        #calculer l indice de correction pour le calcule de la position de la
        #valeur mediane dans la future liste par rapport a la toute premiere
        #liste (orderAll)
        #Si on va a droite, la correction depend de la parite de la
        #longueur de la nouvelle liste mais aussi de la parite de la longueur
        # de la liste actuelle
        if lenorder%2 ==1 and len(order[iddic:])%2 == 0:
            correc = 1
        else:
            correc = 0

        return dichotoSearch(order[iddic:], pos1, id1, d, orderAll, ieme+iddic/2+correc)

    #Si la valeur actuelle est dans un intervale de coord compatible (<2d)
    #commencer la recherche des PS compatibles (a partir de la "ieme" position)
    #dans la liste initiale
    #cas ou la valeur cherchee est a gauche de la position actuelle
    elif cai1-2*d <= pos1 and cai1 >= pos1:
        opp = scanOrder(orderAll, pos1, ieme, d, id1)
        return opp

    #cas ou la valeur cherchee est a droite de la position actuelle
    elif cai1 <= pos1 and cai1+2*d >= pos1:
        opp = scanOrder(orderAll, pos1, ieme, d, id1)
        return opp

#--------------------------------------------------------------------------
def scanOrder(order, pos1, idps, d, id1):
    """Scans "order" list for PS compatible with pos1
    pos1 must be compatible with the initial cai1 value for this function to
    work properly (should be the case for dichotoSearch function"""
    opp=[]

    #coord (cai) and id (cID) of the starting PS for the search to the left
    cai1 = order[idps][1]
    cID = order[idps][0]

    #coord (cai) and id (cID) of the starting PS for the search to the right
    cai1_2 = order[idps][1]+0
    cID_2 = order[idps][0]+0
    idps_2 = idps+0

    #Search to the left of cID (including cID)
    OK = True
    while OK == True:
        if U.parameter_ok(pos1, cai1, 2*d):
            opp.append(cID)

            if idps > 0:
                idps=idps-1
            else:
                OK = False
            cai1 = order[idps][1]
            cID = order[idps][0]
        else:
            OK = False

    #Search to the right of cID (including cID)
    OK = True
    while OK == True:
        if U.parameter_ok(pos1, cai1_2, 2*d):
            opp.append(cID_2)
            if idps_2 < len(order)-1:
                idps_2=idps_2+1
            else:
                OK = False
            cai1_2 = order[idps_2][1]
            cID_2 = order[idps_2][0]
        else:
            OK = False

    return opp

#--------------------------------------------------------------------------
def cleanSV(SV, min_minPS):
    """Before Classif, removes elements from SV dict that are too small
	and that dont need to go to the Classif function"""
    SV = dict(SV)
    for ori, SVs in SV.iteritems():
        for sv in SVs:
            for sv2 in sv[1]:
                if len(sv2)<min_minPS:
                    sv[1].remove(sv2)
#--------------------------------------------------------------------------

def updateSubtelcap(subtelcap, inter, ps_min_tn, xsome):
    """ Verifies that no PS from subtelcap are found in inter """
    #1: remove SV if 1PS is already in INTER 
    subtelcapNew = {}
    for key, allSV in subtelcap.iteritems(): #for each ori
        if key[1]!='NA':
            key = tuple([key[0], 'NA'])
        
        l =[]
        for sv in allSV: #for each sv
        #    print "--- working sv", key, sv
#            if len(sv[0])==0:
#                print "11aa", sv
#                nsv = U.tuples2list(sv)
#                nsv[0], nsv[1] = nsv[1], nsv[0]                
#                sv = nsv

            present = False
            sv0 = []

            for PS in sv[0]: #for each PS in sv
                #print "-------working PS", PS, PS not in inter, xsome[PS][0], key
                if PS not in inter: #if PS not in the inter (INS and TR)
                    sv0.append(PS) #append new liste
                else:
                    present = True
                    
            if not present: #if everything ok, just keep as it was
                l.append(sv)
            else:
                if len(sv0)>= ps_min_tn: #else, only keep undetected PS if len still ok
                    nsv = U.tuples2list(sv)
                    nsv[0] = sv0
                    l.append(nsv)

        
        if key not in subtelcapNew:
            subtelcapNew[key] = l
        else:
            subtelcapNew[key].extend(l)



    #2: remove duplicates due to "dual" orientation search (should not happen)
    subtelcapNew2 = {}
    for key, val in subtelcapNew.iteritems():
        #print val
        newVal = list(set(map(U.list2tuples, val)))
        subtelcapNew2[key]=newVal
    
    return subtelcapNew


#--------------------------------------------------------------------------

   
def rmDupRmSub(a,b):
    """ Remove duplicqtes qnd sub-duplicqtes"""
    if a and b:
        #remove duplicates and merge into 1 list
        a_set = set(map(U.list2tuples, a))
        b_set = set(map(U.list2tuples, b))
        a_b = list(a_set.union(b_set))
        #print a_set, "----", b_set, "----", a_b
        #remove sub-svs (svs included in others)
        svs1 = [set(x[0]) for x in a_b] #left part of all svs
        svs2 = [set(x[1]) for x in a_b] #right part of all svs
        #
        a_b2 = []
        for j, sv in enumerate(a_b): #for each sv
            ok = True 
            
            for i,sv1 in enumerate(svs1): 
                if i!=j:
                    if set(sv[0]).issubset(sv1) and set(sv[1]).issubset(svs2[i]) : #check if sv is a subset of another one
                        ok = False
                        break
                    elif set(sv[1]).issubset(sv1) and set(sv[0]).issubset(svs2[i]) : #check if sv is a subset of another one
                        if set(sv[1]) == sv1 and set(sv[0]) == svs2[i]  and j > i: #if same but opposite orientation sv
                            ok = False
                            break
                        elif set(sv[1]) != sv1 or set(sv[0]) != svs2[i]:
                            ok = False
                            break
            if ok  : #if not, append to list
                a_b2.append(sv)
    else:
        if a:
            a_b2 = a
        elif b:
            a_b2 = b
        else:
            a_b2 = []   
    return a_b2


#--------------------------------------------------------------------------    
def rmSubNRT(dicos, d, xsome, type2sv):
    """remouve NRT (finaly i wrote a generic function, should work with INS and
    TR, not tested though) contained in each other """
    #1: identify compatible NRT, creat all 2 by 2 comparisons  
    toFuse = {}    
    for key in dicos.keys():
        #print "KIKI", key
        toFuse[key] = []
        for ieme, i in enumerate(dicos[key]):
            ori = i[2]
            #print "iiiiii", i
            d1,f1,d2,f2 = i[4],i[5],i[6],i[7]
            m1, m2 = (d1+f1)/2, (d2+f2)/2
            for jeme, j in enumerate(dicos[key]):
                if jeme>ieme:
                    orj = j[2]
                    dj1,fj1,dj2,fj2 = j[4],j[5],j[6],j[7]
                    mj1, mj2 = (dj1+fj1)/2, (dj2+fj2)/2
                    #print "VALUES", d1,f1,d2,f2, ori, dj1,fj1,dj2,fj2, orj, d, d1-fj1, f1-dj1
                    #if abs(d1-fj1)<d and abs(f1-dj1)<d and ori == orj:
                    if abs(m1-mj1)<d and abs(m2-mj2)<d and ori == orj:
                        #print "FOUND1"
                        toFuse[key].append([ieme, jeme])
    #print "STARTING FUSION"
    #2: fuse compatible subtel as maximal cliques
    nDic = {}    
    for key,groups in toFuse.iteritems():
        
        G=nx.Graph()
        #print "KIKI1", key, len(groups)
        if len(groups)>2000000:
            idxg = random.sample(range(len(groups)), 2000000)
            #print "idxg", idxg
            G.add_edges_from([groups[x] for x in idxg])
            print "Warning: High coverage, dumped some RP", key, len([groups[x] for x in idxg])
        else:
            G.add_edges_from(groups)
            #print "KIKI2", key, len(groups)
        
        
        
        
        #print "size of graph", sys.getsizeof(G)
        #print "KIKI2.0.1", key
        fusedG = nx.find_cliques_recursive(G)
        #print "KIKI2.0.2", key, sys.getsizeof(fusedG)
        f = set((val for sublist in fusedG for val in sublist)) #flatten list by position of the NRT in the list
        #print "KIKI2.0.3", key
        singletonsNRT = list(set(range(0,len(dicos[key]))).difference(f)) #search singletons (NRT compatible with no other one)
        #print "KIKI2.0.4", key
        fusedG = fusedG + [[x] for x in singletonsNRT] #add singletons
        #print "KIKI2.2", key
        fusedListe = []
        for elt in fusedG:
            #print "KIKI2.3", key
            tmp = []
            mm = []
            pp = [] #will be empty for NRT
            for SV in elt:
                mm.extend(dicos[key][SV][0])
                pp.extend(dicos[key][SV][1])
                mm = list(set(mm))
                pp = list(set(pp))
            tmp.extend([mm, pp])
            tmp.extend(dicos[key][SV][2:4])
            deb1, fin1, deb2, fin2 = U.Extremites(xsome, mm)
            tmp.extend([deb1, fin1, deb2, fin2])
            tmp.append(dicos[key][SV][8])
            define_SV_borders(tmp, tmp[0], tmp[1], xsome, type2sv, d)
            fusedListe.append(tmp)
        nDic[key] = fusedListe
    return nDic
        

    
#--------------------------------------------------------------------------
def rmDupSV(dicos):
    """removes duplicated SVs from dict after Etap4_quick"""
    
    #clean signs
    opposite = {'++':'--', '--':'++', '+-':'-+', '-+':'+-'}
    flag = False
    for key in dicos.keys():
        if key[1] == 'NA':
            key2 = (key[0], opposite[key[0]])
            flag = True
            if key2 in dicos.keys():
                dicos[key2] = dicos[key] + dicos[key2]
            else:        
                dicos[key2] = dicos[key]
        if flag:
            pope = dicos.pop(key, None)
            flag = False
    
    #remove dups
    ppmm = dicos.get(('++','--'))
    mmpp = dicos.get(('--','++'))
    pmmp = dicos.get(('+-','-+'))
    mppm = dicos.get(('-+','+-'))

    #compare ('++','--') with ('--','++')
    ppmm_mmpp2 = rmDupRmSub(ppmm, mmpp)
    pmmp_mppm2 = rmDupRmSub(pmmp, mppm)
    
#    if ppmm and mmpp:
#        #remove duplicates and put in single orientation
#        ppmm_set = set(map(U.list2tuples, ppmm))
#        mmpp_set = set(map(U.list2tuples, mmpp))
#        ppmm_mmpp = list(ppmm_set.union(mmpp_set))
#        
#        #remove sub-svs (svs included in others)
#        svs1 = [set(x[0]) for x in ppmm_mmpp] #left part of all svs
#        svs2 = [set(x[1]) for x in ppmm_mmpp] #right part of all svs
#        #print "sv1[0]", svs1[0]
#        ppmm_mmpp2 = []
#        for j, sv in enumerate(ppmm_mmpp): #for all sv
#            ok = True 
#            
#            for i,sv1 in enumerate(svs1): #check if sv is a subset of another one
#                if i!=j:
#                    if set(sv[0]).issubset(sv1) and set(sv[1]).issubset(svs2[i]) :
#                        ok = False
#                        break
#                    elif set(sv[1]).issubset(sv1) and set(sv[0]).issubset(svs2[i]) :
#                        if set(sv[1]) == sv1 and set(sv[0]) == svs2[i]  and j > i:
#                            ok = False
#                            break
#                        elif set(sv[1]) != sv1 or set(sv[0]) != svs2[i]:
#                            ok = False
#                            break
#            if ok  : #if not, append to list
#                ppmm_mmpp2.append(sv)
#                
#    else:
#        if ppmm:
#            ppmm_mmpp2 = ppmm
#        elif mmpp:
#            ppmm_mmpp2 = mmpp
#        else:
#            ppmm_mmpp2 = []


    #compare ('+-','-+') with ('-+','+-')
#    if pmmp and mppm:
#        #remove duplicates and put in single orientation
#        pmmp_set = set(map(U.list2tuples, pmmp))
#        mppm_set = set(map(U.list2tuples, mppm))
#        pmmp_mppm = list(pmmp_set.union(mppm_set))
#        
#        #remove sub-svs (svs included in others)
#        svs1 = [set(x[0]) for x in pmmp_mppm]
#        svs2 = [set(x[1]) for x in pmmp_mppm]
#        
#        pmmp_mppm2 = []
#        for j, sv in enumerate(pmmp_mppm): #for all sv
#            ok = True 
#            
#            for i,sv1 in enumerate(svs1): #check if sv is a subset of another one
#                if i!=j:
#                    if set(sv[0]).issubset(sv1) and set(sv[1]).issubset(svs2[i]) :
#                        ok = False
#                        break
#                    elif set(sv[1]).issubset(sv1) and set(sv[0]).issubset(svs2[i]) :
#                        if set(sv[1]) == sv1 and set(sv[0]) == svs2[i]  and j > i:
#                            ok = False
#                            break
#                        elif set(sv[1]) != sv1 or set(sv[0]) != svs2[i]:
#                            ok = False
#                            break
#            if ok  : #if not, append to list
#                pmmp_mppm2.append(sv)     
#        
#    else:
#        if pmmp:
#            pmmp_mppm2 = pmmp
#        elif mppm:
#            pmmp_mppm2 = mppm
#        else:
#            pmmp_mppm2 = []

    return {('++','--'):ppmm_mmpp2, ('+-','-+'):pmmp_mppm2 }


#--------------------------------------------------------------------------
def Etape1(xsome, candidats):
    """Classify PS with respect to orientation of the reads."""
    for i, cai in enumerate(xsome):
        #print cai
        oric = cai[1] + cai[4]
        candidats[oric].append(i)

#--------------------------------------------------------------------------
def Etape2(xsome, candidats, oris, d, oppose, ps_type, min_MinPS):
    """Identifie les paires compatibles entre orientations opposees
        Si pas de PS opposees,  met liste vide"""
    pairopp = {}
    pairdiff = {}

    for cand in candidats[oris]:
        #print xsome[cand][3]
        overall, overopps = U.overlap(xsome, candidats, cand, oris,
                                                 d, oppose)
        if overopps:
            # Candidats ok pour les geneconversions et les inversions.
            # Pour les translocations, il faut refiltrer avec candidatOpposeTransloc
            candidatOppose(xsome, cand, overopps, d, ps_type)
        pairopp[cand] = overopps
        pairdiff[cand] = overall
        if U.ps_fictive(xsome[cand]):
            print "Etape2IN", cand, "*", xsome[cand][0], overopps, overall
            for j in overopps:
                print "Etape2IN overopps",j, xsome[j][0]
            for j in overall:
                print "Etape2IN overall",j, xsome[j]

    return pairopp, pairdiff

#--------------------------------------------------------------------------
def Etape3(xsome, liste, d, min_MinPS):
    """Make homogenous sub-groups of PS in given orientation."""
    i = 0
    L = len(liste)
    stop = 0
    sel = []
    while i < L and stop == 0:
        tri = [liste[i]]
        cai = xsome[liste[i]]
        if U.ps_fictive(cai) :
            print "Etape3 #1", liste[i], cai[0]

        j = 0
        while j < L:
            if i != j:
                caj = xsome[liste[j]]
                coord_left_ok = (abs(cai[3]-caj[3]) < d)
                coord_right_ok = (abs(cai[6]-caj[6]) < d)
                if coord_left_ok and coord_right_ok:
                    if U.ps_fictive(cai) and U.ps_fictive(caj):
                        print "Etape3 #2", liste[i], cai[0], liste[j], caj[0], tri
                    tri.append(liste[j])
                if U.ps_fictive(cai) and U.ps_fictive(caj) and i != j:
                    print "Etape3 #3", liste[i], cai[0], liste[j], caj[0], coord_left_ok, coord_right_ok, tri

            j = j+1


        #ne garde que les groupes qui ne sont pas deja formes
        U.inde_groups(sel, tri)
        if len(tri) == L:
            stop = 1
        else:
            i = i+1

#    for j in liste:
#        print "Etape3 #4", j, xsome[j][0]

    #Alex: only keep homo_same sub_groups that are big enough
    sel = [x for x in sel if len(x)>(min_MinPS-1)]

    return sel

#--------------------------------------------------------------------------
def Etape1_quick(xsome, candidats):
    """Classify PS with respect to orientation of the reads.
    Keeps the left and right coord of each PS
    Returns 2 additionnal lists:
        order1 is ordered by left coord
        order2 is ordered by right coord"""
    order1 = {"--":[], "++":[], "+-":[], "-+":[]}
    order2 = {"--":[], "++":[], "+-":[], "-+":[]}
    for i, cai in enumerate(xsome):
        #print cai
        oric = cai[1] + cai[4]
        #candidats[oric].append([i,cai[3], cai[6]])
        candidats[oric].append(".".join([str(i),str(cai[3]), str(cai[6])]))

        order1[oric].append([i,cai[3]])
        order2[oric].append([i,cai[6]])
    order1 = U.sortListsInDicos(order1, 1)
    order2 = U.sortListsInDicos(order2, 1)
    return order1, order2


#--------------------------------------------------------------------------
def Etape2_quick(order1, order2, candidats, oris, d, oppose, ps_type, min_MinPS):
    """Identifie les paires compatibles entre orientations opposees
        Si pas de PS opposees,  retourne liste vide"""
    pairopp = {}
    pairopp2 = {}
    signeop = oppose[oris]

    #pour chaque PS, trouver les PS compatibles en cherchant par dichotomie
    #sa coord left dans order1 et sa coord right dans order2.
    #candidats: [id, left, right]
    for cand in candidats[oris]:
        #search for PS compatible with left coord
        pairoppList1 = dichotoSearch(order1[signeop], cand[1], cand[0], d, order1[signeop],"NA")
        #search for PS compatible with right coord
        pairoppList2 = dichotoSearch(order2[signeop], cand[2], cand[0], d, order2[signeop],"NA")

        #mettre dans le dicos en virant les doublons
        pairopp[cand[0]] = list(set(set(pairoppList1).intersection(pairoppList2)))
        pairopp2[".".join([str(x) for x in cand])] = list(set(set(pairoppList1).intersection(pairoppList2)))



    
    order1Dic = {}
    order2Dic = {}
    for val in order1[signeop]:
        order1Dic[val[0]]= val[1]
    for val in order2[signeop]:
        order2Dic[val[0]]= val[1]
        
    pairOppCoord = {}
    for PS, g in pairopp.iteritems():
        ng = []
        for oppPS in g:
            Left = order1Dic[oppPS]
            Right = order2Dic[oppPS]
            ng.append(".".join( [str(oppPS), str(Left), str(Right)] ) )
        pairOppCoord[PS] = ng
        

    return pairopp, pairopp2, pairOppCoord



#--------------------------------------------------------------------------
def Etape3_quick(xsome, liste, d, min_MinPS, oris):
    """Make homogenous sub-groups of PS in given orientation."""

    #get the PS IDs
    liste = map(U.conv, liste)
    #liste_light = [x[0] for x in liste]

    #sort the list by 1) coord left and 2) right coord
    liste = sorted(liste, key=operator.itemgetter(1, 2))

    sel =[]
    sel2 =[]
    #for each PS, search for the compatible ones in the same orientation
    for i, psi in enumerate(liste):
        #print "+++", i, psi
        #a ps is obviously compatible with at least itself...
        tri = [psi[0]]
        tri2 = [psi]
        #PS are sorted primarily by left coord and secondarily by right coord
        #Therefore, for a given PS, compatible PS are the ones right before
        #(to the left) and right after (to the right). We search for compatible
        #PS to the left (and right) until we find an uncompatible PS on the
        #left coord.
        #initialize for left and right search
        i1 = i+0
        i2 = i+0

        #left search
        while 1:
            if i1<len(liste)-1:
                i1=i1+1
            else:
                break
            coord_left_ok = (abs(psi[1]-liste[i1][1]) < d)
            if coord_left_ok:
                coord_right_ok = (abs(psi[2]-liste[i1][2]) < d)
            else:
                break
            if coord_left_ok and coord_right_ok and psi[0] != liste[i1][0]:
               #tri.append(liste[i1])
               tri.append(liste[i1][0])
               tri2.append(liste[i1])

        #right search
        while 1:
            if i2>0:
                i2=i2-1
            else:
                break
            coord_left_ok = (abs(psi[1]-liste[i2][1]) < d)
            if coord_left_ok:
                coord_right_ok = (abs(psi[2]-liste[i2][2]) < d)
            else:
                break
            if coord_left_ok and coord_right_ok and psi[0] != liste[i2][0]:
               #tri.append(liste[i2])
               tri.append(liste[i2][0])
               tri2.append(liste[i2])

        #Sort PS IDs nicely
        tri2.sort(key=lambda x: x[0])
        tri.sort()

        #only append non-already found groups that are big enough
        if ( tri not in sel ) and len(tri) >= min_MinPS:
            sel.append(tri)

        #only append non-already found groups that are big enough
        if ( tri2 not in sel2 ) and len([x[0] for x in tri2]) >= min_MinPS:
            sel2.append(tri2)


    #Remove sub groups
    sel = U.filter_list(sel)
    sel2 = U.filter_listOfLists(sel2)
    


    #Creat a subset of candidate PS for which to search opposit candidates
    #only PS that can form a group of at leat 2 PS in at leat 1 orientation are valid
    candidats2 = {oris:[item for sublist in sel2 for item in sublist]}

    return sel, candidats2

#--------------------------------------------------------------------------
def Etape4_quick(xsome, SV, homo_same, pairopp, pairdiff, ori, d, minPS, min_MinPS, pairOppCoord):
    """Rempli la liste des SV en rangeant toujours la liste dans l'ordre:
        "groupe ++ puis groupe --" ou "groupe +- puis groupe -+
        1ere liste: PS dans un sens - 2eme liste: PS dans l'autre sens
    """
    #print "ggggg", pairOppCoord
    ieme = 1
    V = 500
    t = datetime.datetime.now()
    SV[ori] = []
    PSdiff = sorted(pairdiff.keys())
    
    #print homo_same
    for g in homo_same:
        #print homo_same
        #print "ggg", g
        if ieme%V == 0:
            tnow = datetime.datetime.now() - t
            print "\t\tSpeed is", round(V/tnow.total_seconds(),2), "cands/s"
            t = datetime.datetime.now()

        #Verifie s'il y a des PS dans autres orientations
        none = True
        Fictive = False
        for i in g:
            if U.binary_search(PSdiff, i) != -1: #cherche des PS opposee pour i (s il n y en a pas, ca sert a rien de continuer...)
                none = none and (not pairdiff[i])
                if not none: #Il suffit d une seule PS opp pour s arreter
                    break
            if U.ps_fictive(xsome[i]):
                print "Etap4 first ", xsome[i][0],i, pairdiff, none
                Fictive = True


        #listes des nouveaux groupes de PS a rajouter ensemble dans SV
        g_homo = [] #celles comprises dans g
        new = []    #celles dans le sens oppose de g
        g_homo2 = [] #celles comprises dans g
        new2 = []    #celles dans le sens oppose de g

        #Si pas de PS dans autre orientation, ajoute le groupe de PS homogene
        # g a la liste des SV avec
        if none:
            if g:
                g.sort()
            if Fictive:
                print "Etap4 second 1 seul groupe ", xsome[g[0]][0],g
            SV[ori].append([g, [[]]])
        else:
             # S'il y a des PS dans orientation opposee, pour chaque PS de g,
             # dresse les groupes homogenes de PS orientation opposee
             # dictionnaire g_homo_opp: cle = PS de g, valeur = liste des PS homogenes opposees
            #print "debut", datetime.datetime.now()
            g_homo_opp = homoOpp(xsome, g, pairopp, d, min_MinPS, pairOppCoord) #THIS IS THE BOTTLENECK !!!
            #print "g_homo_opp", g_homo_opp
            #print "PAIROPP", pairopp
            #print "fin", datetime.datetime.now()

            if Fictive:
                print "Etap4 second 2 groupes ", xsome[g[0]][0],g
                print "Etap4 second 2 groupes g_homo_opp", g_homo_opp


            #ATTENTION, cas particulier, si une seul PS de g possede 
            #des PS opposées ==> c'est une NRT car l'algo ne fonctionne 
            #que pour des groupes >1.
            if len(g_homo_opp.keys())<len(g):
                #print "WARNING, this should not happen in homoOpp"
                if Fictive:
                    print "Etap4 second 1 seul groupe ", xsome[g[0]][0],g
                SV[ori].append([g, [[]]])                

            
            else:
                #listes des nouveaux groupes de PS a rajouter ensemble dans SV
    #            g_homo = [] #celles comprises dans g
    #            new = []    #celles dans le sens oppose de g
                #CEST CA !!!new2 = reduce(lambda x, y: U.intersectall(x, y, minPS), newD.values())  #this is like below but with a reduce to limite the number of combinations
                #listDesGroupes = reduce(lambda x,y: tuple(set([tuple(i) for i in x]+[tuple(j) for j in y])), newD.values())
                #print "g_homo_opp LA", len(g_homo_opp.keys()), [len(val) for key,val in g_homo_opp.iteritems()]
                complexiteD = {}
                for keyy,vall in g_homo_opp.iteritems():
                    cle = tuple(tuple(i) for i in vall)                   
                    if cle in complexiteD:
                        complexiteD[cle].append(keyy)
                    else:
                        complexiteD[cle]=[keyy]
                complexiteD2 = {}
                for k,v in complexiteD.iteritems():
                    for oppG in k:
                        if oppG not in complexiteD2:
                            complexiteD2[oppG]=v
                        else:
                            complexiteD2[oppG].extend(v)
                #print "complexiteD2", complexiteD2
                g_homo_opp_simple = {}
                for k,v in complexiteD.iteritems():
                    #print "V", v
                    g_homo_opp_simple[v[0]] = [list(i) for i in k]
                #print "complexiteD", g_homo_opp_simple
                ###-----------------------------------------------------------------------------------------------
                
                #liste des PS de g avec des groupes ds g_homo_opp
                lisg = g_homo_opp_simple.keys()
                lisnul = [l for l in g if l not in lisg] #ps de g sans ps opposees forme subtelo_cap/transloc non recip
                if lisnul :
                    SV[ori].append([lisnul, [[]]])
                lisg.sort()
                Lh = len(lisg)

                # s'il existe au moins 2 PS de g avec des groupes g_homo_opp
                if Lh > 1:
                    # Regarde si les PS copines signe opposes sont partagees par
                    # d'autres PS de meme signe
                    for j in range(Lh-1):
                        for k in range(j+1, Lh): #inter_kj = PS OPP compatibles avec k et j
                            inter_kj = [l for l in g_homo_opp_simple[lisg[k]] if\
                                     l in g_homo_opp_simple[lisg[j]] ]

                            if inter_kj:

                            # S'il existe au moins 2 PS de g avec des PS opposes
                            # regarde si elles ne font pas partie d'un groupe deja definit
                            # dans new
    
                            #C'est là que ça merde : Il faut rajouter toutes les PS de g qui ont des PS opp communes
                            #or ici ça ne marche pas
                                inter_kj.sort()
    
                                already, inter_kj2, inter = sumCom(g_homo, new, lisg[k], lisg[j], inter_kj) #ca ca append g_homo
    
                                if not already:
                                    g_homo.append([lisg[j], lisg[k]])
                                    if not inter:
                                        new.append(inter_kj2)
                                    else:
                                        new.append(inter_kj)
                #print "DEB", g_homo, new
                
                    for k,v in complexiteD.iteritems(): #complexiteD
                        #if len(v) >= min_MinPS:
                        g_homo.append(sorted(v))
                        new.append([list(i) for i in k])
                    #print "COMPA", g_homo, new

                ###-----------------------------------------------------------------------------------------------                
                
                
                
                
                
                
                
                
                
#                ###ORI-----------------------------------------------------------------------------------------------
#                #liste des PS de g avec des groupes ds g_homo_opp
#                lisg = g_homo_opp.keys()
#                lisnul = [l for l in g if l not in lisg] #ps de g sans ps opposees forme subtelo_cap/transloc non recip
#                if lisnul :
#                    SV[ori].append([lisnul, [[]]])
#                lisg.sort()
#                Lh = len(lisg)
#                if Fictive:
#                    print "Etap4 trois.a 2 groupes ", xsome[g[0]][0],g, lisg, Lh
#                # s'il existe au moins 2 PS de g avec des groupes g_homo_opp
#                if Lh > 1:
#                    # Regarde si les PS copines signe opposes sont partagees par
#                    # d'autres PS de meme signe
#                    for j in range(Lh-1):
#                        for k in range(j+1, Lh): #inter_kj = PS OPP compatibles avec k et j
#                            inter_kj = [l for l in g_homo_opp[lisg[k]] if\
#                                     l in g_homo_opp[lisg[j]] ]
#                            #if 434 in lisg:
#                            #print "FLAG2:", g, "inter_kj", inter_kj
#                            if Fictive:
#                                print "Etap4 trois.b 2 groupes ", g[0], \
#                                xsome[g[0]][0], lisg[k], lisg[j], "* inter ",\
#                                inter_kj
#                                print "Etap4 trois.b.a 2 groupes ", xsome[g[0]][0],\
#                                g_homo_opp[lisg[k]], g_homo_opp[lisg[j]]
#                            if inter_kj:
#    
#                            # S'il existe au moins 2 PS de g avec des PS opposes
#                            # regarde si elles ne font pas partie d'un groupe deja definit
#                            # dans new
#    
#                            #C'est là que ça merde : Il faut rajouter toutes les PS de g qui ont des PS opp communes
#                            #or ici ça ne marche pas
#                                inter_kj.sort()
#                                already, inter_kj2, inter = sumCom(g_homo, new, lisg[k], lisg[j], inter_kj)
#
#                                if not already:
#                                    #print "G_HOMO", g_homo, [lisg[j], lisg[k]]
#                                    g_homo.append([lisg[j], lisg[k]])
#                                    if not inter:
#                                        new.append(inter_kj2)
#                                    else:
#                                        new.append(inter_kj)
#                                #print "sumCom ", lisg[j], lisg[k],  g_homo, "new", new, inter_kj
#                                if Fictive:
#                                    print "After sumCom ",xsome[g[0]][0], g_homo, "*",inter_kj,"*already",already,"*new ",new
#                    if Fictive:
#                        print "Etap4 trois.e 2 groupes ", xsome[g[0]][0],g_homo,new
##                    print "COMPA", g_homo, g_homo2
#                    print "COMPA", g_homo, new
#                    print "SECOND", g_homo2, new2
#                    #sys.exit()
                ###-----------------------------------------------------------------------------------------------
                elif Lh == 1:
                    #s'il existe 1 seul PS de g avec des groupes g_homo_opp
                    g_homo = [lisg]
                    new = [g_homo_opp.values()[0]]
                    if Fictive:
                        print "Etap4 trois.f 1 PS g / 1 PS g_opp", lisg, new[0]
    #                if len(pp) > 0:
    #                    new.append([o, opp])
                else:
                    #si aucun PS de g n'a de groupes g_homo_opp
                    if len(g) > 1:
                        g.sort()
                        SV[ori].append([g, [[]]])
    
    
                nb_min_ps = U.getMinPSvalue(minPS)
                #new = [x for x in new if len(x)>= nb_min_ps]
                for c in range(len(new)):
                    # pour tous les groupes formes de g_homo et new
                    # (groupes de PS associes formant SV)
                    #if ori == "++" or ori == "+-":
                    if len(new[c]) + len(g_homo) >= nb_min_ps:
                        new[c].sort()
                        g_homo[c].sort()

                        if Fictive:
                            #print "Etap4 2 groupes compatibles ", xsome[new[c][0]][0],g_homo,new
                            print "ETap4 Boucle new ",c, g_homo[c], new[c]
    
                        orip = oriPair(g_homo[c][0],new[c][0][0], xsome)
                        if not sumSV( SV,  g_homo[c], new[c][0], orip):
                            #if ori== '++' or ori=='+-':
                            SV[ori].append([g_homo[c], new[c]])
                            #else:
                            #    SV[ori].append([new[c], g_homo[c]])
    
                        if Fictive:
                            print "Etap4 last Fictive ", SV[ori][-1]
        ieme = ieme + 1
        
    #Remove groups that don't pass min_ps
    #cleanSV(SV, min_MinPS)
    cleanSV(SV, 1)
    


#--------------------------------------------------------------------------

def Etape6_quick(transloc, insert1, insert2, subtelcap, ps_min_tn, xsome, d):
    """ Etape4_quick leads to duplicates because we search both orientations (-- and ++ )
	and (-+, +-) and do not classify SV given the orientation. """
     
    transloc = rmDupSV(transloc) #rmDupSV only remove SV that are present in duplicates or subSV of already present SV
    #print "after tr"    
    insert1 = rmDupSV(insert1)
    #print "after ins1"
    insert2 = rmDupSV(insert2)
    #print "after ins2"
    
    transloc = rmSubNRT(transloc, d, xsome, 'tn')
    #print "after trN"    
    insert1 = rmSubNRT(insert1, d, xsome, 'ins') 
    #print "after ins1N"
    insert2 = rmSubNRT(insert2, d, xsome, 'ins')  
    #print "after ins2N"
    subtelcap = rmSubNRT(subtelcap, d, xsome, 'tn') #rmSubNRT fuses closely related SV
    
    #print "SUBTELCAAAP", subtelcap
    
    tr = U.getPSSetListe(transloc)
    ins1 = U.getPSSetListe(insert1)
    ins2 = U.getPSSetListe(insert2)
    
    inter = set(tr+ins1+ins2)
    subtelcap = updateSubtelcap(subtelcap, inter, ps_min_tn, xsome) #remove dup and inter conta

    return transloc, insert1, insert2, subtelcap

#--------------------------------------------------------------------------
def Etape4(xsome, SV, homo_same, pairopp, pairdiff, ori, d, minPS, min_MinPS):
    """Rempli la liste des SV en rangeant toujours la liste dans l'ordre:
        "groupe ++ puis groupe --" ou "groupe +- puis groupe -+
        1ere liste: PS dans un sens - 2eme liste: PS dans l'autre sens
    """
    SV[ori] = []
    for g in homo_same:

        #Verifie s'il y a des PS dans autres orientations
        none = True
        Fictive = False
        for i in g:
            if i in pairdiff:
                none = none and (not pairdiff[i])

            if U.ps_fictive(xsome[i]):
                print "Etap4 first ", xsome[i][0],i, pairdiff, none
                Fictive = True


        #listes des nouveaux groupes de PS a rajouter ensemble dans SV
        g_homo = [] #celles comprises dans g
        new = []    #celles dans le sens oppose de g


        #Si pas de PS dans autre orientation, ajoute le groupe de PS homogene
        # g a la liste des SV avec
        if none:
            if g:
                g.sort()
            if Fictive:
                print "Etap4 second 1 seul groupe ", xsome[g[0]][0],g
            SV[ori].append([g, [[]]])
        else:
             # S'il y a des PS dans autre orientation, pour chaque PS de g,
             # dresse les groupes homogenes de PS orientation opposee
             # donc dictionnaire g_homo_opp: cle = PS de g, valeur = liste des
             # PS homogenes opposees
            g_homo_opp = homoOpp(xsome, g, pairopp, d, homoOpp)
            if Fictive:
                print "Etap4 second 2 groupes ", xsome[g[0]][0],g
                print "Etap4 second 2 groupes g_homo_opp", g_homo_opp

            #listes des nouveaux groupes de PS a rajouter ensemble dans SV
#            g_homo = [] #celles comprises dans g
#            new = []    #celles dans le sens oppose de g

            #liste des PS de g avec des groupes ds g_homo_opp
            lisg = g_homo_opp.keys()
            lisnul = [l for l in g if l not in lisg] #ps de g sans ps opposees forme subtelo_cap/transloc non recip
            if lisnul :
                #SV.append([lisnul, []])
                SV[ori].append([lisnul, [[]]])
            lisg.sort()
            Lh = len(lisg)
            if Fictive:
                print "Etap4 trois.a 2 groupes ", xsome[g[0]][0],g, lisg, Lh
            # s'il existe au moins 2 PS de g avec des groupes g_homo_opp
            if Lh > 1:
                # Regarde si les PS copines signe opposes sont partagees par
                # d'autres PS de meme signe
                for j in range(Lh-1):
                    for k in range(j+1, Lh):
                        inter_kj = [l for l in g_homo_opp[lisg[k]] if\
                                 l in g_homo_opp[lisg[j]] ]
                        if Fictive:
                            print "Etap4 trois.b 2 groupes ", g[0], \
                            xsome[g[0]][0], lisg[k], lisg[j], "* inter ",\
                            inter_kj
                            print "Etap4 trois.b.a 2 groupes ", xsome[g[0]][0],\
                            g_homo_opp[lisg[k]], g_homo_opp[lisg[j]]
                        if inter_kj:

                        # S'il existe au moins 2 PS de g avec des PS opposes
                        # regarde si elles ne font pas partie d'un groupe deja definit
                        # dans new

                        #C'est là que ça merde : Il faut rajouter toutes les PS de g qui ont des PS opp communes
                        #or ici ça ne marche pas
                            inter_kj.sort()
#                            print "sumCom ", g_homo, lisg[k], lisg[j], new, inter_kj
                            already = sumCom(g_homo, new, lisg[k], lisg[j], inter_kj)

                            if not already:
                                g_homo.append([lisg[k], lisg[j]])
                                new.append(inter_kj)
                            if Fictive:
                                print "After sumCom ",xsome[g[0]][0], g_homo, "*",inter_kj,"*already",already,"*new ",new


                if Fictive:
                    print "Etap4 trois.e 2 groupes ", xsome[g[0]][0],g_homo,new
            elif Lh == 1:
                #s'il existe 1 seul PS de g avec des groupes g_homo_opp
                g_homo = [lisg]
                new = [g_homo_opp.values()[0]]
                if Fictive:
                    print "Etap4 trois.f 1 PS g / 1 PS g_opp", lisg, new[0]
#                if len(pp) > 0:
#                    new.append([o, opp])
            else:
                #si aucun PS de g n'a de groupes g_homo_opp
                if len(g) > 1:
                    g.sort()
                    SV[ori].append([g, [[]]])


            nb_min_ps = U.getMinPSvalue(minPS)
            for c in range(len(new)):
                # pour tous les groupes formes de g_homo et new
                # (groupes de PS associes formant SV)
                #if ori == "++" or ori == "+-":
                if len(new[c]) + len(g_homo) >= nb_min_ps:
                    new[c].sort()
                    g_homo[c].sort()
                    if Fictive:
                        #print "Etap4 2 groupes compatibles ", xsome[new[c][0]][0],g_homo,new
                        print "ETap4 Boucle new ",c, g_homo[c], new[c]

                    orip = oriPair(g_homo[c][0],new[c][0][0], xsome)
                    if not sumSV( SV,  g_homo[c], new[c][0], orip):
                        SV[ori].append([g_homo[c], new[c]])

                    if Fictive:
                        print "Etap4 last Fictive ", SV[ori][-1]

    #Remove singletons
    cleanSV(SV, min_MinPS)


#--------------------------------------------------------------------------
def ClassifIndiv(sv, sv1, sv2, xsome, centrom, subtelo, SV, d, allowed_transloc,
                 interne, diff_max_geneconv, diff_max_other, diff_min_other,
                 subtelcap, geneconv, transloc, insert1, insert2, snp, outf,
                 minPS, ps_type ):
    #print sv
    Fictive = False
    if U.ps_fictive(xsome[sv1[0]]):
        Fictive = True
#

    nb_min_ps = U.getMinPSvalue(minPS)
    nb1 = len(sv1)
    nb2 = len(sv2)
    nPS = nb1 + nb2
    typesv = "rejected"
    if (nPS >= nb_min_ps):
    #if (nPS >= 1):    
        # Si une seule jonction
        #if (not sv2) and (len(sv1) >= minPS["tn"]):
        if (not sv2):
            if len(sv1) >= minPS["tn"]:
                #print "y en a ", sv
                isSubtel(xsome, centrom, subtelo, interne, sv, diff_max_geneconv,
                         diff_max_other, diff_min_other, subtelcap, geneconv,
                         insert1, snp, d, minPS)
                #print insert1, "====", subtelcap
                if Fictive:
                    print "Classif isNRT", xsome[sv1[0]][0],sv
                typesv = "tn"
        else:

        # Si PS opposes, classe compatible et non compatible avec transloc

            oktransloc = []
            put_ins = []
            #print "RRRRRRRRRRRRR sv", sv, "sv1", sv1, "sv2", sv2
            compat, notcompat = U.TestArmOri(xsome, centrom, allowed_transloc,
                                                        sv1, sv2)
            #print "\t Armorie:", time.time() - t1

            if Fictive:
                print "Classif Compat ", xsome[sv1[0]][0]," compat ",compat,\
                " not compat ", notcompat, "ps_type:", ps_type

            if compat :
                # Test coordonnees sur meme chromosome pour ok avec transloc
                oktransloc, put_ins = candidatOpposeTransloc(xsome, compat,
                                                             d, ps_type)

                if Fictive:
                    print "Classif Before liste ", compat[0], "*" , compat[1], "oktransloc:", oktransloc, "putins", put_ins, "sv", sv
                    for toto in compat[0]:
                        print "Classif Before list compat 1 ",toto, xsome[toto]
                    for toto in compat[1]:
                        print "Classif Before list compat 2 ",toto, xsome[toto]


            # Si le nb de PS ok pour transloc est egal au nb de PS de la SV,
            # classe en RT
            if sumPSinSV(oktransloc) == nPS :
                if (not put_ins) and (len(oktransloc[0]) >= minPS["tr"]) and \
                   (len(oktransloc[1]) >= minPS["tr"]) :
                    isTransloc(xsome, oktransloc, diff_max_geneconv, diff_max_other,
                           diff_min_other, geneconv, transloc, snp, subtelcap, d,
                           centrom, subtelo, interne)
                    if Fictive:
                        print "Classif isTransloc", xsome[sv1[0]][0], compat, \
                        oktransloc
                    typesv = "tr"
                    ##########################################################################
                    #ELIF BY ALEX: lorsque put_ins n est pas vide et que toutes
                    #les PS passent sumPSinSV ==> il faut garder l'INS ou la TN !!
                elif put_ins and (len(oktransloc[0]) >= minPS["ins"]) and \
                   (len(oktransloc[1]) >= minPS["ins"]) :
                    if Fictive:
                        print "Classif isInsert from tr", xsome[sv1[0]][0], compat, notcompat
                    isInsert(xsome, notcompat+compat, diff_max_geneconv,
                             diff_max_other, diff_min_other, geneconv,
                             insert2, snp, subtelcap, d, centrom, subtelo, interne)
                    typesv = "ins"

                elif put_ins and (put_ins not in oktransloc) and  \
                     ((len(oktransloc[0]) < minPS["tr"]) or (len(oktransloc[1]) < minPS["tr"])) :
                    
                    if U.ps_fictive(xsome[sv[0][0]]):
                        Fictive = True
                    #if 845 in put_ins:
                    #    Fictive = True
                    if Fictive:
                        print "°°°°Classif isSubtel from tr", xsome[sv1[0]][0], "compat", compat, "notcompat", notcompat, "put_ins", put_ins, "oktransloc", oktransloc
                        print "°°°°put_ins", put_ins, "oktransloc", oktransloc, "put_ins not in oktransloc", put_ins not in oktransloc
                    isSubtel(xsome, centrom, subtelo, interne, [put_ins,[]], diff_max_geneconv,
                             diff_max_other, diff_min_other, subtelcap, geneconv,
				       insert1, snp, d, minPS)
                    typesv = "tn"

                    ##########################################################################

            else:
                # Si le nb de PS ok pour transloc est nul ou inferieur au nb
                # de PS de la SV
                # Si 1 jonction n'est supportee que par une seule PS, non-recip transloc
                #print "attention0", sv
                if ((len(sv2) == 1) or (len(sv1) == 1)) and (nPS >= minPS["tn"]) :
                     #print "attention:", sv
                     if Fictive:
                        print "Classif isPutativIns", xsome[sv1[0]][0],compat, notcompat
                     isPutativIns(xsome, notcompat+compat,  diff_max_geneconv,
                                 diff_max_other, diff_min_other, geneconv,
                                 insert1, snp, subtelcap, d, centrom, subtelo, interne)
                     typesv = "ins"


                elif (len(sv1) >= minPS["ins"]) and (len(sv2) >= minPS["ins"]):
                    if Fictive:
                        print "Classif isInsert", xsome[sv1[0]][0], compat, notcompat
                    isInsert(xsome, notcompat+compat, diff_max_geneconv,
                             diff_max_other, diff_min_other, geneconv,
                             insert2, snp, subtelcap, d, centrom, subtelo, interne)
                    typesv = "ins"


#    if typesv != "rejected":
#        print "ClassifIndiv ", typesv, sv


    return nPS, typesv

#--------------------------------------------------------------------------
def Classif(xsome, centrom, subtelo, SV, d, allowed_transloc, interne,
            diff_max_geneconv, diff_max_other, diff_min_other, subtelcap,
            geneconv, transloc, snp, insert1, insert2, outf, minPS, ps_type):
    """ Classify SV."""
    ntot = 0
    #llprint SV
    for sv in SV:
        #print ntot, len(SV)
        ntot += 1
        #TEST TO BE REMOVED
        #if ntot == 1000000000 :
        #    break
        #if 449 in sv[0] or 450 in sv[0]:
        #    print "Classif pre0 ",sv, xsome[sv[0][0]][0], d


        if U.ps_fictive(xsome[sv[0][0]]):
            print "Classif 0 ",sv, xsome[sv[0][0]][0]
        if not sv[1]:
            #print "ntot", ntot
            nPS, typesv = ClassifIndiv(sv, sv[0], sv[1], xsome, centrom, subtelo,
                     SV, d, allowed_transloc, interne, diff_max_geneconv,
                     diff_max_other, diff_min_other, subtelcap, geneconv,
                     transloc, insert1, insert2, snp, outf, minPS, ps_type)
            if U.ps_fictive(xsome[sv[0][0]]):
                print ntot, "Classif 1", nPS, typesv, sv[0]

        else:
            #print "ntot", ntot
            #print "insert1", insert1
            for sv2 in sv[1]:
                if U.ps_fictive(xsome[sv[0][0]]):
                    print "Classif 2 ", sv[0], sv2

                nPS, typesv = ClassifIndiv(sv, sv[0], sv2, xsome, centrom, subtelo,
                     SV, d, allowed_transloc, interne, diff_max_geneconv,
                     diff_max_other, diff_min_other, subtelcap, geneconv,
                     transloc, insert1, insert2, snp, outf, minPS, ps_type)
                     
#    tr = U.getPSSetListe(transloc)
#    ins1 = U.getPSSetListe(insert1)
#    ins2 = U.getPSSetListe(insert2)
#    inter = set(tr+ins1+ins2)
#    #temp = tr+ins1+ins2
#    #print "jjj", temp[0]
#    subtelcap = updateSubtelcap(subtelcap, inter, minPS["tn"], xsome) #remove dup and inter conta

#--------------------------------------------------------------------------
def extremites(liste, sens):

    if sens == "max":
        return max(liste)
    else:
        return min(liste)

#--------------------------------------------------------------------------
def bordersTransloc(xsome, l1, l2, cen_str1, cen_str2) :

    combi = {("L-L+", "L+L-"):["max","L-","min","L+","min","L+","max","L-"],
         ("L-R-", "L+R+"):["min","L+", "max","L-","max","R-","min","R+"],
         ("R-R+", "R+R-"):["max","R-","min","R+","max","R+","max","R-"],
         ("R-L-", "R+L+"):["max","R-","min","R+", "max","L-","min","L+"]}


    coordsign = {}
    borders = {}
    chr1 = xsome[l1[0]][2]
    chr2 = xsome[l1[0]][5]

    if l1 != l2 and cen_str1 != cen_str2:
        signe = []
        comb = []
        for c in combi:
#            print "COMBI ", cen_str1, cen_str2, c
            if cen_str1 in c and cen_str2 in c:
                for elem in (1,3,5,7):
                    signe.append(combi[c][elem])
                for elem in (0,2,4,6):
                    comb.append(combi[c][elem])
        for i in l2:
            ps = xsome[i]
            U.add2Coord(coordsign, (chr1, cen_str2[:2]), ps[3])
            U.add2Coord(coordsign, (chr2, cen_str2[2:]),  ps[6])

#        print xsome[l2[0]]


    else: #Pour TN
        signe = [cen_str1[:2], cen_str1[:2], cen_str1[2:], cen_str1[2:]]
        comb = ["min", "min", "max", "max"]


    for i in l1:
        ps = xsome[i]
        #cle = chromosome, orientation sur le chromosome
        #valeur = coordonnee sur le chromosome

        U.add2Coord(coordsign, (chr1, cen_str1[:2]), ps[3])
        U.add2Coord(coordsign, (chr2, cen_str1[2:]),  ps[6])



#    for c in coordsign:
#        print c, coordsign[c]
#
#    print "SIGNE TRANSLOC", signe
#    print chr1, signe[0]
#    print chr1, signe[1]
#    print chr2, signe[2]
#    print chr2, signe[3]
    borders[(chr1, signe[0])] = extremites(coordsign[(chr1, signe[0])], comb[0])
    borders[(chr1, signe[1])] = extremites(coordsign[(chr1, signe[1])], comb[1])
    borders[(chr2, signe[2])] = extremites(coordsign[(chr2, signe[2])], comb[2])
    borders[(chr2, signe[3])] = extremites(coordsign[(chr2, signe[3])], comb[3])

    return borders
#------------------------------------------------------------------------------
def extrInv(xsome, liste):
    """determine le min et max des coordonnees"""

    minleft1 = xsome[liste[0]][3]
    maxleft1 = xsome[liste[0]][3]
    maxright2 = xsome[liste[0]][6]
    minright2 = xsome[liste[0]][6]

    for i in range(1, len(liste)):
        if xsome[liste[i]][3] <= minleft1 :
            minleft1 = xsome[liste[i]][3]
        elif xsome[liste[i]][3] > maxleft1:
            maxleft1 = xsome[liste[i]][3]

        if xsome[liste[i]][6] > maxright2 :
            maxright2 = xsome[liste[i]][6]
        if xsome[liste[i]][6] <= minright2 :
            minright2 = xsome[liste[i]][6]

    return minleft1, maxleft1, minright2, maxright2

#--------------------------------------------------------------------------
def getBorders(xsome, l1, l2, cen_str1, cen_str2, d, typesv):

    ### ATTENTION LE GETBORDERS NE MARCHE QUE POUR TRANSLOC ET INSERTIONS
    ### QUID DES SNP GENECONVERSIONS

    borders = {}
    chr1 = xsome[l1[0]][2]
    chr2 = xsome[l1[0]][5]

#    print "getborders ",cen_str1, cen_str2
#    print "getborders ", chr1, cen_str1[:2], cen_str2[:2]
#    print "getborders ", chr2, cen_str1[2:], cen_str2[2:]

    borders[(chr1, cen_str1[:2])] = 0
    borders[(chr1,cen_str2[:2] )] = 0
    borders[(chr2, cen_str2[2:])] = 0
    borders[(chr2, cen_str1[2:])] = 0


#Lists of signs and min/max pour TR ou INS
    #print "typeSVS", typesv
    if "t" in typesv :
        borders = bordersTransloc(xsome, l1, l2, cen_str1, cen_str2)
    elif "ins" in typesv :

        signe=[xsome[l1[0]][1], xsome[l1[0]][4],xsome[l2[0]][1],
               xsome[l2[0]][4] ]

        minleft1a, maxleft1a, minright2a, maxright2a  = extrInv(xsome, l1)
        minleft1b, maxleft1b, minright2b, maxright2b  = extrInv(xsome, l2)

        startchr1_left = min(minleft1a, minleft1b)
        endchr1_left = max(maxleft1a, maxleft1b)
        startchr2_left = min(minright2a, minright2b)
        endchr2_left = max(maxright2a, maxright2b)

        if endchr1_left - startchr1_left < 2*d:
            #Sur le chromosome hôte qui reçoit l'INS (celui ou les coord respectent 2d)
#            borders[(chr1, signe[0])] = endchr1_left
#            borders[(chr1, signe[2])] = startchr1_left
            borders[(chr1, cen_str1[:2])] = endchr1_left
            borders[(chr1, cen_str2[:2])] = startchr1_left
            #Pour le fragment inséré (ne respecte pas 2d ou a des coord chevauchantes)
            if xsome[l2[0]][1]+xsome[l2[0]][4] == "--" or  xsome[l2[0]][1]+xsome[l2[0]][4] == "++":
                #Si --/++
                borders[(chr2, cen_str1[2:])] = endchr2_left
                borders[(chr2, cen_str2[2:])] = startchr2_left
            else:
                #Si -+/+-
                borders[(chr2, cen_str1[2:])]  = startchr2_left
                borders[(chr2, cen_str2[2:])] = endchr2_left

        else:
            borders[(chr2, cen_str1[2:])] = endchr2_left
            borders[(chr2, cen_str2[2:])] = startchr2_left
            if signe[0]+signe[1] == "--" or  signe[0]+signe[1]== "++":
                borders[(chr1, cen_str1[:2])] = endchr1_left
                borders[(chr1,cen_str2[:2] )] = startchr1_left
            else:
                borders[(chr1, cen_str1[:2])] = startchr1_left
                borders[(chr1, cen_str2[:2])] = endchr1_left


    return borders
#--------------------------------------------------------------------------
def define_SV_borders(sv, sv1, sv2, xsome, typesv, d):
    """Define SV borders and add borders on A and borders on B at the
    end of the sv list"""
    #print "SV_BOARDERS", sv, sv1, sv2, typesv
    ### ATTENTION LE GETBORDERS NE MARCHE QUE POUR TRANSLOC ET INSERTIONS
    ### QUID DES SNP GENECONVERSIONS

    nb2 = len(sv2)
    cen_strA = sv[2]
    cen_strB = sv[3]
    #print "ddd", sv, typesv
    chr1 = xsome[sv[0][0]][2]
    chr2 = xsome[sv[0][0]][5]
#    print "Define SV border ", xsome[sv[0][0]], sv[0], sv1
#    print "Define SV border ", sv[1], sv2

    signe = [xsome[sv[0][0]][1], xsome[sv[0][0]][4]]

    if nb2 > 0:
        signe = signe + [ xsome[sv[1][0]][1], xsome[sv[1][0]][4] ]
        borders = getBorders(xsome, sv[0], sv[1], cen_strA,
                             cen_strB, d, typesv)
        left1 = borders[(chr1, cen_strA[:2])]
        right1 = borders[(chr1, cen_strB[:2])]
        left2 = borders[(chr2, cen_strA[2:])]
        right2 = borders[(chr2, cen_strB[2:])]

    else:
        signe = signe + [ xsome[sv[0][0]][1], xsome[sv[0][0]][4] ]
        borders = getBorders(xsome, sv[0], sv[0], cen_strA,
                             cen_strA, d, typesv)
        left1 = borders[(chr1, cen_strA[:2])]
        right1 = borders[(chr1, cen_strA[:2])]
        left2 = borders[(chr2, cen_strA[2:])]
        right2 = borders[(chr2, cen_strA[2:])]

    sv.append([left1, right1, left2, right2])

#--------------------------------------------------------------------------
def retrieve_borders_from_sv(sv):
    """ return borders left A, rightA, leftB, right B"""
    #print "totototototo", sv[-1][0], sv[-1][1], sv[-1][2], sv[-1][3]
    return sv[-1][0], sv[-1][1], sv[-1][2], sv[-1][3]
#--------------------------------------------------------------------------
def overlap_SV(xsome, ori, listsv, d):
    """Group the SV with overlapping borders on A and overlapping borders on B
    """

    ### ATTENTION LE GETBORDERS NE MARCHE QUE POUR TRANSLOC ET INSERTIONS
    ### QUID DES SNP GENECONVERSIONS

    new = []
    i = 0
    stop = False
    L = len(listsv[ori])
    toremove = []
    while i < L-1 and not stop :
        if i not in toremove:
            svi = listsv[ori][i]
            merged = [svi]
            svileft1, sviright1, svileft2, sviright2 = retrieve_borders_from_sv(svi)
            flag = 0

            #group the SV overlapping svi
            j = i + 1
            while j < L :

                if j not in toremove:
                    svj = listsv[ori][j]
                    svjleft1, svjright1, svjleft2, svjright2 = retrieve_borders_from_sv(svj)
                    coo_left1_ok = U.parameter_ok(svileft1, svjleft1, d)
                    coo_right1_ok = U.parameter_ok(sviright1, svjright1, d)
                    coo_left2_ok = U.parameter_ok(svileft2, svjleft2, d)
                    coo_right2_ok = U.parameter_ok(sviright2, svjright2, d)
                    if coo_left1_ok and coo_right1_ok and coo_left2_ok and coo_right2_ok:
                        merged.append(svj)
                        #rajoute i a toremove pour ne pas mettre svi
                        if flag == 0:
                            toremove.append(i)
                        #flag à j pour repartir de la 1ere sv non suivante
                        flag = 1
                        #rajoute j a toremove pour ne pas mettre svj
                        toremove.append(j)

                j = j + 1

            #add the merged group of SV to the new list
            if len(merged) > 1:
                new.append(merge_SV(merged))

        i = i + 1


    #Add the independant initial SV

    for i, sv in enumerate(listsv[ori]):
        if i not in toremove :
            new.append(sv)

    return new
#--------------------------------------------------------------------------
def merge_SV(liste, xsome, d):
    """define a new SV with the common PS of the sv in the liste"""
    #intersection of PS from the first groups of PS (sv[0])
    sv0 = liste[0]
    s = set(liste[0][0])
    for j in liste[1:]:
        s = s.intersection(j[0])
    gps1 = list(s)

    s = set(liste[0][1])
    for j in liste[1:]:
        s = s.intersection(j[1])
    gps2 = list(s)

    #reconstruit la sv avec les deux nouveaux groupes de PS
    #vire les borders et la p-avlue
    mergedsv = [gps1, gps2] + sv0[2:-2]

    define_SV_borders(mergedsv, gps1, gps2, xsome, d)
    return mergedsv

#--------------------------------------------------------------------------
def remove_Reuse(liste, reuse, xsome):
    """Remove SV with >= 1/3 PS reused"""

    nbreused = {}
    for ori in liste:
        nbreused[ori]={}

    for ps in reuse:
        if reuse[ps] > 1:
            for ori in liste:
                for i, sv in enumerate(liste[ori]):
                    if (ps in sv[0]) or (sv[1] and ps in sv[1]):
                        nbreused[i] = nbreused.get(i, 0) + 1

    removedsv = 0
    for ori in liste:
        for i, sv in enumerate(liste[ori]):
            nbPS = len(sv[0]) + len(sv[1])
            if (i in nbreused) and (nbreused[i] >= nbPS/3.0):
                removedsv += 1
                del liste[ori][i]


#--------------------------------------------------------------------------
def reuse(liste, xsome):
    """count the number of reuse of each PS"""
    reuse = {}
    for ori in liste:
        for sv in liste[ori]:
            for ps in sv[0]:
                reuse[ps] = reuse.get(ps, 0)  + 1
            if sv[1]:
                for ps in sv[1]:
                    reuse[ps] = reuse.get(ps, 0) + 1
    remove_Reuse(liste, reuse, xsome)

#--------------------------------------------------------------------------
def Affiche(xsome, centrom, liste, pair, fich, typesv, nb):
    """Write SV to output files."""

    ### ATTENTION LE GETBORDERS NE MARCHE QUE POUR TRANSLOC ET INSERTIONS
    ### QUID DES SNP GENECONVERSIONS

    out = open(fich+"_"+typesv+"_byRP.csv", "a")
    outp = open(fich+"_"+typesv+"_bySV.csv", "a")
    for ori in liste:
        for sv in liste[ori]:
            nb1 = len(sv[0])
            nb2 = len(sv[1])
            nbPS = nb1 + nb2

            cen_strA = sv[2]
            cen_strB = sv[3]
            pbal = str(sv[8])
           
            #print "---", typesv, sv
            left1, right1, left2, right2 = retrieve_borders_from_sv(sv)
            if typesv == "non_reciprocal_translocations":
                right1 = -1
                left2 = -1

            nb += 1
            manip = os.path.split(fich)[1].split("_")[0]
            U.write_bySV(outp, (manip, pair, nb, nbPS, nb1, nb2, left1,
                                right1, abs(right1-left1),cen_strA, left2, right2,
                                abs(right2 - left2), cen_strB, -1, -1, pbal, "NA"))


            for s in sv[0]:
                si = xsome[s]
                U.write_byPS(out, (manip, pair, nb, nbPS, nb1, nb2, left1,
                                   right1, abs(right1-left1), cen_strA, left2,
                                   right2, abs(right2-left2), cen_strB, -1, -1,
                                   pbal, "NA", si[0], si[1], si[2], si[3], si[4],
                                   si[5], si[6], 0))


                if sv[2] != U.Centro_Inter(centrom, si[2], si[1],
                            si[3], si[5], si[4], si[6]):
                    #print "PB Affiche ", pair, fich, nb
                    toto = 1
            for s in sv[1]:
                si = xsome[s]
                U.write_byPS(out, (manip, pair, nb, nbPS, nb1, nb2, left1,
                                   right1, abs(right1-left1), cen_strA, left2,
                                   right2, abs(right2-left2), cen_strB, -1, -1,
                                   pbal, "NA", si[0], si[1], si[2], si[3], si[4],
                                   si[5], si[6], 0))



    out.close()
    outp.close()

    return nb
