# -*- coding: utf-8 -*-
#__doc__       = """ texte doc """

import sys
import os
import gc
dirprog = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.append(dirprog)
import sv_inter as sv
import Ulysse_utils as U
import Ulysse_stats as Ualex
import datetime


#import linecache

#--------------------------------------------------------------------------
def runDetectionInter(params, stats, chrDicos):

    outf = params["out"]
    #Il manque le chargement du reformat_dist

#    minPS = int(params["min_nb_ps"])
    #for p in params:
    #    print p, params[p]


    ############################Cree les fichiers de sorties
    outs = ["non_reciprocal_translocations", "reciprocal_translocations",
           "insertions"]
#    if stats["median"] > 500:
#        outs = outs + ["snp", "geneconversions"]

    for fico in outs:
        U.create_files(fico, outf)

    subtelo=U.get_subtelo_limits(params)
    centrom=U.get_centro_coords(params)
    #print "++++++", params

############################Define parameters
    #d = mediane + 3 * MAD
    d = stats["median"] + 3 * stats["mad"]
    #ecart_max entre coordonnees sur meme chromosome pour gene conversion
#    print "d Inter ", d
    #Detection of geneconversions and snp only possible with mate-pairs library
    diff_max_geneconv = 0
    diff_max_other = 0
    diff_min_other = 0
    if stats["median"] > 500 and (2 * stats["mad"] <  stats["median"]):
    	diff_max_geneconv = 1000 #longueur max tract de conversion genique en meiose
    	diff_max_other = (2 * stats["median"]) + (3 * stats["mad"])
    	diff_min_other = (2 * stats["median"]) - (3 * stats ["mad"])
#    print "Detect Inter ", stats["median"], diff_max_geneconv, diff_max_other, diff_min_other

    #PS paired by opposite orientations
    opposite = {"++":"--", "--":"++", "+-":"-+", "-+":"+-"}


    if stats["ps_type"] == "MP":
        # Combinaisons autorisees position centromere/orientation read des PS
        # formant une transloc:
        allowed_transloc = [("L-L+", "L+L-"), ("L-R-", "L+R+"), ("R-R+", "R+R-"),
                            ("R-L-", "R+L+")]

        # Definition des regions les plus internes pour les subtelomeric
        # translocations 1er champ:  signe combine pour un PS:
        # bras-orientation_read sur chrX ,  bras-orientation_read sur chrY,
        # chr le plus interne,  formule
        interne = [("L-L+", 2, "L+"), ("R-R+", 1, "R-"), ("L-R-", 2, "R-"),
                   ("R-L-", 1, "R-"), ("R+L+", 2, "L+"), ("L+R+", 1, "L+"),
                   ("L+L-", 1, "L+"), ("R+R-", 2, "R-")]
    else:
        allowed_transloc = [("L+L-", "L-L+"), ("L+R+", "L-R-"), ("R+R-", "R-R+"),
                            ("R+L+", "R-L-"), ("NiNi", "NiNi")]
        interne = [("L+L-", 2, "L-"), ("R+R-", 1, "R+"), ("L+R+", 2, "R+"),
                   ("R+L+", 1, "R+"), ("R-L-", 2, "L-"), ("L-R-", 1, "L-"),
                   ("L-L+", 1, "L-"), ("R-R+", 2, "R+")]


    with open(params["out"]+".INTER.report.out", "w") as fileo:
        fileo.write("\n-----------\tDetection of Inter Chromosomal SV\n")
        fileo.write("\t\tMaximum distance between PS extremities (d) = %.2f\n\n"\
        % (d))
        if stats["median"] > 500 and (2 * stats["mad"] <  stats["median"]):
            fileo.write("\tLibrary of mate-pairs : detection of SNP and gene\
conversion is on\n\n")


    clasdifx = []
    dicQual = sv.ReadFilesBAM(params, clasdifx, stats["ps_type"])
    print "PS loaded from BAM file"
    os.system("date")

#TEMPORARY TO BE REMOVED
#    names = sv.ReadFilesPaired(params["in"], clasdifx)


    Xsome = {}
    sv.Filter(clasdifx, Xsome, params["range"])
    clasdifx = []

#    for pair, PSS in Xsome.iteritems():
#        #print "pair", pair
#        clasx = list(set([tuple(i) for i in PSS]))
#        Xsome[pair]=clasx


# START ALEX PARAMETERS #######################################################
#   numberOfInterChromosomalPSFile = params["in"]+"_PScandidate_per_Xsome_pairs.csv"
#    list_chr_names = "SEP".join(params["range"])
    #list_chr_length = "SEP".join([ str(chrDicos[x]) for x in params["range"] ])
    list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    list_chr_names = "SEP".join(chrDicos.keys())
    #print list_chr_names
    #print
    #print list_chr_length

    #list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    #list_chr_names = "SEP".join(chrDicos.keys())

    p = U.getcommonstart(params["range"])
    prefix = ""
    for i in p:
        if i.isalpha():
            prefix += i

    if prefix != "" and centrom and not centrom[0][0].startswith(prefix):
        for c in centrom:
            c[0] = prefix+c[0]


    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + float(params["n"]) * stats["mad"]
    sval = stats["median"] + float(3) * stats["mad"]

    #TO REMOVE
    #ps_min_ins, ps_min_tr, ps_min_tn = 1,1,2


    #Determine le nb minimum de ps pour chaque type de SV
    detectionFile = "ClusterLimitSize"
    ps_min_ins, ps_min_tr, ps_min_tn, infoMsg = Ualex.define_MinPS("INTER",
        stats, params["nsv"], list_chr_length, list_chr_names, detectionFile,
        params["in"], params["in"]+".dist.table", sval, params["fdr"], params["out"], stats["rl"],
                                 params["n"])
# END ALEX PARAMETERS #######################################################
    #ps_min_ins = 10
    with open(params["out"]+".INTER.report.out","a") as fileo:
        fileo.write("Detection processed\n")
        fileo.write("results written to "+params["out"]+"_[NRT/INS/RT]_by[PS/SV].[stats.]csv\n")
        fileo.write("\t\tMinimum number of PS to define an Insertion = %d\n\n"\
        % (ps_min_ins))
        fileo.write("\t\tMinimum number of PS to define a reciprocal translocation = %d\n\n"\
        % (ps_min_ins))
        fileo.write("\t\tMinimum number of PS to define a non reciprocal translocation = %d\n\n"\
        % (ps_min_ins))

###################### START DETECTION
    #print Xsome
    for pair in Xsome:
    #for pair in ["chr11-chr5"]:
        #ori_pair = U.get_X_original(pair, names)
        os.system("date")
        print "*** Processing", pair, len(Xsome[pair]), "***"
        candidats = {"--":[], "++":[], "+-":[], "-+":[]}
        #print candidats
        #Classe les ps en fonction de leur orientation
        print "    Starting Etape 1: PS sorting by orientation", datetime.datetime.now()
        #sv.Etape1(Xsome[pair], candidats)
        order1, order2 = sv.Etape1_quick(Xsome[pair], candidats)
        nb = 0

        SV = {}
        subtelcap = {}
        geneconv = {}
        transloc = {}
        insert1 = {}
        insert2 = {}
        snp = {}

        #TODO: temporary ???                
        minPS = {"ins" : ps_min_ins, "tn" : ps_min_tn, "tr" : ps_min_tr}
        if int(stats["median"])<1000:
            minPS = {"ins" : 1, "tn" : 2, "tr" : 1}
        #minPS = {"ins" : 2, "tn" : 2, "tr" : 2}
        #min_MinPS = min(ps_min_ins, ps_min_tn, ps_min_tr)
        
        min_MinPS = min(minPS.values())
        #min_MinPS = 3
        #print "minPS is equal to", minPS

        print "    Starting Detection", datetime.datetime.now()
        for oris in candidats:
        #for oris in ["-+"]:
            print "\tOrientation \""+oris+"\""
            # Pour chaque pair orientation "--" ou "+-" recupere les PS
            # chevauchants dans une autre orientation et classe a part les
            # copains chevauchants dans orientation opposee qui vont bien
            # pour subtelo_cap, insertions (pour transloc il faut refiltrer
            # dans candidatOpposeTransloc).
            # Si pas de PS opposees qui vont bien,  met liste vide

#            for jj in candidats[oris]:
#                print "BeforeEtape2 ", pair, Xsome[pair][jj]


            #Fait des groupes homogene avec les PS orientation opposes a oris
            print "\t    Starting Etape 2: Making groups of PS in the same orientation", datetime.datetime.now()
            homo_same, candidats2 = sv.Etape3_quick(Xsome[pair], candidats[oris], d, min_MinPS, oris)

            
                
            print "\t    Starting Etape 3: Individual opposite PS collection", datetime.datetime.now()
            pairopp, pairopp2, pairOppCoord = sv.Etape2_quick(order1, order2, candidats2, \
                                          oris, d, opposite, stats["ps_type"],\
                                          min_MinPS)


            #
            print "\t    Starting Etape 4: Makes the list of potential SVs from" , len(homo_same), "left groups", datetime.datetime.now()
            #Rempli la liste des SV en rangeant toujours la liste dans
            #l'ordre: "groupe ++ puis groupe --" ou "groupe +- puis groupe -+"

            #appending a big dictionnary (SV) is very very slow
            #so we creat an intermidiate empty one (SVtemp) and then merge them
            SVtemp={}
            #--------------
#            import json
#            dats = [Xsome[pair], SVtemp, homo_same, pairopp, oris, d, minPS, min_MinPS, pairOppCoord]
#            nams = ["Xsome[pair]", "SVtemp", "homo_same", "pairopp", "oris", "d", "minPS", "min_MinPS", "pairOppCoord"]
#            
#            for i,out in enumerate(nams):
#                with open(out+'.pkl', 'w') as outfile:
#                    json.dump(dats[i], outfile)


            #-------
            #print pairopp
            sv.Etape4_quick(Xsome[pair], SVtemp, homo_same, pairopp, pairopp, oris, d,
                      minPS, min_MinPS, pairOppCoord)
            
            SV = dict(SV.items() + SVtemp.items())
            #print SV, d

            
            #print "debut fusion avec update", datetime.datetime.now()
            #SV.update(SVtemp)
            #print "fin fusion avec udate", datetime.datetime.now()
            ###import pickle
            ###output = open('dict_2.pkl', 'wb')
            ###pickle.dump(SV, output)
            ###output.close()

            ###pkl_f = open('dict_2.pkl', 'rb')
            ###SV = pickle.load(pkl_f)
            ###pkl_f.close()
     
            print "\t    Starting Etape 5: Classification of SVs", datetime.datetime.now()
            if SV[oris] != []:
                sv.Classif(Xsome[pair], centrom, subtelo, SV[oris], d,
                           allowed_transloc, interne, diff_max_geneconv,
                           diff_max_other, diff_min_other, subtelcap,
                           geneconv, transloc, snp, insert1, insert2, outf,
                           minPS, stats["ps_type"])
                           
#        import pickle
#        output = open('transloc-2.pkl', 'wb')
#        pickle.dump(transloc, output)
#        output.close()
#
#        output = open('subtelcap-2.pkl', 'wb')
#        pickle.dump(subtelcap, output)
#        output.close()
#        
#        output = open('insert1-2.pkl', 'wb')
#        pickle.dump(insert1, output)
#        output.close()
#        
#        output = open('insert2-2.pkl', 'wb')
#        pickle.dump(insert2, output)
#        output.close()


#        output = open('transloc-2.pkl', 'rb')
#        transloc = pickle.load(output)
#        output.close()
#        
#        output = open('subtelcap-2.pkl', 'rb')
#        subtelcap = pickle.load(output)
#        output.close()
#        
#        output = open('insert1-2.pkl', 'rb')
#        insert1 = pickle.load(output)
#        output.close()
#        
#        output = open('insert2-2.pkl', 'rb')
#        insert2 = pickle.load(output)
#        output.close()        
#        
#        print insert2

        #print transloc
        print "    Starting Etape 6: Remove SV duplicates", datetime.datetime.now()
        transloc, insert1, insert2, subtelcap = sv.Etape6_quick(transloc, insert1, insert2, subtelcap, ps_min_tn, Xsome[pair])
        #print "CHECK2", l, "VS", len(subtelcap)
       

        #print "TN", subtelcap
        #print "TR", transloc

        print "    Start affichage", datetime.datetime.now()
#        sv.reuse(transloc, Xsome[pair])
        nb = sv.Affiche(Xsome[pair], centrom, transloc, pair, outf,
                        "reciprocal_translocations", nb)
#        sv.reuse(subtelcap, Xsome[pair])
        nb = sv.Affiche(Xsome[pair], centrom, subtelcap, pair, outf,
                        "non_reciprocal_translocations", nb)
#        sv.reuse(insert1, Xsome[pair])
        nb = sv.Affiche(Xsome[pair], centrom, insert1, pair, outf,
                        "insertions", nb)
#        sv.reuse(insert2, Xsome[pair])
        nb = sv.Affiche(Xsome[pair], centrom, insert2, pair, outf,
                        "insertions", nb)

#        if stats["median"] > 500:
#            sv.reuse(snp, Xsome[pair])
#            nb = sv.Affiche(Xsome[pair], centrom, snp, pair, outf, "snp", nb)

#            sv.reuse(snp, Xsome[pair])
#            nb = sv.Affiche(Xsome[pair], centrom, geneconv, pair, outf,
#                            "geneconversions", nb)

#Rajoute les stats aux fichiers resultats
    del Xsome
    Xsome = None
    gc.collect()
    
    
    #Do a final dup removal
    U.cleanAnySVFilePair(params["out"]+"_insertions_bySV.csv", params["out"]+"_insertions_byRP.csv")
    U.cleanAnySVFilePair(params["out"]+"_reciprocal_translocations_bySV.csv", params["out"]+"_reciprocal_translocations_byRP.csv")
    U.cleanAnySVFilePair(params["out"]+"_non_reciprocal_translocations_bySV.csv", params["out"]+"_non_reciprocal_translocations_byRP.csv")

    #add qualities
    U.addMeanSVQuality(params["out"]+"_insertions_bySV.csv", params["out"]+"_insertions_byRP.csv", dicQual)
    U.addMeanSVQuality(params["out"]+"_reciprocal_translocations_bySV.csv", params["out"]+"_reciprocal_translocations_byRP.csv", dicQual)
    U.addMeanSVQuality(params["out"]+"_non_reciprocal_translocations_bySV.csv", params["out"]+"_non_reciprocal_translocations_byRP.csv", dicQual)
    
    #reorder insertions
    U.reOrder(params["out"]+"_insertions_bySV.csv", params["out"]+"_insertions_byRP.csv", d)

    toto, msg, pval_seuil_ins = Ualex.functstats("INS", stats, stats["ninter"], params["nsv"],
                                 list_chr_length, list_chr_names,
                                 params["out"]+"_insertions_bySV.csv",
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"], stats["rl"],
                                 params["n"])
    print msg
    print "\n"
    toto, msg, pval_seuil_tr = Ualex.functstats("TR", stats, stats["ninter"], params["nsv"],
                                 list_chr_length, list_chr_names,
                                 params["out"]+"_reciprocal_translocations_bySV.csv",
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"],stats["rl"],
                                 params["n"])
    print msg
    print "\n"
    toto, msg, pval_seuil_tn = Ualex.functstats("TN", stats, stats["ninter"], params["nsv"],
                                 list_chr_length, list_chr_names,
                                 params["out"]+"_non_reciprocal_translocations_bySV.csv",
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"],stats["rl"],
                                 params["n"])
    print msg
    print "\n"

    print "Detection done\nReport file : "+params["out"]+".INTER.report.out"
    print "Results written to "+params["out"]+"_[NRT/INS/RT]_by[PS/SV].[stats.]csv"

    return pval_seuil_ins, pval_seuil_tr, pval_seuil_tn
    os.system("date")


#------------------------------------------------------------------------------

def runStatsInter(params, stats, chrDicos):
    n = int(params["n"])

# START ALEX PARAMETERS #######################################################
    list_chr_names = "SEP".join(params["range"])
    list_chr_length = "SEP".join([ str(chrDicos[x]) for x in params["range"] ])
    #list_chr_length = "SEP".join([str(x) for x in chrDicos.values()])
    #list_chr_names = "SEP".join(chrDicos.keys())
    
    
    #ATTENTION sval que pour del mais le donner tout le temps
    sval = stats["median"] + n * stats["mad"]
    

    print "Statistics on Interchromosomal SV : Dps = ", stats["ninter"]
    print "False Discovery Rate (fdr) = ", params["fdr"]
    
    toto, msg, pval_seuil_ins = Ualex.functstats("INS", stats, stats["ninter"], params["nsv"],
                                 list_chr_length, list_chr_names,
                                 params["out"]+"_insertions_bySV.csv",
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"], stats["rl"],
                                 params["n"])
    print msg
    toto, msg, pval_seuil_tr = Ualex.functstats("TR", stats, stats["ninter"], params["nsv"],
                                 list_chr_length, list_chr_names,
                                 params["out"]+"_reciprocal_translocations_bySV.csv",
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"],stats["rl"],
                                 params["n"])
    print msg
    toto, msg, pval_seuil_tn = Ualex.functstats("TN", stats, stats["ninter"], params["nsv"],
                                 list_chr_length, list_chr_names,
                                 params["out"]+"_non_reciprocal_translocations_bySV.csv",
                                 params["in"], params["in"]+".dist.table",
                                 sval, params["fdr"], params["out"],stats["rl"],
                                 params["n"])
    print msg
    return pval_seuil_ins, pval_seuil_tr, pval_seuil_tn
                             

#--------------------------------------------------------------------------
def launch(params, stats, chrDicos):

    print "**************** Detection of Inter chromosomal SV"

    if params["only_stats"]:
        pval_seuil_ins, pval_seuil_tr, pval_seuil_tn = runStatsInter(params, stats, chrDicos)
    else:
        pval_seuil_ins, pval_seuil_tr, pval_seuil_tn = runDetectionInter(params, stats, chrDicos)        

    return pval_seuil_ins, pval_seuil_tr, pval_seuil_tn
