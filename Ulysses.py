#!/usr/bin/python
"""
Created on Wed Jul 31 16:58:42 2013

@author: ingrid
"""

import os, sys
import argparse
import deletion
import duplication
import inversions
import detect_inter
import textwrap
import gc
import Ulysse_utils as U

#-------------------------------------------------------------------------
def parser():
    
    """Parse command line argument"""
    svname = {"DUP": "Duplications", "DEL":"Deletions", 
              "INV":"Inversions",
              "INTER":"inter chromosomal variants"}
    
    parser = argparse.ArgumentParser(prog='Ulysses.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
 
 ############## ULYSSES V1.0 Detection of Strural Variations ################                                     

 You can set the parameters by command line arguments, 
 or define them within a parameter file. 
 Otherwise, a default file will be created during the preprocessing step

 ######################################
 
 '''), 
                    epilog=textwrap.dedent('''                  
                    
                                            Enjoy!!!                    
                    
                    '''))
    
    parser.add_argument("-in", metavar='input_file_name', 
                        help='name of the library file (BAM/SAM format)')                        
    
    parser.add_argument("-p", metavar='parameter_file_name', 
                        help='file containing all parameter values. default \
                        ulysses_params', default='ulysses_params')

    parser.add_argument("-n", metavar='mulitplicative factor for MAD and d(n)', type=int,
                        default=6, help='Interdistance and insert size consistency factor')


    parser.add_argument("-stats", metavar='statistics_file_name', 
                        help='name of statistics file created during library processing')

    parser.add_argument("-range", metavar='chromosomes used for detection', 
                        help='chromosome numbers (1-10,X or 1,2,Y) to be screened\
                        or \'all\' for all chromosomes')

    parser.add_argument("-nsv", metavar='Maximum number of expected SV candidate', \
                        help='Necessary to define the minimum\
                        number of PS defining an SV', type=int) #, default=10000)

    parser.add_argument("-fdr", metavar='threshold of deletion FDR', 
                        help='Upper threshold for False discovery rate of a\
                        deletion', type=float, default=0.01)
    
    parser.add_argument("-mapq", metavar='mapping_quality_threshold',
                        type=int, default=20, help='Average minimal mapping quality of reads for DUP, DEL, sINS and INV')

    parser.add_argument("-mapqx", metavar='mapping_quality_threshold',
                        type=int, default=20, help='Average minimal mapping quality of reads for INS, RT and NRT')

    parser.add_argument("-typesv", metavar='type_of_SV_detection',default='all',\
                        help='all for all SV types. DUP for duplication, \
                        DEL for deletion, INV for inversion, INTER for \
                        inter-chromosomal SV, default=all.\
                        Add "--stats" to SV type (e.g. DEL--stats) to launch \
                        the statistical analysis only or use statsmod')

    parser.add_argument("-statsmod", metavar='statistical_tests', type=bool,\
                        default = False, help='If TRUE, only the statistical \
                        validation module is performed on candidate detected SV')

    parser.add_argument("-out", metavar='prefix_of_output_file', 
                        help='if not given, library_file_name_Ulysses')
                        
    parser.add_argument("-vcf", metavar='create VCF 4.2 output file', 
                        help='Name of the sample in VCF file (e.g. Tumor', default = False)                        

    parser.add_argument("-a", "--annotation", metavar='name_of_annotation_file', 
                        help='if not given, no distinction of putative reciporcal translocations',
                        default = "NA")
    parser.add_argument("-field_chr", metavar='column_in annotation_file_for_chromosome number', 
                        help='column containing chromosome number in annotation file', 
                        type = int, default = 1)
 
    parser.add_argument("-field_type", metavar='field_nb_for_element_type', 
                        help='column containing feature type in annotation file',
                        type = int, default = 3)
    parser.add_argument("-field_start", metavar='field_nb_for_element_start', 
                        help='column containing start coordinate in annotation file',
                        type = int, default = 4)
    parser.add_argument("-field_end", metavar='field_nb__for_element_end', 
                        help='column containing end coordinate in annotation file',
                        type = int, default = 5)
    parser.add_argument("-field_sep", metavar='field_separator_in_annotation_file', 
                        help='field separator in annotation file', default = "tb")

                    
                        
    
    args = parser.parse_args()
    
    
    print "\n\n\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n"
    print "\t\t\t\t---  ULYSSES v1.0 20140827  ---     \n\n"
    print "   Ulysses: Accurate detection of rare structural variations from high coverage genome sequencing\n"
    print "\t\t Alexandre Gillet, Hugues Richard, Gilles Fischer and Ingrid Lafontaine\n"
    print "\t\t\t http://www.lcqb.upmc.fr/ulysses\n\n"
    print "\t\t\t\t Copyright UPMC - CNRS\n\n"
    print "\n\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n"

    if os.path.isfile(args.p):
        print "\n\nParameter file :", args.p,"\n\n"
        if args.statsmod == True:
            print "Statistics will be be performed for", args.typesv, "\n"
            return vars(args), os.path.abspath(args.p), args.typesv+"--stats"
        else:
            if args.typesv == "all":
                print "All SV type will be detected\n"
                for val in svname:
                    if "stats" not in val:
                        print svname[val]
            else: 
                if "stats" not in args.typesv.lower():
                    print svname[args.typesv.upper()],"will be detected\n"
                else:
                    print "Statistical significance for", args.typesv, "candidates will \
be evaluated\n"
            return vars(args), os.path.abspath(args.p), args.typesv

    else:
        #print "\n\n\t\t Error :", args.p, "doesn't exist\n\n"
        return vars(args), 0, "0"

#-------------------------------------------------------------------------
def completeParams(params, paramfile):
    
    try:
        paramsForUp = U.get_run_info(paramfile)
    except:
        print "\tError: parameter file does not exist or is not well formatted\n\n"
        sys.exit()
    list_chr_real = paramsForUp["range"]
    
    
    updatable = ['stats', 'range', 'nsv', 'fdr', 'annotation', 'field_chr',\
    'field_type','field_start','field_end','field_sep', 'n', 'mapq', 'mapqx'  ]

    for categorie in paramsForUp:
        if categorie in updatable and categorie in params and params[categorie] != None:
            paramsForUp[categorie]=params[categorie]
    
    order = ['in', 'mapq', 'mapqx', 'out', 'vcf', 'stats', 'nsv', 'fdr', 'range', 'n', 'annotation', \
    'field_chr', 'field_type','field_start','field_end','field_sep'  ]
    with open(paramfile,"w") as pfile:
        for cate in order:
            if cate == "range" and "range" in paramsForUp.keys():
                if isinstance( paramsForUp[cate], list):
                    paramsForUp[cate] = ','.join(map(str, paramsForUp[cate]))
            try:
                pfile.write("%s%s%s\n" % (cate,"=", paramsForUp[cate]))  
            except KeyError:
                pfile.write("%s%s%s\n" % (cate,"=", params[cate])) 
 
    #Obselete
    check = {
    "in":"name of input file. example human for human_X X name of chromosome",\
    "out":"NA","stats":"NA","annotation":"NA", "field_chr":1,"field_type":3,\
    "field_start":4,"field_end":5,"sep":";",\
    "range":"all","n":6,"fdr":0.01,"nsv":100000}
    
    for val in check:
        if val not in params:
            if val == "in" :
                params["in"] = raw_input("You must give the name of the Library file") 
            else:
                params[val] = check[val]
        elif not params[val]:
                params[val] = check[val]
    if not params["out"] :
        params["out"] = params["in"]+"_Ulyssse"
    if not params["stats"]:
        params["stats"] = params["in"]+"_stats.txt"
    
    return list_chr_real
    
#-------------------------------------------------------------------------

params, paramfile, typesv = parser()


typesv= typesv.split('--')
if len(typesv)==1:
    typesv=typesv[0]
    onlyStatPerform = False
elif typesv[1] == 'stats':
    typesv=typesv[0]
    onlyStatPerform = True
else:
    typesv=typesv[0]
    onlyStatPerform = False
    

list_chr_real = completeParams(params, paramfile)







#print "typesv", typesv, onlyStatPerform
params = U.get_run_info(paramfile)   
dicovcf = {}

if typesv.lower() == "all":
    pval_seuil_del, pval_seuil_sins = deletion.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of deletion detection"
    gc.collect()
    os.system("date")    
    pval_seuil_dup = duplication.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of duplication detection"
    os.system("date")    
    gc.collect()

    pval_seuil_inv = inversions.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of inversions detection"
    gc.collect()

    pval_seuil_ins, pval_seuil_rt, pval_seuil_nrt = detect_inter.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of inter-chromosomal SV detection"
    gc.collect()
    os.system("date")    
    
    if params["vcf"] != False:
        
        dicovcf["DEL"] = (params["out"]+"_deletions_bySV.csv", pval_seuil_del)
        dicovcf["sINS"] = (params["out"]+"_small_insertions_bySV.csv", pval_seuil_sins)
        dicovcf["DUP"] = (params["out"]+"_duplications_bySV.csv", pval_seuil_dup)
        dicovcf["INV"] = (params["out"]+"_inversions_bySV.csv", pval_seuil_inv)
        dicovcf["RT"] = (params["out"]+"_reciprocal_translocations_bySV.csv", pval_seuil_rt)
        dicovcf["NRT"] = (params["out"]+"_non_reciprocal_translocations_bySV.csv", pval_seuil_nrt)
        dicovcf["INS"] = (params["out"]+"_insertions_bySV.csv", pval_seuil_ins)
        U.write_VCF(params["out"]+".vcf",dicovcf, params["vcf"], params)
    
elif typesv.upper() == "DUP":
    pval_seuil_dup = duplication.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of duplication detection"
    gc.collect()
    os.system("date")  
    if params["vcf"] != False:
        #print "params", params
        dicovcf["DUP"] = (params["out"]+"_duplications_bySV.csv", pval_seuil_dup)
        U.write_VCF(params["out"]+"."+typesv.upper()+".vcf",dicovcf, params["vcf"], params)
    

elif typesv.upper() == "DEL":
    
    
    pval_seuil_del, pval_seuil_sins = deletion.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of deletion detection"
    gc.collect()
    os.system("date")
    if params["vcf"] != False:
        dicovcf["DEL"] = (params["out"]+"_deletions_bySV.csv", pval_seuil_del)
        dicovcf["sINS"] = (params["out"]+"_small_insertions_bySV.csv", pval_seuil_sins)        
        U.write_VCF(params["out"]+"."+typesv.upper()+".vcf",dicovcf, params["vcf"], params)

elif typesv.upper() == "INV":
    pval_seuil_inv = inversions.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of inversions detection"
    gc.collect()
    os.system("date")
    if params["vcf"] != False:
        dicovcf["INV"] = (params["out"]+"_inversions_bySV.csv", pval_seuil_inv)
        U.write_VCF(params["out"]+"."+typesv.upper()+".vcf",dicovcf, params["vcf"], params)

elif typesv.upper() == "INTER":
    pval_seuil_ins, pval_seuil_rt, pval_seuil_nrt = detect_inter.launch(paramfile, onlyStatPerform, list_chr_real)
    print "End of inter-chromosomal SV detection"
    gc.collect()
    os.system("date")
    if params["vcf"] != False:
        dicovcf["RT"] = (params["out"]+"_reciprocal_translocations_bySV.csv", pval_seuil_rt)
        dicovcf["NRT"] = (params["out"]+"_non_reciprocal_translocations_bySV.csv", pval_seuil_nrt)
        dicovcf["INS"] = (params["out"]+"_insertions_bySV.csv", pval_seuil_ins)
        #print dicovcf
        U.write_VCF(params["out"]+"."+typesv.upper()+".vcf",dicovcf, params["vcf"], params)
    

