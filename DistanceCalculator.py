##Type -h to see usage instructions
import os
import utilities
import align_utilities
import pandas as pd
import mauve_utilities
import phylo_utilities
import distance_matrix_utilities
from collections import defaultdict, OrderedDict

import sys
from Bio import AlignIO
from Bio import Phylo

has_plt = False
try:
    import matplotlib.pyplot as plt
    has_plt = True
except ImportError:
    print("Cannot produce plots without the python module 'matplotlib'")
# import shlex

SCRIPT_VERSION = 1
SCRIPT_SUBVERSION = 0

script_base = os.path.basename(__file__)
_outputBase = '{}_v{}.{}'.format(os.path.splitext(script_base)[0],SCRIPT_VERSION,SCRIPT_SUBVERSION)

kSNP_all = 'SNPs_all_matrix.fasta'
kSNP_core = 'core_SNPs_matrix.fasta'
kSNP_core_tree = 'tree.core.tre'
kSNP_all_tree = 'tree.parsimony.tre'
snpExt = '.snps.tab'
analysisExt = '.snp.analysis.txt'
mapper_headers = ['SampleId','TreeLabel','GroupId']
import argparse

###df is a similarity matrix
def convertToFractionalDistance(df):
    my_max = df.max().max()
    if my_max > 1:
        raise ValueError("Similarity matrix cannot have value greater than 1. Max is {}".format(my_max))
    my_min = df.min().min()
    if my_min < 0:
        raise ValueError("Similarity matrix cannot have value less than 0. Min is {}".format(my_min))
    return 1-df

## df is a similarity matrix
def convertToPercentDifference(df):
    my_max = df.max().max()
    if my_max > 100:
        raise ValueError("Similarity matrix cannot have value greater than 100. Max is {}".format(my_max))
    my_min = df.min().min()
    if my_min < 0:
        raise ValueError("Similarity matrix cannot have value less than 0. Min is {}".format(my_min))
    return 100-df

## ROC for how well thresholds distinguing between in-group comparisons and between-group comparisons
    # output_basename: This will be the root name for several output files
    # dist_dict: dictionary of distance matrices. True positive rates for each matrix will be listed under column with the same name as the keys
    # groups is a dict of lists; each list is a group of isolates that should be considered in-group comparisons. All isolates must be in a group to be analyzed -- if a group should only be used for between-group comparisons it needs to be listed as 'nongroup'
def exportROCanalysis(output_basename,dist_dict,groups,exclude_groups=None,nongroup=None):
    if groups is None:
        print("No groups defined. Not performing ROC calculations.")
    else:
        try:
            for n,dist in dist_dict.items():
                if (not isinstance(dist,pd.DataFrame)) or (set(dist.columns.tolist()) != set(dist.index.tolist())):
                    raise ValueError("Distance matrix {} not properly formatted".format(n))
            list_file = utilities.appendToFilename(output_basename, '.pairlist')
            list_file = utilities.setExt(list_file,'.tab')
            roc_file = utilities.appendToFilename(output_basename, '.roc')
            roc_file = utilities.setExt(roc_file,'.xlsx')
            paired_list = distance_matrix_utilities.convertRecipricolTableToList(dist_dict) ##Makes a non-redundant list from a dictionary of reciprocal tables. Columns are name_1, name_2, and the keys from the dicts
            paired_list.to_csv(list_file,sep="\t")    
            classification_header = 'SameGroup'  
            grouped_list = distance_matrix_utilities.applyGroupingToReciprocalList(paired_list, groups, exclude_groups=exclude_groups, classification=classification_header, nongroup=nongroup) #appends 'classification' column, with True/False for whether samples are in same/different group (acconding to the groups dict).
            if len(grouped_list) > 0:
                grouped_list.to_csv(list_file,sep="\t")
                print("Grouped {}/{} pairs, from original set of {}".format(sum(grouped_list[classification_header]),len(grouped_list),len(paired_list)))
                data_cols = [x for x in dist_dict.keys()]
                ROC_results = distance_matrix_utilities.ROCanalysis(grouped_list,classification=classification_header,SimilarityMeasure=False,TestSet=data_cols)##perform ROC for distance thresholds for each item in TestSet, using classification as Positive/Negative result 
                try:
                    ROC_results['ROC'].fillna('=NA()').to_excel(roc_file)
                except IOError as e:
                    print("Failed to save ROC file:")
                    print(e)
                            
                if has_plt:
                    try:
                        image_file = utilities.appendToFilename(output_basename,'.roc')
                        image_file = utilities.setExt(image_file,'.png')                                 
                        fig = ROC_results['ROC'].set_index('fpr')[data_cols].plot(title='AUC={}'.format(ROC_results['AUC']))
                        fig = fig.get_figure()
                        fig.savefig(image_file)    
                    except Exception as e:
                        print('Failed to ROC scatterplot at '+image_file)
                        print(e) 
                    else:
                        try:
                            plt.close(fig)
                        except Exception as e:
                            print("Failed to close image...")
                            print(e)   
                print("AUC values")
                for k,v in ROC_results['AUC'].items():
                    print("\t{}:\t{}".format(k,v)) 
            else:
                print("The grouping list did not identify any isolate pairs with distance data. No ROC.")
        except Exception as e:
            print("Failed to perform ROC analysis")
            utilities.printExceptionDetails(e)

##Identify overlap between group lists. Groups is a dict of lists (Lab_ID) -- or is it a set (or either)
def groupsOverlap(groups):
    overlaps = defaultdict(OrderedDict)
    ##This is the clade (row in final report)
    for g, lids in groups.items():
        for g2, lids2 in groups.items():
            overlaps[g][g2] = sum([x in lids for x in lids2])        
    groupsFrame = pd.DataFrame([pd.Series(v) for v in overlaps.values()]).set_index('Group').astype(int)            
    return groupsFrame

##Identify how groups overlap on tree. Groups is a dict of lists (Lab_ID)
def groupsOnTree(groups,tree):
    clade_description = defaultdict(OrderedDict)
    ##This is the clade (row in final report)
    for g, lids in groups.items():
        nodes = []
        missing = []
#         print(g)
        for x in lids:
            node = tree.find_any(x)
            if node is None:
                missing.append(x)
            else:
                nodes.append(node)
        clade_description[g]['Group'] = g
        clade_description[g]['GroupIsolates'] = len(lids)
        clade_description[g]['GroupInTree'] = len(nodes)
#         print("\tFound {}/{}. Missing {}".format(len(nodes),len(lids),','.join(missing)))
        if len(nodes) > 0:
            anc = tree.common_ancestor(*nodes)
            term_names = [x.name for x in anc.get_terminals(order='level')]
#             print("\tClade has {} isolates".format(len(term_names)))
        else:
            term_names = []
        clade_description[g]['CladeInTree'] = len(term_names)
        for g2, lids2 in groups.items():
            clade_description[g][g2] = sum([x in term_names for x in lids2])        
    groupsFrame = pd.DataFrame([pd.Series(v) for v in clade_description.values()]).set_index('Group').astype(int)            
    return groupsFrame
    
def distance_matrix(snp_args,output_dir,groups=None):
    snp_matrix = snp_args.distance_file 
    if  os.path.isfile(snp_matrix):
        b = os.path.basename(snp_matrix)
        SNPs = pd.read_table(snp_matrix,index_col=0)
        if SNPs.columns.tolist() != SNPs.index.tolist():
            print("Columns do not match index. Assumed TAB delimited. Attempting CSV delimited.")
            SNPs = pd.read_csv(snp_matrix,index_col=0)
            if SNPs.columns.tolist() != SNPs.index.tolist():
                print("Still no luck. Exiting..")
                raise ValueError("Table rows and columns do not match")
            else:
                print("Success!")
        if snp_args.fractional_similarity:
            SNPs = convertToFractionalDistance(SNPs)
        elif snp_args.percent_similarity:
            SNPs = convertToPercentDifference(SNPs)
        SNPFile = os.path.join(output_dir,b+snpExt)
        distance_matrix_utilities.snpMatrixAnalysis(SNPs,analysis_file=os.path.join(output_dir,b+analysisExt),groups=groups,groupDistFile=utilities.appendToFilename(SNPFile, '.groups'))
        ##ROC
        if snp_args.do_roc:
            exportROCanalysis(SNPFile,{'tpr_Distances':SNPs},groups,nongroup=snp_args.ignore_group)
        ##
        if snp_args.distance_cluster:
            tree_file = utilities.setExt(SNPFile,'.tre')
            SNPs = distance_matrix_utilities.distanceMatrixCluster(SNPs,tree_file)
            cluster_file = utilities.appendToFilename(SNPFile, '.clustered')
            SNPs.to_csv(cluster_file,sep='\t')
        ###ROC
        
def alignmentBySet(aln_args,output_dir,groups=None):
    alignment_file = aln_args.alignment_file
    if not os.path.isfile(alignment_file):
        print("Error. Cannot find alignment file: {}".format(alignment_file))
    else:
        b = os.path.basename(alignment_file)
        a = AlignIO.read(alignment_file,'fasta')
        if aln_args.core_only:
            aln = align_utilities.stripDownToGATC(a)
            print("Count only core positions -- where all sequences have cannonical bases.")
        else:
            aln = a
        aln_len = aln.get_alignment_length()
        print('Opened alignment with {} genomes and {} positions.'.format(len(aln),aln_len))
        if aln_args.cluster_ignore:
            print("Using cluster_ignore: {}".format(aln_args.cluster_ignore))
            SNPs = align_utilities.alignment2snpBySetCluster(aln,ignore='-',cluster_ignore=aln_args.cluster_ignore,unambig_only=aln_args.unambig_only)
        else:
            SNPs = align_utilities.alignment2snpBySet(aln,ignore='-',verbose=True,unambig_only=aln_args.unambig_only)
        print('Maximum distance is {}'.format(SNPs.max().max()))
        SNPFile = os.path.join(output_dir,b+snpExt)
        SNPs.to_csv(SNPFile,sep='\t') 
        print("Wrote SNP distance to file: "+SNPFile)
        distance_matrix_utilities.snpMatrixAnalysis(SNPs,analysis_file=os.path.join(output_dir,b+analysisExt),groups=groups,groupDistFile=utilities.appendToFilename(SNPFile, '.groups'))
        SNPs_frac = SNPs/aln_len
        print('Maximum distance (fractional) is {}'.format(SNPs_frac.max().max()))     
        SNPs_frac.to_csv(utilities.appendToFilename(SNPFile, '.fractional'),sep='\t')
        distance_matrix_utilities.snpMatrixAnalysis(SNPs_frac,analysis_file=os.path.join(output_dir,b+'.fractional'+analysisExt),groups=groups,groupDistFile=utilities.appendToFilename(SNPFile, '.fractional.groups'))
        if aln_args.do_roc:
            exportROCanalysis(SNPFile,{'tpr_SNPs':SNPs},groups,nongroup=aln_args.ignore_group)  
        clusteredSNPs = distance_matrix_utilities.distanceMatrixCluster(SNPs,SNPFile)   
        cluster_file = utilities.appendToFilename(SNPFile, '.clustered')
        clusteredSNPs.to_csv(cluster_file,sep='\t') 

def XMFA(xmfa_args,output_dir,groups=None):
    xmfa_file = xmfa_args.xmfa_file
    if  os.path.isfile(xmfa_file):
        b = os.path.basename(xmfa_file)
        _, aln_list = mauve_utilities.parseXMFA(xmfa_file) ### parsnp XMFA does not parse filenames
        print("Read XMFA file. Found {} regions.".format(len(aln_list)))
        xmfaKey = mauve_utilities.getXMFAnameKey(xmfa_file)
        print("Identified the following {} samples:".format(len(xmfaKey)))
        aln_list = [a for a in aln_list if len(a) == len(xmfaKey)]
        print("Found {} regions containing all samples".format(len(aln_list)))        
        for k,v in xmfaKey.items():
            print('{}:\t{}'.format(k,v))
        total_snps = None
        if xmfa_args.core_only:
            print("Count only core positions -- where all sequences have cannonical bases.")
        aln_len = 0
        for a in aln_list:
            aln = align_utilities.stripDownToGATC(a) if xmfa_args.core_only else a##Core, no ambiguous
            aln_len += aln.get_alignment_length()
#             SNPs = align_utilities.alignment2snpBySet(aln,ignore='-')
            if xmfa_args.cluster_ignore:
                print("Using cluster_ignore: {}".format(xmfa_args.cluster_ignore))
                SNPs = align_utilities.alignment2snpBySetCluster(aln,ignore='-',cluster_ignore=xmfa_args.cluster_ignore,unambig_only=xmfa_args.unambig_only)
            else:
                SNPs = align_utilities.alignment2snpBySet(aln,ignore='-',verbose=True,unambig_only=xmfa_args.unambig_only)            
            SNPs_core = {x:x.split(':')[0] for x in SNPs.index} # cant use ParseMauveSeqName because the full string is not here        
            SNPs.rename(index=SNPs_core,columns=SNPs_core,inplace=True)
#             print(SNPs.mean().mean())
            if total_snps is None:
                total_snps = SNPs.copy()
            else:
                total_snps += SNPs ##This fails because of filename parsing
#             print(total_snps.mean().mean())
        print("Aligment is {} nucleotides".format(aln_len))
        total_snps.rename(index=xmfaKey,columns=xmfaKey,inplace=True)
        SNPFile = os.path.join(output_dir,b+snpExt)
        total_snps.to_csv(SNPFile,sep='\t')  
        print('Maximum distance is {}'.format(total_snps.max().max()))      
        print("Wrote SNP distance to file: "+SNPFile)      
        distance_matrix_utilities.snpMatrixAnalysis(total_snps,analysis_file=os.path.join(output_dir,b+analysisExt),groups=groups,groupDistFile=utilities.appendToFilename(SNPFile, '.groups'))
        SNPs_frac = total_snps/aln_len
        print('Maximum distance (fractional) is {}'.format(SNPs_frac.max().max()))
        SNPs_frac.to_csv(utilities.appendToFilename(SNPFile, '.fractional'),sep='\t')         
        distance_matrix_utilities.snpMatrixAnalysis(SNPs_frac,analysis_file=os.path.join(output_dir,b+'.fractional'+analysisExt),groups=groups,groupDistFile=utilities.appendToFilename(SNPFile, '.fractional.groups'))
        if xmfa_args.do_roc:
            exportROCanalysis(SNPFile,{'tpr_SNPS':SNPs},groups,nongroup=xmfa_args.ignore_group)    
        distance_matrix_utilities.distanceMatrixCluster(SNPs,SNPFile)  
            ##Look at    
            
##Outgroup not implemented -- would be name of a group in group dict.
def kSNP(kSNP_args,output_dir,groups=None,outgroup=None):
    kSNP_dir = kSNP_args.kSNP_directory
    if not os.path.isdir(kSNP_dir):
        raise ValueError("kSNP dir is not a directory")
    kSNP_allFile = os.path.join(kSNP_dir,kSNP_all)
    if  os.path.isfile(kSNP_allFile):
        all_aln = AlignIO.read(kSNP_allFile,'fasta')
        try:
            fasta_list = pd.read_table(os.path.join(kSNP_dir,'fasta_list'),names=['Filename','isolate'])
            print("Evaluating kSNP directory. Found 'fasta_list' file with {} entries".format(len(fasta_list)))
            invalid_files = sum(~fasta_list.Filename.apply(os.path.isfile))
            dup_names = sum(fasta_list.isolate.duplicated())
            if (dup_names > 0) or (invalid_files > 0):
                print("\t...containing {} invalid files and {} duplicated isolate names".format(invalid_files,dup_names))
        except IOError:
            pass            
        print("Analyzing 'all SNP' file. Ignoring gaps (which includes ambiguous characters due to kSNP algorithm)")
        print('Opened alignment of length {}.'.format(all_aln.get_alignment_length()))
        print('\twith {} sequences'.format(len(all_aln)))
        all_SNP = align_utilities.alignment2snpBySet(all_aln,ignore='-')
        print('Maximum distance is {}'.format(all_SNP.max().max()))
        all_SNPFile = os.path.join(output_dir,kSNP_all+snpExt)
        all_SNP.to_csv(all_SNPFile,sep='\t')        
        print("Wrote SNP distance (all) to file: "+all_SNPFile)
        ##groups
        groupDistFile = utilities.appendToFilename(all_SNPFile, '.groups')
        distance_matrix_utilities.snpMatrixAnalysis(all_SNP,analysis_file=os.path.join(output_dir,kSNP_all+analysisExt),groups=groups,groupDistFile=groupDistFile)
        ##clustering 
#         all_tree_file = utilities.setExt(all_SNPFile,'.tre')
        SNPs = distance_matrix_utilities.distanceMatrixCluster(all_SNP)
        cluster_file = utilities.appendToFilename(all_SNPFile, '.clustered')
        SNPs.to_csv(cluster_file,sep='\t') 
        dist_dict = {'tpr_kSNP_all':SNPs.copy()}
        ##Look at the core genome file
        kSNP_coreFile = os.path.join(kSNP_dir,kSNP_core)
        if os.path.isfile(kSNP_coreFile):
            print("Analyzing core SNP file")
            try:
                core_aln = AlignIO.read(kSNP_coreFile,'fasta')
            except (ValueError, IOError) as e:
                print("Failed to open core kSNP file")
                print(e)
            else:               
                print('Opened alignment of length {}.'.format(core_aln.get_alignment_length()))
                print('\twith {} sequences'.format(len(core_aln)))
                core_SNP = align_utilities.alignment2snpBySet(core_aln)
                print('Maximum distance is {}'.format(core_SNP.max().max()))
                core_SNPFile = os.path.join(output_dir,kSNP_core+snpExt)
                core_SNP.to_csv(core_SNPFile,sep='\t')
                print("Wrote SNP distance (core) to file: "+core_SNPFile)
                #grouping
                groupDistFile=utilities.appendToFilename(core_SNPFile, '.groups')
                distance_matrix_utilities.snpMatrixAnalysis(core_SNP,analysis_file=os.path.join(output_dir,kSNP_core+analysisExt),groups=groups,groupDistFile=groupDistFile)
                #clustering
                ##clustering 
                SNPs = distance_matrix_utilities.distanceMatrixCluster(core_SNP,core_SNPFile)
                cluster_file = utilities.appendToFilename(core_SNPFile, '.clustered')
                SNPs.to_csv(cluster_file,sep='\t')   
                dist_dict['tpr_kSNP_core'] =SNPs.copy()
        if kSNP_args.do_roc:
            exportROCanalysis(os.path.join(output_dir,'kSNP_ROC'),dist_dict,groups,nongroup=kSNP_args.ignore_group)        
        try: ##midpoint root trees
            intree = os.path.join(kSNP_dir,kSNP_all_tree)
            outtree = os.path.join(output_dir,kSNP_all_tree+'.midpoint.phylo.xml') ##Note: Biopython does not correctly write newick short branches
            tree = Phylo.read(intree,'newick')
            tree.root_at_midpoint()
            Phylo.write(tree,outtree,'phyloxml')
        except Exception as e:
            print("Failed to reroot all SNP tree")
            utilities.printExceptionDetails(e)
        else:
            if isinstance(groups,dict):
                try:
                    allTreeGroups = groupsOnTree(groups,tree)
                except:
                    pass
                else:
                    allTreeGroupsFile=utilities.appendToFilename(outtree.rstrip('.phylo.xml')+'.tab', '.groups')
                    allTreeGroups.to_csv(allTreeGroupsFile,sep='\t')
        try: ##midpoint root trees
            intree = os.path.join(kSNP_dir,kSNP_core_tree)
            outtree = os.path.join(output_dir,kSNP_core_tree+'.midpoint.phylo.xml') ##Note: Biopython does not correctly write newick short branches
            tree = Phylo.read(intree,'newick')
            tree.root_at_midpoint()
            Phylo.write(tree,outtree,'phyloxml')            
        except Exception as e:
            print("Failed to reroot core SNP tree")
            utilities.printExceptionDetails(e)         
        else:
            if isinstance(groups,dict):
                try:
                    coreTreeGroups = groupsOnTree(groups,tree)
                except:
                    pass
                else:
                    coreTreeGroupsFile=utilities.appendToFilename(outtree.rstrip('.phylo.xml')+'.tab', '.groups')            
                    coreTreeGroups.to_csv(coreTreeGroupsFile,sep='\t')   
            
    else:
        raise ValueError("kSNP directory does not have result file: {}".format(kSNP_all))     
    return None

tree_formats = ['newick','nexus','nexml','phyloxml','cdao']

def tree_dist(tree_args,output_dir,groups=None):
    tree_file = tree_args.tree_file
    if  os.path.isfile(tree_file):
        b = os.path.basename(tree_file)
        if tree_args.tree_format not in tree_formats:
            raise ValueError("illegal tree format: {}. Options are {}".format(tree_args.tree_format,",".join(tree_formats)))
        tree = Phylo.read(tree_file,tree_args.tree_format)
        dist = phylo_utilities.phylogeneticDistanceTable(tree)
        print('Maximum distance is {}'.format(dist.max().max()))
        distFile = os.path.join(output_dir,b+snpExt)
        dist.to_csv(distFile,sep='\t')        
        print("Wrote SNP distance to file: "+distFile)
        distance_matrix_utilities.snpMatrixAnalysis(dist,analysis_file=os.path.join(output_dir,b+analysisExt),groups=groups,groupDistFile=utilities.appendToFilename(distFile, '.groups'))
        if tree_args.do_roc:
            exportROCanalysis(distFile,{'tpr_TreeDist':dist},groups,nongroup=tree_args.ignore_group)        
                          
def main():
    parser = argparse.ArgumentParser(description='A program to calculate distances between genomes, starting from either alignments or assemblies.',
                                     epilog = 'Note: alignment distances are based on any nucleotide difference, including ambiguious characters unless otherwise specified (e.g.not case sensitive, gaps not counted (but other bases in the column are))')
    ### general info
    parser.add_argument('--version','-V',action='version',version='%(prog)s {}.{}'.format(SCRIPT_VERSION,SCRIPT_SUBVERSION))
    parser.add_argument('--GT_mapper','-gm',help='Mapper file from Changayil''s Genotyper. For specifying groups')
    parser.add_argument('--SampleId',action='store_true',help='Identify isolates in GT_mapper by SampleId rather than TreeLabel')
    parser.add_argument('--ignore_group','-ig',help='For ROC analysis, do not attempt to classify members of this group. Group name must be in group mapper file.')
    parser.add_argument('--group_list','-gl',help='Specify groups with a text file. Each row has format "group_name:member1,member2,..."')
    parser.add_argument('--do_roc',help='Perform ROC analysis to test if isolates within group are more similar than isolates between groups')
    parser.add_argument('-o','--output',help='Directory to write summary statistics to. If not specified, will make directory starting with: '+_outputBase)
#     parser.add_argument('--debug',action='store_true',help="Preserve intermediate files and do not update reference files")
    subparsers = parser.add_subparsers(description="Select one of the following commands",dest='subcommand')
    subparsers.required = True
    
    ##Extract  
    kSNP_parser = subparsers.add_parser('kSNP',description="Produce two SNP distances from kSNP results (all SNPs and core SNPS)")
    kSNP_parser.set_defaults(func=kSNP)
    kSNP_parser.add_argument('kSNP_directory',help='Location of kSNP results files, containing {} and optionally {}'.format(kSNP_all,kSNP_core))
     
    alignmentSet_parser = subparsers.add_parser('aln',description="Produce SNP distances from an alignment file. Currently limited to FASTA format")
    alignmentSet_parser.set_defaults(func=alignmentBySet)
    alignmentSet_parser.add_argument('alignment_file',help='Location of FASTA alignment file.')
    alignmentSet_parser.add_argument('--cluster_ignore','-c',help='Size of SNP clusters to ignore',type=int)  
    alignmentSet_parser.add_argument('--unambig_only',action='store_true',default=False,help='Only count differences if both isolates have unambiguous nucleotides (GATC/gatc)')  
    alignmentSet_parser.add_argument('--core_only',action='store_true',default=False,help='Limit analysis to positions where all genomes have unambigous nucleotide')
    
    xmfa_parser = subparsers.add_parser('xmfa',description="Produce SNP distances from an XMFA file (tested with Mauve)")
    xmfa_parser.set_defaults(func=XMFA)
    xmfa_parser.add_argument('xmfa_file',help='Location of XMFA file')
    xmfa_parser.add_argument('--cluster_ignore','-c',help='Size of SNP clusters to ignore',type=int,default=0)
    xmfa_parser.add_argument('--unambig_only',action='store_true',default=False)
    xmfa_parser.add_argument('--core_only',action='store_true',default=False,help='Limit analysis to positions where all genomes have unambigous nucleotide')    

    tree_parser = subparsers.add_parser('tree',description="Produce SNP distances from an phylogenetic tree. Currently limited to Newick format")
    tree_parser.set_defaults(func=tree_dist)
    tree_parser.add_argument('tree_file',help='Location of  tree file')
    tree_parser.add_argument('--tree_format',help='Format of tree file. Biopython options.',default='newick')

    SNP_parser = subparsers.add_parser('dist',description="Apply grouping analysis to a distance matrix (tab delimited)")
    SNP_parser.set_defaults(func=distance_matrix)
    SNP_parser.add_argument('--distance_cluster', action='store_true',help='Use a NJ tree to reorder the SNP distance matrix. Save to output directory with "clustered" extension')
    SNP_parser.add_argument('--fractional_similarity',action='store_true',help="Matrix is fractional similarity. Convert to fractional difference before analyzing (D = 1-S")
    SNP_parser.add_argument('--percent_similarity',action='store_true',help="Matrix is percent similarity. Convert to percent difference before analyzing (D = 100-S")
    SNP_parser.add_argument('distance_file',help='Location of distance file (tab delimited)')
    
    mash_parser = subparsers.add_parser('mash',description="Calculacte MASH distances for a group of assemblies")
    mash_parser.set_defaults(func=mash_matrix)
    mash_parser.add_argument('assemblies',help='Location of either a directory with assemblies, or a tab-delimited file listing assemblies (with header "Filename"')

    
    args = parser.parse_args()

    try:        
        output = args.output
        output_dir = None
        if isinstance(output,str):
            try:
                utilities.safeMakeDir(output)
                output_dir = output
                print("Output directory is: "+os.path.realpath(output_dir))
            except OSError:
                print ("error making output dir")
        if output_dir is None:
            output_dir =  utilities.safeMakeOutputFolder(_outputBase + '_' + args.subcommand)   
            print("Created output directory: "+output_dir)
               
        ##Run analysis      
        sys.stdout = utilities.Logger(os.path.join(output_dir,'DistanceCalculator.log'))
        print('Version {}.{}'.format(SCRIPT_VERSION,SCRIPT_SUBVERSION))
        print("Options are:")
        for arg in vars(args):
            print (arg, getattr(args, arg))        
        print()
                ##Read group information if it is provided
        groups = {}
        GT_mapper = args.GT_mapper
        if GT_mapper is not None:
            if not os.path.isfile(GT_mapper):
                raise ValueError("The listed mapper file does not exist")
            print("Reading mapper file...")
            mapFrame = pd.read_table(GT_mapper)
            missingHeaders = set(mapper_headers).difference(set(mapFrame.columns.tolist()))
            if len(missingHeaders) > 0:
                raise ValueError("The mapper file is missing the following headers: {}".format(', '.join(missingHeaders)))
            for n,g in mapFrame.groupby('GroupId'):
                groups[n] = g['SampleId'].tolist() if args.SampleId else g['TreeLabel'].tolist()
        groupFile = args.group_list
        if groupFile is not None:
            if not os.path.isfile(groupFile):
                raise ValueError("The listed group file does not exist")
            print("Reading group file...")
            with open(groupFile,'rt') as fin:
                for line in fin:
                    l = line.strip()
#                     print(l)
                    if l != '':
                        try:
                            gname, gcomma = l.split(':',maxsplit=1)
                        except ValueError:
                            print("Failed to split on colon. Trying comma...")
                            gname, gcomma = l.split(',',maxsplit=1)
                        glist = [x for x in gcomma.split(',') if x != '']
                        groups[gname] = glist
        ##Report group information
        if len(groups) > 0:
            print("Parsed the following {} groups from input file (number of isolates in parentheses)".format(len(groups)))
            for n,g in groups.items():
                print('{} ({})'.format(n,len(g)))
        if args.ignore_group:
            if args.ignore_group in groups.keys():
                print("Ignoring group: {}".format(args.ignore_group))
            else:
                print("Failed to identify {} among the groups. Not ignoring any group".format(args.ignore_group))
        if len(groups) == 0:
            groups = None ##All of the functions expect "None" for absent group, but this routine has an empty dict   
        try:
            args.func(args,output_dir,groups)
        except (ValueError, IOError) as e:
            print('Exception in main function')
            raise 
    except (ValueError, IOError) as e:
        utilities.printExceptionDetails(e)
        parser.print_usage()   
              
    
    
if __name__ == "__main__":
    if not utilities.has_preferred_python():
        raise Exception("Upgrade your python version")
    main()
