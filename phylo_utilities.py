##Calculate stats on a tree


import os
import sys
sys.path.append('/home/adam/Software/scripts/Utility/')
import align_utilities
import align_table_utilities
import pandas as pd
from io import StringIO
# import numpy as np
# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from Bio import Phylo
# from collections import defaultdict
from Bio.Alphabet import IUPAC
unambig = IUPAC.IUPACUnambiguousDNA.letters
unambig_set = set(unambig)
ambig_all_set = set(IUPAC.IUPACAmbiguousDNA.letters)
ambig_only_set = ambig_all_set.difference(unambig_set)

SCRIPT_VERSION = 1.1 #add 'fq' extension

script_base = os.path.basename(__file__)
_outputBase = '{}_v{}'.format(os.path.splitext(script_base)[0],SCRIPT_VERSION)

# has_ete3 = False
# try:
#     from ete3 import Phyloxml    
#     has_ete3 = True
# except ImportError:
#     print("Cannot provide advanced tree manipulation without ETE3")
#     print("\tNote for Adam -- source activate py344 ")

###This dscribes the diversity within each of the clades
##Can pass SNP matrix to save calculation time
def ParseMatrixWithTree(aln,tree,snpFrame=None):
    if snpFrame is None:
        snpFrame = align_utilities.alignment2snp(aln)
    nodes =  tree.get_nonterminals(order='level')
    n = 0
#     cladeInfoFile = os.path.join(workingDir,'clade_info_test.tab')
#     labeledTreeFile = os.path.join(workingDir,'clade_labeled.tree_test.phylo.xml')
    result_list = []
    for node in nodes:
        subclade = [x.name for x in node.get_terminals(order='level')]
        n+=1
        snp_stats = align_utilities.snpSubmatrixStats(snpFrame,subclade)
        site_stats = align_utilities.alignmentStats(aln,subclade)
        result = snp_stats.append(pd.Series(site_stats))
        node.name = str(n)
        result['Isolates'] = len(subclade)
        result['CladeID'] = n
        result['MemberList'] = ', '.join(subclade)
        result_list.append(result)
    cladeFrame = pd.DataFrame(result_list).set_index('CladeID')
    return cladeFrame

##seqFrame is a sequence alignment in the form of a pd.DataFrame, with position information. (Lyve-set out.filteredbcftoolsquery.tsv)
##ReconcileDict converts tree names to alignment names if necessary (subgroupingKey)
##Get seqFrame with openFilteredBCFToolsQuery, possibly stripping ambigs
##reconcileDict maps node.name to the column headers in the seqFrame
##Table should be indexed on position and have no columns other than genome names (no chr, no REF)
## THis will ignore ambiguous nucleotides, so a substitution could be marked to two different branches if the intervening branch is ambiguous.

##This provides detailed, site-specific comparisons of the clades.
def ParseTableWithTree(seqFrame,tree,reconcileDict=None):
    def reconcileNames(node_list,valid_names):
        result = [x.name for x in node_list]
        if reconcileDict is not None:
            for i in range(len(result)):
                if result[i] in reconcileDict:
                    old = result[i] 
                    result[i] = reconcileDict[old]
        result = [x for x in result if x in valid_names]
        return result
    
    full_tree = reconcileNames(tree.get_terminals(order='level'),seqFrame.columns.tolist())
    term_nodes = tree.get_terminals(order='level')
    nodes =  tree.get_nonterminals(order='level') + term_nodes
    n = 0
    summaries1 = []
    details1 = []
    ###### first identify ambiguous sequences so they can be removed from 
    
    ###### iterate over tree
    for node in nodes:
        if node not in term_nodes:
            n += 1
            node.name = n
        subclade = reconcileNames(node.get_terminals(order='level'),seqFrame.columns.tolist())
        outgroup = [x for x in full_tree if x not in subclade] 
        node_summary, clade_details = align_table_utilities.compareCladeSNPS(seqFrame,subclade,node.name,outgroup)
        details1.append(clade_details)   
        summaries1.append(node_summary)     
    return details1, summaries1    

# def relabelPhylogeny(tree,keyFile):
#     for x in tree.get_terminals():
#         old_name = x.name
#     #     print("rename {} with {}".format(old_name,renamed[old_name]))
#         x.name = renamed[old_name]
    
def phylogeneticDistanceTable(tree):
    nodes = tree.get_terminals() #should be clustered
    labels = [n.name for n in nodes]
    d_table = [[0 for i in range(len(nodes))] for j in range(len(nodes))]
    for i in range(len(nodes)):
        n1 = nodes[i]
        for j in range(i):
            n2 = nodes[j]
            d = tree.distance(n1,n2)
            d_table[i][j] += d
            d_table[j][i] += d    
    return pd.DataFrame(d_table,columns=labels,index=labels)

##Identifies proper output file (for parsing bootstrap), and then dose some default reorientation/collapsing
##Directory is the directory where RAxML put the result files
##Project ID is the parameter that was sent to RAxML
##renameIsolates is a dictionary with lookup information (key is in tree, value is the new label you want on the tree)
def openRAxMLtree(directory,projectID,renameIsolates=None):
    treeFile = os.path.join(directory,'RAxML_bipartitions.'+projectID)
    try:
        MLtree = Phylo.read(treeFile,'newick')
    except:
        print('Unable to open RAxML tree: '+treeFile)
        raise
    else:    
        ##RAxML removes duplicate sequences and adds them back after the fact with no boostrap information, and branch length of 0
        print("Opened RAxML bipartitions tree with {} isolates".format(len(MLtree.get_terminals(order='level'))))
        c = 0
        for n in MLtree.get_nonterminals(order='level'):
            if (n.branch_length == 0.0):
                if (n != MLtree.root):
                    MLtree.collapse(n)
                    c += 1
        print("\tCollapsed {} branches with length 0.".format(c))
        c = 0
        for n in MLtree.get_nonterminals(order='level'):
            if (n.confidence == 0):
                if (n != MLtree.root):
                    MLtree.collapse(n)
                    c += 1
        print("\tCollapsed {} branches with confidence = 0.".format(c))        
        if isinstance(renameIsolates,dict):
            renamed = []
            failed = []
            for iso in MLtree.get_terminals(order='level'):
                name = iso.name
                if name in renameIsolates:
                    iso.name = renameIsolates[name]
                    renamed.append(iso.name)
                else:
                    failed.append(iso.name)
            print("Renamed {} isolates; omitted {}".format(len(renamed),len(failed)))
        elif renameIsolates is not None:
            print("Can only rename isolates with a dict")
        print("Remember that BioPython does not reassign confidence values upon rerooting. Confidence is 0-100")
    return MLtree    


# def biopythonTree2ete3(biopython_tree):
#     if not has_ete3:
#         print("Did not detect ETE3 on this installation. Note for Adam -- source activate py35 ")
#     project = Phyloxml()
#     project.build_from_file(StringIO(biopython_tree.format('phyloxml')))
#     project.export()
#     return project.get_phylogeny()[0]