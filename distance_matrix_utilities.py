##Testing
import pandas as pd
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
import utilities
from sklearn import metrics
import sys
from Bio import Phylo

####### Formatting distance matrix  #####################

### Order Based on NJ tree
##Default is to use the matirx itself, alternative is alphabetical
## df is a distance matrix (dataframe); output is a optional output filename
def distanceMatrixCluster(df,tree_output=None):
#     print("Clustering for {} items".format(len(df)))
    if len(df) > 1:
        names = [str(x) for x in df.index.tolist()]
        name_set = set(names)
        cols = df.columns.tolist()
        col_set = set(cols)
        assert name_set <= col_set, "Some of the listed records do not have corresponding columns: {}".format(name_set.difference(col_set))
        extra_cols = [c for c in cols if c not in names] ##Maintain order; these will go first
        matrix = df[names].copy()
        assert matrix.size == np.sum(matrix.count()), 'Failure. There appears to be empty cells in SNP matrix'
        dm = _DistanceMatrix(names)
        for i in range(len(dm.matrix)):
            for j in range(i):
                mean = (matrix.iloc[i,j] + matrix.iloc[j,i]) / 2
                dm[i,j] =  mean
        tree = DistanceTreeConstructor().nj(dm)
        tree.root_at_midpoint()
        try:
            if isinstance(tree_output,str):
                phylo_out = utilities.setExt(tree_output,'.phyloxml')
                try:
                    Phylo.write(tree,phylo_out,'phyloxml')
                except:
                    pass
                newick_out = utilities.setExt(tree_output,'.nwk')
                try:
                    Phylo.write(tree,newick_out,'newick')
                except:
                    pass
        except Exception as e:
            utilities.printExceptionDetails(e)
        ordered = [x.name for x in tree.get_terminals()]
        df = df.reindex(ordered)
        matrix = matrix.reindex(ordered,ordered)
        result = df[extra_cols].merge(matrix,left_index=True,right_index=True)
    else:
        result = df.copy()
    return result

### Order by sequence name 
def alphabetizeMatirx(df):
    names = df.index.tolist()
    name_set = set(names)
    cols = df.columns.tolist()
    col_set = set(cols)
    assert name_set <= col_set, "Some of the listed records do not have corresponding columns: {}".format(name_set.difference(col_set))
    extra_cols = [c for c in cols if c not in names] ##Maintain order; these will go first
    matrix = df[names].copy()
    assert matrix.size == np.sum(matrix.count()), 'Failure. There appears to be empty cells in SNP matrix'
    ordered = sorted(names)
    df = df.reindex(ordered)
    matrix = matrix.reindex(ordered,ordered)
    return df[extra_cols].merge(matrix,left_index=True,right_index=True)



##Groups is a dictionary of group names with lists of isolate names
def conformGroupsToNames(names,groups,remove_extras=False):
    new_groups = {}
    for key,values in groups.items(): ##Values is a list of names -- assume abbreviated versions of "names"
        new_values = []
        for v in values:
            new_v = None
            for n in names:
                if n.startswith(v):
                    if new_v is not None:
                        print("Warning: group value {} matches multiple names: {}, {}...".format(v,new_v,n))
                    new_v = n
            if new_v is not None:
                new_values.append(new_v)
            else:
                if not remove_extras:
                    new_values.append(v)                
        new_groups[key] = new_values
    return new_groups

############## Complex Statistical descriptions of alignments ######################

def snpSubmatrixStats(snpFrame,group_list=None,out_list=None,verbose=False): ##untested
    if group_list is not None:
        old_len = len(group_list)
        group_list = [x for x in group_list if x in snpFrame.columns.tolist()]
        new_len = len(group_list)
        if (new_len < old_len) and verbose:
            print("Only identified {}/{} requested isolates in the distance matrix".format(old_len,new_len))
    minOut = None
    if (group_list is None) or (sorted(group_list) == sorted(snpFrame.columns.tolist())):
        subMatrix = snpFrame
        minOut = None        
    else:
        assert isinstance(group_list,list)
        ###generate the cross matrix to find closest outgroup
        if out_list is None:
            if verbose:
                print("All other sequences in alignment are treated as outgroup for SNP counting")
            out_list = [x for x in snpFrame.columns if x not in group_list]
        if len(out_list) > 0:
            crossMatrix = snpFrame.reindex(index=group_list,columns=out_list)  
            minOut = crossMatrix.values.min()
        ##Could check that group_list is a subset of index and columns
        subMatrix = snpFrame.reindex(index=group_list,columns=group_list) 
    ##Get distribution stats by converting to an array and dropping self-comparisons
    mask = np.ones(subMatrix.shape, dtype=bool)
    np.fill_diagonal(mask, 0)
    maskedFrame = subMatrix.where(mask) 
    mySeries = pd.Series(maskedFrame.values.reshape(-1)).dropna()
    result = mySeries.describe(percentiles= [.05,.10,.25, .5, .75,.90,.95])
    result['MinOut'] = minOut
    return result

##For use with SNP matrices
def findClosestInMatrix(focus,distance_matrix,reference_list=None):
    focus_row = distance_matrix.loc[focus]
    del focus_row[focus]
    if isinstance(reference_list,list):
        focus_row = focus_row[focus_row.index.isin(reference_list)]
    row_min = focus_row.min()
    min_indicies = focus_row[focus_row == row_min].index.tolist()
    return {'Min_value': row_min, 'Closest': min_indicies}

##For use with alignment2InsDelSnp
def reportDistancesInDict(focus,comparison_list,distance_matrices):
    result = {}
    for key,matrix in distance_matrices.items():
        this = {}
        for item in comparison_list:
            this[item] = matrix.loc[focus,item]
        result[key] = this
    return result

##Analysis file is output for a text report; otherwise stdout
## Groups is a dict of lists of isolates
## Group dist file is output for a matrix of minimum distance between groups        
def snpMatrixAnalysis(snpFrame,analysis_file=None,groups=None,groupDistFile=None):


    def snpDiversity(snpMatrix):
        values = snpSubmatrixStats(snpMatrix)
        lprint("Maximum distance: {}".format(values["max"]))
        snp_total = snpMatrix.values.sum()
        lprint("Mean distance: {}".format(values['mean']))  
        return values
    
    def lprint(text='',output=None):
        print(text)
        try:
            print(text,file=output)
        except:
            pass    
            
    maxD = 'Max_Divergence'
    meanD = 'Mean_Divergence'
    if groups is None:
        groups = {}
    else:
        groups = {str(x):y for x,y in groups.items()} ##untested
    assert isinstance(groups,dict)
    groups = conformGroupsToNames(snpFrame.columns,groups)
    
    assert isinstance(snpFrame, pd.DataFrame), 'Need a matrix'
    assert snpFrame.index.tolist() == snpFrame.columns.tolist(), 'Matrix columns do not match rows'
    ##Consider forcing the matrix into a symetrical square below... maybe make an option to force into square
    output = open(analysis_file,'wt') if isinstance(analysis_file,str) else sys.stdout

        
    ### Invert the groups so that I can lookup membership:
    isolate2group = {}
    for name,group in groups.items():
        for g in group:
            isolate2group[g] = name
        
    #Find if there are identical sequences
    mask = np.ones(snpFrame.shape, dtype=bool)
    np.fill_diagonal(mask, 0)
    identicals = (sum(snpFrame.values[mask] == 0))//2 ##Counting each comparison twice; Cannot find anything more elegant
    comparisons = (snpFrame.size - len(snpFrame)) // 2
    lprint("There are {}/{} comparisons with zero SNPs in this dataset".format(identicals,comparisons))
    if identicals > 0:
        ##list them
        maskedFrame = snpFrame.where(mask)   
        bs = (maskedFrame == 0).any() ##boolean_series
        # snpFrame.iloc[0]
        id_set = set(bs[bs].keys().tolist())
        while len(id_set) > 0:
            key = id_set.pop()
            this_list = [key]
            bs2 = (maskedFrame.loc[key] == 0)
            this_list += bs2[bs2].keys().tolist()
            lprint("Identical sequences: "+', '.join(this_list))
            if len(groups) > 0:
                these_groups = set()
                for i in this_list:
                    my_group = isolate2group[i] if i in isolate2group else "NotAssigned"
                    these_groups.add(my_group)
                if len(these_groups) > 1:
                    lprint("\tWarning: multiple groups have identical sequences in them: {}".format(these_groups))
                else:
                    lprint("\tAll are members of the following group: {}".format(these_groups.pop()))
            id_set.difference_update(this_list)
    
    snpDiversity(snpFrame)
    
    if len(groups) > 0:
        
        ##Check membership
        full_set = set(snpFrame.index.tolist())
        ungrouped = full_set.copy()
        ##Identify minimum distance between each
        columns = ['Member_Count',maxD,meanD]+[x for x in groups.keys()]
        groupMin = pd.DataFrame(index=groups.keys(),columns=columns)
        groupMin.index.name = 'GroupID'
#         groupMin_temp = 
        lprint("Group count: {}".format(len(groups)))
        emptyGroups = set()
        for name,group in groups.items():
            ##Check membership -- both that all group names are valid and that all sequences are included in a group (no test for multiple grouping)
            group_set = set(group)
            extras = group_set.difference(full_set)
            if len(extras) > 0:
                lprint("Warning: The sequence groups included sequences that were not in the alignment. They are...\n\t"+
                     ", ".join(extras))                
            ungrouped.difference_update(group_set)
            group_set.difference_update(extras) ##Can't include the isolates that are not present
            groups[name] = list(group_set)
            if len(group_set) == 0:
                lprint("Warning: group {} no longer has any members.".format(name))
                emptyGroups.add(name)
        for empty in emptyGroups:
            if empty in groups:
                del groups[empty]
        for name,group in groups.items():
            if len(group) > 0:
                ### Group data
                lprint()
                lprint('group "{}"; length: {}'.format(name,len(group)))
                lprint('members: '+ ", ".join(group))
                groupMin.loc[name,'Member_Count'] = len(group)
                if len(group) > 1:
                    diversityWithin = snpDiversity(snpFrame.reindex(index=group,columns=group))
                    groupMin.loc[name,maxD] = diversityWithin["max"]
                    groupMin.loc[name,meanD] = diversityWithin["mean"]
                for name2,group2 in groups.items():
                    assert len(group2) > 0, "Empty group in dictionary"
                    snp_min = snpFrame.reindex(index=group,columns=group2).values.min()
                    groupMin.loc[name,name2] = snp_min
        lprint()
        lprint("The following sequences were not grouped: "+", ".join(ungrouped))
        lprint()
        ##This maybe should go into  distanceMatrixCluster
        samples = groupMin.index.tolist()
        groupMin.dropna(axis=0,how="all",subset=samples,inplace=True) ##THis may be bug prone...but need to mask
        groupMin.dropna(axis=1,how="all",inplace=True)
        if len(groupMin) > 2:
            try:
                groupMin = distanceMatrixCluster(groupMin)
            except Exception as e:
                print("Failed to reorient matrix of group minimum distances")
                utilities.printExceptionDetails(e)
        groupMinAlpha = alphabetizeMatirx(groupMin)
        if groupDistFile is not None:
            try:
                groupMin.to_csv(groupDistFile,sep='\t')
            except Exception as e:
                print('Failure to save group distance matrix to file: '+groupDistFile)
                lprint('Minimum distance between groups:')
                lprint(groupMin.to_csv(sep='\t'))
                utilities.printExceptionDetails(e)
            groupDistFileAlpha = groupDistFile + ".alpha.csv"
            groupMinAlpha.to_csv(groupDistFileAlpha)
        else:
            lprint('Minimum distance between groups:')
            lprint(groupMin.to_csv(sep='\t'))

##Makes a non-redundant list from a dictionary of reciprocal tables. 
def convertRecipricolTableToList(df_dict):
    result = pd.DataFrame()
    for name,df in df_dict.items():
        for i in range(len(df.index)):
            idx = df.index[i]
            for j in range(i):
                jdx = df.index[j]
                labels = ','.join(sorted([idx,jdx]))
                result.loc[labels,name] = df.loc[idx,jdx]
                result.loc[labels,'name_1'] = idx
                result.loc[labels,'name_2'] = jdx
    return result

##Uses the output of convertRecipricolTableToList
##Group_dict is a dict of lists, should match items in list frame. Samples that are not in the group_dict will be dropped from the result table. Groups are assumed to be mutually exclusive 
##Result with be a copy of list frame with a new colum using 'classification' as header (e.g. SameGroup) marked True/False according to whether they were in same/different group.
    # Default behavior is to mark True if two samples are in the same group list and False if they are in different group lits (dropped if one is not in a group list)
    # 'nongroup' specifies that all samples in the named group list should be treated as singletons (i.e. not matching anything)
    #  'exclude_groups' is a dict of lists with the same keys as group_dict. A pair is only considered to NOT Match if one is in group_dict and the other is in the matching exclude_groups list

def applyGroupingToReciprocalList(list_frame,group_dict,exclude_groups=None,classification = 'SameGroup',nongroup=None):
    #Validate group dict data structure
    for name,group in group_dict.items():
        if not isinstance(group,list):
            raise ValueError("Groups must be lists")
    ##Make a data frame with GroupID for each sample. Groups are assumed to be mutually exclusive.
    groupFrame = pd.DataFrame(columns=['SampleId','GroupId']).set_index('SampleId')
    for name,group in group_dict.items():
        for g in group:
            groupFrame.loc[g,'GroupId']=name #If ID is listed twice, first instance will be over-written
    print("Found {} grouped samples:".format(len(groupFrame)))
    ##Validate that nongroup is one of the groups
    if nongroup is not None:
        if nongroup not in group_dict:
            raise ValueError("The 'nongroup' identifier must be a key in the 'group_dict'")
        if exclude_groups is not None:
            raise ValueError("The use of 'nongroup' is incompatible with 'exclude_group'") 
    ## Confirm that exclude_groups matches groups
    if exclude_groups is not None:
        if isinstance(exclude_groups,dict):
            gkeys = set(group_dict.keys())
            ekeys = set(exclude_groups.keys())
            if not ekeys == gkeys:
                raise ValueError("Excluded group keys must match the group dict keys")
        else:
            raise ValueError("Excluded group must be dict")
    ## Validate list frame
    if ('name_1' in list_frame.columns) and ('name_2' in list_frame.columns):
        count = sum(list_frame['name_1'].isin(groupFrame.index)) +  sum(list_frame['name_2'].isin(groupFrame.index))
        print("Found {} pairs where both members are in groups.".format(count))
        if count == 0:
            print("No matches. Look at them...")
            print("...from paired distance frame:")
            print(list_frame.name_1.tolist())
            print("...from grouping file:")
            print(groupFrame.index.tolist())
    else:
        print("Error. Could not find name columns in list_frame... continuing anyway without validation")
                    


    ##For earch pair, check if isolates are in same group
    keepers = [] ##this list will keep each row in which both samples are in groupFrame
    for i in list_frame.index:
        (n1,n2) = i.split(',')
        ref1 = groupFrame.loc[n1,'GroupId'] if (n1 in groupFrame.index) else None
        ref2 = groupFrame.loc[n2,'GroupId'] if (n2 in groupFrame.index) else None        
        if exclude_groups is None: ##Groups are tested on exclusion of other groups
            if (ref1 is not None) and (ref2 is not None):
                same_group = (ref1 == ref2) and (ref1 != nongroup) #nongroup items cannot match any group
                r = list_frame.loc[i].copy()
                r[classification] = same_group
                keepers.append(r)
        else: ##Groups are matched against specific test cases
            if (ref1 is not None) or (ref2 is not None):
                same_group = (ref1 == ref2)
                excluded_pair = (n1 in exclude_groups[ref2]) or (n2 in exclude_groups[ref1])
                if same_group or excluded_pair:
                    r = list_frame.loc[i].copy()
                    r[classification] = same_group
                    keepers.append(r)    
    print("Keeping {} pairs where at least one is in a group, and the other is included in analysis.".format(len(keepers)))            
    return pd.DataFrame(keepers)

def applyGroupingThresholdToReciprocalList(list_frame,group_dict,exclude_groups=None,classification = 'SameGroup',nongroup=None):
    groupFrame = pd.DataFrame(['SampleId','GroupId'])
#     sample_list = []
    for name,group in group_dict.items():
        if not isinstance(group,list):
            raise ValueError("Groups must be lists")
#         sample_list += group
    ##Can sample ID be used twice? I guess.
    for name,group in group_dict.items():
        for g in group:
            groupFrame.loc[g,'GroupId']=name
    if nongroup is not None:
        if nongroup not in group_dict:
            raise ValueError("The 'nongroup' identifier must be a key in the 'group_dict'")
        if exclude_groups is not None:
            raise ValueError("The use of 'nongroup' is incompatible with 'exclude_group'") 
    if exclude_groups is not None:
        if isinstance(exclude_groups,dict):
            gkeys = set(group_dict.keys())
            ekeys = set(exclude_groups.keys())
            if not ekeys == gkeys:
                raise ValueError("Excluded group keys must match the group dict keys")
        else:
            raise ValueError("Excluded group must be dict")


    df_clade = groupFrame.set_index('SampleId')
    keepers = []
    for i in list_frame.index:
        (n1,n2) = i.split(',')
        ref1 = df_clade.loc[n1,'GroupId'] if (n1 in df_clade.index) else None
        ref2 = df_clade.loc[n2,'GroupId'] if (n2 in df_clade.index) else None        
        if exclude_groups is None: ##Groups are tested on exclusion of other groups
            if (ref1 is not None) and (ref2 is not None):
                same_group = (ref1 == ref2) and (ref1 != nongroup) 
                r = list_frame.loc[i].copy()
                r[classification] = same_group
                keepers.append(r)
        else: ##Groups are matched against specific test cases
            if (ref1 is not None) or (ref2 is not None):
                same_group = (ref1 == ref2)
                excluded_pair = (n1 in exclude_groups[ref2]) or (n2 in exclude_groups[ref1])
                if same_group or excluded_pair:
                    r = list_frame.loc[i].copy()
                    r[classification] = same_group
                    keepers.append(r)                
    return pd.DataFrame(keepers)

##Perform Reciever Operator Characteristic analysis for each statistic in list_frame
    ## Returns a dict with two dicts: AUC for each stat; and ROC dataframe for each stat.
##SimilarityMeasure describes the statistics in the list_frame. If False, they are assumed to be distance measures (so that we want to test for differences...
####### This is a bit convoluted...but it works for now
def ROCanalysis(list_frame,classification='SameGroup',SimilarityMeasure=True,TestSet=None):
    test_array = list_frame[classification].astype(bool) if SimilarityMeasure else ~list_frame[classification].astype(bool) #If we're looking at a distance measure, then 'true' means beyond the threshold (exclusion)
    AUC = {}
    test_columns = [x for x in TestSet] if (TestSet is not None) else [x for x in list_frame.columns.tolist() if x != classification]
    ROCFrame = pd.DataFrame()
    for tc in test_columns:
        try:
            fpr, tpr, thresh = metrics.roc_curve(test_array,list_frame[tc])
            if tc.startswith('tpr_'):
                stat_name = tc[4:]
            else:
                stat_name = tc
            fpr_header = '{}_fpr'.format(stat_name)
            tpr_header = '{}_tpr'.format(stat_name)
            thresh_header = '{}_thresh'.format(stat_name)
            ROCFrame_temp = pd.DataFrame({fpr_header:fpr,tpr_header:tpr, thresh_header:thresh}).sort(tpr_header)        
            AUC[stat_name] = metrics.auc(fpr,tpr)
            ROCFrame = ROCFrame.append(ROCFrame_temp)
        except:
            print(tc)
            print(test_array)
            print(list_frame[tc])
            print(list_frame.columns)
            raise

    return {'AUC':AUC,'ROC':ROCFrame}