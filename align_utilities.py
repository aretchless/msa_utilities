###Reformat and calculate statistics on alignments

# import os
import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from Bio import AlignIO
from Bio.Alphabet import IUPAC
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
import utilities
# from Bio import AlignIO

# import seq_utilities
unambig = IUPAC.IUPACUnambiguousDNA.letters
unambig_set = set(unambig)
ambig_all_set = set(IUPAC.IUPACAmbiguousDNA.letters)
ambig_only_set = ambig_all_set.difference(unambig_set)
 
base_sets = {
    "DNA4": set([i.lower() for i in unambig_set] + [i.upper() for i in unambig_set]),
    "DNA4_upper": unambig_set,
    "DNA_ambig":set([i.lower() for i in ambig_all_set] + [i.upper() for i in ambig_all_set]),
    "Gaps":set('-')
}

script_version = 1.4 #speed up alignment stats
mode_options = ['fasta_alignment','mauve_xmfa','fasta']





########## Remove information from the alignment (sequences or positions) ##############

### Limit alignment to "name_list" ######
def reduceAlignmentToList(aln_file_in,aln_file_out,name_list):
    seqs = [x for x in SeqIO.parse(aln_file_in,'fasta')]
    reduced = [x for x in seqs if x.id in name_list]    
    SeqIO.write(reduced,aln_file_out,'fasta')
    
#### Limit alignment to Unambiguous (and possibly polymorphic) positions #######
#Note: needs to be tested; poly only    
def stripDownToGATC(aln,upperOnly=False,polyOnly=False):
    keep_set = base_sets['DNA4_upper'] if upperOnly else base_sets['DNA4'] 
    start = stop = 0
    result = aln[:,start:stop] ##Initialize
    keep_region = False
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]            
        col_let = set(col)
        drop_col = len(col_let.difference(keep_set)) > 0 ##Ambiguous
        drop_col |= (polyOnly and len(col_let) == 1) ##Monomorphic
        if drop_col: #has illegit characters
            if keep_region:
                stop = i ##
                result += aln[:,start:stop]
                keep_region = False
        else: #has good characters only
            if keep_region == False:
                start = i
            keep_region = True
    if keep_region == True:
        stop = aln.get_alignment_length()
        result += aln[:,start:stop]
    return result

##This is very slow on large alignments. The way to do it may be to reconstruct the sequences at the end, rather than extending an aln object
#     for s in aln:
#         fh_out.write(">%s\n" % s.id)
#         fh_out.write('%s\n' % ''.join([s.seq[i] for i in keep_cols]))
## Returns a dict with a non-gapped alignment and an integer sets describing the position of gaps;
##  one reports the gap in the alignment index. Follow with convertGapPosToFlankingPos
def removeGapPositions(aln):
    gaps = '-'
    assert set(gaps) == base_sets['Gaps'] ##For now, I am assuming a single gap character. 
    start  = 0
    result = aln[:,start:0] ##Initialize
    gap_positions = set()
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]          
        drop_col = gaps in col  
        if drop_col: #has illegit characters
            result += aln[:,start:i] ## everything before i
            start = i + 1
            ##validate
            expected_len = i - len(gap_positions) ##rsubtract gaps before this position
            assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)
            ### record
            gap_positions.add(i)
    result += aln[:,start:aln.get_alignment_length()]
    ##validate
    expected_len = aln.get_alignment_length() - len(gap_positions) ##rsubtract gaps before this position
    assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)        
    return {'alignment':result,'gap_positions':gap_positions}

def removeGapPositions_set(aln):
    gaps = '-'
    assert set(gaps) == base_sets['Gaps'] ##For now, I am assuming a single gap character. 
    start = 0
    result = aln[:,start:0] ##Initialize
    gap_positions = set()
    for i in range(aln.get_alignment_length()):
        col = set(aln[:,i])          
        drop_col = gaps in col  
        if drop_col: #has illegit characters
            result += aln[:,start:i] ## everything before i
            start = i + 1
            ##validate
            expected_len = i - len(gap_positions) ##rsubtract gaps before this position
            assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)
            ### record
            gap_positions.add(i)
    result += aln[:,start:aln.get_alignment_length()]
    ##validate
    expected_len = aln.get_alignment_length() - len(gap_positions) ##rsubtract gaps before this position
    assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)        
    return {'alignment':result,'gap_positions':gap_positions}

def removePositions(aln,pos_set):
    drop_pos = sorted(list(pos_set))
    for i in drop_pos:
        if not (i <= aln.get_alignment_length() -1) and (i >= 0):
            raise IndexError("Gap position is not within legitimate alignment index") 
    start = stop = 0
    result = aln[:,start:stop] ##Initialize            
    dropped = 0
    for i in drop_pos:
        if i > start:
            result += aln[:,start:i] ## everything before i
        start = i + 1 ## start after i
        ##validate
        expected_len = i - dropped ##rsubtract gaps before this position
        assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)
        ##record
        dropped += 1
    result += aln[:,start:aln.get_alignment_length()]
    ##validate
    expected_len = aln.get_alignment_length() - dropped ##rsubtract gaps before this position
    assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)        
    return result

## Returns a dict with a alignment containing no gaps in the reference, and an integer set describing the position of gaps;
##  one reports the gap in the alignment index, use convertGapPosToFlankingPos to convert that to final 
def removeGapPositionsReference(aln,ref_id):
    gaps = '-'
    assert set(gaps) == base_sets['Gaps'] ##For now, I am assuming a single gap character. 
    start = stop = 0
    result = aln[:,start:stop] ##Initialize
    gap_positions = set()
    ref_index = -1
    for s in range(len(aln)):
        if aln[s].id == ref_id:
            ref_index = s
    for i in range(aln.get_alignment_length()):
        drop_col = aln[ref_index,i] == gaps        
        if drop_col: 
            result += aln[:,start:i] ## everything before i
            start = i + 1
            ##validate
            expected_len = i - len(gap_positions) ##rsubtract gaps before this position
            assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)
            ### record
            gap_positions.add(i)
    result += aln[:,start:aln.get_alignment_length()]
    ##validate
    expected_len = aln.get_alignment_length() - len(gap_positions) ##rsubtract gaps before this position
    assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)        
    return {'alignment':result,'gap_positions':gap_positions}           

##Returns a set reporting any column containing a gap in aln
def reportGapPositions(aln):
    gaps = '-'
    assert set(gaps) == base_sets['Gaps'] ##For now, I am assuming a single gap character. 
    gap_positions = set()
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]          
        if gaps in col:
            gap_positions.add(i)        
    return gap_positions

## scan alignment for positions with characters in charset. If inverse, report positions with characters that are not in charset.
## Note: "inverse" is inverse by character, not by position. i.e. it does not report positions that lack all characters in charset.
def reportPositionsWithChars(aln,charset,inverse=False):
#     gaps = '-'
#     assert set(gaps) == base_sets['Gaps'] ##For now, I am assuming a single gap character. 
    char_positions = set()
    for i in range(aln.get_alignment_length()):
        col_chars = set(aln[:,i])          
        if inverse:
            ##report position with anything that is not in charset
            if len(col_chars.difference(charset)) > 0:
                char_positions.add(i)
        else:
            ##report position with anything that is in charset
            if len(charset.intersection(col_chars)) > 0:
                char_positions.add(i)        
    return char_positions
                    
## Reports gap positions in de-gapped alignment by providing the prior and next base in the final index system   
# aln_length is provided to validate the gap positions and avoid returning the subsequent base if there is a gap in the last position                      
def convertGapPosToFlankingPos(gap_pos,aln_length):
    flanking = set()
    gap_list = sorted(list(gap_pos))
    final_len = aln_length - len(gap_list)
    r = 0 ## removed gaps
    for i in gap_list:
        if not (i <= aln_length -1) and (i >= 0):
            raise IndexError("Gap position is not within legitimate alignment index")
        prior_base = i-r-1
        if prior_base >= 0:
            flanking.add(prior_base)
        ##Catch if there exists a string of gaps at the end of the alignment.
        end_gaps = final_len + r ## if last position is gap, end_gaps == align_len - 1 == i, and so on if there is a run of gaps
        if i < end_gaps: ##If the gap is the last position in the alignment, then there is no next base
            next_base = i-r 
            flanking.add(next_base) 
        r += 1 
    assert r == len(gap_list), "Failed to count all gaps"
    return flanking

##Returns a set reporting any position within flank_size of a gap_pos
def reportFlanks(gap_pos,flank_size,aln_length):
    flanking = set()
    for i in gap_pos:
        if not (isinstance(i,int) and (i <= aln_length -1) and (i >= 0)):
            raise IndexError("Gap position is not within legitimate alignment index")
        start_base = i-flank_size
        stop_base = i+flank_size
        for n in range(start_base,stop_base+1):
            if (n >= 0) and (n < aln_length-1):
                flanking.add(n)
    return flanking

def reportShortSpacers(positions,max_distance,aln_length):
    last = 0 ##By default, treat position 0 as flanking and alignment break
    inter_gap = set()
    for i in sorted(list(positions)):
        if not (isinstance(i,int) and (i <= aln_length -1) and (i >= 0)):
            raise IndexError("Gap position is not within legitimate alignment index")        
        if i - last <= max_distance:
            inter_gap.update(range(last,i+1))
        last = i
    return inter_gap

##Examines alignment (aln) in fragments defined by start_stop list, and returns a set of coordinates that 
### covers any fragment that does not have at least "core_sites" count of gappless positions
#start_stop is a list of pairs (tuple) with first and last coordinate (not what you'd use for slicing)
def maskSmallSegments(aln,start_stop,core_sites,verbose=False):
    masked = set()
    if verbose:
        print("Testing {} fragments".format(len(start_stop)))   
    prev = -1
    for pair in start_stop:
        start = pair[0]
        stop = pair[1]
        if (start <= prev) or (stop < start): ##start may equal stop (1 base fragment)
            raise ValueError("Borders of subalignments must be monotonically increasing")
        if stop >= aln.get_alignment_length():
            raise ValueError("Coordinates must be within the alignment")
    if verbose:
        print("Alignment is {} base pairs".format(aln.get_alignment_length()))
    ### Examine fragments defined by start_stop list coordinates
    kept = 0
    for pair in start_stop:
        start = pair[0]
        stop = pair[1]
        if verbose:
            print("Testing fragment from {} to {}".format(start,stop))        
        core = countCoreInAlignment(aln[:,start:stop+1])
        if verbose:
            print("\tFound {} core sites in fragment of {} basepairs".format(core,stop-start+1))
        if core < core_sites:
            masked.update(range(start,stop+1))
            if verbose:
                print("\tDropping fragment")
        else:
            kept += 1
    if verbose:
        print("Kept is {} fragments".format(kept))            
    return masked

def countCoreInAlignment(aln):
    gaps = '-'
#     assert set(gaps) == base_sets['Gaps'] ##For now, I am assuming a single gap character. 
    core_count = 0
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]          
        if not gaps in col:
            core_count += 1
    return core_count



##Provide two integer sets, both containing indexes in the original alignment. The gap positions will be used to convert the aln_pos into the position in the ungapped alignment
## aln_len is provided to validate the gap positions and avoid returning the subsequent base if there is a gap in the last position and we report flanks
## if aln_pos contains any gap position, it will be dropped unless gaps_to_flanks is true, in which case the flanking positions will be reported
def convertPosToUngapped(aln_pos,gap_pos,aln_length,gaps_to_flanks=False):
    converted = set()
    gap_list = sorted(list(gap_pos),reverse=True) ##use pop to get the lowest item
    aln_list = sorted(list(aln_pos))
    ##Validate
    for x in gap_list + aln_list:
        if not (isinstance(x, int) and (x >=0) and (x <= aln_length - 1)):
            raise IndexError("Positions are not within legitimate alignment index")
    ## count gaps
    prior_gaps = 0
    g = gap_list.pop() ##lowest value
    for a in aln_list:
        while a > g:
            prior_gaps += 1
            try:
                g = gap_list.pop()
            except IndexError: ##No more gaps
                g = aln_list[-1] + 1 ##This is past the last position of interest.
        if a == g:
            if gaps_to_flanks:
                prior_base = a-prior_gaps-1
                if prior_base >= 0:
                    converted.add(prior_base)
                if a < aln_length -1: ##If the gap is the last position in the alignment, then there is no next base
                    next_base = a-prior_gaps 
                    converted.add(next_base)                 
            else:
                print("Warning: position provided is a gap; cannot place it in the ungapped alignment")
        else:       
            converted.add(a-prior_gaps)              
    return converted

########## Provide basic statistics on alingment ########################## 

### Returns a dict with record for each sequence name, and counts for each character
def countBasesInAlignment(aln_list):
    total_counts = defaultdict(defaultdict(int))
    for a in aln_list:
        for s in a:
            name = s.id.split(':')[0]
            counts = Counter(str(s.seq))
            for k,v in counts.items():
                total_counts[name][k] += v
    return total_counts
                
### counts bases (above), then summs counts for each set in the set_dict 
def sumBasesInAlignment(aln_list,set_dict=None):
    if set_dict is None:
        set_dict = base_sets.copy()
    ### Returns a dict with record for each sequence name, and counts for each character
    total_counts = countBasesInAlignment(aln_list)
    
    ### Now sum the counts for each sequence according to character sets
    summed_counts = {}
    for seq_name, seq_counts in total_counts.items():
        ### Setup result dict if needed
        if seq_name not in summed_counts:
            summed_counts[seq_name] = defaultdict(int)
        ### Track uncategorized characters
        remainder_set = set(seq_counts.keys())
        ##Each character set
        for alphabet_name, alphabet_set in set_dict.items():
            ### Count for each set. Sets can be overlapping
            for item in alphabet_set:
                if item in seq_counts:
                    summed_counts[seq_name][alphabet_name] += seq_counts[item]
#                 else:
#                     print("Character {} is not in sequence for {}".format(item,seq_name))
            remainder_set -= alphabet_set
        for item in remainder_set:
            assert item in seq_counts, "Remainder set acquired non-existant item"
            summed_counts['remainder'][remainder_set] += seq_counts[item]
    return summed_counts                
            
# returns the occurances of each variant (e.g. ACGT) in a column of an alignment
### converts all to uppercase
def getVariantSets(col):
    result = defaultdict(set)
    for j in range(len(col)):
        result[col[j].upper()].add(j)
    return result
        
# converts a list of variantSets (from above method) into a SNP table for n isolates
def variantSetList2SNPs(variantSets,n):
    snp_table = [[0 for i in range(n)] for j in range(n)]
    for variants in variantSets:
        keys = [x for x in variants.keys()]
        for i in range(len(keys)):
            set1 = variants[keys[i]]
            for j in range(i):
                set2 = variants[keys[j]]
                for a in set1:
                    for b in set2:
                        snp_table[a][b] += 1
                        snp_table[b][a] += 1   
    return snp_table
    
def alignment2snpBySet(aln,ignore=None,verbose=False,unambig_only=False):
    if ignore == None:
        ignore = []    
    ##Get labels
    labels = [r.id for r in aln]
    variantSets = []
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]
        variants = getVariantSets(col)
        assert len(variants) > 0
        for ex in ignore:#exclude
            if ex in variants.keys():
                del variants[ex]
        if unambig_only:#include
            new_dict = {x:variants[x] for x in variants.keys() if x in unambig_set} ##Uppercase only, but getVariants casts all to upper        
            variants = new_dict 
        if len(variants) > 1:
            variantSets.append(variants)
    if verbose:
        print("Identified {} variant sites.".format(len(variantSets)))
    snp_table = variantSetList2SNPs(variantSets,len(aln))                        
    return pd.DataFrame(snp_table,columns=labels,index=labels,dtype=int)    

# converts a list of variantSets (from above method) into a SNP table for n isolates
def variantSetList2SNPsCluster(variantSets,n,cluster_ignore=0,ignore_edge=False):
    start_site = -1 if ignore_edge else -1 - cluster_ignore ##
    last_site = [[start_site for i in range(n)] for j in range(n)] ##last site is the previous position with a SNP    
    snp_table = [[0 for i in range(n)] for j in range(n)]
    positions = sorted(variantSets.keys())
    for pos in positions:
        variants = variantSets[pos]
        keys = [x for x in variants.keys()]
        for i in range(len(keys)):
            set1 = variants[keys[i]]
            for j in range(i):
                set2 = variants[keys[j]]
                for a in set1:
                    for b in set2:
                        if last_site[a][b] + cluster_ignore < pos:
                            snp_table[a][b] += 1
                            snp_table[b][a] += 1
                        last_site[a][b] = pos #I don't know whether a < b, so fill both
                        last_site[b][a] = pos
    return snp_table
# snp_table = variantSetList2SNPsCluster(variantSets,len(aln),cluster_ignore=cluster_ignore,ignore_edge=ignore_edge)    

def alignment2snpBySetCluster(aln,ignore=None,verbose=False,cluster_ignore=0,ignore_edge=False,unambig_only=False):
    if ignore == None:
        ignore = []    
    ##Get labels
    labels = [r.id for r in aln]
    variantSets = {} ##Store the position of variant sites
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]
        variants = getVariantSets(col)
        assert len(variants) > 0
        for ex in ignore:
            if ex in variants.keys():
                del variants[ex]
        if unambig_only:
            new_dict = {x:variants[x] for x in variants.keys() if x in unambig_set} ##Uppercase only, but getVariants casts all to upper        
            variants = new_dict    
        if len(variants) > 1:
            variantSets[i] = variants
    if verbose:
        print("Identified {} variant sites.".format(len(variantSets)))
    snp_table = variantSetList2SNPsCluster(variantSets,len(aln),cluster_ignore=cluster_ignore,ignore_edge=ignore_edge)                        
    return pd.DataFrame(snp_table,columns=labels,index=labels,dtype=int)  
######### SNP matrix derived from alignment (counts of differences) ##############3
##Count snps in alignment (any character disagreement, not case sensitive), return a dataframe with labels
## pass ['-'] to ignore to ignore gaps

    
def alignment2snpCluster(aln,ignore=None,cluster_ignore=0,ignore_edge=False):
    if ignore == None:
        ignore = []    
    ##Get labels
    labels = [r.id for r in aln]
    ##Count snps
    snp_table = [[0 for i in range(len(aln))] for _ in range(len(aln))]
    start_site = -1 if ignore_edge else -1 - cluster_ignore
    last_site = [[start_site for i in range(len(aln))] for j in range(len(aln))]
    for i in range(aln.get_alignment_length()): #i is position
        col = aln[:,i]
        for j in range(len(col)):
            for k in range(j):
                if not col[j] in ignore and not col[k] in ignore:
                    if col[j].upper() != col[k].upper():
                        if last_site[j][k] + cluster_ignore < i: #i is position, k is always < j
                            snp_table[j][k] += 1
                            snp_table[k][j] += 1  
                        last_site[j][k] = i  # k < j, always
    return pd.DataFrame(snp_table,columns=labels,index=labels,dtype=int)

def alignment2snp(aln,ignore=None):
    if ignore == None:
        ignore = []    
    ##Get labels
    labels = [r.id for r in aln]
    ##Count snps
    snp_table = [[0 for i in range(len(aln))] for j in range(len(aln))]
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]
        for j in range(len(col)):
            for k in range(j):
                if not col[j] in ignore and not col[k] in ignore:
                    if col[j].upper() != col[k].upper():
                        snp_table[j][k] += 1
                        snp_table[k][j] += 1    
    return pd.DataFrame(snp_table,columns=labels,index=labels,dtype=int)

def alignment2InsDelSnp(aln):
    ##Get labels
    labels = [r.id for r in aln]
    ##Count snps
    snp_table = [[0 for i in range(len(aln))] for j in range(len(aln))]
    ins_table = [[0 for i in range(len(aln))] for j in range(len(aln))]
    del_table = [[0 for i in range(len(aln))] for j in range(len(aln))]
    for i in range(aln.get_alignment_length()):
        col = aln[:,i]
        for j in range(len(col)):
            for k in range(j):
                if col[j].upper() != col[k].upper():
                    if col[k] == '-':
                        ins_table[j][k] += 1  ## read row to see that j has an insertion relative to k
                        del_table[k][j] += 1                         
                    elif col[j] == '-':
                        del_table[j][k] += 1
                        ins_table[k][j] += 1
                        pass 
                    else:
                        snp_table[j][k] += 1
                        snp_table[k][j] += 1    
    SNPS = pd.DataFrame(snp_table,columns=labels,index=labels,dtype=int)
    INS = pd.DataFrame(ins_table,columns=labels,index=labels,dtype=int)
    DEL = pd.DataFrame(del_table,columns=labels,index=labels,dtype=int)
    return {'insertions':INS,'deletions':DEL,'SNPs':SNPS}



########### Renaming items ##################            

def findIndiciesInGroup(aln,group=None,useDescription=False):
    if group is None:
        return [x for x in range(len(aln))]
    indicies = []
    for i in range(len(aln)):
        test_value = aln[i].description if useDescription else aln[i].id
        if test_value in group:
            indicies.append(i)        
    return indicies

def alignmentStats(aln,group=None,verbose=False,IsDNA=True):
    if not IsDNA:
        print('Protein comparison not implemented (everything would seem to be ambiguous)')
        return None
    good_indices = findIndiciesInGroup(aln, group)
    gap='Gap'
    ##Basic stats
    result = {'Sites':aln.get_alignment_length(),'Sequences':len(good_indices)}
    ##
    seqs = [aln[i] for i in good_indices]
    if len(seqs) < len(aln):
        newaln = AlignIO.MultipleSeqAlignment(seqs)
    else:
        newaln = aln
    ##Get list of characters
    characters = Counter()
    for s in seqs:
        characters.update(str(s.seq))
    result['Characters'] =characters
    ## identify the gap
    letters = set(characters.keys())
    if '-' in characters:
        result[gap] = '-'
        letters.remove('-')
    else:
        if verbose:
            print("Unable to identify any gap characters")
    ## identify the rest of the alphabet
    if letters <= base_sets["DNA4_upper"]:
        result['Alphabet'] =IUPAC.IUPACUnambiguousDNA
    elif letters <= base_sets["DNA_ambig"]:
        result['Alphabet'] =IUPAC.IUPACAmbiguousDNA
    else:
        print('Failed to identify alphabet. There are probably ambiguous characters, this will affect the outcome of tests')
        print("Casting all characters to uppercase...")
    ambig_letters = letters.difference(base_sets["DNA4_upper"])
    ## characterize site types
    m = g = p = a = 0 #monomorphic, gapped, polymorphic, ambiguous (all sites)
    d = go =  0 #, dimorphic, gap only
    ## gapped sites and ambiguous sites have "unknown" characters, so they preceed assesment of monomorphic and polymorphic
    mono_bases = defaultdict(int)
    for i in range(newaln.get_alignment_length()):
        col_let = set([x for x in aln[:,i].upper()])
        if gap in result and result[gap] in col_let:
            g += 1
            if len(col_let) == 1:
                go += 1
        elif ambig_letters.intersection(col_let) > set():
            a += 1
        elif len(col_let) == 1:
            m += 1
            mono_bases[list(col_let)[0]] += 1
        else:
            p += 1
            if len(col_let) == 2:
                d += 1
    assert (m+g+p+a) == aln.get_alignment_length(), 'Failed to categorize all sites'
    result['Sites-mono'] = m
    result['Mono_bases'] = mono_bases
    result['Sites-gapped'] = g
    result['Sites-poly'] = p
    result['Sites-ambig'] = a
    #
    result['Sites-gapped_only'] = go
    result['Sites-dimorphic'] = d
    return result

def reportFractionWithChars(aln,charset,inverse=False,case_insensitive=False):
    if case_insensitive: #convert all to upper
        charset = set([x.upper() for x in charset])
    result = defaultdict(int)
    for i in range(aln.get_alignment_length()):
        col_chars = Counter(aln[:,i])
        if case_insensitive:
            char_count = defaultdict(int)
            for c in col_chars:
                char_count[c.upper()] += col_chars[c]
        else:
            char_count = col_chars
        count = 0
        if inverse:
            ##report position with anything that is not in charset
            count = sum([v for k,v in char_count.items() if k not in charset])  
        else:
            ##report position with anything that is in charset
            count = sum([v for k,v in char_count.items() if k in charset])  
        result[count] += 1
    return result

##Provide alignment, wiwdow size, and minimum fraction core for each window (0-1)
def slidingWindowTrees(aln,size,min_core):
    start = 0
    stop = size
    final_list = []
    while stop < aln.get_alignment_length():
        w_result = {'start':start,'stop':stop}
        window = aln[start:stop]
        gapless = stripDownToGATC(window)
        gapless_stats = alignmentStats(gapless)
#     result['Sites-gapped'] = g
#     result['Sites-poly'] = p        
        assert gapless_stats['Sites-gapped'] == 0
        aln_length = gapless.get_alignment_length()
        w_result['core_length'] = aln_length
        core_frac = aln_length/size
        if core_frac > min_core:
            poly_sites = gapless_stats['Sites-poly']/aln_length
            w_result['polymorphic'] = poly_sites
            snp_dist = alignment2snp(gapless)
            portion_dist = snp_dist/aln_length
            names = portion_dist.index.tolist()
            dm = _DistanceMatrix(names)
            for i in range(len(dm.matrix)):
                for j in range(i):
                    dm[i,j] =  snp_dist.iloc[i,j]
            tree = DistanceTreeConstructor().nj(dm)           
            w_result['NJ_tree'] = tree            
        else:
            w_result['Note'] = 'Core < minimum ({})'.format(core_frac)
        final_list.append(w_result)
    return pd.DataFrame(final_list)
            
# def writeAlignmentWithCFMask(aln_filename,aln,mask):
#     mask_file = aln_filename + '_masked_sites.txt'
#     nomask_file = aln_filename + '_valid_sites.txt'
#     with open(mask_file,'wt') as mask_out:
#         with open(nomask_file,'wt') as nomask_out:
#             for i in range(aln.get_alignment_length()):
#                 if i in mask:
#                     print(i+1,file=mask_out)
#                 else:
#                     print(i+1,file=nomask_out)
#     print('Masked sites saved to: '+mask_file)
# 
#     print('Valid sites saved to: '+nomask_file )
#     print("Saving full alignment to "+aln_filename)
#     AlignIO.write(aln,aln_filename,'fasta')
#     reduced_aln_file = utilities.appendToFilename(aln_filename, '_MaskRemoved')
#     print("Saving masked alignment to "+reduced_aln_file)
#     reduced = align_utilities.removePositions(aln,mask)
#     AlignIO.write(reduced,reduced_aln_file,'fasta')     
#     return reduced      

def readAligmentMask(mask_file):
    with open(mask_file) as fin:
        mask_list = [int(x)-1 for x in fin.readlines()]
    return set(mask_list)
            
class Alignment_Cleaner:
    def __init__(self,alignment_file):
        self.aln = AlignIO.read(alignment_file,'fasta')
        self.fastafile = alignment_file
        self.phyfile = None
        self.corefastafile = None
        self.stats = None
        self.mask = None
    
    def printStatus(self):
        if isinstance(self.aln,MultipleSeqAlignment):
            print("Loaded alignment of {} sequences with {}bp from {}".format(len(self.aln),self.aln.get_alignment_length(),self.fastafile))
        
    def getAlignmentStats(self):
        if self.stats is None:
            self.stats = alignmentStats(self.aln)
        return self.stats.copy() 
    
    def getMask(self,verbose=False):
        if self.mask is None:
            if verbose:
                print("Creating mask.")
            self.mask = reportPositionsWithChars(self.aln,unambig_set,True)
            if verbose:
                print("Finished creating mask with {} sites".format(len(self.mask)))
        else:
            if verbose:
                print("Found existing mask with {} sites".format(len(self.mask)))
        return self.mask.copy()    
    
    def AlignToPhy(self,phyfile=None):
        if phyfile is None:
            phyfile = utilities.setExt(self.fastafile,'.phy')
        AlignIO.write(self.aln,phyfile,'phylip')
        self.phyfile = phyfile
        return self.phyfile
    
    def WriteCoreUnambigousPositions(self,corefatafile=None,verbose=False):
        if corefatafile is None:
            corefastafile = utilities.appendToFilename(self.fastafile, '_core_unambiguous')
#         if self.corefastafile is None:
        if verbose:
            print("creating core alignment at "+self.corefastafile)
        core_aln = removePositions(self.aln, self.getMask(verbose))
        if verbose:
            print("\tCore alignment has {}bp".format(core_aln.get_alignment_length()))
        AlignIO.write(core_aln,corefastafile,'fasta')
        if verbose:
            print("\tCore alignment written to file.")
#         else:
#             if verbose:
#                 print("Found core alignment at: "+self.corefastafile)
        self.corefastafile = corefastafile
    ##Return dict 
    
    
    def MaskToCFML_file(self,cfml_mask=None,cfml_valid=None):
        if cfml_mask is None:
            cfml_mask = utilities.setExt(utilities.appendToFilename(self.fastafile, '_masked_sites'),'.txt')
        if cfml_valid is None:
            cfml_valid = utilities.setExt(utilities.appendToFilename(self.fastafile,'_valid_sites'),'.txt')
        local_mask = self.getMask()
        with open(cfml_mask,'wt') as mask_out:
            with open(cfml_valid,'wt') as nomask_out:
                for i in range(self.aln.get_alignment_length()):
                    if i in local_mask:
                        print(i+1,file=mask_out)
                    else:
                        print(i+1,file=nomask_out)    
        if local_mask != readAligmentMask(cfml_mask):
            raise IOError("Failed to accurately record mask to mask_file:"+cfml_mask)        
            
                                   
    def reportCompositionOfEach(self):
        values = []
        for s in self.aln:
            bases = Counter(str(s.seq).upper())
            bases['id'] = s.id
            values.append(bases)
        return pd.DataFrame(values).fillna(0).set_index('id')
    
    def makePhymlQsub(self): 
        pass ##TODO
#         qsub_name = 'phyml_qsub_{}.sh'.format(i)
#         print(qsub_name)
#         phylo_dir = os.path.dirname(f)
#         with open(os.path.join(phylo_dir,qsub_name),'wt') as qsub_out:
#             print('#!/bin/bash -l',file=qsub_out)
#             print('#$ -cwd',file=qsub_out)
#             print('module load phyml/3.0',file=qsub_out)
#             print("phyml -i {} -b 500 -s BEST --n_rand_starts 10 --rand_start --run_id {}".format(os.path.basename(f),i),file=qsub_out)
        


# import argparse
def main():
#     print("")
#     print("Running {} from {} at {}".format(SCRIPT_NAME,os.getcwd(),time.ctime()))
#     print("...script and settings are found in {}\n".format(SCRIPT_DIR)) 
#     
#     postCFML_parse = argparse.ArgumentpostCFML_parse(description='A program to consolidate and standardize genome data files.',
#                                      epilog='Settings are defined by the file {}. Each line must have the parameter and then a space. The options are: {}'.format(SETTING_PATH,",".join(setting_options)),
# )
# #                                      ,formatter_class=argparse.MetavarTypeHelpFormatter)
#     ##Info
#     postCFML_parse.add_argument('--version','-V',action='version',version='%(prog)s {}.{}'.format(script_version,script_subversion))
#     postCFML_parse.add_argument('--debug',action='store_true',help='Create a temporary repository in current directory')
#     postCFML_parse.add_argument('--repository','-R',help='Directory containing genome Repository (has e.g. "assemblies" and "reads" as possible subdirectories')
#     subpostCFML_parses = postCFML_parse.add_subpostCFML_parses(description="Select one of the following commands",dest='subcommand')
#     subpostCFML_parses.required = True
#     
#     ##Extract  
#     extract_postCFML_parse = subpostCFML_parses.add_postCFML_parse('extract',description="Retrieve data files from the archive")
#     extract_postCFML_parse.set_defaults(func=extract)
        
#     if mode == mode_options[0]:
#     print('Opening alignment file: ' + aln_file)
#     aln = AlignIO.read(aln_file,'fasta')
#     print('Alignment length: {}'.format(aln.get_alignment_length()))
#     aln_stats = alignmentStats(aln)
#     snpFrame = distanceMatrixCluster(alignment2snp(aln))   
# elif mode == mode_options[1]:
#     pass
# elif mode == mode_options[2]:

    ##Run the program
#     result = args.func(args)
#     if result != 0:
#         postCFML_parse.print_usage()
    pass

if __name__ == '__main__':
    #Get options
    main()
