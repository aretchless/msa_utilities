##### XMFA manipulation methods #########
########### Use "MauveHelper" to queue the program ##########

import pandas as pd
import os
import re
import utilities 
import align_utilities
import io
from copy import deepcopy
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment 
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



# import utilities
script_version = 0.4 ##using string for MauveIndex (was int)
script_subversion = 2

#########################################################
####### Extract and parse headers from XMFA file ########

def extractMauveHeader(xmfa_file):
    with open(xmfa_file,'rt') as xmfa:
        header = []
        line = xmfa.readline()
        while line.startswith('#'):
            header.append(line.split())
            line = xmfa.readline()
    return header

###This is for renaming the mauve output; reports file basename rather than path unless fullPath is requested
def getMauveFileID(mauve_header,fullPath=False):
    fileID = {}
    for head in mauve_header:
        seqIDmatch = re.match(r'#Sequence(\d+)File',head[0])
        if seqIDmatch:
            ID = seqIDmatch.group(1)
            assert ID.isdigit(),'Mauve IDs are integers, but not {}'.format(ID)
            fileID[ID] = head[1] if fullPath else os.path.basename(head[1]) 
    return fileID


def getFileIDFromXMFA(XMFA):
    return getMauveFileID(extractMauveHeader(XMFA))

################################################################
########### parse and organize the Mauve XMFA ################
## regex
# example: >2:1-377 
Mauve_name_re = re.compile(r"^(?P<id>\d+):(?P<begin>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<comments>.*)")

### Parse XMFA as header and list of MSAs
def parseXMFA(xmfa_file):
    aln_stream_list = [io.StringIO()]
    with open(xmfa_file,'rt') as xmfa:
        ### Get header with #
        header = []
        line = xmfa.readline()
        while line.startswith('#'):
            header.append(line.split())
            line = xmfa.readline()
        ###Parse alignment sections after header
        while line != '':
            s = aln_stream_list[-1]
            if line.startswith('='):
                if s.tell() > 0: ##Create new stream
                    aln_stream_list.append(io.StringIO()) ##This will append one blank StringIO to end           
            else:
                if (':0-0 + ' in line):##This is the genome placeholder in the first block
                    skip = True
                elif line.startswith('>'): ##Start of a sequence that is NOT a placeholder (already tested for placeholder)
                    skip = False
                if not skip and (len(line.strip()) > 0):   
                    s.write(line)
            line = xmfa.readline()      
    aln_list = []
    for s in aln_stream_list:
        if s.tell() > 0:
            s.seek(0)        
            aln_list.append(AlignIO.read(s,'fasta'))
    return header, aln_list   

### Parse the name given to each sequence in alignment: sequenc id, position, strand
###   position values are adjusted to the python index -- [begin:end] gives the appropriate section of sequence
def parseMauveSeqName(name_line):
    mn_match = Mauve_name_re.match(name_line.strip())
    if mn_match:
        id_str = mn_match.group('id')
        begin_str = mn_match.group('begin')
        end_str = mn_match.group('end')
        strand_char = mn_match.group('strand')
        comments = mn_match.group('comments')
        if strand_char == '+':
            strand = True
        elif strand_char == '-':
            strand = False
        else:
            raise ValueError("Cannot find strand in mauve name: "+name_line)
        try:
            _ = int(id_str) ##just confirming format
            begin = int(begin_str) - 1 ##Mauve index starts with 1; python 0
            end = int(end_str)
        except ValueError:
            print("Could not convert string to int")
            raise
        else:
            seq_id = id_str ##Assign as string
    else:
        raise ValueError("Mauve name cannot be parsed: "+name_line)
    return {'seq_id':seq_id, 'begin':begin, 'end':end, 'strand':strand, 'comments':comments}

### Orients all alignments to reference sequence; discard chunks that do not have reference or have less than "genomes" sequences in alignment
##. Returns a list with COPIES of alignment segments

def orient_to_reference(aln_list,ref_id,verbose=False,genomes=0):
    oriented_aln_list = []
    reoriented = 0 ## If aligened fragment needs to be reoriented
    core = 0 ##If all genomes are in aligned fragment
    for a in aln_list:
        if len(a) == genomes:
            core += 1
        aln_strand = None
        for s in a:
            try:
                name_parts = parseMauveSeqName(s.description)
                seq_id = name_parts['seq_id']
                strand = name_parts['strand']
            except ValueError:
                print("Failure to parse")
            else:      
                if seq_id == ref_id:
                    assert aln_strand == None, "Found two sequences from reference genome"
                    aln_strand = strand
                    if verbose:
                        print(s.name)
        if aln_strand is not None: ###Found the refenece sequence
            if aln_strand:
                new_aln = deepcopy(a)
            else: ###Reverse Complement all strands; makes new objects
                reoriented += 1
                rc_list = []
                for s in a:
                    name_parts = parseMauveSeqName(s.description)
                    strand_char = '-' if name_parts['strand'] else '+'
                    description_line = '{} {} {}'.format(s.name,strand_char,name_parts['comments'])
                    rc_list.append(s.reverse_complement(id=s.id, name=s.name, description=description_line))#new object
                new_aln = MultipleSeqAlignment(rc_list) #this is a deep "copy"
            oriented_aln_list.append(new_aln)
    if verbose:
        print("Reoriented {} out of {} segments with reference genome ({} to start)".format(reoriented,len(oriented_aln_list),len(aln_list)))
        if genomes > 0:
            print("\tthere are {} segments with at least {} sequences".format(core,genomes))
            print("\tFinal list contains {} segments.".format(len(oriented_aln_list)))
    return oriented_aln_list
              
 
### Orders all alignments to reference sequence; discard chunks that do not have reference or have less than "genomes" sequences in alignment
## Returns alignment segments in order of reference (if used after "orient", then all are fresh copies of alignments)


def order_to_reference(aln_list,ref_id,verbose=False,genomes=0):
    core = 0 ##If all genomes are in aligned fragment
    begin_dict = {}
    for a in aln_list:
        if len(a) == genomes:
            core += 1
        aln_index = None
        for s in a:
            try:
                name_parts = parseMauveSeqName(s.description)
                seq_id = name_parts['seq_id']
                begin = name_parts['begin']
            except ValueError:
                print("Failure to parse")
            else:      
                if seq_id == ref_id:
                    assert aln_index == None, "Found two sequences from reference genome"
                    aln_index = begin
                    if verbose:
                        print(s.name)
        if aln_index is not None: ###Found the refenece sequence
            begin_dict[aln_index] = a
    sorted_keys = sorted(begin_dict.keys())
    ordered_aln_list = [begin_dict[k] for k in sorted_keys]
    if verbose:
        if genomes > 0:
            print("\tthere are {} segments with at least {} sequences".format(core,genomes))
            print("\tFinal list contains {} segments.".format(len(ordered_aln_list)))
    return ordered_aln_list  


def concatenate_on_reference(aln_list,ref_id,genome_set,verbose=False):
    if  ref_id not in genome_set:
        raise ValueError("Reference ID is not in the set of genome IDs")
    if verbose:
        print("Reorienting...")
    reoriented = orient_to_reference(aln_list,ref_id, verbose, len(genome_set)) 
    if verbose:
        print("Reordering...")
    reordered = order_to_reference(reoriented, ref_id, verbose, len(genome_set))
    ### Standardize name and concatenate, record LCB edges, remove other gaps and record
    if verbose:
        print("Appending...")
    align_ends = [] ##list of two-tuples
    reference_ends = [] ##list of two-tuples
    full_align = None
    core = 0
    for a in reordered:
        AppendSegment = False
        for s in a:
            try:
                name_parts = parseMauveSeqName(s.description)
                seq_id = name_parts['seq_id']
                if seq_id not in genome_set:
                    raise ValueError("Sequence ID is not in the set of genome IDs")
                begin = name_parts['begin']
                end = name_parts['end']
            except ValueError:
                print("Failure to parse")
            else:   
                ##Rename for concatenation
                s.id = seq_id
                s.name = seq_id   
                ### Test if we have the reference sequence
                if seq_id == ref_id:
#                     if verbose:
#                         print("Parsed seq_id {} from description: {}".format(seq_id,s.description))
                    assert not AppendSegment, "Found two sequences from reference genome"
                    AppendSegment = True
                    if verbose:
                        print("\tAdding {}".format(s.description))
                    reference_ends.append((begin,end-1)) ## Record the exact positions that are excluded, not the python range
        if AppendSegment: ###Found the reference sequence
            ###Make sure the full set is present
            to_concat = deepcopy(a)
            align_set = set([x.id for x in to_concat])
            if align_set == genome_set:
                core += 1
            else:
                gap_seqs = genome_set - align_set
                segment_length = to_concat.get_alignment_length()
                for g in gap_seqs:
                    to_concat.append(SeqRecord(id=str(g),name=str(g),seq=Seq('-'*segment_length)))
            ### Now append   
            to_concat.sort()     
            if full_align is None:
                aln_start = 0
#                 align_ends.add(0)
                full_align = to_concat
            else:
                aln_start = full_align.get_alignment_length()
#                 align_ends.add(full_align.get_alignment_length()) ##Begining of this segment
                full_align += to_concat
#             align_ends.add() ##End of this segment
            aln_stop = full_align.get_alignment_length()-1
            align_ends.append((aln_start,aln_stop))
    result = {
        'alignment':full_align,
        'align_breaks':align_ends,
        'reference_breaks':reference_ends
              }
    return result

################# Convenience functions working from files ###############

### Produces an single MSA from the XMFA file, oriented according to the reference genome (ref_id), but with all gaps. Validates with reference sequence (somewhat).
# Returns a dict with concatenate alignment, 
###   and lists reporting first and last base of each alignment segment -- one indexed to reference, one indexed to the alignment    
#         result = {
#         'alignment':full_align,
#         'align_breaks':align_ends,
#         'reference_breaks':reference_ends
#               }
## Results may be passed to MaskAlginment
def XMFA_concatenate_on_reference_file(xmfa_file,reference_file,verbose=False):
    header, aln_list = parseXMFA(xmfa_file)
    fileID = getMauveFileID(header)
    ref_basename = os.path.basename(reference_file)
    file_lookup = {os.path.basename(filename):i for i,filename in fileID.items()}
    try:
        ref_id = file_lookup[ref_basename]
    except KeyError:
        raise ValueError('XMFA_concatenate: Reference file not found in XMFA header')
    results =  concatenate_on_reference(aln_list,ref_id,set(fileID.keys()),verbose=verbose)
    ##Validate by checking that reference in alignment has same counts of each nucleotide 
    if verbose:
        print("Validating...")
    aln = results['alignment']
    for s in aln:
        if s.name == ref_id:
            aln_count = Counter(str(s.seq))
    seq = SeqIO.read(reference_file,'fasta') 
    count = Counter(str(seq.seq))
    for k,v in count.items(): ##The alignment will have gaps also
        if not v == aln_count[k]:
            print("{} has {} in reference but {} in alignment".format(k,v,aln_count[k]))
            raise ValueError
    if verbose:
        print("Valid")
    return results
        
        
### Produces an single MSA from the XMFA file, using the reference genome (ref_id) to orient, but with gaps. Returns a dict with concatenate alignment, 
###   and lists reporting first and last base of each alignment segment -- one indexed to reference, one indexed to the alignment    
#         result = {
#         'alignment':full_align,
#         'align_breaks':align_ends,
#         'reference_breaks':reference_ends
#               }
def XMFA_concatenate_on_reference_id(xmfa_file,ref_id,verbose=False,reference_file=None):
    header, aln_list = parseXMFA(xmfa_file)
    fileID = getMauveFileID(header)
    if reference_file: ###confirm that reference file matches ref_Id
        ref_basename = os.path.basename(reference_file)
        file_lookup = {os.path.basename(filename):i for i,filename in fileID.items()}
        if ref_basename in file_lookup:
            if ref_id != file_lookup[ref_basename]:
                raise ValueError('XMFA_concatenate: Reference ID does not match the reference file')   
    return concatenate_on_reference(aln_list,ref_id,fileID,verbose=verbose) 

### Produces an single MSA from the XMFA file, using the reference genome (ref_id) as the index. Validates with reference sequence.
def XMFA_index_on_reference_file(xmfa_file,reference_file,verbose=False):
    header, aln_list = parseXMFA(xmfa_file)
    fileID = getMauveFileID(header)
    ref_basename = os.path.basename(reference_file)
    file_lookup = {os.path.basename(filename):i for i,filename in fileID.items()}
    try:
        ref_id = file_lookup[ref_basename]
    except KeyError:
        raise ValueError('XMFA_concatenate: Reference file not found in XMFA header')  
    concatenated =  concatenate_on_reference(aln_list,ref_id,set(fileID.keys()),verbose=verbose)
    aln_raw = concatenated['alignment']   
    ref_break_flanks = concatenated['reference_breaks'] #From Mauve. TODO: I could validate by converting the aln_breaks to ref index, but some of the break flanks are gaps in reference
    aln_break_flanks = concatenated['align_breaks']    
    degapped = align_utilities.removeGapPositionsReference(aln_raw,ref_id)
    aln_ungapped = degapped['alignment']
    aln_removed_gaps = degapped['gap_positions']
    aln_gap_flanks = align_utilities.convertGapPosToFlankingPos(aln_removed_gaps,aln_raw.get_alignment_length()) ##TODO: I would like a way to validate this
    ### Validate
    if verbose:
        print("Validating...")
    ##Validate by checking that reference in alignment matches the raw sequence
    for s in aln_ungapped:
        if s.name == ref_id:
            aln_str = str(s.seq)
    seq = SeqIO.read(reference_file,'fasta') 
    seq_str = str(seq.seq)
    if not aln_str == seq_str:
        raise ValueError("Reference Sequence Reconstructed from XMFA does not match original file ")
    elif verbose:
        print("\tMatch between aln_str ({}bp) and seq_str ({}bp)".format(len(aln_str),len(seq_str)))
    ###Check that the ref_break_flanks could be derived from aln_break_flanks and aln_removed_gaps
    aln_break_set = set([x for pair in aln_break_flanks for x in pair])
    ref_break_set = set([x for pair in ref_break_flanks for x in pair])
    ref_break_converted = align_utilities.convertPosToUngapped(aln_break_set,aln_removed_gaps,aln_raw.get_alignment_length(),True)
    if not (ref_break_converted == ref_break_set):
        raise ValueError("Failed to confirm Mauve's segment break positions with actual XMFA alignment") 
    elif verbose:
        print("\tMatch between alignment break indicies ({} positions) and mauve's indicies ({} positions)".format(len(ref_break_converted),len(ref_break_flanks)))    
    if verbose:
        print("Valid")    
    return {
            'alignment':aln_ungapped,
            'break_flanks':ref_break_flanks,
            'gap_flanks':aln_gap_flanks
            }
        
### Takes the alignment list from an XMFA and concatenates on the provided reference ID. Returns a dict with concatenate alignment, 
###   and lists reporting first and last base of each alignment segment -- one indexed to reference, one indexed to the alignment


##Provide results of XMFA_index_on_reference_file. Result can be used for removing positions from alignment (align_utilities.removePositions(aln,mask))
### or for making a ClonalFrame mask (writeCFmask)
##From Mauve utilities
def maskAlignment(concat_results,verbose=False,segmentMin=5000,break_flank=50,gap_flank=5,inter_gap=30):
    ##Get basic info about gap positions
    aln = concat_results['alignment'] ## 
    b_f = concat_results['break_flanks'] ## report positions flanking edges of alignment fragments
    break_set = set([x for pair in b_f for x in pair])
    g_f = concat_results['gap_flanks']    ##gaps removed from reference; report flanking positions
    gaps = align_utilities.reportGapPositions(aln)
    gaps_1 = align_utilities.reportFlanks(gaps,1,aln.get_alignment_length())
    all_gap_flanks = gaps_1.union(g_f).union(break_set)
    ##Run
    frag_masking = align_utilities.maskSmallSegments(aln,b_f,segmentMin,True)
    break_masking = align_utilities.reportFlanks(break_set,break_flank,aln.get_alignment_length())
    gap_masking_1 = align_utilities.reportFlanks(g_f,gap_flank-1,aln.get_alignment_length()) ##Prior masking of one on each side
    gap_masking_0 = align_utilities.reportFlanks(gaps,gap_flank,aln.get_alignment_length()) ##No prior masking
    gap_masking = gap_masking_1.union(gap_masking_0)
    between_gaps = align_utilities.reportShortSpacers(all_gap_flanks,inter_gap,aln.get_alignment_length())
    ##Final
    result = gap_masking.union(break_masking).union(frag_masking).union(between_gaps)
    if verbose:
        print("Masking a total of {} sites".format(len(result)))
        print("\t{} gaps".format(len(gaps)))
        prior_masks = gaps.copy()
        print("\t{} on short fragments (minimum length of {})".format(len(frag_masking.difference(prior_masks)),segmentMin))
        prior_masks.update(frag_masking)
        print("\t{} near edges of fragments ({}bp)".format(len(break_masking.difference(prior_masks)),break_flank))
        prior_masks.update(break_masking)
        print("\t{} near gaps ({}bp)".format(len(gap_masking.difference(prior_masks)),gap_flank))
        prior_masks.update(gap_masking)
        print("\t{} in clusters of gaps (within {}bp of each other)".format(len(between_gaps.difference(prior_masks)),inter_gap))
              
    return result

##Replaces "pos_set" positions in aln with mask_char
def maskPositions(aln,pos_set,mask_char='-'):
    mask_pos = sorted(list(pos_set))
    for i in mask_pos:
        if not (i <= aln.get_alignment_length() -1) and (i >= 0):
            raise IndexError("Gap position is not within legitimate alignment index") 
    start = stop = 0
    result = aln[:,start:stop] ##Initialize     
    gap_size = 0
    for i in mask_pos:
        if i > start:
            ##Add gap before adding next chunk
            gap_aln = aln[:,start:stop] 
            mask_str = gap_size * mask_char
            for s in gap_aln:
                s.seq = Seq(mask_str) 
            result += gap_aln
#             print('added gap of length {},{}'.format(gap_aln.get_alignment_length(),gap_size))
            gap_size = 1 ##since i is a gap position
            ##Add next chunk
            result += aln[:,start:i] ## everything before i
        else:
            gap_size += 1
#         result += gap_aln
        start = i + 1 ## start after i
        ##validate
        expected_len = start - gap_size ##rsubtract gaps before this position
        assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)

    ##Add final gap:
    gap_aln = aln[:,start:stop] 
    mask_str = gap_size * mask_char
    for s in gap_aln:
        s.seq = Seq(mask_str) 
    result += gap_aln    
    ##Add final sequence
    result += aln[:,start:aln.get_alignment_length()]
    ##validate
    expected_len = aln.get_alignment_length()  ##rsubtract gaps before this position
    assert result.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(result.get_alignment_length(),expected_len)        
    return result    
# Write masking file for CF
##Takes mask from maskAlignment and makes file for ClonalFrameML
#The mask file is indexed on 1
#### full_aln_file = os.path.join(phylo_dir,'Alignment_IndexedOn_{}.aln'.format(refM))
##Note: this should probably be deprecated; replace with method in Align_Utilities (or Align_Cleaner class)
from warnings import warn
def writeAlignmentWithCFMask(aln_filename,aln,mask):
    warn("Use Align_utilities, not Mauve_utilities",FutureWarning)    
    mask_file = aln_filename + '_masked_sites.txt'
    nomask_file = aln_filename + '_valid_sites.txt'
    with open(mask_file,'wt') as mask_out:
        with open(nomask_file,'wt') as nomask_out:
            for i in range(aln.get_alignment_length()):
                if i in mask:
                    print(i+1,file=mask_out)
                else:
                    print(i+1,file=nomask_out)
    print('Masked sites saved to: '+mask_file)
    if mask != readAligmentMask(mask_file):
        raise IOError("Failed to accurately record mask to mask_file:"+mask_file)
    print('Valid sites saved to: '+nomask_file )
    print("Saving full alignment to "+aln_filename)
    AlignIO.write(aln,aln_filename,'fasta')
    reduced_aln_file = utilities.appendToFilename(aln_filename, '_MaskRemoved')
    print("Saving masked alignment to "+reduced_aln_file)
    reduced = align_utilities.removePositions(aln,mask)
    AlignIO.write(reduced,reduced_aln_file,'fasta')     
    return reduced               

def readAligmentMask(mask_file):
    with open(mask_file) as fin:
        mask_list = [int(x)-1 for x in fin.readlines()]
    return set(mask_list)
        
###This will suggest how to rename the sequences. Gives priority to 
#### PubMLST identifier
#### CDC BML identifier
#### fall back on file basename (which is provided by first line)
# This should be used before renaming files for phyml and RAxML -- if sequences were renamed with "seq" prefix, then 
def getXMFAnameKey(XMFA):
    fileID = getFileIDFromXMFA(XMFA)
    mauve2lookup = {}
    for k,v in fileID.items():
        newkey = k
        pub_match = re.search(r'_(PubMLST\d+)',v)
        m_match = re.match(r'(M\d{5})',v)
        if pub_match:
            lookup = pub_match.group(1)
        elif m_match:
            lookup = m_match.group(1)
        else:
            lookup = v ##filename
            print("Warning: failed to infer lookup name for file:"+v) 
        mauve2lookup[newkey] = lookup
    return mauve2lookup

    
###########Parse data files #####3############3                        

##Converts Mauve SNP file to matirx with headers named after each file
##Note: fileID=None is not tested
def parseMauveSNPs(snp_file,fileID=None):
    if fileID is None:
        fileNames = None
    else:
        fileNames = len(fileID) * [None]
        for x,y in fileID.items():
            fileNames[x-1] = y
    snps = pd.read_table(snp_file)
    snpsRows = [[c for c in x] for x in snps['SNP pattern'].tolist()]
    snpsFrame = pd.DataFrame(snpsRows,columns=fileNames)
    genomeCols = [x for x in snps.columns if "GenWidePos" in x]
    genomeColsReplace = {}
    for x in genomeCols:
        genomeMatch = re.match(r'sequence_(\d+)_GenWidePos\d+',x)
        if genomeMatch:
            ID = int(genomeMatch.group(1))
            genomeColsReplace[x] = fileID[ID]
    snpPos = snps[genomeCols].rename(columns=genomeColsReplace)
    return snpsFrame, snpPos

##Returns "left" and "right" position for each backbone fragment in BBone_file; uses names given by fileID
## if position is 0, the fragment is missing from that genome
## The length of the fragment is right-left +1 (but only for those where right and left are not zero)
def parseMauveBBone(BBone_file,fileID=None):
    bbone = pd.read_table(BBone_file)
    r_cols = [c for c in bbone.columns if c.endswith('rightend')]
    l_cols = [c for c in bbone.columns if c.endswith('leftend')]
    bbone_right = bbone[r_cols].copy()
    bbone_left = bbone[l_cols].copy()
    ##
    if fileID is not None:
        r_rename ={}
        for c in r_cols:
            rightMatch = re.match(r'seq(\d+)_rightend',c)
            if rightMatch:
                rID = int(rightMatch.group(1))
                r_rename[c] = fileID[rID+1]
        bbone_right.rename(columns=r_rename,inplace=True)
        ##
        l_rename ={}
        for c in l_cols:
            leftMatch = re.match(r'seq(\d+)_leftend',c)
            if leftMatch:
                lID = int(leftMatch.group(1))
                l_rename[c] = fileID[lID+1]
        bbone_left.rename(columns=l_rename,inplace=True)
    return bbone_left, bbone_right


############ SNP filtering and stats ################# 
##Parsing the raw matrix to identify ambiguities and gaps
def MauveSNPStoGapsAndFlanks(snpsFrame,snpPos):
    charGaps = (snpsFrame == '-').any(axis=1)
    posGaps = (snpPos == 0).any(axis=1)
    assert (charGaps == posGaps).all(), "Mauve DataFrames do not agree on position of gaps"
    near_indel = pd.rolling_sum(posGaps, window=3, center=True ) >= 1
    return snpPos[near_indel]

# def reportGapSizes(snpPos):
#     posGaps = (snpPos == 0).any(axis=1)
                          
######### BBone stats #################


########## Assign reference sequences ###############
def MauveSNPStoRefSNPS(snpsFrame,snpPos,refID,padGaps=0):
    snpRef = snpsFrame.copy()
    ##Check that frames agree on position of gaps
    charGaps = (snpRef == '-').any(axis=1)
    posGaps = (snpPos == 0).any(axis=1)
    assert (charGaps == posGaps).all(), "Mauve DataFrames do not agree on position of gaps"
    near_indel = pd.rolling_sum(posGaps, window=2*padGaps+1, center=True ) >= 1 ##Search padGaps distance on each side    
    ##This drops rows with ambiguous data and gaps
    snpRefChar = snpRef.isin(['G','A','T','C']).all(axis=1)
    snpRef['position'] = snpPos[refID]
    snpRef = snpRef[snpRefChar & ~near_indel].copy()
    snpRef.set_index('position',inplace=True)
    return snpRef

##Treat gaps like any other character. Still drops ambiguities
def MauveSNPStoRefSNPSKeepGaps(snpsFrame,snpPos,refID):
    snpRef = snpsFrame.copy()
    ##This drops rows with ambiguous data and gaps
    snpRefChar = snpRef.isin(['G','A','T','C','-']).all(axis=1)
    snpRef['position'] = snpPos[refID]
    snpRef = snpRef[snpRefChar].copy()
    snpRef.set_index('position',inplace=True)
    return snpRef

##Todo: this needs better results
def compareCladeBBone(bbone_left,bbone_right,ingroup,outgroup=None): #ingroup_name (dropped argument)
    ##Calculate size of backbone fragments
    bbone_present = bbone_right > 0
    bbone_size = bbone_right - bbone_left + bbone_present
    if outgroup == None:
        outgroup = list(set(bbone_left.columns.tolist()).difference(set(ingroup)))
    presenceF = bbone_right != 0
    in_all = presenceF[ingroup].all(axis=1)
#     in_none = ~presenceF[ingroup].any(axis=1)
    out_all = presenceF[outgroup].all(axis=1)
    out_none = ~presenceF[outgroup].any(axis=1)   
    in_only = in_all & out_none     
#     out_only = out_all & in_none
    return bbone_size[in_only][ingroup + outgroup]
    
###The way Mauve uses gaps in SNP file undermines these statistics: Mauve only reports gaps if there is a SNP at that location
##Modify this to treat gaps like any other character. Also keeps ambiguities
def compareCladeSNPSAndGaps(seqFrame,ingroup,ingroup_name,outgroup=None):
    if outgroup == None:
        outgroup = list(set(seqFrame.columns.tolist()).difference(set(ingroup)))
    discriminatory_snps = set() ##monomorphic within clade; absent outside
    substitution_snp = set() ## site distinguishes ingroup from outgroup
    hypervariable = set()
    homoplasy_snps = set() ##multiple alleles shared by both in and out
    dimorphic = set()
    for p,row in seqFrame.iterrows():
        clade_states = set(row[ingroup].unique().tolist())
        assert len(clade_states) > 0
#         clade_unamb = clade_states.intersection(unambig_set)
        out_states = set(row[outgroup].unique().tolist())
        assert len(out_states) > 0
#         out_unamb = out_states.intersection(unambig_set)
        overlap_states = clade_states.intersection(out_states)
        all_states = clade_states.union(out_states)
        if (len(all_states) == 2) and (len(all_states.intersection(align_utilities.unambig_set)) == 2):
            dimorphic.add(p)
        if len(overlap_states) == 0:
            substitution_snp.add(p)
            if len(clade_states) == 1:
                discriminatory_snps.add(p)
            else:
                hypervariable.add(p)
        elif len(overlap_states) > 1:
            homoplasy_snps.add(p)
    node_summary = {
        'CladeID':ingroup_name,
        'MemberCount':len(ingroup),
        'DiscriminatorySNPs':len(discriminatory_snps),
        'Substitutions':len(substitution_snp),
        'Hypervariable':len(hypervariable),
        'Homoplasies':len(homoplasy_snps),
        'OutgroupSize':len(outgroup),
        'DimorphicDicsriminatory':len(discriminatory_snps.intersection(dimorphic)),
        'DimorphicHomoplasy':len(homoplasy_snps.intersection(dimorphic)),
        'DimorphicTotal':len(dimorphic)
    }
    clade_details = {
        'CladeID':ingroup_name,
        'Disc':discriminatory_snps,
        'Homo':homoplasy_snps
    }
    return node_summary, clade_details




### One-step convenience functions ########
###Returns SNPS and Positions with proper names
def parseMauveSNPsAndHeader(SNPS,XMFA):
    FileIDs = getMauveFileID(extractMauveHeader(XMFA))
    return parseMauveSNPs(SNPS,fileID = FileIDs)

def createRefSNPsFromMauve(XMFA,SNPS,refFile):
    snps,pos = parseMauveSNPs(SNPS,getFileIDFromXMFA(XMFA))
    return MauveSNPStoRefSNPS(snps,pos,refFile) 

##Returns left and right with appropirate names
def parseMauveBBoneWithHeader(BBONE,XMFA):
    FileIDs = getFileIDFromXMFA(XMFA)
    return parseMauveBBone(BBONE,FileIDs)
    

import argparse
def main():
    parser = argparse.ArgumentParser(description='A program to convert Mauve SNP output into a reference-oriented table.')
    ##Info
    parser.add_argument('--version','-V',action='version',version='%(prog)s {}.{}'.format(script_version,script_subversion))
    parser.add_argument('xmfa_file',help='The the XMFA alignment from Mauve (will use the header)')
    parser.add_argument('snp_file',help='The the SNP file exported from Mauve')
    parser.add_argument('reffile',help='The file from the XMFA alignment to use as the reference coordinate')
    parser.add_argument('outfile',help='Location to write result file.')
    args = parser.parse_args()    
    refFile = os.path.basename(args.reffile)
    refSNPFrame = createRefSNPsFromMauve(args.xmfa_file,args.snp_file,refFile)
    refSNPFrame.to_csv(args.outfile,sep='\t')
     
if __name__ == "__main__":
    if not utilities.has_preferred_python():
        raise Exception("Upgrade your python version")
    main()