## requires python 3
## Type -h to see usage instructions (or jump to argparse)

import os
import utilities
from Bio.Alphabet import IUPAC
import sys
unambig = IUPAC.IUPACUnambiguousDNA.letters ##Note: these are all uppercase
unambig_set = set(unambig)
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
import numpy as np

SCRIPT_VERSION = 1.0
SCRIPT_SUBVERSION = 0

script_base = os.path.basename(__file__)
_outputBase = '{}_v{}.{}'.format(os.path.splitext(script_base)[0],SCRIPT_VERSION,SCRIPT_SUBVERSION)
default_logfile = os.path.splitext(script_base)[0]+'.log'

def appendToFilename(filename,insert):
    (base, ext) = os.path.splitext(filename)
    if ext == '.gz':
        (base,ext2) = os.path.splitext(base)
        ext = ext2 + ext
    return base + insert + ext 
    
##New_ext should include the dot
def setExt(filename,new_ext,decompress=False):
    if decompress:
        (name,gz) = os.path.splitext(filename)
        if gz == '.gz':
            filename = name
        else:
            print("Warning: supposedly compressed file needing new extension does not end with gz")
    (base, _) = os.path.splitext(filename)
    if not new_ext[0] == '.':
        new_ext = '.'+new_ext
    return base + new_ext 

##These should probably be lists -- it preserves order and there's no risk of adding duplicates in this routine.
### If I want to apply set logic downstream, I would have to convert to set, perform logic, then convert back to sorted list
## aln is the alignment to evaluate
## maskset are the characters you want to mask; use inverse to mask columns without anything other than these characters
## use case_sensitive you you want the items in maskset to be case_sensitive. Otherwise, everythingwill be cast to upper.
def reportPositionsMaskMonoPoly(aln,maskset,inverse=False,case_sensitive=False):
    masked_sites = defaultdict(int)
    weighted_mask = defaultdict(float)
#     gaps = '-'
#     assert set(gaps) == base_sets['Gaps'] ##For now, I am assuming a single gap character. 
    if not case_sensitive:
        maskset = set([x.upper() for x in maskset])
    if inverse: ##Need to flip around maskset so that it includes everything that is NOT in maskset.
        character_set = set()
        for s in aln:
            character_set.update(set(s.seq))
        if not case_sensitive:
            character_set = set([x.upper() for x in character_set])
        final_maskset = character_set.difference(maskset)
    else:
        final_maskset = maskset
    mask_positions = set()
    mono_positions = set()
    poly_positions = set()
    for i in range(aln.get_alignment_length()):
        col_chars = set(aln[:,i])  
        if not case_sensitive:
            col_chars = set([x.upper() for x in col_chars])
        if (len(final_maskset.intersection(col_chars)) > 0):##report position with anything that is in maskset
            mask_positions.add(i)      
            masked_sequences = [j for j in range(len(aln)) if aln[j,i] in final_maskset]
            for j in masked_sequences:
                masked_sites[j] += 1
                weighted_mask[j] += 1/len(masked_sequences)
        elif len(col_chars) > 1:
            poly_positions.add(i)
        elif len(col_chars) == 1:
            mono_positions.add(i)
        else:
            raise ValueError("Column {} has no characters!".format(i))
    total = len(mask_positions) + len(mono_positions) + len(poly_positions)
    if total != aln.get_alignment_length():
        print("Error")
    return {'mask':mask_positions,'mono':mono_positions,'poly':poly_positions, 'masked_sites':masked_sites,'weighted_mask':weighted_mask}

#position_set is the position to be masked (unless inverse is True, in which case they are kept). Remove means to remove the positions rather than replace them with Ns
def applyMaskToaln(aln,position_set,inverse=False,remove=False):                  
    aln_len = aln.get_alignment_length()
    for i in position_set:
        if not (i <= aln_len -1) and (i >= 0):
            raise IndexError("Gap position is not within legitimate alignment index")
    if inverse:
        drop_list = [x for x in range(aln_len) if x not in position_set] #naturally sorted
    else:
        drop_list = sorted([x for x in position_set])
    start = stop = 0 
    result = {s.id:s[start:stop] for s in aln} ##Initialize 
    #reference = 'Reference'#aln[-1].id
    #dropped = 0
    for i in drop_list: 
        if i > start: 
            gap_len = start - stop
            buffer = '' if remove else 'N'*gap_len 
            stop = i
            for s in aln:
                result[s.id] += buffer + s[start:stop] ## start is beginning of non-gap; everything before i
        start = i + 1 ## start after i (non-gapped position maybe); this will prevent writing of consequtive dropped positions
        ##validate
        #expected_len = i - dropped if remove else i ##rsubtract gaps before this position
#         assert len(result[reference]) == expected_len, "Alignment is {}bp when it should be {}".format(len(result[reference]),expected_len)
        ##record
        #dropped += 1
    for s in aln:#finish by adding terminal sequences
        gap_len = start - stop
        buffer = '' if remove else 'N'*gap_len         
        result[s.id] += buffer + s[start:len(s)]
    ##validate
    #expected_len = aln_len - dropped ##rsubtract gaps before this position
    final_aln = MultipleSeqAlignment(result.values())
    print("Final alignment is {} nucleotides".format(final_aln.get_alignment_length()))
#     assert final_aln.get_alignment_length() == expected_len, "Alignment is {}bp when it should be {}".format(final_aln.get_alignment_length(),expected_len)        
    return final_aln
    
def countNucleotidesAtPositions(aln,seq_id,position_set):
    aln_len = aln.get_alignment_length()
    for i in position_set:
        if not (i <= aln_len -1) and (i >= 0):
            raise IndexError("Gap position is not within legitimate alignment index")
    ref_seq = None
    for s in aln:
        if s.id == seq_id:
            ref_seq = s
    result = defaultdict(int)
    if isinstance(ref_seq,SeqRecord):
        for p in position_set:
            result[ref_seq.seq[p]] += 1
    else:
        print('Error: failed to identify sequence for counting nucleotides')
    return result
                
    

import argparse
def main():
    parser = argparse.ArgumentParser(description='Identify alignment positions that are known in all genomes (core), and replace others with N. Also produce RAxML ascertainment file.')
    ### general info
    parser.add_argument('--version','-V',action='version',version='%(prog)s {}.{}'.format(SCRIPT_VERSION,SCRIPT_SUBVERSION))
#     parser.add_argument('--debug',action='store_true',help="Preserve intermediate files and do not update reference files")

    ### controls

    ### required
    parser.add_argument('input_alignment',help='Alignment that will be masked. The base name will be used for output files')
    parser.add_argument('-o','--out_dir',help='Location to place results. Filenames will be derived from the input alignment')
    parser.add_argument("--no_output",action='store_true',help='Refrain from modifying alignment. Quickly provide alignment description and exit.')

    
    
    args = parser.parse_args()
    try:
        if args.out_dir:
            if os.path.isdir(args.out_dir):
                outdir = args.out_dir
            else:
                outdir = utilities.safeMakeOutputFolder(args.out_dir)
#                 raise ValueError("Draft dir is not a directory")
        else:
            outdir = utilities.safeMakeOutputFolder(_outputBase)
            
        sys.stdout = utilities.Logger('runlog.log')##see locusExtractor for a well-developed example

        logFile = default_logfile  #utilities.appendToFilename(default_logfile, args.projectID) if args.projectID else 
        logFile = os.path.join(outdir,logFile)
        print("LogFile is : "+logFile)

    
        sys.stdout = utilities.Logger(os.path.join(logFile))
    
        print(_outputBase)
        print("Options are:")
        for arg in vars(args):
            print (arg, getattr(args, arg))        
        print()        
#         raise Exception("Test internal")
        aln = AlignIO.read(args.input_alignment,'fasta')
        print("Opened alignment with {} sequences and {} positions".format(len(aln),aln.get_alignment_length()))
        position_types = reportPositionsMaskMonoPoly(aln,unambig_set,inverse=True) ##Flag positions that do not match the unambigous set
        print("Identified the following number of sites:")
        for k in ['mask','mono','poly']:
            print("\t{}:{} ({:.0%})".format(k,len(position_types[k]),100*len(position_types[k])/aln.get_alignment_length()))
        mask_count_index = [position_types['masked_sites'][i] for i in range(len(aln))]
        mask_count = sorted(mask_count_index)
        print("Distribution of masking sites (0, 50, 75, 90, 100): {}, {}, {}, {}, {} ({:.0%})".format(mask_count[1],np.percentile(mask_count,50),np.percentile(mask_count,75),np.percentile(mask_count,90),mask_count[-1],100*mask_count[-1]/aln.get_alignment_length()))#The minimum will  be 0, the reference.
        mask_weight_index = [position_types['weighted_mask'][i] for i in range(len(aln))]
        weight_count = sorted(mask_weight_index)
        print("\tMost masking sites: {}".format(aln[mask_count_index.index(mask_count[-1])].name))
        print("Distribution of weighted masking sites (0, 50, 75, 90, 100): {}, {}, {}, {}, {} ({:.0%})".format(weight_count[1],np.percentile(weight_count,50),np.percentile(weight_count,75),np.percentile(weight_count,90),weight_count[-1],100*weight_count[-1]/len(position_types['mask'])))#The minimum will  be 0, the reference.        
        print("\tMost masking sites (by weight): {}".format(aln[mask_weight_index.index(weight_count[-1])].name))
        #'masked_sites':masked_sites,'weighted_mask'
        if not args.no_output:
            masked_file = os.path.join(outdir,os.path.basename(appendToFilename(args.input_alignment,'_masked')))
            masked_aln = applyMaskToaln(aln,position_types['mask'])
            AlignIO.write(masked_aln,masked_file,'fasta')
            AlignIO.write(masked_aln,setExt(masked_file,'.phy'),'phylip-relaxed')
            reduced_file = os.path.join(outdir,os.path.basename(appendToFilename(args.input_alignment,'_reduced')))
            reduced_aln = applyMaskToaln(aln,position_types['mask'],remove=True)
            AlignIO.write(reduced_aln,reduced_file,'fasta')
            AlignIO.write(reduced_aln,setExt(reduced_file,'.phy'),'phylip-relaxed')    
            partition_base = os.path.join(outdir,os.path.basename(setExt(args.input_alignment,'.txt')))
            ###TODO: format this for RAxML-NG
            partition_guide = appendToFilename(partition_base,'_partition') ##Pass this to RAxML
            counts = countNucleotidesAtPositions(aln,aln[0].id,position_types['mono']) ##Any sequence will do -- these are monomorphic sites 
            partition_data =  appendToFilename(partition_guide,'_data')
            with open(partition_guide,'wt') as partition_out:
                print("[asc~{}], ASC_DNA, p1=1-{}".format(os.path.abspath(partition_data),aln.get_alignment_length()),file=partition_out)
            with open(partition_data,'wt') as data_out:
                print("{} {} {} {}".format(counts['A'],counts['C'],counts['G'],counts['T']),file=data_out)
            partition_NG = appendToFilename(partition_guide,'_NG')
            with open(partition_NG,'wt') as data_out:
                print("GTR+G+ASC_STAM{{{}/{}/{}/{}}}, ALL=1-{}".format(counts['A'],counts['C'],counts['G'],counts['T'],aln.get_alignment_length()),file=data_out)            

        
    

    except (ValueError, IOError) as e:
        print(e)
        parser.print_usage()
#     raise Exception("Test2")
    
    
if __name__ == "__main__":
    if not utilities.has_preferred_python():
        raise Exception("Upgrade your python version")
    main()
