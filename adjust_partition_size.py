## requires python 3
## Type -h to see usage instructions (or jump to argparse)

import os
import utilities
from Bio import AlignIO

SCRIPT_VERSION = 0.1 
SCRIPT_SUBVERSION = 2

script_base = os.path.basename(__file__)
_outputBase = '{}_v{}.{}'.format(os.path.splitext(script_base)[0],SCRIPT_VERSION,SCRIPT_SUBVERSION)
default_logfile = os.path.splitext(script_base)[0]+'.log'
                
def appendToFilename(filename,insert):
    (base, ext) = os.path.splitext(filename)
    if ext == '.gz':
        (base,ext2) = os.path.splitext(base)
        ext = ext2 + ext
    return base + insert + ext     

import argparse
def main():
    parser = argparse.ArgumentParser(description='Modifies (in place) a RAxML partition file if using a different alignment than what was used my mask_mapped_aln.')
    ### general info
    parser.add_argument('--version','-V',action='version',version='%(prog)s {}.{}'.format(SCRIPT_VERSION,SCRIPT_SUBVERSION))
#     parser.add_argument('--debug',action='store_true',help="Preserve intermediate files and do not update reference files")

    ### controls
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--is_NG',action='store_true',help='partition_file is formatted for RAxML NG; otherwise it is formatted for v8')
    group.add_argument('--find_NG',action='store_true',help='partition_file is formatted for RaxML v8; also modify a NG-formatted file with "_NG" at the end of the filename')
    ### required
    parser.add_argument('input_alignment',help='Alignment that will be used by RAxML. Expects phylip format')
    parser.add_argument('partition_file',help='Partition file will be modified in place to reflect the length of the given alignment.')

    
    
    args = parser.parse_args()
    try:
    
        print(_outputBase)
        print("Options are:")
        for arg in vars(args):
            print (arg, getattr(args, arg))        
        print()        
#         raise Exception("Test internal")
        aln = AlignIO.read(args.input_alignment,'phylip-relaxed')
        print("Opened alignment with {} sequences and {} positions".format(len(aln),aln.get_alignment_length()))
        NG_file = None
        if not args.is_NG:
            with open(args.partition_file) as fin:
                partition_lines = fin.read().splitlines()
            if len(partition_lines) > 1:
                print("Error: this file should have only one line")
                print(partition_lines[0])
                print(partition_lines[1])
            else:
                print("Starting with {}".format(partition_lines[0]))
                parts = partition_lines[0].split(',') ##  created with "[asc~{}], ASC_DNA, p1=1-{}"
                if len(parts) != 3:
                    print("Error: the first line should have three comma-delimited parts")
                else:
                    parts[2] = 'p1=1-{}'.format(aln.get_alignment_length())
                    new_line = ', '.join(parts)
                    with open(args.partition_file,'wt') as fout:
                        print(new_line,file=fout)
                    print("Updated the partition file: {}".format(args.partition_file))
                    print("Current text: {}".format(new_line))
        ##RAxML NG file            
        if args.find_NG:
            NG_file = appendToFilename(args.partition_file,'_NG')
        if args.is_NG:
            NG_file = args.partition_file
        if os.path.isfile(NG_file):
            try:
                with open(NG_file) as fin:
                    ng_lines = fin.read().splitlines()
                if len(ng_lines) > 1:
                    print("Error: NG file should have only one line")
                    print(ng_lines[0])
                    print(ng_lines[1])
                else:
                    print("Starting with {}".format(ng_lines[0]))
                parts = ng_lines[0].split(',') ## created with "GTR+G+ASC_STAM{{}/{}/{}/{}}, ALL=1-{}"
                if len(parts) != 2:
                    print("Error: the first line of the NG file should have two comma-delimited parts")
                else:
                    parts[1] = 'ALL=1-{}'.format(aln.get_alignment_length())
                    new_line = ', '.join(parts)
                    with open(NG_file,'wt') as fout:
                        print(new_line,file=fout)
                    print("Updated the NG partiton file: {}".format(NG_file))
                    print("Current text: {}".format(new_line))
            except IOError:
                print("Failure to open NG partition file")
                    
                        
    except (ValueError, IOError) as e:
        print(e)
        parser.print_usage()
#     raise Exception("Test2")
        
                #print(.format(counts['A'],counts['C'],counts['G'],counts['T'],aln.get_alignment_length()),file=data_out)   
    
    
if __name__ == "__main__":
    if not utilities.has_preferred_python():
        raise Exception("Upgrade your python version")
    main()
