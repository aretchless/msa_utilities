###For working with alignment tables (Lyve-set out.filteredbcftoolsquery.tsv) or out.pooled.vcf.gv

from Bio.Alphabet import IUPAC
unambig = IUPAC.IUPACUnambiguousDNA.letters
unambig_set = set(unambig)
ambig_all_set = set(IUPAC.IUPACAmbiguousDNA.letters)
ambig_only_set = ambig_all_set.difference(unambig_set)
from collections import defaultdict
import pandas as pd

#(Lyve-set out.filteredbcftoolsquery.tsv)
P = 'POS'
R = 'REF'
    
def openFilteredBCFToolsQuery(tsv_file):
    tsvFrame = pd.read_table(tsv_file)
    core_cols = [x for x in tsvFrame.columns if ']' in x]
    tsvFrame = tsvFrame[core_cols].copy()
    cols = [x.split(':')[0].strip('#') for x in core_cols]
    cols2 = [x.split("]")[1] for x in cols]
    tsvFrame.columns = cols2
    return tsvFrame

def indexByPosition(tsvFrame):
    names = tsvFrame.columns.tolist()    
#     indicies = ['CHROM', 'POS', 'REF']
    names = names [3:]
    seqFrame = tsvFrame.set_index('POS')[names]
    return seqFrame
    
##Columns can only be sequences
def removeAmbigSites(tsvFrame,Max=0):
    keepers = []
    for _,row in tsvFrame.iterrows():
        counts = row.value_counts()
        c = 0
        for key, value in counts.items():
            if key in ambig_only_set:
                c += value
        if c <= Max:
            keepers.append(row)
    return pd.DataFrame(keepers)    
    
def countAmbiguousSitesByGenome(tsvFrame):
    counts = []
    for c in tsvFrame.columns:
        nucCounts = tsvFrame[c].value_counts()
        nucCounts['ID'] = c
        counts.append(nucCounts)
    return pd.DataFrame(counts).set_index('ID')

###All columns should be sequence columns
def countSiteTypes(tsvFrame):
    counts = defaultdict(int)
    for p,row in tsvFrame.iterrows():
        values = set(row.unique().tolist())
        if len(values) == 1:
            counts['Mono'] += 1
        else:
            unambig = values.difference(ambig_only_set)
            if len(unambig) == 1:
                counts['Mono-Amb'] += 1
                counts['Mono'] += 1
            else:
                counts['Poly'] += 1
                ambiguous = values.difference(unambig_set)
                if len(ambiguous) > 0:
                    counts['Poly-Amb'] += 1
    return counts    

##This is designed for a DataFrame that includes nucleotide characters and has a single reference position.
##Looks for SNPS that distinguish the ingroup from the outgroup. 
##   some results will be nonsensical if the ingroup is not monophyletic with respect to the outgroup
def compareCladeSNPS(seqFrame,ingroup,ingroup_name,outgroup=None):
    if outgroup == None:
        outgroup = list(set(seqFrame.columns.tolist()).difference(set(ingroup)))
    discriminatory_snps = [] ##monomorphic within clade; absent outside
    substitution_snp = [] ## site distinguishes ingroup from outgroup
    hypervariable = []
    homoplasy_snps = [] ##multiple alleles shared by both in and out
    for p,row in seqFrame.iterrows():
        clade_states = set(row[ingroup].unique().tolist())
        clade_unamb = clade_states.intersection(unambig_set)
        out_states = set(row[outgroup].unique().tolist())
        out_unamb = out_states.intersection(unambig_set)
        overlap_states = clade_unamb.intersection(out_unamb)
        if len(overlap_states) == 0:
            substitution_snp.append(p)
            if len(clade_unamb) == 1:
                discriminatory_snps.append(p)
            else:
                hypervariable.append(p)
        elif len(overlap_states) > 1:
            homoplasy_snps.append(p)
    node_summary = {
        'CladeID':ingroup_name,
        'MemberCount':len(ingroup),
        'DiscriminatorySNPs':len(discriminatory_snps),
        'Substitutions':len(substitution_snp),
        'Hypervariable':len(hypervariable),
        'Homoplasies':len(homoplasy_snps),
        'OutgroupSize':len(outgroup)
    }
    clade_details = {
        'CladeID':ingroup_name,
        'Disc':discriminatory_snps,
        'Homo':homoplasy_snps
    }
    return node_summary, clade_details

##Eliminate anything that is only polymorphic because of ambiguities
def reduceToTruePoly(tsvFrame):
    keepers = []
    for p,row in tsvFrame.iterrows():
        values = set(row.unique().tolist())
        unambig = values.difference(ambig_only_set)
        if len(unambig) > 1:
            keepers.append(row)
    return pd.DataFrame(keepers)

    

### Lyve-set out.pooled.vcf.gz


def openVCFfile(vcf_file,returnHeader = False):
    import gzip
    vcf_in = gzip.open(vcf_file,'rt')
    header = []
    cline = vcf_in.readline()
    while cline.startswith('##'):
        header.append(cline)
        cline = vcf_in.readline() ##will stop on header, which only has single #
    names = cline.lstrip('#').split()
    total_vcf = pd.read_table(vcf_in,names=names)
    if returnHeader:
        result = (total_vcf,header)
    else:
        result = total_vcf
    return result

##TODO: return results
def classifySitesInVCF(vcfFrame):
    ###stats from a vcf_file ... just counting bases for ascertainment correction in raxML
    mono_base_count = defaultdict(int)
    a = m = p = c = g = 0
    e = 0
    match_set = set('.') ##What VCF4.1 spec says should be in the ALT column to specify the reference alllele
    tempFrame = vcfFrame#[0:1000]
    for pos,group in tempFrame.groupby('POS'):
        e += 1
        if len(group) == 1:
            for _, row in group.iterrows():
                ref = row['REF'] #.upper?
                alt = row['ALT'] # .upper?
                if (len(ref) == 1):
                    ref_set = set(ref.split(','))
                    assert len(ref_set) == 1, 'Reference should only have a single allele at pos {}'.format(pos)
                    alt_set = set(alt.split(','))
                    ##Note: ALT may have multi-character states... that's a bit ambiguous to interpret
                    if len(alt_set) > 1:
                        temp_set = set()
                        for alt_allele in alt_set:
                            if len(alt_allele) > 1:
                                temp_set.add(alt_allele[0])
    #                            print("Multi-character ALT at pos {}: A={} R={}".format(pos,alt_allele,ref_set.copy().pop())) 
                            else:
                                temp_set.add(alt_allele)
                        alt_set = temp_set
                    if len(alt_set.difference(unambig_set)) > 0: #Test if ALT has anything other than GACT
                        a += 1 #Ambiguous
                    elif alt_set == ref_set: ##Ref_set can only have one allele (see conditionals and assertions above)
                        m += 1 #Monomorphic 
                        mono_base_count[ref] += 1
                    else: ##Polymorphic
                        p += 1
                        
                else:
                    g += 1 #indel
        else:
            c += 1 ##Complex variantion
    
    total = (a + m + p + g + c) 
    print('Final: a {}, m {}, p {}, g {}, c {}'.format(a,m,p,g,c))  
    expected = sum(~tempFrame.POS.duplicated()) 
    assert total == expected ,'Failed to count everything: Total {}; Expected {}; e {}'.format(total,expected,e)
    return mono_base_count
  