# msa_utilities 

This software is provided as documentation for the manuscript "Using Neisseria meningitidis genomic diversity to inform outbreak strain identification". It is provided as-is, with no guarantee of usability.

To prepare alignments for phylogenetic analysis, the following two commands were called:

mask_mapped_aln.py: To mask non-core sites in the snippy alignment. This generates an ascertainment bias correction file for RAxML.

adjust_partition_size.py: To adjust the partition size in the ascertainment bias correction file for RAXML, so that it maps properly to the Gubbins filtered SNPs file.

DistanceCalculator.py: To generate distance matrices from kSNP results and from phylogenetic trees produced by RAxML.


