# msa_utilities 

This software is provided as documentation for the manuscript "Using Neisseria meningitidis genomic diversity to inform outbreak strain identification". It is provided as-is, with no guarantee of usability.

To prepare alignments for phylogenetic analysis, the following two scripts were called:

mask_mapped_aln.py: This generates a masked alignment file from the alignment file produced by snippy-clean_full_aln. It also produces an partition file for RAxML providing informatino ascertainment bias correction.

mask_mapped_aln.py alignment_file

adjust_partition_size.py: To adjust the partition size in the ascertainment bias correction file for RAXML, so that it maps properly to the Gubbins filtered polymorphic sites file.

adjust_partition_size.py gubbins_filtered_polymorphic_sites partition_file

To create distance matricies from kSNP results and from phylogenetic trees produced by RAxML.

DistanceCalculator.py kSNP kSNP_result_directory
or
DistanceCalculator.py tree tree_file



