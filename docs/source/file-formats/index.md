# File formats

TODO: for each format give a description, a small worked example, column
definitions, units, ordering rules, and any assumptions the tools make.

## Copy-number input

A file containing integer absolute/relative copy numbers for all the patient samples and/or the normal sample (*-cn.txt.gz or *-haplotype-cn.txt.gz).

Either compressed file or uncompressed file is fine. There need to be at least four columns, separated by space, in this file: sample_ID, chr_ID, site_ID, CN. Each column is an integer. Note that there should be __no header names__ in this file.

The sample_ID has to be __ordered from 1 to the number of patient samples__.

The chr_ID and site_ID together determine a unique site along the genome of a sample, ordering from 1 to the largest number (1, 2, 3, ...).

The site_ID can be __consecutive numbers__ from 1 to the total number of sites along the genome, or __consecutive numbers__ from 1 to the total number of sites along each chromosome of the genome.

For __haplotype-specific__ CN, there need to be at least five columns, with the last two being cnA, cnB.

If the total CN is larger than the specified maximum CN allowed by the program, the total CN will be automatically decreased to the maximum CN when the input is total CN and the program will exit when the input is haplotype-specific CN.

When the input copy numbers are relative with normal copy being 0 as those output by CGHcall, please specify it with option "--is_rcn 1".

When the input copy numbers are haplotype-specific which have been scaled relative to ploidy or not, please specify it with option "--is_total 0 --is_rcn 0".

When the input copy numbers are in bins of fixed size, please specify it with option "--bin 0" to use original data. By default "--bin 1" is used to get segment-level data by merging consecutive bins with the same copy number in a sample with change points aligned across all the samples.

## Sample timing files

TODO.

## Tree formats

TODO: Newick conventions, branch length units, rooting.

## Mutation outputs

TODO.

## MCMC configuration and trace files

See `mcmc.cfg` in the repository root. TODO.

## Identifiers and conventions

TODO: sample naming, chromosome naming, 0- vs 1-based coordinates.
