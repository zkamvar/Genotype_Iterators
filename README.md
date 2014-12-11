## Genotype Iteration

This repository contains scripts I used to test implementation of Bruvo's 
distance in C.

It also contains a script called `genome_addition_filler.c`, which is a scirpt
that demonstrates the imputation method for the genome model for Bruvo's 
distance. It also demonstrates a generalization of the algorithm to be able to
create DNA models for genotypes of any ploidy. 

These are fully implemented in the [poppr R package](https://github.com/grunwaldlab/poppr)

## Running

You can test out `genome_addition_filler.c` by modifying the `test_DNA()` 
function and compiling it with gcc:

```
$ gcc genome_addition_filler.c -o gam
$ ./gam
```

Here is the output for triploids:

```
$ ./gam

CREATING TRIPLOID GENOTYPE MODEL
20 POSSIBLE COMBINATIONS
===============
GENOTYPE MODEL:	AAA
GENOTYPE MODEL:	AAC
GENOTYPE MODEL:	AAG
GENOTYPE MODEL:	AAT
GENOTYPE MODEL:	ACC
GENOTYPE MODEL:	ACG
GENOTYPE MODEL:	ACT
GENOTYPE MODEL:	AGG
GENOTYPE MODEL:	AGT
GENOTYPE MODEL:	ATT
GENOTYPE MODEL:	CCC
GENOTYPE MODEL:	CCG
GENOTYPE MODEL:	CCT
GENOTYPE MODEL:	CGG
GENOTYPE MODEL:	CGT
GENOTYPE MODEL:	CTT
GENOTYPE MODEL:	GGG
GENOTYPE MODEL:	GGT
GENOTYPE MODEL:	GTT
GENOTYPE MODEL:	TTT

```

If you wanted to specify a higher ploidy, that option is there, too:

```
$ ./gam 12 | head

CREATING DODECAPLOID GENOTYPE MODEL
455 POSSIBLE COMBINATIONS
===============
GENOTYPE MODEL:	AAAAAAAAAAAA
GENOTYPE MODEL:	AAAAAAAAAAAC
GENOTYPE MODEL:	AAAAAAAAAAAG
GENOTYPE MODEL:	AAAAAAAAAAAT
GENOTYPE MODEL:	AAAAAAAAAACC
GENOTYPE MODEL:	AAAAAAAAAACG

$ ./gam 12 | grep MODEL: | wc
     455    1365   13195
```