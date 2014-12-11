#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*

This C script will act as supplementary documentation for the process of the
genome addition model calculation of Bruvo's distance in the R package "poppr".
This was written by Zhian N. Kamvar.

Bruvo's distance was proposed by Bruvo in 2004 with three imputation methods for
dealing with organisms of differing ploidies (DOI:
10.1111/j.1365-294X.2004.02209.x).

When comparing two samples at a single locus, the genotype with fewer alleles
must have the extra alleles imputed in order to make a fair comparison. The
genome addition model assumes that the genotype with more alleles arose due to a
genome duplication event. That being so, a fair imputation would be to take the
average distance across all possible combinations of observed alleles at that
locus.

The R package polysat asserts that there are n^k possible combinations where n
is the number of observed alleles and k is the number of "missing" alleles. This
method assumes that phase is important for the calculation. As the calculation
of Bruvo's distance requires a minimum path be drawn through a matrix of
genotypes, phase does not matter as a genotype of 1221 and 1212 will give the
same distance compared to genotype x. Assuming that there are n^k possible
combinations biases the results to heterozygotic genotypes that would be
overrepresented.

Since phase is not important in the calculation, we should give each possible
genotype equal representation. This is done by implementing an algorithm that
iterates through choose(n+k-1, k) possibilities. This is known as the multiset
coefficient (http://en.wikipedia.org/wiki/Multiset#Counting_multisets). This
also conveniently reduced the number of calculations necessary.

This example contains a function for recursion: genome_addition_iteration() and
a function to set up the data and run it: workhorse().

I will leave inline notes within each to explain why things are happening.
*/

void genome_addition_iteration(int* genotype, int zeroes, int inds, 
                               int* zero_ind, int curr_zero, 
                               int* replacement, int curr_ind, int verbose);
void workhorse();

int main(void)
{
	workhorse();
	return 0;
}


/*
Inputs: genotype a one dimensional array. (in practice, it's a distance matrix).
		zeroes the number of missing alleles (k)
		inds   the number of potential replacements (n)
		zero_ind an array giving the indices for each missing allele
		curr_zero an index for zero_ind
		replacement an array giving the indices for each replacement
		curr_ind an index for replacement
		verbose indicator for status
*/
void genome_addition_iteration(int* genotype, int zeroes, int inds, 
                               int* zero_ind, int curr_zero, 
                               int* replacement, int curr_ind, int verbose)
{

	// I call this the DEBUGGERNAUT
	if (verbose)
	{
		printf("STATUS-------------------\n");
		printf("ZEROES & INDS: %d\t%d\n", zeroes, inds);
		printf("CURRENT ZERO:\t\t%d\n", curr_zero);
		printf("CURRENT REPLACEMENT:\t%d\n", curr_ind);
	}
	

	// First step: replace a missing allele with an observed allele.
	genotype[zero_ind[curr_zero]] = genotype[replacement[curr_ind]];

	int i;
	int j;

	// Second step: Loop through all observed alleles BEGINNING WITH THE 
	// OBSERVED ALLELE entering the function.
	for (i = curr_ind; i < inds; i++)
	{
		// If you haven't filled the genotype, call the function again.
		if (curr_zero < zeroes - 1)
		{
			// note the increment in the curr_zero and i instead of curr_ind.
			genome_addition_iteration(genotype, zeroes, inds, zero_ind, 
									  ++curr_zero, replacement, i, verbose);

			// BASE CASE: imputed genotype.
			if (curr_zero == zeroes - 1)
			{
				return;
			}
		}
		else // BASE CASE: imputed genotype.
		{
			// CALCULATION....
			printf("GENOTYPE:\t");
			for (j = 0; j < zeroes+inds; j++)
			{
				printf(" %d", genotype[j]);
			}
			printf("\n");

			// BASE CASE: one zero or no more replacements.
			if (zeroes == 1 || i == inds - 1)
			{
				return;
			}			
		}
		// IMPORTANT: decrementing the zero index allows you to refill the 
		// missing alleles with different replacements.
		curr_zero--;
		if (verbose) printf("\t~~~\n");
	}

	return;
}

void workhorse()
{
	int i;
	int j;
	int nall = 6;
	int zeroes = 3;
	int inds = 3;
	int curr_zero = 0;
	int curr_ind = 0;
	
	int* genotype;
	int* zero_ind;
	int* replacement;
	genotype = (int*) malloc(nall * sizeof(int));
	zero_ind = (int*) malloc(zeroes * sizeof(int));
	replacement = (int*) malloc(inds * sizeof(int));
	zero_ind[0] = 3; zero_ind[1] = 4; zero_ind[2] = 5;
	replacement[0] = 0; replacement[1] = 1; replacement[2] = 2; 
	genotype[0] = 1;
	genotype[1] = 2;
	genotype[2] = 3;
	genotype[3] = 0;
	genotype[4] = 0;
	genotype[5] = 0;
	
	printf("genotype:");
	for (i = 0; i < nall; i++)
	{
		printf(" %d", genotype[i]);
	}
	printf("\n====================\n");
	for (i = 0; i < inds; i++)
	{
		genome_addition_iteration(genotype, zeroes, inds, zero_ind,  
								  curr_zero, replacement, i, 1);		
	}
	free(genotype);
	free(zero_ind);
	free(replacement);
	return;
}




