#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "adesub.h"
#include <R.h>


void bruvo_dist(int *in, double *out, int *nall, int *perm, int *woo)
{
	int i, j, k, counter=0, n = 2, p = *nall, w = *woo, **genos;
	double dist[p][p], da, res, minn=100;
/*

	This will calculate bruvo's distance between two individuals. 
	All that needs to be done from here is to have it do the pairwise
	calculations. 

    NOTE: The input needs to be divided by the repeat length beforehand for this
    to work. 

	in: a matrix of two individuals
	out: a double value that will be the output from bruvo's distance.
	n: number of individuals(2)
	nall / p: number of alleles
	perm: a vector from the permn function in R
	woo: p * p!
	minn: is a rolling counter of the minimum between allele compairsons.

*/

	/* allocate memory for genotype table.
	tabintalloc(&dist, p, p); */	
	tabintalloc(&genos, n, p);

	/* reconstruct the genotype table */
	for(j=1; j<=p; j++){
		for(i=1; i<=n; i++){
            
            /* Missing data will return with distance of 100 */
            if(in[counter] == 0)
            {
                return;
            }
            else
            {
    			genos[i][j] = in[counter++];
            }        
		}
	}

    /* Construct distance matrix ((THIS WORKS))*/
	for(j=1; j<=p; j++)
	{
		for(i=1; i<=p; i++)
		{
			da = 1- pow(2 ,-abs(genos[1][i]-genos[2][j]));
			dist[i-1][j-1] = da;

            /* DEBUG
			printf("Geno1: %d, Geno2: %d, 1-2^-%d, %11f\n", genos[1][i], 
                genos[2][j], abs(genos[1][i]-genos[2][j]), dist[i-1][j-1]);
			printf("%f\n", da); */
		}
	}

	/* Calculate the smallest s, which is the minimum distance among alleles */

	for(i=0; i < w; i += p)
  	{
		/* DEBUG 
		printf("COUNT: %d\n", i); */
	  	for(j=0; j < p; j++)
	  	{
			if (j == 0)
			{	/* DEBUG 
				printf("point! %d, %d\n", j, *perm); */
				res = dist[j][*perm++];
			}
			else
			{	/* DEBUG 
				printf("point! %d, %d\n", j, *perm); */
				res += dist[j][*perm++];
			}
		}
        /* Checking if the new calculated distance is smaller than the smallest
           distance seen. */
		if ( res < minn )
		{
			minn = res;
        }
		/* DEBUG 
		printf("RES: %11f MIN: %11f\n ", res, minn); */
	}
		
	*out = minn/p;
	/* Free allocated memory 
	freeinttab(dist);*/	
	freeinttab(genos);
}
