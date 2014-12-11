#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

/* The global counter variable */
int count;
void bruvo_dist(int *in, double *out, int *nall, int *perm, int *woo);
void permute(int *a, int i, int n, int *c);
double mindist(int perms, int alleles, int *perm, double *dist);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	The permutation function algorithm is modified from:
	http://www.geeksforgeeks.org/archives/767 
	
	This is meant to run on 4 tetraploid individuals at one locus with the
	following matrix of alleles:

	     A0   A1   A2   A3

	G0   69   65   70   65
	
	G1   68   73   69   69
	
	G2   63   67   74   66
	
	G3   67   74   63   71
	
	To change how the program runs, change the value of numm to change the
	ploidy, and make sure to change testmat as well (though, if it's a ploidy
	less than 4, it can stay the same).
	
	To change the number of individuals, you have to change the following lines:
	line 42
	line 45
	...or you could do a specific grep with \D4 and be careful.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int main(void)
{
	int numm, i, j, *nummptr, ind;
	
	numm = 5;//5
	ind = 4;
	/* Start of the important stuff. */
	//int testmat[4][4] = {{69,70,71,72},{71,69,72,70},{72,71,70,69},{70,72,69,71}};
	int testmat[4][5] = {69,65,70,65,68,
						68,73,69,69,71,
   						63,67,74,66,69,
   						67,74,63,71,69};
   	double *out;
	nummptr = &numm; //Assign the pointer for the number of alleles.
	int allele_array[numm];
	int permutations = fact(numm)*numm;
   	int *perm_array; // Array to store the permutations.
   	int *permptr; // Pointer for the number of permutations. 
   	permptr = &permutations;
   	out = (double *) malloc((ind*(ind-1)/2) * sizeof(double));
	int bruvomat[1+numm*2]; // Temp matrix for comparisons
	int *bruvoptr;
   	/* Allocating memory for the permuation array. */
   	perm_array = (int *) malloc(permutations * sizeof(int));
   	int mem = permutations * sizeof(int);
   	//printf("MEMORY REQUIRED: %d Kb (about %d Mb)\n", mem/1024, mem/1049000 );
   	/* Populating the array to be permuted. */
	for(i=0; i < numm; i++)
	{
		allele_array[i] = i;
	}
	/* Permuting. */
   	permute(allele_array, 0, numm-1, perm_array);
   	//free(perm_array);
   	
   	/* End of the important stuff.   	
	for (i=0; i < permutations; i++)
	{
		printf("%d\t", perm_array[i]);
		int n = i;
		n++;
		if(n%numm == 0 && i != 0)
		{
			printf("\n");
		}
	}
	*/
	printf("\n     GENOTYPE MATRIX:\n\n");
	printf("     A0   A1   A2   A3   A4\n\n"); //A4
	//printf("     A0   A1\n\n");
	for(i=0; i < ind; i++) // pairwise comparisons.
	{
		printf("G%d", i);
		for(j=0; j < numm; j++)
		{
			printf("   %d", testmat[i][j]);
		}
		printf("\n\n");
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		This is where we begin doing the pairwise comparisons over the "alleles"
		that are represented by testmat, or the matrix of alleles that will be
		analyzed. Since this is a crude version dealing with diploids, I've 
		hardcoded the alleles. This will change, however. It shouldn't be too
		hard to code another for loop in there. 
		NOTES:
		count is set to zero here so that bruvo_dist can add sequentially to 
		the array "out".
		bruvoptr points to bruvomat after all of the values have been acquired
		for comparison. 
		nummptr simply points to the number of alleles.
		perm_array is the permutation array that was generated earlier. 
		permptr is the pointer to the number of permutations. 
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	count = 0;
	for(i=0; i < ind; i++) // pairwise comparisons.
	{
		int z;
		for (z=0; z < numm; z++)
		{
			bruvomat[z] = testmat[i][z];
			//printf("testmat[%d][%d]: %d\n", i, z, bruvomat[z]);
		}
		for(j = i+1; j > i && j < ind; j++)
		{	
			for (z=numm; z < numm*2; z++)
			{
				bruvomat[z] = testmat[j][z-numm];
				//printf("testmat[%d][%d]: %d\n", j, z-numm, bruvomat[z]);
			}
			bruvoptr = bruvomat;
			bruvo_dist(bruvoptr, out, nummptr, perm_array, permptr);
		}
	}
	/* DON'T FORGET TO FREE THE MEMORY YOU STOLE! */
	free(perm_array);
	count=0;
	printf("\nBRUVO'S DISTANCES FOR EACH PAIRWISE COMPARISON:\n\n");
	for(i=0; i < ind; i++)
	{
		for(j = i+1; j > i && j < ind; j++)
		{
			printf("G%d to G%d: %9f\n", i, j, *(out+count));
			count++;
		}
	}
	printf("\nNumber of Permutations: %d\n\n", permutations); 
   	return 0;
}

/* Function to swap values at two pointers */
void swap (int *x, int *y)
{
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
}
  
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	Function to print permutations of string
   	This function takes four parameters:
   	1. String
   	2. Starting index of the string
   	3. Ending index of the string. 
   	4. pointer to array of size n*n! 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void permute(int *a, int i, int n, int *c) 
{
	int j;
    if (i == n)
    {
    	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    		'a' will be the array containing the numeric sequence to be
    		permuted. It will be reshuffled into a new pattern each
    		time it reaches this control structure. To place the value
    		into the array 'c', the pointer for a needs to be incremented
    		over all its elements.
    	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    	count += n+1;
    	//DEBUG: printf("\t\tCOUNT: %d\n", count);
    	int ind = count;
		for(j = n; j >= 0; j--)
		{
			c[--ind] = *(a+j);
		}
	}
    else
    {
        for (j = i; j <= n; j++)
        {
          	swap((a+i), (a+j));
 			permute(a, i+1, n, c);
			swap((a+i), (a+j)); //backtrack
       	}
   	}
} 

/* A factorial function for calculating permutations */
int fact(int x)
{
    int f=1;
    int u;
    for (u=x; u>1; u--)
    {
        f*=u;
    }
    return f;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void bruvo_dist(int *in, double *out, int *nall, int *perm, int *woo)
{
	int i, j, k, counter=0, n = 2, p = *nall, w = *woo, genos[2][p];
	double dist[p][p], da, res, minn=100, *distp;	
	/* reconstruct the genotype table */
	for(i=0; i < n; i++){
		for(j=0; j < p; j++){
            
            /* Missing data will return with distance of 100 */
            if(in[counter] == 0)
            {
                return;
            }
            else
            {
            	//printf("genos[%d][%d]: %d\n", i, j, in[counter]);
    			genos[i][j] = in[counter++];
            }        
		}
	}

    /* Construct distance matrix ((THIS WORKS)) */
	for(j=0; j < p; j++)
	{
		for(i=0; i < p; i++)
		{
			da = 1- pow(2 ,-abs(genos[0][i]-genos[1][j]));
			dist[i][j] = da;

			printf("Geno1: %d Geno2: %d | 1-2^-%d =%9f\n", genos[0][i], 
                genos[1][j], abs(genos[0][i]-genos[1][j]), dist[i][j]);
		}
	}

	/* Calculate the smallest s, which is the minimum distance among alleles */
	printf("\n    DISTANCE MATRIX:\n\n");
	printf("       %d        %d        %d        %d        %d\n\n", genos[1][0], genos[1][1], genos[1][2], genos[1][3], genos[1][4]);
	//printf("       %d        %d\n\n", genos[1][0], genos[1][1]);

	for(i=0; i < p; i++)
	{
		printf("%d ", genos[0][i]);
		for(j=0; j < p; j++)
	  	{
	  		printf("%9f ", dist[i][j]);
	  	}
	  	printf("\n\n");
	}
	distp = (double *) &dist;
	out[count++] = mindist(w, p , perm, distp)/p;
/*	for(i=0; i < w; i += p)
  	{
		//printf("COUNT: %d\n", i); 
	  	for(j=0; j < p; j++)
	  	{
			if (j == 0)
			{	 
				printf("dist[%d][%d] : %9f\n", *perm, j, dist[*perm][j]); 
				res = dist[*perm++][j];
				//printf("dist[%d][%d] : %9f\n", j, *perm, dist[j][*perm]); 
				//res = dist[j][*perm++];
			}
			else
			{	
				printf("dist[%d][%d] : %9f\n", *perm, j, dist[*perm][j]); 
				res += dist[*perm++][j];			
				//printf("dist[%d][%d] : %9f\n", j, *perm, dist[j][*perm]); 
				//res += dist[j][*perm++];
			}
		}
        /* Checking if the new calculated distance is smaller than the smallest
           distance seen. 
		if ( res < minn )
		{
			minn = res;
        }
		printf("AVG: %9f MIN: %9f\n", res/p, minn/p);
	}
	out[count++] = minn/p;*/
}

double mindist(int perms, int alleles, int *perm, double *dist)
{
	int i, j, w = perms, p = alleles, counter = 0;
	double res = 0, minn = 100;
	for(i = 0; i < w; i += p)
	{
		for(j = 0; j < p; j++)
		{
			if (j == 0)
			{
				printf("dist[%d][%d] : %9f\n", *(perm + counter), j, dist[*(perm + counter) + p*j]); 
				res = dist[*(perm + counter++) + p*j];
				if(res > minn)
				{
					//printf("\n===\nCurrent count: %d, i: %d, Modified: %d, Total: %d\n===\n", counter, i, i + (w/p), w);
					printf("\t======\n\tBOUND!\n\t======\n");
					j = p;
					counter = i + w/p;
					i = counter;
				}
			}
			else
			{
				printf("dist[%d][%d] : %9f\n", *(perm + counter), j, dist[*(perm + counter) + p*j]); 
				res += dist[*(perm + counter++) + p*j];
				if(j < p-1 && res > minn)
				{
					//printf("DERP!!! Counter: %d, Counter Update: %d, j: %d\n", counter, counter + (p-j-1), j);
					printf("\t======\n\t hop!\n\t======\n");					
					counter += (p-j-1);
					j = p;
				}		
			}
		}
		/*	Checking if the new calculated distance is smaller than the smallest
		 distance seen. */
		if ( res < minn )
		{
			minn = res;
		}
		printf("AVG: %9f MIN: %9f\n", res/p, minn/p);
	}
	return minn;
}

