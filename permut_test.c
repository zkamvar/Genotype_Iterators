#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	This permutation function will take a number from the user asking for the
	length of the array to be permuted. It prints to screen the permuted array.
	The goal is to modify the main function to take input from R and return
	Bruvo's distance for a single locus. The initial algorithm is from:
	http://www.geeksforgeeks.org/archives/767 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* The global counter variable */
int count;

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
/*		for(j=0; j<=n; j++)
		{
			c[--ind] = *(a+j);

			/* DEBUG:
			printf("%d ", c[ind]);
			if (j == n)
				printf("\n");
			
		}*/
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
	Driver program to test above functions. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int main(void)
{
	int numm, i, j, *nummptr;
	printf("Enter a number for permutation.\n");
	scanf("%d", &numm);
	printf("\nPERMUTING %d...\n", numm);
	nummptr = &numm; //Assign the pointer for the number of alleles.
	int permutations = fact(numm)*numm;
	int *permptr; // Pointer for the number of permutations. 
   	permptr = &permutations;
   	int *perm_array;
   	int allele_array[numm];
   	/* Allocating memory for the permuation array. */
   	perm_array = (int *) malloc(permutations * sizeof(int));
   	int mem = permutations * sizeof(int);
   	/* Populating the array to be permuted. */
	for(i=0; i < numm; i++)
	{
		allele_array[i] = i;
	}
	/* Permuting. */
   	permute(allele_array, 0, numm-1, perm_array);
   	//free(perm_array);
   	/* End of the important stuff.*/   	
	for (i=0; i < permutations; i++)
	{
		printf("%d ", perm_array[i]);
		int n = i;
		n++;
		if(n%numm == 0 && i != 0)
		{
			printf("\n");
		}
	}
	/* DON'T FORGET TO FREE THE MEMORY YOU STOLE! */
	free(perm_array);
   	printf("\nMEMORY USED: %d Kb (about %d Mb)\n", mem/1024, mem/1049000 );	
	printf("\nNumber of Permutations: %d\n\n", permutations); 
	
   	return 0;
}