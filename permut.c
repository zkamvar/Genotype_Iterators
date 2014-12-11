# include <stdio.h>
 
/* Function to swap values at two pointers */
void swap (int *x, int *y)
{
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
}
  
/* Function to print permutations of string
   This function takes three parameters:
   1. String
   2. Starting index of the string
   3. Ending index of the string. */
void permute(int *a, int i, int n) 
{
    int j;
    if (i == n)
    {
		for(j=0; j<=n; j++)
		{
			printf("%d ", *(a+j));
			if (j == n)
				printf("\n");
		}
	}
    else
    {
        for (j = i; j <= n; j++)
        {
          	swap((a+i), (a+j));
 			permute(a, i+1, n);
			swap((a+i), (a+j)); //backtrack
       	}
   	}
} 
 
/* Driver program to test above functions */
int main()
{
   	int a[] = {1,2,3,4,5};  
   	permute(a, 0, 4);
   	/*getchar();*/
   	return 0;
}
