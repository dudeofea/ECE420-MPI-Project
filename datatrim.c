/*
 * The script to fetch the subset of the original data
 *
 * There is one command line input which indicate the desired upper bound for the node index.
 * The output will be stored in the file "data_input" with all the links of which the nodes have the index not greater than the value you passed.
 * The total output nodes and edges will be displayed in the screen.
 *
 *
 * Note that the original directed graph is quiet sparse and if you give a number too small, there will be no existing edges and nodes. If the output edges is 0, try larger input value.
 * Also note that due to the sparsity of the original data, the engaged nodes index might not be continues. However, they are less or equal than the passed value.  
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char* argv[])
{
	int N, src, dst;
	FILE * fp_ori, * fp_dest;
	int* flag;
	int Ecount=0;
	int i,j;
	char* tempstore;
	if (argc==1)
	{
		printf("Please assign the number of nodes.\n");
		return 1;
	}
	N=strtol(argv[1], NULL, 10);

	if ((fp_ori=fopen("web-Stanford.txt","r"))==NULL)
	{
		printf("Fail to open a file. \n");
		return 2;
	}
	
	if ((fp_dest=fopen("data_input","w"))==NULL)
	{
		printf("Fail to open a file. \n");
		return 3;
	}
	
	flag=malloc(N*sizeof(int));
	tempstore=malloc(100*sizeof(char));
	for (i=0; i<N; ++i)
		flag[i]=0;	
	
	for (j=0; j<4; ++j)
	{
		fgets(tempstore, 100, fp_ori);
		//puts(tempstore);
	}
	free(tempstore);

	while(!feof(fp_ori))
	{
		fscanf(fp_ori, "%d\t%d\n", &src, &dst);
		if (src<N && dst <N)
		{
			fprintf(fp_dest, "%d\t%d\n", src, dst);
			++Ecount;
			flag[src]=1;
			flag[dst]=1;

		}
	}
	printf("There are %d edges in the sub dataset. \n", Ecount);
	Ecount=0;
	for (i=0; i<N; ++i)
		Ecount+=flag[i];
	printf("There are %d nodes in the sub dataset. \n", Ecount);
	free(flag);
	fclose(fp_ori);
	fclose(fp_dest);
	return 0;
}
