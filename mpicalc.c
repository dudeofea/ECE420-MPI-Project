/*
*	MPI version of PageRank algorithm
*	using an Iterative approach
*
*	Denis Lachance, 1244466
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

//Node structure
typedef struct{
	int id;				//id of the Node

	int links_to_len;	//How many links go to me
	int* links_to;		//Node links to me
	int links_fr_len;	//How many links come from me 

	double prob;		//PageRank value
	double n_prob;
} Node;

//Pre-declarations
void print_nodes(int l_bound, int u_bound, Node* nodes);
Node* load_data(const char* filename, int* n_size, int* n_count);
void save_data(const char* filename, Node* nodes, int size);
double MPI_Get_Sum(double val, int mpi_rank, int mpi_size);
void MPI_Node_Allgather(int size, Node* nodebuf, int mpi_rank, int mpi_size);
void MPI_Node_Bcast(int size, Node* sendbuf, Node* recvbuf, int send_rank, int mpi_rank, int mpi_size);

int main(int argc, char *argv[])
{
	//Init MPI
	int mpi_rank, mpi_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);

	printf("I am process %d of %d\n", mpi_rank, mpi_size);
	
	int n_count = 0;	//total number of nodes
	int n_size = 0;		//size of nodes array
	Node* nodes;		//individual array of nodes
	Node* all_nodes = NULL;	//Total array of all nodes
	//Load
	//TODO: figure out a way to send links_to along with nodes
	//to other processes
	if(mpi_rank == 0){
		all_nodes = load_data("data_input", &n_size, &n_count);
		printf("Input size is %d with %d processes\n", n_size, mpi_size);
		int ind_size = n_size / mpi_size;
		if(ind_size*mpi_size != n_size){
			printf("Please select a proper number of processes for an input size of %d\n", n_size);
			MPI_Finalize();
			return 0;
		}
		n_size = ind_size;
	}
	MPI_Bcast(&n_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	printf("%d %d\n", n_size, n_count);
	if(mpi_rank == 0)
		nodes = all_nodes;	//no need to keep the original data intact
	else
		nodes = (Node*) calloc(n_size*mpi_size, sizeof(Node));

	MPI_Node_Bcast(n_size, all_nodes, nodes, 0, mpi_rank, mpi_size);
	
	//Initialize PR values
	double frac = 1.0 / (double) n_count;
	int ind_count = 0;		//number of non-zero nodes in individual array
	int l_bound = mpi_rank*n_size;		//lower bound of nodes
	int u_bound = (mpi_rank+1)*n_size;	//upper bounds of nodes
	for (int i = 0; i < n_size*mpi_size; ++i)
	{
		if(nodes[i].id > 0){
			nodes[i].prob = frac;
			if(i >= l_bound && i < u_bound){
				ind_count++;
			}
		}
	}
	//Continue until threshold reached
	double diff = 1.0;
	//while (diff > 0.1)
	for (int h = 0; h < 10; ++h)
	{
		Node* n;
		diff = 0.0;
		//Iterate over all nodes
		for (int i = l_bound; i < u_bound; ++i)
		{
			if (nodes[i].id > 0)
			{
				n = &nodes[i];
				double new = 0.0;
				int l;
				for (int i = 0; i < n->links_to_len; ++i)
				{
					l = nodes[n->links_to[i]-1].links_fr_len;
					//If node has no other links
					if(l == 0){
						l = n_count - 1;
					}
					new += nodes[n->links_to[i]-1].prob / l;
				}
				n->n_prob = 0.15 * frac + 0.85 * new;
			}
		}
		//Move n_prob to prob
		for (int i = l_bound; i < u_bound; ++i)
		{
			if(nodes[i].id > 0){
				diff += fabs(nodes[i].n_prob - nodes[i].prob);
				nodes[i].prob = nodes[i].n_prob;
			}
		}
		//normalize diff
		diff /= ind_count;
		//get total diff
		diff = MPI_Get_Sum(diff, mpi_rank, mpi_size);
		MPI_Node_Allgather(n_size, nodes, mpi_rank, mpi_size);
		if(mpi_rank == 0)
			print_nodes(0, n_size*mpi_size, nodes);
	}
	//save_data("data_output", nodes, n_size);
	//Free the nodes
	if(nodes != NULL){
		for (int i = l_bound; i < u_bound; ++i)
		{
			if(nodes[i].id > 0){
				free(nodes[i].links_to);
				nodes[i].links_to = NULL;
			}
		}
	}
	free(nodes);
	MPI_Finalize();
	return 0;
}

double MPI_Get_Sum(double val, int mpi_rank, int mpi_size){
	int tag = 3;
	double total = val;
	//Send your piece to everyone else
	for (int i = 0; i < mpi_size; ++i)
	{
		if(i != mpi_rank){
			MPI_Send(&val, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
		}
	}
	//Receive everyone elses pieces
	for (int i = 0; i < mpi_size; ++i)
	{
		if(i != mpi_rank){
			MPI_Recv(&val, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total += val;
		}
	}
	return total;
}

//Special function to all gather nodes
void MPI_Node_Allgather(int size, Node* nodebuf, int mpi_rank, int mpi_size){
	int tag = 2;
	//Send your piece to everyone else
	for (int i = 0; i < mpi_size; ++i)
	{
		if(i != mpi_rank){
			MPI_Send(nodebuf+mpi_rank*size, size*sizeof(Node), MPI_BYTE, i, tag, MPI_COMM_WORLD);
		}
	}
	//Receive everyone elses pieces
	for (int i = 0; i < mpi_size; ++i)
	{
		if(i != mpi_rank){
			MPI_Recv(nodebuf+i*size, size*sizeof(Node), MPI_BYTE, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
}

//Special function to scatter Node arrays
void MPI_Node_Bcast(int size, Node* sendbuf, Node* recvbuf, int send_rank, int mpi_rank, int mpi_size){
	int tag = 1;
	if(mpi_rank == send_rank){
		//Send whole array chunk
		for (int i = 0; i < mpi_size; ++i)
		{
			if(i != send_rank){
				MPI_Send(sendbuf, size*mpi_size*sizeof(Node), MPI_BYTE, i, tag, MPI_COMM_WORLD);
			}
		}
		//Send data in struct arrays
		for (int i = 0; i < mpi_size; ++i)
		{
			if(i != send_rank){
				for (int j = i*size; j < (i+1)*size; ++j)
				{
					if(sendbuf[j].id > 0){
						MPI_Send(sendbuf[j].links_to, sendbuf[j].links_to_len, MPI_INT, i, tag, MPI_COMM_WORLD);
					}
				}
			}
		}
	}else{
		//Receive whole array chunk
		MPI_Recv(recvbuf, size*mpi_size*sizeof(Node), MPI_BYTE, send_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//Allocate buffers for struct arrays
		for (int i = mpi_rank*size; i < (mpi_rank+1)*size; ++i)
		{
			if(recvbuf[i].id > 0){
				recvbuf[i].links_to = (int*) malloc(recvbuf[i].links_to_len*sizeof(int));
			}
		}
		//Receive values in allocated arrays
		for (int i = mpi_rank*size; i < (mpi_rank+1)*size; ++i)
		{
			if(recvbuf[i].id > 0){
				MPI_Recv(recvbuf[i].links_to, recvbuf[i].links_to_len, MPI_INT, send_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}
}

void print_nodes(int l_bound, int u_bound, Node* nodes){
	//For testing
	Node* n;
	double prob = 0.0;
	printf("----------------\n");
	for (int i = l_bound; i < u_bound; ++i)
	{
		n = &nodes[i];
		if(n->id > 0){
			prob += n->prob;
			printf("(%lf) %d -> ", n->prob, n->id);
			/*for (int j = 0; j < n->links_to_len; ++j)
			{
				printf("%d ", n->links_to[j]);
			}*/
			printf("\n");
		}
	}
	printf("Total: %lf\n", prob);
}

Node* load_data(const char* filename, int* n_size, int* n_count){
	//Open file
	FILE* in = fopen(filename, "r");
	if(in == NULL){
		printf("Could not open file to read\n");
		return NULL;
	}

	//Declare Node array
	int size = 100, count = 0;
	Node* nodes = (Node*) calloc(size, sizeof(Node));

	//Read to end
	int src, dst;
	Node* n;
	while(!feof(in)){
		//Assumes the smallest value is 1
		if(fscanf(in, "%d\t%d\n", &src, &dst) < 2){
			//error in the input file
			return NULL;
		}
		//Check array size
		if(src > size || dst > size){
			while(src > size || dst > size)
				size *= 2;
			nodes = (Node*) realloc(nodes, size * sizeof(Node));
			if(nodes == NULL){
				printf("realloc() failed!\n");
			}
		}
		//Load destination node
		n = &nodes[dst-1];
		if (n->id != dst){	//New Node
			n->id = dst;
			n->links_to_len = 1;
			n->links_to = (int*) malloc(n->links_to_len * sizeof(int));
			n->links_to[0] = src;
			n->links_fr_len = 0;
			n->prob = 0.0;
			count++;
		}else{				//Old Node
			n->links_to_len++;
			n->links_to = (int*) realloc(n->links_to, n->links_to_len * sizeof(int));
			n->links_to[n->links_to_len-1] = src;
		}

		//Load source node
		n = &nodes[src-1];
		if(n->id != src){	//New Node
			n->id = src;
			n->links_to_len = 0;
			n->links_to = NULL;
			n->links_fr_len = 1;
			n->prob = 0.0;
			count++;
		}else{				//Old Node
			n->links_fr_len++;
		}
	}

	fclose(in);
	*n_size = size;
	*n_count = count;
	return nodes;
}

void save_data(const char* filename, Node* nodes, int size){
	//Open file
	FILE* out = fopen(filename, "w");
	if(out == NULL){
		printf("Could not open file to write\n");
		return;
	}

	//Write data
	for (int i = 0; i < size; ++i)
	{
		if(nodes[i].id > 0){
			fprintf(out, "%d\t%lf\n", nodes[i].id, nodes[i].prob);
		}
	}

	fclose(out);
}