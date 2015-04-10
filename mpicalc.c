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
	int* links_fr;		//Node links from me 

	double prob;		//PageRank value
	double n_prob;
} Node;

//Pre-declarations
void print_nodes(int l_bound, int u_bound, Node* nodes);
Node* load_data(const char* filename, int* n_size, int* n_count);
void sort_data(Node *nodes, int n_size, int mpi_size);
void save_data(const char* filename, Node* nodes, int size);
int* get_send_array(int size, Node* nodes, int mpi_rank, int *send_s);
int* get_dest_array(int size, Node* nodes, int mpi_rank, int *send_arr, int send_s);
int* get_recv_array(int size, Node* nodes, int mpi_rank, int *recv_s);
double MPI_Get_Sum(double val, int mpi_rank, int mpi_size);
void MPI_Node_Alltoall(int size, Node* nodebuf, int mpi_rank, int mpi_size);
void MPI_Node_Alltoall2(int size, Node* nodes, int *send_arr, int *dest_arr, int send_s, int *recv_arr, int recv_s);
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

	//get dependency arrays
	/*int send_s, recv_s;					//size of arrays
	int* send_arr = get_send_array(n_size, nodes, mpi_rank, &send_s);
	int* dest_arr = get_dest_array(n_size, nodes, mpi_rank, send_arr, send_s);
	int* recv_arr = get_recv_array(n_size, nodes, mpi_rank, &recv_s);
	printf("[%d] send: %d recv: %d\n", mpi_rank, send_s, recv_s);*/

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
	double thresh = pow(10, -16);
	while (diff > thresh)
	//for (int h = 0; h < 10; ++h)
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
		//get total diff
		diff = MPI_Get_Sum(diff, mpi_rank, mpi_size);
		MPI_Node_Alltoall(n_size, nodes, mpi_rank, mpi_size);
		//MPI_Node_Alltoall2(n_size, nodes, send_arr, dest_arr, send_s, recv_arr, recv_s);
		if(mpi_rank == 0)
			printf("[%d] diff: %e\n", mpi_rank, diff);
			//print_nodes(0, n_size*mpi_size, nodes);
	}
	if(mpi_rank == 0){
		//print_nodes(0, n_size*mpi_size, nodes);
		save_data("data_output", nodes, n_size);
	}
	//Free the nodes
	if(nodes != NULL){
		for (int i = l_bound; i < u_bound; ++i)
		{
			if(nodes[i].id > 0){
				free(nodes[i].links_to);
				nodes[i].links_to = NULL;
				free(nodes[i].links_fr);
				nodes[i].links_fr = NULL;
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
	double tmp;
	for (int i = 0; i < mpi_size; ++i)	//send index
	{
		if(i != mpi_rank){
			//Receive everyone elses pieces
			MPI_Recv(&tmp, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total += tmp;
		}else{
			for (int j = 0; j < mpi_size; ++j)
			{
				if(j != i){	//don't send to yourself
					//Send your piece to everyone else
					MPI_Send(&val, 1, MPI_DOUBLE, j, tag, MPI_COMM_WORLD);
				}
			}
		}
	}
	return total;
}

//Special function for all to all
void MPI_Node_Alltoall(int size, Node* nodebuf, int mpi_rank, int mpi_size){
	int tag = 2;
	for (int i = 0; i < mpi_size; ++i)	//send index
	{
		if(i != mpi_rank){
			//Receive everyone elses pieces
			//printf("%d receiving from %d\n", mpi_rank, i);
			MPI_Recv(nodebuf+i*size, size*sizeof(Node), MPI_BYTE, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}else{
			for (int j = 0; j < mpi_size; ++j)
			{
				if(j != i){	//don't send to yourself
					//printf("%d sending to %d\n", mpi_rank, j);
					//Send your piece to everyone else
					MPI_Send(nodebuf+mpi_rank*size, size*sizeof(Node), MPI_BYTE, j, tag, MPI_COMM_WORLD);
				}
			}
		}
	}
}

//Get arrays telling us what nodes to send
int* get_send_array(int size, Node* nodes, int mpi_rank, int *send_s){
	int send_i = 0;
	int *send_arr = NULL;
	for (int i = mpi_rank*size; i < (mpi_rank+1)*size; ++i)
	{
		if(nodes[i].id > 0){
			Node *n = &nodes[i];
			//Send to those you link
			for (int j = 0; j < n->links_fr_len; ++j)
			{
				int dest = (n->links_fr[j]-1) / size;
				if(dest != mpi_rank){	//don't send to yourself
					send_i++;
					send_arr = (int*) realloc(send_arr, send_i * sizeof(int));
					send_arr[send_i-1] = n->id;
				}
			}
		}
	}
	*send_s = send_i;
	return send_arr;
}

int* get_dest_array(int size, Node* nodes, int mpi_rank, int *send_arr, int send_s){
	int *dest_arr = (int*) malloc(send_s * sizeof(int));
	for (int i = 0; i < send_s; ++i)
	{
		for (int j = 0; j < nodes[send_arr[i]-1].links_fr_len; ++j)
		{
			int dest = (nodes[send_arr[i]-1].links_fr[j]-1) / size;
			if(dest != mpi_rank){	//don't send to yourself
				dest_arr[i] = dest;
			}
		}
	}
	return dest_arr;
}

//Get arrays telling us what nodes to recv
int* get_recv_array(int size, Node* nodes, int mpi_rank, int *recv_s){
	int recv_i = 0;
	int *recv_arr = NULL;
	for (int i = mpi_rank*size; i < (mpi_rank+1)*size; ++i)
	{
		if(nodes[i].id > 0){
			Node *n = &nodes[i];
			//Receive from those who link to you
			for (int j = 0; j < n->links_to_len; ++j)
			{
				int src = (n->links_to[j]-1) / size;
				if(src != mpi_rank){	//don't receive from yourself
					recv_i++;
					recv_arr = (int*) realloc(recv_arr, recv_i * sizeof(int));
					recv_arr[recv_i-1] = n->links_to[j];
				}
			}
		}
	}
	*recv_s = recv_i;
	return recv_arr;
}

//All gather function sending only to nodes which need updating
void MPI_Node_Alltoall2(int size, Node* nodes, int *send_arr, int *dest_arr, int send_s, int *recv_arr, int recv_s){
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	int src;
	MPI_Request *reqs = malloc((send_s+recv_s) * sizeof(MPI_Request));
	for (int i = 0; i < send_s; ++i)
	{
		//printf("[%d] sending %d\n", mpi_rank, send_arr[i]);
		//printf("[%d] dst=%d\n", mpi_rank, dest_arr[i]);
		MPI_Isend(&nodes[send_arr[i]-1].prob, 1, MPI_DOUBLE, dest_arr[i], send_arr[i], MPI_COMM_WORLD, &reqs[i]);
	}
	for (int i = 0; i < recv_s; ++i)
	{
		//printf("[%d] receiving %d\n", mpi_rank, recv_arr[i]);
		src = (recv_arr[i]-1) / size;
		//printf("[%d] src=%d\n", mpi_rank, src);
		MPI_Irecv(&nodes[recv_arr[i]-1].prob, 1, MPI_DOUBLE, src, recv_arr[i], MPI_COMM_WORLD, &reqs[send_s+i]);
	}
	//Wait for termination
	MPI_Status *stats = (MPI_Status*) malloc((send_s+recv_s) * sizeof(MPI_Status));
	MPI_Waitall(send_s + recv_s, reqs, stats);
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
				for (int j = 0; j < mpi_size*size; ++j)
				{
					if(sendbuf[j].id > 0){
						MPI_Send(sendbuf[j].links_to, sendbuf[j].links_to_len, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(sendbuf[j].links_fr, sendbuf[j].links_fr_len, MPI_INT, i, tag, MPI_COMM_WORLD);
					}
				}
			}
		}
	}else{
		//Receive whole array chunk
		MPI_Recv(recvbuf, size*mpi_size*sizeof(Node), MPI_BYTE, send_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//Allocate buffers for struct arrays
		for (int i = 0; i < mpi_size*size; ++i)
		{
			if(recvbuf[i].id > 0){
				recvbuf[i].links_to = (int*) malloc(recvbuf[i].links_to_len*sizeof(int));
				recvbuf[i].links_fr = (int*) malloc(recvbuf[i].links_fr_len*sizeof(int));
			}
		}
		//Receive values in allocated arrays
		for (int i = 0; i < mpi_size*size; ++i)
		{
			if(recvbuf[i].id > 0){
				MPI_Recv(recvbuf[i].links_to, recvbuf[i].links_to_len, MPI_INT, send_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(recvbuf[i].links_fr, recvbuf[i].links_fr_len, MPI_INT, send_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
			printf("(%lf) %d <- ", n->prob, n->id);
			for (int j = 0; j < n->links_to_len; ++j)
			{
				printf("%d ", n->links_to[j]);
			}
			printf(" | -> ");
			for (int j = 0; j < n->links_fr_len; ++j)
			{
				printf("%d ", n->links_fr[j]);
			}
			printf("\n");
		}
	}
	printf("Total: %e\n", prob);
}

Node* load_data(const char* filename, int* n_size, int* n_count){
	//Open file
	FILE* in = fopen(filename, "r");
	if(in == NULL){
		printf("Could not open file to read\n");
		return NULL;
	}

	//Declare Node array
	int size = 1, count = 0;
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
			n->links_fr = NULL;
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
			n->links_fr_len = 1;
			n->links_fr = (int*) malloc(n->links_fr_len * sizeof(int));
			n->links_fr[0] = dst;
			n->links_to_len = 0;
			n->links_to = NULL;
			n->prob = 0.0;
			count++;
		}else{				//Old Node
			n->links_fr_len++;
			n->links_fr = (int*) realloc(n->links_fr, n->links_fr_len * sizeof(int));
			n->links_fr[n->links_fr_len-1] = dst;
		}
	}

	fclose(in);
	*n_size = size;
	*n_count = count;
	return nodes;
}

//Sort node array data to improve performance
void sort_data(Node *nodes, int n_size, int mpi_size){
	return;
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
			fprintf(out, "%d\t%1.10f\n", nodes[i].id, nodes[i].prob);
		}
	}

	fclose(out);
}