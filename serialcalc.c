/*
*	Serial version of PageRank algorithm
*	using an Iterative approach
*
*	Denis Lachance, 1244466
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
void print_nodes(int n_size, Node* nodes);
Node* load_data(const char* filename, int* n_size, int* n_count);
void save_data(const char* filename, Node* nodes, int size);

int main(int argc, char const *argv[])
{
	int n_count = 0;	//total number of nodes
	int n_size = 0;		//size of nodes array
	//Load
	Node* nodes = load_data("data_input", &n_size, &n_count);
	printf("%d %d\n", n_size, n_count);
	//Initialize PR values
	double frac = 1.0 / (double) n_count;
	for (int i = 0; i < n_size; ++i)
	{
		if(nodes[i].id > 0){
			nodes[i].prob = frac;
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
		for (int i = 0; i < n_size; ++i)
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
		for (int i = 0; i < n_size; ++i)
		{
			if(nodes[i].id > 0){
				diff += fabs(nodes[i].n_prob - nodes[i].prob);
				nodes[i].prob = nodes[i].n_prob;
			}
		}
		//print_nodes(n_size, nodes);
		printf("diff: %e\n", diff);
	}
	print_nodes(n_size, nodes);
	save_data("data_output", nodes, n_size);
	return 0;
}

void print_nodes(int n_size, Node* nodes){
	//For testing
	Node* n;
	double prob = 0.0;
	printf("----------------\n");
	for (int i = 0; i < n_size; ++i)
	{
		n = &nodes[i];
		if(n->id > 0){
			prob += n->prob;
			/*printf("(%lf) %d -> ", n->prob, n->id);
			for (int j = 0; j < n->links_to_len; ++j)
			{
				printf("%d ", n->links_to[j]);
			}
			printf("\n");*/
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
	int size = 100, count = 0;
	Node* nodes = (Node*) calloc(size, sizeof(Node));

	//Read to end
	int src, dst;
	Node* n;
	while(!feof(in)){
		//Assumes the smallest value is 1
		fscanf(in, "%d\t%d\n", &src, &dst);
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
			fprintf(out, "%d\t%1.10f\n", nodes[i].id, nodes[i].prob);
		}
	}

	fclose(out);
}