/*
This section shows how to use the library to compute
a minimum cut on the following graph:

		        SOURCE
		       /       \
		     1/         \2
		     /      3    \
		   node0 -----> node1
		     |   <-----   |
		     |      4     |
		     \            /
		     5\          /6
		       \        /
		          SINK

///////////////////////////////////////////////////
*/

#include <stdio.h>
#include "graph.h"

int main()
{
	typedef Graph<double,double,double> GraphType;
	GraphType *g = new GraphType(/*estimated # of nodes*/ 2, /*estimated # of edges*/ 1); 

	g -> add_node(); 
	g -> add_node(); 

	g -> add_tweights( 0,   /* capacities */  2, 0);
	g -> add_tweights( 1,   /* capacities */  0, 2 );
	g -> add_edge( 0, 1,    /* capacities */ 10, 10 );

	int flow = g -> maxflow();

	printf("Flow = %d\n", flow);
	printf("Minimum cut:\n");
	printf("%d \n", g->what_segment(0));
	printf("%d \n",GraphType::SOURCE);
	printf("%d \n",GraphType::SINK);
	if (g->what_segment(0) == GraphType::SOURCE)
		printf("node0 is in the SOURCE set\n");
	else
		printf("node0 is in the SINK set\n");
	if (g->what_segment(1) == GraphType::SOURCE)
		printf("node1 is in the SOURCE set\n");
	else
		printf("node1 is in the SINK set\n");

	delete g;

	return 0;
}
