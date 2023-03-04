#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <queue>
#include <map>
#include <string>
#include <algorithm>


using namespace std;
mt19937 rng(random_device{}());

#include "graph.h"
#include "flow_graph.h"
#include "read_file.cpp"
#include "algorithms/maximum_matching.cpp"

int match_size(const vector<int> &res)
{
    int match = 0;
    for (int i = 0; i < res.size(); i++)
        if (res[i] != -1) match++;
    return match;
}


int main()
{
    graph g = generate_from_file("real_world/bio-CE-GN/bio-CE-GN.txt");
	
    for (int i = 0; i < 1000; i++){
        g.realize(1000);
        vector<int> match = g.maximum_matching();
        cout<<match_size(match)<<" ";
    }
	return 0;
}
