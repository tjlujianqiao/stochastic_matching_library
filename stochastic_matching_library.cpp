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
#include "algorithms/algorithms.h"


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
    // preprocess
    map<pair<int, int>, double> typeProb;
	
    for (int i = 0; i < 1000; i++){
        g.realize(1000);
        vector<int> match = g.maximum_matching();
        vector<int> SWR= g.sampling_without_replacement(typeProb);
        vector<int> ranking= g.ranking();
        vector<int> balance_SWR= g.balance_SWR();
        vector<int> balance_OCS= g.balance_OCS();
        cout<<match_size(match)<<" ";
    }
	return 0;
}
