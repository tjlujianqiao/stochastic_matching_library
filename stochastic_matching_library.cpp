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


void compute_mean_std(const vector<double> &res, double &mean, double &std){
    mean = 0;
    for (auto item: res) mean += item;
    mean /= res.size();
    
    std = 0;
    for (auto item: res) std += (item - mean) * (item - mean);
    
    std = sqrt(std / (res.size() - 1));
}


int main()
{
    graph g = generate_from_file("real_world/bio-CE-GN/bio-CE-GN.txt");
    // preprocess
    
    int numSample = 1000;
    int realSize = g.online_size();
    
    
    map<pair<int, int>, double> typeProb = g.optimal_matching_prob(numSample, realSize);
    
    vector<double> offMass = g.poisson_offline_mass(typeProb);
	
    vector<double> OPT, SWR, ranking, balanceSWR, balanceOCS, poissonOCS, topHalf, minDegree;
    
    for (int i = 0; i < numSample; i++){
        if(i%10==0)cout<<"WORKING on "<<i<<endl;
        g.realize(realSize);
        
        OPT.push_back(match_size(g.maximum_matching()));
        
        SWR.push_back(match_size(g.sampling_without_replacement(typeProb)));
        
        ranking.push_back(match_size(g.ranking()));
        
        balanceSWR.push_back(match_size(g.balance_SWR()));
        
        balanceOCS.push_back(match_size(g.balance_OCS()));
        
        poissonOCS.push_back(match_size(g.poisson_OCS(offMass, typeProb)));
        
        topHalf.push_back(match_size(g.top_half_sampling(typeProb)));
        
        minDegree.push_back(match_size(g.min_degree()));
    }
    
    
    
    double mean, std;
    compute_mean_std(OPT, mean, std);
    cout<<"OPT:        "<<mean<<" "<<std<<endl;
    
    double x, y;
    compute_mean_std(SWR, x, y);
    cout<<"SWR:        "<<x / mean<<" "<<y / mean<<endl;
    
    compute_mean_std(ranking, x, y);
    cout<<"Ranking:    "<<x / mean<<" "<<y / mean<<endl;
    
    compute_mean_std(balanceSWR, x, y);
    cout<<"BalanceSWR: "<<x / mean<<" "<<y / mean<<endl;
    
    compute_mean_std(balanceOCS, x, y);
    cout<<"BalanceOCS: "<<x / mean<<" "<<y / mean<<endl;
    
    compute_mean_std(poissonOCS, x, y);
    cout<<"PoissonOCS: "<<x / mean<<" "<<y / mean<<endl;
    
    compute_mean_std(topHalf, x, y);
    cout<<"TopHalf:    "<<x / mean<<" "<<y / mean<<endl;
    
    compute_mean_std(minDegree, x, y);
    cout<<"MinDegree:  "<<x / mean<<" "<<y / mean<<endl;
    
	return 0;
}
