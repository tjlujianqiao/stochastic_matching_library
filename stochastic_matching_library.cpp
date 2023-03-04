#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <queue>
#include <map>
#include <string>
#include <algorithm>
#include <cstdio>



using namespace std;
mt19937 rng(random_device{}());

#include "graph.h"
#include "flow_graph.h"
#include "decomposite_graph.h"
#include "read_file.cpp"
#include "algorithms/algorithms.h"


int match_size(const vector<int> &res)
{
    int match = 0;
    for (int i = 0; i < res.size(); i++)
        if (res[i] != -1) match++;
    return match;
}


pair<double,double> compute_mean_std(const vector<double> &res){
    double mean = 0;
    for (auto item: res) mean += item;
    mean /= res.size();
    
    double std = 0;
    for (auto item: res) std += (item - mean) * (item - mean);
    
    std = sqrt(std / (res.size() - 1));
    
    return make_pair(mean, std);
}


struct resAlg{
    string name;
    vector<double> resRun;
    vector<double> resSample;
    vector<pair<double, double>> resGraph;
    
    
    resAlg(string s){
        name = s;
        resRun = {};
        resSample = {};
        resGraph = {};
    }
    
    void add_run(double item){
        resRun.push_back(item);
    }
    
    void summary_run(double base = 1){
        resSample.push_back(compute_mean_std(resRun).first / base);
        resRun.clear();
    }
    
    void summary_sample(){
        resGraph.push_back(compute_mean_std(resSample));
        resSample.clear();
    }
    
    
}OPT("OPT"),
 SWR("SWR"),
 ranking("Ranking"),
 minDegree("MinDegree"),
 feldmanMMM("FeldManEtAl"),
 manshadiGS("ManshadiEtAl"),
 jailletLu("JailletLu"),
 topHalf("TopHalfSampling"),
 poissonOCS("PoissonOCS"),
 balanceSWR("BalanceSWR"),
 balanceOCS("BalanceOCS");

int main()
{
    resAlg* resPointer[] = {&OPT, &SWR, &ranking, &balanceSWR, &balanceOCS, &poissonOCS, &topHalf, &minDegree, &jailletLu, &manshadiGS, &feldmanMMM};
    
    
    int numGraph = 1;
    int numSample = 1000;
    
    cout << "Rep";
    for(int i = 1; i <= numGraph; i++)
    {
        cout << " " << i;
        graph g = generate_from_file("real_world/bio-CE-GN/bio-CE-GN.txt");
        
        //Preprocessing
        int realSize = g.online_size();
        map<pair<int, int>, double> typeProb = g.optimal_matching_prob(numSample, realSize);
        vector<double> offMass = g.poisson_offline_mass(typeProb);
        vector<vector<int>> jlList = g.jaillet_lu_list();
        
        vector<int> blueF, redF;
        tie(blueF, redF) = g.feldman_et_al_color();
        
        
        for (int j = 0; j < numSample; j++)
        {
            g.realize(realSize);
            
            OPT.add_run(match_size(g.maximum_matching()));
            SWR.add_run(match_size(g.sampling_without_replacement(typeProb)));
            ranking.add_run(match_size(g.ranking()));
            balanceSWR.add_run(match_size(g.balance_swr()));
            balanceOCS.add_run(match_size(g.balance_ocs()));
            poissonOCS.add_run(match_size(g.poisson_ocs(offMass, typeProb)));
            topHalf.add_run(match_size(g.top_half_sampling(typeProb)));
            minDegree.add_run(match_size(g.min_degree()));
            jailletLu.add_run(match_size(g.jaillet_lu(jlList)));
            manshadiGS.add_run(match_size(g.manshadi_et_al(typeProb)));
            feldmanMMM.add_run(match_size(g.feldman_et_al(blueF, redF)));
        }
        
        double opt = compute_mean_std(OPT.resRun).first;
        for (auto j : resPointer) (*j).summary_run(opt);
    }
    
    //If you want to summarize the data on one graph, uncomment the following line
    //for (auto j : resPointer) (*j).summary_sample();
    
    cout << endl;
    for (auto i : resPointer) {
        
        cout << (*i).name << ":";
        for(auto j : (*i).resSample) 
            cout<<" "<<j;
        
        cout<<endl;
        
    }
    
    cout<<endl;
    
	return 0;
}
