// Main program 

#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <queue>
#include <map>
#include <string>
#include <algorithm>
#include <set>
#include "glpk.h" // For Brubach et al. (2016)

using namespace std;
mt19937 rng(random_device{}());

#include "graph.h"
#include "flow_graph.h"
#include "cycle_break_graph.h"
#include "decomposite_graph.h"
#include "natural_lp.h"
#include "read_file.cpp"
#include "algorithms/algorithms.h"

// Return size of matching
// -1 represents for not matched
int match_size(const vector<int> &res)
{
    int match = 0;
    for (int i = 0; i < (int)res.size(); i++)
        if (res[i] != -1)
            match++;
    return match;
}

//Compute mean and standard deviation
pair<double, double> compute_mean_std(const vector<double> &res)
{
    double mean = 0;
    for (auto item : res)
        mean += item;
    mean /= res.size();

    double std = 0;
    for (auto item : res)
        std += (item - mean) * (item - mean);

    std = sqrt(std / (res.size() - 1));

    return make_pair(mean, std);
}


//Store the results
struct resAlg
{
    // Algorithm name
    string name;
    
    // Results of runs on one type graph (sampling online vertices)
    vector<double> resRun;
    
    // Results of runs on one dataset (sampling type graphs)
    vector<pair<double, double>> resSample;
    
    // Results of runs on all datasets
    vector<pair<double, double>> resDataset;

    resAlg(string s)
    {
        name = s;
        resRun = {};
        resSample = {};
        resDataset = {};
    }
    
    // Add result of one run
    void add_run(double item)
    {
        resRun.push_back(item);
    }
    
    // Summarize runs on one type graph (sampling online vertices)
    void summary_run(double base = 1)
    {
        auto mean_std = compute_mean_std(resRun);
        resSample.push_back(make_pair(mean_std.first / base, mean_std.second / base));
        resRun.clear();
    }

    // Summarize runs on one dataset (sampling type graphs)
    void summary_sample()
    {
        if (resSample.size() == 1)
        {
            resDataset.push_back(resSample[0]);
            resSample.clear();
            return;
        }
        
        vector<double> res;
        for (auto item : resSample)
            res.push_back(item.first);
        
        resDataset.push_back(compute_mean_std(res));
        resSample.clear();
    }

}   OPT("OPT"),
    SWR("SWR"),
    SWRType("SWRType"),
    ranking("Ranking"),
    balanceSWR("BalanceSWR"),
    balanceOCS("BalanceOCS"),
    minDegree("MinDegree"),
    feldmanMMM("FeldManEtAl"),
    BahmaniKapralov("BahmaniKapralov"),
    manshadiGS("ManshadiEtAl"),
    jailletLu("JailletLu"),
    brubachSSX("BrubachEtAl"),
    topHalf("TopHalfSampling"),
    poissonOCS("PoissonOCS");   //All algorithms


//Presentation order in output
vector<resAlg*> resPointer = {
    &OPT,
    &SWR,
    &SWRType,
    &poissonOCS,
    &minDegree,
    &balanceSWR,
    &balanceOCS,
    &ranking,
    &topHalf,
    &brubachSSX,
    &jailletLu,
    &manshadiGS,
    &BahmaniKapralov,
    &feldmanMMM
};

vector <string> datasetName;

// Save results to directory
void save_results_to_files(string directory)
{
    ofstream fileResMean, fileResStd;
    cout << "Save results into file " << directory + "/" + "(resMean.txt,resStd.txt)" << endl;
    fileResMean.open(directory + "/" + "resMean.txt");
    fileResStd.open(directory + "/" + "resStd.txt");
    
    fileResMean << "Algorithm";
    for (auto item : datasetName)
        fileResMean << " " << item;
    fileResMean << endl;
    
    fileResStd << "Algorithm";
    for (auto item : datasetName)
        fileResStd << " " << item;
    fileResStd << endl;
    
    for (auto i : resPointer)
        if (i != &OPT)
        {
            fileResMean << (*i).name << " ";
            fileResStd << (*i).name << " ";
            
            for (auto j : (*i). resDataset)
            {
                fileResMean << j.first  << " ";
                fileResStd << j.second << " ";
            }
            
            fileResMean << endl;
            fileResStd << endl;
            
        }
    fileResMean.close();
    fileResStd.close();

    cout << "Output Results Done!" << endl;
}

// Run experiments on graphs generated from file
void work_from_file(string name)
{
    cout << "Working on file " << name << endl;

    int numGraph = 1;
    int numSample = 1000;

    cout << "Rep";
    for (int i = 1; i <= numGraph; i++)
    {
        cout << " " << i;
        graph g = generate_from_file(name, true, 0);
        
        // Preprocessing
        int realSize = g.online_size();
        
        map<pair<int, int>, double> typeProb  = g.optimal_matching_prob(numSample, realSize);

        // Extremely slow to compute natural LP
        // natural_lp lp(g.get_adj(),g.online_size());
        // map<pair<int, int>, double> SWRProb = lp.solve_lp();
        
        //vector<double> offMass = g.poisson_offline_mass(SWRProb);
        vector<double> offMass = g.poisson_offline_mass(typeProb);
        

        vector<int> blueF, redF;
        tie(blueF, redF) = g.feldman_et_al_color();

        vector<int> blueB, redB;
        tie(blueB, redB) = g.bahmani_kapralov_color();
        
        
        vector<vector<int>> jlList = g.jaillet_lu_list();
        
        map<pair<int, int>, double> brubachLp = g.brubach_et_al_lp();
        vector<vector<pair<int, double>>> brubachSSXH = g.brubach_et_al_h(brubachLp);

        for (int j = 0; j < numSample; j++)
        {
            g.realize(realSize);

            OPT.add_run(match_size(g.maximum_matching()));
            //SWR.add_run(match_size(g.sampling_without_replacement(SWRProb)));
            SWRType.add_run(match_size(g.sampling_without_replacement(typeProb)));
            ranking.add_run(match_size(g.ranking()));
            balanceSWR.add_run(match_size(g.balance_swr()));
            balanceOCS.add_run(match_size(g.balance_ocs()));
            //poissonOCS.add_run(match_size(g.poisson_ocs(offMass, SWRProb)));
            poissonOCS.add_run(match_size(g.poisson_ocs(offMass, typeProb)));
            //topHalf.add_run(match_size(g.top_half_sampling(SWRProb)));
            topHalf.add_run(match_size(g.top_half_sampling(typeProb)));
            minDegree.add_run(match_size(g.min_degree()));
            
            feldmanMMM.add_run(match_size(g.feldman_et_al(blueF, redF)));
            BahmaniKapralov.add_run(match_size(g.bahmani_kapralov(blueB, redB)));
            manshadiGS.add_run(match_size(g.manshadi_et_al(typeProb)));
            jailletLu.add_run(match_size(g.jaillet_lu(jlList)));
            brubachSSX.add_run(match_size(g.brubach_et_al(brubachSSXH)));
        }

        double opt = compute_mean_std(OPT.resRun).first;
        for (auto j : resPointer)
            (*j).summary_run(opt);
    }
    cout << endl;

    for (auto j : resPointer) (*j).summary_sample();
}


// Generate an Erdos-Renyi bipartite graph, with
// n online types, m offline vertices and p probability of each edge
graph bipartite_erdos_renyi(int n, int m, double p)
{
    graph g(n, m);
    binomial_distribution<int> dist(m, p);
    for (int i = 0; i < n; i++)
    {
        int d = dist(rng);

        set<int> nbs;
        uniform_int_distribution<int> uni(0, m - 1);

        while ((int)nbs.size() < d)
            nbs.insert(uni(rng));

        for (auto it = nbs.begin(); it != nbs.end(); it++)
            g.add_edge(i, n + *it);
    }
    return g;
}

// Run experiments on Erdos-Renyi type bipartite graphs with parameters (n, n, c)
void work_from_erdos_renyi(int n, double c)
{
    double p = c / n;
    cout << "Working on Erdos " << n << " " << n << " " << c << endl;

    int numGraph = 10;
    int numSample = 1000;

    cout << "Rep";
    for (int i = 1; i <= numGraph; i++)
    {
        cerr << " " << i;
        graph g = bipartite_erdos_renyi(n, n, p);

        // Preprocessing
        int realSize = g.online_size();
        
        
        map<pair<int, int>, double> typeProb  = g.optimal_matching_prob(numSample, realSize);

        natural_lp lp(g.get_adj(),g.online_size());
        map<pair<int, int>, double> SWRProb = lp.solve_lp();
        
        vector<double> offMass = g.poisson_offline_mass(SWRProb);
        

        vector<int> blueF, redF;
        tie(blueF, redF) = g.feldman_et_al_color();

        vector<int> blueB, redB;
        tie(blueB, redB) = g.bahmani_kapralov_color();
        
        
        vector<vector<int>> jlList = g.jaillet_lu_list();
        
        map<pair<int, int>, double> brubachLp = g.brubach_et_al_lp();
        vector<vector<pair<int, double>>> brubachSSXH = g.brubach_et_al_h(brubachLp);

        for (int j = 0; j < numSample; j++)
        {
            g.realize(realSize);

            OPT.add_run(match_size(g.maximum_matching()));
            SWR.add_run(match_size(g.sampling_without_replacement(SWRProb)));
            SWRType.add_run(match_size(g.sampling_without_replacement(typeProb)));
            ranking.add_run(match_size(g.ranking()));
            balanceSWR.add_run(match_size(g.balance_swr()));
            balanceOCS.add_run(match_size(g.balance_ocs()));
            poissonOCS.add_run(match_size(g.poisson_ocs(offMass, SWRProb)));
            topHalf.add_run(match_size(g.top_half_sampling(SWRProb)));
            minDegree.add_run(match_size(g.min_degree()));
            
            feldmanMMM.add_run(match_size(g.feldman_et_al(blueF, redF)));
            BahmaniKapralov.add_run(match_size(g.bahmani_kapralov(blueB, redB)));
            manshadiGS.add_run(match_size(g.manshadi_et_al(typeProb)));
            jailletLu.add_run(match_size(g.jaillet_lu(jlList)));
            brubachSSX.add_run(match_size(g.brubach_et_al(brubachSSXH)));
        }

        double opt = compute_mean_std(OPT.resRun).first;
        for (auto j : resPointer)
            (*j).summary_run(opt);
    }
    
    cout << endl;

    for (auto j : resPointer) (*j).summary_sample();
}



int main()
{
    // Work on Small Erdos-Renyi type graph
    /*for (double c = 0.2; c <= 2.0; c += 0.2)
    {
        work_from_erdos_renyi(100, c);
        datasetName.push_back("c=" + to_string(c));
    }*/
    
    // Path and name of real-world datasets
    vector<pair<string, string>> file_name = 
    {
        make_pair("real_world/socfb-Caltech36/socfb-Caltech36.txt", "Caltech36"),
        make_pair("real_world/socfb-Reed98/socfb-Reed98.txt", "Reed98"),
        make_pair("real_world/bio-CE-GN/bio-CE-GN.txt", "CE-GN"),
        make_pair("real_world/bio-CE-PG/bio-CE-PG.txt", "CE-PG"),
        make_pair("real_world/econ-beause/econ-beause.txt", "beause"),
        make_pair("real_world/econ-mbeaflw/econ-mbeaflw.txt", "mbeaflw")
    };
    for (auto item : file_name)
    {
        work_from_file(item.first);
        datasetName.push_back(item.second);
    }

    
    save_results_to_files("real_world_result");
    return 0;
}