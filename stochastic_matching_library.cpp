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
#include "cycle_break_graph.h"
#include "decomposite_graph.h"
#include "read_file.cpp"
#include "algorithms/algorithms.h"

int match_size(const vector<int> &res)
{
    int match = 0;
    for (int i = 0; i < res.size(); i++)
        if (res[i] != -1)
            match++;
    return match;
}

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

struct resAlg
{
    string name;
    vector<double> resRun;
    vector<pair<double, double>> resSample;
    vector<pair<double, double>> resGraph;

    resAlg(string s)
    {
        name = s;
        resRun = {};
        resSample = {};
        resGraph = {};
    }

    void add_run(double item)
    {
        resRun.push_back(item);
    }

    void summary_run(double base = 1)
    {
        auto mean_std = compute_mean_std(resRun);
        resSample.push_back(make_pair(mean_std.first / base, mean_std.second / base));
        resRun.clear();
    }

    void summary_sample()
    {
        if (resSample.size() == 1)
        {
            resGraph.push_back(resGraph[0]);
            resSample.clear();
            return;
        }
        
        
        vector<double> res;
        for (auto item : resSample)
            res.push_back(item.first);
        
        resGraph.push_back(compute_mean_std(res));
        resSample.clear();
    }

}   OPT("OPT"),
    SWR("SWR"),
    ranking("Ranking"),
    balanceSWR("BalanceSWR"),
    balanceOCS("BalanceOCS"),
    minDegree("MinDegree"),
    feldmanMMM("FeldManEtAl"),
    BahmaniKapralov("BahmaniKapralov"),
    manshadiGS("ManshadiEtAl"),
    jailletLu("JailletLu"),
    topHalf("TopHalfSampling"),
    poissonOCS("PoissonOCS");

const vector<resAlg*> resPointer = {&OPT, &SWR, &poissonOCS, &minDegree, &balanceSWR, &balanceOCS, &ranking, &topHalf, &jailletLu, &manshadiGS, &BahmaniKapralov, &feldmanMMM};


const vector<pair<string, string>> file_name = 
{
    make_pair("real_world/socfb-Caltech36/socfb-Caltech36.txt", "Caltech36"),
    make_pair("real_world/socfb-Reed98/socfb-Reed98.txt", "Reed98"),
    make_pair("real_world/bio-CE-GN/bio-CE-GN.txt", "CE-GN"),
    make_pair("real_world/bio-CE-PG/bio-CE-PG.txt", "CE-PG"),
    make_pair("real_world/econ-beause/econ-beause.txt", "beause"),
    make_pair("real_world/econ-mbeaflw/econ-mbeaflw.txt", "mbeaflw")
};


void write_to_files(string directory)
{
    ofstream fileResMean, fileResStd;
    cout << "Save Results into the file:" << directory + "/" + "(resMean.txt,resStd.txt)" << endl;
    fileResMean.open(directory + "/" + "resMean.txt");
    fileResStd.open(directory + "/" + "resStd.txt");
    
    fileResMean << "Algorithm";
    for (auto item : file_name)
        fileResMean << " " << item.second;
    fileResMean << endl;
    
    fileResStd << "Algorithm";
    for (auto item : file_name)
        fileResStd << " " << item.second;
    fileResStd << endl;
    
    for (auto i : resPointer)
        if (i != &OPT)
        {

            double avgMean = 0.0, avgStd = 0.0;

            fileResMean << (*i).name << " ";
            fileResStd << (*i).name << " ";
            
            for (auto j : (*i). resGraph)
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

void work_from_file(string name)
{
    cout << "Working on file " << name << endl;

    int numGraph = 10;
    int numSample = 100;

    cout << "Rep";
    for (int i = 1; i <= numGraph; i++)
    {
        cout << " " << i;
        graph g = generate_from_file(name, true);

        // Preprocessing
        int realSize = g.online_size();
        map<pair<int, int>, double> typeProb = g.optimal_matching_prob(numSample, realSize);
        vector<double> offMass = g.poisson_offline_mass(typeProb);
        vector<vector<int>> jlList = g.jaillet_lu_list();

        vector<int> blueF, redF;
        tie(blueF, redF) = g.feldman_et_al_color();

        vector<int> blueB, redB;
        tie(blueB, redB) = g.bahmani_kapralov_color();

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
            BahmaniKapralov.add_run(match_size(g.bahmani_kapralov(blueB, redB)));
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
    for (auto item : file_name)
        work_from_file(item.first);
    
    write_to_files("real_world_result");
    
    return 0;
}
