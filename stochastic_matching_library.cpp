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
    vector<double> resSample;
    vector<double> resStd;
    vector<pair<double, double>> resGraph;

    resAlg(string s)
    {
        name = s;
        resRun = {};
        resSample = {};
        resGraph = {};
        resStd = {};
    }

    void add_run(double item)
    {
        resRun.push_back(item);
    }

    void summary_run(double base = 1)
    {
        auto mean_std = compute_mean_std(resRun);
        resSample.push_back(mean_std.first / base);
        resStd.push_back(mean_std.second / base);
        resRun.clear();
    }

    void summary_sample()
    {
        resGraph.push_back(compute_mean_std(resSample));
        resSample.clear();
    }

} OPT("OPT"),
    SWR("SWR"),
    ranking("Ranking"),
    minDegree("MinDegree"),
    feldmanMMM("FeldManEtAl"),
    BahmaniKapralov("BahmaniKapralov"),
    manshadiGS("ManshadiEtAl"),
    jailletLu("JailletLu"),
    topHalf("TopHalfSampling"),
    poissonOCS("PoissonOCS"),
    balanceSWR("BalanceSWR"),
    balanceOCS("BalanceOCS");

void write_to_files(resAlg *resPointer[], int resLength, string base)
{
    ofstream fileResMean, fileResStd;
    cout << "Save Results into the file:" << base + "/" + "(resMean.txt,resStd.txt)" << endl;
    fileResMean.open(base + "/" + "resMean.txt");
    fileResStd.open(base + "/" + "resStd.txt");
    fileResMean << "Algorithm" << endl;
    fileResStd << "Algorithm" << endl;
    for (int i = 0; i < resLength; i++)
    {

        double avgMean = 0.0, avgStd = 0.0;

        fileResMean << resPointer[i]->name<< ":";
        fileResMean << accumulate(resPointer[i]->resSample.begin(), resPointer[i]->resSample.end(), 0.0)/ resPointer[i]->resSample.size()  << " ";

        fileResStd << resPointer[i] -> name << ":";
        fileResStd << accumulate(resPointer[i]->resStd.begin(), resPointer[i]->resStd.end(), 0.0) / resPointer[i]->resStd.size() << " ";
    }
    fileResMean.close();
    fileResStd.close();

    cout << "Output Results Done!" << endl;
}

int main()
{
    resAlg *resPointer[] = {&OPT, &SWR, &ranking, &balanceSWR, &balanceOCS, &poissonOCS, &topHalf, &minDegree, &jailletLu, &manshadiGS, &feldmanMMM, &BahmaniKapralov};

    int numGraph = 2;
    int numSample = 100;
    int numAlg ;

    cout << "Rep";
    for (int i = 1; i <= numGraph; i++)
    {
        cout << " " << i;
        graph g = generate_from_file("real_world/bio-CE-GN/bio-CE-GN.txt");

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
        numAlg = 0 ;
        for (auto j : resPointer)
            (*j).summary_run(opt), numAlg ++;
    }
    cout << endl;

    // If you want to summarize the data on one graph, uncomment the following line
    // for (auto j : resPointer) (*j).summary_sample();

    cout << endl;
    for (auto i : resPointer) {

        cout << (*i).name << ":";
        for(auto j : (*i).resSample)
            cout<<" "<<j;

        cout<<endl;
    }

    cout<<endl;
    write_to_files(resPointer, numAlg, "real_world_result");

    return 0;
}

