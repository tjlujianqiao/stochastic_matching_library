// Match with the sampled matching probability of each edge
vector<int> graph::manshadi_et_al(map<pair<int, int>, double> &typeProb)
{
    vector<int> res(realSize, -1);
    vector<int> offLine(onSize + offSize, -1);
    vector<vector<pair<double, int>>> Prob;
    vector<pair<double, int>> p1;
    double mass;
    for (int i = 0; i < realSize; i++)
    {
        p1.clear();
        for (auto j : adj[types[i]])
        {
            mass = typeProb[make_pair(types[i], j)];
            p1.push_back(make_pair(mass, j));
        }
        Prob.push_back(p1);
    }

    for (int i = 0; i < realSize; i++)
    {
        double total = 0.0;
        for (int j = 0; j < (int)Prob[i].size(); j++)
            total += Prob[i][j].first;

        if (total < 1.1)
            Prob[i].push_back(make_pair(1.1 - total, -1));

        total = 0.0;
        for (int j = 0; j < (int)Prob[i].size(); j++)
        {
            total += Prob[i][j].first;
            Prob[i][j].first = total;
        }
    }
    uniform_real_distribution<double> rand_sample(0, 1);
    int c1, c2;
    double r1, r2;

    for (int i = 0; i < realSize; i++)
    {
        r1 = rand_sample(rng);
        r2 = r1 > 0.5 ? r1 - 0.5 : r1 + 0.5;
        pair<double, int> val1 = make_pair(r1, -1);
        pair<double, int> val2 = make_pair(r2, -1);
        // upperbound help return the first element bigger than a input val, if element is a pair, default to compare the first element of the pair
        c1 = (*upper_bound(Prob[i].begin(),
                           Prob[i].end(), val1))
                 .second;
        c2 = (*upper_bound(Prob[i].begin(),
                           Prob[i].end(), val2))
                 .second;
        if (c1 != -1 && offLine[c1] == -1)
        {
            res[i] = c1;
            offLine[c1] = i;
        }
        else if (c2 != -1 && offLine[c2] == -1)
        {
            res[i] = c2;
            offLine[c2] = i;
        }
    }
    return res;
}