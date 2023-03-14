// Match with the sampled matching probability of each edge
vector<int> graph::correlated_sampling(map<pair<int, int>, double> &typeProb)
{
    vector<int> res(realSize, -1);
    vector<int> offLine(onSize + offSize, -1);
    vector<vector<pair<double, int>>> Prob;
    vector<pair<double, int>> p1;
    double mass;
    int jStar = -1;
    for (int i = 0; i < realSize; i++)
    {
        p1.clear();
        for (auto j : adj[types[i]])
        {
            mass = typeProb[make_pair(types[i], j)];
            if (mass > 0.5 + 1e-10)
                jStar = j;
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
    int c1 = -1, c2 = -1;
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
        if (c1 != jStar || jStar == -1)
        {
            c2 = (*upper_bound(Prob[i].begin(),
                               Prob[i].end(), val2))
                     .second;
        }
        else
        {
            double massJStar = typeProb[make_pair(types[i], jStar)];
            vector<pair<double, int>> Prob2;
            for (auto j : adj[types[i]])
            {
                mass = typeProb[make_pair(types[i], j)];
                if (mass < 0.5 - 1e-10)
                    Prob2.push_back(make_pair(mass / (1.0 - massJStar), j));
            }
            double rr = rand_sample(rng), sum = 0;
            for (auto item : Prob2)
            {
                sum += item.first;
                if (sum >= rr)
                {
                    c2 = item.second;
                    break;
                }
            }
        }
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