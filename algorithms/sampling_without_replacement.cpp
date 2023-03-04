vector<int> graph::sampling_without_replacement(map<pair<int, int>, double> &typeProb)
{
    vector<int> res(realSize, -1);
    vector<bool> matched(onSize + offSize, false);
    // ith means ith arrival online vertex while i means type  i
    for (int i = 0; i < realSize; i++)
    {
        double totalMass = 0.0, sum = 0.0, mass;
        vector<pair<int, double>> validMass;
        int index = -1, flag = -1;
        for (int j : adj[types[i]])
        {
            mass = typeProb[make_pair(types[i], j)];
            if (not matched[j] and mass > 0.0)
            {
                totalMass += mass;
                validMass.push_back(make_pair(j, mass));
                flag = 1;
            }
        }
        if (flag == 1)
        {
            std::uniform_real_distribution<double> distr(0, totalMass);
            double rr = distr(rng);
            for (auto item : validMass)
            {
                sum += item.second;
                if (sum >= rr)
                {
                    index = item.first;
                    break;
                }
            }
            if (index != -1)
            {
                res[i] = index;
                matched[index] = true;
            }
        }
    }
    return res;
}




map<pair<int, int>, double> graph::optimal_matching_prob(int numSample, int realSize)
{
    uniform_int_distribution<int> uni(0, onSize - 1);
    map<pair<int, int>, double> Prob;
    vector<int> res;
    for (int count = 0; count < numSample; count++)
    {
        realize(realSize);
        res = maximum_matching();
        for (int i = 0; i < realSize; i++)
        {
            if (res[i] != -1)
                Prob[make_pair(types[i], res[i])] += 1.0 / numSample;
        }
    }
    return Prob;
}

