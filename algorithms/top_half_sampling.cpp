// Match by top-half sampling
vector<int> graph::top_half_sampling(map<pair<int, int>, double> &typeProb)
{
    vector<int> res(realSize, -1);
    vector<bool> matched(onSize + offSize, false);
    // ith means ith arrival online vertex while i means type  i
    for (int i = 0; i < realSize; i++)
    {
        double totalMass = 0.0,  mass;
        vector<pair<int, double>> validMass;
        int index = -1;
        for (int j : adj[types[i]])
        {
            mass = typeProb[make_pair(types[i], j)];
            if (not matched[j] and mass > 0.0)
            {
                totalMass += mass;
                validMass.push_back(make_pair(j, mass));
            }
        }
        if (validMass.size())
        {
            std::uniform_real_distribution<double> distr(0, 1.0 / 2);
            double rr = distr(rng), sum = 0;
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
