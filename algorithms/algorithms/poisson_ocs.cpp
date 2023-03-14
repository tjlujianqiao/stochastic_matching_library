// Match with offline mass and weight x_{ij} for each edge
vector<int> graph::poisson_ocs(const vector<double> &offMass, map<pair<int, int>, double> &typeProb)
{
    vector<int> res(realSize, -1);
    vector<bool> matched(onSize + offSize, false);
     for (int i = 0; i < realSize; i++)
    {
        double totalMass = 0.0, mass;
        vector<pair<int, double>> validMass;
        int index = -1;
        for (int j : adj[types[i]])
        {
            mass = exp((1.0 * i / types.size()) * offMass[j]) * typeProb[make_pair(types[i], j)];
            if (not matched[j] and mass > 0.0)
            {
                totalMass += mass;
                validMass.push_back(make_pair(j, mass));
            }
        }
        if (validMass.size())
        {
            std::uniform_real_distribution<double> distr(0, totalMass);
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

// Compute mass of each offline vertex
vector<double> graph::poisson_offline_mass(map<pair<int, int>, double> &typeProb)
{
    vector<double> offMass = {};
    for (int i = 0; i < onSize + offSize; i++)
        offMass.push_back(0.0);
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            offMass[j] += typeProb[make_pair(i, j)];
    return offMass;
}