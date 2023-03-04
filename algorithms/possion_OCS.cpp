vector<int> graph::possion_OCS(const vector<double> &offMass, map<pair<int, int>, double> &typeProb)
{
    vector<int> res(realSize, -1);
    vector<bool> matched(onSize + offSize, false);
     for (int i = 0; i < realSize; i++)
    {
        double totalMass = 0.0,  mass,mass1;
        vector<pair<int, double>> validMass;
        int index = -1;
        for (int j : adj[types[i]])
        {
            mass = exp((1.0 * i / types.size()) * offMass[j]) * typeProb[make_pair(types[i], j)];
            mass1 = typeProb[make_pair(types[i],j)];
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
void graph::poisson_offline_mass(map<pair<int, int>, double> &typeProb, vector<double> &offMass)
{
    offMass.clear();
    for (int i = 0; i < onSize + offSize; i++)
        offMass.push_back(0.0);
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            offMass[j] += typeProb[make_pair(i, j)];
    return;
}