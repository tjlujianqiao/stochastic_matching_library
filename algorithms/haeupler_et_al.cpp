tuple<vector<int>, vector<int>, vector<vector<pair<double, int>>>> graph::haeupler_et_al_advice(map<pair<int, int>, double> lpPseudo)
{
    vector<int> M1(onSize + offSize, -1), M2(onSize + offSize, -1);
    int s = onSize + offSize, t = s + 1;
    flow_graph g1(s, t);

    for (int i = 0; i < onSize; i++)
    {
        g1.add_edge(s, i, 1);
        for (int j : adj[i])
            g1.add_edge(i, j, 1);
    }
    for (int j = onSize; j < onSize + offSize; j++)
        g1.add_edge(j, t, 1);

    g1.max_flow();
    for (int i = 0; i < onSize; i++)
    {
        for (auto e : g1.adj[i])
            if (e.flow > 0)
                M1[i] = e.v, M1[e.v] = i, lpPseudo[make_pair(i, e.v)] = 0;
    }
    flow_graph g2(s, t);
    for (int i = 0; i < onSize; i++)
    {
        g2.add_edge(s, i, 1);
        for (int j : adj[i])
            if (M1[i] != j)
                g2.add_edge(i, j, 1);
    }
    for (int j = onSize; j < onSize + offSize; j++)
    {
        g2.add_edge(j, t, 1);
    }
    g2.max_flow();
    for (int i = 0; i < onSize; i++)
    {
        for (auto e : g2.adj[i])
        {
            if (e.flow > 0)
                M2[i] = e.v, M2[e.v] = i, lpPseudo[make_pair(i, e.v)] = 0;
        }
    }
    for (auto &e : lpPseudo)
    {
        e.second *= 2;
    }

    vector<vector<pair<double, int>>> Prob(onSize);
    for (int j = onSize; j < onSize + offSize; j++)
    {
        double totalMass = 0.0;
        for (int i : adj[j])
        {
            double mass = lpPseudo[make_pair(i, j)], addedMass;
            if (mass > 1e-5 and totalMass < 1.0 + 1e-9)
            {
                addedMass = min(1.0 - totalMass, mass);
                Prob[i].push_back(make_pair(addedMass, j)), totalMass += addedMass;
            }
        }
    }

    for (int i = 0; i < onSize; i++)
    {

        double totalMass = 0.0;
        for (int j = 0; j < (int)Prob[i].size(); j++)
        {
            totalMass += Prob[i][j].first;
            Prob[i][j].first = totalMass;
        }
        Prob[i].push_back(make_pair(2.1, -1));
    }
    return make_tuple(M1, M2, Prob);
}

vector<int> graph::haeupler_et_al(const vector<int> &M1, const vector<int> &M2, const vector<vector<pair<double, int>>> &Prob)
{
    int offlineDegree[onSize + offSize];
    memset(offlineDegree, 0, sizeof(offlineDegree));
    for (int j = onSize; j < onSize + offSize; j++)
    {
        if (M1[j] != -1)
            offlineDegree[j]++;
        if (M2[j] != -1)
            offlineDegree[j]++;
    }


    vector<int> res(realSize, -1);
    vector<bool> matched(onSize + offSize, false);
    uniform_real_distribution<double> rand_sample(0, 2);

    int c1, c2;
    double r1, r2;
    for (int i = 0; i < realSize; i++)
    {
        r1 = rand_sample(rng);
        r2 = rand_sample(rng);
        pair<double, int> val1 = make_pair(r1, -1);
        pair<double, int> val2 = make_pair(r2, -1);
        c1 = (*upper_bound(Prob[types[i]].begin(),
                           Prob[types[i]].end(), val1))
                 .second;
        c2 = (*upper_bound(Prob[types[i]].begin(),
                           Prob[types[i]].end(), val2))
                 .second;

        if (M1[i] != -1 and not matched[M1[i]])
            res[i] = M1[i];
        else if (M2[i] != -1 and not matched[M2[i]])
            res[i] = M2[i];

        else
        {
            if ((c1 != -1 and c2 != -1) and offlineDegree[c1] > offlineDegree[c2])
                swap(c1, c2);


            if (c1 != -1 and not matched[c1])
                res[i] = c1;
            else if (c2 != -1 and not matched[c2])
                res[i] = c1;
        }
    }
    return res;
}
