vector<vector<int>> graph::jaillet_lu_prob()
{
    map<pair<int, int>, int> rGraph;
    pair<int, int> e;
    int s = adj.size(), t = adj.size() + 1;

    for (int i = 0; i < onSize; i++)
        rGraph[make_pair(s,i)] = 3;

    for (int j = onSize; j < onSize + offSize; j++)
        rGraph[make_pair(j,t)] = 3;

    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            rGraph[make_pair(i,j)] = 2;

    vector<vector<int>> adj2 = adj;
    adj2.push_back(vector<int>());
    adj2.push_back(vector<int>());
    for (int i = 0; i < onSize; i++)
    {
        adj2[i].push_back(s);
        adj2[s].push_back(i);
    }
    for (int i = onSize; i < onSize + offSize; i++)
    {
        adj2[i].push_back(t);
        adj2[t].push_back(i);
    }
    // fordFulkerson(rGraph, adj2, s, t);

    vector<vector<int>> res(onSize, vector<int>());
    for (int i = 0; i < onSize; i++)
    {
        for (auto j : adj[i])
        {
            e = make_pair(i, j);
            if (rGraph[e] < 2)
            {
                for (int k = 0; k < 2 - rGraph[e]; k++)
                    res[i].push_back(j);
            }
        }
        while (res[i].size() < 3)
            res[i].push_back(-1);
    }
    return res;
}



vector<int> graph::jaillet_lu_matching( vector<vector<int>> &jlProb)
{
    vector<int> res(types.size(), -1);
    vector<bool> offLine(adj.size(), false);
    mt19937 rng(random_device{}());

    for (int i = 0; i < onSize; i++)
    {
        shuffle(jlProb[types[i]].begin(), jlProb[types[i]].end(), rng);
        for (int j : jlProb[types[i]])
        {
            if (j != -1 && not offLine[j])
            {
                res[i] = i;
                offLine[j] = true;
                break;
            }
        }
    }
    return res;
}