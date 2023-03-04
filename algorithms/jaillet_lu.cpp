vector<vector<int>> graph::jaillet_lu_list()
{
    int s = onSize + offSize, t = s + 1;
    flow_graph g(s, t);
    
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            g.add_edge(i, j, 2);

    for (int i = 0; i < onSize; i++)
        g.add_edge(s, i, 3);

    for (int j = onSize; j < onSize + offSize; j++)
        g.add_edge(j, t, 3);

    g.maxflow();

    vector<vector<int>> res(onSize, vector<int>());
    
    for (int i = 0; i < onSize; i++) 
        for (auto e : g.adj[i])
            for (int j = 0; j < e.flow; j++)
                res[i].push_back(e.v);
    return res;
}



vector<int> graph::jaillet_lu(vector<vector<int>> &jlList)
{
    vector<int> res(types.size(), -1);
    vector<bool> offLine(adj.size(), false);
    mt19937 rng(random_device{}());

    for (int i = 0; i < onSize; i++)
    {
        shuffle(jlList[types[i]].begin(), jlList[types[i]].end(), rng);
        for (int j : jlList[types[i]])
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