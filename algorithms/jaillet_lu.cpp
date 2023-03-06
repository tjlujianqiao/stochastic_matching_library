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

    g.max_flow();

    cycle_break_graph gCycle(onSize, onSize + offSize);
    for (int i = 0; i < onSize; i ++)
        for (auto e : g.adj[i])
            if (e.flow > 0)
                gCycle.add_edge(i, e.v, e.flow);
            
    gCycle.frac_to_int();
    gCycle.cycle_break();

    vector<vector<int>> res(onSize, vector<int>());
    
    for (int i = 0; i < onSize; i++) 
        for (auto e : gCycle.vInt[i])
            for (int j = 0; j < e.second; j++)
                res[i].push_back(e.first);
    return res;
}



vector<int> graph::jaillet_lu(vector<vector<int>> &jlList)
{
    vector<int> res(realSize, -1);
    vector<bool> matched(offSize + onSize, false);

    for (int i = 0; i < realSize; i++)
    {
        shuffle(jlList[types[i]].begin(), jlList[types[i]].end(), rng);
        for (int j : jlList[types[i]])
        {
            if (j != -1 && not matched[j])
            {
                res[i] = j;
                matched[j] = true;
                break;
            }
        }
    }
    return res;
}