// Compute lists in Jaillet and Lu (2013)
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
    for (int i = 0; i < onSize; i++)
        for (auto e : g.adj[i])
            if (e.flow > 0)
                gCycle.add_edge(i, e.v, e.flow);

    // Break cycles of type C2 and C3 required for theoretical guarantee
    gCycle.frac_to_int();
    gCycle.cycle_break();

    vector<vector<int>> res(onSize, vector<int>());

    for (int i = 0; i < onSize; i++)
        for (auto e : gCycle.vInt[i])
            for (int j = 0; j < e.second; j++)
                res[i].push_back(e.first);
    return res;
}

// Match with list of each online types
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

// Compute LP solution for non-integral algorithm in Jaillet and Lu (2013)
map<pair<int, int>, double> graph::jaillet_lu_non_integral()
{
    // adding jb for each offline vertex j
    int s = onSize + 2 * offSize, t = s + 1;
    int mul = 1e9, ln = 306852819;
    flow_graph g(s, t);

    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
        {
            int jb = j + offSize;
            g.add_edge(i, j, mul / 2);
            g.add_edge(i, jb, mul);
            g.add_edge(jb, j, ln);
        }

    for (int i = 0; i < onSize; i++)
        g.add_edge(s, i, mul);

    for (int j = onSize; j < onSize + offSize; j++)
        g.add_edge(j, t, mul);
    g.max_flow();
    map<pair<int, int>, double> jlProb;
    for (int i = 0; i < onSize; i++)
        for (auto e : g.adj[i])
            if (e.flow > 0)
            {
                int j = (e.v < onSize + offSize) ? e.v : e.v - offSize;
                jlProb[make_pair(i, j)] += (double)e.flow / mul;
            }

    return jlProb;
}