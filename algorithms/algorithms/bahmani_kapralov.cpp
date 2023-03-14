// Compute edge colors in Bahmani and Kapralov (2010)
pair<vector<int>, vector<int>> graph::bahmani_kapralov_color()
{
    vector<int> blue(onSize, -1);
    vector<int> red(onSize, -1);
    
    int s = onSize + offSize, t = s + 1;
    flow_graph g(s, t);
    
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            g.add_edge(i, j, 1);

    for (int i = 0; i < onSize; i++)
        g.add_edge(s, i, 2);

    for (int j = onSize; j < onSize + offSize; j++)
        g.add_edge(j, t, 2);
    
    
    g.max_flow();
    
    vector<int> inS = g.min_cut();
    
    decomposite_graph gDecom(onSize + offSize, onSize);
            
    //Build graph Gs and Gt
    flow_graph gS(s, t);
    flow_graph gT(s, t);
    for (int i = 0; i < onSize; i++)
        for (auto e : g.adj[i])
        {
            int j = e.v;    //Edge (i, j)
            if (j < onSize + offSize && inS[i] == 1 && inS[j] == 1)
            {
                if (e.flow > 0) gS.add_edge(j, i, 1);
                else gS.add_edge(i, j, 1);
                    
            }
            if (j < onSize + offSize && inS[i] == 0 && inS[j] == 0)
            {
                if (e.flow > 0) gT.add_edge(i, j, 1);
                else gT.add_edge(j, i, 1);
            }
        }
    
    for (auto e : g.adj[s])
    {
        int i = e.v;
        if (inS[i] == 1)
        {
            if (e.flow == 0) gS.add_edge(s, i, 1);
            if (e.flow == 2) gS.add_edge(i, t, 1);
        }
    }
    
    for (auto e : g.adj[t])
    {
        int j = e.v;
        if (inS[j] == 1)
        {
            if (e.flow == 0) gT.add_edge(s, j, 1);
            if (e.flow == -2) gT.add_edge(j, t, 1);
        }
    }
    
    gS.max_flow();
    gT.max_flow();
    
    vector<map<int, int>> flowSum(t + 1);
    for (int i = 0; i <= t; i++)
    {
        flowSum[i] = {};
        
        for (auto e : g.adj[i])
            flowSum[i][e.v] += e.flow;
        
        for (auto e : gS.adj[i])
            flowSum[i][e.v] += e.flow;
        
        for (auto e : gT.adj[i])
            flowSum[i][e.v] += e.flow;
    }
    
    
    for (int i = 0; i < onSize; i++)
        for (auto e : flowSum[i])
        {
            int j = e.first;
            
            if (j < onSize + offSize && e.second)
                gDecom.add_edge(i, j);
        }
                
    gDecom.set_color();
    
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
        {
            if (gDecom.color[i][j] == 1)
                blue[i] = j;
            if (gDecom.color[i][j] == 2)
                red[i] = j;
        }
        
    return make_pair(blue, red);
}

// Match online vertices with advice of blue and red edges
vector<int> graph::bahmani_kapralov(vector<int> &blue, vector<int> &red)
{
    vector<int> res(realSize, -1);
    vector<int> num(onSize, 0);
    vector<bool> matched(offSize + onSize, false);

    for (int i = 0; i < realSize; i++)
    {
        int x = blue[types[i]], y = red[types[i]];
        
        num[types[i]]++;
        
        if (num[types[i]] == 1 && x != -1 && !matched[x])
        {
            res[i] = x;
            matched[x] = true;
        }
        
        if (num[types[i]] == 2 && y != -1 && !matched[y])
        {
            res[i] = y;
            matched[y] = true;
        }
    }
    return res;
}