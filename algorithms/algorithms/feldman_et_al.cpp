// Compute edge colors in Feldman et al. (2009)
pair<vector<int>, vector<int>> graph::feldman_et_al_color()
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
    
    decomposite_graph gDecom(onSize + offSize, onSize);
    for (int i = 0; i < onSize; i++)
        for (auto e : g.adj[i])
            if (e.flow > 0)
                gDecom.add_edge(i, e.v);
                
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
vector<int> graph::feldman_et_al(vector<int> &blue, vector<int> &red)
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