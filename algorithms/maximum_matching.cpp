// Compute offline optimal matching in realization graph
vector<int> graph::maximum_matching()
{
    int s = realSize + offSize, t = s + 1;
    flow_graph g(s, t);
    
    for (int i = 0; i < realSize; i++)
        for (int j : adj[types[i]])
            g.add_edge(i, j - onSize + realSize, 1);
        
    for (int i = 0; i < realSize; i++)
        g.add_edge(s, i, 1);
    
    for (int j = realSize; j < realSize + offSize; j++)
        g.add_edge(j, t, 1);
    
    g.max_flow();
    
    vector<int> res(realSize, -1);
    for(int i = 0; i < realSize; i++)
        for (auto e : g.adj[i])
            if (e.flow > 0)
                res[i] = e.v - realSize + onSize;
            
    return res;
}
