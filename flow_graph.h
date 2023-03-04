//Compute maximum flow by Dinic's algorithm
//Assume t is the vertex with largest label 

struct flow_graph{
    
    struct edge
    {
        int v, cap, flow, rev;
        edge(int v, int cap, int flow, int rev) : v(v), cap(cap), flow(flow), rev(rev){}
    };
    

    
    vector<vector<edge>> adj;
    
    
    vector<int> dep, cur;
    
    
    int s, t;
    
    const int inf = 1e9;
    
    

    flow_graph(int S, int T)
    {
        s = S, t = T;
        adj.resize(t + 1);
        for (int i = 0; i <= t; i++)
            adj[i] = {};
    }

    void add_edge(int x, int y, int cap)
    {
        int x_e = adj[x].size(), y_e = adj[y].size();
        adj[x].push_back(edge(y, cap, 0, y_e));
        adj[y].push_back(edge(x, 0, 0, x_e));
    }

    bool bfs()
    {
        queue<int> q;
        fill(dep.begin(), dep.end(), 0);
        dep[s] = 1;
        q.push(s);
        do
        {
            int u = q.front();
            q.pop();
            for (int u_e = 0; u_e < adj[u].size(); u_e++){
                int v = adj[u][u_e].v;
                if (!dep[v] && adj[u][u_e].cap){
                    dep[v] = dep[u] + 1;
                    q.push(v);
                }
            }
        }
        while(!q.empty());
        
        return dep[t];
    }
    
    int dfs(int u, int num)
    {
        if ( (u == t) || !num ) return num;
        for (int &u_e = cur[u]; u_e < adj[u].size(); u_e++)
        {
            int v = adj[u][u_e].v, v_e = adj[u][u_e].rev;
            if ((dep[v] == dep[u] + 1) && adj[u][u_e].cap)
            {
                int d = dfs(v, min(num, adj[u][u_e].cap));
                if (d) {
                    adj[u][u_e].cap -= d;
                    adj[u][u_e].flow += d;
                    
                    adj[v][v_e].cap += d;
                    adj[v][v_e].flow -= d;
                    return d;
                }
            }
        }
        return 0;
    }
    
    
    
    void max_flow()
    {
        dep.resize(t + 1);
        cur.resize(t + 1);
        while (bfs()) {
            fill(cur.begin(), cur.end(), 0);
            while (int d = dfs(s, inf));
        }
    }
    
    
    
    vector<int> min_cut()
    {
        vector<int> inS(t + 1, 0);
        queue<int> Q;
        Q.push(s);
        inS[s] = 1;
        
        while(!Q.empty())
        {
            int u = Q.front();
            Q.pop();
            
            for (int u_e = 0; u_e < adj[u].size(); u_e++){
                int v = adj[u][u_e].v;
                if (!inS[v] && adj[u][u_e].cap){
                    inS[v] = 1;
                    Q.push(v);
                }
            }
        }
        return inS;
    }
    
};
