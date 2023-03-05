struct cycle_break_graph{
    
    vector<map<int, int>> val;

    int size;

    cycle_break_graph(int n)
    {
        size = n;
        
        val.resize(size);
        for (int i = 0; i < size; i++)
            val[i] = {};
    }
    
    void add_edge(int x, int y, int z)
    {
        val[x][y] = z;
    }
    
    inline int cycle_type(int u1, int u2, int v1, int v2)
    {
        return 7 - val[u1][v1] - val[u1][v2] - val[u2][v1] - val[u2][v2];
    }
    
    void break_one_cycle(int u1, int u2, int v1, int v2)
    {
        if (val[u1][v1] + val[u1][v2] < val[u2][v1] + val[u2][v2])
            swap(u1, u2);
        
        if (val[u1][v1] < val[u1][v2])
            swap(v1, v2);
        
        if (cycle_type(u1, u2, v1, v2) == 2)    //Type C2
        {
            val[u1][v1] = 1;
            val[u1][v2] = 2;
            val[u2][v1] = 1;
            val[u2].erase(val[u2].find(v2));
        }
        else                                    //Type C3
        {
            val[u1][v1] = 2;
            val[u1].erase(val[u1].find(v2));
            val[u2].erase(val[u2].find(v1));
            val[u2][v2] = 2;
        }
    }
    
    int detect_cycle(int &U1, int &U2, int &V1, int &V2)
    {
        map<pair<int, int>, vector<int>> M = {};
        
        int flag = 0;
        
        for (int u1 = 0; u1 < size; u1 ++)
            for (auto e1 : val[u1]) for (auto e2 : val[u1])
            {
                int v1 = e1.first, v2 = e2.first; 
                if (v1 < v2)
                {
                    if (M.count(make_pair(v1, v2)))
                    {
                        for (int u2 : M[make_pair(v1, v2)])
                        {
                            int type = cycle_type(u1, u2, v1, v2);
                            if (type == 2 || type == 3)
                            {
                                U1 = u1, U2 = u2, V1 = v1, V2 = v2;
                                if (type == 2) return 1;
                                flag = 1;
                            }
                        }
                    
                        M[make_pair(v1, v2)].push_back(u1);
                    }
                    else 
                        M[make_pair(v1, v2)] = {u1};
                }
            }
        
        return flag;
    }
    
    void cycle_break()
    {
        int u1, u2, v1, v2;
        while (detect_cycle(u1, u2, v1, v2))
        {
            break_one_cycle(u1, u2, v1, v2);
        }
    }
};
