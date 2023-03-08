// Perform cycle break for type graph in Jaillet and Lu (2013) and Brubach et al. (2016)

struct cycle_break_graph{
    
    // Fractional value of each edge
    vector<map<int, double>> vFrac;
    
    // Rounded value of each edge
    vector<map<int, int>> vInt;

    // Total number of vertices
    int size;
    // Number of online types
    int onSize;

    // Construct an empty type graph with:
    // n online types and m offline vertices
    cycle_break_graph(int n, int m)
    {
        onSize = n, size = m;
        
        vFrac.resize(size);
        vInt.resize(size);
        
        for (int i = 0; i < size; i++)
        {
            vFrac[i] = {};
            vInt[i] = {};
        }
    }
    
    // Add an edge (x, y) with value z
    void add_edge(int x, int y, double z)
    {
        vFrac[x][y] = z;
        vFrac[y][x] = z;
    }
    
    // Decide whether the number is fractional
    bool is_fractional(double z){
        return abs(z - floor(z + 0.5)) > 1e-10;
    }
    
    // Rounding all fractional values to integrals
    void frac_to_int()
    {
        for (int i = 0; i < size; i++)
            for (auto e : vFrac[i])
                if (e.second > 1e-10)
                {
                    int j = e.first; double x = e.second;
                    vInt[i][j] = (int)floor(x + 0.5);
                }
    }
    
    // Find a cycle or maximal path with all edges fractional
    vector<int> find_fractional_cycle_or_path()
    {
        int start = -1;
        for (int i = 0; i < size; i++)
        {
            int num_frac = 0;
            for (auto e : vFrac[i])
                if (is_fractional(e.second))
                    num_frac++;
            
            if (num_frac >= 1)
                start = i;
            if (num_frac == 1)
                break;
        }
    
        if (start == -1) return {};
        
        vector<int> visit(size, 0);
        
        vector<int> vList(1, start);
        visit[start] = 1;
        int cur = start, pre = -1;
        
        
        bool flag;
        do
        {
            flag = false;
            for (auto e : vFrac[cur])
                if (e.first != pre && (is_fractional(e.second)))
                {
                    int next = e.first;
                    if (!visit[next])
                    {
                        visit[next] = 1;
                        vList.push_back(next);
                        pre = cur, cur = next;
                        flag = true;
                        break;
                    }
                    else
                    {
                        vList.push_back(next);
                        while (*vList.begin() != next)
                            vList.erase(vList.begin());
                        break;
                    }
                }
        }
        while(flag);
        
        return vList;
    }
    
    // Apply Gandhi et. al (2006)'s dependent rounding
    void gandhi_et_al_rounding()
    {
        do
        {
            vector<int> vList = find_fractional_cycle_or_path();
            
            if (vList.empty())
                break;
            
            double alpha = 1, beta = 1;
            
            for (int i = 0; i < (int)vList.size() - 1; i++)
            {
                int x = vList[i], y = vList[i + 1];
                double e = vFrac[x][y];
                if (i % 2 == 0)
                {
                    alpha = min(alpha, ceil(e) - e);
                    beta = min(beta, e - floor(e));
                }
                else
                {
                    alpha = min(alpha, e - floor(e));
                    beta = min(beta, ceil(e) - e);
                }
            }
            
            uniform_real_distribution<double> curRand(0.0, alpha + beta);
            
            if (curRand(rng) < beta)
                for (int i = 0; i < (int)vList.size() - 1; i++)
                {
                    int x = vList[i], y = vList[i + 1];
                    if (i % 2 == 0)
                    {
                        vFrac[x][y] += alpha;
                        vFrac[y][x] += alpha;
                    }
                    else
                    {
                        vFrac[x][y] -= alpha;
                        vFrac[y][x] -= alpha;
                    }
                }
            else
                for (int i = 0; i < (int)vList.size() - 1; i++)
                {
                    int x = vList[i], y = vList[i + 1];
                    if (i % 2 == 0)
                    {
                        vFrac[x][y] -= beta;
                        vFrac[y][x] -= beta;
                    }
                    else
                    {
                        vFrac[x][y] += beta;
                        vFrac[y][x] += beta;
                    }
                }
            
        } while(true);
    }
    
    //Decide a cycle's type
    inline int cycle_type(int u1, int u2, int v1, int v2)
    {
        return 7 - vInt[u1][v1] - vInt[u1][v2] - vInt[u2][v1] - vInt[u2][v2];
    }
    
    //Break a cycle
    void break_one_cycle(int u1, int u2, int v1, int v2)
    {
        if (vInt[u1][v1] + vInt[u1][v2] < vInt[u2][v1] + vInt[u2][v2])
            swap(u1, u2);
        
        if (vInt[u1][v1] < vInt[u1][v2])
            swap(v1, v2);
        
        if (cycle_type(u1, u2, v1, v2) == 2)    //Type C2
        {
            vInt[u1][v1] = vInt[v1][u1] = 1;
            vInt[u1][v2] = vInt[v2][u1] = 2;
            vInt[u2][v1] = vInt[v1][u2] = 1;
            vInt[u2].erase(vInt[u2].find(v2));
            vInt[v2].erase(vInt[v2].find(u2));
        }
        else                                    //Type C3
        {
            vInt[u1][v1] = vInt[v1][u1] = 2;
            vInt[u1].erase(vInt[u1].find(v2));
            vInt[v2].erase(vInt[v2].find(u1));
            vInt[u2].erase(vInt[u2].find(v1));
            vInt[v1].erase(vInt[v1].find(u2));
            vInt[u2][v2] = vInt[v2][u2] = 2;
        }
    }
    
    //Detect cycles of type C2 and C3
    int detect_cycle(int &U1, int &U2, int &V1, int &V2)
    {
        map<pair<int, int>, vector<int>> M = {};
        
        int flag = 0;
        
        for (int u1 = 0; u1 < onSize; u1 ++)
            for (auto e1 : vInt[u1]) for (auto e2 : vInt[u1])
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
    
    // Break all cycles of type C2 and C3
    void cycle_break()
    {
        int u1, u2, v1, v2;
        while (detect_cycle(u1, u2, v1, v2))
            break_one_cycle(u1, u2, v1, v2);
    }
    

};
