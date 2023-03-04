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

    g.maxflow();
    
    vector<int> visit(t + 1, 0);
    
    for (int i = 0; i < onSize + offSize; i++)
        if (visit[i] < 2)
        {
            int cur = i, pre = -1, flag;
            visit[i] = 1;
            do
            {
                flag = false;
                for (auto e : g.adj[cur])
                    if (e.flow && e.v < onSize + offSize && e.v != pre)
                    {
                        int next = e.v;
                        if (!visit[next])
                        {
                            visit[next] = 1;
                            pre = cur, cur = next;
                            flag = true;
                            break;
                        }
                    }
            }while (flag);
            
            int start = cur;
            vector<int> vList;
            vList.push_back(start);
            visit[start] = 2;
            pre = -1;
            bool isCycle = false;
            
            do
            {
                flag = false;
                for (auto e : g.adj[cur])
                    if (e.flow && e.v < onSize + offSize && e.v != pre)
                    {
                        int next = e.v;
                        if (visit[next] < 2)
                        {
                            visit[next] = 2;
                            vList.push_back(next);
                            pre = cur, cur = next;
                            flag = true;
                            break;
                        }
                        else if (next == start)
                        {
                            vList.push_back(next);
                            isCycle = true;
                            break;
                        }
                    }
            }
            while(flag);
            
            for (int j = 0; j < vList.size() - 1; j++)
            {
                int x = vList[j], y = vList[j + 1];
                if (x > y) swap(x, y);
                
                if (isCycle) //Cycle
                    if (j % 2 == 0) blue[x] = y; else red[x] = y;
                else if (vList.size() % 2 == 0) //Odd-length paths
                    if (j % 2 == 0) blue[x] = y; else red[x] = y;
                else if (vList[0] >= onSize) //Even-length paths starting from an offline vertex
                    if (j % 2 == 0) blue[x] = y; else red[x] = y;
                else                         //Even-length paths starting from an online vertex
                    if (j % 2 == 1 || j == 0) blue[x] = y; else red[x] = y;
            }
        }        
        
    return make_pair(blue, red);
}



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