//Compute maximum flow by Dinic's algorithm
//Assume t is the vertex with largest label 

struct decomposite_graph{
    
    vector<map<int, int>> color;

    int size;
    int onSize;

    decomposite_graph(int n, int m)
    {
        size = n;
        onSize = m;
        
        color.resize(n);
        for (int i = 0; i < n; i++)
            color[i] = {};
    }
    
    void add_edge(int x, int y)
    {
        color[x][y] = color[y][x] = -1;
    }
    
    void set_color()
    {
        vector<int> visit(size, 0);
        
        for (int i = 0; i < size; i++)
            if (visit[i] < 2)
            {
                int cur = i, pre = -1, flag;
                visit[i] = 1;
                do
                {
                    flag = false;
                    for (auto e : color[cur])
                        if (e.first != pre)
                        {
                            int next = e.first;
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
                vector<int> vList(1, start);
                visit[start] = 2;
                pre = -1;
                bool isCycle = false;
                
                do
                {
                    flag = false;
                    for (auto e : color[cur])
                        if (e.first != pre)
                        {
                            int next = e.first;
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
                    
                    if (isCycle) //Cycle
                        if (j % 2 == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                    else if (vList.size() % 2 == 0) //Odd-length paths
                        if (j % 2 == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                    else if (vList[0] >= onSize) //Even-length paths starting from an offline vertex
                        if (j % 2 == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                    else                         //Even-length paths starting from an online vertex
                        if (j % 2 == 1 || j == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                }
            }
    }
};
