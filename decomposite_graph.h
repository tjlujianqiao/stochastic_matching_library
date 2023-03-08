//Input: A graph with only paths and cycles
//Decomposite into blue and red edges in Feldman et al. (2009) and Bahmani and Kapralov (2010)


struct decomposite_graph{
    
    // Color of each edge
    vector<map<int, int>> color;

    // Total number of vertices 
    int size;
    
    // Number of online types
    int onSize;

    // Construct an empty type graph with:
    // n vertices in total, and m online types 
    decomposite_graph(int n, int m)
    {
        size = n;
        onSize = m;
        
        color.resize(n);
        for (int i = 0; i < n; i++)
            color[i] = {};
    }
    
    // Add an edge (x, y)
    void add_edge(int x, int y)
    {
        color[x][y] = color[y][x] = -1;
    }
    
    // Assign colors for each edge
    void set_color()
    {
        vector<int> visit(size, 0);
        
        for (int i = 0; i < size; i++)
            if (visit[i] < 2)
            {
                // Start from any vertex to find the end of a path (or a cycle)
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
                
                
                // Start from one end of a path (or a cycle) to find the other end
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
                
                for (int j = 0; j < (int)vList.size() - 1; j++)
                {
                    int x = vList[j], y = vList[j + 1];
                    
                    if (isCycle)                    //Cycle
                        if (j % 2 == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                    else if (vList.size() % 2 == 0) //Odd-length paths
                        if (j % 2 == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                    else if (vList[0] >= onSize)    //Even-length paths starting from an offline vertex
                        if (j % 2 == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                    else                            //Even-length paths starting from an online vertex
                        if (j % 2 == 1 || j == 0) color[x][y] = color[y][x] = 1; else color[x][y] = color[y][x] = 2;
                }
            }
    }
};
