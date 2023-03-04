graph generate_from_file(string path)
{
    ifstream fin(path);
    
    //Ignore line 1
    string line;
    getline(fin, line);
    
    //Ignore the first character % in line 2
    char ch;
    int n, m;
    fin >> ch >> m >> n;
    
    graph g(n/2, n/2);
    
    vector<int> id(n);
    iota(id.begin(), id.end(), 0);
    shuffle(id.begin(), id.end(), rng);
    
    for (int i = 0; i < m; i++)
    {
        int x, y;
        fin >> x >> y;
        
        //Index starts from 1 in input file
        x--, y--;
        
        //if (id[x] > id[y])
            //swap(x, y);
        
        if (id[x] < n/2 && id[y] >= n/2 && id[y] < n/2 + n/2)
            g.add_edge(id[x], id[y]);
    }
    
    fin.close();
    
    return g;
}