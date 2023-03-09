graph generate_from_file(string path, bool dup = false, int subSample = 0)
{
    ifstream fin(path);

    // Ignore line 1
    string line;
    getline(fin, line);

    // Ignore the first character % in line 2
    char ch;
    int n, m;
    
    fin >> ch >> m >> n;

    if (not dup)
    {

        vector<int> id(n);
        iota(id.begin(), id.end(), 0);
        shuffle(id.begin(), id.end(), rng);
        if (subSample != 0)
            n = min(n, subSample);
        
        graph g(n / 2, n / 2);
        
        for (int i = 0; i < m; i++)
        {
            int x, y;
            fin >> x >> y;
            
            // Ignore the third number in line
            string line;
            getline(fin, line);

            // Index starts from 1 in input file
            x--, y--;

            if (id[x] < n / 2 && id[y] >= n / 2 && id[y] < n / 2 + n / 2)
                g.add_edge(id[x], id[y]);
        }
        fin.close();
        return g;
    }
    else
    {
        graph g(n, n);
        for (int i = 0; i < m; i++)
        {
            int x, y;
            fin >> x >> y;
            
            // Ignore the third number in line
            string line;
            getline(fin, line);

            // Index starts from 1 in input file
            x--, y--;

            g.add_edge(x, n + y);
        }
        fin.close();
        return g;
    }
}
