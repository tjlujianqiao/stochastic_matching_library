// Match by RANKING algorithm
vector<int> graph::ranking()
{
    vector<int> rank(onSize + offSize, 0);
    vector<int> res(realSize, -1);
    vector<int> matched(onSize + offSize, -1);

    iota(rank.begin(), rank.end(), 0);
    shuffle(rank.begin(), rank.end(), rng);

    for (int i = 0; i < realSize; i++)
    {
        int match = -1;
        
        for (int j : adj[types[i]])
            if (matched[j] == -1)
                if (match == -1 || rank[j] > rank[match])
                    match = j;

        res[i] = match;
        if (match != -1)
            matched[match] = i;
        
    }
    return res;
}

