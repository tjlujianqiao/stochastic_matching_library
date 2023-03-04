vector<int> graph::ranking()
{
    vector<int> rank(realSize, 0);
    vector<int> res(realSize, -1);
    vector<int> matched(onSize + offSize, -1);

    for (int j = 0; j < onSize + offSize; j++)
        rank[j] = j;

    shuffle(rank.begin() + onSize, rank.end(), rng);

    vector<int> sortedNb;

    for (int i = 0; i < realSize; i++)
    {
        sortedNb = adj[types[i]];
        sort(sortedNb.begin(), sortedNb.end(), [&rank](const int& j1, const int & j2) {return rank[j1] > rank[j2];});
        for (int j : sortedNb)
        {
            if (matched[j] == -1)
            {
                res[i] = j;
                matched[j] = i;
                break;
            }
        }
    }
    return res;
}

