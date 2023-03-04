vector<int> graph::min_degree()
{
    vector<int> res(realSize, -1);
    vector<int> degree(onSize + offSize, 0);
    vector<bool> matched(onSize + offSize, 0);
    vector<pair<int, int>> unmatchedNeighbors;
    int index;
    for (int i = 0; i < realSize; i++)
    {
        unmatchedNeighbors.clear();
        for (auto j : adj[types[i]])
            if (not matched[j])
                degree[j]++, unmatchedNeighbors.push_back(make_pair(degree[j], j));
        if (unmatchedNeighbors.size())
        {
            sort(unmatchedNeighbors.begin(), unmatchedNeighbors.end(), [](const pair<int, int> &p1, const pair<int, int> &p2)
                 { return (p1.first == p2.first) ? p1.second < p2.second : p1.first < p2.first; });
            index = unmatchedNeighbors[0].second;
            matched[index] = true;
            res[i] = index;
        }
    }
    return res;
}