// Match by min-degree algorithm
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
            // By default the sort function sorts the vector elements on basis of first element of pairs.
            sort(unmatchedNeighbors.begin(), unmatchedNeighbors.end());
            index = unmatchedNeighbors[0].second;
            matched[index] = true;
            res[i] = index;
        }
    }
    return res;
}