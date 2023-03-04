double fillWater(vector<double> &level, const double water)
{
    double l = 0, r = water, eps = 1e-8;
    while ((r - l) > eps)
    {
        double mid = (l + r) / 2, tt = 0;
        for (auto item : level)
            if (item < mid)
                tt += mid - item;
        if (tt >= water)
            r = mid;
        else
            l = mid;
    }
    return (l + r) / 2;
}

vector<int> graph::balance_SWR()
{
    vector<double> currentLevel(onSize + offSize, 0);
    vector<bool> selected(onSize + offSize, false);
    vector<int> res(realSize, -1);

    for (int i = 0; i < realSize; i++)
    {
        vector<double> level;

        for (int j : adj[types[i]])
            if (not selected[j] and currentLevel[j] < 1.0)
                level.push_back(currentLevel[j]);

        if (not level.size())
            continue;

        double newLevel = fillWater(level, 1);

        double mass = 0, chosen = 0;
        vector<pair<int, double>> validMass;
        for (int j : adj[types[i]])
            if (not selected[j] and currentLevel[j] < 1.0)
            {
                mass += max(newLevel - currentLevel[j], 0.0);
            }

        std::uniform_real_distribution<double> dist(0, mass);
        double sample = dist(rng);

        for (int j : adj[types[i]])
        {
            if (not selected[j] and currentLevel[j] < 1.0)
            {
                chosen += max(newLevel - currentLevel[j], 0.0);
                if (chosen >= sample)
                {
                    selected[j] = true;
                    res[i] = j;
                    break;
                }
            }
        }
        for (int j : adj[types[i]])
            if (not selected[j] and currentLevel[j] < 1.0)
                currentLevel[j] = max(newLevel, currentLevel[j]);
    }
    return res;
}

vector<int> graph::balance_OCS()
{
    vector<double> currentLevel(onSize + offSize, 0);
    vector<bool> selected(onSize + offSize, false);
    vector<int> res(realSize, -1);
    auto w = [](double y)
    {
        double c = (4 - 2 * sqrt(3)) / 3;
        return exp(1.0 * y + y * y / 2.0 + c * y * y * y);
    };
    for (int i = 0; i < realSize; i++)
    {
        vector<double> level;

        for (int j : adj[types[i]])
            if (not selected[j])
                level.push_back(currentLevel[j]);

        if (not level.size())
            continue;

        double newLevel = fillWater(level, 1);

        double mass = 0, chosen = 0;
        for (int j : adj[types[i]])
            if (not selected[j])
                mass += max((newLevel - currentLevel[j]), 0.0) * w(currentLevel[j]);

        std::uniform_real_distribution<double> dist(0, mass);
        double sample = dist(rng);

        for (int j : adj[types[i]])
        {
            if (not selected[j])
            {
                chosen += max(newLevel - currentLevel[j], 0.0) * w(currentLevel[j]);
                if (chosen >= sample)
                {
                    selected[j] = true;
                    res[i] = j;
                    break;
                }
            }
        }
        for (int j : adj[types[i]])
            if (not selected[j])
                currentLevel[j] = max(newLevel, currentLevel[j]);
    }
    return res;
}
