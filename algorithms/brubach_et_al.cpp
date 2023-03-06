vector<vector<pair<int, double>>> graph::brubach_et_al_h(map<pair<int, int>, double> lpSol)
{
    cycle_break_graph gCycle(onSize, onSize + offSize);
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            if (lpSol[make_pair(i, j)] > 1e-10)
                gCycle.add_edge(i, j, lpSol[make_pair(i,j)]);

            
    gCycle.gandhi_et_al_rounding();
    gCycle.frac_to_int();
    gCycle.cycle_break();
    
    vector<int> onX(onSize + offSize, 0), offX(onSize + offSize, 0);
    for (int i = 0; i < onSize; i++)
        for (auto e : gCycle.vInt[i])
            onX[i] += e.second, offX[e.first] += e.second;
        
    for (int i = 0; i < onSize + offSize; i++)
        if(onX[i] > 3 || offX[i] > 3)cout<<"WRONG "<<onX[i]<<" "<<offX[i]<<endl;
    
    
    vector<vector<pair<int, double>>> res(onSize, vector<pair<int, double>>());
    
    for (int i = 0; i < onSize; i++) 
        for (auto e : gCycle.vInt[i])
                res[i].push_back(make_pair(e.first, e.second));
    return res;
}


vector<int> graph::brubach_et_al(vector<vector<pair<int, double>>> h)
{
    vector<int> res(realSize, -1);
    vector<bool> matched(offSize + onSize, false);

    for (int i = 0; i < realSize; i++)
    {
        /*shuffle(jlList[types[i]].begin(), jlList[types[i]].end(), rng);
        for (int j : jlList[types[i]])
        {
            if (j != -1 && not matched[j])
            {
                res[i] = j;
                matched[j] = true;
                break;
            }
        }*/
    }
    return res;
}