// Compute LP in Brubach et al. (2016)

map<pair<int, int>, double> graph::brubach_et_al_lp()
{
    vector<int> ia, ja;
    vector<double> ar;
    vector<map<int, int>> eID(onSize);

    ia.push_back(-1); ja.push_back(-1); ar.push_back(-1.0);

    glp_prob *lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_term_out(GLP_OFF);

    int nRow = onSize + offSize;
    glp_add_rows(lp, nRow);
    for (int i = 1; i <= nRow; i++)
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, 1.0);

    for (int i = 0; i < onSize; i++)
        if (adj[i].size())
            glp_add_cols(lp, adj[i].size());
    
    int num = 0;
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
        {
            num++;
            eID[i][j] = num;
            
            ia.push_back(i + 1);
            ja.push_back(num);
            ar.push_back(1.0);

            ia.push_back(j + 1);
            ja.push_back(num);
            ar.push_back(1.0);
            
            glp_set_col_bnds(lp, num, GLP_DB, 0.0, 1.0 - exp(-1.0));
            glp_set_obj_coef(lp, num, 1.0);
        }

    for (int j = onSize; j < onSize + offSize; j++)
        for (int i1 : adj[j]) for (int i2 : adj[j])
            if (i1 < i2)
            {
                nRow++;
                glp_add_rows(lp, 1);
                glp_set_row_bnds(lp, nRow, GLP_UP, 0.0, 1.0 - exp(-2.0));
                
                ia.push_back(nRow);
                ja.push_back(eID[i1][j]);
                ar.push_back(1.0);

                ia.push_back(nRow);
                ja.push_back(eID[i2][j]);
                ar.push_back(1.0);
            }

    glp_load_matrix(lp, (int)ia.size() - 1, &ia[0], &ja[0], &ar[0]);

    glp_simplex(lp, NULL);

    map<pair<int, int>, double> res;

    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            res[make_pair(i, j)] = glp_get_col_prim(lp, eID[i][j]);

    glp_delete_prob(lp);
    return res;
}

// Compute H' in Brubach et al. (2016)
vector<vector<pair<int, double>>> graph::brubach_et_al_h(map<pair<int, int>, double> &lpSol)
{
    cycle_break_graph gCycle(onSize, onSize + offSize);
    for (int i = 0; i < onSize; i++)
        for (int j : adj[i])
            if (lpSol[make_pair(i, j)] > 1e-10)
                gCycle.add_edge(i, j, lpSol[make_pair(i,j)] * 3);

            
    gCycle.gandhi_et_al_rounding();
    gCycle.frac_to_int();
    gCycle.cycle_break();
    
    vector<int> offX(onSize + offSize, 0);
    for (int i = 0; i < onSize; i++)
        for (auto e : gCycle.vFrac[i])
            offX[e.first] += e.second;
    
    vector<vector<pair<int, double>>> res(onSize, vector<pair<int, double>>());
    
    for (int i = 0; i < onSize; i++)
        if (gCycle.vInt[i].size() == 1)
        {
            auto iter = gCycle.vInt[i].begin();
            
            int v1 = (*iter).first;
            res[i].push_back(make_pair(v1, 1));
            
        }
        else if (gCycle.vInt[i].size() == 2)
        {
            auto iter = gCycle.vInt[i].begin();
            map<int, int> &val = gCycle.vInt[i];
            
            int v1 = (*iter).first;
            iter++;
            int v2 = (*iter).first;
            
            if (offX[v1] > offX[v2]) swap(v1, v2);
            
            double x1 = val[v1] / 3.0, x2 = val[v2] / 3.0;
            
            if (offX[v1] == 1 && offX[v2] == 3 && val[v2] == 2)       // Case 1
                x1 = 0.1, x2 = 0.9;
            else if (offX[v1] == 2 && offX[v2] == 3 && val[v2] == 2)  // Case 2
                x1 = 0.15, x2 = 0.85;
            else if (offX[v1] == 2 && offX[v2] == 3 && val[v1] == 2)  // Case 3
                x1 = 0.6, x2 = 0.4;
            else if (offX[v1] == 1 && offX[v2] == 2)                  // Case 8
                x1 = 0.25, x2 = 0.75;
            else if (offX[v1] == 2 && offX[v2] == 2 && val[v1] == 1)  // Case 9
                x1 = 0.3, x2 = 0.7;
            else if (offX[v1] == 3 && offX[v2] == 3)
            {
                if (val[v1] < val[v2])
                {
                    swap(v1, v2);
                    x1 = val[v1] / 3.0, x2 = val[v2] / 3.0;
                }
                if (val[v1] == 2 && val[v2] == 1)
                {
                    if (gCycle.vInt[v2].size() == 2)                  // Case 10
                        x1 = 1 - 0.2744, x2 = 0.2744;
                    else                                              // Case 11
                        x1 = 1 - 0.15877, x2 = 0.15877;
                }
            }
            
            res[i].push_back(make_pair(v1, x1));
            res[i].push_back(make_pair(v2, x2));
            
        }
        else if (gCycle.vInt[i].size() == 3)
        {
            auto iter = gCycle.vInt[i].begin();
            map<int, int> &val = gCycle.vInt[i];
            
            int v1 = (*iter).first;
            *iter++;
            int v2 = (*iter).first;
            *iter++;
            int v3 = (*iter).first;
            
            if (offX[v1] > offX[v2]) swap(v1, v2);
            if (offX[v1] > offX[v3]) swap(v1, v3);
            if (offX[v2] > offX[v3]) swap(v2, v3);
            
            double x1 = val[v1] / 3.0, x2 = val[v2] / 3.0, x3 = val[v3] / 3.0;
            
            if (offX[v1] == 1 && offX[v2] == 3 && offX[v3] == 3)      // Case 4
                x1 = 0.1, x2 = 0.45, x3 = 0.45;
            else if (offX[v1] == 2 && offX[v2] == 3 && offX[v3] == 3) // Case 5
                x1 = 0.2, x2 = 0.4, x3 = 0.4;
            else if (offX[v1] == 1 && offX[v2] == 2 && offX[v3] == 3) // Case 6
                x1 = 0.15, x2 = 0.2, x3 = 0.65;
            else if (offX[v1] == 1 && offX[v2] == 1 && offX[v3] == 3) // Case 7
                x1 = 0.1, x2 = 0.1, x3 = 0.8;
            else if (offX[v1] == 2 && offX[v2] == 2 && offX[v3] == 3) // Case 12
                x1 = 0.25, x2 = 0.25, x3 = 0.5;
            
            
            res[i].push_back(make_pair(v1, x1));
            res[i].push_back(make_pair(v2, x2));
            res[i].push_back(make_pair(v3, x3));
            
        }
    return res;
}

// Match online vertices with weight H'
vector<int> graph::brubach_et_al(vector<vector<pair<int, double>>> h)
{
    vector<int> res(realSize, -1);
    vector<bool> matched(offSize + onSize, false);

    

    for (int i = 0; i < realSize; i++)
    {
        auto hI = h[types[i]];
        vector<int> list = {};
        if (hI.size() == 1)
        {
            int v1 = hI[0].first;
            list.push_back(v1);
        }
        else if (hI.size() == 2)
        {
            int v1 = hI[0].first, v2 = hI[1].first;
            double x1 = hI[0].second, x2 = hI[1].second;
            
            uniform_real_distribution<double> curRand(0.0, x1 + x2);
            if (curRand(rng) <= x1)
            {
                list.push_back(v1);
                list.push_back(v2);
            }
            else
            {
                list.push_back(v2);
                list.push_back(v1);
            }
        }
        else if (hI.size() == 3)
        {
            int v1 = hI[0].first, v2 = hI[1].first, v3 = hI[2].first;
            double x1 = hI[0].second, x2 = hI[1].second, x3 = hI[2].second;
            
            vector<pair<vector<int>, double>> sampleList = {};
            sampleList.push_back(make_pair(vector<int>{v1, v2, v3}, x1 * x2 / (x2 + x3)));
            sampleList.push_back(make_pair(vector<int>{v1, v3, v2}, x1 * x3 / (x3 + x2)));
            sampleList.push_back(make_pair(vector<int>{v2, v1, v3}, x2 * x1 / (x1 + x3)));
            sampleList.push_back(make_pair(vector<int>{v2, v3, v1}, x2 * x3 / (x3 + x1)));
            sampleList.push_back(make_pair(vector<int>{v3, v1, v2}, x3 * x1 / (x1 + x2)));
            sampleList.push_back(make_pair(vector<int>{v3, v2, v1}, x3 * x2 / (x2 + x1)));
            
            uniform_real_distribution<double> curRand(0.0, 1.0);
            double sample = curRand(rng), sum = 0;
            
            for (auto j : sampleList)
            {
                sum += j.second;
                if (sum >= sample)
                {
                    list = j.first;
                    break;
                }
            }
        }
        
        for (int j : list)
        {
            if (j != -1 && not matched[j])
            {
                res[i] = j;
                matched[j] = true;
                break;
            }
        }
    }
    return res;
}
