map<pair<int, int>, double> graph::brubach_et_al_lp()
{
    vector<int> ia, ja;
    vector<double> ea;
    map<pair<int, int>, int> varnames;
    vector<pair<int, int> > vars;

    ia.push_back(-1); ja.push_back(-1); ea.push_back(-1.0);
    vars.push_back(make_pair(-1, -1));

    glp_prob *lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);

    glp_term_out(GLP_OFF);

    //cerr << "************ SET ROW BOUNDS" << endl;

    int nrows = onSize + offSize;
    glp_add_rows(lp, nrows);
    for (int i = 0; i < nrows; i++)
        glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, 1.0);

    //cerr << "************ SET COLUMN BOUNDS (1-1/e) AND OBJECTIVE" << endl;

    int ncols = 0;
    for (int i = 0; i < onSize; i++)
        ncols += (int)adj[i].size();
    
    
    glp_add_cols(lp, ncols);
    for (int i = 0; i < ncols; i++) {
        glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0 - 1.0 / exp(1.0));
        glp_set_obj_coef(lp, i + 1, 1.0);
    }

    //cerr << "************ SET UP ADJACENCY MATRIX CONSTRAINTS" << endl;

    int var = 1;
    for (int i = 0; i < onSize; i++) {
        for (int j = 0; j < (int)adj[i].size(); j++) {
            ia.push_back(i + 1);
            ja.push_back(var);
            ea.push_back(1.0);

            ia.push_back(adj[i][j] + 1);
            ja.push_back(var);
            ea.push_back(1.0);

            varnames[make_pair(i, adj[i][j])] = var;
            vars.push_back(make_pair(i, adj[i][j]));
            var++;
        }
    }

    //cerr << "***************** SET UP 1-1/e^2 CONSTRAINTS" << endl;

    int currow = nrows;
    for (int j = onSize; j < (int)onSize + offSize; j++) {
        for (int i1 = 0; i1 < (int)adj[j].size(); i1++) {
            for (int i2 = i1 + 1; i2 < (int)adj[j].size(); i2++) {
                pair<int, int> e1 = make_pair(adj[j][i1], j);
                pair<int, int> e2 = make_pair(adj[j][i2], j);
                glp_add_rows(lp, 1);
                currow++;
                glp_set_row_bnds(lp, currow, GLP_UP, 0.0, 1.0 - 1.0 / 
                    (exp(1.0)*exp(1.0)));
                
                ia.push_back(currow);
                ja.push_back(varnames[e1]);
                ea.push_back(1.0);

                ia.push_back(currow);
                ja.push_back(varnames[e2]);
                ea.push_back(1.0);
            }
        }
    }

    //cerr << "************ LOAD MATRIX OF CONSTRAINTS" << endl;

    glp_load_matrix(lp, (int)ia.size() - 1, &ia[0], &ja[0], &ea[0]);

    //cerr << "************ RUN SIMPLEX" << endl;
    
    glp_simplex(lp, NULL);

    //cerr << "************ DONE" << endl;
    
    map<pair<int, int>, double> res;

    for (int i = 1; i <= ncols; i++) {
        double val = glp_get_col_prim(lp, i);
        res[vars[i]] = val;
    }

    glp_delete_prob(lp);
    glp_free_env();
    return res;
}


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