// Compute natural LP in Huang, Shu (2021) by ellipsoid method
// Extremely slow, only works on very small graph
// Reference: https://web.stanford.edu/class/ee364b/lectures/ellipsoid_method_notes.pdf

class natural_lp
{
private:
    map<pair<int, int>, int> eID;
    vector<vector<int>> adjLP;
    const double eps_obj = 1e-3;
    const double eps_feas = 1e-3;
    double f_best = 0;
    int onSize, n;
    double *x;
    double *g_k;
    double **P;
    double *lambda;

public:
    // Initialize
    natural_lp(const vector<vector<int>> &adj, int onsize)
    {
        adjLP = adj;
        onSize = onsize;
        n = setid();
    };

    // Assign ids to each edge
    int setid()
    {
        int id = 1;
        for (int i = 0; i < onSize; i++)
            for (int j : adjLP[i])
            {
                eID[make_pair(i, j)] = id;
                id++;
            }
        return --id;
    }

    map<pair<int, int>, double> solve_lp()
    {

        map<pair<int, int>, double> typeProb;
        if (n == 1)
        {
            for (int i = 0; i < onSize; i++)
                for (int j : adjLP[i])
                    typeProb[make_pair(i, j)] = 1 - exp(-1);
            return typeProb;
        }
        // Initial point (0, 0, ..., 0) is feasible
        P = (double **)malloc(sizeof(double *) * (n + 1));
        g_k = (double *)malloc(sizeof(double) * (n + 1));
        x = (double *)malloc(sizeof(double) * (n + 1));
        lambda = (double *)malloc(sizeof(double) * (onSize));
        

        for (int i = 0; i < n + 1; i++)
        {
            *(P + i) = (double *)malloc(sizeof(double) * (n + 1));
        }

        for (int i = 0; i < n + 1; i++)
            for (int j = 0; j < n + 1; j++)
                P[i][j] = 0;

        for (int i = 1; i < n + 1; i++)
        {
            P[i][i] = n / 4.0;
            g_k[i] = 0;
            x[i] = 0.5;
        }
        for (int i = 0; i < onSize; i++)
            lambda[i] = 1;

        f_best = 0;
        while ( not iterate_ellipsoid() );

        for (int i = 0; i < onSize; i++)
            for (int j : adjLP[i])
                typeProb[make_pair(i, j)] = x[eID[make_pair(i, j)]];

        free(P);
        P = nullptr;
        free(g_k);
        g_k = nullptr;
        free(x);
        x = nullptr;

        return typeProb;
    }

    // Print the optimal solution
    void print()
    {
        for (int i = 0; i < onSize; i++)
        {
            for (int j : adjLP[i])
                cout << "X[" << i << ", " << j << "] = " << x[eID[make_pair(i, j)]] << endl;
        }
        return;
    }

    // Iterating in ellipsoid method
    bool iterate_ellipsoid()
    {
        double stop_value, sum;

        fill(g_k, g_k + n + 1, 0.0);

        for (int i = 0; i < onSize; i++)
        {
            // Check constraint of each online type
            sum = -lambda[i];
            for (int j : adjLP[i])
            {
                int id = eID[make_pair(i, j)];
                sum += x[id];
            }

            if (sum > eps_feas)
            {
                for (int j : adjLP[i])
                {
                    int id = eID[make_pair(i, j)];
                    g_k[id] = (double)1.0;
                }
                stop_value = cal_stop_value();
                if (stop_value < sum)
                {
                    return true;
                }
                update_ellipsoid(sum, stop_value);
                return false;
            }
        }

        // Check the natural constraint of each offline type
        for (int j = onSize; j < (int)adjLP.size(); j++)
        {
            vector<pair<double, int>> x_j;
            for (int i : adjLP[j])
            {
                int id = eID[make_pair(i, j)];
                x_j.push_back(make_pair(x[id], i));
            }
            sort(x_j.begin(), x_j.end(), [](const pair<double, int> &p1, const pair<double, int> &p2)
                 { return p1.first > p2.first; });

            double sum_l = 0, sum_x = 0;

            for (int s = 0; s < (int)x_j.size(); s++)
            {
                sum_l += lambda[x_j[s].second];
                sum_x += x_j[s].first;
                sum = sum_x + expl(-sum_l) - 1;
                if (sum > eps_feas)
                {
                    for (int e = 0; e <= s; e++)
                    {
                        int id = eID[make_pair(x_j[s].second, j)];
                        g_k[id] = (double)1.0;
                    }
                    stop_value = cal_stop_value();

                    if (stop_value < sum)
                    {
                        return true;
                    }
                    update_ellipsoid(sum, stop_value);
                    return false;
                }
            }
        }

        // Check constraint x_{ij} >= 0 of each edge
        for (int i = 0; i < onSize; i++)
        {
            for (int j : adjLP[i])
            {
                int id = eID[make_pair(i, j)];
                sum = -x[id];
                if (sum > eps_feas)
                {
                    g_k[id] = (double) -1.0;
                    stop_value = cal_stop_value();
                    if (stop_value < sum)
                    {
                        return true;
                    }
                    update_ellipsoid(sum, stop_value);
                    return false;
                }
            }
        }

        for (int i = 0; i < n + 1; i++) g_k[i] = (double) -1.0;
        stop_value = cal_stop_value();
        if (stop_value <= eps_obj)
        {
            return true;
        }
        sum = get_obj();
        update_ellipsoid(sum, stop_value, false);
        return false;
    }

    // Compute stop value
    double cal_stop_value()
    {
        double stop_value = 0;
        for (int i = 1; i < n + 1; i++)
            for (int j = 1; j < n + 1; j++)
                stop_value += P[i][j] * g_k[i] * g_k[j];

        return sqrt(stop_value);
    }

    // Compute current value of objective function
    double get_obj()
    {
        double fv = 0;
        for (int i = 1; i < n + 1; i++)
            fv -= x[i];
        return fv;
    }

    // Update x and P in one iteration
    void update_ellipsoid(double sum, double stop_value, bool flag = true)
    {
        double g_n[n + 1];
        for (int i = 0; i < n + 1; i++)
        {
            g_n[i] = g_k[i];
        }

        for (int i = 1; i < n + 1; i++)
            g_n[i] /= stop_value;

        double alpha = 0;

        if (flag)
        {
            alpha = 1.0 * sum / stop_value;
        }
        else
        {
            f_best = min(f_best, sum);
            alpha = 1.0 * (sum - f_best) / stop_value;
        }

        double term1[n + 1];
        double term2[n + 1];
        for (int i = 1; i < n +1 ; i++)
        term1[i] = term2[i] = 0;
        

        for (int i = 1; i < n + 1; i++)
            for (int j = 1; j < n + 1; j++)
            {
                term1[i] += P[i][j] * g_n[j];
                term2[j] += g_n[i] * P[i][j];
            }
        for (int i = 1; i < n + 1; i++)
        {
            x[i] -= (double)(1.0 + n * alpha) / (n + 1) * term1[i];
            for (int j = 1; j < n + 1; j++)
            {
                P[i][j] = 1.0 * n * n *(1.0 - alpha * alpha) / (n * n - 1) *
                          (P[i][j] - 2.0 * (1.0 + n * alpha) / ((n + 1) * (1 + alpha)) * term1[i] * term2[j]);
            }
        }
        
    }
};
