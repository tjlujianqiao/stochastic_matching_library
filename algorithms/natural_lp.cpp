
class natural_lp
{
private:
    map<pair<int, int>, int> eID;
    vector<vector<int>> adjLP;
    const double eps = 0.001;
    double f_best = 0;
    int onSize, n;
    vector<double> x, g_k, lambda;
    vector<vector<double>> P;

public:
    natural_lp(const vector<vector<int>> &adj, int onsize)
    {
        adjLP = adj;
        onSize = onsize;
        lambda = vector<double>(onsize, 1);
        n = setid();
    };
    int setid()
    {
        int id = 1;
        for (int i = 0; i < onSize; i++)
        {
            for (int j : adjLP[i])
            {
                eID[make_pair(i, j)] = id;
                id++;
            }
        }
        return --id;
    }
    map<pair<int, int>, double> solve_lp()
    {
        vector<double> x_init(n + 1, 0);
        vector<vector<double>> P_init(n + 1, vector<double>(n + 1, 0));
        lambda = vector<double>(onSize, 1);
        for (int i = 1; i <= n; i++)
        {
            P_init[i][i] = 1.0 * n;
        }
        x = x_init;
        P = P_init;
        // the initial point is feasible
        f_best = get_obj();
        while (not iterate_ellipsoid())
            ;
        map<pair<int, int>, double> typeProb;
        for (int i = 0; i < onSize; i++)
        {
            for (int j : adjLP[i])
            {
                typeProb[make_pair(i, j)] = x[eID[make_pair(i, j)]];
            }
        }
        return typeProb;
    }
    void print()
    {
        for (int i = 0; i < onSize; i++)
        {
            for (int j : adjLP[i])
                printf("i:%d,j:%d, %f \n", i, j, x[eID[make_pair(i, j)]]);
        }
        return;
    }

    bool iterate_ellipsoid()
    {
        double stop_value, sum;
        g_k = vector<double>(x.size(), 0);

        for (int i = 0; i < onSize; i++)
        {
            sum = -lambda[i];
            for (int j : adjLP[i])
            {
                int id = eID[make_pair(i, j)];
                sum += x[id];
            }

            if (sum > 0)
            {
                for (int j : adjLP[i])
                {
                    int id = eID[make_pair(i, j)];
                    g_k[id] = 1;
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
        for (int j = onSize; j < adjLP.size(); j++)
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

            for (int s = 0; s < x_j.size(); s++)
            {
                sum_l += lambda[x_j[s].second];
                sum_x += x_j[s].first;
                int id = eID[make_pair(x_j[s].second, j)];
                sum = sum_x + exp(-sum_l) - 1;
                if (sum > 0)
                {
                    for (int e = 0; e <= s; e++)
                    {
                        int id = eID[make_pair(x_j[s].second, j)];
                        g_k[id] = 1;
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
        for (int i = 0; i < onSize; i++)
        {
            for (int j : adjLP[i])
            {
                int id = eID[make_pair(i, j)];
                sum = -x[id];
                if (sum > 0)
                {
                    g_k[id] = -1;
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

        fill(g_k.begin(), g_k.end(), -1);
        stop_value = cal_stop_value();
        sum = get_obj();
        update_ellipsoid(sum, stop_value, false);
        if (stop_value <= eps)
        {
            return true;
        }
        return false;
    }
    double cal_stop_value()
    {
        double stop_value = 0;
        for (int i = 1; i < P.size(); i++)
        {
            for (int j = 1; j < P[i].size(); j++)
            {
                stop_value += P[i][j] * g_k[i] * g_k[j];
            }
        }
        // cout << stop_value << endl;
        // assert(stop_value > 0);
        return sqrt(stop_value);
    }
    double get_obj()
    {
        double fv = 0;
        for (int i = 1; i < x.size(); i++)
            fv -= x[i];
        return fv;
    }
    void update_ellipsoid(double sum, double stop_value, bool flag = true)
    {
        vector<double> g_n = g_k;
        for (auto &e : g_n)
        {
            e /= stop_value;
        }
        double alpha, n = x.size();
        if (flag)
        {
            alpha = 1.0 * sum / stop_value;
        }
        else
        {
            f_best = min(f_best, get_obj());
            alpha = 1.0 * (sum - f_best) / stop_value;
        }
        vector<double> x_n(x.begin(), x.end());
        vector<vector<double>> P_n(n, vector<double>(n, 0));
        vector<double> midterm1(n, 0);
        vector<double> midterm2(n, 0);
        vector<vector<double>> midterm3(n, vector<double>(n, 0));
        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < n; j++)
            {
                midterm1[i] += P[i][j] * g_n[j];
                midterm2[i] += g_n[j] * P[j][i];
            }
        }

        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < n; j++)
            {
                midterm3[i][j] = midterm1[i] * midterm2[j];
            }
        }
        for (int i = 1; i < n; i++)
        {
            x_n[i] -= (1.0 + n * alpha) / (n + 1) * midterm1[i];
            for (int j = 1; j < n; j++)
            {
                P_n[i][j] = 1.0 * n * n / (n * n - 1) *
                            (1 - alpha * alpha) *
                            (P[i][j] - 2.0 * (1 + n * alpha) / ((n + 1) * (1 + alpha)) * midterm3[i][j]);
            }
        }
        x = x_n;
        P = P_n;
    }
};
