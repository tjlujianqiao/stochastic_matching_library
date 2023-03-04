// Type graph stored in adjacency list representation, with online types and offline vertices
// Realization graph stored by type array of online vertices
class graph
{

public:
    // Construct an empty graph with n online types and m offline vertices
    graph(int n, int m)
    {
        onSize = n;
        offSize = m;
        
        adj.resize(n + m);
        for (int i = 0; i < onSize + offSize; i++)
            adj[i] = {};
        
        realSize = 0;
        types.clear();
    }
    
    // Add an edge (i, j)
    void add_edge(int i, int j)
    {
        adj[i].push_back(j);
        adj[j].push_back(i);
    }
    
    int online_size()
    {
        return onSize;
    }
    
    // Print graph
    void print()
    {
        cout << "Print type graph. Each line contains offline neighbors of an online vertex." << endl;
        for (int i = 0; i < onSize; i++)
        {
            cout << "Online vertex " << i << ":";
            for (int j : adj[i]) cout << " " << j;
            cout << endl;
        }
        cout << "Finished" << endl;
    }
    
    //Randomly construct a realization graph with n online vertices
    void realize(int n)
    {
        realSize = n;
        types.resize(realSize);
        uniform_int_distribution<int> typeDist(0, onSize - 1);
        for (int i = 0; i < realSize; i++)
            types[i] = typeDist(rng);
    }
    
    //Randomly construct a realization graph with n online vertices
    void print_type()
    {
        cout << "Type list:";
        for (int i : types) cout << " " << i;
        cout << endl;
    }
    
    
//NOTE: All following functions compute matchings in realization graph
    
    vector<int> maximum_matching();
    
    vector<int> sampling_without_replacement(map<pair<int, int>, double> &typeProb);
    map<pair<int, int>, double> optimal_matching_prob(int n_samples, int onSizeSample);

    vector<int> poisson_ocs(const vector<double> &offMass, map<pair<int, int>, double> &typeProb);
    vector<double> poisson_offline_mass(map<pair<int, int>, double> &typeProb);

    vector<int> top_half_sampling(map<pair<int, int>, double> &typeProb);
    
    pair<vector<int>, vector<int>> feldman_et_al_color();
    vector<int> feldman_et_al(vector<int> &blue, vector<int> &red);
    
    vector<int> manshadi_et_al(map<pair<int, int>, double> &typeProb);

    vector<vector<int>> jaillet_lu_list();
    vector<int> jaillet_lu(vector<vector<int>> &jlList);

    vector<int> min_degree();

    vector<int> ranking();

    vector<int> balance_swr();
    vector<int> balance_ocs();
    
private:
    // Adjacency list representation
    vector<vector<int>> adj;
    
    // Type of each online vertex in realization graph
    vector<int> types;
    
    // Number of online type
    int onSize;

    // Number of offline vertices
    int offSize;
    
    // Number of online vertices
    int realSize;
    
};
