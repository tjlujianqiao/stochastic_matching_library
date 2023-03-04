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
        type.clear();
    }
    
    // Add an edge (i, j)
    void add_edge(int i, int j)
    {
        adj[i].push_back(j);
        adj[j].push_back(i);
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
        type.resize(realSize);
        uniform_int_distribution<int> typeDist(0, onSize - 1);
        for (int i = 0; i < realSize; i++)
            type[i] = typeDist(rng);
    }
    
    //Randomly construct a realization graph with n online vertices
    void print_type()
    {
        cout << "Type list:";
        for (int i : type) cout << " " << i;
        cout << endl;
    }
    
    
//NOTE: All following functions compute matchings in realization graph
    
    vector<int> maximum_matching();
    
private:
    // Adjacency list representation
    vector<vector<int>> adj;
    
    // Type of each online vertex in realization graph
    vector<int> type;
    
    // Number of online type
    int onSize;

    // Number of offline vertices
    int offSize;
    
    // Number of online vertices
    int realSize;
    
};
