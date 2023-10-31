#pragma GCC target("avx")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

#include <vector>
#include <set>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include "vec_op.hpp"
#include <eigen-3.3.8/Eigen/Dense>

using namespace std;
using namespace Eigen;

#define INF 1e12
struct edge
{
    long long to, cost;
};
typedef pair<long long, long long> P;
struct graph
{
    long long V;
    vector<vector<edge>> G;
    vector<long long> d;

    graph(long long n)
    {
        init(n);
    }

    void init(long long n)
    {
        V = n;
        G.resize(V);
        d.resize(V);
        for (int i = 0; i < V; i++)
        {
            d[i] = INF;
        }
    }

    void add_edge(long long s, long long t, long long cost)
    {
        edge e;
        e.to = t, e.cost = cost;
        G[s].push_back(e);
    }

    void dijkstra(long long s)
    {
        for (int i = 0; i < V; i++)
        {
            d[i] = INF;
        }
        d[s] = 0;
        priority_queue<P, vector<P>, greater<P>> que;
        que.push(P(0, s));
        while (!que.empty())
        {
            P p = que.top();
            que.pop();
            long long v = p.second;
            if (d[v] < p.first)
                continue;
            for (auto e : G[v])
            {
                if (d[e.to] > d[v] + e.cost)
                {
                    d[e.to] = d[v] + e.cost;
                    que.push(P(d[e.to], e.to));
                }
            }
        }
    }
};

struct gr_system
{
public:
    using state = vector<double>;
    int genes;
    vector<vector<double>> GRN;
    vector<double> theta;
    int num_final_attractors = 0;
    int depth = 0;
    vector<int> diff_list{};
    vector<vector<gr_system::state>> attractors = vector<vector<gr_system::state>>(1, vector<gr_system::state>(0));
    vector<vector<int>> goes_to = vector<vector<int>>(1, vector<int>(1, 0));

public:
    gr_system(int genes_, vector<vector<double>> GRN_, vector<double> theta_)
        : genes(genes_), GRN(GRN_), theta(theta_) {}
};

bool compare_by_fitness(gr_system idv_a, gr_system idv_b)
{
    return idv_a.num_final_attractors > idv_b.num_final_attractors;
}

double dot(const vector<double> &a, const vector<double> &b)
{
    double result = 0;
    for (int i = 0; i < a.size(); i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

int move(int genes, const gr_system &System, gr_system::state &current_state, vector<double> &dxdt, gr_system::state &destination, double &c, int &c_id)
{
    double infinitesimal = 1e-12;
    double c_;

    for (int i = 0; i < genes; i++)
    {
        c_ = (System.theta[i] - dot(System.GRN[i], current_state)) / dot(System.GRN[i], dxdt);
        if (c_ > 0 && c_ < c)
        {
            c = c_;
            c_id = i;
        }
    }
    if (c > 1)
    {
        current_state = destination;
        return 1;
    }
    else
    {
        current_state += (c + infinitesimal) * dxdt;
        // destination[c_id] = !destination[c_id];
        for (int i = 0; i < genes; i++)
        {
            destination[i] = (double)(dot(current_state, System.GRN[i]) - System.theta[i] > 0);
        }
        dxdt = destination - current_state;
        return 0;
    }
}

void check(const int genes, gr_system::state &current_state, bool &is_fp, const gr_system &System)
{
    int count = 0, itr = 0; // lc 判定（こんなんで良いのか？？）
    // gr_system::state current_state(genes), destination(genes);
    gr_system::state destination(genes);
    // current_state = initial_state;
    for (int i = 0; i < genes; i++)
    {
        destination[i] = (double)(dot(current_state, System.GRN[i]) - System.theta[i] > 0);
    }
    gr_system::state former_state = current_state, dx = former_state - current_state;
    vector<double> dxdt = destination - current_state;
    while (true)
    {
        former_state = current_state;
        int c_id = -1;
        double c = 1e12;
        if (move(genes, System, current_state, dxdt, destination, c, c_id) == 1)
        {
            break;
        };
        dx = former_state - current_state;
        itr++;
        if (itr > 5 && all_of(dx.begin(), dx.end(), [](double dx_i)
                              { return abs(dx_i) < 1e-4; }))
        {
            itr = 0;
            vector<int> area(genes, 0), area_init(genes), difference(genes, 1);
            for (int i = 0; i < genes; i++)
            {
                if (dot(System.GRN[i], current_state) - System.theta[i] > 0)
                {
                    area[i] = 1; // is plus
                }
            }
            area_init = area;
            set<int> which_plane_list{};
            bool goes_to_lattice = false, error = false;
            while (!all_of(difference.begin(), difference.end(), [](bool difference_i)
                           { return difference_i == 0; }))
            {
                if (itr % 100 == 0)
                {
                    area_init = area; // 定義し直し
                    which_plane_list = {};
                }
                double c_add = 1e12;
                int c_add_id = -1;
                if (move(genes, System, current_state, dxdt, destination, c_add, c_add_id) == 1)
                {
                    goes_to_lattice = true;
                    break;
                }
                area[c_add_id] = !area[c_add_id];
                which_plane_list.insert(c_add_id);
                difference = area - area_init;
                itr++;
                if (itr > 1000)
                {
                    error = true;
                    break;
                }
            }
            if (goes_to_lattice)
            {
                break;
            }
            else if (error)
            {
                is_fp = false; // リミットサイクル関係ないけど、強制終了として
                break;
            }
            // calculate_hyperplane_intersection;
            MatrixXd A(genes, genes);
            VectorXd b(genes);
            for (int i = 0; i < genes; i++)
            {
                if (which_plane_list.find(i) != which_plane_list.end())
                {
                    for (int j = 0; j < genes; j++)
                    {
                        A(i, j) = System.GRN[i][j];
                    }
                    b(i) = System.theta[i];
                }
                else
                {
                    for (int j = 0; j < genes; j++)
                    {
                        A(i, j) = (i == j) ? 1. : 0.;
                    }
                    b(i) = area[i];
                }
            }

            VectorXd x = A.colPivHouseholderQr().solve(b);
            // move()
            double c_eff = 1e12;
            int c_eff_id = -1;
            for (int i = 0; i < genes; i++)
            {
                if (which_plane_list.find(i) != which_plane_list.end())
                {
                    destination[i] = x(i);
                }
                else
                {
                    destination[i] = area[i];
                }
            }
            dxdt = destination - current_state;

            if (move(genes, System, current_state, dxdt, destination, c_eff, c_eff_id) == 1)
            {
                break;
            }
            if (which_plane_list.find(c_eff_id) != which_plane_list.end())
            {
                for (int i = 0; i < genes; i++)
                {
                    current_state[i] = x(i);
                }
                break;
            }
            itr = 0;
        }
        count++;
        if (count > 100)
        {
            is_fp = false;
            break;
        }
    }
}

void check_atr_property_step(const int genes, int mother, gr_system &System, gr_system::state &state, bool &is_fp, bool &is_already, bool &has_lc)
{
    is_fp = true;
    int count = 0, itr = 0; // lc 判定（こんなんで良いのか？？）
    gr_system::state current_state = state;
    check(genes, current_state, is_fp, System);

    if (is_fp)
    {
        for (int a = 0; a < System.attractors.size(); a++)
        {
            bool is_same = true;
            for (int i = 0; i < genes; i++)
            {
                if (abs(current_state[i] - System.attractors[a][0][i]) < 1e-2)
                {
                    continue;
                }
                else
                {
                    is_same = false;
                    break;
                }
            }
            if (is_same)
            {
                System.goes_to[mother][a] += 1;
                is_already = true;
                break;
            }
        }
        if (is_already)
        {
            return;
        }
        else // new attractor
        {
            vector<gr_system::state> v(1, current_state);
            System.attractors.push_back(v);
            System.goes_to.push_back(vector<int>(System.attractors.size() - 1, 0));
            for (int a = 0; a < System.attractors.size(); a++)
            {
                System.goes_to[a].push_back(0);
            }
            int id = System.attractors.size() - 1;
            System.goes_to[mother][id] += 1;
        }
    }
    else
    {
        has_lc = true;
    }
}

auto is_minus = [](double x)
{
    return (x < 0);
};
auto is_greater_than_one = [](double x)
{
    return (x > 1);
};

void dfs_step(const int genes, const vector<gr_system::state> &noise_arr, int mother, gr_system &System, bool &has_lc)
{
    for (auto noise : noise_arr)
    {
        gr_system::state noise_added_state(genes);
        for (int i = 0; i < genes; i++)
        {
            noise_added_state[i] = System.attractors[mother][0][i] + noise[i];
            if (noise_added_state[i] > 1)
            {
                noise_added_state[i] -= 2 * noise[i];
            }
            else if (noise_added_state[i] < 0)
            {
                noise_added_state[i] -= 2 * noise[i];
            }
        }
        bool is_fp = true, is_already = false;
        check_atr_property_step(genes, mother, System, noise_added_state, is_fp, is_already, has_lc);
        if (has_lc)
        {
            return;
        }
        else if (!is_already)
        {
            int id = System.attractors.size() - 1;
            dfs_step(genes, noise_arr, id, System, has_lc);
        }
    }
}

void calculate_fitness_step(const int genes, vector<gr_system::state> &noise_arr, const gr_system::state &initial_state, gr_system &System, const int num_noise_type)
{
    gr_system::state current_state = initial_state;
    bool is_fp = true;
    check(genes, current_state, is_fp, System);

    if (is_fp)
    {
        System.attractors[0].push_back(current_state);
    }
    else
    {
        System.num_final_attractors = 0;
        return;
    }

    bool has_lc = false;
    dfs_step(genes, noise_arr, 0, System, has_lc);
    if (has_lc)
    {
        System.num_final_attractors = 0;
        System.depth = 0;
    }
    else
    {
        vector<set<int>> st(System.goes_to.size(), set<int>{});
        for (int atr = 0; atr < System.goes_to.size(); atr++)
        {
            for (int ch = 0; ch < System.goes_to[atr].size(); ch++)
            {
                st[atr].insert(System.goes_to[atr][ch]);
            }
        }
        int differentiated = 0;
        for (int i = 0; i < System.goes_to.size(); i++)
        {
            if (System.goes_to[i][i] == num_noise_type)
            {
                differentiated++;
                System.diff_list.push_back(i);
            }
        }
        System.num_final_attractors = differentiated;

        graph g(System.goes_to.size());
        for (int atr = 0; atr < System.goes_to.size(); atr++)
        {
            for (int dest = 0; dest < System.goes_to[atr].size(); dest++)
            {
                if (System.goes_to[atr][dest] != 0)
                {
                    g.add_edge(atr, dest, 1);
                }
            }
        }
        g.dijkstra(0);
        for (auto diff_type : System.diff_list)
        {
            if (g.d[diff_type] > System.depth)
            {
                System.depth = g.d[diff_type];
            }
        }
    }
}

int powint(int a, int b)
{
    int res = 1;
    for (int i = 0; i < b; i++)
    {
        res *= a;
    }
    return res;
}

int main(int argc, const char **argv)
{
    int genes = stoi(argv[1]);
    int num_noise_type = stoi(argv[2]);
    int seed = stoi(argv[3]);
    double delta = stod(argv[4]);
    string str_delta = argv[4];
    mt19937_64 mt(seed);

    string dir_name = "../data_random_noise" + str_delta + "/genes_" + to_string(genes) + "_step/num_noise_type_" + to_string(num_noise_type) + "/seed_" + to_string(seed);
    ofstream sysinfo_writing_file(dir_name + "/fitness.txt");
    ofstream GRN_writing_file(dir_name + "/GRN.txt");
    ofstream theta_writing_file(dir_name + "/theta.txt");
    ofstream depth_writing_file(dir_name + "/depth.txt");
    ofstream differentiate_writing_file(dir_name + "/diffenrentiate_network.txt");
    ofstream attractors_writing_file(dir_name + "/attractors.txt");
    ofstream goes_to_writing_file(dir_name + "/atr_adj_mat.txt");

    ifstream N_noise_read("../perturbation_lib/N_noise_genes" + to_string(genes) + "_delta_" + str_delta + "_seed_" + to_string(seed) + "_noise_set_N=" + to_string(num_noise_type) + "_including_e.txt");
    vector<gr_system::state> noise_arr(num_noise_type, gr_system::state(genes, 0));
    for (int type = 0; type < num_noise_type; type++)
    {
        for (int i = 0; i < genes; i++)
        {
            N_noise_read >> noise_arr[type][i];
        }
    }

    uniform_real_distribution<double> rand_for_theta(-1, 1);
    vector<vector<double>> GRN(genes, vector<double>(genes));
    vector<double> theta(genes);
    double mutation_rate = 0.05;
    normal_distribution<double> mutation(0, 0.5);

    for (int i = 0; i < genes; i++)
    {
        for (int j = 0; j < genes; j++)
        {
            GRN[i][j] = rand_for_theta(mt);
        }
    }
    for (int i = 0; i < genes; i++)
    {
        theta[i] = rand_for_theta(mt);
    }

    gr_system System{genes, GRN, theta};
    gr_system::state initial_state(genes);
    for (int i = 0; i < genes; i++)
    {
        initial_state[i] = 0;
    }

    int max_gen = 1e6, population = 100;
    uniform_real_distribution<double> rand(0, 1);

    vector<gr_system> idv_vec(population, gr_system{genes, GRN, theta});
    for (auto idv : idv_vec)
    {
        idv.depth = 1;
    }
    double dice;
    for (long long generation = 0; generation < max_gen; generation++)
    {
        for (long long good_one = 0; good_one < population / 5; good_one++)
        {
            for (int k = 1; k < 5; k++)
            {
                vector<vector<double>> mutGRN(genes, vector<double>(genes));
                for (int i = 0; i < genes; i++)
                {
                    for (int j = 0; j < genes; j++)
                    {
                        dice = rand(mt);
                        if (dice < mutation_rate)
                        {
                            mutGRN[i][j] = rand_for_theta(mt);
                        }
                        else
                        {
                            mutGRN[i][j] = idv_vec[good_one].GRN[i][j];
                        }
                    }
                }

                vector<double> muttheta(genes);
                for (int i = 0; i < genes; i++)
                {
                    dice = rand(mt);
                    if (dice < mutation_rate)
                    {
                        muttheta[i] = rand_for_theta(mt);
                    }
                    else
                    {
                        muttheta[i] = idv_vec[good_one].theta[i];
                    }
                }
                idv_vec[good_one * 5 + k] = gr_system{genes, mutGRN, muttheta};
                calculate_fitness_step(genes, noise_arr, initial_state, idv_vec[good_one * 5 + k], num_noise_type);
            }
        }
        sort(idv_vec.begin(), idv_vec.end(), compare_by_fitness);

        depth_writing_file << idv_vec[0].depth << endl;
        sysinfo_writing_file << idv_vec[0].num_final_attractors << endl;

        cout << generation << " " << idv_vec[0].depth << " " << idv_vec[0].num_final_attractors << endl;
        for (int i = 0; i < genes; i++)
        {
            theta_writing_file << idv_vec[0].theta[i] << " ";
            for (int j = 0; j < genes; j++)
            {
                GRN_writing_file << idv_vec[0].GRN[i][j] << " ";
            }
            GRN_writing_file << endl;
        }
        theta_writing_file << endl;

        for (int atr = 0; atr < idv_vec[0].attractors.size(); atr++)
        {
            for (auto ch : idv_vec[0].attractors[atr])
            {
                for (int i = 0; i < genes; i++)
                {
                    attractors_writing_file << ch[i] << " ";
                }
            }
            attractors_writing_file << "; ";
        }
        attractors_writing_file << endl;

        for (int atr = 0; atr < idv_vec[0].goes_to.size(); atr++)
        {
            for (int dest = 0; dest < idv_vec[0].goes_to[atr].size(); dest++)
            {
                goes_to_writing_file << idv_vec[0].goes_to[atr][dest] << " ";
                if (idv_vec[0].goes_to[atr][dest] != 0)
                {
                    differentiate_writing_file << dest << " ";
                }
            }
            differentiate_writing_file << "; ";
            goes_to_writing_file << "; ";
        }
        differentiate_writing_file << endl;
        goes_to_writing_file << endl;
    }
}
