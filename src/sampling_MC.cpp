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
    int count = 0, itr = 0;
    gr_system::state destination(genes);
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
                is_fp = false;
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
    int count = 0, itr = 0;
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

void calculate_fitness_step(const int genes, vector<gr_system::state> &noise_arr, const gr_system::state &initial_state, gr_system &System)
{
    int num_noise_type = noise_arr.size();
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

bool is_flat(const vector<long long> &Histogram)
{
    double ave = accumulate(Histogram.begin(), Histogram.end(), 0LL) * 1.0 / Histogram.size();
    for (int h = 0, H = Histogram.size(); h < H; h++)
    {
        if (Histogram[h] > ave * 0.9)
        {
            continue;
        }
        else
        {
            return false;
        }
    }
    return true;
}

int main(int argc, const char **argv)
{
    int genes = stoi(argv[1]);
    int num_noise_type = stoi(argv[2]);
    int max_fitness = stoi(argv[3]);
    double delta = stod(argv[4]);
    string str_delta = argv[4];

    int seed_for_perturbation = 1;
    int seed = stoi(argv[5]);
    mt19937_64 mt(seed);

    uniform_real_distribution<double> rand_for_theta(-1, 1);
    gr_system::state initial_state(genes);
    for (int i = 0; i < genes; i++)
    {
        initial_state[i] = 0;
    }

    ifstream N_noise_read("../perturbation_lib/N_noise_genes" + to_string(genes) + "_delta_" + str_delta + "_seed_" + to_string(seed_for_perturbation) + "_noise_set_N=" + to_string(num_noise_type) + "_including_e.txt");
    vector<gr_system::state> noise_arr(num_noise_type, gr_system::state(genes, 0));
    for (int type = 0; type < num_noise_type; type++)
    {
        for (int i = 0; i < genes; i++)
        {
            N_noise_read >> noise_arr[type][i];
        }
    }

    string input_dir = "../data_random_noise" + str_delta + "/genes_" + to_string(genes) + "_MC_step/num_noise_type_" + to_string(num_noise_type) + "/max=" + to_string(max_fitness);

    vector<double> S(max_fitness + 1, 0);
    vector<long long> H(max_fitness + 1, 0);

    double value;
    int renewed = 0;
    ifstream H_reading_file(input_dir + "/hist_f=2^-0.txt");
    while (H_reading_file.is_open())
    {
        for (int fit = 0; fit < max_fitness + 1; fit++)
        {
            H_reading_file >> value;
            S[fit] += value * pow(2, -renewed);
        }
        renewed++;
        H_reading_file.close();
        H_reading_file.open(input_dir + "/hist_f=2^-" + to_string(renewed) + ".txt");
    }
    if (renewed < 20)
    {
        return 0;
    }

    string output_dir = input_dir + "/samples/seed_" + to_string(seed);

    vector<ofstream> matrix_output_vec(max_fitness + 1);
    vector<ofstream> theta_output_vec(max_fitness + 1);
    vector<ofstream> differentiate_network_output_vec(max_fitness + 1);
    vector<ofstream> goes_to_output_vec(max_fitness + 1);
    vector<ofstream> attractors_output_vec(max_fitness + 1);
    vector<ofstream> depth_output_vec(max_fitness + 1);
    for (int fit = 0; fit < max_fitness + 1; fit++)
    {
        matrix_output_vec[fit].open(output_dir + "/matrix_fit=" + to_string(fit) + ".txt", ios::out);
        theta_output_vec[fit].open(output_dir + "/theta_fit=" + to_string(fit) + ".txt", ios::out);
        differentiate_network_output_vec[fit].open(output_dir + "/differentiate_network_fit=" + to_string(fit) + ".txt", ios::out);
        attractors_output_vec[fit].open(output_dir + "/attractors_fit=" + to_string(fit) + ".txt", ios::out);
        goes_to_output_vec[fit].open(output_dir + "/goes_to_fit=" + to_string(fit) + ".txt", ios::out);
        depth_output_vec[fit].open(output_dir + "/depth_fit=" + to_string(fit) + ".txt", ios::out);
    }

    uniform_real_distribution<double> rand(0, 1);
    double mutation_rate = 0.05;

    vector<vector<double>> GRN_init(genes, vector<double>(genes));
    for (int i = 0; i < genes; i++)
    {
        for (int j = 0; j < genes; j++)
        {
            GRN_init[i][j] = rand_for_theta(mt);
        }
    }

    vector<double> theta_init(genes);
    for (int i = 0; i < genes; i++)
    {
        theta_init[i] = rand_for_theta(mt);
    }

    gr_system System_former(genes, GRN_init, theta_init);
    calculate_fitness_step(genes, noise_arr, initial_state, System_former);

    long long step = 0;
    for (long long sample = 0; sample < 1e6; sample++)
    {
        vector<vector<double>> GRN(genes, vector<double>(genes));
        for (int i = 0; i < genes; i++)
        {
            for (int j = 0; j < genes; j++)
            {
                if (rand(mt) < mutation_rate)
                {
                    GRN[i][j] = rand_for_theta(mt);
                }
                else
                {
                    GRN[i][j] = System_former.GRN[i][j];
                }
            }
        }
        vector<double> theta(genes);
        for (int i = 0; i < genes; i++)
        {
            if (rand(mt) < mutation_rate)
            {
                theta[i] = rand_for_theta(mt);
            }
            else
            {
                theta[i] = System_former.theta[i];
            }
        }
        gr_system System_current(genes, GRN, theta);
        calculate_fitness_step(genes, noise_arr, initial_state, System_current);

        double prob = min(1., exp(S[System_former.num_final_attractors] - S[System_current.num_final_attractors]));
        if (rand(mt) < prob && System_current.num_final_attractors <= max_fitness)
        {
            System_former = System_current;
        }

        if (step % 100 == 0)
        {
            int fit_former = System_former.num_final_attractors;
            for (int i = 0; i < genes; i++)
            {
                theta_output_vec[fit_former] << System_former.theta[i] << " ";
                for (int j = 0; j < genes; j++)
                {
                    matrix_output_vec[fit_former] << System_former.GRN[i][j] << " ";
                }
                matrix_output_vec[fit_former] << endl;
            }
            theta_output_vec[fit_former] << endl;

            vector<set<int>> st(System_former.goes_to.size(), set<int>{});
            int SIZE = System_former.goes_to.size();
            attractors_output_vec[fit_former] << System_former.attractors.size() << " ";
            for (int atr = 0; atr < SIZE; atr++)
            {
                for (int ch = 0; ch < SIZE; ch++)
                {
                    goes_to_output_vec[fit_former] << System_former.goes_to[atr][ch] << " ";
                    if (System_former.goes_to[atr][ch] != 0)
                    {
                        differentiate_network_output_vec[fit_former] << ch << " ";
                    }
                }
                goes_to_output_vec[fit_former] << "; ";
                differentiate_network_output_vec[fit_former] << "; ";
                for (auto ch_atr : System_former.attractors[atr])
                {
                    for (auto v : ch_atr)
                    {
                        attractors_output_vec[fit_former] << v << " ";
                    }
                }
                attractors_output_vec[fit_former] << "; ";
            }
            goes_to_output_vec[fit_former] << endl;
            attractors_output_vec[fit_former] << endl;
            differentiate_network_output_vec[fit_former] << endl;

            depth_output_vec[fit_former] << System_former.depth << endl;
        }
        step++;
    }
}