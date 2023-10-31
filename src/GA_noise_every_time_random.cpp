#pragma GCC target("avx")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

#include <vector>
#include <set>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
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
    int num_all_attractors = 0;
    int num_progenitors = 0;
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

int move_not_on_plane(const gr_system &System, gr_system::state &current_state, vector<vector<double>> &constraint_A, vector<double> &constraint_b, vector<int> &constraint_idnum, vector<bool> &on_plane, gr_system::state &destination_effective, double &c, int &c_id)
{
    double infinitesimal = 1e-10;
    double c_;
    bool has_op = false;
    vector<double> dxdt = destination_effective - current_state;
    for (int i = 0; i < System.genes; i++)
    {
        c_ = (System.theta[i] - dot(System.GRN[i], current_state)) / dot(System.GRN[i], dxdt);
        if (c_ > 0 && c_ < c)
        {
            c = c_;
            c_id = i;
        }
    }
    if (c > 1 - 1e-2)
    {
        current_state = destination_effective;
        return 1;
    }
    else
    {
        current_state += (c + infinitesimal) * dxdt;
        for (int i = 0; i < System.genes; i++)
        {
            if (current_state[i] > 1)
            {
                current_state[i] = 1 - infinitesimal;
            }
            else if (current_state[i] < 0)
            {
                current_state[i] = infinitesimal;
            }
        }
        for (int i = 0; i < System.genes; i++)
        {
            destination_effective[i] = (double)(dot(current_state, System.GRN[i]) - System.theta[i] > 0);
            if (destination_effective[i] >= 1)
            {
                destination_effective[i] = 1 - infinitesimal;
            }
            else if (destination_effective[i] <= 0)
            {
                destination_effective[i] = infinitesimal;
            }
        }
        return 0;
    }
}

int move_on_plane(const gr_system &System, gr_system::state &current_state, vector<vector<double>> &constraint_A, vector<double> &constraint_b, vector<int> &constraint_idnum, vector<bool> &on_plane, gr_system::state &destination_effective, double &c, int &c_id)
{
    double infinitesimal = 1e-10;
    double c_;
    bool has_op = false;
    vector<double> dxdt = destination_effective - current_state;
    for (int i = 0; i < System.genes; i++)
    {
        if (on_plane[i])
        {
            continue;
        }
        c_ = (System.theta[i] - dot(System.GRN[i], current_state)) / dot(System.GRN[i], dxdt);
        if (c_ > 0 && c_ < c)
        {
            c = c_;
            c_id = i;
        }
    }

    double a;
    for (int i = 0; i < System.genes; i++)
    {
        if (on_plane[i + System.genes])
        {
            continue;
        }
        if (abs(dxdt[i]) < 1e-4)
        {
            continue;
        }
        c_ = (-current_state[i]) / dxdt[i];
        if (c_ > 0 && c_ < c)
        {
            c = c_;
            c_id = i + System.genes;
            a = infinitesimal; // i.e. 0
        }
        c_ = (1 - current_state[i]) / dxdt[i];
        if (c_ > 0 && c_ < c)
        {
            c = c_;
            c_id = i + System.genes;
            a = 1 - infinitesimal; // i.e. 1
        }
    }

    if (c > 1 - 1e-2)
    {
        current_state = destination_effective;
        return 1;
    }
    else
    {
        current_state += (c + infinitesimal) * dxdt;
        if (c_id > System.genes)
        {
            on_plane[c_id] = true;
        }
        else
        {
            on_plane[c_id + System.genes] = false;
        }

        for (int i = 0; i < System.genes; i++)
        {
            if (current_state[i] > 1)
            {
                current_state[i] = 1 - infinitesimal;
            }
            else if (current_state[i] < 0)
            {
                current_state[i] = infinitesimal;
            }
        }
        for (int i = 0; i < System.genes; i++)
        {
            on_plane[i] = false;
            destination_effective[i] = (double)(dot(current_state, System.GRN[i]) - System.theta[i] > 0);
            if (destination_effective[i] >= 1)
            {
                destination_effective[i] = 1 - infinitesimal;
            }
            else if (destination_effective[i] <= 0)
            {
                destination_effective[i] = infinitesimal;
            }
        }
        constraint_A = vector<vector<double>>(0, vector<double>(System.genes));
        constraint_b = vector<double>(0);
        constraint_idnum = vector<int>(System.genes * 2, -1);
        return 0;
    }
}

void check(const int genes, gr_system::state &current_state, bool &is_fp, const gr_system &System)
{
    int count = 0;
    double infinitesimal = 1e-10;
    gr_system::state destination(genes), destination_effective(genes);
    for (int i = 0; i < genes; i++)
    {
        destination_effective[i] = (double)(dot(current_state, System.GRN[i]) - System.theta[i] > 0);
        if (destination_effective[i] >= 1)
        {
            destination_effective[i] = 1 - infinitesimal;
        }
        else if (destination_effective[i] <= 0)
        {
            destination_effective[i] = infinitesimal;
        }
    }
    gr_system::state former_state = current_state, dx = former_state - current_state;
    vector<double> dxdt = destination - current_state;
    bool done = false;
    vector<vector<double>> constraint_A(0, vector<double>(System.genes));
    vector<double> constraint_b(0);
    vector<int> constraint_idnum(System.genes * 2, -1);
    vector<bool> on_plane(genes * 2, false), new_plane(genes, false);
    for (int iteration = 0; iteration < 5; iteration++)
    {
        int c_id = -1;
        double c = 1e12;
        if (move_not_on_plane(System, current_state, constraint_A, constraint_b, constraint_idnum, on_plane, destination_effective, c, c_id) == 1)
        {
            done = true;
            break;
        };
    }

    while (!done)
    {
        former_state = current_state;
        int c_id = -1;
        double c = 1e12;
        if (any_of(on_plane.begin(), on_plane.begin() + System.genes, [](double on_plane_i)
                   { return on_plane_i; }))
        {
        }
        else
        {
            if (move_not_on_plane(System, current_state, constraint_A, constraint_b, constraint_idnum, on_plane, destination_effective, c, c_id) == 1)
            {
                break;
            }
        }
        vector<double> dxdt = destination_effective - current_state;
        dx = former_state - current_state;
        if (all_of(dx.begin(), dx.end(), [](double dx_i)
                   { return abs(dx_i) < 1e-4; }))
        {
            for (int c_add_id = 0; c_add_id < genes; c_add_id++)
            {
                if (System.GRN[c_add_id][c_add_id] < 0 && abs(System.theta[c_add_id] - dot(System.GRN[c_add_id], current_state)) < 1e-8 && !on_plane[c_add_id + System.genes])
                {
                    on_plane[c_add_id] = true;
                    new_plane[c_add_id] = true;
                }
            }
            // calculate_hyperplane_intersection;
            MatrixXd A(genes, genes);
            VectorXd b(genes);
            for (int i = 0; i < genes; i++)
            {
                if (on_plane[i])
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
                    b(i) = double(dot(System.GRN[i], current_state) - System.theta[i] > 0);
                    if (b(i) > 1)
                    {
                        b(i) -= 2 * (b(i) - 1);
                    }
                    if (b(i) < 0)
                    {
                        b(i) = -b(i);
                    }
                }
            }

            VectorXd x = A.colPivHouseholderQr().solve(b);
            // move()
            for (int i = 0; i < genes; i++)
            {
                if (on_plane[i])
                {
                    destination_effective[i] = x(i);
                    if (new_plane[i])
                    {
                        constraint_A.push_back(System.GRN[i]);
                        constraint_b.push_back(System.theta[i]);
                        constraint_idnum[i] = constraint_A.size() - 1;
                    }
                }
                else
                {
                    destination_effective[i] = double(dot(System.GRN[i], current_state) - System.theta[i] > 0);
                    if (destination_effective[i] >= 1)
                    {
                        destination_effective[i] = 1 - infinitesimal;
                    }
                    else if (destination_effective[i] <= 0)
                    {
                        destination_effective[i] = infinitesimal;
                    }
                }
            }

            int c_id = -1;
            double c = 1e12;
            if (move_on_plane(System, current_state, constraint_A, constraint_b, constraint_idnum, on_plane, destination_effective, c, c_id) == 1)
            {
                break;
            }
            new_plane = vector<bool>(genes, false);
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

void dfs_step(const int genes, const int num_noise_type, const double delta, normal_distribution<double> &normal_dist, mt19937_64 &mt, long long &random_count, int mother, gr_system &System, bool &has_lc)
{
    vector<gr_system::state> noise_arr(num_noise_type, gr_system::state(genes, 0));
    double norm;
    for (int n = 0; n < num_noise_type; n++)
    {
        norm = 0;
        for (int i = 0; i < genes; i++)
        {
            noise_arr[n][i] = normal_dist(mt);
            random_count++;
            norm += pow(noise_arr[n][i], 2);
        }
        for (int i = 0; i < genes; i++)
        {
            noise_arr[n][i] = noise_arr[n][i] / sqrt(norm) * delta;
        }
    }
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
            dfs_step(genes, num_noise_type, delta, normal_dist, mt, random_count, id, System, has_lc);
        }
    }
}

void calculate_fitness_step_every_time_random(const int genes, const int num_noise_type, const double delta, const gr_system::state &initial_state, gr_system &System, normal_distribution<double> &normal_dist, mt19937_64 &mt, long long &random_count)
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
        System.num_all_attractors = 0;
        System.num_progenitors = System.num_all_attractors - System.num_final_attractors;
        return;
    }

    bool has_lc = false;
    dfs_step(genes, num_noise_type, delta, normal_dist, mt, random_count, 0, System, has_lc);
    if (has_lc)
    {
        System.num_final_attractors = 0;
        System.num_all_attractors = 0;
        System.num_progenitors = System.num_all_attractors - System.num_final_attractors;
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
        System.num_all_attractors = System.attractors.size();
        System.num_progenitors = System.num_all_attractors - System.num_final_attractors;

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

int main(int argc, const char **argv)
{
    int genes = stoi(argv[1]);
    int num_noise_type = stoi(argv[2]);
    int seed = stoi(argv[3]);
    double delta = stod(argv[4]);
    string str_delta = argv[4];
    mt19937_64 mt(seed);

    string dir_name = "../data_random_noise" + str_delta + "/genes_" + to_string(genes) + "_step/num_noise_type_" + to_string(num_noise_type) + "_every_time_random/seed_" + to_string(seed);
    ofstream fitness_writing_file(dir_name + "/fitness.txt");
    ofstream num_attractors_writing_file(dir_name + "/num_attractors.txt");
    ofstream random_counts_writing_file(dir_name + "/random_counts.txt");
    ofstream GRN_writing_file(dir_name + "/GRN.txt");
    ofstream theta_writing_file(dir_name + "/theta.txt");
    ofstream depth_writing_file(dir_name + "/depth.txt");
    ofstream differentiate_writing_file(dir_name + "/diffenrentiate_network.txt");
    ofstream attractors_writing_file(dir_name + "/attractors.txt");
    ofstream goes_to_writing_file(dir_name + "/atr_adj_mat.txt");
    ofstream in_process_GRN_population_file;
    ofstream in_process_theta_population_file;

    uniform_real_distribution<double> rand_for_theta(-1, 1);
    normal_distribution<double> normal_dist(0., 1.);
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
    long long random_count = 0;
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
                        random_count++;
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
                    random_count++;
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
                calculate_fitness_step_every_time_random(genes, num_noise_type, delta, initial_state, idv_vec[good_one * 5 + k], normal_dist, mt, random_count);
            }
        }
        sort(idv_vec.begin(), idv_vec.end(), compare_by_fitness);

        depth_writing_file << idv_vec[0].depth << endl;
        fitness_writing_file << idv_vec[0].num_final_attractors << endl;
        num_attractors_writing_file << idv_vec[0].attractors.size() << endl;
        random_counts_writing_file << random_count << endl;
        if (generation % 100 == 0)
        {
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
        if (generation % 10000 == 0)
        {
            in_process_GRN_population_file.open(dir_name + "/in_process_GRN.txt");
            in_process_theta_population_file.open(dir_name + "/in_process_theta.txt");
            in_process_GRN_population_file << generation << " " << random_count << endl;
            in_process_theta_population_file << generation << " " << random_count << endl;
            for (long long idv = 0; idv < population; idv++)
            {
                for (int i = 0; i < genes; i++)
                {
                    for (int j = 0; j < genes; j++)
                    {
                        in_process_GRN_population_file << setprecision(15) << idv_vec[idv].GRN[i][j] << " ";
                    }
                    in_process_GRN_population_file << endl;
                    in_process_theta_population_file << setprecision(15) << idv_vec[idv].theta[i] << " ";
                }
                in_process_theta_population_file << endl;
            }
            in_process_GRN_population_file.close();
            in_process_theta_population_file.close();
        }
    }
}
