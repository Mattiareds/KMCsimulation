#include <vector>
#include <fstream>
#include <random>
#include <array>
#include <tuple>


#ifndef METROPOLIS_H
#define METROPOLIS_H


class metropolis{
    private:
    coordinates& crd;
    coordinates& c_crd;
    geometry& s; //sites
    std::random_device rd;
    std::mt19937 rng;
    std::uniform_real_distribution<double> initial_pos;
    std::uniform_real_distribution<double> dist_time;

    std::vector<bool> atoms;
    std::vector<int> interested_sites; //sites interested in the last iteration
    std::vector<int> nnn_atoms; //number of nn for each atom
    int movements_per_class[100]; // sarà ridimensionato dopo aver letto n_class
    std::vector<std::vector<int>> table_of_processes; //map initial site -> type of possible processes
    std::vector<std::vector<int>> table_of_end_pos; //map initial site -> possible end sites 
    std::vector<std::vector<int>> table_of_initial_pos; //map finale site -> initial site
    std::vector<int> deactivated_sites;
    std::vector<double> P; // sarà ridimensionato dopo aver letto n_class
    std::vector<std::tuple<int,int,int>> dynamic_processes;
    //std::vector<std::vector<std::vector<int>>>  map_class_processes;//map class -> process
    int map_class_processes[100][1000][2];

    std::vector<float> barriere;
    static constexpr double Kb = 8.617333262e-5;
    int n_deposited;
    float E_b;
    float temperatura;
    float nu_0;
    double F; //deposition rate
    float filling;
    int n_class;
    float total_time;
    int pos; //current pos
    int next;
    double time=0.0;
    bool deposition=false;
    std::string process_name;
    std::ofstream output;

    //int classes [11]={0, 6, 10, 15, 19, 24, 28, 31, 36, 41, 44}; //classes
    
    void interested_sites_calc(bool diff); 
    void nn_updater(bool b);
    double probability_calculator(int a);
    void table_of_processes_filler(); 
    void algorithm();
    void second_layer_updates(int i);
    void second_layer_activation();
    void classification();
    void probability_filler();
    void start_of_the_sim();
    void time_prob_calc();
    void map_of_class_eraser(int c, int i, int j);
    void map_of_class_position_eraser(int i, int nn); 
    void map_of_class_next_eraser(int i); 
    void map_of_class_filler(int c, int i , int j);
    int get_MCP_size(int c);
    void deposition_func();
    void print_output();
    void print_configuration();
    void print_final_output();
    void barrier_chooser(int cl, int &i, int& j);
    void get_process_name(int cl);
    void shifter(int c, int ind, int tg);
    bool is_deactivated(int i);
    bool is_border(int i);
    bool are_different_border(int i, int j);
    std::array<int,2> activate_edges(int upper_site, int pv);
    void table_of_processes_updater(int i);
    void clear_dynamic_processes();

    public:
    //builder
    metropolis(coordinates& structure, coordinates& core, geometry& sites)
        : crd(structure), c_crd(core), s(sites),
          rng(rd()), initial_pos(0.0,1.0) , dist_time(0.0,1.0)
    { 
    }

    //std::ofstream debug_file;

    void file_reader(std::ifstream& ifile);
    void simulation(); //questa unica accessibile da esterno

};

#endif