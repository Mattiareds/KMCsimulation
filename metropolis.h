#include <vector>
#include <fstream>
#include <random>
#include <array>
#include <tuple>
#include <chrono>
#include <string>


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
    static constexpr int MAX_CLASSES = 150; //size of the class arrays, must be > highest base class + 12
    static constexpr int MAX_PROCESSES = 1000; //max processes stored per class
    int movements_per_class[MAX_CLASSES];
    std::vector<std::vector<int>> table_of_processes; //map initial site -> type of possible processes
    std::vector<std::vector<int>> table_of_end_pos; //map initial site -> possible end sites 
    std::vector<std::vector<int>> table_of_initial_pos; //map finale site -> initial site
    std::vector<int> deactivated_sites;
    std::vector<double> P; 
    std::vector<std::tuple<int,int,int>> dynamic_processes;
    //std::vector<std::vector<std::vector<int>>>  map_class_processes;//map class -> process
    int map_class_processes[MAX_CLASSES][MAX_PROCESSES][2];

    std::vector<float> barriere;
    int additional_cost_corner[2];
    static constexpr double Kb = 8.617333262e-5;
    int n_deposited = 0;
    float E_b;
    float temperatura;
    float nu_0;
    double F; //deposition rate
    float filling;
    int n_class = 0;
    float total_time;
    int pos = 0; //current pos
    int next = 0;
    double time=0.0;
    bool deposition=false;
    bool output_file;
    bool fill_first_layer;
    bool print_bond_file;
    bool differentiate_corner_sites = false;
    double max_runtime_hours = 0.0; //wall clock limit of the program, 0 = no limit
    std::chrono::steady_clock::time_point start_runtime;
    std::string stop_reason = "physical simulation time reached";
    bool final_report_written = false;
    int dropped_from = -1;
    int dropped_target = -1;
    std::string process_name;
    std::ofstream output;
    std::ofstream occupation_file;
    std::ofstream bond_file;
    std::ofstream run_info_file; //always written: settings of the run and why it stopped

    
    void interested_sites_calc(bool diff); 
    void nn_updater(bool b);
    double probability_calculator(int a);
    void table_of_processes_filler(); 
    void algorithm();
    void second_layer_updates(int i);
    void second_layer_activation();
    void activate_all();
    void second_layer_deactivation(int deserted_site);
    void classification();
    void probability_filler();
    void first_layer_filler();
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
    void print_final_configuration();
    void print_final_output();
    void print_occupation();
    void write_run_info();
    void write_final_report(const std::string& reason);
    [[noreturn]] void stop_simulation(const std::string& reason);
    double elapsed_runtime() const; //seconds since the start of the simulation
    bool runtime_exceeded() const;
    void barrier_chooser(int cl, int &i, int& j);
    void get_process_name(int cl);
    void shifter(int c, int ind, int tg);
    bool is_deactivated(int i);
    bool is_border(int i);
    bool is_corner(int i);
    bool are_different_border(int i, int j);
    bool are_different_corner(int i, int j);
    std::array<int,2> activate_edges(int upper_site, int pv);
    void table_of_processes_updater(int i);
    void clear_dynamic_processes();
    int bond_breaking_counter = 0;
    int bond_breaking_size_counter[8] = {0}; 
    int barrier_crossing_counter = 0;

    public:
    //builder
    metropolis(coordinates& structure, coordinates& core, geometry& sites)
        : crd(structure), c_crd(core), s(sites),
          rng(rd()), initial_pos(0.0,1.0) , dist_time(0.0,1.0)
    { 
    }

    //std::ofstream debug_file;

    void file_reader(std::ifstream& ifile);
    void simulation();

    //getter
    bool get_differentiate_corner_sites() const { return differentiate_corner_sites; }

};

#endif