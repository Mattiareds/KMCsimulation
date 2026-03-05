#include <vector>
#include <fstream>
#include <random>

#ifndef METROPOLIS_H
#define METROPOLIS_H


class metropolis{
    private:
    coordinates& crd;
    geometry& s; //sites
    geometry& h_s; //high sites
    std::random_device rd;
    std::mt19937 rng;
    std::uniform_real_distribution<double> initial_pos;
    std::uniform_real_distribution<double> dist_time;

    std::vector<bool> atoms;
    std::vector<int> interested_sites; //sites interested in the last iteration
    std::vector<int> nnn_atoms; //number of nn for each atom
    std::vector<int> atoms_per_class; // sarà ridimensionato dopo aver letto n_class
    std::vector<std::vector<int>> table_of_processes; //map initial site -> type of possible processes
    std::vector<std::vector<int>> table_of_end_pos; //map initial site -> possible end sites 
    std::vector<double> P; // sarà ridimensionato dopo aver letto n_class
    std::vector<std::vector<std::vector<int>>>  map_class_processes;//map class -> process

    std::vector<double> barriere;
    double E_b; 
    double Kb=0.00008617;
    int mc_step;
    int metodo;
    double temperatura;
    double nu_0;
    double F; //deposition rate
    double filling;
    int n_class;
    double total_time;
    int pos; //current pos
    int next;
    double time=0.0;
    bool deposition=false;
    std::string process_name;
    std::ofstream output;

    int* classes = new int[11]{0, 6, 10, 15, 19, 24, 28, 31, 36, 41, 44}; //classes
    
    void interested_sites_calc(bool diff); 
    void nn_updater(bool b);
    double probability_calculator(int a);
    void table_of_processes_filler(); 
    void algorithm();
    void second_layer_updates();
    void classification();
    void probability_filler();
    void start_of_the_sim();
    void time_prob_calc();
    void map_of_class_eraser(int c, int i, int j);
    void deposition_func();
    void print_output();
    int * barrier_chooser(int cl);

    public:
    //builder
    metropolis(coordinates& structure, geometry& sites, geometry& high_sites)
        : crd(structure), s(sites), h_s(high_sites),
          rng(rd()), initial_pos(0.0,1.0) , dist_time(0.0,1.0)
    { 
    }

    void file_reader(std::ifstream& ifile);
    void simulation(); //questa unica accessibile da esterno

};

#endif