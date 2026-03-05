#include <vector>
#include <fstream>

#ifndef COORDINATES_H
#define COORDINATES_H

class coordinates{
    private:
    std::string second_line;
    std::vector<std::string> specie_chimica; //chemical species
    std::vector<std::vector<double>> siti; //coodrdinates of diffusion sites (vector Nx3)

    
    int N; //number of coordinates
    int M; //number of sites of the upper level
    double d_pv;
    int spigolo;
    int* n_cut = new int[4];
    
    public:
    
    //read some coordinates
    void initialize(std::ifstream& ifile);
    //read settings file
    void settings_reader(std::ifstream& ifile);
    //set the output file 
    void output_writer(std::ofstream& ofile, int n_pos);

    //setters and getters
    const std::string get_second_line() const {return second_line;}

    const std::string get_chem_sp(int i) const {return specie_chimica.at(i);}
    //get the components of the position of a particle
    const double site_x(int i) const {return siti[i][0];}
    const double site_y(int i) const {return siti[i][1];}
    const double site_z(int i) const {return siti[i][2];}
    //get the number of coordinates
    const int sites_size() const {return siti.size();}
    //get the number of particles
    const int get_N() const {return N;}
    //get the length of the edges wihout cuts
    const int get_spigolo() const {return spigolo;}
    //get the number of cutted planes
    const int get_Ncut(int indx) const {
        if(indx < 4) return n_cut[indx];
        else return 0;
    }
    const double get_dpv() const {return d_pv;}


};

#endif