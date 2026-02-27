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

    void initialize(std::ifstream& ifile);
    void settings_reader(std::ifstream& ifile);
    void output_writer(std::ofstream& ofile, int n_pos);

    //ragionare su quali funzioni siano davvero accessibili dall'esterno e quali no, probabilmente basta che alcuni membri siano accessibili così da poterli utilizzare in altre classi 
    //invece tutte le altre funzioni vanno messe private
    const std::string get_second_line() const {return second_line;}
    const std::string get_chem_sp(int i) const {return specie_chimica.at(i);}
    const int site_x(int i) const {return siti[i][0];}
    const int site_y(int i) const {return siti[i][1];}
    const int site_z(int i) const {return siti[i][2];}
    const int sites_size() const {return siti.size();}
    const int get_N() const {return N;}
    const int get_spigolo() const {return spigolo;}
    const int get_Ncut(int indx) const {
        if(indx < 4) return n_cut[indx];
    }
    const int get_dpv() const {return d_pv;}


};

#endif