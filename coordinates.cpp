#include "coordinates.h"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <iostream>

//reader of the coordinates
void coordinates::initialize(std::ifstream& ifile1){
    if(ifile1){
        ifile1>>N;
        ifile1.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //move pointer to second line
        std::getline(ifile1, second_line);
        siti.resize(N, std::vector<double>(3,0.0));
        specie_chimica.resize(N);
        for(int i=0;i<N;i++){
            ifile1>>specie_chimica[i]>>siti[i][0]>>siti[i][1]>>siti[i][2];
        }
    }else std::cout<< "From site.h: no main input file!!! " << std::endl;
}

//reader of the settings
void coordinates::settings_reader(std::ifstream& ifile){
    if(ifile){
        ifile>>d_pv>>spigolo;
        for(int i=0;i<4;i++){
            ifile>>n_cut[i];
        }
    }else std::cout<< "From site.h: no settings file!!!" <<std::endl;
}

//output of coordinates
void coordinates::output_writer(std::ofstream& ofile, int n_pos){
    ofile<<n_pos<<std::endl;
    ofile << second_line << std::endl;
    for(int i=0;i<siti.size();i++){
        ofile<<specie_chimica[i]<<" "<<siti[i][0]<<" "<<siti[i][1]<<" "<<siti[i][2]<<std::endl;
    }
}
