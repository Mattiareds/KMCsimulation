#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"
#include <iostream>

using namespace std;

int main(){
    //define objects of the classes
    coordinates core;
    coordinates shell;
    coordinates high_shell;
    
    geometry g_core;
    geometry g_shell;
    geometry g_high_shell;

    metropolis alg(core, g_core, g_shell);

    //set the core and the shell
    g_core.set(core,alg);
    g_shell.set(shell, alg);
    
    //input files
    ifstream core_coord("piramide.xyz");
    ifstream shell_coord("siti_jmol.xyz");
    ifstream high_shell_coord("high_sites.xyz");
    ifstream settings("mc_settings.in");
    ifstream coordinates_settings("coordinates_settings.in");
    ifstream shell_settings("shell_settings.in");
    ofstream output_file("out.out");

    //read input files
    core.initialize(core_coord);
    shell.initialize(shell_coord);
    core.settings_reader(coordinates_settings);
    shell.settings_reader(shell_settings);

    alg.file_reader(settings);

    //set and calculate list of nn, geometry, planes etc
    g_core.starter();
    g_shell.starter();
    
    std::cout<<"e fin qui..."<<std::endl;
    alg.simulation();

}

//compiler: g++ geometry.cpp coordinates.cpp Nanopyramid.cpp -o Nanopyramid