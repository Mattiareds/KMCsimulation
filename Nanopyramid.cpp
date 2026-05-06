#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
    //define objects of the classes
    coordinates core;
    coordinates shell;
    
    geometry g_core;
    geometry g_shell;

    metropolis alg(shell,core, g_shell); 

    //set the core and the shell
    g_core.set(core,alg);
    g_shell.set(shell, alg);
    
    //input files
    ifstream core_coord("piramide.xyz");
    if (!core_coord.is_open()) {
        cerr << "ERROR: Could not open file 'piramide.xyz'" << endl;
        return 1;
    }
    
    ifstream shell_coord("siti_jmol.xyz");
    if (!shell_coord.is_open()) {
        cerr << "ERROR: Could not open file 'siti_jmol.xyz'" << endl;
        return 1;
    }
    
    ifstream settings("mc_settings.in");
    if (!settings.is_open()) {
        cerr << "ERROR: Could not open file 'mc_settings.in'" << endl;
        return 1;
    }
    
    ifstream coordinates_settings("coordinates_settings.in");
    if (!coordinates_settings.is_open()) {
        cerr << "ERROR: Could not open file 'coordinates_settings.in'" << endl;
        return 1;
    }
    
    ifstream shell_settings("shell_settings.in");
    if (!shell_settings.is_open()) {
        cerr << "ERROR: Could not open file 'shell_settings.in'" << endl;
        return 1;
    }
    
    ofstream output_file("out.out");
    if (!output_file.is_open()) {
        cerr << "ERROR: Could not create output file 'out.out'" << endl;
        return 1;
    }

    cout << "All input files successfully opened." << endl;

    //read input files
    core.initialize(core_coord);
    shell.initialize(shell_coord);
    core.settings_reader(coordinates_settings);
    shell.settings_reader(shell_settings);

    alg.file_reader(settings);

    //set and calculate list of nn, geometry, planes etc
    g_core.starter();
    g_shell.starter();
    
    g_shell.pv_substitution();
    g_shell.test_border(g_core);
   
    //simulation
    alg.simulation();
    
    return 0;
}