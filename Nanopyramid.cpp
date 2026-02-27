#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"

using namespace std;

int main(){
    coordinates core;
    coordinates shell;
    coordinates high_shell;
    
    geometry g_core;
    geometry g_shell;
    geometry g_high_shell;

    metropolis alg(core, g_core, g_shell);

    g_core.set(core,alg);
    g_shell.set(shell, alg);
    
    ifstream core_coord("piramide.xyz");
    ifstream shell_coord("siti_jmol.xyz");
    ifstream high_shell_coord("high_sites.xyz");
    ifstream settings("mc_settings.in");
    ofstream output_file("out.out");

    core.initialize(core_coord);
    shell.initialize(shell_coord);
    //DA FINIRE DI SETTARE E FARE LE PRIME PROVE 
}

//compiler: g++ geometry.cpp coordinates.cpp Nanopyramid.cpp -o Nanopyramid