#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <new>
#include <unordered_set>

/*
 * Reads the simulation parameters from the input file.
 * Initializes energy barriers, temperature, flux (F), and sets up the internal 
 * data structures for process classification .
 */
void metropolis::file_reader(std::ifstream& ifile){
    if(ifile){
        int n_barriers;
        ifile>>n_barriers;
        barriere.resize(n_barriers);
        for(int i=0 ; i<n_barriers ; i++){
            ifile>>barriere[i];
        }
        ifile>>E_b;
        ifile>>F;
        ifile>>filling;
        ifile>>temperatura;
        ifile>>nu_0;
        ifile>>n_class;
        ifile>>total_time;

        std::string bool_tracker;
        ifile>>bool_tracker;
        if(bool_tracker=="T") output_file=true;
        else if(bool_tracker=="F") output_file=false;
        ifile>>bool_tracker;
        if(bool_tracker=="T") fill_first_layer=true;
        else if(bool_tracker=="F") fill_first_layer=false;

        for (int i=0; i<100; i++) movements_per_class[i]=0;
        for (int i=0 ; i< 100; i++){ for(int j=0 ; j<1000 ; j++){for (int k=0; k<2; k++) map_class_processes[i][j][k]=-1;}}
        P.resize(n_class, 0.0);
        n_deposited=0;
    }
}

/*
 * Identifies the lattice sites that require a recalculation of their 
 * possible processes after a move or deposition event. It optimizes performance 
 * by targeting only the relevant neighborhood .
 */
void metropolis::interested_sites_calc(bool dep){
    auto nn_pos = s.get_table_of_nn(pos);
    interested_sites.clear();
    interested_sites.push_back(pos);
    for(size_t i=0 ; i<nn_pos.size() ; i++){
        if(nn_pos[i]!=pos) interested_sites.push_back(nn_pos[i]);
    }
    if(!dep){
        auto nn_next = s.get_table_of_nn(next);
        for(size_t i=0 ; i<nn_next.size() ; i++){
            bool already_here=false;
            for(size_t l=0 ; l< interested_sites.size() ; l++) if(nn_next[i]==interested_sites[l]) { already_here=true; break; }
            if(nn_next[i]!=pos && nn_next[i]!=next && !already_here) interested_sites.push_back(nn_next[i]);
        }
    }

    if (dropped_from != -1) {
        auto nn_drop = s.get_table_of_nn(dropped_from);
        for(size_t i = 0; i < nn_drop.size(); i++){
            bool already_here = false;
            for(size_t l = 0; l < interested_sites.size(); l++) {
                if(nn_drop[i] == interested_sites[l]) { already_here = true; break; }
            }
            if(nn_drop[i] != pos && nn_drop[i] != next && !already_here) {
                interested_sites.push_back(nn_drop[i]);
            }
        }
        interested_sites.push_back(dropped_from);
    }

    if (dropped_target != -1) {
        auto nn_tgt = s.get_table_of_nn(dropped_target);
        for(size_t i=0 ; i<nn_tgt.size() ; i++){
            bool already_here=false;
            for(size_t l=0 ; l< interested_sites.size() ; l++) if(nn_tgt[i]==interested_sites[l]) { already_here=true; break; }
            if(nn_tgt[i]!=pos && nn_tgt[i]!=next && !already_here) interested_sites.push_back(nn_tgt[i]);
        }
        interested_sites.push_back(dropped_target);
    }
}

/*
 * Updates the number of nearest neighbors (coordination) for sites 
 * affected by the recent atomic movement. This count is fundamental 
 * for calculating transition energy barriers .
 */
void metropolis::nn_updater(bool dep){
    bool pos_checked=false;
    bool next_checked=false;
    std::unordered_set<int> sites_to_update;

    // Neighbors of the departure position
    auto nn_pos = s.get_table_of_nn(pos);
    sites_to_update.insert(nn_pos.begin(), nn_pos.end());

    // If diffusion, also add neighbors of the arrival position
    if(!dep){
        auto nn_next = s.get_table_of_nn(next);
        sites_to_update.insert(nn_next.begin(), nn_next.end());
    }

    //id dropped
    if (dropped_from != -1) {
        auto nn_drop = s.get_table_of_nn(dropped_from);
        sites_to_update.insert(nn_drop.begin(), nn_drop.end());
        sites_to_update.insert(dropped_from);
    }

    if (dropped_target != -1) {
        auto nn_tgt = s.get_table_of_nn(dropped_target);
        sites_to_update.insert(nn_tgt.begin(), nn_tgt.end());
        sites_to_update.insert(dropped_target);
    }

    // Update each site coordination precisely once
    for(int site : sites_to_update){
        
        int count = 0;
        auto nn_of_nn = s.get_table_of_nn(site);
        auto edge = s.get_edge_map(site);

        for(size_t i=0 ; i< nn_of_nn.size() ; i++){
            int nearest=nn_of_nn[i];
            int edges= edge[i];
                if(atoms[nearest] && (edges < 0)){ 
                    if(nearest == pos && !pos_checked){ 
                        pos_checked = true;
                        count++;
                        continue;
                    }else if(nearest == pos && pos_checked){
                        continue;
                    }
                    if(!dep && nearest == next && !next_checked){
                        next_checked = true;
                        count++;
                        continue;
                    }else if(!dep && nearest == next && next_checked){
                        continue;
                    }
                    count++;
                }
        }
        if(s.get_TOPl(site)==14 || s.get_TOPl(site)==15) count = count - 4;
        if(s.get_TOPl(site)==10) count = count -1;
        nnn_atoms[site] = count;
    }
}

/*
 * Determines if a site index belongs to a facet boundary (border) based 
 * on the geometric mapping (TOPl values) .
 */
bool metropolis::is_border(int i){
    if(s.get_TOPl(i)==8 || s.get_TOPl(i)==7 || s.get_TOPl(i)==15){
        return true;
    }
    return false;
}

/*
 * Checks if two border sites belong to different facets of the nanoparticle .
 */
bool metropolis::are_different_border(int i, int j){
    if(is_border(i) && is_border(j)){
        if(s.get_TOPl(i)!=s.get_TOPl(j)){
            return true;
        }
    }
    return false;
}

/*
 * Populates the static tables of possible physical processes based on 
 * the crystal geometry. It identifies all valid jumps between centers 
 * and borders of different facets .
 */
void metropolis::table_of_processes_filler(){
    table_of_processes.resize(crd.sites_size(), std::vector<int>(0,0));
    table_of_end_pos.resize(crd.sites_size(), std::vector<int>(0,0));
    table_of_initial_pos.resize(crd.sites_size(), std::vector<int>(0,0));
    
    for(long unsigned int i=0 ; i<crd.sites_size() ; i++){
        std::vector<int> nns= s.get_table_of_nn(i);
        for(size_t j=0 ; j<nns.size() ; j++){
            int pv=nns[j];
            if(s.get_TOPl(i)==1 || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)==s.get_plane(pv,1))){
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 0){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(0); // 111c to 111c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            else if((s.get_TOPl(i)==0 && s.get_TOPl(pv)==0) || (s.get_TOPl(i)==7 && s.get_TOPl(pv)==7)) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 36){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(36); // 100c to 100c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)!=s.get_plane(pv,1)) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 15){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(15); // 111 inter-facet jump
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==7) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 29){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(29); // 111 jump to 100
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            else if(s.get_TOPl(i)==7 && s.get_TOPl(pv)==8) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 22){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(22); // 100 jump to 111
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            else if((s.get_TOPl(i)==7 && s.get_TOPl(pv)==0)) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 45){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(45); // 100b to 100c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==1) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 8){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(8); // 111b to 111c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            else if((s.get_TOPl(i)==0 && s.get_TOPl(pv)==7)) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 53){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(53); // 100c to 100b
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
        }
    }
}

/*
 * Dynamically updates the process table during the simulation. It handles 
 * the activation of new sites (like "mixed" sites) that appear as the 
 * cluster grows .
 */
void metropolis::table_of_processes_updater(int i){    
    std::vector<int> nns= s.get_table_of_nn(i);
    for(size_t j=0 ; j<nns.size() ; j++){
        int pv=nns[j];
        if(s.get_TOPl(i)==1 || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)==s.get_plane(pv,1)) || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==10) || (s.get_TOPl(i)==10 && s.get_TOPl(pv)==10) || (s.get_TOPl(i)==10 && s.get_TOPl(pv)==8) ){
            bool already_present = false;
            for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 0){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(0); 
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
            }
        }
        else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)!=s.get_plane(pv,1)) {
            bool already_present = false;
            for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 15){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(15); 
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
            }
        }
        else if((s.get_TOPl(i)==8 || s.get_TOPl(i)==10) && s.get_TOPl(pv)==7) {
            bool already_present = false;
            for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 29){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(29); 
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
            }
        }
        else if(s.get_TOPl(i)==7 && (s.get_TOPl(pv)==8 || s.get_TOPl(pv)==10)) {
            bool already_present = false;
            for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 22){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(22); 
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
            }
        }
        else if((s.get_TOPl(i)==8 || s.get_TOPl(i)==10) && s.get_TOPl(pv)==1) {
            bool already_present = false;
            for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 8){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(8); 
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
            }
        }
        else if(s.get_TOPl(i)==10 && s.get_TOPl(pv)>10){ 
            bool already_present = false;
            for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 66){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(66);
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
            }
        }
    }
}

/*
 * Checks if a site is currently in the "deactivated" state. This usually 
 * applies to the 100-facet second layer sites before they are 
 * geometrically reachable .
 */
bool metropolis::is_deactivated(int i){
    for(int j : deactivated_sites) if (j==i) return true;
    return false;
}

/*
 * Activates the geometric edges between facets. When an atom is nearby, 
 * it "unlocks" previously hidden coordination sites to allow inter-facet diffusion .
 */
std::array<int,2> metropolis::activate_edges(int upper_site, int pv){
    std::array<int,2> back={-1,-1};
    auto edge = s.get_edge_map(pv);
    int index=0;
    for(size_t i=0 ; i<edge.size() ; i++){
        if(edge[i]>=0){
            s.set_table_of_nn(pv , i, edge[i]);
            s.append_nn(edge[i], upper_site);
            s.set_TOPl(edge[i],10);  
            for(auto old_nn : s.get_table_of_nn(edge[i])){
                if(s.get_TOPl(old_nn)==8){
                    int oppc = 0;
                    for (int p : table_of_processes[old_nn]){
                        if(p==0) oppc++;
                    }
                    s.append_nn(old_nn, edge[i]);
                    if(oppc > 2) s.set_TOPl(old_nn, 1);
                    s.set_plane(edge[i], 20);
                    table_of_processes_updater(old_nn); 
                }
            }
            table_of_processes_updater(edge[i]);
            s.reset_edge_map(pv, i);

            back[index]=edge[i];
            index++;
        }
        if(index>1) return back;
    }
    return back;
}

/*
 * Specific update for the second layer of the nanoparticle. It defines how 
 * atoms on newly created surfaces can move based on local coordination .
 */
void metropolis::second_layer_updates(int upper_site){
    auto nn=s.get_table_of_nn(upper_site);
    for (int pv : nn){
        if(s.get_TOPl(upper_site)==s.get_TOPl(pv)){
            if(is_deactivated(pv)==false){
                //std::cout<< "process " << 36 << " from "<< upper_site << " to " << pv <<std::endl; 
                table_of_processes[upper_site].push_back(36); 
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
            }
        }
        else if((s.get_TOPl(upper_site)==14 && s.get_TOPl(pv)==15) ){
            if(is_deactivated(pv)==false){
                //std::cout<< "process " << 53 << " from "<< upper_site << " to " << pv <<std::endl; 
                table_of_processes[upper_site].push_back(53); 
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
            }
        }
        else if((s.get_TOPl(upper_site)==15 && s.get_TOPl(pv)==14)){
            if(is_deactivated(pv)==false){
                //std::cout<< "process " << 45 << " from "<< upper_site << " to " << pv <<std::endl; 
                table_of_processes[upper_site].push_back(45); 
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
            }
        }
        else if(s.get_TOPl(pv)==0 || s.get_TOPl(pv)==7){ 
            table_of_processes[upper_site].push_back(72);
            table_of_end_pos[upper_site].push_back(pv);
            table_of_initial_pos[pv].push_back(upper_site);
        }
        if(s.get_TOPl(pv)==7){ 
            auto edge_sites = activate_edges(upper_site,pv);
            for (auto edge_site : edge_sites){
                if(edge_site>=0){
                    table_of_processes[upper_site].push_back(22);
                    table_of_end_pos[upper_site].push_back(edge_site);
                    table_of_initial_pos[edge_site].push_back(upper_site);
                }
            }
        }
    }
}

/*
 * Controls the activation of the second atomic layer. When a site 
 * reaches a coordination of 4, it is no longer considered "virtual" 
 * and enters the simulation dynamics .
 */
 void metropolis::second_layer_activation(){
    for (size_t i=0 ; i<deactivated_sites.size() ; i++){
        int counter = 0;
        int upper_site = deactivated_sites[i];
        for(int j : s.get_table_of_nn(upper_site)){
            if(atoms[j]==true) counter++;
        }
        if(counter >= 4){
            second_layer_updates(upper_site); 
            deactivated_sites.erase(deactivated_sites.begin()+i);
        }
    }
}

void metropolis::second_layer_deactivation(int deserted_site) {
    auto nn = s.get_table_of_nn(deserted_site);
    bool fallen_down = false;
    for (int upper_site : nn) {
        // Verifica se è un sito del secondo strato (assumendo TOPl 14 o 15)
        if (s.get_TOPl(upper_site) == 14 || s.get_TOPl(upper_site) == 15) {
            // Se il sito è correntemente attivo
            if (!is_deactivated(upper_site)) {
                
                // Ricalcola il supporto
                int counter = 0;
                for (int j : s.get_table_of_nn(upper_site)) {
                    if (atoms[j]) counter++;
                }
                
                // Condizione non più fisica: meno di 4 supporti
                if (counter < 4) {
                    
                    // Se il sito è occupato, l'atomo cade giù
                    if (atoms[upper_site]) {
                        // Cerca un sito di supporto vuoto (almeno deserted_site lo sarà)
                        atoms[upper_site] = false;
                        atoms[deserted_site] = true;
                        fallen_down = true;
                        
                        if (output_file) {
                            output << "Drop down: da " << upper_site << " a " << deserted_site 
                                   << " al tempo: " << time << std::endl;
                        }
//                        std::cout << "Drop down: da " << upper_site << " a " << deserted_site 
//                                   << " al tempo: " << time << "con processo "<< pos << " to "<< next<< std::endl;

                        // Salva traccia per il ricalcolo degli interessati
                        dropped_from = upper_site;
                        dropped_target = deserted_site;
                    }
                    
                    // Pulizia totale dei processi per il sito disattivato
                    //for (int t = 0; t <= 12; t++) map_of_class_position_eraser(upper_site, t);
                    //map_of_class_next_eraser(dropped_target);
                    //map_of_class_next_eraser(upper_site);
                    
                    // Rimuoviamo il sito dalle tabelle statiche dei vicini
                    for (int pv : table_of_end_pos[upper_site]) {
                        auto& init_pos = table_of_initial_pos[pv];
                        auto it = init_pos.begin();
                        while (it != init_pos.end()) {
                            if (*it == upper_site) it = init_pos.erase(it);
                            else ++it;
                        }
                    }
                    table_of_processes[upper_site].clear();
                    table_of_end_pos[upper_site].clear();
                }
            }
        }
        if(fallen_down) break;
    }
}

/*
void metropolis::second_layer_activation() {
    std::vector<int> to_remove;

    for (size_t i = 0; i < deactivated_sites.size(); i++) {
        int counter = 0;
        int upper_site = deactivated_sites[i];
        for (int j : s.get_table_of_nn(upper_site)) {
            if (atoms[j] == true) counter++;
        }
        if (counter >= 4) {
            second_layer_updates(upper_site);
            to_remove.push_back(upper_site); 
        }
    }

    //re-calculate all the processes if more than 1 site have been activated
    for (auto site : to_remove) {
        second_layer_updates(site);
    }

    deactivated_sites.erase(
        std::remove_if(deactivated_sites.begin(), deactivated_sites.end(),
            [&](int site) {
                return std::find(to_remove.begin(), to_remove.end(), site) != to_remove.end();
            }),
        deactivated_sites.end()
    );
}
*/

/*
 * Maps a process class to the correct energy barrier ID and neighbor multiplier. 
 * This is used to compute the Arrhenius probability of the event .
 */
void metropolis::barrier_chooser(int cl, int &bar_id , int &neighbours ){
    std::string num;
    int bar[2];

    if (cl >= 8 && cl <= 14) { bar[0] = 6; bar[1] = cl - 8; }
    else if (cl >= 15 && cl <= 21) { bar[0] = 2; bar[1] = cl - 15; }
    else if (cl >= 22 && cl <= 28) { bar[0] = 4; bar[1] = cl - 22; }
    else if (cl >= 29 && cl <= 35) { bar[0] = 3; bar[1] = cl - 29; }
    else if (cl >= 36 && cl <= 44) { bar[0] = 1; bar[1] = cl - 36; }
    else if (cl >= 45 && cl <= 52) { bar[0] = 5; bar[1] = cl - 45; }
    else if (cl >= 53 && cl <= 59) { bar[0] = 7; bar[1] = cl - 53; }
    else if (cl >= 60 && cl <= 65) { bar[0] = 8; bar[1] = cl - 60; }
    else if (cl >= 66 && cl <= 71) { bar[0] = 9; bar[1] = cl - 66; }
    else if (cl >= 72) { bar[0] = 10; bar[1] = cl - 72; }
    else { bar[0] = 0; bar[1] = cl; }

    bar_id=bar[0];
    neighbours=bar[1];
}

/*
 * Assigns a string name to each process class based on the geometric 
 * transition (e.g., 111border to 111center) for readable output logging .
 */
void metropolis::get_process_name(int cl){
    std::string num;
    int nn;

    if (cl >= 8 && cl <=14)
    {
        nn = cl - 8; num=std::to_string(nn);
        process_name = "111border to 111center same facet with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 15 && cl <= 21)
    {
        nn = cl - 15; num=std::to_string(nn);
        process_name = "111border to 111border different facets with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 22 && cl <= 28)
    {
        nn = cl - 22; num=std::to_string(nn);
        process_name = "100 to 111 with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 29 && cl <= 35)
    {
        nn = cl - 29; num=std::to_string(nn);
        process_name = "111 to 100 with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 36 && cl <= 44)
    {
        nn = cl - 36; num=std::to_string(nn);
        process_name = "100center to 100center same facet with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 45 && cl <= 52)
    {
        nn = cl - 45; num=std::to_string(nn);
        process_name = "100border to 100center same facet with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 53 && cl <= 59)
    {
        nn = cl - 53; num=std::to_string(nn);
        process_name = "100center to 100border with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 60 && cl <= 65)
    {
        nn = cl - 60; num=std::to_string(nn);
        process_name = "111border-push-100  with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 66 && cl <= 71)
    {
        nn = cl - 66; num=std::to_string(nn);
        process_name = "111 mixed to 100 double with: " + num + " nearest neighbors";
        return;
    }
    else if (cl >= 72)
    {
        nn = cl - 72; num=std::to_string(nn);
        process_name = "kink: 100shell to 100core with: " + num + " nearest neighbors";
        return;
    }
    else
    {
        nn = cl; num=std::to_string(nn);
        process_name = "111center to 111center same facet with: " + num + " nearest neighbors";
        return;
    }
}

/*
 * Calculates the jump frequency (probability rate) using the Arrhenius law: 
 * $Rate = \nu_0 \exp\left(-\frac{E_{barrier} + n \cdot E_b}{K_b T}\right)$ .
 */
double metropolis::probability_calculator(int cl){
    int b_id = 0;
    int neighbours = 0;
    barrier_chooser(cl, b_id, neighbours);
    double rt_value = nu_0 * (exp((-1.) *(barriere[b_id] + (neighbours*E_b)) / (Kb*temperatura) )) ;
    return (rt_value);
}

/*
 * Maintenance function that shifts the processes list to ensure there 
 * are no "holes" in the memory structure after an event .
 */
void metropolis::shifter(int c, int ind, int tg){
    for(int i = ind ; i<tg ; i++){
        map_class_processes[c][i][0] = map_class_processes[c][i+1][0];
        map_class_processes[c][i][1] = map_class_processes[c][i+1][1];
    }
}

/*
 * Removes a single specific jump from the process catalog when it 
 * is no longer physically possible (e.g., neighbor coordination changed) .
 */
void metropolis::map_of_class_eraser(int c, int i, int j) {
    auto& processes= map_class_processes[c];
    int tg = get_MCP_size(c);
    for(int ind=0; ind< tg ; ind++){
        auto& couple = processes[ind];
        if(couple[0]==-1) break;
        if(couple[0]==i && couple[1]==j){
            couple[0]=-1;
            couple[1]=-1;
            movements_per_class[c]--;
            shifter(c,ind,tg);
            break;
        }
    }
}

/*
 * Erases all processes originating from a specific site. Used when 
 * an atom moves away from that site, making all previous outgoing 
 * jumps invalid .
 */
void metropolis::map_of_class_position_eraser(int i, int n_nn) {
    auto processes = table_of_processes[i];
    for(size_t prc=0 ; prc<processes.size() ; prc++){
        int end_p = table_of_end_pos[i][prc];
        map_of_class_eraser(processes[prc] + n_nn ,i,end_p);
    }
}

/*
 * Erases all processes targeting a specific site. Used when a site 
 * becomes occupied, making it an invalid destination for other atoms .
 */
void metropolis::map_of_class_next_eraser(int i) {
    auto copy_tip=table_of_initial_pos[i];
    for(int i_site : copy_tip){
        for(int n=0; n<= 12; n++){
            for(size_t i_prc=0 ; i_prc < table_of_end_pos[i_site].size() ; i_prc++) if(table_of_end_pos[i_site][i_prc]==i){
                map_of_class_eraser(table_of_processes[i_site][i_prc] + n, i_site , i);
            }
        }
    }
}

/*
 * Adds a physical jump to the global catalog for a given class .
 */
void metropolis::map_of_class_filler(int c, int i , int j){
    auto& processes = map_class_processes[c];
    for (auto& couple: processes){
        if(couple[0]<0){
            couple[0]=i;
            couple[1]=j;
            movements_per_class[c]++;
            break;
        }
    }
}

/*
 * Returns the current number of valid processes stored in a specific class .
 */
int metropolis::get_MCP_size(int c){
    int counter=0;
    auto& processes = map_class_processes[c];
    for (auto& couple: processes){
        if(couple[0]<0){
            break;
        }
        counter++;
    }
    return counter;
}

/*
 * Clears conditional processes like "Push" or "Kink" which depend on 
 * transient local configurations and must be rerealculated every step .
 */
void metropolis::clear_dynamic_processes(){
    for(auto& [site, pv, cl] : dynamic_processes){
        map_of_class_eraser(cl, site, pv);
    }
    dynamic_processes.clear();
}

/*
 * The core logic for Kinetic Monte Carlo classification. For each 
 * interested site, it checks which physical moves are valid based 
 * on current occupancy and coordination .
 */
void metropolis::classification(){
    clear_dynamic_processes();
    for(int site : interested_sites){
        if(atoms[site]){
            int nn_count = nnn_atoms[site];
            
            for(int t=0 ; t <= 12 ; t++) map_of_class_position_eraser(site, t);

            std::vector<int>& top = table_of_processes[site];

            for(size_t j=0 ; j < top.size() ; j++){  
                int base_process = top[j]; 
                int pv = table_of_end_pos[site][j];
                
                if(!atoms[pv]){ 
                    int actual_class = base_process + nn_count; 
                    if(actual_class >= 0 && actual_class < 100) {
                        //if(actual_class>35 && actual_class < 60) std::cout<< "process " << actual_class << "site    " << site<< " neighbors "<<nn_count << std::endl;
                        map_of_class_filler(actual_class, site, pv);                        
                    } else {
                        std::cerr << "Error: Class " << actual_class << " exceeds array size!" << std::endl;
                    }
                }
                // Push logic
                else if(atoms[pv] && s.get_TOPl(pv)==7 && s.get_TOPl(site)==8){
                    auto nn = s.get_table_of_nn(pv);
                    for(int n : nn){
                        if(!atoms[n] && s.get_TOPl(n)==0){
                            int actual_class = 60 + nn_count;
                            dynamic_processes.emplace_back(site, n, actual_class);
                            map_of_class_filler(actual_class, site, n);
                        }
                    }
                }
                // Kink logic
                else if(atoms[pv] && base_process==72){
                    auto nn = s.get_table_of_nn(pv);
                    for(int n : nn){
                        if(!atoms[n] && (s.get_TOPl(n)==0 || s.get_TOPl(n)==7)){
                            int actual_class = 72 + nn_count;
                            dynamic_processes.emplace_back(site, n, actual_class);
                            map_of_class_filler(actual_class, site, n);
                        }
                    }
                }
            }
        }
    }
}

/*
 * Calculates the total probability for each class by multiplying the 
 * individual process rate by the number of active processes in that class .
 */
void metropolis::probability_filler(){
    std::fill(P.begin(), P.end(), 0.0); 
    for(int i=0 ; i < n_class ; i++){
        if(movements_per_class[i] > 0) {
            double p = probability_calculator(i);
            if(movements_per_class[i]> 6*n_deposited) std::exit(0);
            P[i] = movements_per_class[i] * p;
        }
    }
}

/*
 * Handles the deposition of a new atom. It randomly selects a valid 
 * surface site on the nanoparticle and updates the system occupancy .
 */
void metropolis::deposition_func(){
    int p; 
    int n_nodep = deactivated_sites.size();
    int max_sites = crd.sites_size();

    while(true) {
        p = static_cast<int>(initial_pos(rng) * (max_sites - n_nodep));
        if(p < 0) p = 0;
        if(p >= max_sites) p = max_sites - 1;
        if(s.get_TOPl(p) < 2 && !atoms[p]) {
            break; 
        }
    }
    atoms[p] = true;
    pos = p;
    n_deposited++;
    output << "New atom: " << p << "    at time: "<< time <<std::endl;
}

/*
 * Calculates the time interval between events and chooses the next 
 * event to occur (Deposition vs Diffusion) using a weighted 
 * random selection based on total probabilities .
 */
void metropolis::time_prob_calc(){
    double r_time=dist_time(rng);
    double P_tot=0.0;
    double P_dep=crd.sites_size()*F;
    for(size_t i=0 ; i<P.size() ; i++) P_tot+=P[i];
    double P_diff=P_tot;
    P_tot+=P_dep;
    
    // Time increment based on the sum of all possible event rates
    time += ( -(log(r_time)) / P_tot);
    if(time>total_time) {
        std::cerr<<" Too much long time increment: "<<time<<" with a total available time for the simulation of: "<<total_time<<std::endl;
        print_final_configuration();
        std::exit(0);
    }

    double random = dist_time(rng);
    double r_p = random * P_tot;
    if(r_p>P_diff){
        deposition=true;
    } else {
        deposition=false;
        int chosen_class=-1;
        double p_i = 0;
        for(size_t i = 0 ; i < P.size() ; i++){
            if( (p_i < r_p) && (r_p < (p_i+P[i])) ) chosen_class = i;
            p_i += P[i];
        }

        int mcount = movements_per_class[chosen_class];
        if(mcount <= 0){
            std::cerr << "Warning: chosen_class has no movements: " << chosen_class << std::endl;
            std::exit(0);
        }
        if(chosen_class<0){
            std::cerr << "Warning: chosen_class has no size "<<std::endl;
            std::exit(0);
        }
        std::uniform_int_distribution<> dist(0, mcount - 1);
        int movement = dist(rng);
        if(atoms[map_class_processes[chosen_class][movement][0]]==false){
            std::cerr<<"Atomo inesistente...  classe "<<chosen_class << " movement "<<movement<<std::endl;
            print_configuration();
            /*
            for (size_t i=0 ; i<100 ; i++){
                bool print=false;
                for (size_t j=0 ; j<1000 ; j++){
                    for (size_t t=0 ; t<2 ; t++){
                        if (map_class_processes[i][j][t] > 0 ){
                             std::cout<< "class "<< i <<" movement  "<<j<<" : "<< map_class_processes[i][j][t]<<" ";
                             print=true;
                        }
                    }
                    if(print) std::cout<<"        ";
                }
                if (print) std::cout<<std::endl;
            }
            */

            std::exit(0);
        }
        if(movement < 0 || movement >= get_MCP_size(chosen_class)){
            std::cerr << "Error: movement index out of range: " << movement << " >= " << get_MCP_size(chosen_class) << std::endl;
            return;
        }
        const auto& mv = map_class_processes[chosen_class][movement];
        next = mv[1];
        pos = mv[0];
        get_process_name(chosen_class);
        if(output_file) print_output();
        
        for(int t=0 ; t <= 12 ; t++){
            map_of_class_position_eraser(pos, t);
        }
        map_of_class_next_eraser(next);
    }
}

/*
 * Logs the specific process that just occurred to the output file .
 */
void metropolis::print_output(){
    output<< process_name << "  from    " << pos <<"  to    "<< next <<std::endl;
}

/*
 * Exports the current spatial configuration of the nanoparticle into 
 * a file for visualization and later analysis .
 */
void metropolis::print_configuration(){
    std::string str1 = std::to_string(n_deposited);
    int numZeros = std::max(0, 4 - static_cast<int>(str1.length()));
    str1.insert(0, numZeros, '0');
    
    std::string fileName = "occ" + str1 + ".out";
    std::ofstream intermed_config(fileName);
    
    c_crd.output_writer(intermed_config, c_crd.sites_size() + n_deposited, time);
    crd.output_writer_partial(intermed_config, n_deposited, atoms);
    
    if(output_file) output << "Printed configuration with: " << n_deposited 
           << " occupated sites, at time: " << time * pow(10, 12) 
           << " picoseconds, = " << time << " seconds " << std::endl;
}

void metropolis::print_final_configuration(){
    std::string str1 = std::to_string(n_deposited);
    int numZeros = std::max(0, 4 - static_cast<int>(str1.length()));
    str1.insert(0, numZeros, '0');
    
    std::string fileName = "occ" + str1 + "FINAL.out";
    std::ofstream intermed_config(fileName);
    
    c_crd.output_writer(intermed_config, c_crd.sites_size() + n_deposited, time);
    crd.output_writer_partial(intermed_config, n_deposited, atoms);
    
    if(output_file) output << "Printed configuration with: " << n_deposited 
           << " occupated sites, at time: " << time * pow(10, 12) 
           << " picoseconds, = " << time << " seconds " << std::endl;
}

void metropolis::first_layer_filler(){
    for(int i=0 ; i<crd.get_N() ; i++){
        if(s.get_TOPl(i)==0 || s.get_TOPl(i)==7){
            atoms[i] = true;
            n_deposited++;
        }
        second_layer_activation();
    }
}


/*
 * Sets up the initial simulation state, initializes the nanoparticle, 
 * and triggers the first deposition event .
 */
void metropolis::start_of_the_sim(){
    if(output_file) output.open("MC_processes.out");
    table_of_processes_filler(); 
    deactivated_sites=s.get_upper_sites();
    atoms.resize(crd.get_N(), false);
    if(fill_first_layer) first_layer_filler();
    deposition_func();
    interested_sites_calc(true);
    nnn_atoms.resize(atoms.size());
    nn_updater(deposition);
}


/*
 * The main Kinetic Monte Carlo loop. It iterates through process 
 * classification, probability calculation, and event execution 
 * until the total simulation time is reached .
 */
void metropolis::algorithm(){
    while(time < total_time){
        dropped_from = -1;
        dropped_target = -1;
        
        classification();
        probability_filler();
        time_prob_calc();
        
        if(deposition){
            deposition_func();
            print_configuration();
        } else {
            atoms[pos] = false;
            atoms[next]= true; 
        }
        
        second_layer_deactivation(pos);
        second_layer_activation();
        interested_sites_calc(deposition);
        nn_updater(deposition);
        pos = next;
        
        if(n_deposited == (int)(filling * (float) crd.sites_size())) {
            F = 0; 
            if(output_file) output << "End of the deposition, there are " << n_deposited << " filled sites " << std::endl;
        }
    }
}

/*
 * The public entry point to start the entire simulation workflow .
 */
void metropolis::simulation(){
    start_of_the_sim();
    algorithm();
    if(output_file) print_output();
    print_final_configuration();
}