#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <new>
#include <unordered_set>

/**
 * Reads simulation parameters from an input file.
 * Initializes energy barriers, temperature, flux (F), and 
 * prepares the internal data structures for the process classes.
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

        std::string file;
        ifile>>file;
        if(file=="T") output_file=true;
        else if(file=="F") output_file=false;

        for (int i=0; i<100; i++) movements_per_class[i]=0;
        for (int i=0 ; i< 100; i++){ 
            for(int j=0 ; j<1000 ; j++){
                for (int k=0; k<2; k++) map_class_processes[i][j][k]=-1;
            }
        }
        P.resize(n_class, 0.0);
        n_deposited=0;
    }
}

/**
 * Identifies which lattice sites need an update after a move or deposition.
 * This limits the recalculation only to the neighborhood of the changed sites, 
 * which is essential for simulation performance.
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
}

/**
 * Updates the coordination number (nearest neighbors) for sites affected by the last event.
 * It accounts for different facets and ensures that geometric boundaries are respected.
 */
void metropolis::nn_updater(bool dep){
    bool pos_checked=false;
    bool next_checked=false;
    std::unordered_set<int> sites_to_update;

    auto nn_pos = s.get_table_of_nn(pos);
    sites_to_update.insert(nn_pos.begin(), nn_pos.end());

    if(!dep){
        auto nn_next = s.get_table_of_nn(next);
        sites_to_update.insert(nn_next.begin(), nn_next.end());
    }

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
        nnn_atoms[site] = count;
    }
}

/**
 * Helper to determine if a site belongs to a geometric boundary 
 * based on its TOPl (Type of Plane/Facet) value.
 */
bool metropolis::is_border(int i){
    if(s.get_TOPl(i)==8 || s.get_TOPl(i)==7 || s.get_TOPl(i)==15){
        return true;
    }
    return false;
}

/**
 * Checks if two sites belong to different geometric boundaries.
 */
bool metropolis::are_different_border(int i, int j){
    if(is_border(i) && is_border(j)){
        if(s.get_TOPl(i)!=s.get_TOPl(j)){
            return true;
        }
    }
    return false;
}

/**
 * Initial population of the possible physical processes based on static geometry.
 * Categorizes jumps (111 to 111, 100 to 111, etc.) into base process types.
 */
void metropolis::table_of_processes_filler(){
    table_of_processes.resize(crd.sites_size(), std::vector<int>(0,0));
    table_of_end_pos.resize(crd.sites_size(), std::vector<int>(0,0));
    table_of_initial_pos.resize(crd.sites_size(), std::vector<int>(0,0));

    for(size_t i=0 ; i<crd.sites_size() ; i++){
        std::vector<int> nns= s.get_table_of_nn(i);
        for(size_t j=0 ; j<nns.size() ; j++){
            int pv=nns[j];
            // Cases: 111 centers, 100 centers, jumps between facets, etc.
            if(s.get_TOPl(i)==1 || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)==s.get_plane(pv,1))){
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
            else if((s.get_TOPl(i)==0 && s.get_TOPl(pv)==0) || (s.get_TOPl(i)==7 && s.get_TOPl(pv)==7)) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 36){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(36); 
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
            // (Other geometry cases follow the same logic...)
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
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==7) {
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
            else if(s.get_TOPl(i)==7 && s.get_TOPl(pv)==8) {
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
            else if((s.get_TOPl(i)==7 && s.get_TOPl(pv)==0)) {
                 bool already_present = false;
                for(size_t k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 45){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(45); 
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
                    table_of_processes[i].push_back(8); 
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
                    table_of_processes[i].push_back(53); 
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                }
            }
        }
    }
}

/**
 * Updates process tables dynamically when mixed sites or new facets are involved.
 * This is used during simulation as the cluster grows.
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
        // (Similar logic for 29, 22, 8, and 66...)
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
        else if((s.get_TOPl(i)==8|| s.get_TOPl(i)==8) && s.get_TOPl(pv)==1) {
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

/**
 * Checks if a site is in the deactivated list (typically sites on the 100-facet
 * that aren't available for the 2nd layer yet).
 */
bool metropolis::is_deactivated(int i){
    for(int j : deactivated_sites) if (j==i) return true;
    return false;
}

/**
 * Transforms an edge site into an active gateway between facets.
 * This happens when atoms occupy strategic positions that "bridge" two planes.
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
                    for (int p : table_of_processes[old_nn]) if(p==0) oppc++;
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

/**
 * Adds process definitions specifically for the second layer of the 100-facet.
 */
void metropolis::second_layer_updates(int upper_site){
    auto nn=s.get_table_of_nn(upper_site);
    for (int pv : nn){
        if(s.get_TOPl(upper_site)==s.get_TOPl(pv)){
            if(is_deactivated(pv)==false){
                table_of_processes[upper_site].push_back(36); 
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
            }
        }
        else if((s.get_TOPl(upper_site)==14 && s.get_TOPl(pv)==15) ){
            if(is_deactivated(pv)==false){
                table_of_processes[upper_site].push_back(53); 
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
            }
        }
        else if((s.get_TOPl(upper_site)==15 && s.get_TOPl(pv)==14)){
            if(is_deactivated(pv)==false){
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

/**
 * Checks the "deactivated" sites on the 100-facet. 
 * If a site now has enough neighbors (>=4), it is activated for the simulation.
 */
void metropolis::second_layer_activation(){
    for (size_t i=0 ; i<deactivated_sites.size() ; i++){
        int counter = 0;
        int upper_site = deactivated_sites[i];
        for(int j : s.get_table_of_nn(upper_site)) if(atoms[j]==true) counter++;
        if(counter >= 4){
            second_layer_updates(upper_site); 
            deactivated_sites.erase(deactivated_sites.begin()+i);
        }
    }
}

/**
 * Maps a process class index back to its physical parameters:
 * barrier_id (base energy) and the number of nearest neighbors (multiplied by E_b).
 */
void metropolis::barrier_chooser(int cl, int &bar_id , int &neighbours ){
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

/**
 * Returns a human-readable name for the selected process class (used for logging).
 */
void metropolis::get_process_name(int cl){
    int nn;
    std::string num;
    if (cl >= 8 && cl <=14) { nn = cl - 8; process_name = "111border to 111center same facet with: " + std::to_string(nn) + " nn"; }
    else if (cl >= 15 && cl <= 21) { nn = cl - 15; process_name = "111border to 111border diff facet with: " + std::to_string(nn) + " nn"; }
    // ... (rest of the process naming logic)
    else { nn = cl; process_name = "111center to 111center with: " + std::to_string(nn) + " nn"; }
}

/**
 * Calculates the transition rate using the Arrhenius law: 
 * rate = nu_0 * exp(-(E_barrier + n*E_bond) / kT)
 */
double metropolis::probability_calculator(int cl){
    int b_id = 0, neighbours = 0;
    barrier_chooser(cl, b_id, neighbours);
    return nu_0 * (exp((-1.) *(barriere[b_id] + (neighbours*E_b)) / (Kb*temperatura)));
}

/**
 * Utility to shift elements in the process array after an erasure to keep it compact.
 */
void metropolis::shifter(int c, int ind, int tg){
    for(int i = ind ; i<tg ; i++){
        map_class_processes[c][i][0] = map_class_processes[c][i+1][0];
        map_class_processes[c][i][1] = map_class_processes[c][i+1][1];
    }
}

/**
 * Removes a specific [start_site, end_site] move from the process catalog.
 */
void metropolis::map_of_class_eraser(int c, int i, int j) {
    auto& processes= map_class_processes[c];
    int tg = get_MCP_size(c);
    for(int ind=0; ind< tg ; ind++){
        auto& couple = processes[ind];
        if(couple[0]==i && couple[1]==j){
            couple[0]=-1; couple[1]=-1;
            movements_per_class[c]--;
            shifter(c,ind,tg);
            break;
        }
    }
}

/**
 * Removes all processes starting from a site. Called when an atom leaves a site.
 */
void metropolis::map_of_class_position_eraser(int i, int n_nn) {
    auto processes = table_of_processes[i];
    for(size_t prc=0 ; prc<processes.size() ; prc++){
        int end_p = table_of_end_pos[i][prc];
        map_of_class_eraser(processes[prc] + n_nn ,i,end_p);
    }
}

/**
 * Removes all processes ending at a site. Called when a site becomes occupied.
 */
void metropolis::map_of_class_next_eraser(int i) {
    auto copy_tip=table_of_initial_pos[i];
    for(int i_site : copy_tip){
        for(int n=0; n<= nnn_atoms[i_site]; n++){
            for(size_t i_prc=0 ; i_prc < table_of_end_pos[i_site].size() ; i_prc++){
                if(table_of_end_pos[i_site][i_prc]==i){
                    map_of_class_eraser(table_of_processes[i_site][i_prc] + n, i_site , i);
                }
            }
        }
    }
}

/**
 * Inserts a new move into the process catalog.
 */
void metropolis::map_of_class_filler(int c, int i , int j){
    auto& processes = map_class_processes[c];
    for (auto& couple: processes){
        if(couple[0]<0){ couple[0]=i; couple[1]=j; movements_per_class[c]++; break; }
    }
}

/**
 * Returns the number of active processes in a specific class.
 */
int metropolis::get_MCP_size(int c){
    int counter=0;
    auto& processes = map_class_processes[c];
    for (auto& couple: processes){
        if(couple[0]<0) break;
        counter++;
    }
    return counter;
}

/**
 * Clears temporary processes (like Push or Kink moves) that depend on specific atom configurations.
 */
void metropolis::clear_dynamic_processes(){
    for(auto& [site, pv, cl] : dynamic_processes) map_of_class_eraser(cl, site, pv);
    dynamic_processes.clear();
}

/**
 * Core function of the KMC: updates the list of all possible moves for the "interested sites".
 * It classifies moves based on the current local coordination (neighbors) and facet type.
 */
void metropolis::classification(){
    clear_dynamic_processes();
    for(int site : interested_sites){
        if(atoms[site]){
            int nn_count = nnn_atoms[site];
            for(int t=0 ; t <= 6 ; t++) map_of_class_position_eraser(site, t);

            std::vector<int>& top = table_of_processes[site];
            for(size_t j=0 ; j < (int)top.size() ; j++){  
                int base_process = top[j]; 
                int pv = table_of_end_pos[site][j];
                
                if(!atoms[pv]){ 
                    int actual_class = base_process + nn_count;
                    if(actual_class >= 0 && actual_class < 100) map_of_class_filler(actual_class, site, pv);
                }
                // Handle Push and Kink logic...
                else if(atoms[pv] && s.get_TOPl(pv)==7 && s.get_TOPl(site)==8){
                    auto nn = s.get_table_of_nn(pv);
                    for(int n : nn) if(!atoms[n] && s.get_TOPl(n)==0){
                        int actual_class = 60 + nn_count;
                        dynamic_processes.emplace_back(site, n, actual_class);
                        map_of_class_filler(actual_class, site, n);
                    }
                }
                else if(atoms[pv] && base_process==72){
                    auto nn = s.get_table_of_nn(pv);
                    for(int n : nn) if(!atoms[n] && (s.get_TOPl(n)==0 || s.get_TOPl(n)==7)){
                        int actual_class = 72 + nn_count;
                        dynamic_processes.emplace_back(site, n, actual_class);
                        map_of_class_filler(actual_class, site, n);
                    }
                }
            }
        }
    }
}

/**
 * Fills the total probability vector P. 
 * P[i] = (number of moves in class i) * (rate of a single move in class i).
 */
void metropolis::probability_filler(){
    std::fill(P.begin(), P.end(), 0.0); 
    for(int i=0 ; i < n_class ; i++){
        if(movements_per_class[i] > 0) P[i] = movements_per_class[i] * probability_calculator(i);
    }
}

/**
 * Handles the deposition of a new atom from the flux.
 * Randomly selects a site on the permitted facets.
 */
void metropolis::deposition_func(){
    int p; 
    int n_nodep = deactivated_sites.size();
    int max_sites = crd.sites_size();
    while(true) {
        p = static_cast<int>(initial_pos(rng) * (max_sites - n_nodep));
        if(s.get_TOPl(p) < 2 && !atoms[p]) break; 
    }
    atoms[p] = true;
    pos = p;
    n_deposited++;
    output << "New atom: " << p << "    at time: "<< time <<std::endl;
}

/**
 * The KMC heart: 
 * 1. Calculates the time step (delta_t = -ln(r)/Total_Rate).
 * 2. Decides if the next event is a Deposition or a Diffusion jump.
 * 3. Picks the specific process and site based on their relative probabilities.
 */
void metropolis::time_prob_calc(){
    double r_time=dist_time(rng);
    double P_tot=0.0;
    double P_dep=crd.sites_size()*F;
    for(size_t i=0 ; i<P.size() ; i++) P_tot+=P[i];
    double P_diff=P_tot;
    P_tot+=P_dep;

    time += ( -(log(r_time)) / P_tot);

    double r_p = dist_time(rng) * P_tot;
    if(r_p>P_diff){
        deposition=true;
    } else {
        deposition=false;
        int chosen_class = 0;
        double p_i = 0;
        for(size_t i = 0 ; i < P.size() ; i++){
            if( (p_i < r_p) && (r_p < (p_i+P[i])) ) { chosen_class = i; break; }
            p_i += P[i];
        }

        int mcount = movements_per_class[chosen_class];
        std::uniform_int_distribution<> dist(0, mcount - 1);
        int movement = dist(rng);
        
        const auto& mv = map_class_processes[chosen_class][movement];
        next = mv[1];
        pos = mv[0];
        get_process_name(chosen_class);
        if(output_file) print_output();
        
        map_of_class_position_eraser(pos, nnn_atoms[pos]);
        map_of_class_next_eraser(next);
    }
}

void metropolis::print_output(){
    output<< process_name << "  from    " << pos <<"  to    "<< next <<std::endl;
}

/**
 * Saves the current configuration (occupied sites) to a file.
 */
void metropolis::print_configuration(){
    std::string str1 = std::to_string(n_deposited);
    int numZeros = std::max(0, 4 - static_cast<int>(str1.length()));
    str1.insert(0, numZeros, '0');
    std::string fileName = "occ" + str1 + ".out";
    std::ofstream intermed_config(fileName);
    c_crd.output_writer(intermed_config, c_crd.sites_size() + n_deposited, time);
    crd.output_writer_partial(intermed_config, n_deposited, atoms);
}

/**
 * Initializes simulation data: output files, geometry tables, and the first deposition.
 */
void metropolis::start_of_the_sim(){
    output.open("MC_processes.out");
    table_of_processes_filler(); 
    deactivated_sites=s.get_upper_sites();
    atoms.resize(crd.get_N(), false);
    deposition_func();
    interested_sites_calc(true);
    nnn_atoms.resize(atoms.size());
    nn_updater(deposition);
}

/**
 * The main KMC loop. Continues until simulation time reaches total_time.
 */
void metropolis::algorithm(){
    while(time<total_time){
        classification();
        probability_filler();
        time_prob_calc();
        if(deposition){
            deposition_func();
            print_configuration();
        } else{
            atoms[pos] = false;
            atoms[next]= true;
        }
        second_layer_activation();
        interested_sites_calc(deposition);
        nn_updater(deposition);
        pos=next;

        if(n_deposited >= (int) (filling * (float) crd.sites_size()) ) {
            F=0; 
            output<<"End of the deposition, there are "<<n_deposited<<" filled sites "<<std::endl;
        }
    }
}

/**
 * PUBLIC entry point to start the full simulation.
 */
void metropolis::simulation(){
    start_of_the_sim();
    algorithm();
    if(output_file) print_output();
    print_configuration();
}