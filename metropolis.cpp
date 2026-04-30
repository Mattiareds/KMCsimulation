#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"
#include <cstdlib>
#include <iostream>
#include <unordered_set>


/*per il calcolo delle probabilità:
1) si ha un vettore con tutte le probabilità associate ad ogni sito, a prescindere dall'occupazione o meno da parte di un atomo (probabilità di base, no pv)
2) chiaramente se il sito non è occupato non si tiene conto di quelle probabilità
3) si ha un altro vettore che ha al suo interno il numero di pv per ogni sito
4) quando si effettua una mossa si vanno a ricalcolare il numero dei pv solo dei siti pv coinvolti nella mossa (quindi tutti i pv del sito di partenza e del sito di arrivo, a parte i siti
coinvolti che non avranno nessun tipo di modifica )
5) quando si calcola la probailità finale, si tiene conto del numero dei pv che sono in quel sito
6) la probabilità va divisa per tipo di processo coinvolto!!! quindi una per ogni barriera!!!, quindi quando si effettua una mossa bisogna ricalcolare la probabilità complessiva di tutti
i processi di quel tipo
7) uno sa che processo avviene, questo andr� a variare la probabilità di alcuni processi generali. se io conosco, per ogni atomo, che tipo di processi mi fa, mi vado ad aggiungere/togliere
probabilità da quelle tipologie di processo, andando a ridurre drasticamente il tempo di calcolo di quel processi. La funzione per calcolare i nn la ho già, mi serve un vettore che mi dica
per ogni tipo di processo quali sono i siti, così se un sito cambia vado a calcolare solo la total probability di quel sito (moltiplicata per il riempimento)
8)assumendo che la probabilità sia calcolata, ad una mossa che modifica il numero di pv di un atomo vado a ridurre/aumentare la probabilità totale di quel tipo di processi in base alla differenza
tra i primi vicini di prima e di ora

*/

void metropolis::file_reader(std::ifstream& ifile){
    if(ifile){
        ifile>>mc_step;
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

        for (int i=0; i<100; i++) movements_per_class[i]=0;
        for (int i=0 ; i< 100; i++){ for(int j=0 ; j<1000 ; j++){for (int k=0; k<2; k++) map_class_processes[i][j][k]=-1;}}
        P.resize(n_class, 0.0);
        n_deposited=0;
    }
}

void metropolis::interested_sites_calc(bool dep){
    auto nn_pos = s.get_table_of_nn(pos);
    interested_sites.clear();
    interested_sites.push_back(pos);
    for(int i=0 ; i<nn_pos.size() ; i++){
        if(nn_pos[i]!=pos) interested_sites.push_back(nn_pos[i]);
    }
    if(!dep){
        auto nn_next = s.get_table_of_nn(next);
        for(int i=0 ; i<nn_next.size() ; i++){
            bool already_here=false;
            for(size_t l=0 ; l< interested_sites.size() ; l++) if(nn_next[i]==interested_sites[l]) { already_here=true; break; }
            if(nn_next[i]!=pos && nn_next[i]!=next && !already_here) interested_sites.push_back(nn_next[i]);
        }
    }
}

//function for update the number of nn after a jump
//pos is the initial and, after, the final position of the jump
void metropolis::nn_updater(bool dep){
    bool pos_checked=false;
    bool next_checked=false;
    std::unordered_set<int> sites_to_update;

    // Vicini di pos
    auto nn_pos = s.get_table_of_nn(pos);
    sites_to_update.insert(nn_pos.begin(), nn_pos.end());

    // Se diffusione, aggiungi anche i vicini di next
    if(!dep){
        auto nn_next = s.get_table_of_nn(next);
        sites_to_update.insert(nn_next.begin(), nn_next.end());
    }
    //debug_file<<"sites to update "<<sites_to_update.size()<<std::endl;
    // Aggiorna ogni sito una sola volta
    for(int site : sites_to_update){
        
        int count = 0;
        auto nn_of_nn = s.get_table_of_nn(site);
        auto edge = s.get_edge_map(site);

        for(int i=0 ; i< nn_of_nn.size() ; i++){
            int nearest=nn_of_nn[i];
            int edges= edge[i];
            //if( (s.get_plane(site, 1)==s.get_plane(nearest,1)) || ((s.get_plane(site, 1)==14||s.get_plane(site, 1)==15) || (s.get_plane(nearest, 1)==14 || s.get_plane(nearest, 1)==15) ) || (s.get_TOPl(site)==10 || s.get_TOPl(nearest)==10 ) ){
                if(atoms[nearest] && (edges < 0)){ //if lt 0 not changed: same facet
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
            //}
        }
        nnn_atoms[site] = count;
    }
}

bool metropolis::is_border(int i){
    if(s.get_TOPl(i)==8 || s.get_TOPl(i)==7 || s.get_TOPl(i)==15){
        return true;
    }
    return false;
}

bool metropolis::are_different_border(int i, int j){
    if(is_border(i) && is_border(j)){
        if(s.get_TOPl(i)!=s.get_TOPl(j)){
            return true;
        }
    }
    return false;
}


//fill the processes based on the geometry (initialize)
void metropolis::table_of_processes_filler(){
    table_of_processes.resize(crd.sites_size(), std::vector<int>(0,0));
    table_of_end_pos.resize(crd.sites_size(), std::vector<int>(0,0));
    table_of_initial_pos.resize(crd.sites_size(), std::vector<int>(0,0));
    //run over all the sites
    for(int i=0 ; i<crd.sites_size() ; i++){
        std::vector<int> nns= s.get_table_of_nn(i);
        //run over all the pv
        for(int j=0 ; j<nns.size() ; j++){
            int pv=nns[j];
            if(s.get_TOPl(i)==1 || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)==s.get_plane(pv,1))){
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 0){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(0); //111c 111c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<0<<std::endl;
                }
            }
            else if((s.get_TOPl(i)==0 && s.get_TOPl(pv)==0) || (s.get_TOPl(i)==7 && s.get_TOPl(pv)==7)) {
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 25){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(25); //100c 100c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<25<<std::endl;
                }
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)!=s.get_plane(pv,1)) {
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 10){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(10); //111 jump 111
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<10<<std::endl;
                }
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==7) {
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 20){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(20); //111 jump 100
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<20<<std::endl;
                }
            }
            else if(s.get_TOPl(i)==7 && s.get_TOPl(pv)==8) {
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 15){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(15); //100 jump 111
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<15<<std::endl;
                }
            }
            else if((s.get_TOPl(i)==7 && s.get_TOPl(pv)==0)) {
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 33){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(33); //100b 100c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<33<<std::endl;
                }
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==1) {
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 6){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(6); //111b 111c
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<6<<std::endl;
                }
            }
            else if((s.get_TOPl(i)==0 && s.get_TOPl(pv)==7)) {
                 bool already_present = false;
                for(int k = 0; k < table_of_end_pos[i].size(); k++){
                    if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 40){
                        already_present = true; break;
                    }
                }
                if(!already_present){
                    table_of_processes[i].push_back(40); //100c 100b
                    table_of_end_pos[i].push_back(pv);
                    table_of_initial_pos[pv].push_back(i);
                    //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<40<<std::endl;
                }
            }
        }
    }
}


//fill the processes based on the geometry (update)
void metropolis::table_of_processes_updater(int i){    
    std::vector<int> nns= s.get_table_of_nn(i);
    //run over all the pv
    for(int j=0 ; j<nns.size() ; j++){
        int pv=nns[j];
        if(s.get_TOPl(i)==1 || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)==s.get_plane(pv,1)) || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==10) || (s.get_TOPl(i)==10 && s.get_TOPl(pv)==10) || (s.get_TOPl(i)==10 && s.get_TOPl(pv)==8) ){
            bool already_present = false;
            for(int k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 0){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(0); //111c 111c
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<0<<std::endl;
            }
        }
        else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)!=s.get_plane(pv,1)) {
            bool already_present = false;
            for(int k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 10){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(10); //111 jump 111
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<10<<std::endl;
            }
        }
        else if((s.get_TOPl(i)==8 || s.get_TOPl(i)==10) && s.get_TOPl(pv)==7) {
            bool already_present = false;
            for(int k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 20){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(20); //111 jump 100
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<20<<std::endl;
            }
        }
        else if(s.get_TOPl(i)==7 && (s.get_TOPl(pv)==8 || s.get_TOPl(pv)==10)) {
            bool already_present = false;
            for(int k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 15){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(15); //100 jump 111
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<15<<std::endl;
            }
        }
        else if((s.get_TOPl(i)==8|| s.get_TOPl(i)==8) && s.get_TOPl(pv)==1) {
            bool already_present = false;
            for(int k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 6){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(6); //111b 111c
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<6<<std::endl;
            }
        }
        else if(s.get_TOPl(i)==10 && s.get_TOPl(pv)>10){ //jump mixed site->new site
            bool already_present = false;
            for(int k = 0; k < table_of_end_pos[i].size(); k++){
                if(table_of_end_pos[i][k] == pv && table_of_processes[i][k] == 52){
                    already_present = true; break;
                }
            }
            if(!already_present){
                table_of_processes[i].push_back(52);
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //debug_file<<"from "<<i<<"to "<<pv<<" type: "<<52<<" mixed-new "<<std::endl;
            }
        }
        //here we don't have 100 processes since are not useful 
    }
    
}

bool metropolis::is_deactivated(int i){
    for(int j : deactivated_sites) if (j==i) return true;
    return false;
}


//the old edge site is now a site that allows to pass to the new facet
std::array<int,2> metropolis::activate_edges(int upper_site, int pv){
    std::array<int,2> back={-1,-1};
    auto edge = s.get_edge_map(pv);
    //debug_file<<"size "<<edge.size()<<std::endl;
    int index=0;
    for(int i=0 ; i<edge.size() ; i++){
        //debug_file<<edge[i]<<std::endl;
        if(edge[i]>=0){
            //update the edge situation: new nn which is the high site
            s.set_table_of_nn(pv , i, edge[i]);
            s.append_nn(edge[i], upper_site);
            s.set_TOPl(edge[i],10);  //10=mixed
            //debug_file<<" update: new nn for ex-edge guy: "<< edge[i] <<" "<<upper_site <<std::endl;
            for(auto old_nn : s.get_table_of_nn(edge[i])){
                if(s.get_TOPl(old_nn)==8){
                    //update old nn informations
                    int oppc = 0;
                    for (int p : table_of_processes[old_nn]){
                        if(p==0) oppc++;
                    }
                    s.append_nn(old_nn, edge[i]);
                    if(oppc > 2) s.set_TOPl(old_nn, 1);
                    s.set_plane(edge[i], 20);
                    //debug_file<<" update: new nn for ex_border 111 guy: "<< old_nn <<" "<<edge[i]<<std::endl;
                    table_of_processes_updater(old_nn); //old_nn should be an old 111b which becomes a 111c
                }
            }
            table_of_processes_updater(edge[i]);
            //reset the edge map 
            s.reset_edge_map(pv, i);

            //return now since it is only one
            back[index]=edge[i];
            index++;
        }
        if(index>1) return back;
    }
    return back;
}



//ci vuole una funzione che aggiorni quella qui sopra quando viene introdotto anche il secondo strato
//sarà una funzione che gira solo su siti di bordo faccetta 100
//function for update the atoms position
void metropolis::second_layer_updates(int upper_site){
//bisogna mettere qua quello che sta sopra: controllo ad ogni ciclo (girando solo sugli ultimi) e aggiorno la tabella 
    auto nn=s.get_table_of_nn(upper_site);
    for (int pv : nn){
        //debug_file<<pv<<std::endl;
        if(s.get_TOPl(upper_site)==s.get_TOPl(pv)){
            if(is_deactivated(pv)==false){
                table_of_processes[upper_site].push_back(25); //100c 100c
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
                //debug_file<<"from "<<upper_site<<"to "<<pv<<" type: "<<25<<" 100c to 100c "<<std::endl;
            }
        }
        else if((s.get_TOPl(upper_site)==14 && s.get_TOPl(pv)==15) ){
            if(is_deactivated(pv)==false){
                table_of_processes[upper_site].push_back(40); //100c 100b
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
                //debug_file<<"from "<<upper_site<<"to "<<pv<<" type: "<<40<<" 100c to 100b "<<std::endl;
            }
        }
        else if((s.get_TOPl(upper_site)==15 && s.get_TOPl(pv)==14)){
            if(is_deactivated(pv)==false){
                table_of_processes[upper_site].push_back(33); //100b 100c
                table_of_end_pos[upper_site].push_back(pv);
                table_of_initial_pos[pv].push_back(upper_site);
                //debug_file<<"from "<<upper_site<<"to "<<pv<<" type: "<<33<<" 100b to 100c "<<std::endl;
            }
        }
        else if(s.get_TOPl(pv)==0 || s.get_TOPl(pv)==7){ // fake kink (fix after)
            table_of_processes[upper_site].push_back(56);
            table_of_end_pos[upper_site].push_back(pv);
            table_of_initial_pos[pv].push_back(upper_site);
            //debug_file<<"from "<<upper_site<<"to "<<pv<<" type: "<<56<<" kink "<<std::endl;
        }
        if(s.get_TOPl(pv)==7){ //100 new to 111 edge site
            auto edge_sites = activate_edges(upper_site,pv);
            for (auto edge_site : edge_sites){
                if(edge_site>=0){
                    table_of_processes[upper_site].push_back(15);
                    table_of_end_pos[upper_site].push_back(edge_site);
                    table_of_initial_pos[edge_site].push_back(upper_site);
                    //debug_file<<"from "<<upper_site<<"to "<<edge_site<<" type: "<<15<<" new 111 edge "<<std::endl;
                }
            }
        }
    }
}


void metropolis::second_layer_activation(){
    //debug_file<<"second layer test "<<std::endl;
    for (int i=0 ; i<deactivated_sites.size() ; i++){
        //debug_file<<" second layer loop "<<std::endl;
        int counter = 0;
        int upper_site = deactivated_sites[i];
        for(int j : s.get_table_of_nn(upper_site)){
            if(atoms[j]==true) counter++;
        }
        if(counter >= 4){
            //debug_file<<upper_site<<" updating the second layer           --------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
            second_layer_updates(upper_site); 
            deactivated_sites.erase(deactivated_sites.begin()+i);
        }
    }
}

//second layer deactivation is not implemented since it's rare... extremely
void metropolis::barrier_chooser(int cl, int &bar_id , int &neighbours ){
    std::string num;
    int bar[2];

    if (cl >= 6 && cl <= 9)
    {
        bar[0] = 6;
        bar[1] = cl - 6;
        num=std::to_string(bar[1]);
        process_name = "111border to 111center same facet with: " + num + " nearest neighbors";
    }
    else if (cl >= 10 && cl <= 14)
    {
        bar[0] = 2;
        bar[1] = cl - 10;
        num=std::to_string(bar[1]);
        process_name = "111border to 111border different facets with: " + num + " nearest neighbors";
    }
    else if (cl >= 15 && cl <= 19)
    {
        bar[0] = 4;
        bar[1] = cl - 15;
        num=std::to_string(bar[1]);
        process_name = "100 to 111 with: " + num + " nearest neighbors";
    }
    else if (cl >= 20 && cl <= 24)
    {
        bar[0] = 3;
        bar[1] = cl - 20;
        num=std::to_string(bar[1]);
        process_name = "111 to 100 with: " + num + " nearest neighbors";
    }
    else if (cl >= 25 && cl <= 32)
    {
        bar[0] = 1;
        bar[1] = cl - 25;
        num=std::to_string(bar[1]);
        process_name = "100center to 100center same facet with: " + num + " nearest neighbors";
    }
    else if (cl >= 33 && cl <= 39)
    {
        bar[0] = 5;
        bar[1] = cl - 33;
        num=std::to_string(bar[1]);
        process_name = "100border to 100center same facet with: " + num + " nearest neighbors";
    }
    else if (cl >= 40 && cl <= 46)
    {
        bar[0] = 7;
        bar[1] = cl - 40;
        num=std::to_string(bar[1]);
        process_name = "100center to 100border with: " + num + " nearest neighbors";
    }
    else if (cl >= 47 && cl <= 51)
    {
        bar[0] = 8;
        bar[1] = cl - 47;
        num=std::to_string(bar[1]);
        process_name = "111border-push-100  with: " + num + " nearest neighbors";
    }
    else if (cl >= 52 && cl <= 55)
    {
        bar[0] = 9;
        bar[1] = cl - 52;
        num=std::to_string(bar[1]);
        process_name = "111 mixed to 100 double with: " + num + " nearest neighbors";
    }
    else if (cl >= 56)
    {
        bar[0] = 10;
        bar[1] = cl - 56;
        num=std::to_string(bar[1]);
        process_name = "kink: 100shell to 100core with: " + num + " nearest neighbors";
    }
    else
    {
        bar[0] = 0;
        bar[1] = cl;
        num=std::to_string(bar[1]);
        process_name = "111center to 111center same facet with: " + num + " nearest neighbors";
    }
    bar_id=bar[0];
    neighbours=bar[1];
}


//calculate
double metropolis::probability_calculator(int cl){
    int b_id = 0;
    int neighbours = 0;
    barrier_chooser(cl, b_id, neighbours);
    double rt_value = nu_0 * (exp((-1.) *(barriere[b_id] + (neighbours*E_b)) / (Kb*temperatura) )) ;
    //debug_file<<"prob for process "<<b_id << " "<< rt_value<<std::endl;
    return (rt_value);
}

void metropolis::shifter(int c, int ind, int tg){
    for(int i = ind ; i<tg ; i++){
        map_class_processes[c][i][0] = map_class_processes[c][i+1][0];
        map_class_processes[c][i][1] = map_class_processes[c][i+1][1];
    }
}

//eraser for the single process
void metropolis::map_of_class_eraser(int c, int i, int j) {
    auto& processes= map_class_processes[c];
    int tg = get_MCP_size(c);
    for(int ind=0; ind< tg ; ind++){
        auto& couple = processes[ind];
        if(couple[0]==-1) break;
        if(couple[0]==i && couple[1]==j){
            //debug_file <<"eraser     " <<couple[0] << "    "<<couple[1]<<"      "<<c<<std::endl;
            couple[0]=-1;
            couple[1]=-1;
            movements_per_class[c]--;
            shifter(c,ind,tg);
            break;
        }
        //ind++;
    }
}

//when a process happens, the site is empty, so remove all possible processes that start in that place
void metropolis::map_of_class_position_eraser(int i, int n_nn) {
    //get the list of processes
    auto processes = table_of_processes[i];
    for(size_t prc=0 ; prc<processes.size() ; prc++){
        //get the final position of each process and erase
        int end_p = table_of_end_pos[i][prc];
        map_of_class_eraser(processes[prc] + n_nn ,i,end_p);
    }
}

void metropolis::map_of_class_next_eraser(int i) {
    auto copy_tip=table_of_initial_pos[i];
    for(int i_site : copy_tip){
        for(int n=0; n<= nnn_atoms[i_site]; n++){
            for(int i_prc=0 ; i_prc < table_of_end_pos[i_site].size() ; i_prc++) if(table_of_end_pos[i_site][i_prc]==i){
                map_of_class_eraser(table_of_processes[i_site][i_prc] + n, i_site , i);
            }
        }
    }
}

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

void metropolis::clear_dynamic_processes(){
    for(auto& [site, pv, cl] : dynamic_processes){
        map_of_class_eraser(cl, site, pv);
    }
    dynamic_processes.clear();
}

/* update the table_of_processes for each site nn to what happened
update how many atoms per class
the only advantage of having a table of processes is that you don't have to re-calculate the possible move for each atom
*/
void metropolis::classification(){
    clear_dynamic_processes();
    for(int site : interested_sites){
        if(atoms[site]){
            int nn_count = nnn_atoms[site];
            
            for(int t=0 ; t <= 6 ; t++) map_of_class_position_eraser(site, t);

            //debug_file<<"eraser called for site                                            "<< site<<std::endl;
            
            // Get the list of BASE processes for this site
            std::vector<int>& top = table_of_processes[site];

            for(int j=0 ; j < (int)top.size() ; j++){  
                int base_process = top[j]; 
                int pv = table_of_end_pos[site][j];
                
                if(!atoms[pv]){ 
                    int actual_class = base_process + nn_count; // Calculate temp index
                    //debug_file<<"The number of nn of process  "<< base_process << " from site: "<<site<< " to site:  "<<pv<< " is " <<nn_count<<std::endl;

                    if(actual_class >= 0 && actual_class < 100) {
                        map_of_class_filler(actual_class, site, pv);                        
                    } else {
                        std::cerr << "Error: Class " << actual_class << " exceeds array size!" << std::endl;
                    }
                }
                // push
                else if(atoms[pv] && s.get_TOPl(pv)==7 && s.get_TOPl(site)==8){
                    auto nn = s.get_table_of_nn(pv);
                    for(int n : nn){
                        if(!atoms[n] && s.get_TOPl(n)==0){
                            int actual_class = 47 + nn_count;
                            dynamic_processes.emplace_back(site, n, actual_class);
                            map_of_class_filler(actual_class, site, n);
                            //debug_file<<"Push possible from: "<<site<<" to "<<n<<std::endl;
                        }
                    }
                }
                // kink
                else if(atoms[pv] && base_process==56){
                    auto nn = s.get_table_of_nn(pv);
                    for(int n : nn){
                        if(!atoms[n] && (s.get_TOPl(n)==0 || s.get_TOPl(n)==7)){
                            int actual_class = 56 + nn_count;
                            dynamic_processes.emplace_back(site, n, actual_class);
                            map_of_class_filler(actual_class, site, n);
                            //debug_file<<"Kink possible from: "<<site<<" to "<<n<<std::endl;
                        }
                    }
                }
            }
        }
    }
}

void metropolis::probability_filler(){
    // Reset the vector so we don't add to old values
    std::fill(P.begin(), P.end(), 0.0); 
    for(int i=0 ; i < n_class ; i++){
        if(movements_per_class[i] > 0) {
            double p = probability_calculator(i);
            //debug_file<<"class "<<i<<" with movements: "<<movements_per_class[i]<<std::endl;
            if(movements_per_class[i]> 6*n_deposited) std::exit(0);
            P[i] = movements_per_class[i] * p;
        }
    }
}

//commento: avrei anche potuto eliminare i siti di edge tra quelli papabili per la diffusione (e tra i siti in generale, ma se uno vuole fare altri strati poi quelli potrebbero servire)
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
    //debug_file<< "New atom: " << p << "    at time: "<< time <<std::endl;
}


void metropolis::time_prob_calc(){
    //update the time
    double r_time=dist_time(rng);
    double P_tot=0.0;
    double P_dep=crd.sites_size()*F;
    for(int i=0 ; i<P.size() ; i++) P_tot+=P[i];
    double P_diff=P_tot;
    P_tot+=P_dep;
    //std::cout<<"prob dep "<<P_dep<<" prob diffusion "<<P_diff <<" prob total "<<P_tot<<std::endl;
    time += ( -(log(r_time)) / P_tot);

    //choice of the class
    double random = dist_time(rng);
    double r_p = random * P_tot;
    //std::cout<<"random number"<<r_p<<std::endl;
    if(r_p>P_diff){
        deposition=true;
    } else {
        //std::cout<<"diffusion"<<std::endl;
        deposition=false;
        int chosen_class;
        double p_i = 0;
        for(int i = 0 ; i < P.size() ; i++){
            if( (p_i < r_p) && (r_p < (p_i+P[i])) ) chosen_class = i;
            p_i += P[i];
        }

        //choice of the process
        int mcount = movements_per_class[chosen_class];
        if(mcount <= 0){
            std::cerr << "Warning: chosen_class has no movements: " << chosen_class << std::endl;
            std::exit(0);
        }
        std::uniform_int_distribution<> dist(0, mcount - 1);
        int movement = dist(rng);
        if(atoms[map_class_processes[chosen_class][movement][0]]==false){
            std::cerr<<"Atomo inesistente... "<<std::endl;
            std::exit(0);
        }
        //std::cout << "Selected class " << chosen_class << " That has number of elements " << get_MCP_size(chosen_class) << "      "<<movements_per_class[chosen_class]<<" Selected movement " << movement << std::endl;
        if(movement < 0 || movement >= get_MCP_size(chosen_class)){
            std::cerr << "Error: movement index out of range: " << movement << " >= " << get_MCP_size(chosen_class) << std::endl;
            return;
        }
        const auto& mv = map_class_processes[chosen_class][movement];
        next = mv[1];
        pos = mv[0];
        print_output();
        //std::cout<<"ora pos "<<pos<<" mentre next "<< next<<std::endl;
        //remove the already executed process to the list of possible processes 
        map_of_class_position_eraser(pos, nnn_atoms[pos]);
        map_of_class_next_eraser(next);
    }
}


void metropolis::print_output(){
    output<< process_name << "  from    " << pos <<"  to    "<< next <<std::endl;
    //debug_file<<process_name << " PROCESS THAT HAPPENED: from    " << pos <<"  to    "<< next <<std::endl;
}

void metropolis::print_final_output(){
    output<<"Final configuration obtained after "<<time<<" seconds; "<<mc_step<<" Monte Carlo steps"<<std::endl;
    crd.output_writer(output,crd.sites_size());
}

//CHECK ALL
void metropolis::start_of_the_sim(){
    output.open("MC_processes.txt");
    //debug_file.open("debug.out");
    table_of_processes_filler(); //initialize
    deactivated_sites=s.get_upper_sites();
    atoms.resize(crd.get_N(), false);
    deposition_func();
    interested_sites_calc(true);
    nnn_atoms.resize(atoms.size());
    nn_updater(deposition);
    //debug
    //debug_file<<"start with this config:"<<std::endl;
    //for(int i=0; i<atoms.size() ; i++) debug_file<<i<<"  "<<atoms[i]<<" in plane "<<s.get_TOPl(i)<<std::endl;
    //end debug
}

void metropolis::algorithm(){
    //compute all the possible processes
    while(time<total_time){
        classification();
        probability_filler();
        time_prob_calc();
        if(deposition){
            deposition_func();
            crd.output_writer(output,crd.sites_size());
        } else{
            atoms[pos] = false;
            atoms[next]= true;
        }
        second_layer_activation();
        interested_sites_calc(deposition);
        nn_updater(deposition);
        pos=next;
        //stopper
        if(n_deposited >= (int) (filling * (float) crd.sites_size()) ) {
            F=0; 
            output<<"End of the deposition, there are "<<n_deposited<<" filled sites "<<std::endl;
        }

        /*  FOR DEBUGGING
        int cnt=0;
         for(int i=0; i<atoms.size() ; i++){
            debug_file<<"atoms: "<<i<<"  "<<atoms[i]<<std::endl;
            if(atoms[i]) cnt++;
        }
        if(cnt!=n_deposited) {std::cout<<"what "<<cnt<<" "<<n_deposited<<std::endl; std::exit(0);}
        for(int c=0 ; c<100 ; c++){
            for(int i=0; i<1000; i++){
                if(map_class_processes[c][i][0]>0){
                    debug_file<<"process "<<c<<" from "<<map_class_processes[c][i][0]<<" to "<<map_class_processes[c][i][1]<<std::endl;
                }
            }
        }
        debug_file<<n_deposited<<" -----------------------------------end mc loop------------------------------------"<<std::endl;
        */
    }
}

//the PUBLIC function for the simulation, the only one that has to be called by the user
//in addiction with something for outputs
void metropolis::simulation(){
    start_of_the_sim();
    algorithm();
    print_output();
}
