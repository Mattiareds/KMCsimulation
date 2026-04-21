#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"
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
        ifile>>metodo;
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

    // Aggiorna ogni sito una sola volta
    for(int site : sites_to_update){
        int count = 0;
        auto nn_of_nn = s.get_table_of_nn(site);

        for(int nearest : nn_of_nn){
            // stessa logica del tuo codice:
            // conta solo atomi occupati, ma escludi pos e (se serve) next
            if(atoms[nearest]){
                if(nearest == pos && !pos_checked){ 
                    count++;
                    continue;
                }else if(nearest == pos && pos_checked){
                    continue;
                }
                if(!dep && nearest == next && !next_checked){
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


//OBSOLETO
//void metropolis::nn_updater(bool dep){
//    auto nn_pos = s.get_table_of_nn(pos);
//    std::cout<<nn_pos.size()<<std::endl;
//    for(int j=0 ; j< nn_pos.size() ; j++){ //nn of the jump
//        std::cout<<"nn of what happened " << nn_pos[j]<<std::endl;
//        int count=0;
//        auto nn_of_nn = s.get_table_of_nn(nn_pos[j]);
//        for(int nearest : nn_of_nn) { //for each nn update its number of nn
//            if(atoms[nearest]) count++; //count how many nn
//            if(nearest==pos) count--;
//        }
//        nnn_atoms[nn_pos[j]]=count; //update the number of nn
//        std::cout << "Updating site: " << nn_pos[j] << std::endl;
//    }
//    if(!dep){ //if diffusion
//        auto nn_next = s.get_table_of_nn(next);
//        for(int j=0 ; j< nn_next.size() ; j++){ //nn of the jump
//            int count=0;
//            auto nn_of_nn_1 = s.get_table_of_nn(nn_next[j]);
//            for(int nearest_1 : nn_of_nn_1) { //for each nn update its number of nn
//                if(atoms[nearest_1]) count++; //count how many nn
//                if(nearest_1 == next) count--;
//            }
//            nnn_atoms[nn_next[j]]=count; //update the number of nn
//             std::cout << "Updating again: " << nn_pos[j] << std::endl;
//        }
//    }
//}
//END OBSOLETO


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
                 table_of_processes[i].push_back(0); //111c 111c
                 table_of_end_pos[i].push_back(pv);
                 table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<0<<std::endl;
            }
            else if((s.get_TOPl(i)==0 && s.get_TOPl(pv)==0) || (s.get_TOPl(i)==7 && s.get_TOPl(pv)==7)) {
                table_of_processes[i].push_back(24); //100c 100c
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<24<<std::endl;
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)!=s.get_plane(pv,1)) {
                table_of_processes[i].push_back(10); //111 jump 111
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<10<<std::endl;
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==7) {
                table_of_processes[i].push_back(19); //111 jump 100
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<19<<std::endl;
            }
            else if(s.get_TOPl(i)==7 && s.get_TOPl(pv)==8) {
                table_of_processes[i].push_back(15); //100 jump 111
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<15<<std::endl;
            }
            else if(s.get_TOPl(i)==7 && s.get_TOPl(pv)==0) {
                table_of_processes[i].push_back(28); //100b 100c
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<28<<std::endl;
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==1) {
                table_of_processes[i].push_back(6); //111b 111c
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<6<<std::endl;
            }
            else if(s.get_TOPl(i)==0 && s.get_TOPl(pv)==7) {
                table_of_processes[i].push_back(31); //100c 100b
                table_of_end_pos[i].push_back(pv);
                table_of_initial_pos[pv].push_back(i);
                //std::cout<<"from "<<i<<"to "<<pv<<" type: "<<31<<std::endl;
            }
        }
    }
}

//ci vuole una funzione che aggiorni quella qui sopra quando viene introdotto anche il secondo strato
//sarà una funzione che gira solo su siti di bordo faccetta 100
//function for update the atoms position
void metropolis::second_layer_updates(){

}

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
    else if (cl >= 15 && cl <= 18)
    {
        bar[0] = 4;
        bar[1] = cl - 15;
        num=std::to_string(bar[1]);
        process_name = "100 to 111 with: " + num + " nearest neighbors";
    }
    else if (cl >= 19 && cl <= 23)
    {
        bar[0] = 3;
        bar[1] = cl - 19;
        num=std::to_string(bar[1]);
        process_name = "111 to 100 with: " + num + " nearest neighbors";
    }
    else if (cl >= 24 && cl <= 27)
    {
        bar[0] = 1;
        bar[1] = cl - 24;
        num=std::to_string(bar[1]);
        process_name = "100center to 100center same facet with: " + num + " nearest neighbors";
    }
    else if (cl >= 28 && cl <= 30)
    {
        bar[0] = 5;
        bar[1] = cl - 28;
        num=std::to_string(bar[1]);
        process_name = "100border to 100center same facet with: " + num + " nearest neighbors";
    }
    else if (cl >= 31 && cl <= 35)
    {
        bar[0] = 7;
        bar[1] = cl - 31;
        num=std::to_string(bar[1]);
        process_name = "100center to 100border with: " + num + " nearest neighbors";
    }
    else if (cl >= 36 && cl <= 40)
    {
        bar[0] = 8;
        bar[1] = cl - 36;
        num=std::to_string(bar[1]);
        process_name = "111border-push-100  with: " + num + " nearest neighbors";
    }
    else if (cl >= 41 && cl <= 43)
    {
        bar[0] = 9;
        bar[1] = cl - 41;
        num=std::to_string(bar[1]);
        process_name = "111 mixed to 100 with: " + num + " nearest neighbors";
    }
    else if (cl >= 44)
    {
        bar[0] = 10;
        bar[1] = cl - 44;
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
    std::cout<<"prob for process "<<b_id << " "<< rt_value<<std::endl;
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
            //std::cout <<"eraser     " <<couple[0] << "    "<<couple[1]<<std::endl;
            couple[0]=-1;
            couple[1]=-1;
            movements_per_class[c]--;
            shifter(c,ind,tg);
            break;
        }
        ind++;
    }
}

//when a process happens, the site is empty, so remove all possible processes that start in that place
void metropolis::map_of_class_position_eraser(int i) {
    //get the list of processes
    auto processes = table_of_processes[i];
    for(size_t prc=0 ; prc<processes.size() ; prc++){
        //get the final position of each process and erase
        int end_p = table_of_end_pos[i][prc];
        map_of_class_eraser(processes[prc],i,end_p);
    }
}

void metropolis::map_of_class_next_eraser(int i) {
    auto copy_tip=table_of_initial_pos[i];
    for(int i_site : copy_tip){
        for(int i_prc=0 ; i_prc < table_of_end_pos[i_site].size() ; i_prc++) if(table_of_end_pos[i_site][i_prc]==i){
            map_of_class_eraser(table_of_processes[i_site][i_prc], i_site , i);
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

/* update the table_of_processes for each site nn to what happened
update how many atoms per class
the only advantage of having a table of processes is that you don't have to re-calculate the possible move for each atom
*/
void metropolis::classification(){
    for(int site : interested_sites){
        if(atoms[site]){
            if(n_deposited > 1) map_of_class_position_eraser(site);
            
            // Get the list of BASE processes for this site
            std::vector<int>& top = table_of_processes[site];

            for(int j=0 ; j < (int)top.size() ; j++){  
                int base_process = top[j]; 
                int pv = table_of_end_pos[site][j];

                if(!atoms[pv]){ 
                    int nn_count = nnn_atoms[site];
                    int actual_class = base_process + nn_count; // Calculate temp index
                    std::cout<<"The number of nn of process: "<<site<< " is " <<nn_count<<std::endl;


                    if(actual_class >= 0 && actual_class < 100) {
                        map_of_class_filler(actual_class, site, pv);                        
                    } else {
                        std::cerr << "Error: Class " << actual_class << " exceeds array size!" << std::endl;
                    }

                    //aggiungere eraser dei processi con numero nullo di pv
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
            std::cout<<"class "<<i<<" with movements: "<<movements_per_class[i]<<std::endl;
            P[i] = movements_per_class[i] * p;
        }
    }
}

//commento: avrei anche potuto eliminare i siti di edge tra quelli papabili per la diffusione (e tra i siti in generale, ma se uno vuole fare altri strati poi quelli potrebbero servire)
void metropolis::deposition_func(){
    bool no_edge=false;
    int p; 
    while(no_edge==false){
        p = static_cast<int>(initial_pos(rng) * crd.sites_size());
        if(p < 0) p = 0;
        if(p >= crd.sites_size()) p = crd.sites_size() - 1;
        if(s.get_TOPl(p)<2) no_edge=true;
    }
    atoms[p] = true;
    pos=p;
    n_deposited++;
    output << "New atom: " << p << "    at time: "<< time <<std::endl;
    //std::cout << "New atom: " << p << "    at time: "<< time <<std::endl;
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
            return;
        }
        std::uniform_int_distribution<> dist(0, mcount - 1);
        int movement = dist(rng);
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
        map_of_class_position_eraser(pos);
        map_of_class_next_eraser(next);
    }
}


void metropolis::print_output(){
    //pensare a come fare il file di output
    output<< process_name << "  from    " << pos <<"  to    "<< next <<std::endl;
    //std::cout<<process_name << "  from    " << pos <<"  to    "<< next <<std::endl;
}

//CHECK ALL
void metropolis::start_of_the_sim(){
    output.open("MC_processes.txt");
    table_of_processes_filler(); //initialize
    atoms.resize(crd.get_N(), false);
    deposition_func();
    interested_sites_calc(true);
    nnn_atoms.resize(atoms.size());
    nn_updater(deposition);
    //std::cout<<"start with this config:"<<std::endl;
    for(int i=0; i<atoms.size() ; i++) std::cout<<i<<"  "<<atoms[i]<<std::endl;

}

void metropolis::algorithm(){
    //compute all the possible processes
    while(time<total_time){
        //std::cout<<"bug3" << std::endl;
        classification();
        probability_filler();
        time_prob_calc();
        //from now next is imporant 
        //std::cout<<"bug4" << std::endl;
        if(deposition){
            deposition_func();
            //crd.output_writer(output,crd.sites_size());
        } else{
            atoms[pos]= false;
            atoms[next]= true;
        }
        interested_sites_calc(deposition);
        nn_updater(deposition);
        pos=next;
        //here next is useful
        if(n_deposited>1) for(int i=0; i<atoms.size() ; i++) std::cout<<"atoms: "<<i<<"  "<<atoms[i]<<std::endl;
        //std::cout<<"-----------------------------------end mc loop------------------------------------"<<std::endl;
    }
}

//the PUBLIC function for the simulation, the only one that has to be called by the user
//in addiction with something for outputs
void metropolis::simulation(){
    //std::cout<<"bug0000" << std::endl;
    start_of_the_sim();
    algorithm();
    //std::cout<<"bug4" << std::endl;
    print_output();
}
