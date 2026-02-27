#include "coordinates.h"
#include "geometry.h"
#include "metropolis.h"

/*per il calcolo delle probabilità:
1) si ha un vettore con tutte le probabilità associate ad ogni sito, a prescindere dall'occupazione o meno da parte di un atomo (probabilità di base, no pv)
2) chiaramente se il sito non è occupato non si tiene conto di quelle probabilità 
3) si ha un altro vettore che ha al suo interno il numero di pv per ogni sito 
4) quando si effettua una mossa si vanno a ricalcolare il numero dei pv solo dei siti pv coinvolti nella mossa (quindi tutti i pv del sito di partenza e del sito di arrivo, a parte i siti 
coinvolti che non avranno nessun tipo di modifica )
5) quando si calcola la probailità finale, si tiene conto del numero dei pv che sono in quel sito 
6) la probabilità va divisa per tipo di processo coinvolto!!! quindi una per ogni barriera!!!, quindi quando si effettua una mossa bisogna ricalcolare la probabilità complessiva di tutti
i processi di quel tipo
7) uno sa che processo avviene, questo andrà a variare la probabilità di alcuni processi generali. se io conosco, per ogni atomo, che tipo di processi mi fa, mi vado ad aggiungere/togliere 
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

        atoms_per_class.resize(n_class, 0);
        map_class_processes.resize(n_class,std::vector<std::vector<int>>());
        P.resize(n_class, 0.0);
    }
}

void metropolis::interested_sites_calc(bool dep){
    std::vector<int> nn_pos = s.get_table_of_nn(pos);
    interested_sites.push_back(pos);
    for(int i=0 ; i<nn_pos.size() ; i++){
        if(nn_pos[i]!=pos || nn_pos[i]!=next)interested_sites.push_back(nn_pos[i]);
    }
    if(!dep){
        std::vector<int> nn_next = s.get_table_of_nn(next);
        interested_sites.push_back(next);
        for(int i=0 ; i<nn_next.size() ; i++){
            if(nn_next[i]!=pos || nn_next[i]!=next)interested_sites.push_back(nn_next[i]);
        }
    }
}

//function for update the number of nn after a jump
//pos is the initial and, after, the final position of the jump 
void metropolis::nn_updater(bool dep){
    std::vector<int> nn_pos = s.get_table_of_nn(pos);
    for(int j=0 ; j< nn_pos.size() ; j++){ //nn of the jump
        int count=0;
        std::vector<int> nn_of_nn = s.get_table_of_nn(nn_pos[j]);
        for(int i=0 ; i<nn_of_nn.size() ; i++) { //for each nn update its number of nn
            if(atoms[nn_of_nn[i]]) count++; //count how many nn 
        }
        nnn_atoms[nn_pos[j]]=count-1; //update the number of nn 
    }
    if(!dep){
        std::vector<int> nn_next = s.get_table_of_nn(next);
        for(int j=0 ; j< nn_next.size() ; j++){ //nn of the jump
            int count=0;
            std::vector<int> nn_of_nn_1 = s.get_table_of_nn(nn_next[j]);
            for(int i=0 ; i<nn_of_nn_1.size() ; i++) { //for each nn update its number of nn
                if(atoms[nn_of_nn_1[i]]) count++; //count how many nn 
            }
            nnn_atoms[nn_next[j]]=count-1; //update the number of nn 
        }
    }
}

//fill the processes based on the geometry (initialize)
void metropolis::table_of_processes_filler(){
    table_of_processes.resize(crd.sites_size(), std::vector<int>(0,0));
    table_of_end_pos.resize(crd.sites_size(), std::vector<int>(0,0));
    //giro su tutti siti
    for(int i=0 ; i<crd.sites_size() ; i++){
        std::vector<int> nns= s.get_table_of_nn(i);
        for(int j=0 ; j<nns.size() ; j++){
            int pv=nns[j];
            if(s.get_TOPl(i)==1 || (s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)==s.get_plane(pv,1))){
                 table_of_processes[i].push_back(0); //111c 111c
                 table_of_end_pos[i].push_back(pv);
            }
            else if((s.get_TOPl(i)==0 && s.get_TOPl(pv)==0) || (s.get_TOPl(i)==7 && s.get_TOPl(pv)==7)) {
                table_of_processes[i].push_back(24); //100c 100c
                table_of_end_pos[i].push_back(pv);
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==8 && s.get_plane(i,1)!=s.get_plane(pv,1)) {
                table_of_processes[i].push_back(10); //111 jump 111
                table_of_end_pos[i].push_back(pv);
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==7) {
                table_of_processes[i].push_back(19); //111 jump 100
                table_of_end_pos[i].push_back(pv);
            }
            else if(s.get_TOPl(i)==7 && s.get_TOPl(pv)==8) {
                table_of_processes[i].push_back(15); //100 jump 111
                table_of_end_pos[i].push_back(pv);
            }
            else if(s.get_TOPl(i)==7 && s.get_TOPl(pv)==0) {
                table_of_processes[i].push_back(28); //100b 100c
                table_of_end_pos[i].push_back(pv);
            }
            else if(s.get_TOPl(i)==8 && s.get_TOPl(pv)==1) {
                table_of_processes[i].push_back(6); //111b 111c
                table_of_end_pos[i].push_back(pv);
            }
            else if(s.get_TOPl(i)==0 && s.get_TOPl(pv)==0) {
                table_of_processes[i].push_back(31); //100c 100b
                table_of_end_pos[i].push_back(pv);
            }
        }
    }
}
//ci vuole una funzione che aggiorni quella qui sopra quando viene introdotto anche il secondo strato
//sarà una funzione che gira solo su siti di bordo faccetta 100
//function for update the atoms position
void metropolis::second_layer_updates(){

}

int * metropolis::barrier_chooser(int cl){
    std::string num;
    int* bar = new int[2];

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

    return bar;
}


//calculate 
double metropolis::probability_calculator(int cl){
    int* par=barrier_chooser(cl);
    return ( nu_0 * (exp((-1.) *(barriere[par[0]] + (par[1]*E_b)) / (Kb*temperatura) )) );
}

void metropolis::map_of_class_eraser(int c, int i, int j){
    for(int k=0 ; k< map_class_processes[c].size() ; k++){
        if(map_class_processes[c][k][0] == i && map_class_processes[c][k][1] == j){
            map_class_processes[c].erase(map_class_processes[c].begin() + k);
            break;
        }
    }
}

//update the table_of_processes for each site nn to what happened
//update how many atoms per class
//the only advantage of having a table of processes is that you don't have to re-calculate the possible move for each atom
void metropolis::classification(){
    for(int i=0 ; i<interested_sites.size() ; i++){
        std::vector<int> top = table_of_processes[i];
        for(int j=0 ; j<top.size() ; j++){
            int cl=top[j]; //process
            int pv=table_of_end_pos[i][j]; //associated pv
            if(!atoms[pv] && atoms[interested_sites[i]]){ //if pv empty && the site is full
                int n = nnn_atoms[pv]; //nn of associated pv
                table_of_processes[i][j] = cl + n;
                atoms_per_class[cl+n]++;
                std::vector<int> posit={i,j};
                map_class_processes[cl+n].push_back(posit);
                atoms_per_class[cl]--;
                map_of_class_eraser(cl,i,j);
            }
        }
    }
}

void metropolis::probability_filler(){
    for(int i=0 ; i<n_class ; i++){
        double p = probability_calculator(i);
        P[i]+= atoms_per_class[i]*p;
    }
}

void metropolis::deposition_func(){
    int p = initial_pos(rng) * crd.sites_size();
    atoms[p] = true;
    pos=p;
}


void metropolis::time_prob_calc(){
    //update the time
    double r_time=dist_time(rng);
    double P_tot=0.0;
    double P_dep=crd.sites_size()*F;
    for(int i=0 ; i<P.size() ; i++) P_tot+=P[i];
    double P_diff=P_tot;
    P_tot+=P_dep;
    time+= ( -(log(r_time)) / P_tot);

    //choice of the class
    double random = dist_time(rng);
    double r_p = random * P_tot;
    if(r_p>P_diff) deposition=true;
    if(!deposition){
        int chosen_class;
        double p_i = 0;
        for(int i = 0 ; i < P.size() ; i++){
            if( (p_i<r_p) && (r_p < (p_i+P[i])) ) chosen_class = i;
            p_i += P[i];
        }
        //choice of the process
        std::uniform_int_distribution<> dist(0,atoms_per_class[chosen_class]-1);
        int atom = dist(rng);
        pos=map_class_processes[chosen_class][atom][0];
        next=map_class_processes[chosen_class][atom][0];
    }else deposition_func();
}


void metropolis::print_output(){
    //pensare a come fare il file di output
    output<< process_name << "  from    " << pos <<"  to    "<< next <<std::endl;
}

void metropolis::start_of_the_sim(){
    table_of_processes_filler(); //initialize
    atoms.resize(crd.get_N(), false);
    deposition_func();
    interested_sites_calc(false);
    output.open("MC_processes.txt");
}

void metropolis::algorithm(){
    //compute all the possible processes
    while(time<total_time){
        nn_updater(deposition);
        classification();
        probability_filler();
        time_prob_calc();
        if(deposition){ 
            deposition_func;
            crd.output_writer(output,crd.sites_size());
        }
        interested_sites_calc(deposition);
    }
}

void metropolis::simulation(){
    //questa deve richiamare tutte le funzioni sopra, deve essere L'UNICO comando che è usabile dall'esterno, ovvero pubblico 
    //assieme ovviamente a quelli che generano i file di output
    start_of_the_sim();
    algorithm();
    print_output();
}

