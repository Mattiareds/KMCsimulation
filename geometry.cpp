#include "geometry.h"
#include "metropolis.h"
#include "coordinates.h"
#include <iostream>


void geometry::make_planes(){
    //100 tipo 0
    std::cout<<"Making planes with: "<<st->get_Ncut(0)<<" "<<st->get_Ncut(1)<<" "<<st->get_Ncut(2)<<std::endl;
    planes.push_back({0,0,1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1-st->get_Ncut(0))}); //piano sopra
    planes.push_back({0,0,-1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1-st->get_Ncut(1))}); //piano sotto
    planes.push_back({1,0,0,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1-st->get_Ncut(2))});  
    planes.push_back({-1,0,0,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1-st->get_Ncut(2))}); 
    planes.push_back({0,1,0,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1-st->get_Ncut(2))});
    planes.push_back({0,-1,0,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1-st->get_Ncut(2))});
    //111 tipo 1
    planes.push_back({1,1,1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)});  //piani grandi sopra
    planes.push_back({-1,1,1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)});
    planes.push_back({-1,-1,1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)});
    planes.push_back({1,-1,1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)});
    planes.push_back({-1,-1,-1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)}); //piani piccoli sotto
    planes.push_back({-1,1,-1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)});
    planes.push_back({1,-1,-1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)});
    planes.push_back({1,1,-1,st->get_dpv()/(sqrt(2))*(st->get_spigolo()-1)});
}


bool geometry::test_piano(int i,int j){
    if(abs(st->site_x(i)*planes[j][0] + st->site_y(i)*planes[j][1] + st->site_z(i)*planes[j][2] - planes[j][3]) < 0.05){
        return true;
    }
    return false;
}

int plane_type(int i){
    int t;
    if(i<6 || i==14){ //100
        t=0;
    }
    if(i>5 && i<15){ //111
        t=1;
    }
    return t;
}


//to do, for each type of coordinate readed, BEFORE the simulation
void geometry::site_characterisation(){
    //which plane?
    info_plane_sites.resize(st->sites_size(), std::vector<int>(5,0));
    type_of_plane.resize(st->sites_size(), 0);
    std::vector<int> counter(info_plane_sites.size(), 0);
    for(int i=0; i<st->sites_size(); i++){
        for(int j=0;j<planes.size();j++){
            if(test_piano(i,j)==true) {
                if(counter[i]==0){
                    info_plane_sites[i][1] = (double) j;
                    type_of_plane[i] = plane_type(j);
                }
                else if(counter[i]==1){
                    info_plane_sites[i][2]=(double) j;
                    if((plane_type(j) + type_of_plane[i])== 2) type_of_plane[i] = 2; //two 111 
                    else if((plane_type(j) && type_of_plane[i]) ==1) type_of_plane[i] = 3; //one 111 and one 100
                    else if((plane_type(j) && type_of_plane[i]) ==0) type_of_plane[i] = 4; //two 100
                }
                else if(counter[i]==2){
                    info_plane_sites[i][3]=(double) j;
                    type_of_plane[i] = 5; //111 111 100 is the only case
                }
                else if(counter[i]==3){
                    info_plane_sites[i][4]=(double) j;
                    type_of_plane[i] = 6; //111 111 100 100 is the only case
                }
                counter[i]++;
            }
        }
        info_plane_sites[i][0]=counter[i];
    }
}


//test if the site of the lattice are on the planes of the base (= are border sites)
void geometry::test_border(geometry& core){
    for(size_t i=0 ; i<type_of_plane.size() ; i++){
      //mi restringo ai siti che stanno su un solo piano
      if(type_of_plane[i] == 0 || type_of_plane[i] == 1){
    	    for(size_t j=0 ; j<core.planes.size() ; j++){
    	    //testo se appartengo ai piani del più piccolo, e nel caso (ne basta uno) ritorno 1, vale sia per 111 che per 100
    	        if((abs(st->site_x(i)*core.planes[j][0] + st->site_y(i)*core.planes[j][1] + st->site_z(i)*core.planes[j][2] - core.planes[j][3])<0.10)){
                    type_of_plane[i] = type_of_plane[i] + 7; //7if 100b 8 if 111b
    	        }
    	    }
        }
    }   
}

//initialize a table of pv, that is going to be updated when we add higher sites
void geometry::initialize_nn(){
    table_of_nn.resize(st->sites_size(),std::vector<int>(0,0));
    for(int j=0;j<st->sites_size();j++){
        for(int i=0;i<st->sites_size();i++){
            if(i!=j && (sqrt( pow((st->site_x(j) - st->site_x(i)),2) + pow((st->site_y(j) - st->site_y(i)),2) + pow((st->site_z(j)-st->site_z(i)),2) ) - st->get_dpv() ) < 0.08 ) table_of_nn[j].push_back(i); 
        }
    }
}

//initialize a table of pv for higher sites
void geometry::initialize_external_nn(coordinates& core){
    table_of_nn.resize(st->sites_size(),std::vector<int>(0,0));
    for(int j=0;j<st->sites_size();j++){
        for(int i=0;i<core.sites_size();i++){
            if(i!=j) if((sqrt( pow((st->site_x(j) - core.site_x(i)),2) + pow((st->site_y(j) - core.site_y(i)),2) + pow((st->site_z(j)-core.site_z(i)),2) ) - st->get_dpv() ) < 0.08 ) table_of_nn[j].push_back(i);
        }
    }
}

//implementare una funzione che sostituisca tutti i pv di edge con il pv più vicino della faccetta adiacente
void geometry::pv_substitution(){
    
}

void geometry::starter(){
    planes.resize(0,std::vector<double>(0,0.0));
    make_planes();
    site_characterisation();
    initialize_nn();
}
