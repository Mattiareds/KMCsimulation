#include "geometry.h"
#include "metropolis.h"
#include "coordinates.h"
#include <iostream>
#include <cmath>
#include <vector>


void geometry::make_planes(){
    //100 tipo 0
    //std::cout<<"Making planes with: "<<st->get_Ncut(0)<<" "<<st->get_Ncut(1)<<" "<<st->get_Ncut(2)<<std::endl;
    const int num_planes = 15;
    const int components = 4; // x, y, z, offset

    planes = new double*[num_planes];
    for (int i = 0; i < num_planes; ++i) {
        planes[i] = new double[components];
    }
    
    double factor = st->get_dpv() / std::sqrt(2.0);
    double spigolo_minus_1 = st->get_spigolo() - 1.0;

    planes[0][0] = 0.0;  planes[0][1] = 0.0;  planes[0][2] = 1.0;
    planes[0][3] = factor * (spigolo_minus_1 - st->get_Ncut(0));

    // Plane 1: Piano sotto
    planes[1][0] = 0.0;  planes[1][1] = 0.0;  planes[1][2] = -1.0;
    planes[1][3] = factor * (spigolo_minus_1 - st->get_Ncut(1));

    //Plane 14: Secondo strato
    planes[14][0] = 0.0;  planes[14][1] = 0.0;  planes[14][2] = -1.0;
    planes[14][3] = factor * (st->get_spigolo() - st->get_Ncut(1) );

    // Planes 2-5: Lateral planes (X and Y axes)
    // Indices: 2(+x), 3(-x), 4(+y), 5(-y)
    double ncut2_offset = factor * (spigolo_minus_1 - st->get_Ncut(2));
    
    planes[2][0] = 1.0;  planes[2][1] = 0.0;  planes[2][2] = 0.0;  planes[2][3] = ncut2_offset;
    planes[3][0] = -1.0; planes[3][1] = 0.0;  planes[3][2] = 0.0;  planes[3][3] = ncut2_offset;
    planes[4][0] = 0.0;  planes[4][1] = 1.0;  planes[4][2] = 0.0;  planes[4][3] = ncut2_offset;
    planes[5][0] = 0.0;  planes[5][1] = -1.0; planes[5][2] = 0.0;  planes[5][3] = ncut2_offset;

    // Planes 6-13: The "111" family (Diagonals)
    double diag_offset = factor * spigolo_minus_1;
    double signs[8][3] = {
        {1,1,1}, {-1,1,1}, {-1,-1,1}, {1,-1,1},   // Upper diagonals
        {-1,-1,-1}, {-1,1,-1}, {1,-1,-1}, {1,1,-1} // Lower diagonals
    };

    for (int i = 0; i < 8; ++i) {
        int idx = i + 6; // Start from index 6
        planes[idx][0] = signs[i][0];
        planes[idx][1] = signs[i][1];
        planes[idx][2] = signs[i][2];
        planes[idx][3] = diag_offset;
    }
}


bool geometry::test_piano(int i,int j){
    if(!st) return false;
    // compute signed value of plane equation
    double val = st->site_x(i)*planes[j][0] + st->site_y(i)*planes[j][1] + st->site_z(i)*planes[j][2] - planes[j][3];
    double nx = planes[j][0];
    double ny = planes[j][1];
    double nz = planes[j][2];
    double norm = std::sqrt(nx*nx + ny*ny + nz*nz);
    double dist = std::fabs(val) / (norm > 0.0 ? norm : 1.0);
    const double eps = 1e-3; // distance tolerance
    return dist < eps;
}

int plane_type(int i){
    int t;
    if(i<6){ //100
        t=0;
    }
    if(i>5 && i<14){ //111
        t=1;
    }
    if(i==14){
        t=i;
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
        for(int j=0; j<15 ; j++){
            if(test_piano(i,j)==true) {
                if(counter[i]==0){
                    info_plane_sites[i][1] = (double) j;
                    type_of_plane[i] = plane_type(j);
                    if(plane_type(j)==14){
                        upper_sites.push_back(i);
                    }
                }
                else if(counter[i]==1){
                    info_plane_sites[i][2]=(double) j;
                    if((plane_type(j) + type_of_plane[i]) == 2){ type_of_plane[i] = 2; }//two 111 
                    else if((plane_type(j) && type_of_plane[i]) ==1){ type_of_plane[i] = 3; }//one 111 and one 100
                    else if((plane_type(j) && type_of_plane[i]) ==0){ type_of_plane[i] = 4;}//two 100
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
    if(!st){
        std::cerr << "Warning: geometry::st not set in test_border\n";
        return;
    }
    const double eps = 1e-3; // distance tolerance
    for(size_t i=0 ; i<type_of_plane.size() ; i++){ //run on which is the plane for each site
        int info = type_of_plane[i];
      //sites only on one plane
      if(info < 2 || info==14){
            //for the planes of the core
            for(size_t j=0 ; j<14 ; j++){
                double val = st->site_x(i)*core.planes[j][0] + st->site_y(i)*core.planes[j][1] + st->site_z(i)*core.planes[j][2] - core.planes[j][3];
                double nx = core.planes[j][0];
                double ny = core.planes[j][1];
                double nz = core.planes[j][2];
                double norm = std::sqrt(nx*nx + ny*ny + nz*nz);
                double dist = std::fabs(val) / (norm > 0.0 ? norm : 1.0);
                if(dist < eps){
                    bool has_different_nn = false;
                    for(int nn : table_of_nn[i]){
                        if(type_of_plane[nn] != type_of_plane[i]){
                            has_different_nn = true; break;
                        }
                    }
                    if(!has_different_nn) break; // è interno, non border
                    if(info==14){
                        type_of_plane[i] = 15;
                    }else{
                        type_of_plane[i] = type_of_plane[i] + 7; //7 if 100b, 8 if 111b
                    }
                    break; // one plane match is enough
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
            if(i!=j && (sqrt( pow((st->site_x(j) - st->site_x(i)),2) + pow((st->site_y(j) - st->site_y(i)),2) + pow((st->site_z(j)-st->site_z(i)),2) ) - st->get_dpv() ) < 0.08 ){
                table_of_nn[j].push_back(i); 
            } 
        }
    }
}


// bool geometry::is_same_plane(int site_i, int site_f){
//     auto planes_i=info_plane_sites[site_i];
//     auto planes_f=info_plane_sites[site_f];
//     bool same_p=false;
//     for(size_t i=1 ; i<planes_i.size() ; i++){
//         for(size_t j=1 ; j<planes_f.size() ; j++){
//             if(planes_i[i]==planes_f[j]) { same_p=true; return same_p; }
//         }
//     }
//     return same_p;
// }

//implementare una funzione che sostituisca tutti i pv di edge con il pv più vicino della faccetta adiacente
void geometry::pv_substitution(){
    auto table_copy=table_of_nn;
    edge_map.resize(st->sites_size());
    for (size_t i=0 ; i<table_of_nn.size() ; i++){
        edge_map[i].resize(table_of_nn[i].size(), -1);  // Resize each row individually
    }
    //run on sites
    for (size_t i=0 ; i<table_of_nn.size() ; i++){
        //if my site is not on one edge
        if(type_of_plane[i] < 2){
            //get the nn of my site
            for (size_t j=0 ; j< table_of_nn[i].size(); j++){
                int possible_edge_site=table_of_nn[i][j];
                //my site is on the edge?
                if(type_of_plane[possible_edge_site] > 1){
                    //run on all the nearest neighbours of the site that is on the edge
                    for(size_t k=0 ; k< table_of_nn[possible_edge_site].size(); k++){
                        //this is the nearest neighbour of the nearest neighbour of my site
                        int nn_of_the_nn = table_of_nn[possible_edge_site][k];
                        //if is not on the edge, it means that it has only one plane
                        if(type_of_plane[nn_of_the_nn]<2){
                            //the first one that is good and that is not on the edge will be my new nn
                            if(info_plane_sites[i][1]!=info_plane_sites[nn_of_the_nn][1]) {
                                //std::cout<<"site: "<<i<<" SHELL substitution: "<<possible_edge_site<<" "<<table_of_nn[i][j]<< " with "<< nn_of_the_nn<<std::endl;
                                edge_map[i][j] = possible_edge_site;
                                table_of_nn[i][j] = nn_of_the_nn;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    //for (int i=0 ; i<table_of_nn.size(); i++){
    //    std::cout<<i<<" con pv:       ";
    //    for (int j: table_of_nn[i]) std::cout<<j<<" ";
    //    std::cout<<std::endl;
    //}
}


void geometry::starter(){
    make_planes();
    site_characterisation();
    initialize_nn();
}


