#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

/* PLATINUM DATA
*lattice constant: 3.9242 ang
distance between nearest neighbors: 2.7748 ang

THIS VERSION CREATES A TRUNCATED OCTAHEDRON AND RETURNS A SECOND FILE WITH THE SITES WHERE THE ADATOM MOVES
COMPARED TO VERSION 3, IT IS CREATED WITH THE Z AXIS INSTEAD OF THE X AXIS
*/

vector<vector<double>> square_maker(int l,double* c, double d_pv){
    vector<vector<double>> square;

    // Create at the origin and then translate
    for(int i=0;i<l;i++){
        for(int j=0;j<l;j++){
            vector<double> coord={d_pv*j,d_pv*i,0};
            square.push_back(coord);
        }
    }

    // Translate to the origin
    vector <double> cm={0,0,0};

    for(int t=0;t<(l*l);t++){
        cm[0]+=(square[t][0]/double(l*l));
        cm[1]+=(square[t][1]/double(l*l));
        cm[2]+=(square[t][2]/double(l*l));
    }
    for(int i=0;i<(l*l);i++){
        square[i][0]=square[i][0]-cm[0];
        square[i][1]=square[i][1]-cm[1];
        square[i][2]=square[i][2]-cm[2];        
    }

    // Translate to the desired position
    for(int i=0;i<(l*l);i++){
        square[i][0]=square[i][0]+c[0];
        square[i][1]=square[i][1]+c[1];
        square[i][2]=square[i][2]+c[2];        
    }
    return square;
}

vector<vector<double>> cutted_square_maker(int l,double* c, double d_pv,int n_cut){
    vector<vector<double>> square=square_maker(l,c,d_pv);
    vector<vector<double>> cutted;
    for(int i=0;i<square.size();i++){
        if(((fabs(square[i][1])+fabs(square[i][0]))<(double(l-n_cut)*(d_pv-(0.05*d_pv))))){ // Translate the line by a fraction of the step for safety, otherwise truncation caused issues
            cutted.push_back(square[i]);
        }
    }
    return cutted;
}

bool equal_rows(vector<double> a,vector<double> b,double tolerance){
    for(int i=0;i<a.size();i++){
        if(abs(a[i]-b[i])>tolerance){
            return false;
        }
    }
    return true;
}

vector<vector<double>> rotate_xy_45deg(vector<vector<double>> v){
    vector<vector<double>> rotated(v.size(),vector<double>(3,0.0));
    double a=sqrt(2)/2.;

    for(int i=0;i<v.size();i++){
        rotated[i][0]= a*v[i][0]-a*v[i][1];
        rotated[i][1]=a*v[i][0]+a*v[i][1];
        rotated[i][2]=v[i][2];
    }
    return rotated;
}

int main(){

    // DATA
    ifstream ifile("dati_piramide.txt");
    double d_pv;
    int n_l;
    int n_cut[4]; 
    ifile>>d_pv>>n_l;
    for(int i=0;i<4;i++){
        ifile>>n_cut[i];
    }
    double h= d_pv*sin(M_PI/4.);
    // Set the center vector of the square to be positioned
    double* c=new double[3];
    for (int i=0;i<3;i++) c[i]=0;

    // To create the border, the program runs twice; in the second run, two are added
    // At the end, for the border, the two results are compared, and all common atoms are discarded
    
    // Vectors to store everything
    vector<vector<double>> small;
    vector<vector<double>> large;

    //Having 1 or more Au layers on bottom
    int n_layers;
    cout << "How many Au layers do you want on the base?"<< endl;
    cin>>n_layers;

    for(int cycle=0;cycle<2;cycle++){

        if(cycle==1) {
            n_l=n_l+2;
            for(int p=0;p<3;p++) n_cut[p]=n_cut[p]+1;
        }

        vector<vector<double>> atoms;

        // Create the first one centered at the origin
        vector<vector<double>> q=cutted_square_maker(n_l,c,d_pv,n_cut[2]);

        for(int i=0;i<q.size();i++){
            atoms.push_back(q[i]);
        }

        // Create the others in the positive semi-axis THIS TIME FOR Z, NOT X
        int j=1;
        while(j<(n_l-n_cut[0])){
            int l=n_l-j;
            c[2]=j*h;

            // Add truncated squares as long as it makes sense
            if(l>(n_l-n_cut[2])){
                vector<vector<double>> bmb=cutted_square_maker(l,c,d_pv,n_cut[2]-j);
                for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
                j++;
            }else{
                vector<vector<double>> bmb=square_maker(l,c,d_pv);
                for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
                j++;
            }
        }

        // Create the negative axis
        int y=1;
        while(y<(n_l-n_cut[1])){
            int q=n_l-y;
            c[2]=(-1)*y*h;

            // Add truncated squares
            if(q>(n_l-n_cut[2])){
                vector<vector<double>> bmb=cutted_square_maker(q,c,d_pv,n_cut[2]-y);
                for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
                y++;
            }else{
                vector<vector<double>> bmb=square_maker(q,c,d_pv);
                for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
                y++;
            }
        }
        cout<<atoms.size()<<endl;

        if(cycle==0) small=atoms;
        if(cycle==1) large=atoms;

        for (int i=0;i<3;i++) c[i]=0;
    }
    

    // Save the indices to be removed
    vector<int> index;
    for(int j=0;j<large.size();j++){
        for(int i=0;i<small.size();i++){
            if(equal_rows(large[j],small[i],0.001)==true){
                index.push_back(j);
                cout<<j<<endl;
                cout<<i<<endl;
                cout<<large[j][0]<<" "<<large[j][1]<<" "<<large[j][2]<<endl;
                cout<<small[i][0]<<" "<<small[i][1]<<" "<<small[i][2]<<endl;
            }
        }
    }

    // Sort the indices in descending order because it is convenient
    // to remove them from the highest, otherwise the index of interest changes during the operation
    sort(index.rbegin(),index.rend());

    vector<vector<double>>::iterator it = large.begin();
    for(int i=0;i<index.size();i++){
        large.erase(it+index[i]);
    }

    //Add an Au layer on bottom of the pyramid
    int core_size = small.size();
    int l_add = 0;
    for(int i=0 ; i<n_layers; i++){
        l_add = n_cut[1]-1 - i;
        c[2] = (-1.)*(n_l-n_cut[1]-1+i)*h;
        for (int j=0 ; j<(l_add*l_add) ; j++) small.push_back(square_maker(l_add,c,d_pv)[j]);
    }

    //adsorption sites active from t>0 if the growth is 3D

    int n_high_holes;
    int slab = n_cut[3] -1 ;
    for (int i=0 ; i< slab; i++){
        n_high_holes += pow((n_cut[1]+1-i),2);
    }
    
    vector<vector<double>> high_holes;
    for(int i=0 ; i< slab ; i++){
        int q = n_cut[1] - i;
        c[2] = (-1)*(n_l-n_cut[1]+i) *h;

        for(int j=0 ; j<(q*q) ; j++) high_holes.push_back(square_maker(q,c,d_pv)[j]);
    }

    // BEFORE PRINTING, ROTATE BY 45 DEGREES IN THE XY PLANE SO IT HAS THE CORRECT ANGLES
    vector<vector<double>> r_large=rotate_xy_45deg(large);
    vector<vector<double>> r_small=rotate_xy_45deg(small);   
    vector<vector<double>> r_sites=rotate_xy_45deg(high_holes); 

    ofstream ofile;
    ofile.open("piramide.xyz");
    ofile<<r_small.size()<<endl;
    ofile<<"Pt "<<"Pt "<<"1 "<<"1"<<endl;
    for(int i=0;i<core_size;i++){
        //cout<<i<<endl;
        ofile<<"Pt "<<r_small[i][0]<<" "<<r_small[i][1]<<" "<<r_small[i][2]<<endl;
    }
    for(int i=0 ; i<(l_add*l_add); i++){
        ofile<< "Au " <<r_small[i+core_size][0]<<" "<<r_small[i+core_size][1]<<" "<<r_small[i+core_size][2]<<endl;
    }

    ofstream ofile1;
    ofile1.open("siti_jmol.xyz");
    ofile1<<r_large.size()<<endl;
    ofile1<<"Au "<<"Au "<<"1 "<<"1"<<endl;
    for(int i=0;i<r_large.size();i++){
        ofile1<<"Au "<<r_large[i][0]<<" "<<r_large[i][1]<<" "<<r_large[i][2]<<endl;
    }

    ofstream ofile2;
    ofile2.open("high_sites.xyz");
    ofile2<<r_sites.size()<<endl;
    ofile2<<"Au "<<"Au "<<"1 "<<"1"<<endl;
    for(int i=0;i<r_sites.size();i++){
        ofile2<<"Au "<<r_sites[i][0]<<" "<<r_sites[i][1]<<" "<<r_sites[i][2]<<endl;
    }

}
