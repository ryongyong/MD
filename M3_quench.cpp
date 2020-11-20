#include<bits/stdc++.h>
using namespace std ; 


void set_parameter(long long &TIME ,int  &NUM ,long double &DT ,long double &DELTA ,long double &EPS ,long double &SIGMA ,long double &RHO){
    TIME = 20000000 ;
    NUM = 3 ; 
    DT = 0.005 ;
    DELTA = 0.0 ;    
    EPS = 1.0 ;
    SIGMA = 1.0 ; 
    RHO = 3.0 ; 
} 

long double calc_pes_dist(vector<pair<long double , long double> > dist_x,vector<pair<long double , long double> > dist_y, vector<vector<long double > > pos){
    long double dist_x_dist = 0 ;
    long double dist = 0 ; 

    for(int i = 0 ; i < 3 ; i++){
        dist_x_dist += ( dist_x.at(i).first - pos.at(0).at(i) )*( dist_x.at(i).first - pos.at(0).at(i) ) ; 
        dist_x_dist += ( dist_y.at(i).second - pos.at(1).at(i) )*( dist_y.at(i).second - pos.at(1).at(i) ); 
    }
    dist = sqrt(dist_x_dist) ; 
    return dist ; 
}

long double calc_rad(vector<long double>r01_vec  ,vector<long double> r02_vec  ){
    long double X = 0 ; 
    long double Y = 0 ; 
    long double Z = 0 ; 

    for(int i = 0 ; i < 3 ; i++){
        X += r01_vec.at(i)*r02_vec.at(i) ; 
    }
    for(int i = 0 ; i < 3 ; i++){
        Y += r01_vec.at(i)*r01_vec.at(i) ; 
        Z += r02_vec.at(i)*r02_vec.at(i) ; 
    }
    
    X /= sqrt(Y)*sqrt(Z)  ; 
    X = acosl(X)*180.0/M_PI ; 

    return X ; 
}

long double calc_kinetic_energy(vector<vector< long double > > v , vector<long double > MASS , int NUM){
    long double x = 0 ; 
    for(int i = 0 ; i < 3 ; i++)for(int j = 0 ; j < NUM ; j++){
        x += (1.0/2.0) * MASS.at(j)*v.at(i).at(j)*v.at(i).at(j) ;
    }
    return x ; 
}

long double calc_potential_energy(vector<vector< long double > > &engp , int NUM){
    long double x = 0 ; 
    for(int i = 0 ; i < NUM ; i++){
        for(int j = 0 ; j < NUM ; j++){
            x +=  engp.at(i).at(j) ;
            engp.at(i).at(j) = 0.0 ;
        }
    }
    return x ; 
}

long double mean_time(vector<long double> X ){
    long double x = 0 ; 
    for(int i = 0 ; i < X.size() ; i++){
        x += X.at(i) ; 
    }
    x /= X.size() ; 
    return x ; 
}

void calc_pos(vector<vector<long double > > &pos , vector<vector<long double> > v ,vector<vector<long double> > force1 , vector<long double> MASS , long double DT , int NUM ){
    for(int i = 0 ; i < NUM ; i++)for(int j = 0 ; j < 3 ; j++){
      pos.at(j).at(i) = pos.at(j).at(i) + DT*v.at(j).at(i) + (0.5*DT*DT/MASS.at(i)*force1.at(j).at(i));
    }
}

void calc_velocity(vector<vector<long double> > &v,vector<vector<long double> >force1 ,vector<vector<long double> > force2 ,vector<long double> MASS , long double DT , int NUM){
    for(int i = 0 ; i < NUM ; i++)for(int j = 0 ; j < 3 ; j++){
        v.at(j).at(i) = v.at(j).at(i) + ( (DT/MASS.at(i) )*(force1.at(j).at(i) + force2.at(j).at(i) ) ) /2.0 ;
    }
}

void initial_velocity(vector<vector<long double> > &v  , int NUM){
    for(int i = 0 ; i < 3 ; i++)for(int j = 0 ; j < NUM ; j++) v.at(i).at(j) = 0.0 ; 
}

void initial_pos(vector<vector<long double> > &pos ,vector<vector<long double> > D, int NUM , long double SIGMA , long double DELTA ){

    D[0][0] = -1.0 ;
    D[0][1] =  1.0 ;
    D[0][2] =  0.0;
    D[1][0] = -1.0 ;
    D[1][1] = -1.0 ;
    D[1][2] =  1.0 ;
    D[2][0] =  0.0 ;
    D[2][1] =  0.0 ;
    D[2][2] =  0.0 ;

    pos.at(0).at(0) = 0.0 ;//原子１
    pos.at(1).at(0) = 0.0 ;
    pos.at(2).at(0) = 0.0 ;
    
    pos.at(0).at(1) = 1.0 ;//原子2
    pos.at(1).at(1) = 0.00001 ;
    pos.at(2).at(1) = 0.0 ;
    
    pos.at(0).at(2) = 2.0 ; 
    pos.at(1).at(2) = 0.0 ;
    pos.at(2).at(2) = 0.0 ;

    for(int i = 0 ; i < 3 ; i++)for(int j = 0 ; j < NUM ; j++) pos.at(i).at(j) += D.at(i).at(j)*DELTA ; 
    
}

void calc_force(vector<vector<long double> > &force1, vector<vector<long double > > pos, long double RHO , int NUM, bool ok  , vector<vector< long double > > &engp){

    for(int i = 0 ; i < NUM ; i++){
        for(int j = i+1 ; j < NUM ; j++){
                
            long double x = pos.at(0).at(i) - pos.at(0).at(j);//原子間のx軸距離
            long double y = pos.at(1).at(i) - pos.at(1).at(j);
            long double z = pos.at(2).at(i) - pos.at(2).at(j);
            long double rxy  = sqrt(x*x+y*y) ; 
            long double r = sqrt(x*x+y*y+z*z) ;
            long double fc = -(-2.0*RHO*exp(-2.0*RHO*(r-1.0)) + 2.0*RHO*exp(-RHO*(r-1.0))) ; 

            force1.at(0).at(i) += fc * (x/rxy) ;
            force1.at(1).at(i) += fc * (y/rxy) ;
            force1.at(0).at(j) -= fc * (x/rxy) ;
            force1.at(1).at(j) -= fc * (y/rxy) ;

            if(ok){
                double ep = exp(-2*RHO*(r-1.0)) - 2.0*exp(-RHO*(r-1.0)); 
                engp.at(i).at(j) = ep ;
            }
        }
    }
}

int main(void){
    
    long long TIME ;
    int  NUM ; 
    long double DT ;
    long double DELTA ;    
    long double EPS ;
    long double SIGMA ; 
    long double RHO ; 

    set_parameter(TIME , NUM , DT , DELTA , EPS , SIGMA , RHO) ; 

    vector<vector<long double > > D(3 , vector<long double>(NUM)) ;  
    vector<vector<long double > > force1(3 , vector<long double>(NUM)) ; 
    vector<vector<long double > > force2(3 , vector<long double>(NUM)) ; 
    vector<vector<long double > > engp(NUM , vector<long double>(NUM)) ; 
    vector<vector<long double > > v(3 , vector<long double> (NUM)) ;  
    vector<vector<long double > > pos(3 , vector<long double>(NUM)) ;   
    vector<long double> ke(TIME) , U(TIME) , E(TIME);  
    vector<long double> r01_vec(3) , r02_vec(3) ;
    long double T_mean = 0 , E_mean = 0 , ke_mean = 0 ;
    vector<long double> theta(TIME) ;  
    vector<pair<long double , long double> > dist_x(3) , dist_y(3) ; 
    long double dist_pes = 0 ; 

    for(int i = 0 ; i < 3 ; i++)for(int j = 0 ; j < NUM ; j++){
        force1.at(i).at(j) = 0.0 ;
        force2.at(i).at(j) = 0.0 ;
        engp.at(i).at(j) = 0.0 ; 
    } 

    vector<long double> MASS(NUM) ; 
    for(int i = 0 ; i < NUM ; i++) MASS.at(i) = 1.0 ;     


    initial_pos(pos , D , NUM , SIGMA , DELTA) ; 
    initial_velocity(v , NUM) ;     
        
///////////////////////////////////// start verlet algorithm  /////////////////////////////////////////////
        
    for(long long  k = 0 ; k < TIME ; k++){
        bool ok = false ; 

        for(int i = 0 ; i < 3 ; i++){
            dist_x.at(i).first = pos.at(0).at(i) ; //x
            dist_y.at(i).second = pos.at(1).at(i) ; //y
        }



        calc_force(force1 , pos , RHO , NUM ,ok , engp) ; //calc force1 
        ok = true ; 
        calc_pos(pos , v , force1 , MASS , DT  , NUM ) ; //calc position
        calc_force(force2 , pos , RHO , NUM , ok  , engp) ; //calc force2 
        calc_velocity(v ,force1 , force2 , MASS , DT ,NUM) ;  //calc velocity 
        dist_pes += calc_pes_dist(dist_x , dist_y , pos) ; 

        ke.at(k) = calc_kinetic_energy(v , MASS ,NUM) ; 
        U.at(k) = calc_potential_energy(engp , NUM) ; 
        E.at(k) = U.at(k) + ke.at(k) ; // E = Total energy 


        for (int i = 0 ; i < NUM ; i++){
            r01_vec.at(i) = pos.at(i).at(1) - pos.at(i).at(0)   ;
            r02_vec.at(i) = pos.at(i).at(2) - pos.at(i).at(0)   ;
        }        
        theta.at(k) = calc_rad(r01_vec , r02_vec ) ;  
        if(k % 1000 == 0 && k > 30000) cout << fixed << setprecision(20) << dist_pes << "," <<  U.at(k) << endl;            

        for(int i = 0 ; i < 3 ; i++){
            for(int j = 0 ; j < NUM ; j++){
                force1.at(i).at(j) = 0.0 ;
                force2.at(i).at(j) = 0.0 ;
                v.at(i).at(j) = 0.0 ; 
            }
        }

    }

//////////////////////////////////end verlet algorithm //////////////////////////////////////////////
    
    E_mean = mean_time(E) ; 
    ke_mean = mean_time(ke) ; 
        
    E_mean = 0 ;
    T_mean = 0 ;
/*
    for(int i = 0 ; i < 3 ; i++){
        for(int j = 0 ; j < NUM ; j++){
            cout << pos.at(i).at(j) << "," ;
        }
        cout << endl;  
    }
*/
    return (0);
}
