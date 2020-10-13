#include<bits/stdc++.h>
using namespace std ; 

int main(void){

    long double DELTA = 1 ; 

    while(DELTA < 2){    
    DELTA += 0.01 ; 
    long long TIME = 500000 ;
    int  NUM = 3 ; 
    long long B = 500000 ; 
    long double DT = 0.005 ;
   
    vector<vector<long double > > D(3 , vector<long double>(NUM)) ;  
    long double EPS = 1.0 , SIGMA = 1.0 , RHO = 1.0; 
    vector<long double> MASS(NUM) ; 
    for(int i = 0 ; i < NUM ; i++) MASS.at(i) = 3.0 ; 
    vector<vector<long double > > force1(3 , vector<long double>(NUM)), force2(3 , vector<long double>(NUM)) ; 
    vector<vector<long double > > engp(NUM , vector<long double>(NUM)) ; 
    long double x = 0 , y=0  , z=0  , fx=0 , fy=0 , fz=0 , ep=0 , rxy=0 , r = 0 , fc = 0  ;
    vector<long double> ke(TIME) , U(TIME) , E(TIME);  
    //long double r01_vec[3] , r02_vec[3] ,C ;
    long double T_mean = 0 , E_mean = 0 , ke_mean = 0 ;

    for(int i = 0 ; i < 3 ; i++)for(int j = 0 ; j < NUM ; j++){
      force1.at(i).at(j) = 0.0 ;
      force2.at(i).at(j) = 0.0 ;
      engp.at(i).at(j) = 0.0 ; 
    } 


/////////////////////////////////////////////initial position//////////////////////////////////////////////////
    vector<vector<long double > > pos(3 , vector<long double>(NUM)) ;   


    D[0][0] =  0.484507151549531 ;
    D[0][1] = -0.728112699618685 ;
    D[0][2] =  0.339269852827254 ;
    D[1][0] =  0.649266472076588 ;
    D[1][1] = -0.608513580689966 ;
    D[1][2] = -0.880300767367624 ;
    D[2][0] =  0.0 ;
    D[2][1] =  0.0 ;
    D[2][2] =  0.0 ;

    pos.at(0).at(0) = 0.0 ;//原子１
    pos.at(1).at(0) = 0.0 ;
    pos.at(2).at(0) = 0.0 ;
    
    pos.at(0).at(1) = pow(2.0 , (1.0)/(6.0))*SIGMA  ;
    pos.at(1).at(1) = 0.0 ;
    pos.at(2).at(1) = 0.0 ;
    
    pos.at(0).at(2) = pow(2.0,(1.0/6.0))*SIGMA/2.0  ;
    pos.at(1).at(2) = pow(2.0,(1.0/6.0))*SIGMA*sin(M_PI/3.0) ;
    pos.at(2).at(2) = 0.0 ;

    for(int i = 0 ; i < 3 ; i++)for(int j = 0 ; j < NUM ; j++) pos.at(i).at(j) += D.at(i).at(j)*DELTA ; 


//////////////////////////////////////// initial velocity ///////////////////////////////////////////////

    vector<vector<long double > > v(3 , vector<long double> (NUM)) ; 
    for(int i = 0 ; i < 3 ; i++)for(int j = 0 ; j < NUM ; j++) v.at(i).at(j) = 0.0 ; 
        
        
/////////////////////////////////////  verlet algorithm   /////////////////////////////////////////////
        
        
    for (long long  k = 0 ; k < TIME ; k++){

    ////////////////////////////////////// calc force 1 //////////////////////////////

      for(int i = 0 ; i < NUM ; i++){
        for(int j = i+1 ; j < NUM ; j++){

          x = pos.at(0).at(i) - pos.at(0).at(j);//原子間のx軸距離
          y = pos.at(1).at(i) - pos.at(1).at(j);
          //z = pos.at(2).at(i) - pos.at(2).at(j);
          rxy  = sqrt(x*x+y*y) ; 
          r = sqrt(x*x+y*y+z*z) ;
          fc = (long double) -(-2.0*RHO*exp(-2.0*RHO*(r-1.0)) + 2.0*RHO*exp(-RHO*(r-1.0))) ; 
          fx = fc * (x/rxy) ; 
          fy = fc * (y/rxy) ;
          //fz = (long double) 2.0*RHO*exp(-2.0*RHO*(z-1.0)) - 2*RHO*exp(-RHO*(z-1.0)) ;

          force1.at(0).at(i) += fx ;
          force1.at(1).at(i) += fy ;
          //force1.at(2).at(i) += fz ;

          force1.at(0).at(j) -= fx ;
          force1.at(1).at(j) -= fy ;
          //force1.at(2).at(j) -= fz ;
          //cout << x << " " << " " << y << " " << z << endl; 

          x = 0.0 ; 
          fx = 0.0 ; 
          y = 0.0 ; 
          fy = 0.0 ; 
          z = 0.0 ; 
          fz = 0.0 ; 
              //cout << force1.at(0).at(i) << endl ;
              //cout << force1.at(0).at(j) << endl ;
        }
    }
    ////////////////////////////////// 位置の更新 ///////////////////////////
    for(int i = 0 ; i < NUM ; i++){
      pos.at(0).at(i) = pos.at(0).at(i) + DT*v.at(0).at(i) + (0.5*DT*DT/MASS.at(i)*force1.at(0).at(i));
      pos.at(1).at(i) = pos.at(1).at(i) + DT*v.at(1).at(i) + (0.5*DT*DT/MASS.at(i)*force1.at(1).at(i));
      pos.at(2).at(i) = pos.at(2).at(i) + DT*v.at(2).at(i) + (0.5*DT*DT/MASS.at(i)*force1.at(2).at(i));
    }    
///////////////////////////////////////calc force 2 ////////////////////////////////////////

    for(int i = 0 ; i < NUM ; i++){
      for(int j = i+1 ; j < NUM ; j++){
        x = pos.at(0).at(i) - pos.at(0).at(j);
        y = pos.at(1).at(i) - pos.at(1).at(j);
        //z = pos.at(2).at(i) - pos.at(2).at(j);
        r = sqrt(x*x+y*y+z*z) ;
        rxy  = sqrt(x*x+y*y) ; 

        fc = (long double) -(-2.0*RHO*exp(-2.0*RHO*(r-1.0)) + 2.0*RHO*exp(-RHO*(r-1.0))) ; 
        fx = fc * (x/rxy) ; 
        fy = fc * (y/rxy) ;
        
        ep = exp(-2*RHO*(r-1.0)) - 2.0*exp(-RHO*(r-1.0)); 
        engp.at(i).at(j) = ep ;
        force2.at(0).at(i) += fx ;
        force2.at(1).at(i) += fy ;
        //force2.at(2).at(i) += fz ;

        force2.at(0).at(j) -= fx ;
        force2.at(1).at(j) -= fy ;
        //force2.at(2).at(j) -= fz ;

        x = 0 ; 
        fx = 0 ; 
        y = 0 ; 
        fy = 0 ; 
        z = 0 ; 
        fz = 0 ;    
        ep = 0 ;      
      }
    }

////////////////////////////////////////calc velocity /////////////////////////         
    for(int i = 0 ; i < NUM ; i++){
      v.at(0).at(i) = v.at(0).at(i) + ( (DT/MASS.at(i) )*(force1.at(0).at(i) + force2.at(0).at(i) ) ) /2.0 ;
      v.at(1).at(i) = v.at(1).at(i) + ( (DT/MASS.at(i) )*(force1.at(1).at(i) + force2.at(1).at(i) ) ) /2.0 ;
      v.at(2).at(i) = v.at(2).at(i) + ( (DT/MASS.at(i) )*(force1.at(2).at(i) + force2.at(2).at(i) ) ) /2.0 ;
    }    
//////////////////////////////////////////////////////////////////////

/*
        for (i = 0 ; i<=2 ; i++){
        
        r01_vec[i] = pos[i][1] - pos[i][0]  ;
        r02_vec[i] = pos[i][2] - pos[i][0]  ;
        
        }

        i = 0 ;
        j = 0 ;

        C = r01_vec[0]*r02_vec[1] - r01_vec[1]*r02_vec[0] ; 
        
*/        


    for(int i = 0 ; i < 3 ; i++){
      for(int j = 0 ; j < NUM ; j++){
        //cout << ke.at(k) << endl ;
        ke.at(k) += (1.0/2.0) * MASS.at(j)*v.at(i).at(j)*v.at(i).at(j) ;
      }
    }


    for(int i = 0 ; i < NUM ; i++){
      for(int j = 0 ; j < NUM ; j++){
        U.at(k) +=  engp.at(i).at(j) ;
        engp.at(i).at(j) = 0.0 ;
      }
    }
    //cout << U.at(k) << " " << ke.at(k) << endl;
    E.at(k) = U.at(k) + ke.at(k) ;
            
            
    x = 0.0 ;
    y = 0.0 ;
    z = 0.0 ;
    ep = 0.0 ;
            
            
    if(B < TIME ){
        //cout << U.at(k) << " " <<ke.at(k) <<  " " << E.at(k) << endl;
        //cout << E.at(k) << endl ;
    }
            
            
        
    for(int i = 0 ; i < 3 ; i++){
      for(int j = 0 ; j < NUM ; j++){
        force1.at(i).at(j) = 0.0 ;
        force2.at(i).at(j) = 0.0 ;
        //v.at(i).at(j) = 0 ;
      }
    }
            
    //cout << "U = " << U.at(k) << endl; 
    }
    ///////////////////////////////////////ベルレ法終わり/////////////////////////////////
    //cout << E.at(TIME-1) << endl ;
        for(long long i = 0 ; i < TIME ; i++){
            E_mean += E.at(i) ;
            T_mean += ke.at(i) ;
            
        }
        
        E_mean /= TIME ;
        T_mean *= 2 ;
        T_mean /= 3*NUM - 6 ;
        T_mean /= TIME ;
        
        cout <<fixed << setprecision(10) << E_mean << "," << T_mean << endl ;
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
    }
    return (0);
}