//
//  LJ4_variable_DT.cpp
//  
//
//  Created by ryoma Matsumura on 2020/09/22.
//

#include<bits/stdc++.h>
using namespace std ;

int main(void)
{
    double DT = 0.00001 ; //データを取るときは0.005でやること
    
    while(DT <= 0.05){
        double DELTA = 0.00002 ;
        long double E_val = 0;
        //DELTA += 0.0005 ;
        DT *= 2  ;
        long long TIME ;//総ステップ数
        TIME = 1000000;
        int NUM = 4 ;  //粒子数
        long long B ; //平衡化時間
        B = 50000 ;
        //double DT = 0.005 ; //データを取るときは0.005でやること
        vector<long double> E(TIME) ;//エネルギーのを格納する配列
        vector<vector<long double > > D(3 , vector<long double>(NUM)) ; //初期座標からのズレ比率
        double EPS = 1.0 , SIGMA = 1.0 ;
        double CE12 ,CE06 ,CF12 ,CF06 ;
        CE12 = 4.0*EPS*pow(SIGMA,12.0);
        CE06 = 4.0*EPS*pow(SIGMA,6.0) ;
        CF12 = CE12*12.0;
        CF06 = CE06*6.0;
        
        vector<double> MASS(NUM) ;
        for(int i = 0 ; i< NUM ; i++) MASS.at(i) = 1.0 ;
        vector<vector<long double > > force1(3 , vector<long double >(NUM)),force2(3 , vector<long double >(NUM)), engp(NUM , vector<long double >(NUM));
        long double x , y , z , r2 , r2i , r06i , r12i , fc , fx , fy , fz , ep ;
        vector<long double> ke(TIME) ;
        vector<long double> U(TIME) ;
        
        double T_mean = 0 , E_mean = 0 , ke_mean = 0  ;
        
        ////////////////////////////////初期座標の設定////////////////////////////
        
        
        
        vector<vector<double > > pos(3 , vector<double>(NUM)) ;
        
        
        D[0][0] =  0.484507151549531 ;
        D[0][1] = -0.728112699618685 ;
        D[0][2] =  0.339269852827254 ;
        D[1][0] =  0.649266472076588 ;
        D[1][1] = -0.608513580689966 ;
        D[1][2] = -0.880300767367624 ;
        D[2][0] =  0.0 ;
        D[2][1] =  0.0 ;
        D[2][2] =  0.0 ;
        
        pos.at(0).at(0) = -0.0118091 ;//原子１
        pos.at(1).at(0) = 0.0248446 ;
        pos.at(2).at(0) = 0.0696386 ;
        
        pos.at(0).at(1) = 1.11008 ;
        pos.at(1).at(1) = 0.0252804 ;
        pos.at(2).at(1) = 0.105516 ;
        
        pos.at(0).at(2) = 0.550963 ;
        pos.at(1).at(2) = 0.994693 ;
        pos.at(2).at(2) = 0.018625 ;
        
        pos.at(0).at(3) = 0.5205  ;
        pos.at(1).at(3) = 0.413303  ;
        pos.at(2).at(3) = 0.978302  ;
        
        for(int i = 0 ; i < 3 ; i++){
            for(int j = 0 ; j < NUM ; j++){
                pos.at(i).at(j) += D.at(i).at(j)*DELTA  ;
            }
        }
        
        
        //////////////////////////////////// 初期速度//////////////////////////////////
        
        vector<vector<double> > v(3 , vector<double>(NUM)) ;
        
        
        
        //////////////////////////////  verlet algorithm//////////////////////////////////////
        
        for( long long k = 0 ; k < TIME ; k++){
            
            ///////////////////////////////////////////////////////////////////////////////
            for(int i = 0 ; i < NUM ; i++){
                for(int j = i+1 ; j < NUM ; j++){
                    x = pos.at(0).at(i) - pos.at(0).at(j);
                    y = pos.at(1).at(i) - pos.at(1).at(j);
                    z = pos.at(2).at(i) - pos.at(2).at(j);
                    
                    
                    r2 = x*x + y*y + z*z ;
                    r2i = 1.0/r2 ;
                    r06i = r2i*r2i*r2i ;
                    r12i = r06i * r06i ;
                    fc = (CF12 * r12i - CF06 * r06i)*r2i ;
                    fx = fc*x ;
                    fy = fc*y ;
                    fz = fc*z ;
                    
                    force1.at(0).at(i) += fx ;
                    force1.at(1).at(i) += fy ;
                    force1.at(2).at(i) += fz ;
                    
                    force1.at(0).at(j) -= fx ;
                    force1.at(1).at(j) -= fy ;
                    force1.at(2).at(j) -= fz ;
                    
                    //cout << force1.at(0).at(i) << endl ;
                    //cout << force1.at(0).at(j) << endl ;
                }
            }
            ////////////////////////////////////////////////////////////////////////////
            /*
             for(int i = 0 ; i < 3 ; i++){
             for(int j = 0 ; j < NUM ; j++){
             //cout << force1.at(i).at(j) << " " ;
             }
             //cout << endl ;
             }
             */
            ///////////////////////////////////////位置の更新//////////////////////////////////////
            for(int i = 0 ; i < NUM ; i++){
                pos.at(0).at(i) = pos.at(0).at(i) + DT*v.at(0).at(i) + (0.5*DT*DT/MASS.at(i)*force1.at(0).at(i));
                pos.at(1).at(i) = pos.at(1).at(i) + DT*v.at(1).at(i) + (0.5*DT*DT/MASS.at(i)*force1.at(1).at(i));
                pos.at(2).at(i) = pos.at(2).at(i) + DT*v.at(2).at(i) + (0.5*DT*DT/MASS.at(i)*force1.at(2).at(i));
                
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////
            
            for(int i = 0 ; i < NUM ; i++){
                for(int j = i+1 ; j < NUM ; j++){
                    x = pos.at(0).at(i) - pos.at(0).at(j);
                    y = pos.at(1).at(i) - pos.at(1).at(j);
                    z = pos.at(2).at(i) - pos.at(2).at(j);
                    
                    
                    r2 = x*x + y*y + z*z ;
                    r2i = 1.0/r2 ;
                    r06i = r2i*r2i*r2i ;
                    r12i = r06i * r06i ;
                    fc = (CF12 * r12i - CF06 * r06i)*r2i ;
                    fx = fc*x ;
                    fy = fc*y ;
                    fz = fc*z ;
                    
                    ep = CE12*r12i - CE06*r06i ;
                    engp.at(i).at(j) = ep ;
                    
                    force2.at(0).at(i) += fx ;
                    force2.at(1).at(i) += fy ;
                    force2.at(2).at(i) += fz ;
                    
                    force2.at(0).at(j) -= fx ;
                    force2.at(1).at(j) -= fy ;
                    force2.at(2).at(j) -= fz ;
                }
                
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////
            
            for(int i = 0 ; i < NUM ; i++){
                v.at(0).at(i) = v.at(0).at(i) + ( (DT/MASS.at(i) )*(force1.at(0).at(i) + force2.at(0).at(i) ) ) /2.0 ;
                v.at(1).at(i) = v.at(1).at(i) + ( (DT/MASS.at(i) )*(force1.at(1).at(i) + force2.at(1).at(i) ) ) /2.0 ;
                v.at(2).at(i) = v.at(2).at(i) + ( (DT/MASS.at(i) )*(force1.at(2).at(i) + force2.at(2).at(i) ) ) /2.0 ;
            }
            
            ///////////////////////////////////////////////////////////////////////////////////////////////
            
            
            for(int i = 0 ; i < 3 ; i++){
                for(int j = 0 ; j < NUM ; j++){
                    //cout << ke.at(k) << endl ;
                    ke.at(k) += (1.0/2.0) * MASS.at(j)* v.at(i).at(j)*v.at(i).at(j) ;
                }
            }
            
            
            
            for(int i = 0 ; i < NUM ; i++){
                for(int j = 0 ; j < NUM ; j++){
                    U.at(k) +=  engp.at(i).at(j) ;
                    engp.at(i).at(j) = 0.0 ;
                }
            }
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
                    force1.at(i).at(j) = 0 ;
                    force2.at(i).at(j) = 0 ;
                    //v.at(i).at(j) = 0 ;
                }
            }
            
            
            E_val += abs( E.at(0)-E.at(k) );
        
            //cout << E_val << endl ;
            
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
        
        cout <<fixed << setprecision(16) <<  DT << "," << E_val/TIME << endl ;
        E_mean = 0 ;
        T_mean = 0 ;
        //E_val = 0 ;
        ///////////////////////////最後の座標の出力////////////
        /*
         for(int i =  0 ; i < 3 ; i++){
         for(int j = 0 ; j < NUM ; j++){
         //cout << pos.at(i).at(j) << "," ;
         }
         //cout << endl ;
         }
         */
        ///////////////////////////////////////////////////
    }
    
    
    
    return (0);
}
