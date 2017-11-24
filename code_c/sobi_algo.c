// -----------------------------------------------------------------//
// Code property:                                                   //
// ==============                                                   //
// TELECOM BRETAGNE						                            //
// Dpt. Signal et Communications				                    //
// Technopole Brest-Iroise					                        //
// CS 83818 - 29238 Brest Cedex 3 - France			                //
//                                                                  //
// General info:                                                    //
// =============                                                    //
// program: sobi_algo.c                                             //
// last update: 14/10/2017                                          //
// info: yi.xu@imt-atlantique.net                                   //
//                                                                  //
// Code info:  	                                                    //
// ==========                                                       //
// A code for calculating matrix A for X=AS with 2 mixzed signals   //
//------------------------------------------------------------------//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "math_functions.h"

double* sobi_algo(double * input_1, double * input_2, int nb, int size_window){
    int j;
    int k;
    int i;
    static double A[4];
    double R_1[size_window][size_window];
    double R_2[size_window][size_window];
    double R_12[size_window][size_window];
    double T_1 = 0.0;
    double T_2 = 0.0;
    double T_12 = 0.0;
    double F_1 = 0.0;
    double F_2 = 0.0;
    double F_12 = 0.0;
    int nb_blocks;
    center(input_1, nb);
    center(input_2, nb);

    // initialization of R_x, NOT EVERY VALUE is initialized to zero !!! sometimes -nan !!!!
    for(i=0; i<size_window; i++){
        for(j=0; j<size_window; j++){
            R_1[i][j] = 0.0;
            R_2[i][j] = 0.0;
            R_12[i][j] = 0.0;
        }
    }
    nb_blocks = nb/size_window;

    for(k=0; k<nb_blocks; k++){
        int begin = size_window*k;
        for(i=0; i<size_window; i++){
            for(j=0; j<size_window; j++){
                R_1[i][j] =  R_1[i][j] + input_1[begin+i]*input_1[begin+j]/nb_blocks;
                R_2[i][j] = R_2[i][j] + input_2[begin+i]*input_2[begin+j]/nb_blocks;
                R_12[i][j] = R_12[i][j] + input_1[begin+i]*input_2[begin+j]/nb_blocks;
            }
        }
    }
    // calculate parameters
    for(i=0; i<size_window; i++){
                T_1 += R_1[i][i];
                T_2 += R_2[i][i];
                T_12 += R_12[i][i];
        for(j=0; j<size_window; j++){
                F_1 += R_1[i][j];
                F_2 += R_2[i][j];
                F_12 += R_12[i][j];
        }
    }

    F_1  = (F_1 - T_1)/(size_window)/(size_window - 1.0);
    T_1 /= size_window;
    F_2  = (F_2 - T_2)/(size_window)/(size_window - 1.0);
    T_2 /= size_window;
    F_12  = (F_12 - T_12)/(size_window)/(size_window - 1.0);
    T_12 /= size_window;

    printf("F1 = %f\n",  F_1);
    printf("F2 = %f\n",  F_2);
    printf("F12 = %f\n",  F_12);
    printf("T1 = %f\n",  T_1);
    printf("T2 = %f\n",  T_2);
    printf("T12 = %f\n",  T_12);
    double Alpha = 2*F_12*T_12 - (F_1*T_2 + F_2*T_1);
    double Beta = 2*(T_12*T_12 - T_1*T_2);
    double gamma_2 = (F_1*T_2 - F_2*T_1)*(F_1*T_2 - F_2*T_1) + 4*(F_12*T_2-T_12*F_2)*(F_12*T_1 - T_12*F_1);
    double d_1 = Alpha - sqrt(gamma_2);
    double d_2 = Alpha + sqrt(gamma_2);

    A[0] = Beta*F_1-T_1*d_1; // A11
    A[1] = Beta*F_12-T_12*d_2; // A12
    A[2] = Beta*F_12 - T_12*d_1; // A21
    A[3] = Beta*F_2-T_2*d_2; // A22
    return A;

}


