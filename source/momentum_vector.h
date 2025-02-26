#ifndef MOMENTUM_VECTOR_H
#define MOMENTUM_H

#include<bits/stdc++.h>
#include "fastGL_based_functions.h"
#include "functions_mom_based.h"

typedef complex<double> comp;

void flavor_based_momentum_vector(  std::vector<comp> &qvec,
                                    std::vector<comp> &weights, 
                                    comp En, 
                                    double mi, //mass of the spectator
                                    double number_of_points )
{
    comp epsilon_for_kvec = 0;// 1.0e-10; 
    comp kmin = 0.0 + epsilon_for_kvec; 
    comp kmax = pmom(En, 0.0, mi) - epsilon_for_kvec;
    line_maker_with_weights(qvec, weights, kmin, kmax, number_of_points); 
}

void simple_momentum_vector(    std::vector<comp> &qvec,
                                std::vector<comp> &weights, 
                                comp kmin,
                                comp kmax,
                                double number_of_points )
{
    comp epsilon_for_kvec = 0;// 1.0e-10; 
    //comp kmin = 0.0 + epsilon_for_kvec; 
    //comp kmax = pmom(En, 0.0, mi) - epsilon_for_kvec;
    line_maker_with_weights(qvec, weights, kmin, kmax, number_of_points); 
}


/* This contour goes below the real momentum 
axis in such a way that it misses the two bound-state
poles for M2 with m1 and m2 particles */
void contour_1( std::vector<comp> &qvec, 
                std::vector<comp> &weights, 
                comp kmin, 
                comp qb_1, 
                comp qb_2, 
                comp kmax, 
                double eps,
                double number_of_points )
{
    comp ii = {0.0, 1.0}; 
    double size1 = number_of_points/3.0; 
    comp first_point = qb_1 - ii*eps; 
    comp second_point = qb_2 - ii*eps; 

    line_maker_with_weights(qvec, weights, kmin, first_point, size1); 
    line_maker_with_weights(qvec, weights, first_point, second_point, size1);
    line_maker_with_weights(qvec, weights, second_point, kmax, size1); 

}

/* This contour goes below the real momentum axis
below on single two-body bound-state pole, this is 
for M2 with m1 and m1 particles */
void contour_2( std::vector<comp> &qvec, 
    std::vector<comp> &weights, 
    comp kmin, 
    comp qb_1, 
    comp kmax, 
    double eps,
    double number_of_points )
{
comp ii = {0.0, 1.0}; 
double size1 = number_of_points/2.0; 
comp first_point = qb_1 - ii*eps; 

line_maker_with_weights(qvec, weights, kmin, first_point, size1); 
line_maker_with_weights(qvec, weights, first_point, kmax, size1); 

}


#endif