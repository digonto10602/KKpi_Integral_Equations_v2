/* This file is to test the old code 
and check with the 2+1 system code */

#include<bits/stdc++.h>
#include "../old_source/functions_momrep_based.h"
#include "../old_source/momentum_vector_maker.h"
#include "../old_source/integralequation_momrep_based.h"
#include "../old_source/interpolator_momrep.h"
#include "../old_source/solvers.h"
#include <omp.h>
#include <Eigen/Eigenvalues> //for vertex factor calculations
//#include "../old_source/LUSolver_modified.h"
#include "../old_source/fastGL_based_functions.h"
#include "../old_source/SA_method_functions.h"
#include <sys/time.h>

typedef complex<double> comp; 

void check_degen_integral_eqn_component()
{
    double m = 1.0; 
    double a = 2.0; 
    comp sigb = sigmab(a,m); 

    std::cout<<"sig_b = "<<sigb<<std::endl; 

    comp s = 8.2; 
    comp En = std::sqrt(s); 
    comp qb = pmom(s, sigb, m); 

    std::cout<<"q_b = "<<qb<<std::endl; 

    comp kmin = 0.0; 
    comp kmax = pmom(s, 0.0, m); 

    std::cout<<"kmax = "<<kmax<<std::endl; 

    

    std::vector<comp> qvec; 
    std::vector<comp> weights; 
    double number_of_points = 1000; 

    double eps_for_m2k = 0.0; 
    double eps_for_ope = 0.0; 
    double eps_for_cutoff = 0.0; 
    double eps_eta = 25.0; 
    double eps_check = eps_momrep_above_threshold(eps_eta, s, a, m, number_of_points); 
    eps_for_m2k = eps_check; 
    eps_for_ope = eps_for_m2k; 
    std::cout<<"eps_for_m2k = "<<eps_for_m2k<<std::endl; 
    
    line_maker_with_weights(qvec, weights, kmin, kmax, number_of_points); 

    for(int i=0; i<qvec.size(); ++i)
    {
        std::cout<<"p="<<qvec[i]<<'\t'
                 <<"w="<<weights[i]<<std::endl;
    }

    int size = qvec.size(); 
    Eigen::MatrixXcd Bmat(size, size); 
    Eigen::VectorXcd Gvec(size); 
    Eigen::VectorXcd dvec(size); 

    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps_for_ope,eps_for_m2k);
    double relerror;

    std::cout<<"Bmat = "<<std::endl; 
    //std::cout<<Bmat<<std::endl; 
    for(int i=0; i<size; ++i)
    {
        for(int j=0; j<size; ++j)
        {
            std::cout<<i<<'\t'
                     <<j<<'\t'
                     <<Bmat(i,j)<<std::endl;
        }
    }

    Gvec_maker_momrep(Gvec, s, qvec, qb, m, eps_for_ope); 
    Gvec = -1.0*Gvec; 

    std::cout<<"Gvec = "<<std::endl; 
    for(int i=0; i<size; ++i)
    {
        std::cout<<i<<'\t'<<qvec[i]<<'\t'<<Gvec(i)<<std::endl; 
    }

    double relerr; 
    LinearSolver_2(Bmat, dvec, Gvec, relerr); 

    std::cout<<"rel err = "<<relerr<<std::endl; 
    std::cout<<"dmat = "<<std::endl; 
    for(int i=0; i<size; ++i)
    {
        std::cout<<i<<'\t'<<dvec(i)<<std::endl; 
    }

    comp result; 
    interpolator_ds_integraleq_momrep_2eps_with_weights(dvec, qvec, weights, s, qb, qb, a, m, eps_for_ope, eps_for_m2k, result);

    comp gfunc = gfuncConst(a, 0.0, m); 
    comp gsq = gfunc*gfunc; 
    comp ds = result; 
    std::cout<<"g = "<<gfunc<<std::endl;
    comp mphib = gsq*result; 
    std::cout<<"mphib = "<<mphib<<std::endl; 
    comp rhopb = rhophib(s, a, m); 
    double error = abs(((1.0/mphib).imag() + rhopb)/rhopb)*100.0; 
    std::cout<<"im(1/mphib) = "<<(1.0/mphib).imag()<<std::endl; 
    std::cout<<"rho phib = "<<rhopb<<std::endl;
    std::cout<<"error = "<<error<<std::endl; 

}

int main()
{
    check_degen_integral_eqn_component();

    return 0; 

}