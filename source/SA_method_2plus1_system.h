/* SA Method for 2+1 system */

#ifndef SAMETHOD_2PLUS1_H
#define SAMETHOD_2PLUS1_H

#include<bits/stdc++.h>
//#include "momentum_vector.h"
//#include "integral_equations_2plus1_system.h"
#include "linear_solvers_cpu.h"
//#include "gpu_solvers.h"


//delta_kernel made for 2+1 system build with delta_M2k_ERE
void test_delta_kernel_2plus1_system_ERE(   Eigen::MatrixXcd &kern_mat, 
                                            comp En, 
                                            std::vector<comp> &pvec_for_m1m2,
                                            std::vector<comp> &kvec_for_m1m1, 
                                            std::vector<comp> &weights_for_pvec_for_m1m2, 
                                            std::vector<comp> &weights_for_kvec_for_m1m1,
                                            comp sigb1,
                                            comp sigb2, 
                                            double m1, 
                                            double m2, 
                                            double eps_for_m2k,
                                            double eps_for_ope, 
                                            double eps_for_cutoff, 
                                            comp total_P,
                                            double a0_m1, 
                                            double r0_m1, 
                                            double eta_1, 
                                            double a0_m2, 
                                            double r0_m2, 
                                            double eta_2,
                                            char debug  )
{
    comp ii = {0.0, 1.0}; 
    comp pi = std::acos(-1.0); 

    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size(); 
    int total_size = size1 + size2; 

    Eigen::MatrixXcd del_M2k_1(size1, size1); 
    Eigen::MatrixXcd del_M2k_2(size2, size2); 
    del_M2k_1 = Eigen::MatrixXcd::Zero(size1, size1); 
    del_M2k_2 = Eigen::MatrixXcd::Zero(size2, size2); 

    Eigen::MatrixXcd Filler0_12(size1, size2); 
    Eigen::MatrixXcd Filler0_21(size2, size1); 
    Eigen::MatrixXcd Filler0_22(size2, size2); 
    Filler0_12 = Eigen::MatrixXcd::Zero(size1, size2); 
    Filler0_21 = Eigen::MatrixXcd::Zero(size2, size1); 
    Filler0_22 = Eigen::MatrixXcd::Zero(size2, size2); 

    //delta_M2k matrix build
    //when m1 is spectator 
    double mi = m1; 
    double mj = m1;
    double mk = m2; 

    for(int i=0; i<size1; ++i)
    {
        comp mom_p = pvec_for_m1m2[i]; 
        
        comp delta_M2k = delta_M2k_ERE(eta_1, En, mom_p, total_P, sigb1, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);
    
        del_M2k_1(i,i) = delta_M2k; 
    }
    //when m2 is spectator 
    mi = m2; 
    mj = m1; 
    mk = m1; 

    for(int i=0; i<size2; ++i)
    {
        comp mom_k = kvec_for_m1m1[i]; 
        
        comp delta_M2k = delta_M2k_ERE(eta_2, En, mom_k, total_P, sigb2, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);
    
        del_M2k_2(i,i) = delta_M2k; 
    }

    Eigen::MatrixXcd delta_M2k_mat(total_size, total_size); 

    delta_M2k_mat << del_M2k_1, Filler0_12, 
                     Filler0_21, 0.5*del_M2k_2; 

    if(debug=='y')
    {
        std::cout<<"===========M2MAT============="<<std::endl; 
        std::cout<<delta_M2k_mat<<std::endl; 
        std::cout<<"============================="<<std::endl; 

    }

    //Gmat build
    //for (i,j) = 1 1 and k = 2
    mi = m1;
    mj = m1; 
    mk = m2; 
    Eigen::MatrixXcd WG_11(size1, size1);

    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size1; ++j)
        {
            comp p = pvec_for_m1m2[i];
            comp k = pvec_for_m1m2[j];
            comp weights = weights_for_pvec_for_m1m2[j];
            comp omgk = omega_func(k,mj);

            comp G = GS_pk(En, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
            WG_11(i,j) = W*G; 
        }
    }

    //for (i,j) = 1 2 and k = 1
    mi = m1; 
    mj = m2; 
    mk = m1; 
    Eigen::MatrixXcd WG_12(size1, size2); 

    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size2; ++j)
        {
            comp p = pvec_for_m1m2[i];
            comp k = kvec_for_m1m1[j]; 
            comp weights = weights_for_kvec_for_m1m1[j];
            comp omgk = omega_func(k,mj);


            comp G = GS_pk(En, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
            WG_12(i,j) = std::sqrt(2.0)*W*G; 
        }
    }

    //for (i,j) = 2 1 and k = 1
    mi = m2; 
    mj = m1; 
    mk = m1; 
    Eigen::MatrixXcd WG_21(size2, size1); 

    for(int i=0; i<size2; ++i)
    {
        for(int j=0; j<size1; ++j)
        {
            comp p = kvec_for_m1m1[i];
            comp k = pvec_for_m1m2[j]; 
            comp weights = weights_for_pvec_for_m1m2[j];
            comp omgk = omega_func(k,mj);

            comp G = GS_pk(En, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
            WG_21(i,j) = std::sqrt(2.0)*W*G; 
        }
    }


    //WG_12 = std::sqrt(2.0)*WG_12; 
    //WG_21 = std::sqrt(2.0)*WG_21; 

    Eigen::MatrixXcd WG_mat(size1 + size2, size1 + size2); 
    //std::cout<<__PRETTY_FUNCTION__<<std::endl; 

    WG_mat <<   WG_11, WG_12,
                WG_21, Filler0_22; 

    if(debug=='y')
    {
        std::cout<<"===========GMAT============="<<std::endl; 
        std::cout<<WG_mat<<std::endl; 
        std::cout<<"============================"<<std::endl; 
    }

    kern_mat = WG_mat*delta_M2k_mat; 

}

//Delta Bmat for 2+1 system build with delta_M2k_ERE
void test_delta_Bmat_2plus1_system_ERE( Eigen::MatrixXcd &B_mat, 
                                        comp En, 
                                        std::vector<comp> &pvec_for_m1m2,
                                        std::vector<comp> &kvec_for_m1m1, 
                                        std::vector<comp> &weights_for_pvec_for_m1m2,
                                        std::vector<comp> &weights_for_kvec_for_m1m1,
                                        comp sigb1, 
                                        comp sigb2, 
                                        double m1, 
                                        double m2, 
                                        double eps_for_m2k, 
                                        double eps_for_ope, 
                                        double eps_for_cutoff, 
                                        comp total_P, 
                                        double a0_m1, 
                                        double r0_m1, 
                                        double eta_1, 
                                        double a0_m2, 
                                        double r0_m2,
                                        double eta_2,
                                        char debug    )
{
    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size(); 
    int total_size = size1 + size2; 

    Eigen::MatrixXcd delta_kern_mat(total_size, total_size); 
    Eigen::MatrixXcd Imat(total_size, total_size); 
    Imat = Eigen::MatrixXcd::Identity(total_size, total_size); 

    test_delta_kernel_2plus1_system_ERE(delta_kern_mat, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, sigb1, sigb2, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2, debug);

    B_mat = Imat + delta_kern_mat; 
}

void negative_GSpqbi_mat_for_Kphib_ERE( Eigen::MatrixXcd &negGpqbi_mat,
                                        comp En, 
                                        std::vector<comp> &pvec_for_m1m2, 
                                        std::vector<comp> &kvec_for_m1m1, 
                                        std::vector<comp> &weights_for_m1m2, 
                                        std::vector<comp> &weights_for_m1m1, 
                                        comp qb_1, 
                                        comp qb_2, 
                                        double m1, 
                                        double m2, 
                                        double eps_for_ope, 
                                        double eps_for_cutoff, 
                                        comp total_P,
                                        char debug     )
{
    comp ii = {0.0, 1.0};
    comp pi = std::acos(-1.0); 
    
    double mi, mj, mk;

    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd Filler0_22(size2, 1); 
    Filler0_22 = Eigen::MatrixXcd::Zero(size2, 1); 

    //for (i,j) = 1 1 and k = 2
    
    mi = m1; 
    mj = m1;
    mk = m2; 

    Eigen::MatrixXcd G_11(size1, 1);

    for(int i=0; i<size1; ++i)
    {
        comp p = pvec_for_m1m2[i];
        comp G = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
        G_11(i,0) = G;
    } 

    //for (i,j) = 1 2 and k = 1
    
    mi = m1; 
    mj = m2;
    mk = m1; 

    Eigen::MatrixXcd G_12(size1, 1);

    for(int i=0; i<size1; ++i)
    {
        comp p = pvec_for_m1m2[i];
        comp G = GS_pk(En, p, qb_2, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
        G_12(i,0) = G;
    } 
    //for (i,j) = 2 1 and k = 1
    
    mi = m2; 
    mj = m1;
    mk = m1; 

    Eigen::MatrixXcd G_21(size2, 1);

    for(int i=0; i<size2; ++i)
    {
        comp p = kvec_for_m1m1[i];
        comp G = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
        G_21(i,0) = G;
    } 

    G_12 = std::sqrt(2.0)*G_12; 
    G_21 = std::sqrt(2.0)*G_21; 

    Eigen::MatrixXcd Gmat(size1 + size2, 2); 

    Gmat << G_11, G_12, 
            G_21, Filler0_22;

    negGpqbi_mat = -1.0*Gmat; 

    if(debug=='y')
    {
        std::cout<<"======== negative Gmat(p,qb_i) ========="<<std::endl; 
        std::cout<<negGpqbi_mat<<std::endl; 
        std::cout<<"======================================"<<std::endl; 
    }

}

//This is Kphib_pqbi, the definition is a little different
//than setting the outgoing spectator momentum to qb_1 or qb_2
//outgoing spectator momentum qb_i is based on what flavor we 
//have as spectator. 
void Kphib_pqbi_2plus1_system_ERE(  Eigen::MatrixXcd &Kphib_pqbi_mat, 
                                    comp En, 
                                    std::vector<comp> &pvec_for_m1m2, 
                                    std::vector<comp> &kvec_for_m1m1, 
                                    std::vector<comp> &weights_for_pvec_for_m1m2, 
                                    std::vector<comp> &weights_for_kvec_for_m1m1, 
                                    comp sigb1, 
                                    comp sigb2,
                                    double m1, 
                                    double m2, 
                                    double eps_for_m2k, 
                                    double eps_for_ope, 
                                    double eps_for_cutoff, 
                                    comp total_P, 
                                    double a0_m1, 
                                    double r0_m1, 
                                    double eta_1, 
                                    double a0_m2, 
                                    double r0_m2, 
                                    double eta_2, 
                                    comp qb_1, 
                                    comp qb_2,
                                    char debug  )
{
    comp ii = {0.0, 1.0}; 
    comp pi = std::acos(-1.0); 

    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size();
    int total_size = size1 + size2; 

    Eigen::MatrixXcd delta_Bmat; 
    Eigen::MatrixXcd negGpqbi_mat; 
    Eigen::MatrixXcd Kphib_mat; 

    //comp qb_1 = qb_i(En, sigb_1, m1); 
    //comp qb_2 = qb_i(En, sigb_2, m2); 

    test_delta_Bmat_2plus1_system_ERE(delta_Bmat, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, sigb1, sigb2, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2, debug); 
    negative_GSpqbi_mat_for_Kphib_ERE(negGpqbi_mat, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, qb_1, qb_2, m1, m2, eps_for_ope, eps_for_cutoff, total_P, debug); 

    double relerr; 
    LinearSolver_2(delta_Bmat, Kphib_mat, negGpqbi_mat, relerr); 

    //GPU Solver 
    int mat_length = pvec_for_m1m2.size() + kvec_for_m1m1.size(); 
    int mat_width = static_cast<int>(negGpqbi_mat.cols()); 
    //d_mat = Eigen::MatrixXcd(mat_length, mat_width);

    //std::cout<<"here"<<std::endl; 
    //std::cout<<mat_width<<std::endl; 
    //cusolverComplex_mat(B_mat, negG_mat, d_mat, mat_length, mat_width); 

    Kphib_pqbi_mat = Kphib_mat; 
}

//We do the Kphib interpolator in 3 ways
//Kphib(qb_i, qb_i) which modifies the incoming 
//and outgoing spectator binding momentum based on 
//the spectator itself, the other two ways should be 
//based on the arbitrary spectator momentum.

//We first build Kphib(qb_i, qb_i) interpolator
//If we set both of the qb_i to qb_1 or qb_2, this 
//will give me interpolated results for K(q,qb_i)
//for arbitrary q = qb_1 or qb_2 

void test_Kphib_qbi_qbi_interpolator(   Eigen::MatrixXcd &result_Kphib_qq, 
                                        Eigen::MatrixXcd &Kphib_pqbi, //this is input Kphib_p_qbi
                                        comp En, 
                                        std::vector<comp> &pvec_for_m1m2, 
                                        std::vector<comp> &kvec_for_m1m1, 
                                        std::vector<comp> &weights_for_pvec_for_m1m2, 
                                        std::vector<comp> &weights_for_kvec_for_m1m1, 
                                        comp sigb1, 
                                        comp sigb2, 
                                        double m1, 
                                        double m2, 
                                        double eps_for_m2k, 
                                        double eps_for_ope, 
                                        double eps_for_cutoff, 
                                        comp total_P, 
                                        double a0_m1, 
                                        double r0_m1, 
                                        double eta_1, 
                                        double a0_m2, 
                                        double r0_m2, 
                                        double eta_2, 
                                        comp qb_1, 
                                        comp qb_2,
                                        char debug  ) 
{
    double mi, mj, mk;
    comp ii = {0.0, 1.0}; 
    comp pi = std::acos(-1.0); 

    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size(); 
    int total_size = size1 + size2; 

    //G(qb_i, qb_i) build 
    
    Eigen::MatrixXcd Gqbiqbi_mat(2, 2); 

    //(i,j) = 1 1 and k = 2
    mi = m1; 
    mj = m1; 
    mk = m2; 

    comp G11 = GS_pk(En, qb_1, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff);

    //(i,j) = 1 2 and k = 1
    mi = m1; 
    mj = m2;
    mk = m1; 

    comp G12 = GS_pk(En, qb_1, qb_2, mi, mj, mk, eps_for_ope, eps_for_cutoff);

    //(i,j) = 2 1 and k = 1
    mi = m2; 
    mj = m1;
    mk = m1; 

    comp G21 = GS_pk(En, qb_2, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff);

    Gqbiqbi_mat(0,0) = G11; 
    Gqbiqbi_mat(0,1) = std::sqrt(2.0)*G12; 
    Gqbiqbi_mat(1,0) = std::sqrt(2.0)*G21; 
    Gqbiqbi_mat(1,1) = 0.0; 

    //kernel build

    Eigen::MatrixXcd delta_M2k_1(size1, size1); 
    Eigen::MatrixXcd delta_M2k_2(size2, size2); 

    delta_M2k_1 = Eigen::MatrixXcd::Zero(size1, size1); 
    delta_M2k_2 = Eigen::MatrixXcd::Zero(size2, size2); 

    Eigen::MatrixXcd M2k_filler0_12(size1, size2); 
    Eigen::MatrixXcd M2k_filler0_21(size2, size1); 

    M2k_filler0_12 = Eigen::MatrixXcd::Zero(size1, size2); 
    M2k_filler0_21 = Eigen::MatrixXcd::Zero(size2, size1); 

    //spectator = m1
    mi = m1; 
    mj = m1; 
    mk = m2; 

    for(int i=0; i<size1; ++i)
    {
        comp k = pvec_for_m1m2[i]; 
        comp sigb1 = sigma(En, qb_1, mi, total_P);
        comp delta_M2k = delta_M2k_ERE(eta_1, En, k, total_P, sigb1, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);
        delta_M2k_1(i,i) = delta_M2k; 
    }

    //spectator = m2
    mi = m2; 
    mj = m1; 
    mk = m1; 

    for(int i=0; i<size2; ++i)
    {
        comp k = kvec_for_m1m1[i]; 
        comp sigb2 = sigma(En, qb_2, mi, total_P); 
        comp delta_M2k = delta_M2k_ERE(eta_2, En, k, total_P, sigb2, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);
        delta_M2k_2(i,i) = delta_M2k; 
    }

    Eigen::MatrixXcd delta_M2k_mat(total_size, total_size); 
    delta_M2k_mat << delta_M2k_1, M2k_filler0_12,
                     M2k_filler0_21, 0.5*delta_M2k_2; 
    
    
    //G(qb_i, k) build 

    //(i,j) = 1 1 and k = 2
    mi = m1; 
    mj = m1; 
    mk = m2; 
    Eigen::MatrixXcd WG_11(1, size1); 

    for(int i=0; i<size1; ++i)
    {
        comp k = pvec_for_m1m2[i]; 
        comp weights = weights_for_pvec_for_m1m2[i]; 
        comp omgk = omega_func(k, mj); 

        comp G = GS_pk(En, qb_1, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
        comp W = weights*k*k/(pow(2.0*pi, 2.0)*omgk);
        WG_11(0,i) = W*G; 
    }

    //(i,j) = 1 2 and k = 1
    mi = m1; 
    mj = m2; 
    mk = m1; 
    Eigen::MatrixXcd WG_12(1, size2); 

    for(int i=0; i<size2; ++i)
    {
        comp k = kvec_for_m1m1[i]; 
        comp weights = weights_for_kvec_for_m1m1[i]; 
        comp omgk = omega_func(k, mj); 

        comp G = GS_pk(En, qb_1, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
        comp W = weights*k*k/(pow(2.0*pi, 2.0)*omgk);
        WG_12(0,i) = std::sqrt(2.0)*W*G; 
    }

    //(i,j) = 2 1 and k = 1
    mi = m2; 
    mj = m1; 
    mk = m1; 
    Eigen::MatrixXcd WG_21(1, size1); 

    for(int i=0; i<size1; ++i)
    {
        comp k = pvec_for_m1m2[i]; 
        comp weights = weights_for_pvec_for_m1m2[i]; 
        comp omgk = omega_func(k, mj); 

        comp G = GS_pk(En, qb_2, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
        comp W = weights*k*k/(pow(2.0*pi, 2.0)*omgk);
        WG_21(0,i) = std::sqrt(2.0)*W*G; 
    }

    Eigen::MatrixXcd WG_mat(2, total_size); 
    Eigen::MatrixXcd WG_filler0_22(1, size2); 
    WG_filler0_22 = Eigen::MatrixXcd::Zero(1, size2); 

    WG_mat <<   WG_11, WG_12, 
                WG_21, WG_filler0_22; 

    result_Kphib_qq = - Gqbiqbi_mat - WG_mat*delta_M2k_mat*Kphib_pqbi;

    
}

//Rho_phib matrix is a matrix that has
//the kinematic phase space and two-body
//bound-state pole residue 

void Rho_phib_mat(  Eigen::MatrixXcd &Rho_mat, 
                    comp En, 
                    double eta_1, 
                    double eta_2, 
                    double m1, 
                    double m2, 
                    comp qb_1, 
                    comp qb_2,
                    comp total_P   )
{
    comp ii = {0.0, 1.0}; 

    comp sigb1 = sigma(En, qb_1, m1, total_P); 
    comp sigb2 = sigma(En, qb_2, m2, total_P); 

    comp g1 = gfunc_i(eta_1, sigb1, m1, m2);
    comp g2 = gfunc_i(eta_2, sigb2, m1, m1); 

    comp rho1 = rhophib(qb_1, En); 
    comp rho2 = rhophib(qb_2, En); 

    Eigen::MatrixXcd R_mat(2, 2); 

    R_mat(0,0) = ii*g1*g1*rho1; 
    R_mat(0,1) = 0.0; 
    R_mat(1,0) = 0.0;
    R_mat(1,1) = 0.5*ii*g2*g2*rho2; 

    Rho_mat = R_mat; 
}

//dS(qb_i, qb_i) 

void test_dS_qbiqbi_SA_Method_2plus1_system_ERE(    Eigen::MatrixXcd &dqbiqbi_mat, 
                                                    comp En,
                                                    double m1, 
                                                    double m2,
                                                    std::vector<comp> &pvec_for_m1m2, 
                                                    std::vector<comp> &weights_for_pvec_for_m1m2,
                                                    std::vector<comp> &kvec_for_m1m1, 
                                                    std::vector<comp> &weights_for_kvec_for_m1m1, 
                                                    comp sigb1, 
                                                    comp sigb2, 
                                                    comp qb_1, 
                                                    comp qb_2, 
                                                    double eps_for_m2k, 
                                                    double eps_for_ope, 
                                                    double eps_for_cutoff, 
                                                    comp total_P, 
                                                    double a0_m1, 
                                                    double r0_m1, 
                                                    double eta_1, 
                                                    double a0_m2, 
                                                    double r0_m2, 
                                                    double eta_2, 
                                                    char debug      )
{
    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd Kphib_pqbi_mat; 
    Eigen::MatrixXcd Kphib_qbiqbi_mat; 
    Kphib_pqbi_2plus1_system_ERE(Kphib_pqbi_mat, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, sigb1, sigb2, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2, qb_1, qb_2, debug); 
    test_Kphib_qbi_qbi_interpolator(Kphib_qbiqbi_mat, Kphib_pqbi_mat, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, sigb1, sigb2, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2, qb_1, qb_2, debug); 
    Eigen::MatrixXcd I_mat;
    I_mat = Eigen::MatrixXcd::Identity(2, 2); 

    Eigen::MatrixXcd Rho_mat; 
    Rho_phib_mat(Rho_mat, En, eta_1, eta_2, m1, m2, qb_1, qb_2, total_P); 

    Eigen::MatrixXcd num = Kphib_qbiqbi_mat; 
    Eigen::MatrixXcd denom = I_mat - Kphib_qbiqbi_mat*Rho_mat; 

    Eigen::MatrixXcd denom_inv = denom.inverse(); 

    Eigen::MatrixXcd res = denom_inv*num; 

    dqbiqbi_mat = res; 
}

#endif 