#ifndef INTEGRALEQN_MOM_BASED_H
#define INTEGRALEQN_MOM_BASED_H

#include "functions_mom_based.h"
#include<Eigen/Dense> 


/* Here the convention is the 2+1 system is defined
with 2 particles of mass m1 and one particle of mass
m2, */
/* For the analysis of KKpi system, mK is m1 
and mpi is m2 */
/* This function is the matrix form of 
weights*k^2/{(2pi)^2 omega_k} G M2 */
/* We call the weights*integral_measure*G matrix as WG_matrix */
void kernel_2plus1_system_ERE(  Eigen::MatrixXcd &kern_mat,
                                comp E, 
                                std::vector<comp> &pvec_for_m1m2,
                                std::vector<comp> &kvec_for_m1m1, 
                                std::vector<comp> &weights_for_pvec_for_m1m2,
                                std::vector<comp> &weights_for_kvec_for_m1m1,
                                double m1, 
                                double m2,
                                double eps_for_m2k,
                                double eps_for_ope,
                                double eps_for_cutoff,
                                comp total_P,
                                double a0_m1,
                                double r0_m1, 
                                double a0_m2,
                                double r0_m2    )
{
    //Boiler Plate
    char debug = 'n';
    comp ii = {0.0,1.0};
    comp pi = std::acos(-1.0); 
    //double r = 0.0;//setting effective range to 0 for both M2k

    int size_for_M2_1 = pvec_for_m1m2.size(); 
    int size_for_M2_2 = kvec_for_m1m1.size(); 

    //when m1 and m2 are in twobody subchannel
    Eigen::MatrixXcd M2k_for_m1m2(size_for_M2_1, size_for_M2_1);
    //when m1 and m1 are in twobody subchannel 
    Eigen::MatrixXcd M2k_for_m1m1(size_for_M2_2, size_for_M2_2);
    M2k_for_m1m2 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_1); 
    M2k_for_m1m1 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2); 
    //Now M2k for m1 and m2 goes in the 11 position of the M2 matrix
    //and M2k for m1 and m1 goes in the 22 position of the M2 matrix
    //The rest shall be filled with zeros 

    Eigen::MatrixXcd Filler0_12(size_for_M2_1, size_for_M2_2);
    Eigen::MatrixXcd Filler0_21(size_for_M2_2, size_for_M2_1);
    Eigen::MatrixXcd Filler0_22(size_for_M2_2, size_for_M2_2);

    Filler0_12 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_2);
    Filler0_21 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_1);
    Filler0_22 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2);


    //First we setup the matrix for M2
    //if m1 is spectator then the 2body sub-channel
    //have m1 and m2 
    
    double mi = m1;
    double mj = m1;
    double mk = m2; 
    
    

    for(int i=0; i<size_for_M2_1; ++i)
    {
        comp mom_p = pvec_for_m1m2[i];

        comp M2k = M2k_tilde_ERE(E, mom_p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);

        M2k_for_m1m2(i,i) = M2k; 
    }

    //Now m2 is the spectator, so the two body subchannel 
    //have m1 and m1 
    
    mi = m2;
    mj = m1; 
    mk = m1; 
    
    for(int i=0; i<size_for_M2_2; ++i)
    {
        comp mom_k = kvec_for_m1m1[i];

        comp M2k = M2k_tilde_ERE(E, mom_k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);

        M2k_for_m1m1(i,i) = M2k; 
    }

    //M2k_for_m1m1 = 0.5*M2k_for_m1m1; 

    Eigen::MatrixXcd M2k_mat(size_for_M2_1 + size_for_M2_2, size_for_M2_1 + size_for_M2_2);

    M2k_mat <<  M2k_for_m1m2, Filler0_12,
                Filler0_21, M2k_for_m1m1;

    if(debug=='y')
    {
        std::cout<<"===========M2MAT============="<<std::endl; 
        std::cout<<M2k_mat<<std::endl; 
        std::cout<<"============================="<<std::endl; 

    }


    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size();

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

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
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


            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
            WG_12(i,j) = W*G; 
        }
    }

    //for (i,j) = 2 1 and k = 1
    mi = m2; 
    mj = m1; 
    mk = m1; 
    Eigen::MatrixXcd WG_21(size2, size1); 

    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size2; ++j)
        {
            comp p = kvec_for_m1m1[i];
            comp k = pvec_for_m1m2[j]; 
            comp weights = weights_for_pvec_for_m1m2[j];
            comp omgk = omega_func(k,mj);

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
            WG_21(i,j) = W*G; 
        }
    }

    WG_12 = std::sqrt(2.0)*WG_12; 
    WG_21 = std::sqrt(2.0)*WG_21; 

    Eigen::MatrixXcd WG_mat(size1 + size2, size1 + size2); 

    WG_mat <<   WG_11, WG_12,
                WG_21, Filler0_22; 

    if(debug=='y')
    {
        std::cout<<"===========GMAT============="<<std::endl; 
        std::cout<<WG_mat<<std::endl; 
        std::cout<<"============================"<<std::endl; 
    }

    kern_mat = WG_mat*M2k_mat; 

}

//This code works but it has the twobody
//symmetry factor hardcoded in it, we are no
//longer using this code, use the next iteration
//from now on. 
void test_kernel_2plus1_system_ERE(  Eigen::MatrixXcd &kern_mat,
                                comp E, 
                                std::vector<comp> &pvec_for_m1m2,
                                std::vector<comp> &kvec_for_m1m1, 
                                std::vector<comp> &weights_for_pvec_for_m1m2,
                                std::vector<comp> &weights_for_kvec_for_m1m1,
                                double m1, 
                                double m2,
                                double eps_for_m2k,
                                double eps_for_ope,
                                double eps_for_cutoff,
                                comp total_P,
                                double a0_m1,
                                double r0_m1, 
                                double a0_m2,
                                double r0_m2    )
{
    //Boiler Plate
    char debug = 'n';
    comp ii = {0.0,1.0};
    comp pi = std::acos(-1.0); 

    //this is the symmetrization factor
    //we need to comment this later on 
    double eta_1 = 0.5; 
    double eta_2 = 0.5; 
    //double r = 0.0;//setting effective range to 0 for both M2k

    int size_for_M2_1 = pvec_for_m1m2.size(); 
    int size_for_M2_2 = kvec_for_m1m1.size(); 

    //when m1 and m2 are in twobody subchannel
    Eigen::MatrixXcd M2k_for_m1m2(size_for_M2_1, size_for_M2_1);
    //when m1 and m1 are in twobody subchannel 
    Eigen::MatrixXcd M2k_for_m1m1(size_for_M2_2, size_for_M2_2);
    M2k_for_m1m2 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_1); 
    M2k_for_m1m1 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2); 
    //Now M2k for m1 and m2 goes in the 11 position of the M2 matrix
    //and M2k for m1 and m1 goes in the 22 position of the M2 matrix
    //The rest shall be filled with zeros 

    Eigen::MatrixXcd Filler0_12(size_for_M2_1, size_for_M2_2);
    Eigen::MatrixXcd Filler0_21(size_for_M2_2, size_for_M2_1);
    Eigen::MatrixXcd Filler0_22(size_for_M2_2, size_for_M2_2);

    Filler0_12 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_2);
    Filler0_21 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_1);
    Filler0_22 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2);


    //First we setup the matrix for M2
    //if m1 is spectator then the 2body sub-channel
    //have m1 and m2 
    
    double mi = m1;
    double mj = m1;
    double mk = m2; 
    
    

    for(int i=0; i<size_for_M2_1; ++i)
    {
        comp mom_p = pvec_for_m1m2[i];

        //comp M2k = M2k_tilde_ERE(E, mom_p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);
        comp M2k = M2k_ERE(eta_1, E, mom_p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);

        M2k_for_m1m2(i,i) = M2k; 
    }

    //Now m2 is the spectator, so the two body subchannel 
    //have m1 and m1 
    
    mi = m2;
    mj = m1; 
    mk = m1; 
    
    for(int i=0; i<size_for_M2_2; ++i)
    {
        comp mom_k = kvec_for_m1m1[i];

        //comp M2k = M2k_tilde_ERE(E, mom_k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);
        comp M2k = M2k_ERE(eta_2, E, mom_k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);

        M2k_for_m1m1(i,i) = M2k; 
    }

    //M2k_for_m1m1 = 0.5*M2k_for_m1m1; 

    Eigen::MatrixXcd M2k_mat(size_for_M2_1 + size_for_M2_2, size_for_M2_1 + size_for_M2_2);

    /*
    M2k_mat <<  M2k_for_m1m2, Filler0_12,
                Filler0_21, M2k_for_m1m1;
    */

    M2k_mat <<  M2k_for_m1m2, Filler0_12,
                Filler0_21, 0.5*M2k_for_m1m1;
    
    if(debug=='y')
    {
        std::cout<<"===========M2MAT============="<<std::endl; 
        std::cout<<M2k_mat<<std::endl; 
        std::cout<<"============================="<<std::endl; 

    }

    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size();

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

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
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


            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
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

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
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

    kern_mat = WG_mat*M2k_mat; 


}

//the next iteration will be used from now on,
//this has added dependency 
//on the twobody symmetry factor
void test_kernel_2plus1_system_ERE_1(   Eigen::MatrixXcd &kern_mat,
                                        comp E, 
                                        std::vector<comp> &pvec_for_m1m2,
                                        std::vector<comp> &kvec_for_m1m1, 
                                        std::vector<comp> &weights_for_pvec_for_m1m2,
                                        std::vector<comp> &weights_for_kvec_for_m1m1,
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
                                        double eta_2    )
{
    //Boiler Plate
    char debug = 'n';
    comp ii = {0.0,1.0};
    comp pi = std::acos(-1.0); 

    //this is the symmetrization factor
    //we need to comment this later on 
    //double eta_1 = 0.5; 
    //double eta_2 = 0.5; 
    //double r = 0.0;//setting effective range to 0 for both M2k

    int size_for_M2_1 = pvec_for_m1m2.size(); 
    int size_for_M2_2 = kvec_for_m1m1.size(); 

    //when m1 and m2 are in twobody subchannel
    Eigen::MatrixXcd M2k_for_m1m2(size_for_M2_1, size_for_M2_1);
    //when m1 and m1 are in twobody subchannel 
    Eigen::MatrixXcd M2k_for_m1m1(size_for_M2_2, size_for_M2_2);
    M2k_for_m1m2 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_1); 
    M2k_for_m1m1 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2); 
    //Now M2k for m1 and m2 goes in the 11 position of the M2 matrix
    //and M2k for m1 and m1 goes in the 22 position of the M2 matrix
    //The rest shall be filled with zeros 

    Eigen::MatrixXcd Filler0_12(size_for_M2_1, size_for_M2_2);
    Eigen::MatrixXcd Filler0_21(size_for_M2_2, size_for_M2_1);
    Eigen::MatrixXcd Filler0_22(size_for_M2_2, size_for_M2_2);

    Filler0_12 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_2);
    Filler0_21 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_1);
    Filler0_22 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2);


    //First we setup the matrix for M2
    //if m1 is spectator then the 2body sub-channel
    //have m1 and m2 
    
    double mi = m1;
    double mj = m1;
    double mk = m2; 
    
    

    for(int i=0; i<size_for_M2_1; ++i)
    {
        comp mom_p = pvec_for_m1m2[i];

        //comp M2k = M2k_tilde_ERE(E, mom_p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);
        comp M2k = M2k_ERE(eta_1, E, mom_p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);

        M2k_for_m1m2(i,i) = M2k; 
    }

    //Now m2 is the spectator, so the two body subchannel 
    //have m1 and m1 
    
    mi = m2;
    mj = m1; 
    mk = m1; 
    
    for(int i=0; i<size_for_M2_2; ++i)
    {
        comp mom_k = kvec_for_m1m1[i];

        //comp M2k = M2k_tilde_ERE(E, mom_k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);
        comp M2k = M2k_ERE(eta_2, E, mom_k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);

        M2k_for_m1m1(i,i) = M2k; 
    }

    //M2k_for_m1m1 = 0.5*M2k_for_m1m1; 

    Eigen::MatrixXcd M2k_mat(size_for_M2_1 + size_for_M2_2, size_for_M2_1 + size_for_M2_2);

    /*
    M2k_mat <<  M2k_for_m1m2, Filler0_12,
                Filler0_21, M2k_for_m1m1;
    */

    M2k_mat <<  M2k_for_m1m2, Filler0_12,
                Filler0_21, 0.5*M2k_for_m1m1;
    
    if(debug=='y')
    {
        std::cout<<"===========M2MAT============="<<std::endl; 
        std::cout<<M2k_mat<<std::endl; 
        std::cout<<"============================="<<std::endl; 

    }

    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size();

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

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
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


            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
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

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
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

    kern_mat = WG_mat*M2k_mat; 


}

//This has distinct epsilons to be induced for different 
//two body subchannels, ope and cutoff
//the 2 epsilons for M2k are straight forward to use
//for OPE, their usage depends on which M2k is added to 
//the ope and which spectator correspond to i and j val. 
//the two ope epsilons should be same and it should depend 
//on which spectator is latched with j index 
void test_kernel_2plus1_system_ERE_2(   Eigen::MatrixXcd &kern_mat,
                                        comp E, 
                                        std::vector<comp> &pvec_for_m1m2,
                                        std::vector<comp> &kvec_for_m1m1, 
                                        std::vector<comp> &weights_for_pvec_for_m1m2,
                                        std::vector<comp> &weights_for_kvec_for_m1m1,
                                        double m1, 
                                        double m2,
                                        double eps_for_m2k_1,
                                        double eps_for_m2k_2,
                                        double eps_for_ope_1,
                                        double eps_for_ope_2,
                                        double eps_for_cutoff_1,
                                        double eps_for_cutoff_2,
                                        comp total_P,
                                        double a0_m1,
                                        double r0_m1, 
                                        double eta_1, 
                                        double a0_m2,
                                        double r0_m2,
                                        double eta_2    )
{
    //Boiler Plate
    char debug = 'n';
    comp ii = {0.0,1.0};
    comp pi = std::acos(-1.0); 

    //this is the symmetrization factor
    //we need to comment this later on 
    //double eta_1 = 0.5; 
    //double eta_2 = 0.5; 
    //double r = 0.0;//setting effective range to 0 for both M2k

    int size_for_M2_1 = pvec_for_m1m2.size(); 
    int size_for_M2_2 = kvec_for_m1m1.size(); 

    //when m1 and m2 are in twobody subchannel
    Eigen::MatrixXcd M2k_for_m1m2(size_for_M2_1, size_for_M2_1);
    //when m1 and m1 are in twobody subchannel 
    Eigen::MatrixXcd M2k_for_m1m1(size_for_M2_2, size_for_M2_2);
    M2k_for_m1m2 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_1); 
    M2k_for_m1m1 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2); 
    //Now M2k for m1 and m2 goes in the 11 position of the M2 matrix
    //and M2k for m1 and m1 goes in the 22 position of the M2 matrix
    //The rest shall be filled with zeros 

    Eigen::MatrixXcd Filler0_12(size_for_M2_1, size_for_M2_2);
    Eigen::MatrixXcd Filler0_21(size_for_M2_2, size_for_M2_1);
    Eigen::MatrixXcd Filler0_22(size_for_M2_2, size_for_M2_2);

    Filler0_12 = Eigen::MatrixXcd::Zero(size_for_M2_1, size_for_M2_2);
    Filler0_21 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_1);
    Filler0_22 = Eigen::MatrixXcd::Zero(size_for_M2_2, size_for_M2_2);


    //First we setup the matrix for M2
    //if m1 is spectator then the 2body sub-channel
    //have m1 and m2 
    
    double mi = m1;
    double mj = m1;
    double mk = m2; 
    
    

    for(int i=0; i<size_for_M2_1; ++i)
    {
        comp mom_p = pvec_for_m1m2[i];

        //comp M2k = M2k_tilde_ERE(E, mom_p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);
        comp M2k = M2k_ERE(eta_1, E, mom_p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k_1);

        M2k_for_m1m2(i,i) = M2k; 
    }

    //Now m2 is the spectator, so the two body subchannel 
    //have m1 and m1 
    
    mi = m2;
    mj = m1; 
    mk = m1; 
    
    for(int i=0; i<size_for_M2_2; ++i)
    {
        comp mom_k = kvec_for_m1m1[i];

        //comp M2k = M2k_tilde_ERE(E, mom_k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);
        comp M2k = M2k_ERE(eta_2, E, mom_k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k_2);

        M2k_for_m1m1(i,i) = M2k; 
    }

    //M2k_for_m1m1 = 0.5*M2k_for_m1m1; 

    Eigen::MatrixXcd M2k_mat(size_for_M2_1 + size_for_M2_2, size_for_M2_1 + size_for_M2_2);

    /*
    M2k_mat <<  M2k_for_m1m2, Filler0_12,
                Filler0_21, M2k_for_m1m1;
    */

    M2k_mat <<  M2k_for_m1m2, Filler0_12,
                Filler0_21, 0.5*M2k_for_m1m1;
    
    if(debug=='y')
    {
        std::cout<<"===========M2MAT============="<<std::endl; 
        std::cout<<M2k_mat<<std::endl; 
        std::cout<<"============================="<<std::endl; 

    }

    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size();

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

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope_1, eps_for_cutoff_1);
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


            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope_2, eps_for_cutoff_2);
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

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope_1, eps_for_cutoff_1);
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

    kern_mat = WG_mat*M2k_mat; 


}


void Bmat_2plus1_system_ERE(    Eigen::MatrixXcd &B_mat,
                                comp E, 
                                std::vector<comp> &pvec_for_m1m2,
                                std::vector<comp> &kvec_for_m1m1, 
                                std::vector<comp> &weights_for_pvec_for_m1m2,
                                std::vector<comp> &weights_for_kvec_for_m1m1,
                                double m1, 
                                double m2,
                                double eps_for_m2k,
                                double eps_for_ope,
                                double eps_for_cutoff,
                                comp total_P,
                                double a0_m1,
                                double r0_m1,
                                double a0_m2, 
                                double r0_m2    )
{
    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd kern_mat(size1 + size2, size1 + size2); 

    Eigen::MatrixXcd Imat(size1 + size2, size1 + size2); 

    Imat = Eigen::MatrixXcd::Identity(size1 + size2, size1 + size2); 

    kernel_2plus1_system_ERE(kern_mat, E, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, a0_m2, r0_m2); 

    B_mat = Imat + kern_mat; 
}

//This code is redundant, the next iteration will 
//used from now on, which has added eta_i dependency
//(two body symmetry factor)
void test_Bmat_2plus1_system_ERE(    Eigen::MatrixXcd &B_mat,
                                comp E, 
                                std::vector<comp> &pvec_for_m1m2,
                                std::vector<comp> &kvec_for_m1m1, 
                                std::vector<comp> &weights_for_pvec_for_m1m2,
                                std::vector<comp> &weights_for_kvec_for_m1m1,
                                double m1, 
                                double m2,
                                double eps_for_m2k,
                                double eps_for_ope,
                                double eps_for_cutoff,
                                comp total_P,
                                double a0_m1,
                                double r0_m1,
                                double a0_m2, 
                                double r0_m2    )
{
    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd kern_mat(size1 + size2, size1 + size2); 

    Eigen::MatrixXcd Imat(size1 + size2, size1 + size2); 

    Imat = Eigen::MatrixXcd::Identity(size1 + size2, size1 + size2); 

    test_kernel_2plus1_system_ERE(kern_mat, E, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, a0_m2, r0_m2); 
    //std::cout<<"test_Bmat_2plus1_system_ERE"<<std::endl; 
    
    B_mat = Imat + kern_mat; 
}


//The next iteration uses distinct epsilons for the two
//two-body subchannels
//this has eta_i dependency
void test_Bmat_2plus1_system_ERE_1(     Eigen::MatrixXcd &B_mat,
                                        comp E, 
                                        std::vector<comp> &pvec_for_m1m2,
                                        std::vector<comp> &kvec_for_m1m1, 
                                        std::vector<comp> &weights_for_pvec_for_m1m2,
                                        std::vector<comp> &weights_for_kvec_for_m1m1,
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
                                        double eta_2    )
{
    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd kern_mat(size1 + size2, size1 + size2); 

    Eigen::MatrixXcd Imat(size1 + size2, size1 + size2); 

    Imat = Eigen::MatrixXcd::Identity(size1 + size2, size1 + size2); 

    test_kernel_2plus1_system_ERE_1(kern_mat, E, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2); 
    //std::cout<<"test_Bmat_2plus1_system_ERE"<<std::endl; 
    
    B_mat = Imat + kern_mat; 
}

//This code will be used from now on
//this has two distinct epsilons for the two twobody
//subchannels, ope and cutoff 
void test_Bmat_2plus1_system_ERE_2(     Eigen::MatrixXcd &B_mat,
                                        comp E, 
                                        std::vector<comp> &pvec_for_m1m2,
                                        std::vector<comp> &kvec_for_m1m1, 
                                        std::vector<comp> &weights_for_pvec_for_m1m2,
                                        std::vector<comp> &weights_for_kvec_for_m1m1,
                                        double m1, 
                                        double m2,
                                        double eps_for_m2k_1,
                                        double eps_for_m2k_2,
                                        double eps_for_ope_1,
                                        double eps_for_ope_2,
                                        double eps_for_cutoff_1,
                                        double eps_for_cutoff_2,
                                        comp total_P,
                                        double a0_m1,
                                        double r0_m1,
                                        double eta_1, 
                                        double a0_m2, 
                                        double r0_m2,
                                        double eta_2    )
{
    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd kern_mat(size1 + size2, size1 + size2); 

    Eigen::MatrixXcd Imat(size1 + size2, size1 + size2); 

    Imat = Eigen::MatrixXcd::Identity(size1 + size2, size1 + size2); 

    test_kernel_2plus1_system_ERE_2(kern_mat, E, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k_1, eps_for_m2k_2, eps_for_ope_1, eps_for_ope_2, eps_for_cutoff_1, eps_for_cutoff_2, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2); 
    //std::cout<<"test_Bmat_2plus1_system_ERE"<<std::endl; 
    
    B_mat = Imat + kern_mat; 
}


//Here we matrix multiply the Gs(p,k')M2(k') to build 
//the B matrix, this is to check the numbers for the 
//actual B matrix used in the calculation coded above 
//This only depends on the scattering length a0_i 
void test_Bmat_GMtilde_multiplied_2plus1_system_ERE(    Eigen::MatrixXcd &B_mat,
                                                        comp En, 
                                                        std::vector<comp> &pvec_for_m1m2,
                                                        std::vector<comp> &kvec_for_m1m1, 
                                                        std::vector<comp> &weights_for_pvec_for_m1m2,
                                                        std::vector<comp> &weights_for_kvec_for_m1m1,
                                                        double m1, 
                                                        double m2,
                                                        double eps_for_m2k,
                                                        double eps_for_ope,
                                                        double eps_for_cutoff,
                                                        comp total_P,
                                                        double a0_m1,
                                                        double a0_m2    )
{
    comp pi = std::acos(-1.0); 
    comp ii = {0.0, 1.0}; 

    double mi = m1; 
    double mj = m1; 
    double mk = m2; 

    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd WG11M1(size1, size1); 
    Eigen::MatrixXcd WG12M2(size1, size2);
    Eigen::MatrixXcd WG21M1(size2, size1); 
    Eigen::MatrixXcd Filler22(size2, size2); 
    Filler22 = Eigen::MatrixXcd::Zero(size2, size2); 

    //for G11 M1 
    mi = m1; 
    mj = m1; 
    mk = m2; 

    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size1; ++j)
        {
            comp mom_p = pvec_for_m1m2[i]; 
            comp mom_k = pvec_for_m1m2[j]; 
            comp weights = weights_for_pvec_for_m1m2[j]; 
            comp omgk = omega_func(mom_k, mj); 
            comp W = weights*mom_k*mom_k/(std::pow(2.0*pi,2.0)*omgk);
            comp G = GS_pk(En, mom_p, mom_k, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
            comp M = M2k_tilde_ERE(En, mom_k, total_P, a0_m1, 0.0, mi, mj, mk, eps_for_m2k);

            WG11M1(i,j) = W*G*M; 
        }
    } 

    //for G12 M2 
    mi = m1; 
    mj = m2; 
    mk = m1; 

    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size2; ++j)
        {
            comp mom_p = pvec_for_m1m2[i]; 
            comp mom_k = kvec_for_m1m1[j]; 
            comp weights = weights_for_kvec_for_m1m1[j]; 
            comp omgk = omega_func(mom_k, mj); 
            comp W = weights*mom_k*mom_k/(std::pow(2.0*pi,2.0)*omgk);
            comp G = GS_pk(En, mom_p, mom_k, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
            comp M = M2k_tilde_ERE(En, mom_k, total_P, a0_m2, 0.0, mj, mk, mi, eps_for_m2k);

            WG12M2(i,j) = std::sqrt(2)*W*G*M; 
        }
    }

    //for G21 M1 
    mi = m2; 
    mj = m1; 
    mk = m1; 

    for(int i=0; i<size2; ++i)
    {
        for(int j=0; j<size1; ++j)
        {
            comp mom_p = kvec_for_m1m1[i]; 
            comp mom_k = pvec_for_m1m2[j]; 
            comp weights = weights_for_pvec_for_m1m2[j]; 
            comp omgk = omega_func(mom_k, mj); 
            comp W = weights*mom_k*mom_k/(std::pow(2.0*pi,2.0)*omgk);
            comp G = GS_pk(En, mom_p, mom_k, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
            comp M = M2k_tilde_ERE(En, mom_k, total_P, a0_m1, 0.0, mj, mk, mi, eps_for_m2k);

            WG21M1(i,j) = std::sqrt(2)*W*G*M; 
        }
    }

    Eigen::MatrixXcd Imat(size1 + size2, size1 + size2); 
    Imat = Eigen::MatrixXcd::Identity(size1 + size2, size1 + size2); 

    Eigen::MatrixXcd WGM(size1 + size2, size1 + size2); 
    WGM << WG11M1, WG12M2, 
           WG21M1, Filler22; 
    
    B_mat = Imat + WGM; 

}


//This part is mainly used to build the inhomogeneous part of 
//the integral equation B(p,k)d(p,k) = -G(p,k), this negative 
//G matrix is built here. 
void negative_Gmat_2plus1_system(   Eigen::MatrixXcd &neg_G_mat,
                                    comp E, 
                                    std::vector<comp> &pvec_for_m1m2, 
                                    std::vector<comp> &kvec_for_m1m1,
                                    std::vector<comp> &weights_for_pvec_for_m1m2, 
                                    std::vector<comp> &weights_for_kvec_for_m1m1, 
                                    double m1, 
                                    double m2, 
                                    double eps_for_ope,
                                    double eps_for_cutoff,
                                    comp total_P    )
{
    char debug = 'n';
    double mi,mj,mk; 
    comp ii = {0.0,1.0};
    comp pi = std::acos(-1.0); 

    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd Filler0_22(size2, size2);
    Filler0_22 = Eigen::MatrixXcd::Zero(size2, size2);

    //for (i,j) = 1 1 and k = 2
    
    mi = m1; 
    mj = m1;
    mk = m2; 

    Eigen::MatrixXcd G_11(size1, size1); 

    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size1; ++j)
        {
            comp p = pvec_for_m1m2[i];
            comp k = pvec_for_m1m2[j]; 

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
            G_11(i,j) = G; 
        }
    }

    //for (i,j) = 1 2 and k = 1
    mi = m1; 
    mj = m2; 
    mk = m1; 

    Eigen::MatrixXcd G_12(size1, size2); 

    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size2; ++j)
        {
            comp p = pvec_for_m1m2[i];
            comp k = kvec_for_m1m1[j]; 

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
            G_12(i,j) = G; 
        }
    }

    //for (i,j) = 2 1 and k = 1
    mi = m2; 
    mj = m1; 
    mk = m1; 

    Eigen::MatrixXcd G_21(size2, size1); 

    for(int i=0; i<size2; ++i)
    {
        for(int j=0; j<size1; ++j)
        {
            comp p = kvec_for_m1m1[i];
            comp k = pvec_for_m1m2[j]; 

            comp G = GS_pk(E, p, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
            G_21(i,j) = G; 
        }
    }

    G_12 = std::sqrt(2.0)*G_12; 
    G_21 = std::sqrt(2.0)*G_21; 

    Eigen::MatrixXcd Gmat(size1 + size2, size1 + size2); 

    Gmat << G_11, G_12, 
            G_21, Filler0_22;

    neg_G_mat = -Gmat; 

    if(debug=='y')
    {
        std::cout<<"======== negative Gmat ========="<<std::endl; 
        std::cout<<neg_G_mat<<std::endl; 
        std::cout<<"================================"<<std::endl; 
    }


}

/*===========================================================*/

                //M3df Integral Equation

/*===========================================================*/


//The following is based on the Ls rescattering matrix 
void LSmat_2plus1_system(   Eigen::MatrixXcd &Lsmat
                            )
{
    
}



#endif 