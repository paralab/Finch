//
// Created by milinda on 1/19/17.
//

/**
 *
 * @author Milinda Fernando
 * School of Computing University of Utah.
 * @brief Constains lapack routines such as linear system solve, eigen solve to build the interpolation matrices.
 *
 *
 * */


#ifndef SFCSORTBENCH_LAPAC_H
#define SFCSORTBENCH_LAPAC_H

#include <cstring>
#include <iostream>
//#include "lapacke.h"

extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,double* b, int* ldb, int* info );
extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,double* w, double* work, int* lwork, int* info);

namespace lapack
{

/**
 *  @brief: Wrapper for LAPACK DGESV solver for AX=B. Parameters are given below.
 *  @param[in] n : number of rows or columns of linear system
 *  @param[in] nrhs: number of right hand sides.
 *  @param[in] A: matrix A (row major)
 *  @param[in] lda: leading dimention of the array A
 *  @param[in] B: matrix B (row major)
 *  @param[out] X:  matrix X (solution)
 *  @param[in]  ldb:  leading dimention of B
 *  @param[out] info:  returns the status of the solve.
 */
inline void lapack_DGESV(int n , int nrhs, const double * A, int lda, double * B, double * X, int ldb, int info)
    {
        int * ipiv = new int [n] ;
        //memcpy(X,B,sizeof(double)*n*nrhs);

        for(unsigned int i=0;i<nrhs;i++)
            for(unsigned int j=0;j<n;j++)
                X[i*n+j]=B[j*nrhs+i];

        double * L=new double[n*n];

        for(unsigned int i=0;i<n;i++)
            for(unsigned int j=0;j<n;j++)
                L[i*n+j] =A [j*n+i];

        dgesv_(&n,&nrhs,L,&lda,ipiv,X,&lda,&info);

        if( info > 0 ) {
            printf( "The diagonal element of the triangular factor of A,\n" );
            printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
            printf( "the solution could not be computed.\n" );

        }

        memcpy(L,X,sizeof(double)*n*nrhs);

        for(unsigned int i=0;i<nrhs;i++)
            for(unsigned int j=0;j<n;j++)
                X[i*nrhs+j] = L[j*n+i];


        /*lapack_int * ipiv = (lapack_int *)malloc(n*sizeof(lapack_int)) ;
        memcpy(X,B,sizeof(double)*n*nrhs);

        double * L=new double[n*n];
        memcpy(L,A,sizeof(double)*n*n);

        info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, L, lda, ipiv, X,ldb);*/


        if (info <0) {printf(" lapack linear solve failed. \n"); }
        delete [] L;
        delete [] ipiv;

        return ;
    }

/**
 *  @brief: Wrapper for LAPACK DGESV compute eigen values of a square matrix of A. Parameters are given below.
 *  @param[in] n : number of rows or columns of linear system
 *  @param[in] A: matrix A (row major)
 *  @param[in] lda: leading dimention of the array A
 *  @param[out] wr: real part of eigen values
 *  @param[out] vs eigen vectors (row major)
 *  @param[out] info:  returns the status of the solve.
 */


inline void lapack_DSYEV(int n, const double * A, int lda, double * wr,double * vs,int info )
{


    double * L=new double[n*n];

    for(unsigned int i=0;i<n;i++)
        for(unsigned int j=0;j<n;j++)
            L[i*n+j] =A [j*n+i];

    double wkopt;
    double* work;
    int lwork=-1;

    char Vectors[] = "Vectors";
    char Upper[] = "Upper";

    dsyev_( Vectors, Upper, (int*)&n, L, (int*)&lda, wr, &wkopt, &lwork,(int*)&info );
    lwork = (int)wkopt;
    work = new double[lwork];
    /* Solve eigenproblem */
    dsyev_( Vectors, Upper, (int*)&n, L, (int*)&lda, wr, work, &lwork, (int*)&info );


    for(unsigned int i=0;i<n;i++)
        for(unsigned int j=0;j<n;j++)
            vs[i*n+j]=L[j*n+i];

   /* std::cout<<"lapack: "<<std::endl;

    for(unsigned int i=0;i<n;i++)
        std::cout<<" i: "<<i<<" eig : "<<wr[i]<<std::endl;

    for(unsigned int i=0;i<n;i++)
    {
        for(unsigned int j=0;j<n;j++)
        {
            std::cout<<vs[i*(n)+j]<<" ";
        }

        std::cout<<std::endl;
    }

    memcpy(vs,A,sizeof(double)*n*n);
    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',n,vs,lda,wr);
    std::cout<<"lapacke: "<<std::endl;
    for(unsigned int i=0;i<n;i++)
        std::cout<<" i: "<<i<<" eig : "<<wr[i]<<std::endl;

    for(unsigned int i=0;i<n;i++)
    {
        for(unsigned int j=0;j<n;j++)
        {
            std::cout<<vs[i*(n)+j]<<" ";
        }

        std::cout<<std::endl;
    }*/

    delete [] work;
    delete [] L;
    if(info!=0) std::cout<<"lapack eigen solve failed. "<<std::endl;
    return;

}






}// end of namespace


#endif //SFCSORTBENCH_LAPAC_H
