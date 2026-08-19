//
//  SVDCov.h
//  TreeSeg
//
//  Created by Carole Sudre on 15/05/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 07/05/2026.
//

#pragma once

#include "math.h"
#include <iostream>
using namespace std;


class SVD
{
    
    
    float * U;
    float * V ;
    float * s;
    
    int m, n;
    
public:
    ~SVD(){
        if (this->U!=NULL) {
            delete [] this->U;
            this->U=NULL;
        }
        if (this->V!=NULL) {
            delete [] this->V;
            this->V=NULL;
        }
        if (this->s!=NULL) {
            delete [] this->s;
            this->s=NULL;
        }
    }
    
    SVD (const float * MatrixToSVD , int m,int n) {
        
        
        
        int nu = m>=n?n:m;
        int tmpsize=(m+1)>=n?n:(m+1);
        this->m=m;
        this->n=n;
        this->s = new float[tmpsize];//{0};
        for (int i=0; i<tmpsize; i++) {
            this->s[i]=0;
        }
        this->U = new float[m*nu];//{0};
        for (int i=0; i<m*nu; i++) {
            this->U[i]=0;
        }
        this->V = new float[n*n];//{0};
        for (int i=0; i<n*n; i++) {
            this->V[i]=0;
        }
        float * e=new float[n];
        float * work=new float[m];
        float * A=new float[m*n];
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                A[i+m*j]=MatrixToSVD[i+j*m];
            }
        }
        
        int wantu = 1;                                    /* boolean */
        int wantv = 1;                                    /* boolean */
        int i=0, j=0, k=0;
        
        // Reduce A to bidiagonal form, storing the diagonal elements
        // in s and the super-diagonal elements in e.
        
        int nct = m-1>=n?n:(m-1);
        int tmp=(n-2)>=m?m:(n-2);
        int nrt = tmp>=0?tmp:0;
        int maxn=nrt>=nct?nrt:nct;
        for (k = 0; k < maxn; k++) {
            if (k < nct) {
                
                // Compute the transformation for the k-th column and
                // place the k-th diagonal in s[k].
                // Compute 2-norm of k-th column without under/overflow.
                s[k]=0;
                
                
                for (i = k; i < m; i++) {
                    s[k]=sqrt(powf(s[k],2)+powf(A[i+m*k],2));
                    
                    //s[k] = hypot(s[k],A[i][k]);
                }
                
                if (s[k] != 0.0) {
                    if (A[k+m*k] < 0.0) {
                        s[k]=-s[k];
                    }
                    for (i = k; i < m; i++) {
                        A[i+k*m]/=s[k];
                        
                    }
                    
                    A[k+m*k] += 1.0;
                }
                
                s[k] = -s[k];
            }
            for (j = k+1; j < n; j++) {
                if ((k < nct) && (this->s[k] != 0.0))  {
                    
                    // Apply the transformation.
                    
                    float t = 0;
                    for (i = k; i < m; i++) {
                        t+=A[i+k*m]*A[i+m*j];
                        //t += A[i][k]*A[i][j];
                    }
                    
                    t = -t/A[k+k*m];
                    for (i = k; i < m; i++) {
                        
                        A[i+j*m]+= t*A[i+m*k];
                    }
                }
                
                // Place the k-th row of A into e for the
                // subsequent calculation of the row transformation.
                //e.setvalue(j,1,A.getvalue(k,j));
                e[j] = A[k+j*m];
            }
            if (wantu & (k < nct)) {
                
                // Place the transformation in U for subsequent back
                // multiplication.
                
                for (i = k; i < m; i++) {
                    //U.setvalue(i,k,A.getvalue(i,k));
                    U[i+k*m] = A[i+k*m];
                }
            }
            if (k < nrt) {
                
                // Compute the k-th row transformation and place the
                // k-th super-diagonal in e[k].
                // Compute 2-norm without under/overflow.
                //e.setvalue(k,1,0);
                e[k] = 0;
                for (i = k+1; i < n; i++) {
                    //e.setvalue(k,1,sqrt(powf(e.getvalue(k,1),2)+powf(e.getvalue(i,1),2)))
                    e[k] = sqrt(powf(e[k],2)+powf(e[i],2));
                }
                if (e[k] != 0.0) {
                    if (e[k+1] < 0.0) {
                        //e.setvalue(k,1,-e.getvalue(k,1));
                        e[k] = -e[k];
                    }
                    for (i = k+1; i < n; i++) {
                        //e.setvalue(i,1,e.getvalue(i,1)/e.getvalue(k,1));
                        e[i] /= e[k];
                    }
                    
                    e[k+1] += 1.0;
                }
                e[k] = -e[k];
                if ((k+1 < m) & (e[k] != 0.0)) {
                    
                    // Apply the transformation.
                    
                    for (i = k+1; i < m; i++) {
                        work[i] = 0.0;
                    }
                    for (j = k+1; j < n; j++) {
                        for (i = k+1; i < m; i++) {
                            work[i] += e[j]*A[i+j*m];
                        }
                    }
                    for (j = k+1; j < n; j++) {
                        float t = -e[j]/e[k+1];
                        for (i = k+1; i < m; i++) {
                            A[i+j*m] += t*work[i];
                        }
                    }
                }
                if (wantv) {
                    
                    // Place the transformation in V for subsequent
                    // back multiplication.
                    
                    for (i = k+1; i < n; i++) {
                        this->V[i+k*m] = e[i];
                    }
                }
            }
        }
        
        // Set up the final bidiagonal matrix or order p.
        
        int p = n>=(m+1)?(m+1):n;
        if (nct < n) {
            this->s[nct] = A[nct+nct*m];
        }
        if (m < p) {
            this->s[p-1] = 0.0;
        }
        if (nrt+1 < p) {
            e[nrt] = A[nrt+(p-1)*m];
        }
        e[p-1] = 0.0;
        
        // If required, generate U.
        
        if (wantu) {
            for (j = nct; j < nu; j++) {
                for (i = 0; i < m; i++) {
                    U[i+j*m] = 0.0;
                }
                U[j+j*m] = 1.0;
            }
            for (k = nct-1; k >= 0; k--) {
                if (this->s[k] != 0.0) {
                    for (j = k+1; j < nu; j++) {
                        float t = 0;
                        for (i = k; i < m; i++) {
                            t += U[i+k*m]*U[i+j*m];
                        }
                        t = -t/U[k+k*m];
                        for (i = k; i < m; i++) {
                            U[i+j*m] += t*U[i+k*m];
                        }
                    }
                    for (i = k; i < m; i++ ) {
                        U[i+k*m] = -U[i+k*m];
                    }
                    U[k+k*m] = 1.0 + U[k+k*m];
                    for (i = 0; i < k-1; i++) {
                        U[i+k*m] = 0.0;
                    }
                } else {
                    for (i = 0; i < m; i++) {
                        U[i+k*m] = 0.0;
                    }
                    U[k+k*m] = 1.0;
                }
            }
        }
        
        // If required, generate V.
        
        if (wantv) {
            for (k = n-1; k >= 0; k--) {
                if ((k < nrt) & (e[k] != 0.0)) {
                    for (j = k+1; j < nu; j++) {
                        float t = 0;
                        for (i = k+1; i < n; i++) {
                            t += this->V[i+k*n]*this->V[i+j*n];
                        }
                        t = -t/V[k+1+k*n];
                        for (i = k+1; i < n; i++) {
                            V[i+j*n] += t*V[i+k*n];
                        }
                    }
                }
                for (i = 0; i < n; i++) {
                    V[i+k*n] = 0.0;
                }
                V[k+k*n] = 1.0;
            }
        }
        
        // Main iteration loop for the singular values.
        
        int pp = p-1;
        int iter = 0;
        float eps = pow(2.0,-52.0);
        while (p > 0) {
            int k=0;
            int kase=0;
            
            // Here is where a test for too many iterations would go.
            if (iter>10000){
                cout << "Reached 10000 in SVD"<<endl;
                if (m==n && m==3){
                    for (int i=0;i<9;i++) {
                        cout << MatrixToSVD[i] << " "<< endl;
                    }
                }
                break;
            }
            // This section of the program inspects for
            // negligible elements in the s and e arrays.  On
            // completion the variables kase and k are set as follows.
            
            // kase = 1     if s(p) and e[k-1] are negligible and k<p
            // kase = 2     if s(k) is negligible and k<p
            // kase = 3     if e[k-1] is negligible, k<p, and
            //              s(k), ..., s(p) are not negligible (qr step).
            // kase = 4     if e(p-1) is negligible (convergence).
            
            for (k = p-2; k >= -1; k--) {
                if (k == -1) {
                    break;
                }
                if (fabs(e[k]) <= eps*(fabs(s[k]) + fabs(s[k+1]))) {
                    e[k] = 0.0;
                    break;
                }
            }
            if (k == p-2) {
                kase = 4;
            } else {
                int ks;
                for (ks = p-1; ks >= k; ks--) {
                    if (ks == k) {
                        break;
                    }
                    float t = (ks != p ? fabs(e[ks]) : 0.) +
                    (ks != k+1 ? fabs(e[ks-1]) : 0.);
                    if (fabs(s[ks]) <= eps*t)  {
                        s[ks] = 0.0;
                        break;
                    }
                }
                if (ks == k) {
                    kase = 3;
                } else if (ks == p-1) {
                    kase = 1;
                } else {
                    kase = 2;
                    k = ks;
                }
            }
            k++;
            
            // Perform the task indicated by kase.
            
            switch (kase) {
                    
                    // Deflate negligible s(p).
                    
                case 1: {
                    float f = e[p-2];
                    e[p-2] = 0.0;
                    for (j = p-2; j >= k; j--) {
                        float t =sqrt(powf(s[j],2)+powf(f,2));
                        float cs = s[j]/t;
                        float sn = f/t;
                        s[j] = t;
                        if (j != k) {
                            f = -sn*e[j-1];
                            e[j-1] = cs*e[j-1];
                        }
                        if (wantv) {
                            for (i = 0; i < n; i++) {
                                t = cs*V[i+j*n] + sn*V[i+(p-1)*n];
                                V[i+(p-1)*n] = -sn*V[i+j*n] + cs*V[i+(p-1)*n];
                                V[i+j*n] = t;
                            }
                        }
                    }
                }
                    break;
                    
                    // Split at negligible s(k).
                    
                case 2: {
                    float f = e[k-1];
                    e[k-1] = 0.0;
                    for (j = k; j < p; j++) {
                        float t =sqrt(powf(s[j],2)+powf(f,2));
                        float cs = s[j]/t;
                        float sn = f/t;
                        s[j] = t;
                        f = -sn*e[j];
                        e[j] = cs*e[j];
                        if (wantu) {
                            for (i = 0; i < m; i++) {
                                t = cs*U[i+j*m] + sn*U[i+(k-1)*m];
                                U[i+(k-1)*m] = -sn*U[i+j*m] + cs*U[i+(k-1)*m];
                                U[i+j*m] = t;
                            }
                        }
                    }
                }
                    break;
                    
                    // Perform one qr step.
                    
                case 3: {
                    
                    // Calculate the shift.
                    
                    float scale = max(max(max(max(
                                                  fabs(s[p-1]),fabs(s[p-2])),fabs(e[p-2])),
                                          fabs(s[k])),fabs(e[k]));
                    float sp = s[p-1]/scale;
                    float spm1 = s[p-2]/scale;
                    float epm1 = e[p-2]/scale;
                    float sk = s[k]/scale;
                    float ek = e[k]/scale;
                    float b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
                    float c = (sp*epm1)*(sp*epm1);
                    float shift = 0.0;
                    if ((b != 0.0) | (c != 0.0)) {
                        shift = sqrt(b*b + c);
                        if (b < 0.0) {
                            shift = -shift;
                        }
                        shift = c/(b + shift);
                    }
                    float f = (sk + sp)*(sk - sp) + shift;
                    float g = sk*ek;
                    
                    // Chase zeros.
                    
                    for (j = k; j < p-1; j++) {
                        float t = sqrt(powf(f,2)+powf(g,2));
                        float cs = f/t;
                        float sn = g/t;
                        if (j != k) {
                            e[j-1] = t;
                        }
                        f = cs*s[j] + sn*e[j];
                        e[j] = cs*e[j] - sn*s[j];
                        g = sn*s[j+1];
                        s[j+1] = cs*s[j+1];
                        if (wantv) {
                            for (i = 0; i < n; i++) {
                                t = cs*V[i+j*n] + sn*V[i+(j+1)*n];
                                V[i+(j+1)*n] = -sn*V[i+j*n] + cs*V[i+(j+1)*n];
                                V[i+j*n] = t;
                            }
                        }
                        t = sqrt(powf(f,2)+powf(g,2));
                        cs = f/t;
                        sn = g/t;
                        s[j] = t;
                        f = cs*e[j] + sn*s[j+1];
                        s[j+1] = -sn*e[j] + cs*s[j+1];
                        g = sn*e[j+1];
                        e[j+1] = cs*e[j+1];
                        if (wantu && (j < m-1)) {
                            for (i = 0; i < m; i++) {
                                t = cs*U[i+j*m] + sn*U[i+(j+1)*m];
                                U[i+(j+1)*m] = -sn*U[i+j*m] + cs*U[i+(j+1)*m];
                                U[i+j*m] = t;
                            }
                        }
                    }
                    e[p-2] = f;
                    iter = iter + 1;
                }
                    break;
                    
                    // Convergence.
                    
                case 4: {
                    
                    // Make the singular values positive.
                    
                    if (s[k] <= 0.0) {
                        s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                        if (wantv) {
                            for (i = 0; i <= pp; i++) {
                                V[i+k*n] = -V[i+k*n];
                            }
                        }
                    }
                    
                    // Order the singular values.
                    
                    while (k < pp) {
                        if (s[k] >= s[k+1]) {
                            break;
                        }
                        float t = s[k];
                        s[k] = s[k+1];
                        s[k+1] = t;
                        if (wantv && (k < n-1)) {
                            for (i = 0; i < n; i++) {
                                t = V[i+(k+1)*n]; V[i+(k+1)*n] = V[i+k*n]; V[i+k*n] = t;
                            }
                        }
                        if (wantu && (k < m-1)) {
                            for (i = 0; i < m; i++) {
                                t = U[i+(k+1)*m]; U[i+(k+1)*m] = U[i+k*m]; U[i+k*m] = t;
                            }
                        }
                        k++;
                    }
                    iter = 0;
                    p--;
                }
                    break;
            }
        }
        if(e!=NULL){
            delete [] e;
            e=NULL;
        }
        if(work!=NULL){
            delete [] work;
            work=NULL;
        }
        if(A!=NULL){
            delete [] A;
            A=NULL;
        }
    }
    
    float * getU(){
        return this->U;
//        int minm = min(m+1,n);
//        
//        float * A = new float[m* minm];
//        
//        for (int i=0; i<m; i++)
//            for (int j=0; j<minm; j++)
//                A[i+j*m] = U[i+j*m];
//        return A;
    }
    
//    void getU (float & A) 
//    {
//        int minm = min(m+1,n);
//        
//        float * A = new float[m* minm];
//        
//        for (int i=0; i<m; i++)
//            for (int j=0; j<minm; j++)
//                A[i+j*m] = U[i+j*m];
//        
//    }
    
    /* Return the right singular vectors */
    
//    00475    void getV (float & A) 
//    {
//        A = V;
//    }
    
    float * getV(){
        return this->V;
    }
    
    
//    00482    void getSingularValues (float &x) 
//    {
//        x = s;
//    }
    
    float * getSingularValues(){
        return this->s;
    }
    
    float * getS(){
        float *A = new float [n*n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i+j*n] = 0.0;
            }
            A[i+i*n] = s[i];
        }
        return A;
    }
//    00491    void getS (float &A) {
//
//    }
    
    
    double norm2 () {
        return s[0];
    }
    
    
    double cond () {
        return s[0]/s[min(m,n)-1];
    }
    
    
    int rank () 
    {
        double eps = pow(2.0,-52.0);
        double tol = max(m,n)*s[0]*eps;
        int r = 0;
        for (int i = 0; i < min(m,n); i++) {
            if (s[i] > tol) {
                r++;
            }
        }
        return r;
    }
};

