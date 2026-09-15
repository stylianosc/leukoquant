//
//  EigenvaluesCov.h
//  TreeSeg
//
//  Created by Carole Sudre on 15/05/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 07/05/2026.
//
#pragma once

#include "math.h"

class Eigenvalue
{
    
    
    int n;
    
    int issymmetric; /* boolean*/
    
    
    float * d;         /* real part */
    float * e;         /* img part */
    float * V;
    float * H;
    float * ort;
    
    
    // Symmetric Householder reduction to tridiagonal form.
    
    void tred2() {
        
        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        
        for (int j = 0; j < n; j++) {
            d[j] = V[n-1+j*n];
        }
        
        // Householder reduction to tridiagonal form.
        
        for (int i = n-1; i > 0; i--) {
            
            // Scale to avoid under/overflow.
            
            float scale = 0.0;
            float h = 0.0;
            for (int k = 0; k < i; k++) {
                scale = scale + fabs(d[k]);
            }
            if (scale == 0.0) {
                e[i] = d[i-1];
                for (int j = 0; j < i; j++) {
                    d[j] = V[i-1+j*n];
                    V[i+j*n] = 0.0;
                    V[j+i*n] = 0.0;
                }
            } else {
                
                // Generate Householder vector.
                
                for (int k = 0; k < i; k++) {
                    d[k] /= scale;
                    h += d[k] * d[k];
                }
                float f = d[i-1];
                float g = sqrt(h);
                if (f > 0) {
                    g = -g;
                }
                e[i] = scale * g;
                h = h - f * g;
                d[i-1] = f - g;
                for (int j = 0; j < i; j++) {
                    e[j] = 0.0;
                }
                
                // Apply similarity transformation to remaining columns.
                
                for (int j = 0; j < i; j++) {
                    f = d[j];
                    V[j+i*n] = f;
                    g = e[j] + V[j+j*n] * f;
                    for (int k = j+1; k <= i-1; k++) {
                        g += V[k+j*n] * d[k];
                        e[k] += V[k+j*n] * f;
                    }
                    e[j] = g;
                }
                f = 0.0;
                for (int j = 0; j < i; j++) {
                    e[j] /= h;
                    f += e[j] * d[j];
                }
                float hh = f / (h + h);
                for (int j = 0; j < i; j++) {
                    e[j] -= hh * d[j];
                }
                for (int j = 0; j < i; j++) {
                    f = d[j];
                    g = e[j];
                    for (int k = j; k <= i-1; k++) {
                        V[k+j*n] -= (f * e[k] + g * d[k]);
                    }
                    d[j] = V[i-1+j*n];
                    V[i+j*n] = 0.0;
                }
            }
            d[i] = h;
        }
        
        // Accumulate transformations.
        
        for (int i = 0; i < n-1; i++) {
            V[n-1+i*n] = V[i+i*n];
            V[i+i*n] = 1.0;
            float h = d[i+1];
            if (h != 0.0) {
                for (int k = 0; k <= i; k++) {
                    d[k] = V[k+(i+1)*n] / h;
                }
                for (int j = 0; j <= i; j++) {
                    float g = 0.0;
                    for (int k = 0; k <= i; k++) {
                        g += V[k+(i+1)*n] * V[k+j*n];
                    }
                    for (int k = 0; k <= i; k++) {
                        V[k+j*n] -= g * d[k];
                    }
                }
            }
            for (int k = 0; k <= i; k++) {
                V[k+(i+1)*n] = 0.0;
            }
        }
        for (int j = 0; j < n; j++) {
            d[j] = V[n-1+j*n];
            V[n-1+j*n] = 0.0;
        }
        V[n-1+(n-1)*n] = 1.0;
        e[0] = 0.0;
    }
    
    // Symmetric tridiagonal QL algorithm.
    
    void tql2 () {
        
        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        
        for (int i = 1; i < n; i++) {
            e[i-1] = e[i];
        }
        e[n-1] = 0.0;
        
        float f = 0.0;
        float tst1 = 0.0;
        float eps = powf(2.0,-52.0);
        for (int l = 0; l < n; l++) {
            
            // Find small subdiagonal element
            
            tst1 = tst1>=(fabs(d[l]) + fabs(e[l]))?tst1:(fabs(d[l]) + fabs(e[l]));
            int m = l;
            
            // Original while-loop from Java code
            while (m < n) {
                if (fabs(e[m]) <= eps*tst1) {
                    break;
                }
                m++;
            }
            
            
            // If m == l, d[l] is an eigenvalue,
            // otherwise, iterate.
            
            if (m > l) {
                int iter = 0;
                do {
                    iter = iter + 1;  // (Could check iteration count here.)
                    
                    // Compute implicit shift
                    
                    float g = d[l];
                    float p = (d[l+1] - g) / (2.0 * e[l]);
                    float r = sqrt(powf(p,2.0)+powf(1.0,2.0));
                    if (p < 0) {
                        r = -r;
                    }
                    d[l] = e[l] / (p + r);
                    d[l+1] = e[l] * (p + r);
                    float dl1 = d[l+1];
                    float h = g - d[l];
                    for (int i = l+2; i < n; i++) {
                        d[i] -= h;
                    }
                    f = f + h;
                    
                    // Implicit QL transformation.
                    
                    p = d[m];
                    float c = 1.0;
                    float c2 = c;
                    float c3 = c;
                    float el1 = e[l+1];
                    float s = 0.0;
                    float s2 = 0.0;
                    for (int i = m-1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e[i];
                        h = c * p;
                        r = sqrt(powf(p,2)+powf(e[i],2));
                        e[i+1] = s * r;
                        s = e[i] / r;
                        c = p / r;
                        p = c * d[i] - s * g;
                        d[i+1] = h + s * (c * g + s * d[i]);
                        
                        // Accumulate transformation.
                        
                        for (int k = 0; k < n; k++) {
                            h = V[k+(i+1)*n];
                            V[k+(i+1)*n] = s * V[k+i*n] + c * h;
                            V[k+i*n] = c * V[k+i*n] - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * e[l] / dl1;
                    e[l] = s * p;
                    d[l] = c * p;
                    
                    // Check for convergence.
                    
                } while (fabs(e[l]) > eps*tst1);
            }
            d[l] = d[l] + f;
            e[l] = 0.0;
        }
        
        // Sort eigenvalues and corresponding vectors.
        
        for (int i = 0; i < n-1; i++) {
            int k = i;
            float p = d[i];
            for (int j = i+1; j < n; j++) {
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            }
            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                for (int j = 0; j < n; j++) {
                    p = V[j+i*n];
                    V[j+i*n] = V[j+k*n];
                    V[j+k*n] = p;
                }
            }
        }
    }
    
    // Nonsymmetric reduction to Hessenberg form.
    
    void orthes () {
        
        //  This is derived from the Algol procedures orthes and ortran,
        //  by Martin and Wilkinson, Handbook for Auto. Comp.,
        //  Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutines in EISPACK.
        
        int low = 0;
        int high = n-1;
        
        for (int m = low+1; m <= high-1; m++) {
            
            // Scale column.
            
            float scale = 0.0;
            for (int i = m; i <= high; i++) {
                scale = scale + fabs(H[i+(m-1)*n]);
            }
            if (scale != 0.0) {
                
                // Compute Householder transformation.
                
                float h = 0.0;
                for (int i = high; i >= m; i--) {
                    ort[i] = H[i+(m-1)*n]/scale;
                    h += ort[i] * ort[i];
                }
                float g = sqrt(h);
                if (ort[m] > 0) {
                    g = -g;
                }
                h = h - ort[m] * g;
                ort[m] = ort[m] - g;
                
                // Apply Householder similarity transformation
                // H = (I-u*u'/h)*H*(I-u*u')/h)
                
                for (int j = m; j < n; j++) {
                    float f = 0.0;
                    for (int i = high; i >= m; i--) {
                        f += ort[i]*H[i+j*n];
                    }
                    f = f/h;
                    for (int i = m; i <= high; i++) {
                        H[i+j*n] -= f*ort[i];
                    }
                }
                
                for (int i = 0; i <= high; i++) {
                    float f = 0.0;
                    for (int j = high; j >= m; j--) {
                        f += ort[j]*H[i+j*n];
                    }
                    f = f/h;
                    for (int j = m; j <= high; j++) {
                        H[i+j*n] -= f*ort[j];
                    }
                }
                ort[m] = scale*ort[m];
                H[m+(m-1)*n] = scale*g;
            }
        }
        
        // Accumulate transformations (Algol's ortran).
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                V[i+j*n] = (i == j ? 1.0 : 0.0);
            }
        }
        
        for (int m = high-1; m >= low+1; m--) {
            if (H[m+(m-1)*n] != 0.0) {
                for (int i = m+1; i <= high; i++) {
                    ort[i] = H[i+(m-1)*n];
                }
                for (int j = m; j <= high; j++) {
                    float g = 0.0;
                    for (int i = m; i <= high; i++) {
                        g += ort[i] * V[i+j*n];
                    }
                    // Double division avoids possible underflow
                    g = (g / ort[m]) / H[m+(m-1)*n];
                    for (int i = m; i <= high; i++) {
                        V[i+j*n] += g * ort[i];
                    }
                }
            }
        }
    }
    
    
    // Complex scalar division.
    
    float cdivr, cdivi;
    void cdiv(float xr, float xi, float yr, float yi) {
        float r,d;
        if (fabs(yr) > fabs(yi)) {
            r = yi/yr;
            d = yr + r*yi;
            cdivr = (xr + r*xi)/d;
            cdivi = (xi - r*xr)/d;
        } else {
            r = yr/yi;
            d = yi + r*yr;
            cdivr = (r*xr + xi)/d;
            cdivi = (r*xi - xr)/d;
        }
    }
    
    
    // Nonsymmetric reduction from Hessenberg to real Schur form.
    
    void hqr2 () {
        
        //  This is derived from the Algol procedure hqr2,
        //  by Martin and Wilkinson, Handbook for Auto. Comp.,
        //  Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        
        // Initialize
        
        int nn = this->n;
        int n = nn-1;
        int low = 0;
        int high = nn-1;
        float eps = powf(2.0,-52.0);
        float exshift = 0.0;
        float p=0,q=0,r=0,s=0,z=0,t,w,x,y;
        
        // Store roots isolated by balanc and compute matrix norm
        
        float norm = 0.0;
        for (int i = 0; i < nn; i++) {
            if ((i < low) || (i > high)) {
                d[i] = H[i+i*n];
                e[i] = 0.0;
            }
            for (int j = max(i-1,0); j < nn; j++) {
                norm = norm + fabs(H[i+j*n]);
            }
        }
        
        // Outer loop over eigenvalue index
        
        int iter = 0;
        while (n >= low) {
            
            // Look for single small sub-diagonal element
            
            int l = n;
            while (l > low) {
                s = fabs(H[l-1+(l-1)*nn]) + fabs(H[l+l*nn]);
                if (s == 0.0) {
                    s = norm;
                }
                if (fabs(H[l+(l-1)*nn]) < eps * s) {
                    break;
                }
                l--;
            }
            
            // Check for convergence
            // One root found
            
            if (l == n) {
                H[n+n*nn] = H[n+n*nn] + exshift;
                d[n] = H[n+n*nn];
                e[n] = 0.0;
                n--;
                iter = 0;
                
                // Two roots found
                
            } else if (l == n-1) {
                w = H[n+(n-1)*nn] * H[n-1+n*nn];
                p = (H[n-1+(n-1)*nn] - H[n+n*nn]) / 2.0;
                q = p * p + w;
                z = sqrt(fabs(q));
                H[n+n*nn] = H[n+n*nn] + exshift;
                H[n-1+(n-1)*nn] = H[n-1+(n-1)*nn] + exshift;
                x = H[n+n*nn];
                
                // Real pair
                
                if (q >= 0) {
                    if (p >= 0) {
                        z = p + z;
                    } else {
                        z = p - z;
                    }
                    d[n-1] = x + z;
                    d[n] = d[n-1];
                    if (z != 0.0) {
                        d[n] = x - w / z;
                    }
                    e[n-1] = 0.0;
                    e[n] = 0.0;
                    x = H[n+(n-1)*nn];
                    s = fabs(x) + fabs(z);
                    p = x / s;
                    q = z / s;
                    r = sqrt(p * p+q * q);
                    p = p / r;
                    q = q / r;
                    
                    // Row modification
                    
                    for (int j = n-1; j < nn; j++) {
                        z = H[n-1+j*nn];
                        H[n-1+j*nn] = q * z + p * H[n+j*nn];
                        H[n+j*nn] = q * H[n+j*nn] - p * z;
                    }
                    
                    // Column modification
                    
                    for (int i = 0; i <= n; i++) {
                        z = H[i+(n-1)*nn];
                        H[i+(n-1)*nn] = q * z + p * H[i+n*nn];
                        H[i+n*nn] = q * H[i+n*nn] - p * z;
                    }
                    
                    // Accumulate transformations
                    
                    for (int i = low; i <= high; i++) {
                        z = V[i+(n-1)*n];
                        V[i+(n-1)*nn] = q * z + p * V[i+n*nn];
                        V[i+n*nn] = q * V[i+n*nn] - p * z;
                    }
                    
                    // Complex pair
                    
                } else {
                    d[n-1] = x + p;
                    d[n] = x + p;
                    e[n-1] = z;
                    e[n] = -z;
                }
                n = n - 2;
                iter = 0;
                
                // No convergence yet
                
            } else {
                
                // Form shift
                
                x = H[n+n*nn];
                y = 0.0;
                w = 0.0;
                if (l < n) {
                    y = H[n-1+(n-1)*nn];
                    w = H[n+(n-1)*nn] * H[n-1+n*nn];
                }
                
                // Wilkinson's original ad hoc shift
                
                if (iter == 10) {
                    exshift += x;
                    for (int i = low; i <= n; i++) {
                        H[i+i*nn] -= x;
                    }
                    s = fabs(H[n+(n-1)*nn]) + fabs(H[n-1+(n-2)*nn]);
                    x = y = 0.75 * s;
                    w = -0.4375 * s * s;
                }
                
                // MATLAB's new ad hoc shift
                
                if (iter == 30) {
                    s = (y - x) / 2.0;
                    s = s * s + w;
                    if (s > 0) {
                        s = sqrt(s);
                        if (y < x) {
                            s = -s;
                        }
                        s = x - w / ((y - x) / 2.0 + s);
                        for (int i = low; i <= n; i++) {
                            H[i+i*nn] -= s;
                        }
                        exshift += s;
                        x = y = w = 0.964;
                    }
                }
                
                iter = iter + 1;   // (Could check iteration count here.)
                
                // Look for two consecutive small sub-diagonal elements
                
                int m = n-2;
                while (m >= l) {
                    z = H[m+m*nn];
                    r = x - z;
                    s = y - z;
                    p = (r * s - w) / H[m+1+m*nn] + H[m+(m+1)*nn];
                    q = H[m+1+(m+1)*nn] - z - r - s;
                    r = H[m+2+(m+1)*nn];
                    s = fabs(p) + fabs(q) + fabs(r);
                    p = p / s;
                    q = q / s;
                    r = r / s;
                    if (m == l) {
                        break;
                    }
                    if (fabs(H[m+(m-1)*nn]) * (fabs(q) + fabs(r)) <
                        eps * (fabs(p) * (fabs(H[m-1+(m-1)*nn]) + fabs(z) +
                                          fabs(H[m+1+(m+1)*nn])))) {
                        break;
                    }
                    m--;
                }
                
                for (int i = m+2; i <= n; i++) {
                    H[i+(i-2)*nn] = 0.0;
                    if (i > m+2) {
                        H[i+(i-3)*nn] = 0.0;
                    }
                }
                
                // Double QR step involving rows l:n and columns m:n
                
                for (int k = m; k <= n-1; k++) {
                    int notlast = (k != n-1);
                    if (k != m) {
                        p = H[k+(k-1)*nn];
                        q = H[k+1+(k-1)*nn];
                        r = (notlast ? H[k+2+(k-1)*nn] : 0.0);
                        x = fabs(p) + fabs(q) + fabs(r);
                        if (x != 0.0) {
                            p = p / x;
                            q = q / x;
                            r = r / x;
                        }
                    }
                    if (x == 0.0) {
                        break;
                    }
                    s = sqrt(p * p + q * q + r * r);
                    if (p < 0) {
                        s = -s;
                    }
                    if (s != 0) {
                        if (k != m) {
                            H[k+(k-1)*nn] = -s * x;
                        } else if (l != m) {
                            H[k+(k-1)*nn] = -H[k+(k-1)*nn];
                        }
                        p = p + s;
                        x = p / s;
                        y = q / s;
                        z = r / s;
                        q = q / p;
                        r = r / p;
                        
                        // Row modification
                        
                        for (int j = k; j < nn; j++) {
                            p = H[k+j*nn] + q * H[k+1+j*nn];
                            if (notlast) {
                                p = p + r * H[k+2+j*nn];
                                H[k+2+j*nn] = H[k+2+j*nn] - p * z;
                            }
                            H[k+j*nn] = H[k+j*nn] - p * x;
                            H[k+1+j*nn] = H[k+1+j*nn] - p * y;
                        }
                        
                        // Column modification
                        
                        for (int i = 0; i <= min(n,k+3); i++) {
                            p = x * H[i+k*nn] + y * H[i+(k+1)*nn];
                            if (notlast) {
                                p = p + z * H[i+(k+2)*nn];
                                H[i+(k+2)*nn] = H[i+(k+2)*nn] - p * r;
                            }
                            H[i+k*nn] = H[i+k*nn] - p;
                            H[i+(k+1)*nn] = H[i+(k+1)*nn] - p * q;
                        }
                        
                        // Accumulate transformations
                        
                        for (int i = low; i <= high; i++) {
                            p = x * V[i+k*nn] + y * V[i+(k+1)*nn];
                            if (notlast) {
                                p = p + z * V[i+(k+2)*nn];
                                V[i+(k+2)*nn] = V[i+(k+2)*nn] - p * r;
                            }
                            V[i+k*nn] = V[i+k*nn] - p;
                            V[i+(k+1)*nn] = V[i+(k+1)*nn] - p * q;
                        }
                    }  // (s != 0)
                }  // k loop
            }  // check convergence
        }  // while (n >= low)
        
        // Backsubstitute to find vectors of upper triangular form
        
        if (norm == 0.0) {
            return;
        }
        
        for (n = nn-1; n >= 0; n--) {
            p = d[n];
            q = e[n];
            
            // Real vector
            
            if (q == 0) {
                int l = n;
                H[n+n*nn] = 1.0;
                for (int i = n-1; i >= 0; i--) {
                    w = H[i+i*nn] - p;
                    r = 0.0;
                    for (int j = l; j <= n; j++) {
                        r = r + H[i+j*nn] * H[j+n*nn];
                    }
                    if (e[i] < 0.0) {
                        z = w;
                        s = r;
                    } else {
                        l = i;
                        if (e[i] == 0.0) {
                            if (w != 0.0) {
                                H[i+n*nn] = -r / w;
                            } else {
                                H[i+n*nn] = -r / (eps * norm);
                            }
                            
                            // Solve real equations
                            
                        } else {
                            x = H[i+(i+1)*nn];
                            y = H[i+1+i*nn];
                            q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                            t = (x * s - z * r) / q;
                            H[i+n*nn] = t;
                            if (fabs(x) > fabs(z)) {
                                H[i+1+n*nn] = (-r - w * t) / x;
                            } else {
                                H[i+1+n*nn] = (-s - y * t) / z;
                            }
                        }
                        
                        // Overflow control
                        
                        t = fabs(H[i+n*nn]);
                        if ((eps * t) * t > 1) {
                            for (int j = i; j <= n; j++) {
                                H[j+n*nn] = H[j+n*nn] / t;
                            }
                        }
                    }
                }
                
                // Complex vector
                
            } else if (q < 0) {
                int l = n-1;
                
                // Last vector component imaginary so matrix is triangular
                
                if (fabs(H[n+(n-1)*nn]) > fabs(H[n-1+n*nn])) {
                    H[n-1+(n-1)*nn] = q / H[n+(n-1)*nn];
                    H[n-1+n*nn] = -(H[n+n*nn] - p) / H[n+(n-1)*nn];
                } else {
                    cdiv(0.0,-H[n-1+n*nn],H[n-1+(n-1)*nn]-p,q);
                    H[n-1+(n-1)*nn] = cdivr;
                    H[n-1+n*nn] = cdivi;
                }
                H[n+(n-1)*nn] = 0.0;
                H[n+n*nn] = 1.0;
                for (int i = n-2; i >= 0; i--) {
                    float ra,sa,vr,vi;
                    ra = 0.0;
                    sa = 0.0;
                    for (int j = l; j <= n; j++) {
                        ra = ra + H[i+j*nn] * H[j+(n-1)*nn];
                        sa = sa + H[i+j*nn] * H[j+n*nn];
                    }
                    w = H[i+i*nn] - p;
                    
                    if (e[i] < 0.0) {
                        z = w;
                        r = ra;
                        s = sa;
                    } else {
                        l = i;
                        if (e[i] == 0) {
                            cdiv(-ra,-sa,w,q);
                            H[i+(n-1)*nn] = cdivr;
                            H[i+n*nn] = cdivi;
                        } else {
                            
                            // Solve complex equations
                            
                            x = H[i+(i+1)*nn];
                            y = H[i+1+i*nn];
                            vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                            vi = (d[i] - p) * 2.0 * q;
                            if ((vr == 0.0) && (vi == 0.0)) {
                                vr = eps * norm * (fabs(w) + fabs(q) +
                                                   fabs(x) + fabs(y) + fabs(z));
                            }
                            cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                            H[i+(n-1)*nn] = cdivr;
                            H[i+n*nn] = cdivi;
                            if (fabs(x) > (fabs(z) + fabs(q))) {
                                H[i+1+(n-1)*nn] = (-ra - w * H[i+(n-1)*nn] + q * H[i+n*nn]) / x;
                                H[i+1+n*nn] = (-sa - w * H[i+n*nn] - q * H[i+(n-1)*nn]) / x;
                            } else {
                                cdiv(-r-y*H[i+(n-1)*nn],-s-y*H[i+n*nn],z,q);
                                H[i+1+(n-1)*nn] = cdivr;
                                H[i+1+n*nn] = cdivi;
                            }
                        }
                        
                        // Overflow control
                        
                        t = max(fabs(H[i+(n-1)*nn]),fabs(H[i+n*nn]));
                        if ((eps * t) * t > 1) {
                            for (int j = i; j <= n; j++) {
                                H[j+(n-1)*nn] = H[j+(n-1)*nn] / t;
                                H[j+n*nn] = H[j+n*nn] / t;
                            }
                        }
                    }
                }
            }
        }
        
        // Vectors of isolated roots
        
        for (int i = 0; i < nn; i++) {
            if (i < low || i > high) {
                for (int j = i; j < nn; j++) {
                    V[i+j*nn] = H[i+j*nn];
                }
            }
        }
        
        // Back transformation to get eigenvectors of original matrix
        
        for (int j = nn-1; j >= low; j--) {
            for (int i = low; i <= high; i++) {
                z = 0.0;
                for (int k = low; k <= min(j,high); k++) {
                    z = z + V[i+k*nn] * H[k+j*nn];
                }
                V[i+j*nn] = z;
            }
        }
    }
    
public:
    
    
    
    Eigenvalue(const float * A, int n) {
        //n = A.dim2();
        this->n=n;
        V = new float[n*n];//{0};//Array2D<Real>(n,n);
        for (int i=0; i<n*n; i++) {
            V[i]=0;
        }
        d = new float[n];//{0};//Array1D<Real>(n);
        for (int i=0; i<n; i++) {
            d[i]=0;
        }
        e = new float[n];//{0};//Array1D<Real>(n);
        for (int i=0; i<n; i++) {
            e[i]=0;
        }
        
        issymmetric = 1;
        for (int j = 0; (j < n) && issymmetric; j++) {
            for (int i = 0; (i < n) && issymmetric; i++) {
                issymmetric = (A[i+j*n] == A[j+i*n]);
            }
        }
        
        if (issymmetric) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    V[i+j*n] = A[i+j*n];
                }
            }
            
            // Tridiagonalize.
            tred2();
            
            // Diagonalize.
            tql2();
            
        } else {
            H = new float[n*n];//{0};//TNT::Array2D<R(n,n);
            for (int i=0; i<n*n; i++) {
                H[i]=0;
            }
            ort = new float[n];//{0};//TNT::Array1D<Real>(n);
            for (int i=0; i<n; i++) {
                ort[i]=0;
            }
            
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    H[i+j*n] = A[i+j*n];
                }
            }
            
            // Reduce to Hessenberg form.
            orthes();
            
            // Reduce Hessenberg to real Schur form.
            hqr2();
        }
    }
    
    
    
//    void getV (TNT::Array2D<Real> &V_) {
//        V_ = V;
//        return;
//    }
//    
//    
//    void getRealEigenvalues (TNT::Array1D<Real> &d_) {
//        d_ = d;
//        return ;
//    }
//    
//    void getImagEigenvalues (TNT::Array1D<Real> &e_) {
//        e_ = e;
//        return;
//    }
//    
//    
//    void getD (TNT::Array2D<Real> &D) {
//        D = Array2D<Real>(n,n);
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                D[i][j] = 0.0;
//            }
//            D[i+i*n] = d[i];
//            if (e[i] > 0) {
//                D[i+(i+1)*n] = e[i];
//            } else if (e[i] < 0) {
//                D[i+(i-1)*n] = e[i];
//            }
//        }
//    }
    float * getV(){
        return this->V;
    }
    
    float * getRealEigenValues(){
        return this->d;
    }
    
    float * getImagEigenValues(){
        return this->e;
    }
    
    float * getD(){
        int n=this->n;
        float * D = new float[n*n];//{0};// Array2D<Real>(n,n);
        for (int i=0; i<n*n; i++) {
            D[i]=0;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D[i+j*n] = 0.0;
            }
            D[i+i*n] = d[i];
            if (e[i] > 0) {
                D[i+(i+1)*n] = e[i];
            } else if (e[i] < 0) {
                D[i+(i-1)*n] = e[i];
            }
        }
        return D;
    }
};


