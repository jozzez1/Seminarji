#include <iostream>
#include <cmath>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "../dlib-18.10/dlib/optimization.h"

typedef Eigen::Matrix<long double, 3, 1> Vector3Ld;
typedef Eigen::Matrix<long double, 3, 3> Matrix3Ld;
typedef Eigen::Matrix<long double, 12, 1> ChiVecLd;
typedef Eigen::Matrix<long double, 13, 1> parameterVectorLd;

////////////////////////////////////////////
/* This segment creates the mass matrices */
////////////////////////////////////////////

long double
alpha (long double a11, long double a12, long double u1, long double u2)
{
    return 1.0L/(1.0L + a11*u1 + a12*u2);
}

Vector3Ld
vector_x (long double z1, long double z2, long double z3,
        long double w1, long double w2, long double w3,
        long double a11, long double a12,
        long double u1, long double u2)
{
    Vector3Ld z, w;
    z << z1,z2,z3;
    w << w1,w2,w3;

    return z + alpha(a11, a12, u1, u2)*w;
}

long double
theta (Vector3Ld xg, Vector3Ld xh)
{
    long double lambda_inv = sqrtl(1.0L + (xg.transpose() * xh));
    return 1.0L/(lambda_inv * (1.0L + lambda_inv));
}

Matrix3Ld
kronecker_product (Vector3Ld v, Vector3Ld w)
{
    Matrix3Ld M;
    M << v(0)*w(0), v(0)*w(1), v(0)*w(2),
         v(1)*w(0), v(1)*w(1), v(1)*w(2),
         v(2)*w(0), v(2)*w(1), v(2)*w(2);

    return M;
}

Matrix3Ld
matrix_Lambda (Vector3Ld xg, Vector3Ld xh)
{
    Matrix3Ld M, I;
    I.setIdentity(3, 3);
    M = I - theta(xg, xh)*kronecker_product(xg, xh);
    return M;
}

Matrix3Ld
mass_matrix (long double y1, long double y2, long double y3,
        long double z1, long double z2, long double z3,
        long double w1, long double w2, long double w3,
        long double u1, long double u2, long double vev,
        long double a11, long double a12,
        long double a21, long double a22)
{
    Matrix3Ld Y, Lg, Lh, M;
    Vector3Ld xg, xh;

    Y << y1, 0.0L, 0.0L, 0.0L, y2, 0.0L, 0.0L, 0.0L, y3;
    xg = vector_x (z1, z2, z3, w1, w2, w3, a11, a12, u1, u2);
    xh = vector_x (z1, z2, z3, w1, w2, w3, a21, a22, u1, u2);
    Lg = matrix_Lambda (xg, xg);
    Lh = matrix_Lambda (xh, xh);

    M = vev * Lg * (Y + kronecker_product(xg, xh)) * Lh;
    return M;
}

Matrix3Ld
up_mass_matrix_squared (long double y1, long double y2, long double y3,
        long double z1, long double z2, long double z3,
        long double w1, long double w2, long double w3,
        long double u1, long double u2, long double vev)
{
    Matrix3Ld MU;
    MU = mass_matrix (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vev, 1.0L, 1.0L, 0.0L, -1.0L);
    return MU.transpose() * MU;
}

Matrix3Ld
down_mass_matrix_squared (long double y1, long double y2, long double y3,
        long double z1, long double z2, long double z3,
        long double w1, long double w2, long double w3,
        long double u1, long double u2, long double vev)
{
    Matrix3Ld MD;
    MD = mass_matrix (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vev, -1.0L, 1.0L, 0.0L, -1.0L);
    return MD.transpose() * MD;
}

Matrix3Ld
electron_mass_matrix_squared (long double y1, long double y2, long double y3,
        long double z1, long double z2, long double z3,
        long double w1, long double w2, long double w3,
        long double u1, long double u2, long double vev)
{
    Matrix3Ld ME;
    ME = mass_matrix (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vev, -1.0L, -3.0L, 0.0L, 3.0L);
    return ME.transpose() * ME;
}

////////////////////////////////////////////////
/* This part is for the minimization function */
////////////////////////////////////////////////

ChiVecLd
vector_B (long double y1, long double y2, long double y3,
        long double z1, long double z2, long double z3,
        long double w1, long double w2, long double w3,
        long double u1, long double u2,
        long double vu, long double vd)
{
    ChiVecLd B, A, tmp;
    Matrix3Ld MU2, MD2, ME2, LU, LD, VCKM;
    Vector3Ld eig;

    MU2 = up_mass_matrix_squared (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu);
    MD2 = down_mass_matrix_squared (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vd);
    ME2 = electron_mass_matrix_squared (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, (-3)*vd);
    Eigen::SelfAdjointEigenSolver <Matrix3Ld> es;

    // B accepts theoretical values
    es.compute (MU2, Eigen::ComputeEigenvectors);
    eig  = es.eigenvalues();
    LU   = es.eigenvectors();
    assert (fabsl(LU.determinant() - 1.0L) < 1e-6L);
    tmp(0) = eig(0);
    tmp(1) = eig(1);
    tmp(2) = eig(2);

    es.compute (MD2, Eigen::ComputeEigenvectors);
    eig  = es.eigenvalues();
    LD   = es.eigenvectors();
    assert (fabsl(LD.determinant() - 1.0L) < 1e-6L);
    tmp(3) = eig(0);
    tmp(4) = eig(1);
    tmp(5) = eig(2);

    VCKM = LU.transpose() * LD;
    tmp(6) = VCKM(1,2);
    tmp(7) = VCKM(1,3);
    tmp(8) = VCKM(2,3);

    es.compute (ME2, Eigen::ComputeEigenvectors);
    eig  = es.eigenvalues();

    tmp(9) = eig(0);
    tmp(10)= eig(1);
    tmp(11)= eig(2);

    // A are the experimental values, masses are in MeV
    const long double mu   = 0.4565L,
                      mc   = 0.2225e+3L,
                      mt   = 70.5188e+3L,
                      md   = 1.0773L,
                      ms   = 20.4323L,
                      mb   = 0.9321e+3L,
                      me   = 0.4413L,
                      mmu  = 93.116L,
                      mtau = 1.1609L,
                      Vus  = 0.225341L,
                      Vub  = 0.00351001L,
                      Vcb  = 0.0411997L;

    A << mu*mu, mc*mc,   mt*mt,
         md*md, ms*ms,   mb*mb,
         Vus,   Vub,     Vcb,
         me*me, mmu*mmu, mtau*mtau;

    // B returns relative error in units of 1%
    for (int i = 0; i <= 11; i++)
        B(i) = (A(i) - tmp(i))/(0.01L * A(i));
    
    return B;
}

long double
chi2 (long double y1, long double y2, long double y3,
        long double z1, long double z2, long double z3,
        long double w1, long double w2, long double w3,
        long double u1, long double u2,
        long double vu, long double vd)
{
    ChiVecLd B = vector_B (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu, vd);

    long double X2 = B.transpose() * B;
    return X2/(12.0L);
}

//////////////////////////////////////////////////////////////////
/* Possible minimization algorith should come here, but it's ok */
//////////////////////////////////////////////////////////////////

parameterVectorLd
gradient_of_chi2 (long double y1, long double dy1,
        long double y2, long double dy2,
        long double y3, long double dy3,
        long double z1, long double dz1,
        long double z2, long double dz2,
        long double z3, long double dz3,
        long double w1, long double dw1,
        long double w2, long double dw2,
        long double w3, long double dw3,
        long double u1, long double du1,
        long double u2, long double du2,
        long double vu, long double dvu,
        long double vd, long double dvd)
{
    parameterVectorLd grad;

    long double X2 = chi2 (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu, vd);
    grad(0)  = (chi2(y1 + dy1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu, vd) - X2)/dy1;
    grad(1)  = (chi2(y1, y2 + dy2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu, vd) - X2)/dy2;
    grad(2)  = (chi2(y1, y2, y3 + dy3, z1, z2, z3, w1, w2, w3, u1, u2, vu, vd) - X2)/dy3;
    grad(3)  = (chi2(y1, y2, y3, z1 + dz1, z2, z3, w1, w2, w3, u1, u2, vu, vd) - X2)/dz1; 
    grad(4)  = (chi2(y1, y2, y3, z1, z2 + dz2, z3, w1, w2, w3, u1, u2, vu, vd) - X2)/dz2; 
    grad(5)  = (chi2(y1, y2, y3, z1, z2, z3 + dz3, w1, w2, w3, u1, u2, vu, vd) - X2)/dz3; 
    grad(6)  = (chi2(y1, y2, y3, z1, z2, z3, w1 + dw1, w2, w3, u1, u2, vu, vd) - X2)/dw1; 
    grad(7)  = (chi2(y1, y2, y3, z1, z2, z3, w1, w2 + dw2, w3, u1, u2, vu, vd) - X2)/dw2; 
    grad(8)  = (chi2(y1, y2, y3, z1, z2, z3, w1, w2, w3 + dw3, u1, u2, vu, vd) - X2)/dw3; 
    grad(9)  = (chi2(y1, y2, y3, z1, z2, z3, w1, w2, w3, u1 + du1, u2, vu, vd) - X2)/du1; 
    grad(10) = (chi2(y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2 + du2, vu, vd) - X2)/du2; 
    grad(11) = (chi2(y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu + dvu, vd) - X2)/dvu; 
    grad(12) = (chi2(y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu, vd + dvd) - X2)/dvd; 

    return grad;
}



int main (void)
{
    Vector3Ld v;
    v << 1.0L,2.0L,3.0L;

    std::cout << "v*v =\n" << v.transpose() * v << std::endl;

    exit(EXIT_SUCCESS);
}
