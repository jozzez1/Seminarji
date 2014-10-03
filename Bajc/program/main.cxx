#include <iostream>
#include <cmath>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

typedef Eigen::Matrix<long double, 3, 1> Vector3Ld;
typedef Eigen::Matrix<long double, 3, 3> Matrix3Ld;
typedef Eigen::Matrix<long double, 12, 1> ChiVecLd;
typedef Eigen::Matrix<long double, 13, 1> parameterVectorLd;
typedef Eigen::Matrix<long double, 12, 13> parameterGradientLd;
typedef Eigen::Matrix<long double, 13, 13> HessianMatrixLd;

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
vector_B (parameterVectorLd a)
{
    ChiVecLd B;
    Matrix3Ld MU2, MD2, ME2, LU, LD, VCKM;
    Vector3Ld eig;

    long double y1 = a(0),
                y2 = a(1),
                y3 = a(2),
                z1 = a(3),
                z2 = a(4),
                z3 = a(5),
                w1 = a(6),
                w2 = a(7),
                w3 = a(8),
                u1 = a(9),
                u2 = a(10),
                vu = a(11),
                vd = a(12);

    MU2 = up_mass_matrix_squared (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vu);
    MD2 = down_mass_matrix_squared (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, vd);
    ME2 = electron_mass_matrix_squared (y1, y2, y3, z1, z2, z3, w1, w2, w3, u1, u2, (-3)*vd);
    Eigen::SelfAdjointEigenSolver <Matrix3Ld> es;

    // B accepts theoretical values
    es.compute (MU2, Eigen::ComputeEigenvectors);
    eig  = es.eigenvalues();
    LU   = es.eigenvectors();
    assert (fabsl(LU.determinant() - 1.0L) < 1e-6L);
    B(0) = eig(0);
    B(1) = eig(1);
    B(2) = eig(2);

    es.compute (MD2, Eigen::ComputeEigenvectors);
    eig  = es.eigenvalues();
    LD   = es.eigenvectors();
    assert (fabsl(LD.determinant() - 1.0L) < 1e-6L);
    B(3) = eig(0);
    B(4) = eig(1);
    B(5) = eig(2);

    VCKM = LU.transpose() * LD;
    B(6) = VCKM(1,2);
    B(7) = VCKM(1,3);
    B(8) = VCKM(2,3);

    es.compute (ME2, Eigen::ComputeEigenvectors);
    eig  = es.eigenvalues();

    B(9) = eig(0);
    B(10)= eig(1);
    B(11)= eig(2);
   
    return B;
}

ChiVecLd
vector_A (void)
{
    ChiVecLd A;

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

    return A;
}

long double
chi2 (parameterVectorLd a)
{
    ChiVecLd B, A, tmp;
    tmp = vector_B (a);
    A = vector_A ();

    // B now returns relative error in units of 1%
    for (int i = 0; i <= 11; i++)
        B(i) = (A(i) - tmp(i))/(0.01L * A(i));
 

    long double X2 = B.transpose() * B;
    return X2/(12.0L);
}

//////////////////////////////////////////
/* Minimization via Levenberg-Marquardt */
//////////////////////////////////////////

parameterGradientLd
gradient_of_B (parameterVectorLd a, parameterVectorLd da, long double lambda)
{
    ChiVecLd B, B0;
    parameterVectorLd ei;
    parameterGradientLd gradient_tensor;
    long double h = 0.05L * lambda,
                H;

    B0 = vector_B (a);

    for (int i = 0; i <= 12; i++)
    {
        H = h*da(i);
        B = (1.0L/H) * (vector_B (a + H*ei.Unit(i)) - B0);
        for (int j = 0; j <= 11; j++)
            gradient_tensor(j,i) = B(j);
    }

    return gradient_tensor;
}

parameterVectorLd
gradient_of_chi2 (parameterVectorLd a, parameterVectorLd da, long double lambda)
{
    ChiVecLd A, B, tmp1, tmp2;
    parameterVectorLd grad;
    parameterGradientLd gradient_tensor;

    gradient_tensor = gradient_of_B (a, da, lambda);
    A    = vector_A ();
    tmp1 = vector_B (a);
    tmp2 = 0.0001L * A.array() * A.array();
    for (int i = 0; i <= 11; i++)
        B(i) = (A(i) - tmp1(i))/tmp2(i);

    grad = (-0.5L) * gradient_tensor.transpose() * B;

    return grad;
}

HessianMatrixLd
modified_Hessian_of_chi2 (parameterVectorLd a, parameterVectorLd da, long double lambda)
{
    parameterGradientLd gradient_tensor;
    HessianMatrixLd H;
    gradient_tensor = gradient_of_B (a, da, lambda);
    H = (0.5L) * gradient_tensor.transpose() * gradient_tensor;

    for (int i = 0; i <= 12; i++)
        H(i,i) += lambda;

    return H;
}

parameterVectorLd
levenberg_marquardt_step (parameterVectorLd a, parameterVectorLd da,
        long double * lambda, long double * X2)
{
    HessianMatrixLd H = modified_Hessian_of_chi2 (a, da, *lambda);
    parameterVectorLd grad_of_chi2 = gradient_of_chi2 (a, da, *lambda),
                      da_new = H.fullPivLu().solve (grad_of_chi2);
    long double X2n = chi2 (a + da_new);

    while (X2n >= *X2)
    {
        *lambda *= 10.0L;
        H = modified_Hessian_of_chi2 (a, da, *lambda);
        grad_of_chi2 = gradient_of_chi2 (a, da, *lambda);
        da_new = H.fullPivLu().solve(grad_of_chi2);
        X2n = chi2 (a + da_new);
    }

    *X2 = X2n;

    return da_new;
}

parameterVectorLd
levenberg_marquardt (parameterVectorLd a, parameterVectorLd da,
        long double * lambda, long double * X2, int maxIter)
{
    int iter = 0;
    *X2 = chi2 (a);
    parameterVectorLd da_new;

    do
    {
        da_new = levenberg_marquardt_step (a, da, lambda, X2);
        a += da_new;
        da = da_new;
        iter++;
    } while (*X2 > 1 && iter < maxIter);

    return a;
}

HessianMatrixLd
covariance_matrix (HessianMatrixLd H, long double lambda)
{
    for (int i = 0; i <= 12; i++)
        H(i,i) -= lambda;

    return H.inverse();
}

///////////////////
/* Main function */
///////////////////

int main (void)
{
    Vector3Ld v;
    v << 1.0L,2.0L,3.0L;

    std::cout << "v*v =\n" << v.transpose() * v << std::endl;

    ChiVecLd ei;
    for (int i = 0; i <= 5; i++)
        std::cout << "Some unit vector:" << ei.Unit(i).transpose() << std::endl;

    exit (EXIT_SUCCESS);
}
