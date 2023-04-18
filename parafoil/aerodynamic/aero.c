#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979323846

void print_matrix(double *M, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%6.3f", *(M + i * cols + j));
        }
        printf("\n");
    }
}

void v_diff(double *res, double *v1, double *v2)
{
    res[0] = v1[0] - v2[0];
    res[1] = v1[1] - v2[1];
    res[2] = v1[2] - v2[2];
}

void v_sum(double *res, double *v1, double *v2)
{
    res[0] = v1[0] + v2[0];
    res[1] = v1[1] + v2[1];
    res[2] = v1[2] + v2[2];
}

double v_norm(double *v)
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void v_cross(double *res, double *v1, double *v2)
{
    res[0] = v1[1] * v2[2] - v1[2] * v2[1];
    res[1] = v1[2] * v2[0] - v1[0] * v2[2];
    res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void v_mult(double *res, double *v, double num)
{
    for (size_t i = 0; i < 3; i++)
    {
        res[i] = v[i] * num;
    }
}

double v_dot(double *v1, double *v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void vortxl(double *u, double *X1, double *X2, double *XP, double gamma)
{
    // x1 = X1[0]
    // y1 = X1[1]
    // z1 = X1[2]
    // x2 = X2[0]
    // y2 = X2[1]
    // z2 = X2[2]
    // xp = XP[0]
    // yp = XP[1]
    // zp = XP[2]

    double r0[3];
    v_diff(r0, X2, X1);
    double r1[3];
    v_diff(r1, XP, X1);
    double r2[3];
    v_diff(r2, XP, X2);

    double norm_r1 = v_norm(r1);
    double norm_r2 = v_norm(r2);

    double r1xr2[3];
    v_cross(r1xr2, r1, r2);

    double norm_r1xr2 = v_norm(r1xr2);

    double inv_r1xr2 = 0;
    if (norm_r1xr2 != 0)
    {
        inv_r1xr2 = 1.0 / norm_r1xr2;
    }

    double inv_r1 = 0;
    if (norm_r1 != 0)
    {
        inv_r1 = 1.0 / norm_r1;
    }

    double inv_r2 = 0;
    if (norm_r2 != 0)
    {
        inv_r2 = 1.0 / norm_r2;
    }

    double a[3];
    v_mult(a, r0, inv_r1xr2);
    double rr1[3];
    v_mult(rr1, r1, inv_r1);
    double rr2[3];
    v_mult(rr2, r2, inv_r2);
    double b[3];
    v_diff(b, rr1, rr2);
    double c = v_dot(a, b);

    // Записываем результат в переменную u
    v_mult(u, r1xr2, gamma * (0.25 / PI) * c * inv_r1xr2);
}

void biot_savart(int N, double b, double *A, double *B, double *normals, double *xctrl, double *xbound, double *coord, double *Vinf)
{
    double one = 1.0;
    double xa[3] = {0, 0, 0};
    double xb[3] = {0, 0, 0};
    double xc[3] = {0, 0, 0};
    double xd[3] = {0, 0, 0};
    double xp[3] = {0, 0, 0};
    double nunit[3] = {0, 0, 0};
    double coord_diff = 0;
    double uind1[3] = {0, 0, 0};
    double uind2[3] = {0, 0, 0};
    double uind3[3] = {0, 0, 0};
    double uindt[3] = {0, 0, 0};
    double uinf[3] = {0, 0, 0};

    for (size_t j = 0; j < N; j++)
    {
        for (size_t r = 0; r < 3; r++)
        {
            // xp[r] = xctrl[r, j];
            xp[r] = *(xctrl + r * N + j);
            // nunit[r] = normals[r, j];
            nunit[r] = *(normals + r * N + j);
            uinf[r] = -*(Vinf + j * 3 + r);
        }
        // B[j, 0] = -v_dot(uinf, nunit);
        *(B + j * 1 + 0) = -v_dot(uinf, nunit);

        for (size_t k = 0; k < N; k++)
        {
            // xb[0] = xbound[0, k];
            xb[0] = *(xbound + 0 * N + k);
            xa[0] = xb[0] + 20 * b;
            xc[0] = xb[0];
            xd[0] = xa[0];

            // xb[1] = coord[1, k];
            xb[1] = *(coord + 1 * (N + 1) + k);
            xa[1] = xb[1];
            // xc[1] = coord[1, k + 1];
            xc[1] = *(coord + 1 * (N + 1) + k + 1);
            xd[1] = xc[1];

            // (coord[2, k + 1] - coord[2, k]) / 2
            coord_diff = (*(coord + 2 * (N + 1) + k + 1) - *(coord + 2 * (N + 1) + k)) / 2;
            // xbound[2,k] - coord_diff
            xb[2] = *(xbound + 2 * N + k) - coord_diff;
            xa[2] = xb[2];
            // xbound[2,k] + coord_diff
            xc[2] = *(xbound + 2 * N + k) + coord_diff;
            xd[2] = xc[2];

            // First trailing vortex
            vortxl(uind1, xa, xb, xp, one);
            // Bounded vortex
            vortxl(uind2, xb, xc, xp, one);
            // Second trailing vortex
            vortxl(uind3, xc, xd, xp, one);

            v_sum(uindt, uind1, uind2);
            v_sum(uindt, uindt, uind3);
            // A[j, k] = v_dot(uindt, nunit)
            *(A + j * N + k) = v_dot(uindt, nunit);
        }
    }
}

void kutta_joukowsky(int N, double b, double Rho, double *local_force, double *xbound, double *coord, double *Circ, double *Vinf2)
{
    double xa[3] = {0, 0, 0};
    double xb[3] = {0, 0, 0};
    double xc[3] = {0, 0, 0};
    double xd[3] = {0, 0, 0};
    double xp[3] = {0, 0, 0};
    double coord_diff = 0;
    double Wk[3] = {0, 0, 0};
    double uinf[3] = {0, 0, 0};
    double wkind1[3] = {0, 0, 0};
    double wkind3[3] = {0, 0, 0};
    double v_i[3] = {0, 0, 0};
    double d_g[3] = {0, 0, 0};
    double f[3] = {0, 0, 0};

    for (size_t j = 0; j < N; j++)
    {
        for (size_t r = 0; r < 3; r++)
        {
            // xp[r] = xbound[r, j];
            xp[r] = *(xbound + r * N + j);
            Wk[r] = 0;
            uinf[r] = -*(Vinf2 + j * 3 + r);
        }
        for (size_t k = 0; k < N; k++)
        {
            // xb[0] = xbound[0, k];
            xb[0] = *(xbound + 0 * N + k);
            xa[0] = xb[0] + 20 * b;
            xc[0] = xb[0];
            xd[0] = xa[0];

            // xb[1] = coord[1, k];
            xb[1] = *(coord + 1 * (N + 1) + k);
            xa[1] = xb[1];
            // xc[1] = coord[1, k + 1];
            xc[1] = *(coord + 1 * (N + 1) + k + 1);
            xd[1] = xc[1];

            // (coord[2, k + 1] - coord[2, k]) / 2
            coord_diff = (*(coord + 2 * (N + 1) + k + 1) - *(coord + 2 * (N + 1) + k)) / 2;
            // xbound[2,k] - coord_diff
            xb[2] = *(xbound + 2 * N + k) - coord_diff;
            xa[2] = xb[2];
            // xbound[2,k] + coord_diff
            xc[2] = *(xbound + 2 * N + k) + coord_diff;
            xd[2] = xc[2];

            // First trailing vortex
            vortxl(wkind1, xa, xb, xp, *(Circ + k * 1 + 0));
            // Second trailing vortex
            vortxl(wkind3, xc, xd, xp, *(Circ + k * 1 + 0));

            // Down Wash
            v_sum(Wk, Wk, wkind1);
            v_sum(Wk, Wk, wkind3);
        }

        // Induced + Kinematic velocity
        v_sum(v_i, Wk, uinf);

        // Bound vortex
        for (size_t r = 0; r < 3; r++)
        {
            // xb = coord[:, j]
            xb[r] = *(coord + r * (N + 1) + j);
            // xc = coord[:, j+1]
            xc[r] = *(coord + r * (N + 1) + j + 1);
        }
        // if (1)
        // {
        //     printf("xb = ");
        //     print_matrix(xb, 1, 3);
        // }
        // Gamma is per unit lenght
        v_diff(d_g, xc, xb);
        v_mult(d_g, d_g, *(Circ + j * 1 + 0));
        // KJ
        v_cross(f, v_i, d_g);
        v_mult(f, f, Rho);
        for (size_t r = 0; r < 3; r++)
        {
            // local_force[:, j]
            *(local_force + r * N + j) = f[r];
        }
    }
}