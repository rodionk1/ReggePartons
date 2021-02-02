/* qtsvd.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

//#include "f2c.h"

/* Subroutine */ int qtsvd_c(double *rint, double *x, int *n,
        double *f)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static double a, b;
    static int i__;
    static double a1, a2, x1, f12;
 
    /* Parameter adjustments */
    --f;
    --x;
    --rint;

    /* Function Body */
    rint[1] = 0.;
    x1 = (x[2] - x[1]) / 2.;
    f12 = f[1] * ((x1 - x[2]) * (x1 - x[3]) * (x1 - x[4])) / ((x[1] - x[2]) * 
            (x[1] - x[3]) * (x[1] - x[4])) + f[2] * ((x1 - x[1]) * (x1 - x[3])
             * (x1 - x[4])) / ((x[2] - x[1]) * (x[2] - x[3]) * (x[2] - x[4])) 
            + f[3] * ((x1 - x[1]) * (x1 - x[2]) * (x1 - x[4])) /
            ((x[3] - x[1]) * (x[3] - x[2]) * (x[3] - x[4])) + f[4] *
            ((x1 - x[1]) * (x1 - x[2]) * (x1 - x[3])) /
            ((x[4] - x[1]) * (x[4] - x[2]) * (x[4] - x[3]));
    rint[2] = (f[1] + f12 * 4. + f[2]) * (x[2] - x[1]) / 6.;
    a1 = x[2] - x[1];
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
        a = x[i__] - x[i__ - 2];
        b = a * a;
        a2 = x[i__] - x[i__ - 1];
        rint[i__] = rint[i__ - 2] + f[i__ - 2] * (-b / 6. + a * a1 / 2.) /
                a1 + f[i__ - 1] * (a * b / 6.) / (a2 * a1) + f[i__] *
                (b / 3. - a * a1 / 2.) / a2;
        a1 = a2;
/* L1: */
    }
    return 0;
} /* qtsvd_c */
