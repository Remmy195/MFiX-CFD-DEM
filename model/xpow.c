/* xpow.c - replacement pow() function that detects integer and dyadic exponents
 *
 * build with:
 *
 * gcc -O3 -march=native -fPIC -c xpow.c
 * ld -shared -o xpow.so xpow.o -ldl -lm
 *
 * gcc -Wl,wrap=pow
 *
 * Charles G Waldman  2021-10-20
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>

double __real_pow(double x, double y);

// see https://lwn.net/Articles/691932/ for function multi-versioning
//  (we are currently not using this)

// optimized pow function to replace library function
double __wrap_pow(double x, double y) {
    double r, s, x0, y0;
    if (y == 0.0)
        return 1.0;
    if (y == 1.0)
        return x;
    if (y == 2.0)
        return x*x;
    if (x == 0.0)
        return 0.0;
    if (y == -1.0)
        return 1.0/x;
    if (y < 0.0)
        return __real_pow(x, y);
    int iy = y;
    x0 = x;
    y0 = y;
    r = 1.0;

    if (iy) {  /// __builtin_clz may yield garbage for 0, if not using LZCNT instruction.
        int zeros =  __builtin_clz(iy); // find highest bit set
        //printf("%d %d\n", iy, zeros);
        //fflush(stdout);
        if (zeros < 28)
            return __real_pow(x0, y0);
        y -= iy;
        switch(zeros) {
        case 28:
            if (iy & 1) r *= x;
            iy >>= 1;
            x *= x;
        case 29:
            if (iy & 1) r *= x;
            iy >>= 1;
            x *= x;
        case 30:
            if (iy & 1) r *= x;
            iy >>= 1;
            x *= x;
        case 31:
            if (iy & 1) r *= x;
        }
    }
    // handle dyadic fractional part
    if (y == 0)
        return r;
    s = sqrt(x0);
    if (y == 0.5)
        return s*r;
    else if (y == 0.25)
        return sqrt(s)*r;
    else if (y == 0.75)
        return s*sqrt(s)*r;
    // fall back to math library
    return __real_pow(x0, y0);
}

#ifdef BUILD_TEST
#include <stdio.h>
int main(int argc, char **argv)
{
    double x, y, p1, p2;
    for (x=-5; x<=5; x+=0.5) {
        for (y=-5; y<=30; y += 0.125) {
            p1 = __wrap_pow(x,y);
            p2 = pow(x,y);
            if (fabs((p1-p2)/p2) > 1e-15)
                printf("%f %f %.20f %.20f\n", x, y, p1, p2);
        }
    }
}
#endif
