/* https://gcc.gnu.org/onlinedocs/gfortran/Non-Fortran-Main-Program.html */

extern void _gfortran_set_args(int argc, char *argv[]);
extern void _gfortran_set_options(int num, int options[]);
extern void _gfortran_set_fpe(int val);
#ifndef WIN64
#define _GNU_SOURCE // for fedisableexcept
#include <fenv.h>
#endif

int gfortran_init_() {
    int opts[] = {68, 511, 1, 1, 1, 63};
    char *argv[] = {"python", "-m", "mfixgui.pymfix"};
    _gfortran_set_args(3, argv);
    _gfortran_set_options(6, opts);
//  _gfortran_set_fpe(13); // default:  zero, overflow, invalid
    return 57;
}

int gfortran_set_fpe_(int *flags) {
    /*    IEEE exceptions.  Possible values are (bitwise
     *    or-ed) zero (0, default) no trapping,
     *    ‘GFC_FPE_INVALID’ (1), ‘GFC_FPE_DENORMAL’ (2),
     *    ‘GFC_FPE_ZERO’ (4), ‘GFC_FPE_OVERFLOW’ (8),
     *    ‘GFC_FPE_UNDERFLOW’ (16), and ‘GFC_FPE_INEXACT’
     *     (32).
     */
#if !defined(WIN64) && !defined(__APPLE__)
/* fedisableexcept does not exist on Windows
 * which means pausing, disabling FPE and restarting
 * may not work correctly */
    fedisableexcept(FE_ALL_EXCEPT);
#endif
    _gfortran_set_fpe(*flags);
    return 58;
}

int disable_backtrace_() {
    int opts[] = {68, 511, 0, 0, 0};
    _gfortran_set_options(5, opts);
    return 0;
}
