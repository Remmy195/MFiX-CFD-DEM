// Author:  Charles G Waldman 2024-07-31
// Replacement for Fortran "SYSTEM_CLOCK" which suffers from overflow and other issues
// See https://mfix.netl.doe.gov/gitlab/develop/mfix/-/issues/48

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

double etime() {
    struct timeval tv;
    static double t0 = 0.0;

    if (t0 == 0.0) { // first call
        if (gettimeofday(&tv, NULL)) {
                perror("gettimeofday");
                return -1;
        }
        t0 = tv.tv_sec + tv.tv_usec/1000000.0;
        return 0.0;
    }
    if (gettimeofday(&tv, NULL)) {
        perror("gettimeofday");
        return -1;
    }
    return tv.tv_sec + tv.tv_usec/1000000.0 - t0;
}

double _gfortran_cpu_time_8(double x) {
    fprintf(stderr,"Please do not use CPU_TIME.  Use the functions in time_mod instead.\n");
    fprintf(stderr,"See https://mfix.netl.doe.gov/gitlab/develop/mfix/-/issues/48\n");
    exit(1);
}
