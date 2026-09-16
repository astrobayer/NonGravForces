/* ngfHwm14Load.c -- load HWM14's tables from a folder, once, without letting the model STOP.
 *
 * Shared by mex/ngfHwm14Gateway.c and mex/ngfRhs.c so the rule lives in one place:
 * HWM14's FINDANDOPEN looks in the current directory and executes STOP when a file is
 * missing, so the three files are checked first, the folder is made the PROCESS current
 * directory for the one evaluation that loads every table (quiet, disturbed and the QD
 * coordinates), and the previous directory is restored before returning.
 *
 * Returns 0 on success, 1 a file is missing, 2 the directory could not be read or entered,
 * 3 it could not be restored (the tables are loaded all the same).
 *
 * Created 2026-09-16 by Claude Opus 5 (Anthropic).
 */
#ifndef _WIN32
#define _POSIX_C_SOURCE 200112L   /* chdir and getcwd under -std=c99 on glibc */
#endif
#include <stdio.h>
#ifdef _WIN32
#include <direct.h>
#define NGF_CHDIR _chdir
#define NGF_GETCWD _getcwd
#else
#include <unistd.h>
#define NGF_CHDIR chdir
#define NGF_GETCWD getcwd
#endif

extern void ngf_hwm14_calc(int iyd, double sec, double alt, double glat, double glon, double ap2, double *w);

int ngfHwm14Load(const char *dir, char *missing, int missingSize)
{
    static const char *files[3] = {"hwm123114.bin", "dwm07b104i.dat", "gd2qd.dat"};
    char here[1024], full[1400];
    double w[2];
    int i;
    FILE *f;
    for (i = 0; i < 3; i++) {
        snprintf(full, sizeof full, "%s/%s", dir, files[i]);
        f = fopen(full, "rb");
        if (f == NULL) {
            if (missing != NULL && missingSize > 0) {
                snprintf(missing, (size_t)missingSize, "%s", files[i]);
            }
            return 1;
        }
        fclose(f);
    }
    if (NGF_GETCWD(here, sizeof here) == NULL) {
        return 2;
    }
    if (NGF_CHDIR(dir) != 0) {
        return 2;
    }
    ngf_hwm14_calc(95150, 43200.0, 250.0, 0.0, 0.0, 4.0, w);
    if (NGF_CHDIR(here) != 0) {
        return 3;
    }
    return 0;
}
