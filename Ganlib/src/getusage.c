/*
   -----------------------------------------------------------------------
   Copyright (C) 2019 Ecole Polytechnique de Montreal

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   RECOVER USER ALLOCATED MEMORY USED

   -----------------------------------------------------------------------
 */

#include <stdio.h>
#ifndef MSDOS
#include <sys/time.h>
#include <sys/resource.h>

double getusage()
{
   struct rusage r_usage;
   getrusage(RUSAGE_SELF, &r_usage);
   double utime=(double) r_usage.ru_maxrss;
   return utime;
#else
#include <windows.h>
#include <psapi.h>

double getusage()
{
    HANDLE h_process = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS pmc;

    if(GetProcessMemoryInfo(h_process, &pmc, sizeof(pmc))) {
        // The *nix equivalent is expressed in KB, but Windows uses bytes
        SIZE_T mem_kb = pmc.WorkingSetSize / 1024;
        return (double)mem_kb;
    } else {
        return 0.0;
    }
}
#endif
}
