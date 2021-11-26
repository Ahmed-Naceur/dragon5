
/*--------------------------------*/
/* Python3-LCM bindings           */
/* author: A. Hebert (03/07/2020) */
/*--------------------------------*/

/*
Copyright (C) 2020 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/
#include "lcm.h"
#include "cle2000.h"

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    lifo *stack;         /* internal structure */
    int_32 impx_lifo;    /* print flag */
} lifoobject;

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    lcm *iplist;         /* lcm/xsm/file object handle */
    int_32 impx_lcm;     /* print flag */
    int_32 lrda_lcm;     /* da size */
    int_32 iact_lcm;     /* access mode */
    char type_lcm[13];   /* object type */
    char name_lcm[73];   /* object name */
} pylcmobject;
