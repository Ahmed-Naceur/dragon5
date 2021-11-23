
/*--------------------------------*/
/* Python3-Cle2000 bindings       */
/* author: A. Hebert (03/07/2020) */
/*--------------------------------*/

/*
Copyright (C) 2020 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/
#include <stdio.h>
#include <Python.h>
#include <structmember.h>
#include <setjmp.h>
#include "pylcm.h"
jmp_buf buf;

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    lifo *stack;           /* internal structure */
    int_32 impx_cle2000;   /* print flag */
    char name_cle2000[13]; /* procedure name */
} cle2000object;

static char AbortString[132];
static PyObject *PyCle2000Error = NULL;
static PyTypeObject Cle2000Type;     /* shared type-descriptor */

void xabort_c(char *msg){
  fflush(stdout); fflush(stderr);
  PyErr_SetString(PyCle2000Error, msg);
  longjmp(buf, 1);
}

static PyObject *cle2000_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
   char *nomsub="cle2000_new";
   PyObject *stacklist, *item=NULL;
   char *procname;
   long impx = 0;
  
   /* Parse arguments */
   static char *kwlist[] = {"procname", "stack", "impx", NULL};
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO|l", kwlist,
                                      &procname, &stacklist, &impx)) {
     return NULL;
   }

   cle2000object *self = (cle2000object*)PyObject_New(cle2000object, type);
   self->impx_cle2000 = impx;

   if (strlen(procname) > 12) {
     sprintf(AbortString,"%s: character procname overflow",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }
   strcpy(self->name_cle2000,procname);

   PyObject *module = PyImport_ImportModule("lifo");
   if (module == NULL) {
     sprintf(AbortString,"%s: module lifo should be installed first",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }
   PyObject *moduleDict = PyModule_GetDict(module);
   PyObject *protocolClass = PyDict_GetItemString(moduleDict, "new");
   if (protocolClass == NULL) {
     sprintf(AbortString,"%s: unable to find class new in module lifo",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }

   /* Set the Lifo stack */
   lifo *my_lifo;
   lifo_node *my_node;
   if (PyObject_IsInstance(stacklist, protocolClass)) { /* stacklist is a lifo stack object */
     lifoobject *lifo_object = (lifoobject *)stacklist;
     my_lifo = lifo_object->stack;
   } else { /* stacklist is a unique readonly PyObject */
     cleopn(&my_lifo);
     item = stacklist;
     my_node = (lifo_node *) malloc(sizeof(lifo_node));
     strcpy(my_node->name, "inpu_val0001"); 
     if (PyObject_IsInstance(item, (PyObject *)&PyLong_Type)) {
       long my_int = (long)PyLong_AsLong(item);
       if (impx > 0) printf("%s: unique item is an integer --> %ld\n", nomsub, my_int);
       my_node->type = 11; my_node->value.ival = my_int;
     } else if (PyObject_IsInstance(item, (PyObject *)&PyFloat_Type)) {
       double my_double = (double)PyFloat_AsDouble(item);
       if (impx > 0) printf("%s: unique item is a floating point number --> %e\n", nomsub, my_double);
       my_node->type = 14; my_node->value.dval = my_double;
     } else if (PyObject_IsInstance(item, (PyObject *)&PyBool_Type)) {
       long my_int = 1;
       if (item == Py_False) my_int = 0;
       if (impx > 0) printf("%s: unique item is a boolean number --> %ld\n", nomsub, my_int);
       my_node->type = 15; my_node->value.ival = my_int;
     } else if (PyObject_IsInstance(item, (PyObject *)&PyUnicode_Type)) {
       int len = PyUnicode_GET_LENGTH(item);
       if (len > 72) {
         sprintf(AbortString,"%s: character data overflow",nomsub);
         PyErr_SetString(PyCle2000Error, AbortString);
         return NULL;
       }
       int kind = PyUnicode_KIND(item);
       void *data = PyUnicode_DATA(item);
       char* my_string;
       my_string = (char *)malloc(len+1);
       int n;
       for (n=0; n<len; n++) {
         my_string[n] = (Py_UCS4)PyUnicode_READ(kind, data, n);
       }
       my_string[len] = '\0';
       if (impx > 0) printf("%s: unique item is a string --> %s\n", nomsub, my_string);
       my_node->type = 13;
       strcpy(my_node->value.hval, my_string);
       free(my_string);
     } else {
       sprintf(AbortString,"%s: use a list or a Lifo object as second argument to cle2000.new", nomsub);
       PyErr_SetString(PyCle2000Error, AbortString);
       return NULL;
     }
     clepush(&my_lifo, my_node);
   }
   self->stack = my_lifo;
   return (PyObject *)self;
}

static void cle2000_dealloc(cle2000object *self)
{
   char *nomsub="cle2000_dealloc";
   if (self->impx_cle2000 > 0) printf("%s: desallocate a cle2000 object\n", nomsub);
}

static PyObject *cle2000_exec(cle2000object *self) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="cle2000_exec";
     lifo *my_lifo;
     int ipos;
     int_32 ganmod(char *cmodul, int_32 nentry, char (*hentry)[13], int_32 *ientry,
            int_32 *jentry, lcm **kentry, char (*hparam)[73]);
     int_32 trimod(char *cmodul, int_32 nentry, char (*hentry)[13], int_32 *ientry,
            int_32 *jentry, lcm **kentry, char (*hparam)[73]);
     int_32 dramod(char *cmodul, int_32 nentry, char (*hentry)[13], int_32 *ientry,
            int_32 *jentry, lcm **kentry, char (*hparam)[73]);
     int_32 donmod(char *cmodul, int_32 nentry, char (*hentry)[13], int_32 *ientry,
            int_32 *jentry, lcm **kentry, char (*hparam)[73]);

     /* Assign the code entry point */
#if (defined __trivac__)
     int_32 (*dummod)() = &trimod;
     if (self->impx_cle2000 > 0) printf("%s: install trivac library\n", nomsub);
#elif (defined __dragon__)
     int_32 (*dummod)() = &dramod;
     if (self->impx_cle2000 > 0) printf("%s: install dragon library\n", nomsub);
#elif (defined __donjon__)
     int_32 (*dummod)() = &donmod;
     if (self->impx_cle2000 > 0) printf("%s: install donjon library\n", nomsub);
#else
     int_32 (*dummod)() = &ganmod;
     if (self->impx_cle2000 > 0) printf("%s: install ganlib library\n", nomsub);
#endif

     my_lifo = self->stack;
     int nstack = (int)my_lifo->nup;
     /* close the LCM objects */
     for (ipos=0; ipos<nstack; ++ipos) {
       lifo_node *my_node = clepos(&my_lifo, ipos);
         if ((my_node->type == 3) || (my_node->type == 4)) {
           if (my_node->access > 0) lcmcl_c(&(my_node->value.mylcm), 1);
       }
     }

     /* call the parametrized procedure */
     int_32 ier, ilevel = 1;
     ier = cle2000_c(ilevel, dummod, self->name_cle2000, self->impx_cle2000-1, my_lifo);
     if (ier != 0) {
       sprintf(AbortString,"%s: cle2000 failure (%s.c2m). ier=%d", nomsub, self->name_cle2000, ier);
       PyErr_SetString(PyCle2000Error, AbortString);
       return NULL;
     }

     /* reopen the LCM objects */
     for (ipos=0; ipos<nstack; ++ipos) {
       lifo_node *my_node = clepos(&my_lifo, ipos);
       if ((my_node->type == 3) || (my_node->type == 4)) {
         if (my_node->access == 0) my_node->access = 1;
         lcmop_c(&(my_node->value.mylcm), my_node->OSname, my_node->access, my_node->type-2, 0);
       }
     }
     if (self->impx_cle2000 > 0) clelib(&my_lifo);
     clecls(&my_lifo);
     fflush(stdout);
     return Py_None;
   }
}

static PyObject *cle2000_getLifo(cle2000object *self) {
   char *nomsub="getLifo";

   PyObject *module = PyImport_ImportModule("lifo");
   if (module == NULL) {
     sprintf(AbortString,"%s: module lifo should be installed first",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }
   PyObject *moduleDict = PyModule_GetDict(module);
   PyTypeObject *protocolClass = (PyTypeObject *)PyDict_GetItemString(moduleDict, "new");
   if (protocolClass == NULL) {
     sprintf(AbortString,"%s: unable to find class new in module lifo",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }

   lifoobject *rv = (lifoobject*)PyObject_New(lifoobject, protocolClass);
   rv->stack = self->stack;
   rv->impx_lifo = self->impx_cle2000;
   Py_INCREF(rv);
   return (PyObject *)rv;
}

static PyObject *cle2000_putLifo(cle2000object *self, PyObject *args) {
   char *nomsub="putLifo";
   PyObject *item=NULL;
  
   /* Parse arguments */
   if(!PyArg_ParseTuple(args, "O", &item)) {
       return NULL;
   }

   PyObject *module = PyImport_ImportModule("lifo");
   if (module == NULL) {
     sprintf(AbortString,"%s: module lifo should be installed first",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }
   PyObject *moduleDict = PyModule_GetDict(module);
   PyObject *protocolClass = PyDict_GetItemString(moduleDict, "new");
   if (protocolClass == NULL) {
     sprintf(AbortString,"%s: unable to find class new in module lifo",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }
   if (!PyObject_IsInstance(item, protocolClass)) {
     sprintf(AbortString,"%s: the argument is not a lifo instance",nomsub);
     PyErr_SetString(PyCle2000Error, AbortString);
     return NULL;
   }

   lifoobject *rv = (lifoobject*)item;
   Py_DECREF(self->stack);
   self->stack = rv->stack;
   return Py_None;
}

static PyMethodDef cle2000_methods[] = {
    {"exec", (PyCFunction) cle2000_exec, METH_NOARGS, "Execute the cle-2000 procedure"},
    {"getLifo", (PyCFunction) cle2000_getLifo, METH_NOARGS, "Recover the embedded Lifo stack"},
    {"putLifo", (PyCFunction) cle2000_putLifo, METH_VARARGS, "Replace the embedded Lifo stack"},
    {NULL, NULL, 0, NULL}
};

static PyMemberDef cle2000_members[] = {
    {"_impx", T_INT, offsetof(cle2000object, impx_cle2000), 0, "print index"},
    {NULL}  /* Sentinel */
};

static PyTypeObject Cle2000Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "cle2000.new",             /*tp_name*/
    sizeof(cle2000object),     /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)cle2000_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,  /*tp_flags*/
    "Custom objects",          /*tp_doc*/
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    cle2000_methods,           /* tp_methods */
    cle2000_members,           /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    cle2000_new,               /* tp_new */
};

// Module definition
static struct PyModuleDef cle2000module = { 
    PyModuleDef_HEAD_INIT,
    "cle2000",
    "A Python module that provides support to the Ganlib API",
    -1, 
    NULL, NULL, NULL, NULL, NULL
};

// Module initialization
PyMODINIT_FUNC PyInit_cle2000(void) {
    PyObject* m;

    Py_Initialize();
    if (PyType_Ready(&Cle2000Type) < 0) return NULL;

    m = PyModule_Create(&cle2000module);
    if (m == NULL) return NULL;

    PyModule_AddObject(m, "new", (PyObject *)&Cle2000Type);

    /* Initialize new exception object */
    PyCle2000Error = PyErr_NewException("cle2000.PyCle2000Error", NULL, NULL);

    /* Add exception object to your module */
    PyModule_AddObject(m, "PyCle2000Error", PyCle2000Error);
    return m;
}
