
/*--------------------------------*/
/* Python3-Lifo bindings          */
/* author: A. Hebert (03/07/2020) */
/*--------------------------------*/

/*
Copyright (C) 2020 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/
#include <Python.h>
#include <structmember.h>
#include "pylcm.h"

static char AbortString[132];
static PyObject *PyLifoError = NULL;
static PyTypeObject LifoType;     /* shared type-descriptor */

static PyObject *lifo_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
     int_32 impx = 0;
  
     /* Parse arguments */
     static char *kwlist[] = {"impx", NULL};
     if (!PyArg_ParseTupleAndKeywords(args, kwds, "|l", kwlist, &impx)) {
         return NULL;
     }

    lifoobject *self = (lifoobject*)PyObject_New(lifoobject, type);
    self->impx_lifo = impx;
    cleopn(&(self->stack));
    return (PyObject *)self;
}

static void lifo_dealloc(lifoobject *self)
{
   char *nomsub="lifo_dealloc";
   clecls(&(self->stack));
   if (self->impx_lifo > 0) printf("%s: desallocate a lifo object\n",nomsub);
}

static PyObject *lifo_lib(lifoobject *self) {
   fflush(stdout);
   clelib(&(self->stack));
   fflush(stdout);
   return Py_None;
}

static PyObject *lifo_push(lifoobject *self, PyObject *args) {
   char *nomsub="lifo_push";
   PyObject *item=NULL;
   lifo_node *my_node;
  
   /* Parse arguments */
   if(!PyArg_ParseTuple(args, "O", &item)) {
       return NULL;
   }

   PyObject *module = PyImport_ImportModule("lcm");
   if (module == NULL) {
     sprintf(AbortString,"%s: module lcm should be installed first",nomsub);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }
   PyObject *moduleDict = PyModule_GetDict(module);
   PyObject *protocolClass = PyDict_GetItemString(moduleDict, "new");
   if (protocolClass == NULL) {
     sprintf(AbortString,"%s: unable to find class new in module lcm",nomsub);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }

   lifo *my_lifo = self->stack;
   int ipos = my_lifo->nup;
   my_node = (lifo_node *) malloc(sizeof(lifo_node));
   sprintf(my_node->name,"inpu_val%04d", ipos+1);
   if (PyObject_IsInstance(item, (PyObject *)&PyLong_Type)) {
     long my_int = (long)PyLong_AsLong(item);
     if (self->impx_lifo > 1) printf("item %d is an integer --> %ld\n", ipos, my_int);
     my_node->type = 11; my_node->value.ival = my_int;
   } else if (PyObject_IsInstance(item, (PyObject *)&PyFloat_Type)) {
     double my_double = (double)PyFloat_AsDouble(item);
     if (self->impx_lifo > 1) printf("item %d is a double floating point number --> %e\n", ipos, my_double);
     my_node->type = 14; my_node->value.dval = my_double;
   } else if (PyObject_IsInstance(item, (PyObject *)&PyBool_Type)) {
     long my_int = 1;
     if (item == Py_False) my_int = -1;
     if (self->impx_lifo > 0) printf("%s: unique item is a boolean number --> %ld\n", nomsub, my_int);
     my_node->type = 15; my_node->value.ival = my_int;
   } else if (PyObject_IsInstance(item, (PyObject *)&PyUnicode_Type)) {
     int len = PyUnicode_GET_LENGTH(item);
     if (len > 72) {
       sprintf(AbortString,"%s: character data overflow",nomsub);
       PyErr_SetString(PyLifoError, AbortString);
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
     if (self->impx_lifo > 1) printf("item %d is a string --> %s\n", ipos, my_string);
     my_node->type = 13;
     strcpy(my_node->value.hval, my_string);
     free(my_string);
   } else if (PyObject_IsInstance(item, protocolClass)) {
     pylcmobject *lcm_object = (pylcmobject *)item;
     char *pytype = lcm_object->type_lcm;
     if (strcmp(pytype,"LCM") == 0) {
       my_node->type = 3;
     } else if (strcmp(pytype,"XSM") == 0) {
       my_node->type = 4;
     } else if (strcmp(pytype,"BINARY") == 0) {
       my_node->type = 5;
     } else if (strcmp(pytype,"ASCII") == 0) {
       my_node->type = 6;
     } else if (strcmp(pytype,"DA")== 0) {
       my_node->type = 7;
       my_node->lparam = lcm_object->lrda_lcm;
     } else {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLifoError, AbortString);
       return NULL;
     }
     if (self->impx_lifo > 1) printf("item %d is a lcm object named %s (type=%d)\n", ipos, lcm_object->name_lcm, my_node->type);
     strcpy(my_node->OSname, lcm_object->name_lcm);
     my_node->access = lcm_object->iact_lcm;
     my_node->value.mylcm = lcm_object->iplist;
   } else {
     /* Create an empty slot in the lifo */
     sprintf(my_node->name,"outp_val%04d", ipos+1);
     PyObject* objectsRepresentation = PyObject_Repr(item);
     int len = PyUnicode_GET_LENGTH(objectsRepresentation);
     int kind = PyUnicode_KIND(objectsRepresentation);
     void *data = PyUnicode_DATA(objectsRepresentation);
     char* my_string;
     my_string = (char *)malloc(len+1);
     int n;
     for (n=0; n<len; n++) {
       my_string[n] = (Py_UCS4)PyUnicode_READ(kind, data, n);
     }
     my_string[len] = '\0';
     if (self->impx_lifo > 1) printf("item %d is a %s\n", ipos, my_string);
     if (strcmp(my_string, "<class 'int'>") == 0) {
       my_node->type = -11;
     } else if (strcmp(my_string, "<class 'str'>") == 0) {
       my_node->type = -13;
     } else if (strcmp(my_string, "<class 'float'>") == 0) {
       my_node->type = -14;
     } else if (strcmp(my_string, "<class 'bool'>") == 0) {
       my_node->type = -15;
     }
     free(my_string);
   }
   clepush(&my_lifo, my_node);
   return Py_None;
}

static PyObject *lifo_pushEmpty(lifoobject *self, PyObject *args) {
   char *nomsub="lifo_pushEmpty";
   lifo_node *my_node;
   char *OSname = NULL;
   char *htype = NULL;
  
   /* Parse arguments */
   if(!PyArg_ParseTuple(args, "s|s", &OSname, &htype)) {
       return NULL;
   }
   
   lifo *my_lifo = self->stack;
   int ipos = my_lifo->nup;
   my_node = (lifo_node *) malloc(sizeof(lifo_node));
   if (self->impx_lifo > 1) printf("item %d is an empty %s object\n", ipos, htype);
   if (strlen(OSname) > 72) {
      sprintf(AbortString,"%s: character OSname overflow",nomsub);
      PyErr_SetString(PyLifoError, AbortString);
      return NULL;
   }
   strcpy(my_node->OSname, OSname);
   sprintf(my_node->name,"outp_val%04d", ipos+1);
   my_node->access = 0;
   if (htype == NULL) {
       my_node->type = -3;
   } else {
     if (strcmp(htype,"LCM") == 0) {
       my_node->type = -3;
     } else if (strcmp(htype,"XSM") == 0) {
       my_node->type = -4;
     } else if (strcmp(htype,"BINARY") == 0) {
       my_node->type = -5;
     } else if (strcmp(htype,"ASCII") == 0) {
       my_node->type = -6;
     } else if (strcmp(htype,"DA")== 0) {
       my_node->type = -7;
     } else if (strcmp(htype,"I")== 0) {
       my_node->type = -11;
     } else if (strcmp(htype,"S")== 0) {
       my_node->type = -13;
     } else if (strcmp(htype,"D")== 0) {
       my_node->type = -14;
     } else if (strcmp(htype,"B")== 0) {
       my_node->type = -15;
     } else {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,htype);
       PyErr_SetString(PyLifoError, AbortString);
       return NULL;
     }
   }
   clepush(&my_lifo, my_node);
   return Py_None;
}

static PyObject *lifo_pop(lifoobject *self) {
   char *nomsub="lifo_pop";
   PyObject *item=NULL;
   pylcmobject *rv = NULL;
   char lcm_type[8];
   lifo *my_lifo = self->stack;

   PyObject *module = PyImport_ImportModule("lcm");
   if (module == NULL) {
     sprintf(AbortString,"%s: module lcm should be installed first",nomsub);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }
   PyObject *moduleDict = PyModule_GetDict(module);
   PyTypeObject *protocolClass = (PyTypeObject *)PyDict_GetItemString(moduleDict, "new");
   if (protocolClass == NULL) {
     sprintf(AbortString,"%s: unable to find class new in module lcm",nomsub);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }

   lifo_node *my_node = clepop(&my_lifo);
   switch (my_node->type) {
   case 3 :
     strcpy(lcm_type, "LCM");
     goto c37;
   case 4 :
     strcpy(lcm_type, "XSM");
     goto c37;
   case 5 :
     strcpy(lcm_type, "BINARY");
     goto c37;
   case 6 :
     strcpy(lcm_type, "ASCII");
     goto c37;
   case 7 :
     strcpy(lcm_type, "DA");
     c37: 
     rv = (pylcmobject*)PyObject_New(pylcmobject, protocolClass);
     rv->iplist = my_node->value.mylcm;
     rv->impx_lcm = self->impx_lifo;
     strcpy(rv->type_lcm, lcm_type);
     strcpy(rv->name_lcm, my_node->OSname);
     rv->iact_lcm = my_node->access;
     rv->lrda_lcm = my_node->lparam;
     item = (PyObject *) rv;
     Py_INCREF(item);
     break;
   case 11 :
     item = Py_BuildValue("i", my_node->value.ival);
     break;
   case 13 :
     item = Py_BuildValue("s", my_node->value.hval);
     break;
   case 14 :
     item = Py_BuildValue("d", my_node->value.dval);
     break;
   case 15 :
     if (my_node->value.ival == -1) {
       item = Py_False;
     } else {
       item = Py_True;
     }
     break;
   default :
     sprintf(AbortString,"%s: undefined node - use node method",nomsub);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }
   return item;
}

static PyObject *lifo_pos(lifoobject *self, PyObject *args) {
   char *nomsub="lifo_pos";
   PyObject *item=NULL;
   pylcmobject *rv = NULL;
   long ipos = 0;
   char lcm_type[8];
   /* Parse arguments */
   if(!PyArg_ParseTuple(args, "O", &item)) {
       return NULL;
   }
   lifo *my_lifo = self->stack;
   if (PyObject_IsInstance(item, (PyObject *)&PyLong_Type)) {
     ipos = (long)PyLong_AsLong(item);
     if (self->impx_lifo > 0) printf("%s: recover node at position %ld\n", nomsub, ipos);
   } else if (PyObject_IsInstance(item, (PyObject *)&PyUnicode_Type)) {
     int len = PyUnicode_GET_LENGTH(item);
     if (len > 72) {
       sprintf(AbortString,"%s: CHARACTER DATA OVERFLOW.",nomsub);
       PyErr_SetString(PyLifoError, AbortString);
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
     for (n=0; n<my_lifo->nup; n++) {
       lifo_node *my_node = clepos(&my_lifo, n);
       if (strcmp(my_string, my_node->OSname) == 0) {
         ipos = n;
         goto findit;
       }
     }
     sprintf(AbortString,"%s: unable to find a node with OSname=%s",nomsub, my_string);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;

     findit:
     if (self->impx_lifo > 0) printf("%s: recover node %s at position %ld\n", nomsub, my_string, ipos);
     free(my_string);
   }
   if (ipos+1 > my_lifo->nup) {
     sprintf(AbortString,"%s: node index %ld overflow the stack length=%d",nomsub, ipos, my_lifo->nup);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }
       
   PyObject *module = PyImport_ImportModule("lcm");
   if (module == NULL) {
     sprintf(AbortString,"%s: module lcm should be installed first",nomsub);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }
   PyObject *moduleDict = PyModule_GetDict(module);
   PyTypeObject *protocolClass = (PyTypeObject *)PyDict_GetItemString(moduleDict, "new");
   if (protocolClass == NULL) {
     sprintf(AbortString,"%s: unable to find class new in module lcm",nomsub);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }

   lifo_node *my_node = clepos(&my_lifo, ipos);
   switch (my_node->type) {
   case 3 :
     strcpy(lcm_type, "LCM");
     goto c37;
   case 4 :
     strcpy(lcm_type, "XSM");
     goto c37;
   case 5 :
     strcpy(lcm_type, "BINARY");
     goto c37;
   case 6 :
     strcpy(lcm_type, "ASCII");
     goto c37;
   case 7 :
     strcpy(lcm_type, "DA");
     c37: 
     rv = (pylcmobject*)PyObject_New(pylcmobject, protocolClass);
     rv->iplist = my_node->value.mylcm;
     rv->impx_lcm = self->impx_lifo;
     strcpy(rv->type_lcm, lcm_type);
     strcpy(rv->name_lcm, my_node->OSname);
     rv->iact_lcm = my_node->access;
     rv->lrda_lcm = my_node->lparam;
     item = (PyObject *) rv;
     Py_INCREF(item);
     break;
   case 11 :
     item = Py_BuildValue("i", my_node->value.ival);
     break;
   case 13 :
     item = Py_BuildValue("s", my_node->value.hval);
     break;
   case 14 :
     item = Py_BuildValue("d", my_node->value.dval);
     break;
   case 15 :
     if (my_node->value.ival == -1) {
       item = Py_False;
     } else {
       item = Py_True;
     }
     break;
   default :
     sprintf(AbortString,"%s: node %ld is undefined",nomsub, ipos);
     PyErr_SetString(PyLifoError, AbortString);
     return NULL;
   }
   return item;
}

static PyObject *lifo_len(lifoobject *self) {
   lifo *my_lifo = self->stack;
   return Py_BuildValue("i", my_lifo->nup);
}

static PyObject *lifo_osname(lifoobject *self, PyObject *args) {
   char *nomsub="lifo_osname";
   long ipos = 0;
   /* Parse arguments */
   if(!PyArg_ParseTuple(args, "l", &ipos)) {
       return NULL;
   }
   if (self->impx_lifo > 0) printf("%s: recover node OSname at position %ld\n", nomsub, ipos);
       
   lifo *my_lifo = self->stack;
   lifo_node *my_node = clepos(&my_lifo, ipos);
   return Py_BuildValue("s", my_node->OSname);
}

static PyMethodDef lifo_methods[] = {
    {"lib", (PyCFunction)lifo_lib, METH_NOARGS, "Print the table-of-content of the stack"},
    {"push", (PyCFunction)lifo_push, METH_VARARGS, "Push an entry in the stack"},
    {"pushEmpty", (PyCFunction)lifo_pushEmpty, METH_VARARGS, "Push an empty Lcm object in the stack"},
    {"pop", (PyCFunction)lifo_pop, METH_NOARGS, "Pop an entry from the stack"},
    {"node", (PyCFunction)lifo_pos, METH_VARARGS, "Recover an entry from the stack"},
    {"getMax", (PyCFunction)lifo_len, METH_NOARGS, "Return the length of the stack"},
    {"OSName", (PyCFunction)lifo_osname, METH_VARARGS, "Return the OSname of a node"},
    {NULL}  /* Sentinel */
};

static PyMemberDef lifo_members[] = {
    {"_impx", T_INT, offsetof(lifoobject, impx_lifo), 0, "print index"},
    {NULL}  /* Sentinel */
};

static PyTypeObject LifoType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "lifo.new",                /*tp_name*/
    sizeof(lifoobject),        /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)lifo_dealloc,  /*tp_dealloc*/
    (printfunc) lifo_lib,      /*tp_print*/
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
    lifo_methods,              /* tp_methods */
    lifo_members,              /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    lifo_new,                  /* tp_new */
};

static PyModuleDef lifomodule = {
    PyModuleDef_HEAD_INIT,
    "lifo",
    "A Python module for accessing the lifo stack from C code.",
    -1,
    NULL, NULL, NULL, NULL, NULL
};
PyMODINIT_FUNC PyInit_lifo(void) {
    PyObject* m;

    Py_Initialize();
    if (PyType_Ready(&LifoType) < 0) return NULL;

    m = PyModule_Create(&lifomodule);
    if (m == NULL) return NULL;

    PyModule_AddObject(m, "new", (PyObject *)&LifoType);

    /* Initialize new exception object */
    PyLifoError = PyErr_NewException("lifo.PyLifoError", NULL, NULL);

    /* Add exception object to your module */
    PyModule_AddObject(m, "PyLifoError", PyLifoError);
    return m;
}
