
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
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <structmember.h>
#include <setjmp.h>
#include "numpy/arrayobject.h"
#include "pylcm.h"
jmp_buf buf;
enum Tydata {TAB = 0, ENTIER = 1, REEL_SP = 2, STRING_4C = 3, REEL_DP = 4, BOOLEEN = 5,  COMPLEX = 6, LIS = 10, INDEF = 99 };

static char AbortString[132];
static PyObject *PyLcmError = NULL;
static PyTypeObject PyLcmType;     /* shared type-descriptor */

void xabort_c(char *msg){
  fflush(stdout); fflush(stderr);
  PyErr_SetString(PyLcmError, msg);
  longjmp(buf, 1);
}

char *filled_string_with_blank(int nbmots, char *nomchaine) {
  char *chaine_cptee = NULL;
  int i,j;
  chaine_cptee= (char *)malloc(nbmots*4+1);
  strcpy(chaine_cptee, nomchaine);
  for (i = 0; i < (nbmots*4+1); i++) {
    if (chaine_cptee[i] == '\0') {
      for (j = i; j < (nbmots*4+1); j++) {
	chaine_cptee[j] = ' ';
      }
      break;
    }
  }
  return chaine_cptee;
}

static PyObject *lcm_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="lcm_new";
     char *pytype = NULL;
     char *name = NULL;
     char s[73];
     int_32 iact = 0;
     int_32 lrda = 128;
     int_32 impx = 0;
  
     /* Parse arguments */
     static char *kwlist[] = {"pytype", "name", "iact", "lrda", "impx", NULL};
     if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|slll", kwlist,
                                      &pytype, &name, &iact, &lrda, &impx)) {
         return NULL;
     }

     pylcmobject *self = (pylcmobject*)PyObject_New(pylcmobject, type);
     if (name == NULL) {
       long iii=(long)((PyObject_Hash((PyObject *)self)-1)%1000000000000 + 1);
       sprintf(s,"LCM_%12ld",iii);
     } else {
       int len = strlen(name);
       if (len > 72) {
         sprintf(AbortString,"%s: character name overflow",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       strcpy(s,name);
     }	
     if (impx > 0) {
       printf("%s: new lcm type=%s name=%s iact=%d impx=%d\n", nomsub, pytype, s, iact, impx);
       if (strcmp(pytype,"DA") == 0) printf("%s: lrda=%d\n", nomsub, lrda);
     }
     if (iact < 0 || iact > 2) {
       sprintf(AbortString,"%s: invalid iact=%d (0, 1, 2 expected)",nomsub,iact);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     self->impx_lcm = impx;
     self->iact_lcm = iact;
     strcpy(self->type_lcm, pytype);
     if (strcmp(pytype,"LCM") == 0) {
       int_32 medium = 1;
       lcmop_c(&(self->iplist), s, iact, medium, impx);
     } else if (strcmp(pytype,"XSM") == 0) {
       int_32 medium = 2;
       lcmop_c(&(self->iplist), s, iact, medium, impx);
     } else if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0)) {
     } else if (strcmp(pytype,"DA")== 0) {
       self->lrda_lcm = lrda;
     } else if (strcmp(pytype,"LCM_INP") == 0) {
       strcpy(self->type_lcm, "LCM");
       int_32 medium = 1; /* create a lcm file */
       int_32 imode = 2;  /* from ascii file */
       int_32 idir = 2;   /* to import */
       if (iact != 0) {
         sprintf(AbortString,"%s: invalid iact=%d (0 expected)",nomsub,iact);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       char s2[73];
       sprintf(s2,"_%s",s);
       FILE *fp = fopen(s2, "r");
       lcmop_c(&(self->iplist), s, iact, medium, impx);
       lcmexp_c(&(self->iplist), impx, fp, imode, idir);
       fclose(fp);
     } else if (strcmp(pytype,"XSM_INP") == 0) {
       strcpy(self->type_lcm, "XSM");
       int_32 medium = 2;  /* create a xsm file */
       int_32 imode = 2;  /* from ascii file */
       int_32 idir = 2;   /* to import */
       if (iact != 0) {
         sprintf(AbortString,"%s: invalid iact=%d (0 expected)",nomsub,iact);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       char s2[73];
       sprintf(s2,"_%s",s);
       FILE *fp = fopen(s2, "r");
       lcmop_c(&(self->iplist), s, iact, medium, impx);
       lcmexp_c(&(self->iplist), impx, fp, imode, idir);
       fclose(fp);
     } else {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     strcpy(self->name_lcm, s);
     return (PyObject *)self;
   }
}

static void PyLCM_dealloc(pylcmobject *self) {
   if (!setjmp(buf)) {
     char *nomsub="PyLCM_dealloc";
     if ((strncmp(self->type_lcm,"LCM",3)==0) || (strncmp(self->type_lcm,"XSM",3)==0)) {
       char nameobj[73], namedir[13];
       int_32 vide, longueur, nature, access;
       lcminf_c(&(self->iplist), nameobj, namedir, &vide, &longueur, &nature, &access);
       if (strcmp(namedir,"/") == 0) {
         if (self->impx_lcm > 0) printf("%s: call lcmcl_c for object %s of type %s\n",
              nomsub,self->name_lcm,self->type_lcm);
         int_32 iact = 2;
         lcmcl_c(&(self->iplist), iact);
       } else {
         if (self->impx_lcm > 0) printf("%s: desallocate pylcm handle for object %s of type %s\n",
              nomsub,self->name_lcm,self->type_lcm);
       }
     }
   }
}

static PyObject *PyLCM_lib(pylcmobject *self) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="PyLCM_lib";
     char *pytype = self->type_lcm;
     if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
         (strcmp(pytype,"DA") == 0)) {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     fflush(stdout);
     lcmlib_c(&(self->iplist));
     fflush(stdout);
     return Py_None;
   }
}

static PyObject *PyLCM_val(pylcmobject *self) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="PyLCM_val";
     char *pytype = self->type_lcm;
     if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
         (strcmp(pytype,"DA") == 0)) {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     lcmval_c(&(self->iplist)," ");
     return Py_None;
   }
}

static PyObject *PyLCM_keys(pylcmobject *self) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="PyLCM_keys";
     PyObject *ret = NULL;
     PyObject *keys = NULL;
     PyObject *cle = NULL;
     char nom[13], nameobj[73], namedir[13], name[13];
     int_32 vide, longueur, nature, access;
  
     keys=PyList_New(0);
     lcminf_c(&(self->iplist), nameobj, namedir, &vide, &longueur, &nature, &access);
     if (longueur == -1) { /* associative table */
       if(!vide) {
         strcpy(name," ");
         lcmnxt_c(&(self->iplist),name);
         strcpy(nom,name);
         do {
           lcmnxt_c(&(self->iplist),nom);
	   cle = Py_BuildValue("s", nom);
	   PyList_Append(keys,cle);
	   Py_DECREF(cle);
         } while(strcmp(nom,name) != 0);
       }
       ret = keys;
     } else {
       sprintf(AbortString,"%s: associative table expected (nature=%d)",nomsub,nature);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     return ret;
   }
}

static PyObject *PyLCM_len(pylcmobject *self) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="pylcm_len";
     char nameobj[73], namedir[13];
     int_32 vide, longueur, memoire, access;
     char *pytype = self->type_lcm;
     if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
         (strcmp(pytype,"DA") == 0)) {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     lcminf_c(&(self->iplist), nameobj, namedir, &vide, &longueur, &memoire, &access);
     return Py_BuildValue("i", longueur);
   }
}

static PyObject *PyLCM_copy(pylcmobject *self, PyObject *args) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="PyLCM_copy";
     pylcmobject *rv = NULL;
     char *input_name = NULL;
     lcm *lcm_object;
     char s[73];
     int_32 medium = 1;
     char nameobj[73], namedir[13];
     int_32 vide, longueur, memoire, access;
  
     /* Parse arguments */
     if (!PyArg_ParseTuple(args, "|sl", &input_name, &medium)) {
       return NULL;
     }

     if (input_name == NULL) {
       long iii=(long)((PyObject_Hash((PyObject *)self)-1)%1000000000000 + 1);
       sprintf(s,"LCM_%12ld",iii);
     } else {
       int len = strlen(input_name);
       if (len > 72) {
         sprintf(AbortString,"%s: character name overflow",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       strcpy(s,input_name);
     }	
     char *pytype = self->type_lcm;
     if (self->impx_lcm > 0) {
       printf("%s: new lcm type=%s name=%s impx=%d\n", nomsub, pytype, s, self->impx_lcm);
     }
     if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
         (strcmp(pytype,"DA") == 0)) {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     if (medium < 1 || medium > 2) {
       sprintf(AbortString,"%s: invalid medium=%d (1, 2 expected)",nomsub,medium);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     lcminf_c(&(self->iplist), nameobj, namedir, &vide, &longueur, &memoire, &access);
     if (longueur != -1) {
       sprintf(AbortString,"%s: associative table expected for copy",nomsub);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     int_32 imp = 0;
     lcmop_c(&lcm_object, s, imp, medium, self->impx_lcm);
     lcmequ_c(&(self->iplist), &lcm_object);
     rv = (pylcmobject*)PyObject_New(pylcmobject, &PyLcmType);
     rv->iplist = lcm_object;
     rv->impx_lcm = self->impx_lcm;
     rv->iact_lcm = imp;
     if (medium == 1) {
       strcpy(rv->type_lcm, "LCM");
     } else {
       strcpy(rv->type_lcm, "XSM");
     }
     strcpy(rv->name_lcm, s);
     return (PyObject *) rv;
   }
}

static PyObject *PyLCM_rep(pylcmobject *self, PyObject *args) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="pylcm_rep";
     pylcmobject *ret = NULL;
     lcm *lcm_object;
     char nameobj[73], namedir[13];
     int_32 vide, longueur, memoire, access;
     char *pytype = self->type_lcm;
     if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
         (strcmp(pytype,"DA") == 0)) {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     lcminf_c(&(self->iplist), nameobj, namedir, &vide, &longueur, &memoire, &access);
     if (longueur == -1) {
       char *nomcle = NULL;
       if (!PyArg_ParseTuple(args, "s", &nomcle)) {
         return NULL;
       }
       if (strlen(nomcle) > 12) {
         sprintf(AbortString,"%s: character name overflow",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       if (self->impx_lcm > 0) {
         printf("%s: create daughter associative table (key=%s)\n", nomsub, nomcle);
       }
       lcm_object = lcmdid_c(&(self->iplist), nomcle);
     } else {
       int_32 iset;
       if (!PyArg_ParseTuple(args, "l", &iset)) {
         return NULL;
       }
       if (self->impx_lcm > 0) {
         printf("%s: create daughter associative table (key=%d)\n", nomsub, iset);
       }
       lcm_object = lcmdil_c(&(self->iplist), iset);
     }
     ret = (pylcmobject *)PyObject_New(pylcmobject, &PyLcmType);
     ret->iplist = lcm_object;
     ret->impx_lcm = self->impx_lcm;
     strcpy(ret->type_lcm, self->type_lcm);
     strcpy(ret->name_lcm, namedir);
     return (PyObject *)ret;
   }
}

static PyObject *PyLCM_lis(pylcmobject *self, PyObject *args) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="pylcm_lis";
     pylcmobject *ret = NULL;
     lcm *lcm_object;
     char nameobj[73], namedir[13];
     int_32 vide, longueur, memoire, access;
     char *pytype = self->type_lcm;
     long list_length = 0;
     if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
         (strcmp(pytype,"DA") == 0)) {
       sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     lcminf_c(&(self->iplist), nameobj, namedir, &vide, &longueur, &memoire, &access);
     if (longueur == -1) {
       char *nomcle = NULL;
       if (!PyArg_ParseTuple(args, "sl", &nomcle, &list_length)) {
         return NULL;
       }
       if (strlen(nomcle) > 12) {
         sprintf(AbortString,"%s: character name overflow",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       if (self->impx_lcm > 0) {
         printf("%s: create daughter heterogeneous list (key=%s)\n", nomsub, nomcle);
       }
       lcm_object = lcmlid_c(&(self->iplist), nomcle, list_length);
     } else {
       int_32 iset;
       if (!PyArg_ParseTuple(args, "ll", &iset, &list_length)) {
         return NULL;
       }
       if (self->impx_lcm > 0) {
         printf("%s: create daughter heterogeneous list (key=%d)\n", nomsub, iset);
       }
       lcm_object = lcmlil_c(&(self->iplist), iset, list_length);
     }
     ret = (pylcmobject *)PyObject_New(pylcmobject, &PyLcmType);
     ret->iplist = lcm_object;
     ret->impx_lcm = self->impx_lcm;
     strcpy(ret->type_lcm, self->type_lcm);
     strcpy(ret->name_lcm, namedir);
     return (PyObject *)ret;
   }
}

static PyObject *PyLCM_expor(pylcmobject *self, PyObject *args) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="pylcm_expor";
     FILE *file = NULL;
     char * input_name = NULL;
     char namt[74];
  
     if (!PyArg_ParseTuple(args, "|s", &input_name)) {
       return NULL;
     }
     long imode = 2; /* ascii */
     long idir = 1; /* exportation */
     if ((input_name == NULL) || (strcmp(input_name, " ") == 0)) {
       sprintf(namt,"_%s", self->name_lcm);
     } else {
       if (strncmp(input_name,"_",1) != 0) {
         sprintf(AbortString,"%s: leading '_' expected in file name",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       int len = strlen(input_name);
       if (len > 73) {
         sprintf(AbortString,"%s: character name overflow",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return NULL;
       }
       strcpy(namt, input_name);
     }
	
     if (self->impx_lcm > 1) printf("%s: export %s --> %s\n", nomsub, namt, input_name);
     file = fopen(namt, "w");
     lcmexp_c(&(self->iplist), self->impx_lcm, file, imode, idir);
     fclose(file);
     return Py_None;
   }
}

static PyObject *pylcm_dict(pylcmobject *self, PyObject *key) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="pylcm_dict";
     PyObject *ret = NULL;
     char namekey[13], nameobj[73], namedir[13];
     int i;
     lcm *objet;
     int_32 *ndata =NULL;
     int_32 nbmots = 0;
     int_32 tydata = -1;
     char *chaine = NULL;
     int_32 vide, longueur, memoire, access;
     pylcmobject *rv = NULL;
     int D_1 = 1;
     long taille[1];
  
     int len = PyUnicode_GET_LENGTH(key);
     if (len > 12) {
       sprintf(AbortString,"%s: character data overflow",nomsub);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     }
     long kind = PyUnicode_KIND(key);
     void *data = PyUnicode_DATA(key);
     int n;
     for (n=0; n<len; n++) {
       namekey[n] = (Py_UCS4)PyUnicode_READ(kind, data, n);
     }
     namekey[len] = '\0';
     if (self->impx_lcm > 1) printf("%s: dictionary key --> %s\n", nomsub, namekey);

     lcmlen_c(&(self->iplist), namekey, &nbmots, &tydata);
     switch (tydata) {
     case INDEF :
       sprintf(AbortString,"%s: wrong key",nomsub);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     case TAB :
     case LIS : /* creation of a daughter pylcm object */
       objet = lcmgid_c(&(self->iplist), namekey);
       lcminf_c(&objet, nameobj, namedir, &vide, &longueur, &memoire, &access);
       rv = (pylcmobject*)PyObject_New(pylcmobject, &PyLcmType);
       rv->iplist = objet;
       rv->impx_lcm = self->impx_lcm;
       rv->lrda_lcm = 0;
       rv->iact_lcm = self->iact_lcm;
       strcpy(rv->type_lcm,self->type_lcm);
       strcpy(rv->name_lcm,namedir);
       ret = (PyObject *)rv;
       break;
     case STRING_4C : /* string */
       ndata = (int_32 *) malloc(nbmots*(sizeof(*ndata)));
       lcmget_c(&(self->iplist), namekey, ndata);
       chaine = (char *) malloc((int)nbmots*4 + 1); /* +1 pour \0 */
       for (i=0;i<nbmots;i++) strncpy ((chaine+4*i),(char *) (ndata + i), 4);
       chaine[nbmots*4] = '\0';
       ret = Py_BuildValue("s", chaine);
       free(chaine);
       free(ndata);
       ndata = NULL;
       chaine = NULL;
       break;
     default :
       taille[0] = (int)nbmots;
       switch (tydata) {
       case ENTIER :
         ret = PyArray_SimpleNew(D_1, taille, NPY_INT32);
         break;
       case REEL_SP :
         ret = PyArray_SimpleNew(D_1, taille, NPY_FLOAT);
         break;
       case REEL_DP :
         ret = PyArray_SimpleNew(D_1, taille, NPY_DOUBLE);
         break;
       case BOOLEEN :
         ret = PyArray_SimpleNew(D_1, taille, NPY_BOOL);
         break;
       case COMPLEX :
         ret = PyArray_SimpleNew(D_1, taille, NPY_CFLOAT);
         break;
       }
       int_32 *tabdata = (int_32 *)PyArray_DATA((PyArrayObject *)ret);
       lcmget_c(&(self->iplist), namekey, tabdata);
     }
     return ret;
  }
}

static long pylcm_assign_dict(pylcmobject *self, PyObject *key, PyObject *v) {
   if (setjmp(buf)) {
     return -1;
   } else {
     char *nomsub="pylcm_assign_dict";
     long ret = 0;
     int_32 nbmots = 0;
     int_32 tydata = 99; /* undefined */
     int_32 *tabdata = NULL;
     char namekey[13];

     int len = PyUnicode_GET_LENGTH(key);
     if (len > 12) {
       sprintf(AbortString,"%s: character data overflow",nomsub);
       PyErr_SetString(PyLcmError, AbortString);
       return -1;
     }
     long kind = PyUnicode_KIND(key);
     void *data = PyUnicode_DATA(key);
     int n;
     for (n=0; n<len; n++) {
       namekey[n] = (Py_UCS4)PyUnicode_READ(kind, data, n);
     }
     namekey[len] = '\0';
     if (self->impx_lcm > 1) printf("%s: dictionary key --> %s\n", nomsub, namekey);

     if(PyUnicode_Check(v)){
       char *filled_string = NULL;
       tydata = STRING_4C;

       int len = PyUnicode_GET_LENGTH(v);
       int kind = PyUnicode_KIND(v);
       void *data = PyUnicode_DATA(v);
       char* my_string;
       my_string = (char *)malloc(len+1);
       int n;
       for (n=0; n<len; n++) {
         my_string[n] = (Py_UCS4)PyUnicode_READ(kind, data, n);
       }
       my_string[len] = '\0';
       nbmots = len/4 + (len%4+3)/4;

       filled_string = filled_string_with_blank(nbmots,my_string);
       free(my_string);
       tabdata = (int_32 *) filled_string;
       lcmput_c(&(self->iplist), namekey, nbmots, tydata, tabdata);
       free(filled_string);
       filled_string = NULL;
     } else if(PyArray_Check(v)){
       long type_num = PyArray_TYPE((PyArrayObject *)v);
       long nb_elts = PyArray_SIZE((PyArrayObject *)v);
       nbmots = nb_elts;
       tabdata = (int_32 *)PyArray_DATA((PyArrayObject *)v);
       switch (type_num) {
       case NPY_INT32 :
         tydata = ENTIER;
         break;
       case NPY_FLOAT :
         tydata = REEL_SP;
         break;
       case NPY_DOUBLE :
         tydata = REEL_DP;
         break;
       case NPY_BOOL :
         tydata = BOOLEEN;
         break;
       case NPY_CFLOAT :
         tydata = COMPLEX;
         break;
       default :
         sprintf(AbortString,"%s: numerical array type %ld is not implemented",nomsub, type_num);
         PyErr_SetString(PyLcmError, AbortString);
         return -1;
       }
       lcmput_c(&(self->iplist), namekey, nbmots, tydata, tabdata); 	
     } else {
       if (((PyObject *)v)->ob_type == &PyLcmType) {
         lcm *rep_vide;
         int_32 vide, longueur, memoire, access;
         char nameobj[73], namedir[13];
         pylcmobject *new_pylcm = (pylcmobject *)v;
         char *pytype = new_pylcm->type_lcm;
         if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
             (strcmp(pytype,"DA") == 0)) {
           sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
           PyErr_SetString(PyLcmError, AbortString);
           return -1;
         }
         rep_vide = lcmgid_c(&(self->iplist), namekey);
         lcminf_c(&rep_vide,nameobj,namedir,&vide, &longueur, &memoire, &access);
         if (longueur == -1 && memoire) {
           lcmequ_c(&(new_pylcm->iplist),&rep_vide);
         } else {
           sprintf(AbortString,"%s: use rep method before making a copy",nomsub);
           PyErr_SetString(PyLcmError, AbortString);
           return -1;
         }
       } else {
         sprintf(AbortString,"%s: cannot assign this object type",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return -1;
       }
     }
     return ret;
   }
}

static PyObject *pylcm_list(pylcmobject *self, PyObject *indice) {
   if (setjmp(buf)) {
     return NULL;
   } else {
     char *nomsub="pylcm__list";
     char nameobj[73], namedir[13];
     int i;
     int_32 *ndata =NULL;
     PyObject *ret = NULL;
     int_32 ind = -1;
     int_32 nbmots = 0;
     int_32 tydata = -1;
     lcm *objet;
     char *chaine = NULL;
     int_32 vide, longueur, memoire, access;
     pylcmobject *rv = NULL;
     int D_1 = 1;
     long taille[1];
  
     ind = PyLong_AsLong(indice);
     lcmlel_c(&(self->iplist), ind, &nbmots, &tydata);
     switch (tydata) {
     case INDEF :
       sprintf(AbortString,"%s: wrong key",nomsub);
       PyErr_SetString(PyLcmError, AbortString);
       return NULL;
     case TAB :
     case LIS :
       objet = lcmgil_c(&(self->iplist), ind);
       lcminf_c(&objet, nameobj, namedir, &vide, &longueur, &memoire, &access);
       rv = (pylcmobject*)PyObject_New(pylcmobject, &PyLcmType);
       rv->iplist = objet;
       rv->impx_lcm = self->impx_lcm;
       rv->lrda_lcm = 0;
       rv->iact_lcm = self->iact_lcm;
       strcpy(rv->type_lcm,self->type_lcm);
       strcpy(rv->name_lcm,namedir);
       ret = (PyObject *)rv;
       break;
     case STRING_4C : /* type chaine */
       ndata = (int_32 *) malloc(nbmots*(sizeof(*ndata)));
       lcmgdl_c(&(self->iplist), ind, ndata);
       chaine = (char *) malloc((int)nbmots*4 + 1); /* +1 for \0 */
       for (i=0;i<nbmots;i++) strncpy ((chaine+4*i),(char *) (ndata + i), 4);
       chaine[nbmots*4] = '\0';
       ret = Py_BuildValue("s", chaine);
       free(chaine);
       free(ndata);
       ndata = NULL;
       chaine = NULL;
       break;
     default :
       taille[0] = (int)nbmots;
       switch (tydata) {
       case ENTIER :
         ret = PyArray_SimpleNew(D_1, taille, NPY_INT32);
         break;
       case REEL_SP :
         ret = PyArray_SimpleNew(D_1, taille, NPY_FLOAT);
         break;
       case REEL_DP :
         ret = PyArray_SimpleNew(D_1, taille, NPY_DOUBLE);
         break;
       case BOOLEEN :
         ret = PyArray_SimpleNew(D_1, taille, NPY_BOOL);
         break;
       case COMPLEX :
         ret = PyArray_SimpleNew(D_1, taille, NPY_CFLOAT);
         break;
       }
       int_32 *tabdata = (int_32 *)PyArray_DATA((PyArrayObject *)ret);
       lcmgdl_c(&(self->iplist), ind, tabdata);
     }
     return ret;
   }
}

static long pylcm_assign_list(pylcmobject *self, PyObject* i, PyObject *v) {
   if (setjmp(buf)) {
     return -1;
   } else {
     char *nomsub="pylcm_assign_list";
     long ret = 0;
     int_32 nbmots = 0;
     int_32 tydata = 99; /* undefined */
     int_32 *tabdata = NULL;
     int_32 pos = PyLong_AsLong(i);
  
     if(PyUnicode_Check(v)){
       char *filled_string = NULL;
       tydata = STRING_4C;

       int len = PyUnicode_GET_LENGTH(v);
       int kind = PyUnicode_KIND(v);
       void *data = PyUnicode_DATA(v);
       char* my_string;
       my_string = (char *)malloc(len+1);
       int n;
       for (n=0; n<len; n++) {
         my_string[n] = (Py_UCS4)PyUnicode_READ(kind, data, n);
       }
       my_string[len] = '\0';
       nbmots = len/4 + (len%4+3)/4;

       filled_string = filled_string_with_blank(nbmots,my_string);
       free(my_string);
       tabdata = (int_32 *) filled_string;
       lcmpdl_c(&(self->iplist), pos, nbmots, tydata, tabdata);
       free(filled_string);
       filled_string=NULL;
     } else if(PyArray_Check(v)) {
       long type_num = PyArray_TYPE((PyArrayObject *)v);
       long nb_elts = PyArray_SIZE((PyArrayObject *)v);
       nbmots = nb_elts;
       tabdata = (int_32 *)PyArray_DATA((PyArrayObject *)v);
       switch (type_num) {
       case NPY_INT32 :
         tydata = ENTIER;
         break;
       case NPY_FLOAT :
         tydata = REEL_SP;
         break;
       case NPY_DOUBLE :
         tydata = REEL_DP;
         break;
       case NPY_BOOL :
         tydata = BOOLEEN;
         break;
       case NPY_CFLOAT :
         tydata = COMPLEX;
         break;
       default :
         sprintf(AbortString,"%s: numerical array type %ld is not implemented",nomsub, type_num);
         PyErr_SetString(PyLcmError, AbortString);
         return -1;
       }
       lcmpdl_c(&(self->iplist), pos, nbmots, tydata, tabdata); 	
     } else {
       if (((PyObject *)v)->ob_type == &PyLcmType) {
         lcm *rep_vide;
         int_32 vide, longueur, memoire, access;
         char nameobj[73], namedir[13];
         pylcmobject *new_pylcm = (pylcmobject *)v;
         char *pytype = new_pylcm->type_lcm;
         if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
             (strcmp(pytype,"DA") == 0)) {
           sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
           PyErr_SetString(PyLcmError, AbortString);
           return -1;
         }
         rep_vide = lcmgil_c(&(self->iplist), pos);
         lcminf_c(&rep_vide, nameobj, namedir, &vide, &longueur, &memoire, &access);
         if (longueur == -1 && memoire){
           lcmequ_c(&(new_pylcm->iplist),&rep_vide);
         } else {
           sprintf(AbortString,"%s: use rep method before making a copy",nomsub);
           PyErr_SetString(PyLcmError, AbortString);
           return -1;
         }
       } else {
         sprintf(AbortString,"%s: cannot assign this object type",nomsub);
         PyErr_SetString(PyLcmError, AbortString);
         return -1;
       }
     }
     return ret;
   }
}

static PyObject *PyLCM_subscript(pylcmobject *self, PyObject *index) {
   char *nomsub="PyLCM_subscript";
   PyObject *ret = NULL;
   char *pytype = self->type_lcm;
   if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
       (strcmp(pytype,"DA") == 0)) {
     sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
     PyErr_SetString(PyLcmError, AbortString);
     return NULL;
   }
   if(PyUnicode_Check(index)) { /* dictionary */
     ret = pylcm_dict(self, index);
   } else if (PyLong_Check(index)) { /* list */
     ret = pylcm_list(self, index);
   } else {
     sprintf(AbortString,"%s: invalid index type",nomsub);
     PyErr_SetString(PyLcmError, AbortString);
     return NULL;
   }
   return ret;
}

static long PyLCM_assign_subscript(pylcmobject *self, PyObject *index, PyObject *v) {
   char *nomsub="PyLCM_assign_subscript";
   long ret = 0;
   char *pytype = self->type_lcm;
   if ((strcmp(pytype,"BINARY") == 0) || (strcmp(pytype,"ASCII") == 0) ||
       (strcmp(pytype,"DA") == 0)) {
     sprintf(AbortString,"%s: invalid type=%s",nomsub,pytype);
     PyErr_SetString(PyLcmError, AbortString);
     return -1;
   }
   if(PyUnicode_Check(index)) {
     ret = pylcm_assign_dict(self, index, v);
   } else if (PyLong_Check(index)) {
     ret = pylcm_assign_list(self, index, v);
   } else {
     sprintf(AbortString,"%s: invalid index type",nomsub);
     PyErr_SetString(PyLcmError, AbortString);
     return -1;
   }
   return ret;
}

int PyLCM_length(PyObject *self) {
   return 0;
}

static PyMethodDef lcm_methods[] = {
    {"lib", (PyCFunction)PyLCM_lib, METH_NOARGS, "Print the table-of-content of the lcm object"},
    {"val", (PyCFunction)PyLCM_val, METH_NOARGS, "Validate a lcm object"},
    {"keys", (PyCFunction)PyLCM_keys, METH_NOARGS, "Set the keys content of an associative table"},
    {"len", (PyCFunction)PyLCM_len, METH_NOARGS, "Return the length of a lcm object"},
    {"copy", (PyCFunction)PyLCM_copy, METH_VARARGS, "Deep copy of a pylcm object"},
    {"rep", (PyCFunction)PyLCM_rep, METH_VARARGS, "Create a daughter associative table"},
    {"lis", (PyCFunction)PyLCM_lis, METH_VARARGS, "Create a daughter heterogeneous list"},
    {"expor", (PyCFunction)PyLCM_expor, METH_VARARGS, "Export a lcm objects towards an ascii file"},
    {NULL}  /* Sentinel */
};

static PyMemberDef lcm_members[] = {
    {"_impx", T_INT, offsetof(pylcmobject, impx_lcm), 0, "print index"},
    {"_access", T_INT, offsetof(pylcmobject, iact_lcm), READONLY, "access index"},
    {"_type", T_STRING_INPLACE, offsetof(pylcmobject, type_lcm), READONLY, "object type"},
    {"_name", T_STRING_INPLACE, offsetof(pylcmobject, name_lcm), READONLY, "object name"},
    {NULL}  /* Sentinel */
};

static PyMappingMethods PyLCM_as_mapping = {
    (lenfunc)PyLCM_length,                  /* mp_length            */
    (binaryfunc)PyLCM_subscript,            /* mp_subscript         */
    (objobjargproc)PyLCM_assign_subscript,  /* mp_assign_subscript  */
};

static PyTypeObject PyLcmType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "lcm.new",                 /*tp_name*/
    sizeof(pylcmobject),       /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PyLCM_dealloc, /*tp_dealloc*/
    (printfunc) PyLCM_lib,     /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    &PyLCM_as_mapping,         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,  /*tp_flags*/
    "Custom objects",          /*tp_doc*/
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    lcm_methods,               /* tp_methods */
    lcm_members,               /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    lcm_new,                   /* tp_new */
};

static PyModuleDef lcmmodule = {
    PyModuleDef_HEAD_INIT,
    "lcm",
    "A Python module for accessing LCM objects from C code.",
    -1,
    NULL, NULL, NULL, NULL, NULL
};
PyMODINIT_FUNC PyInit_lcm(void) {
    PyObject* m;

    Py_Initialize();
    import_array();
    if (PyType_Ready(&PyLcmType) < 0) return NULL;

    m = PyModule_Create(&lcmmodule);
    if (m == NULL) return NULL;

    PyModule_AddObject(m, "new", (PyObject *)&PyLcmType);

    /* Initialize new exception object */
    PyLcmError = PyErr_NewException("lcm.PyLcmError", NULL, NULL);

    /* Add exception object to your module */
    PyModule_AddObject(m, "PyLcmError", PyLcmError);
    return m;
}
