
/**********************************/
/* C API for kdi file support     */
/* author: A. Hebert (30/04/2002) */
/**********************************/

/*
Copyright (C) 2002 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#if defined(CRAY)
#define lnword 8
#else
#define lnword 4
#endif

#if !defined(MSDOS)
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include "kdi.h"

int_32 offset;

kdi_file * kdiop_c(char *nomC,int_32 iactio)
{
  kdi_file *my_file;
  my_file = (kdi_file *) malloc(sizeof(*my_file));
  strcpy(my_file->nom,nomC);
  my_file->fd = NULL;
  if (iactio == 0) {
     FILE *file ;
     long fd;
     if ( ( file = fopen(nomC,"rb") ) != NULL ) {
        fclose(file) ;
        perror ("open error 0 in kdiop_c ");
        return NULL;
     }
     fd = creat(nomC,0600);
     close(fd);
     my_file->fd = fopen(nomC,"r+b");
  } else if (iactio == 1) {
     FILE *file ;
     if ( ( file = fopen(nomC,"rb") ) == NULL ) {
        perror ("open error 1 in kdiop_c ");
        return NULL;
     }
     fclose(file) ;
     my_file->fd = fopen(nomC,"r+b");
  } else if (iactio == 2) {
     FILE *file ;
     if ( ( file = fopen(nomC,"rb") ) == NULL ) {
        perror ("open error 2 in kdiop_c ");
        return NULL;
     }
     my_file->fd = file;
  } else {
     return NULL;
  }
  if ( my_file->fd == NULL ) {
     perror ("open error 3 in kdiop_c ");
     return NULL;
  }
  return my_file;
}

int_32 kdiput_c(kdi_file *my_file,int_32 *data,int_32 iofset,int_32 length)
{
  int_32 irc=0;
  offset=iofset*lnword;
  if (my_file == NULL) {
     irc = -1;
  } else if (fseek(my_file->fd,offset,0) >= 0) {
     int_32 n, iof=0;
     while ((n = fwrite(&data[iof],lnword,length,my_file->fd)) < length-iof) {
        if (n < 0) return n-1;
        iof+=n;
     }
  } else {
     irc = -3;
  }
  return irc;
}

int_32 kdiget_c(kdi_file *my_file,int_32 *data,int_32 iofset,int_32 length)
{
  int_32 irc=0;
  offset=iofset*lnword;
  if (my_file == NULL) {
     irc = -1;
  } else if (fseek(my_file->fd,offset,0) >= 0) {
     int_32 n, iof=0;
     while ((n = fread(&data[iof],lnword,length,my_file->fd)) < length-iof) {
        if (n < 0) return n-1;
        iof+=n;
     }
  } else {
     irc = -3;
  }
  return irc;
}

int_32 kdicl_c(kdi_file *my_file,int_32 istatu)
{
  long irc;
  if (my_file == NULL) {
     irc = -1;
  } else if (istatu == 1) {
     irc = fclose(my_file->fd);
     free(my_file);
  } else if (istatu == 2) {
     irc = fclose(my_file->fd);
     if (irc != 0 ) {
        perror ("close error 1 in kdicl_c ");
        return irc;
     }
     irc = remove(my_file->nom);
     free(my_file);
  } else {
     irc = -999;
  }
  my_file = NULL;
  if (irc != 0) perror ("close error 2 in kdicl_c ");
  return irc;
}
