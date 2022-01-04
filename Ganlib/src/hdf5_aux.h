
/**********************************/
/* C API for hdf5 file support    */
/* (auxiliary functions)          */
/* author: A. Hebert (30/11/2021) */
/**********************************/

/*
   Copyright (C) 2021 Ecole Polytechnique de Montreal

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
 */

#if defined(HDF5_LIB)
#include "ganlib.h"
#include "hdf5.h"
void hdf5_open_file_c(const char *, hid_t *, int_32);
void hdf5_close_file_c(hid_t *);
void hdf5_list_c(hid_t *, const char *);
void hdf5_get_dimensions_c(hid_t *, const char *, int_32 *);
void hdf5_get_num_group_c(hid_t *, const char *, int_32 *);
void hdf5_list_datasets_c(hid_t *, const char *, int_32 *, char *idata);
void hdf5_list_groups_c(hid_t *, const char *, int_32 *, char *idata);
void hdf5_info_c(hid_t *, const char *, int_32 *, int_32 *, int_32 *, int_32 *);
void hdf5_read_data_int_c(hid_t *, const char *, int_32 *);
void hdf5_read_data_real4_c(hid_t *, const char *, float *);
void hdf5_read_data_real8_c(hid_t *, const char *, double *);
void hdf5_read_data_string_c(hid_t *, const char *, char *);
void hdf5_write_data_int_c(hid_t *, const char *, int_32, int_32 *, int_32 *);
void hdf5_write_data_real4_c(hid_t *, const char *, int_32, int_32 *, float *);
void hdf5_write_data_real8_c(hid_t *, const char *, int_32, int_32 *, double *);
void hdf5_write_data_string_c(hid_t *, const char *, int_32, int_32, int_32 *, char *);
#endif
