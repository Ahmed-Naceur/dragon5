
/**********************************/
/* C API for freesteam support    */
/* author: A. Hebert (27/05/2012) */
/**********************************/

/*
   Copyright (C) 2012 Ecole Polytechnique de Montreal

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
 */

#include "steam.h"
#include "steam_pT.h"
#include "steam_Tx.h"
#include "region.h"

void free_pT(double *P, double *T, double *RHO, double *H, double *ZK, double *ZMU, double *CP){
  SteamState State = freesteam_set_pT((double)*P, (double)*T);
  *RHO = freesteam_rho(State);  /*density*/
  *H = freesteam_h(State);      /*enthalpy*/
  *ZK = freesteam_k(State);     /*thermal conductivity*/
  *ZMU = freesteam_mu(State);   /*dynamic viscosity*/
  *CP = freesteam_cp(State);    /*isobaric heat capacity*/
}

void free_Tx(double *T, double *X, double *RHO, double *H, double *ZK, double *ZMU, double *CP){
  SteamState State = freesteam_set_Tx((double)*T, (double)*X);
  *RHO = freesteam_rho(State);  /*density*/
  *H = freesteam_h(State);      /*enthalpy*/
  *ZK = freesteam_k(State);     /*thermal conductivity*/
  *ZMU = freesteam_mu(State);   /*dynamic viscosity*/
  *CP = freesteam_cp(State);    /*isobaric heat capacity*/
}
