/*
freesteam - IAPWS-IF97 steam tables library
Copyright (C) 2004-2009  John Pye

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*//** @file
	Functions to return SteamState objects for points on the
	boundary of the ...
*/

#define FREESTEAM_BUILDING_LIB
#include "bounds.h"

#include "region.h"
#include "b23.h"
#include "zeroin.h"

#include <stdlib.h>

typedef struct{
	double p, T;
} SteamPTData;

static double pT3_fn(double rho, void *user_data){
#define D ((SteamPTData *)user_data)
	return D->p - freesteam_region3_p_rhoT(rho, D->T);
#undef D
}

SteamState freesteam_bound_pmax_T(double T){
	SteamState S;
	if(T <= REGION1_TMAX){
		S.region = 1;
		S.valueR.R1.p = IAPWS97_PMAX;
		S.valueR.R1.T = T;
	}else if(T > freesteam_region2_s_pT(IAPWS97_PMAX,freesteam_b23_T_p(IAPWS97_PMAX))){
		S.region = 2;
		S.valueR.R2.p = IAPWS97_PMAX;
		S.valueR.R2.T = T;
	}else{
		double tol = 1e-7;
		double sol, err = 0;
                double lb, ub;
		SteamPTData D = {IAPWS97_PMAX, T};
		S.region = 3;
		S.valueR.R3.T = T;
		lb = 1./freesteam_region2_v_pT(freesteam_b23_p_T(T),T);
		ub = 1./freesteam_region1_v_pT(IAPWS97_PMAX,REGION1_TMAX);

		if(zeroin_solve(&pT3_fn, &D, lb, ub, tol, &sol, &err)){
			fprintf(stderr,"%s (%s:%d): failed to solve for rho\n",__func__,__FILE__,__LINE__);
			exit(1);
		}
		S.valueR.R3.rho = sol;
	}
	return S;
}
