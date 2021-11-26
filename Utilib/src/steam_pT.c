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
*/
#define FREESTEAM_BUILDING_LIB
#include "steam_pT.h"
#include "region.h"
#include "zeroin.h"
#include "b23.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef struct{
	double p, T;
} SteamPTData;

static double pT_region3_fn(double rho, void *user_data){
#define D ((SteamPTData *)user_data)
	return D->p - freesteam_region3_p_rhoT(rho, D->T);
#undef D
}

/**
	This function will never return region 4, because it's not possible
	to 'sit on the knife' of saturation. If you need to set saturated states,
	you should use another function such as freesteam_region1_set_Tx.
*/
SteamState freesteam_set_pT(double p, double T){
        double tol, sol, err;
	SteamState S;
	if(T < REGION1_TMAX){
		if(p > freesteam_region4_psat_T(T)){
			S.region = 1;
			S.valueR.R1.T = T;
			S.valueR.R1.p = p;
		}else{
			S.region = 2;
			S.valueR.R2.T = T;
			S.valueR.R2.p = p;
		}
	}else{
		/* FIXME some optimisation possiblxe here with test for lower pressures */
		double T23 = freesteam_b23_T_p(p);
		double p23min = freesteam_b23_p_T(REGION1_TMAX);
		if(p < p23min || T > T23){
			S.region = 2;
			S.valueR.R2.T = T;
			S.valueR.R2.p = p;
		}else{
			/* FIXME the limit values are all wrong here! */
                        double ub, lb;
			SteamPTData D = {p,T};
			ub = 1./freesteam_region1_v_pT(IAPWS97_PMAX,REGION1_TMAX);
			lb = 1./freesteam_region2_v_pT(freesteam_b23_p_T(T),T);
			/* if we're in the little wee area around the critical pt... */
			if(T < IAPWS97_TCRIT){
				double psat = freesteam_region4_psat_T(T);
				if(p < psat){
					ub = freesteam_region4_rhog_T(T);
					assert(lb<ub);
				}else{
					lb = freesteam_region4_rhof_T(T);
					assert(lb<ub);
				}
			}
			tol = 1e-7;
			err = 0;
			if(zeroin_solve(&pT_region3_fn, &D, lb, ub, tol, &sol, &err)){
				fprintf(stderr,"%s (%s:%d): failed to solve for rho\n",__func__,__FILE__,__LINE__);
				exit(1);
			}
			S.region = 3;
			S.valueR.R3.T = T;
			S.valueR.R3.rho = sol;
		}
	}
	return S;
}		
