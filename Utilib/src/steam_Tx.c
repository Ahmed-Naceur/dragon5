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
#include "steam_Tx.h"
#include "region.h"
#include "b23.h"
#include "zeroin.h"

/* .... */

#include <stdlib.h>
#include <assert.h>
#include <math.h>


int freesteam_bounds_Tx(double T, double x, int verbose){

#define BOUND_WARN(MSG) \
	if(verbose){\
		fprintf(stderr,"%s (%s:%d): WARNING " MSG " (T = %g K, x = %g)\n"\
		,__func__,__FILE__,__LINE__,T,x);\
	}

	if(T <= IAPWS97_TMIN){
		BOUND_WARN("T <= TMIN");
		return 1;
	}
	if(T > IAPWS97_TCRIT){
		BOUND_WARN("T > TCRIT");
		return 2;
	}

	if(x < 0){
		BOUND_WARN("x < 0");
		return 3;
	}else if(x > 1){
		BOUND_WARN("x > 1");
		return 4;
	}
	return 0;
#undef BOUND_WARN
}

int freesteam_region_Tx(double T, double x){
	if(T >= IAPWS97_TCRIT)return 3;

	if(x <= 0){
		if(T > REGION1_TMAX)return 3;
		return 1;
	}

	if(x >= 1){
		if(T > REGION1_TMAX)return 3;
		return 2;
	}

	return 4;
}

typedef struct{
	double T, s;
} SolveTSData;

static double Ts_region3_fn(double rho, void *user_data){
#define D ((SolveTSData *)user_data)
	return D->s - freesteam_region3_s_rhoT(rho, D->T);
#undef D
}

/**
	This function will always return saturated mixtures; no negative or >1 
	values of x are being 'understood' here (although one can give them meaning
	based on extrapolated values of u or h of v, for example...)
*/
SteamState freesteam_set_Tx(double T, double x){
	SteamState S;

	if(T >= IAPWS97_TCRIT){
		/* region 3 supercritical. just return a state with the specified
		temperature and the critical point entropy. arbitrary. */
                double ub, lb, tol, sol, err;
		SolveTSData D = {T, freesteam_region3_s_rhoT(IAPWS97_RHOCRIT, IAPWS97_TCRIT)};
		ub = 1./freesteam_region1_v_pT(IAPWS97_PMAX,REGION1_TMAX);
		lb = 1./freesteam_region2_v_pT(freesteam_b23_p_T(T),T);
		tol = 1e-7;
		err = 0;
		if(zeroin_solve(&Ts_region3_fn, &D, lb, ub, tol, &sol, &err)){
			fprintf(stderr,"%s (%s:%d): failed to solve for rho\n",__func__,__FILE__,__LINE__);
			exit(1);
		}
		S.region = 3;
		S.valueR.R3.T = T;
		S.valueR.R3.rho = sol;
	}else if(x <= 0){
		if(T > REGION1_TMAX){
			S.region = 3;
			S.valueR.R3.T = T;
			S.valueR.R3.rho = freesteam_region4_rhof_T(T);	
			/* FIXME iteratively refine the value */
		}else{
			S.region = 1;
			S.valueR.R1.p = freesteam_region4_psat_T(T);
			S.valueR.R1.T = T;
		}
	}else if(x >= 1){
		if(T > REGION1_TMAX){
			S.region = 3;
			S.valueR.R3.T = T;
			S.valueR.R3.rho = freesteam_region4_rhog_T(T);	
			/* FIXME iteratively refine the value */
		}else{
			S.region = 2;
			S.valueR.R1.p = freesteam_region4_psat_T(T);
			S.valueR.R1.T = T;
		}
	}else{
		/* finally! */
		S.region = 4;
		S.valueR.R4.T = T;
		S.valueR.R4.x = x;
	}

	return S;
}
