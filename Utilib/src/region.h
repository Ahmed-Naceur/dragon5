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

#ifndef FREESTEAM_REGION_H
#define FREESTEAM_REGION_H

#include "common.h"

FREESTEAM_DLL double freesteam_region1_u_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_v_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_s_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_h_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_cp_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_cv_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_w_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_a_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_g_pT(double p, double T);

FREESTEAM_DLL double freesteam_region2_v_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_u_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_s_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_h_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_cp_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_cv_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_w_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_a_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_g_pT(double p, double T);

FREESTEAM_DLL double freesteam_region3_p_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_u_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_s_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_h_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_cp_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_cv_rhoT(double rho, double T);

FREESTEAM_DLL double freesteam_region4_psat_T(double T);
FREESTEAM_DLL double freesteam_region4_Tsat_p(double p);

FREESTEAM_DLL double freesteam_region4_rhof_T(double T);
FREESTEAM_DLL double freesteam_region4_rhog_T(double T);

FREESTEAM_DLL double freesteam_region4_v_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_u_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_h_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_s_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_cp_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_cv_Tx(double T, double x);

FREESTEAM_DLL double freesteam_region4_dpsatdT_T(double T);

/* used in calculations of derivatives, see derivs.c */
double freesteam_region1_alphav_pT(double p, double T);
double freesteam_region1_kappaT_pT(double p, double T);

double freesteam_region2_alphav_pT(double p, double T);
double freesteam_region2_kappaT_pT(double p, double T);

double freesteam_region3_alphap_rhoT(double rho, double T);
double freesteam_region3_betap_rhoT(double rho, double T);

FREESTEAM_DLL double freesteam_region1_u_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_v_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_s_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_h_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_cp_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_cv_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_w_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_a_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_g_pT(double p, double T);

/* used in calculations of derivatives, see derivs.c */
double freesteam_region1_alphav_pT(double p, double T);
double freesteam_region1_kappaT_pT(double p, double T);

FREESTEAM_DLL double freesteam_region2_v_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_u_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_s_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_h_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_cp_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_cv_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_w_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_a_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_g_pT(double p, double T);

/* used in calculations of derivatives, see derivs.c */
double freesteam_region2_alphav_pT(double p, double T);
double freesteam_region2_kappaT_pT(double p, double T);


FREESTEAM_DLL double freesteam_region3_p_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_u_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_s_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_h_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_cp_rhoT(double rho, double T);
FREESTEAM_DLL double freesteam_region3_cv_rhoT(double rho, double T);

/* FIXME implement freesteam_region3_w_rhoT */

/* used in calculations of derivatives, see derivs.c */
double freesteam_region3_alphap_rhoT(double rho, double T);
double freesteam_region3_betap_rhoT(double rho, double T);

FREESTEAM_DLL double freesteam_region4_psat_T(double T);
FREESTEAM_DLL double freesteam_region4_Tsat_p(double p);

FREESTEAM_DLL double freesteam_region4_rhof_T(double T);
FREESTEAM_DLL double freesteam_region4_rhog_T(double T);

FREESTEAM_DLL double freesteam_region4_v_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_u_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_h_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_s_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_cp_Tx(double T, double x);
FREESTEAM_DLL double freesteam_region4_cv_Tx(double T, double x);

FREESTEAM_DLL double freesteam_region4_dpsatdT_T(double T);

#define REGION1_TMAX 623.15 /* K */
#define REGION2_TMAX 1073.15

#endif
