#ifndef FREESTEAM_COMMON_H
#define FREESTEAM_COMMON_H

#include "config.h"

#define FREESTEAM_CHAR int

/*
	ASCEND code in base/generic only EXPORTS symbols, no imports.
	The FREESTEAM_DLLSPEC macro will, depending on whether we are
	FREESTEAM_BUILDING_LIBASCEND (building libascend.so aka ascend.dll)
	or FREESTEAM_BUILDING_INTERFACE (building for example _ascpy.dll or
	ascendtcl.dll), act respectively to declare symbols as being
	*exported* or *imported*.

	New versions of GCC are able to make use of these declarations
	as well.
*/
#ifdef __WIN32__
# define FREESTEAM_EXPORT __declspec(dllexport)
# define FREESTEAM_IMPORT __declspec(dllimport)
#else
# ifdef HAVE_GCCVISIBILITY
#  define FREESTEAM_EXPORT __attribute__ ((visibility("default")))
#  define FREESTEAM_IMPORT
# else
#  define FREESTEAM_EXPORT
#  define FREESTEAM_IMPORT
# endif
#endif

#ifdef FREESTEAM_BUILDING_LIB
# define FREESTEAM_DLL extern FREESTEAM_EXPORT
#else
# define FREESTEAM_DLL extern FREESTEAM_IMPORT
#endif

#if !defined(FREESTEAM_DLL) || !defined(FREESTEAM_EXPORT) || !defined(FREESTEAM_IMPORT)
# error "NO FREESTEAM_DLL DEFINED"
#endif

/* Constants used throughout IAPWS-IF97 */

#define IAPWS97_PMAX 100e6 /* Pa */
#define IAPWS97_TMIN 273.15 /* K */
#define IAPWS97_TMAX 1073.15 /* K */

#define IAPWS97_TCRIT 647.096 /* K */
#define IAPWS97_PCRIT 22.064e6 /* Pa */
#define IAPWS97_RHOCRIT 322. /* kg/m3 */

#define IAPWS97_PTRIPLE 611.657 /* Pa */

#define IAPWS97_R 461.526 /* J/kgK */

#ifndef __GNUC__
# define __func__ "<function>"
#endif

#ifdef IAPWS97_WARN_APPROX
# define IAPWS97_APPROXIMATE \
	static char _warn_approx=0; \
	if(!_warn_approx){ \
		_warn_approx = 1; \
		fprintf(stderr \
			,"WARNING: %s (%s:%d): backwards or approximation function used!\n" \
			,__func__,__FILE__,__LINE__ \
		); \
	}
#else
# define IAPWS97_APPROXIMATE
#endif

#define SQ(X) ((X)*(X))

/* Basic math routines, if necesary... */

FREESTEAM_DLL double freesteam_ipow(double x, int n);

#ifdef FREESTEAM_BUILDING_LIB
/* our local ipow implementation */
# define ipow freesteam_ipow
/* 'isnan' function for use with Windows */
# ifdef WIN32
#  include <float.h>
#  define isnan _isnan
# endif
#endif

#endif /* FREESTEAM_COMMON_H */
