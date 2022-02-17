/*
** svn $Id: upwelling.h 1054 2021-03-06 19:47:12Z arango $
*******************************************************************************
** Copyright (c) 2002-2021 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Upwelling Test.
**
** Application flag:   Capitalized name of this folder
** Input script:       liveocean.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG

#define DJ_GRADPS
#define TS_DIF2
#define  MIX_GEO_TS
#define SALINITY
#define SOLVE3D

#undef AVERAGES
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV

/* LiveOcean specific choices */
#define NONLIN_EOS
#define MASKING
#define SPHERICAL
/* undef atm forcing for analytical run */
#undef SOLAR_SOURCE
#undef BULK_FLUXES
#undef LONGWAVE_OUT
/*  */
#define DEFLATE
#define HDF5
#define RADIATION_2D
#undef RAMP_TIDES
#define SSH_TIDES
#define UV_TIDES
#define ADD_FSOBC
#define ADD_M2OBC

#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX

#define GLS_MIXING
#if defined GLS_MIXING
# define CANUTO_A
# define N2S2_HORAVG
# define RI_SPLINES
#endif

#if defined BIO_FENNEL  || defined ECOSIM || \
    defined NPZD_POWELL || defined NEMURO
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_SRFLUX
#endif

#if defined NEMURO
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
#endif

#ifdef BIO_FENNEL
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
#endif

#ifdef PERFECT_RESTART
# undef  AVERAGES
# undef  DIAGNOSTICS_BIO
# undef  DIAGNOSTICS_TS
# undef  DIAGNOSTICS_UV
# define OUT_DOUBLE
#endif
