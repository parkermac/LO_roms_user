/*
** svn $Id: fennel_var.h 523 2011-01-05 03:21:38Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2011 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Fennel et al. (2006) ecosystem   **
**  model variables that are used in input and output NetCDF files.   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iDetr)')
                idTvar(iDetr)=varid
              CASE ('idTvar(iLDetr)')
                idTvar(iLDetr)=varid
!              CASE ('idTbry(ieast,iSDet)')
!                idTbry(ieast,iSDet)=varid
!              CASE ('idTbry(iwest,iSDet)')
!                idTbry(iwest,iSDet)=varid
!              CAS E ('idTbry(inorth,iSDet)')
!                idTbry(inorth,iSDet)=varid
!              CASE ('idTbry(isouth,iSDet)')
!                idTbry(isouth,iSDet)=varid
# ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid

              CASE ('idTvar(iCaCO3)')
                idTvar(iCaCO3)=varid

#endif
# ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
# endif

# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY   || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR || \
     defined SO_SEMI
              CASE ('idTads(iNO3_)')
                idTads(iNO3_)=varid
              CASE ('idTads(iPhyt)')
                idTads(iPhyt)=varid
              CASE ('idTads(iZoop)')
                idTads(iZoop)=varid
              CASE ('idTads(iDetr)')
                idTads(iDetr)=varid
              CASE ('idTads(iLDetr)')
                idTads(iLDetr)=varid
#  ifdef OXYGEN
              CASE ('idTads(iOxyg)')
                idTads(iOxyg)=varid
#  endif
# endif
              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid
              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid
              CASE ('idTbry(iwest,iDetr)')
                idTbry(iwest,iDetr)=varid
              CASE ('idTbry(ieast,iDetr)')
                idTbry(ieast,iDetr)=varid
              CASE ('idTbry(isouth,iDetr)')
                idTbry(isouth,iDetr)=varid
              CASE ('idTbry(inorth,iDetr)')
                idTbry(inorth,iDetr)=varid
              CASE ('idTbry(iwest,iLDetr)')
                idTbry(iwest,iLDetr)=varid
              CASE ('idTbry(ieast,iLDetr)')
                idTbry(ieast,iLDetr)=varid
              CASE ('idTbry(isouth,iLDetr)')
                idTbry(isouth,iLDetr)=varid
              CASE ('idTbry(inorth,iLDetr)')
                idTbry(inorth,iLDetr)=varid
#ifdef CARBON

              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid

              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid

              CASE ('idTbry(iwest,iCaCO3)')
                idTbry(iwest,iCaCO3)=varid
              CASE ('idTbry(ieast,iCaCO3)')
                idTbry(ieast,iCaCO3)=varid
              CASE ('idTbry(isouth,iCaCO3)')
                idTbry(isouth,iCaCO3)=varid
              CASE ('idTbry(inorth,iCaCO3)')
                idTbry(inorth,iCaCO3)=varid

#endif
# ifdef OXYGEN
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid
# endif
# ifdef DIAGNOSTICS_BIO
#  ifdef OXYGEN
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
              CASE ('iDbio2(iDsed)')
                iDbio2(iDsed)=varid
              CASE ('iDbio3(iDNIT)')
                iDbio3(iDNIT)=varid
#  endif
# ifdef CARBON
!              CASE ('iDbio2(iCOfx)')
!                iDbio2(iCOfx)=varid
!              CASE ('iDbio2(ipCO2)')
!                iDbio2(ipCO2)=varid
# endif
              CASE ('iDbio3(iPPro)')
                iDbio3(iPPro)=varid
              CASE ('iDbio3(iGraz)')
                iDbio3(iGraz)=varid
              CASE ('iDbio3(iRem)')
                iDbio3(iRem)=varid

# endif





/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iPhyt)')
                idRtrc(iPhyt)=varid
              CASE ('idRtrc(iZoop)')
                idRtrc(iZoop)=varid
              CASE ('idRtrc(iDetr)')
                idRtrc(iDetr)=varid
              CASE ('idRtrc(iLDetr)')
                idRtrc(iLDetr)=varid

#ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid

              CASE ('idRtrc(iCaCO3)')
                idRtrc(iCaCO3)=varid

#endif
#ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
#endif




