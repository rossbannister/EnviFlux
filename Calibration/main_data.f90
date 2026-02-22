module main_data
        
! Description:
! Contains data for INVICAT
!
! Method :
! ------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      12/2014  Initial version (C. Wilson)
!           10/2020  Modified, Ross Bannister
!
! Code Description:
!   Language:           Fortran 90.

implicit none

integer, parameter               :: res=21
integer, parameter               :: ylon=64, ylat=32, ylev=56
integer, parameter               :: nvtot=1
real   , parameter               :: mod_ts=1800.
integer, parameter               :: nmonth=37
integer, parameter               :: nswitch=1, gswitch=0, iswitch=0
integer, parameter               :: ypick=2012, mpick=1, dpick=1.
integer, parameter               :: rsy=2012, rsm=12, rsd=32.
integer, parameter               :: iniinc=1, invicat_rs=1
real   , parameter               :: fems=2.50, fini=1.00
real   , parameter               :: oberr=3.
!real   , parameter               :: ndays=365.+31.+27.
real   , parameter               :: ndays=31.+31.+27.
!real,    parameter               :: ndays = 2.
integer, parameter               :: ns_stop=30
integer, parameter               :: nsmax=120
integer, parameter               :: crashmax=5
integer, parameter               :: fluxmulti=1
integer, parameter               :: extraloop=80
integer, parameter               :: montherr=1
integer, parameter               :: fluxinterp=1
integer, parameter               :: diagb_on=1
integer, parameter               :: deadobs=0
integer, parameter               :: nfields=1
integer, parameter, dimension(1) :: spacecorr = (/0/)
integer, parameter, dimension(1) :: timecorr = (/0/)

integer                          :: ncy_in
integer                          :: gs
integer                          :: nini
integer                          :: nfl
integer, parameter               :: dim_1d   = ylon * ylat * (nmonth * nfields + ylev)
integer                          :: dim_work
integer                          :: nbias

integer                          :: dim_obs, nonobs, gonobs, isnobs, sunobs
character*90, dimension(nfields) :: fluxfile, fluxfilepre, fluxfilepost
character*90                     :: biasfile1, biasfile2
integer                          :: nitr, optim, nsim, isim

        
end module main_data
