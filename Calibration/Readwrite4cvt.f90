!=================================================================================================
SUBROUTINE read_concentration_fc ( filename, & ! IN  name of netCDF file
                                   fc,       & ! OUT forecast field (last in file)
                                   lons,     & ! OUT (optional) longitudes
                                   lats,     & ! OUT (optional) latitudes
                                   alt,      & ! OUT (optional) altitudes (3d field)
                                   p,        & ! OUT (optional) pressures (3d field)
                                   T )         ! OUT (optional) temperature (3d field)
                                   

! Description:
! Reads in the last concentration field in the specified file.
! Also reads in the longs, lats, level altitudes, pressures, and temperatures if these (optional) arguments are provided.
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      09/2019  New code.  R.N.Bannister
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev


IMPLICIT NONE

INCLUDE 'netcdf.inc'

! Arguments:
CHARACTER(LEN=*),           INTENT(IN)  :: filename
REAL,                       INTENT(OUT) :: fc(1:ylat, 1:ylon, 1:ylev)
REAL,             OPTIONAL, INTENT(OUT) :: lons(1:ylon)
REAL,             OPTIONAL, INTENT(OUT) :: lats(1:ylat)
REAL,             OPTIONAL, INTENT(OUT) :: alt(1:ylat, 1:ylon, 1:ylev)
REAL,             OPTIONAL, INTENT(OUT) :: p(1:ylat, 1:ylon, 1:ylev)
REAL,             OPTIONAL, INTENT(OUT) :: T(1:ylat, 1:ylon, 1:ylev)


! Local variables:
INTEGER                                 :: lon, lat, lev, times, nlev_file
INTEGER                                 :: status, ncid
INTEGER                                 :: dimidtime, dimidsigma
INTEGER                                 :: varidlon, varidlat, varidalt
INTEGER                                 :: varidp, varidT, varidfc
INTEGER                                 :: startA(1), countA(1), startD(4), countD(4)


!Open the netCDF file
!---------------------
status = NF_OPEN(filename, NF_NOWRITE, ncid)
IF ( status .NE. 0 ) THEN
  PRINT*, ' *** Error opening file ***'
  PRINT*, status, NF_STRERROR(status)
  PRINT*, 'FILE :: ', filename
  STOP
ENDIF

!Get the necessary dimension ids
!-------------------------------
status = NF_INQ_DIMID(ncid, 'time', dimidtime)
IF ( status .NE. 0 ) THEN
  PRINT*, ' *** Error obtaining dimid for time ***'
  PRINT*, status, NF_STRERROR(status)
  PRINT*, 'FILE :: ', filename
  STOP
END IF

status = NF_INQ_DIMID(ncid, 'niv', dimidsigma)
IF ( status .NE. 0 ) THEN
  PRINT*, ' *** Error obtaining dimid for niv (sigma level) ***'
  PRINT*, status, NF_STRERROR(status)
  PRINT*, 'FILE :: ', filename
  STOP
END IF


!Get the variable ids
!--------------------
status = NF_INQ_VARID(ncid, 'TR1', varidfc)
IF ( status .NE. 0 ) THEN
  PRINT*, ' *** Error obtaining varid for TR1 (methane forecast) ***'
  PRINT*, status, NF_STRERROR(status)
  PRINT*, 'FILE :: ', filename
  STOP
END IF
IF (PRESENT(lons)) THEN
  status = NF_INQ_VARID(ncid, 'lon', varidlon)
  IF ( status .NE. 0 ) THEN
    PRINT*, ' *** Error obtaining varid for lon ***'
    PRINT*, status, NF_STRERROR(status)
    PRINT*, 'FILE :: ', filename
    STOP
  END IF
END IF
IF (PRESENT(lats)) THEN
  status = NF_INQ_VARID(ncid, 'lat', varidlat)
  IF ( status .NE. 0 ) THEN
    PRINT*, ' *** Error obtaining varid for lat ***'
    PRINT*, status, NF_STRERROR(status)
    PRINT*, 'FILE :: ', filename
    STOP
  END IF
END IF
IF (PRESENT(alt)) THEN
  status = NF_INQ_VARID(ncid, 'g3d', varidalt)
  IF ( status .NE. 0 ) THEN
    PRINT*, ' *** Error obtaining varid for g3d (altitude) ***'
    PRINT*, status, NF_STRERROR(status)
    PRINT*, 'FILE :: ', filename
    STOP
  END IF
END IF
IF (PRESENT(p)) THEN
  status = NF_INQ_VARID(ncid, 'p', varidp)
  IF ( status .NE. 0 ) THEN
    PRINT*, ' *** Error obtaining varid for p (pressure) ***'
    PRINT*, status, NF_STRERROR(status)
    PRINT*, 'FILE :: ', filename
    STOP
  END IF
END IF
IF (PRESENT(T)) THEN
  status = NF_INQ_VARID(ncid, 't3d', varidt)
  IF ( status .NE. 0 ) THEN
    PRINT*, ' *** Error obtaining varid for t3d (temperature) ***'
    PRINT*, status, NF_STRERROR(status)
    PRINT*, 'FILE :: ', filename
    STOP
  END IF
END IF



! Find the number of times present in the file
status = NF_INQ_DIMLEN(ncid, dimidtime, times)
IF ( status .NE. 0 ) THEN
  PRINT*, ' *** Error obtaining the number of times in the file ***'
  PRINT*, status, NF_STRERROR(status)
  PRINT*, 'FILE :: ', filename
  STOP
END IF


! Find the number of vertical levels present in the file
status = NF_INQ_DIMLEN(ncid, dimidsigma, nlev_file)
IF ( status .NE. 0 ) THEN
  PRINT*, ' *** Error obtaining the number of vertical levels in the file ***'
  PRINT*, status, NF_STRERROR(status)
  PRINT*, 'FILE :: ', filename
  STOP
END IF



!Get the main forecast variable
!------------------------------
startD(1) = 1     ! lon
startD(4) = times ! time
countD(1) = ylon  ! lon
countD(2) = 1     ! lat
countD(3) = 1     ! lev
countD(4) = 1     ! time
DO lev = 1, ylev
  DO lat = 1, ylat
    startD(2) = lat
    startD(3) = nlev_file - ylev + lev
    status = NF_GET_VARA_DOUBLE (ncid, varidfc, startD, countD,   &
                                 fc(lat, 1:ylon, lev))
    IF ( status .NE. 0 ) THEN
      PRINT*, ' *** Error reading TR1 (methane forecast) ***'
      PRINT*, 'Level ', lev
      PRINT*, 'Lat   ', lat
      PRINT*, 'Time  ', times
      PRINT*, status, NF_STRERROR(status)
      PRINT*, 'FILE :: ', filename
      STOP
    END IF
  END DO
END DO


!Get the longitudes
!------------------

IF (PRESENT(lons)) THEN
  startA(1) = 1
  countA(1) = ylon
  status = NF_GET_VARA_DOUBLE (ncid, varidlon, startA, countA, lons(1:ylon))
  IF ( status .NE. 0 ) THEN
    PRINT*, ' *** Error reading longitudes ***'
    PRINT*, status, NF_STRERROR(status)
    PRINT*, 'FILE :: ', filename
    STOP
  END IF
END IF

IF (PRESENT(lats)) THEN
  startA(1) = 1
  countA(1) = ylat
  status = NF_GET_VARA_DOUBLE (ncid, varidlat, startA, countA, lats(1:ylat))
  IF ( status .NE. 0 ) THEN
    PRINT*, ' *** Error reading latitudes ***'
    PRINT*, status, NF_STRERROR(status)
    PRINT*, 'FILE :: ', filename
    STOP
  END IF
END IF

IF (PRESENT(lats)) THEN
  startD(1) = 1     ! lon
  startD(4) = times ! time
  countD(1) = ylon  ! lon
  countD(2) = 1     ! lat
  countD(3) = 1     ! lev
  countD(4) = 1     ! time
  DO lev = 1, ylev
    DO lat = 1, ylat
      startD(2) = lat
      startD(3) = nlev_file - ylev + lev
      status = NF_GET_VARA_DOUBLE (ncid, varidalt, startD, countD,   &
                                   alt(lat, 1:ylon, lev))
      IF ( status .NE. 0 ) THEN
        PRINT*, ' *** Error reading g3d (altitudes) ***'
        PRINT*, 'Level ', lev
        PRINT*, 'Lat   ', lat
        PRINT*, 'Time  ', times
        PRINT*, status, NF_STRERROR(status)
        PRINT*, 'FILE :: ', filename
        STOP
      END IF
    END DO
  END DO
END IF

IF (PRESENT(p)) THEN
  startD(1) = 1     ! lon
  startD(4) = times ! time
  countD(1) = ylon  ! lon
  countD(2) = 1     ! lat
  countD(3) = 1     ! lev
  countD(4) = 1     ! time
  DO lev = 1, ylev
    DO lat = 1, ylat
      startD(2) = lat
      startD(3) = nlev_file - ylev + lev
      status = NF_GET_VARA_DOUBLE (ncid, varidp, startD, countD,   &
                                   p(lat, 1:ylon, lev))
      IF ( status .NE. 0 ) THEN
        PRINT*, ' *** Error reading p (pressures) ***'
        PRINT*, 'Level ', lev
        PRINT*, 'Lat   ', lat
        PRINT*, 'Time  ', times
        PRINT*, status, NF_STRERROR(status)
        PRINT*, 'FILE :: ', filename
        STOP
      END IF
    END DO
  END DO
END IF

IF (PRESENT(T)) THEN
  startD(1) = 1     ! lon
  startD(4) = times ! time
  countD(1) = ylon  ! lon
  countD(2) = 1     ! lat
  countD(3) = 1     ! lev
  countD(4) = 1     ! time
  DO lev = 1, ylev
    DO lat = 1, ylat
      startD(2) = lat
      startD(3) = nlev_file - ylev + lev
      status = NF_GET_VARA_DOUBLE (ncid, varidt, startD, countD,   &
                                   T(lat, 1:ylon, lev))
      IF ( status .NE. 0 ) THEN
        PRINT*, ' *** Error reading t3d (temperatures) ***'
        PRINT*, 'Level ', lev
        PRINT*, 'Lat   ', lat
        PRINT*, 'Time  ', times
        PRINT*, status, NF_STRERROR(status)
        PRINT*, 'FILE :: ', filename
        STOP
      END IF
    END DO
  END DO
END IF


!Close the netCDF file
!---------------------
status = NF_CLOSE(ncid)
IF ( status .NE. 0 ) THEN
  PRINT*, ' *** Error closing file ***'
  PRINT*, status, NF_STRERROR(status)
  PRINT*, 'FILE :: ', filename
  STOP
ENDIF


END SUBROUTINE read_concentration_fc






!=================================================================================================
SUBROUTINE Write_time_evolving_fields (filename, Ny, Nx, Nz, Nt, lats, lons, field, fieldname)

!********************************************************
!* Subroutine to write a 4D field                       *
!*                                                      *
!*                                                      *
!* R. Bannister, 11-10-19                               *
!*                                                      *
!********************************************************


IMPLICIT NONE

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE 'netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*), INTENT(IN)  :: filename
INTEGER,          INTENT(IN)  :: Ny, Nx, Nz, Nt
REAL,             INTENT(IN)  :: lats(Ny)
REAL,             INTENT(IN)  :: lons(Nx)
REAL,             INTENT(IN)  :: field(1:Ny, 1:Nx, 1:Nz, 1:Nt)
CHARACTER(LEN=*), INTENT(IN)  :: fieldname


!Declare local variables
!------------------------
INTEGER                       :: ncid, ierr, ddA(1), ddB(4)
INTEGER                       :: dimidx, dimidy, dimidz, dimidt
INTEGER                       :: varidx, varidy, varidz, varidt, varid_field
INTEGER                       :: startA(1), countA(1), startB(4), countB(4)
INTEGER                       :: ierr1, ierr2, ierr3, ierr4
REAL                          :: zs(1:Nz), ts(1:Nt)

INTEGER                       :: z, t

!*****************************************************************************************
!PRINT*, 'Write_time_evolving_fields'
!*****************************************************************************************

  ! Create netCDF file
  !-------------------------------------
  ierr = NF_CREATE(filename, NF_CLOBBER, ncid)
  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error creating file ***'
    PRINT*, filename
    PRINT*, ierr, NF_STRERROR(ierr)
    STOP
    !ELSE
    !PRINT*, 'FILE CREATED'
  ENDIF

  !Define the dimensions
  !---------------------
  ierr1 = NF_DEF_DIM(ncid, 'lat', Ny, dimidy)
  ierr2 = NF_DEF_DIM(ncid, 'lon', Nx, dimidx)
  ierr3 = NF_DEF_DIM(ncid, 'lev', Nz, dimidz)
  ierr4 = NF_DEF_DIM(ncid, 'yr_index', Nt, dimidt)

  IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
     PRINT*, '***Error defining dimension ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     STOP
  !ELSE
  !   PRINT*, 'Dimension ids defined'
  ENDIF

  !Define the variables (include variables giving the dim. values)
  !---------------------------------------------------------------

  ! Dimension variables

  ddA(1) = dimidy
  ierr1  = NF_DEF_VAR(ncid, 'y', NF_DOUBLE, 1, ddA, varidy)
  ddA(1) = dimidx
  ierr2  = NF_DEF_VAR(ncid, 'x', NF_DOUBLE, 1, ddA, varidx)
  ddA(1) = dimidz
  ierr3  = NF_DEF_VAR(ncid, 'z', NF_DOUBLE, 1, ddA, varidz)
  ddA(1) = dimidt
  ierr3  = NF_DEF_VAR(ncid, 't', NF_DOUBLE, 1, ddA, varidt)

  IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
     PRINT*, '***Error defining dimension variable ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     STOP
  !ELSE
  !PRINT*, 'Dimension variable ids defined'
  ENDIF

  ! Main variable

  ddB(1)  = dimidy
  ddB(2)  = dimidx
  ddB(3)  = dimidz
  ddB(4)  = dimidt
  ierr    = NF_DEF_VAR(ncid, fieldname, NF_DOUBLE, 4, ddB, varid_field)

  IF (ierr /= 0) THEN
    PRINT*, '***Error defining main variable id ***'
    PRINT*, 'ierr ',  ierr,  NF_STRERROR(ierr)
    STOP
  ENDIF

  !Change mode of netCDF operation from define to write
  !------------------------------------------------------
  ierr = NF_ENDDEF(ncid)

  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error changing mode of netCDF operation ***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  !ELSE
  !PRINT*, 'Mode Changed'
  END IF

  !---------------------------------------------
  ! Output the values of the dimension variables
  ! --------------------------------------------
  DO z = 1, Nz
    zs(z) = REAL(z)
  END DO
  DO t = 1, Nt
    ts(t) = REAL(t)
  END DO


PRINT *, '-------------------'
PRINT *, 'Lats'
PRINT *, lats(1:Ny)
PRINT *, '-------------------'
PRINT *, 'Lons'
PRINT *, lons(1:Nx)
PRINT *, '-------------------'




  startA(1) = 1
  countA(1) = Ny
  ierr1     = NF_PUT_VARA_DOUBLE(ncid, varidy, startA, countA, lats(1:Ny))
  countA(1) = Nx
  ierr2     = NF_PUT_VARA_DOUBLE(ncid, varidx, startA, countA, lons(1:Nx))
  countA(1) = Nz
  ierr3     = NF_PUT_VARA_DOUBLE(ncid, varidz, startA, countA, zs(1:Nz))
  countA(1) = Nt
  ierr4     = NF_PUT_VARA_DOUBLE(ncid, varidt, startA, countA, ts(1:Nt))


  IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
    PRINT*, '***Error writing dimension data ***'
    PRINT*, 'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*, 'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*, 'ierr3', ierr3, NF_STRERROR(ierr3)
    PRINT*, 'ierr4', ierr4, NF_STRERROR(ierr4)
    STOP
  !ELSE
    !PRINT*, 'Dimension Variables output ok'
  ENDIF

!--------------------------------------------
! Output the values of the main variable
! -------------------------------------------

startB(1) = 1
countB(1) = Ny
startB(2) = 1
countB(2) = Nx
startB(3) = 1
countB(3) = Nz
startB(4) = 1
countB(4) = Nt

ierr      = NF_PUT_VARA_DOUBLE(ncid, varid_field, startB, countB, field(1:Ny, 1:Nx, 1:Nz, 1:Nt))

IF (ierr /= 0) THEN
  PRINT*, '***Error writing main variable data ***'
  PRINT*,'ierr',   ierr,   NF_STRERROR(ierr)
  STOP
!ELSE
  !PRINT*, 'Main data written'
ENDIF

!Close-up the file
!-----------------
ierr = NF_CLOSE(ncid)

IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error closing netCDF file ***'
  PRINT*,'ierr', ierr, NF_STRERROR(ierr)
  PRINT*,'xconv -i ', filename, ' &'
  STOP
!ELSE
  !PRINT*, 'File closed'
ENDIF


END SUBROUTINE Write_time_evolving_fields



!=================================================================================================
SUBROUTINE Write_one_field (filename, Nx, Ny, Nz, field, fieldname, xs, ys, zs)

!********************************************************
!* Subroutine to write a single field                   *
!*                                                      *
!*                                                      *
!* R. Bannister, vn1.4da, 14-12-17                      *
!*                                                      *
!********************************************************


IMPLICIT NONE

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
! INCLUDE '/usr/include/netcdf.inc'
INCLUDE 'netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*), INTENT(IN)  :: filename
INTEGER,          INTENT(IN)  :: Nx, Ny, Nz
REAL,             INTENT(IN)  :: field(1:Nx, 1:Ny, 1:Nz)
CHARACTER(LEN=*), INTENT(IN)  :: fieldname
REAL,             INTENT(IN)  :: xs(1:Nx)
REAL,             INTENT(IN)  :: ys(1:Ny)
REAL,             INTENT(IN)  :: zs(1:Nz)



!Declare local variables
!------------------------
INTEGER                       :: ncid, ierr, ddA(1), ddB(3)
INTEGER                       :: dimidx, dimidy, dimidz
INTEGER                       :: varidx, varidy, varidz, varid_field
INTEGER                       :: startA(1), countA(1), startB(3), countB(3)
INTEGER                       :: ierr1, ierr2, ierr3
INTEGER                       :: x, y, z

!*****************************************************************************************
!PRINT*, 'Write one field'
!*****************************************************************************************

  ! Create netCDF file
  !-------------------------------------
  ierr = NF_CREATE(filename, NF_CLOBBER, ncid)
  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error creating file ***'
    PRINT*, filename
    PRINT*, ierr, NF_STRERROR(ierr)
    STOP
    !ELSE
    !PRINT*, 'FILE CREATED'
  ENDIF

  !Define the dimensions
  !---------------------
  ierr1 = NF_DEF_DIM(ncid, 'x', Nx, dimidx)
  ierr2 = NF_DEF_DIM(ncid, 'y', Ny, dimidy)
  ierr3 = NF_DEF_DIM(ncid, 'z', Nz, dimidz)

  IF ((ierr1 + ierr2 + ierr3) /= 0) THEN
     PRINT*, '***Error defining dimension ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     STOP
  !ELSE
  !   PRINT*, 'Dimension ids defined'
  ENDIF

  !Define the variables (include variables giving the dim. values)
  !---------------------------------------------------------------

  ! Dimension variables

  ddA(1) = dimidx
  ierr1  = NF_DEF_VAR(ncid, 'x', NF_DOUBLE, 1, ddA, varidx)

  ddA(1) = dimidy
  ierr2  = NF_DEF_VAR(ncid, 'y', NF_DOUBLE, 1, ddA, varidy)

  ddA(1) = dimidz
  ierr3  = NF_DEF_VAR(ncid, 'z', NF_DOUBLE, 1, ddA, varidz)

  IF ((ierr1 + ierr2 + ierr3) /= 0) THEN
     PRINT*, '***Error defining dimension variable ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     STOP
  !ELSE
  !PRINT*, 'Dimension variable ids defined'
  ENDIF

  ! Main variable

  ddB(1)  = dimidx
  ddB(2)  = dimidy
  ddB(3)  = dimidz
  ierr    = NF_DEF_VAR(ncid, fieldname, NF_DOUBLE, 3, ddB, varid_field)

  IF (ierr /= 0) THEN
    PRINT*, '***Error defining main variable id ***'
    PRINT*, 'ierr ',  ierr,  NF_STRERROR(ierr)
    STOP
  ENDIF

  !Change mode of netCDF operation from define to write
  !------------------------------------------------------
  ierr = NF_ENDDEF(ncid)

  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error changing mode of netCDF operation ***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  !ELSE
  !PRINT*, 'Mode Changed'
  END IF

  !---------------------------------------------
  ! Output the values of the dimension variables
  ! --------------------------------------------
  startA(1) = 1
  countA(1) = Nx
  ierr1     = NF_PUT_VARA_DOUBLE(ncid, varidx, startA, countA, xs(1:Nx))
  countA(1) = Ny
  ierr2     = NF_PUT_VARA_DOUBLE(ncid, varidy, startA, countA, ys(1:Ny))
  countA(1) = Nz
  ierr3     = NF_PUT_VARA_DOUBLE(ncid, varidz, startA, countA, zs(1:Nz))

  IF ((ierr1 + ierr2 + ierr3) /= 0) THEN
    PRINT*, '***Error writing dimension data ***'
    PRINT*, 'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*, 'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*, 'ierr3', ierr3, NF_STRERROR(ierr3)
    STOP
  !ELSE
    !PRINT*, 'Dimension Variables output ok'
  ENDIF

!--------------------------------------------
! Output the values of the main variable
! -------------------------------------------

startB(1) = 1
countB(1) = Nx
startB(2) = 1
countB(2) = Ny
startB(3) = 1
countB(3) = Nz

ierr      = NF_PUT_VARA_DOUBLE(ncid, varid_field, startB, countB, field(1:Nx, 1:Ny, 1:Nz))

IF (ierr /= 0) THEN
  PRINT*, '***Error writing main variable data ***'
  PRINT*,'ierr',   ierr,   NF_STRERROR(ierr)
  STOP
!ELSE
  !PRINT*, 'Main data written'
ENDIF

!Close-up the file
!-----------------
ierr = NF_CLOSE(ncid)

IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error closing netCDF file ***'
  PRINT*,'ierr', ierr, NF_STRERROR(ierr)
  PRINT*,'xconv -i ', filename, ' &'
  STOP
!ELSE
  !PRINT*, 'File closed'
ENDIF

END SUBROUTINE Write_one_field






!=================================================================================================
SUBROUTINE Write_state_vector (filename, oned, for_adj)

!********************************************************
!* Subroutine to write a state vector (or adjoint)      *
!*                                                      *
!*                                                      *
!* R. Bannister, 13-09-21                               *
!*                                                      *
!********************************************************

USE main_data,    ONLY : ylon, ylat, ylev, nmonth, dim_1d, nfields
USE onedvar_data, ONLY : CVT_HorizCors, CVT_VertCors, CVT_TemporalCors

IMPLICIT NONE

!NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE 'netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*), INTENT(IN)  :: filename         ! Netcdf filename
REAL,             INTENT(IN)  :: oned(1:dim_1d)   ! The state vector/adjoint
INTEGER,          INTENT(IN)  :: for_adj          ! 1 = forward, 2 = adjoint

!Declare local variables
!-----------------------
REAL                          :: flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                          :: tracer(1:ylat, 1:ylon, 1:ylev)
INTEGER                       :: ierr(5), ncid, l
INTEGER                       :: ddA(1), ddC(3), ddD(4)
INTEGER                       :: startA(1), countA(1)
INTEGER                       :: startC(3), countC(3)
INTEGER                       :: startD(4), countD(4)
INTEGER                       :: dimidy, dimidx, dimidz, dimidt, dimidf
INTEGER                       :: varidy, varidx, varidz, varidt, varidf
INTEGER                       :: varidflux, varidtracer
REAL                          :: fields(1:nfields)

!Convert 1-d variables to physical structures
!--------------------------------------------
CALL Physical_linear2multid ( flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                              tracer(1:ylat, 1:ylon, 1:ylev),            &
                              oned(1:dim_1d) )

!Create netCDF file
!------------------
ierr(1) = NF_CREATE(filename, NF_CLOBBER, ncid)
IF ( SUM(ierr(1:1)) /= 0 ) THEN
  PRINT*, ' *** Error creating file ***'
  PRINT*, filename
  DO l = 1, 1
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF

!Define the dimensions
!---------------------
ierr(1) = NF_DEF_DIM(ncid, 'latitude',  ylat,    dimidy)
ierr(2) = NF_DEF_DIM(ncid, 'longitude', ylon,    dimidx)
ierr(3) = NF_DEF_DIM(ncid, 'altitude',  ylev,    dimidz)
ierr(4) = NF_DEF_DIM(ncid, 'month',     nmonth,  dimidt)
ierr(5) = NF_DEF_DIM(ncid, 'field',     nfields, dimidf)
IF ( SUM(ierr(1:5)) /= 0 ) THEN
  PRINT*, ' *** Error defining the dimensions ***'
  DO l = 1, 5
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF

!Define the variables (include variables giving the dim. values)
!---------------------------------------------------------------

! Dimension variables
ddA(1)  = dimidy
ierr(1) = NF_DEF_VAR(ncid, 'latitude',  NF_DOUBLE, 1, ddA, varidy)
ddA(1)  = dimidx
ierr(2) = NF_DEF_VAR(ncid, 'longitude', NF_DOUBLE, 1, ddA, varidx)
ddA(1)  = dimidz
ierr(3) = NF_DEF_VAR(ncid, 'altitude',  NF_DOUBLE, 1, ddA, varidz)
ddA(1)  = dimidt
ierr(4) = NF_DEF_VAR(ncid, 'month',     NF_DOUBLE, 1, ddA, varidt)
ddA(1)  = dimidf
ierr(5) = NF_DEF_VAR(ncid, 'field',     NF_DOUBLE, 1, ddA, varidf)
IF ( SUM(ierr(1:5)) /= 0 ) THEN
  PRINT*, ' *** Error defining the dimension variables ***'
  DO l = 1, 5
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF

! Main variables
ddD(1)  = dimidy
ddD(2)  = dimidx
ddD(3)  = dimidt
ddD(4)  = dimidf
IF (for_adj == 1) THEN
  ierr(1) = NF_DEF_VAR(ncid, 'flux', NF_DOUBLE, 4, ddD, varidflux)
ELSE
  ierr(1) = NF_DEF_VAR(ncid, 'flux_adj', NF_DOUBLE, 4, ddD, varidflux)
END IF
ddC(1)  = dimidy
ddC(2)  = dimidx
ddC(3)  = dimidz
IF (for_adj == 1) THEN
  ierr(2) = NF_DEF_VAR(ncid, 'tracer', NF_DOUBLE, 3, ddC, varidtracer)
ELSE
  ierr(2) = NF_DEF_VAR(ncid, 'tracer_adj', NF_DOUBLE, 3, ddC, varidtracer)
END IF
IF ( SUM(ierr(1:2)) /= 0 ) THEN
  PRINT*, ' *** Error defining the main variables ***'
  DO l = 1, 2
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF


!Change mode of netCDF operation from define to write
!----------------------------------------------------
ierr(1) = NF_ENDDEF(ncid)
IF ( SUM(ierr(1:1)) /= 0 ) THEN
  PRINT*, ' *** Error changing the mode ***'
  DO l = 1, 1
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF

!Output the values of the dimension variables
!--------------------------------------------
startA(1) = 1
countA(1) = ylat
ierr(1)   = NF_PUT_VARA_DOUBLE(ncid, varidy, startA, countA, CVT_HorizCors % TOMCAT_lats(1:ylat))
countA(1) = ylon
ierr(2)   = NF_PUT_VARA_DOUBLE(ncid, varidx, startA, countA, CVT_HorizCors % TOMCAT_longs(1:ylon))
countA(1) = ylev
ierr(3)   = NF_PUT_VARA_DOUBLE(ncid, varidz, startA, countA, CVT_VertCors % alts(1:ylev))
countA(1) = nmonth
ierr(4)   = NF_PUT_VARA_DOUBLE(ncid, varidt, startA, countA, CVT_TemporalCors % months(1:nmonth))
countA(1) = nfields
DO l = 1, nfields
  fields(l) = REAL(l)
END DO
ierr(5)   = NF_PUT_VARA_DOUBLE(ncid, varidf, startA, countA, fields(1:nfields))
IF ( SUM(ierr(1:5)) /= 0 ) THEN
  PRINT*, ' *** Error outputting dimension variables ***'
  DO l = 1, 5
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF

!Output the values of the main variables
!---------------------------------------
startD(1) = 1
countD(1) = ylat
startD(2) = 1
countD(2) = ylon
startD(3) = 1
countD(3) = nmonth
startD(4) = 1
countD(4) = nfields
ierr(1)   = NF_PUT_VARA_DOUBLE(ncid, varidflux, startD, countD, flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
startC(1) = 1
countC(1) = ylat
startC(2) = 1
countC(2) = ylon
startC(3) = 1
countC(3) = ylev
ierr(2)   = NF_PUT_VARA_DOUBLE(ncid, varidtracer, startC, countC, tracer(1:ylat, 1:ylon, 1:ylev))
IF ( SUM(ierr(1:2)) /= 0 ) THEN
  PRINT*, ' *** Error outputting main variables ***'
  DO l = 1, 2
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF

!Close the file
!--------------
ierr(1) = NF_CLOSE(ncid)
IF ( SUM(ierr(1:1)) /= 0 ) THEN
  PRINT*, ' *** Error closing netcdf file ***'
  DO l = 1, 1
    PRINT*, ierr(l), NF_STRERROR(ierr(l))
  END DO
  STOP
ENDIF

END SUBROUTINE Write_state_vector





!=================================================================================================
SUBROUTINE Write_state_vector_control (filename, oned, for_adj)

!********************************************************
!* Subroutine to write a state vector (or adjoint)      *
!*                                                      *
!*                                                      *
!* R. Bannister, 14-09-21                               *
!*                                                      *
!********************************************************


USE onedvar_data, ONLY: n1d, L

USE main_data,    ONLY : ylev, nmonth, nfields
USE onedvar_data, ONLY : CVT_VertCors, CVT_TemporalCors

IMPLICIT NONE

!NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE 'netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*), INTENT(IN)  :: filename         ! Netcdf filename
REAL,             INTENT(IN)  :: oned(1:n1d)      ! The state vector/adjoint
INTEGER,          INTENT(IN)  :: for_adj          ! 1 for forward, 2 for adjoint state

!Declare local variables
!-----------------------
REAL                          :: flux(1:2, 0:L, 0:L, 1:nmonth, 1:nfields)
REAL                          :: tracer(1:2, 0:L, 0:L, 1:ylev)
INTEGER                       :: f, mon, ll, index_1d, mm, lev
INTEGER                       :: ierr(6), ncid, n
INTEGER                       :: ddA(1), ddD(4), ddE(5)
INTEGER                       :: startA(1), countA(1)
INTEGER                       :: startD(4), countD(4)
INTEGER                       :: startE(5), countE(5)
INTEGER                       :: dimidreim, dimidtotwn, dimidmeridwn, dimidvertlev, dimidt, dimidf
INTEGER                       :: varidreim, varidtotwn, varidmeridwn, varidvertlev, varidt, varidf
INTEGER                       :: varidflux, varidtracer
REAL                          :: fields(1:nfields), reim(1:2), wn(0:L), vertlev(1:ylev)

!Convert 1-d variables to physical structures
!--------------------------------------------

index_1d = 0

!Surface flux
! - - - - - -
flux(1:2, 0:L, 0:L, 1:nmonth, 1:nfields) = 0.0
! Loop over each flux type
DO f = 1, nfields
  ! Loop over each month (or timescale, depending on option)
  DO mon = 1, nmonth
    ! Loop over total wavenumber
    DO ll = 0, L
      ! Special case for mm=0 (no Sine contribution)
      ! The cosine term
      index_1d = index_1d + 1
      flux(1,ll,0,mon,f) = oned(index_1d)
      IF (ll > 0) THEN
        ! Loop over remaining meridional wavenumbers
        DO mm = 1, ll
          ! The cosine term
          index_1d = index_1d + 1
          flux(1,ll,mm,mon,f) = oned(index_1d)
          ! The sine term
          index_1d = index_1d + 1
          flux(2,ll,mm,mon,f) = oned(index_1d)
        END DO
      END IF
    END DO
  END DO
END DO

! Initial tracer
!  - - - - - - -
tracer(1:2, 0:L, 0:L, 1:ylev) = 0.0
! Loop over each vertical level
DO lev = 1, ylev
  ! Loop over total wavenumber
  DO ll = 0, L
    ! Special case for mm=0 (no Sine contribution)
    ! The cosine term
    index_1d = index_1d + 1
    tracer(1,ll,0,lev) = oned(index_1d)
    IF (ll > 0) THEN
      ! Loop over remaining meridional wavenumbers
      DO mm = 1, ll
        ! The cosine term
        index_1d = index_1d + 1
        tracer(1,ll,mm,lev) = oned(index_1d)
        ! The sine term
        index_1d = index_1d + 1
        tracer(2,ll,mm,lev) = oned(index_1d)
      END DO
    END IF
  END DO
END DO


! Create the dimension indices
DO n = 1, 2
  reim(n) = REAL(n)
END DO
DO n = 0, L
  wn(n) = REAL(n)
END DO
DO n = 1, nfields
  fields(n) = REAL(n)
END DO
DO n = 1, ylev
  vertlev(n) = REAL(n)
END DO

!Create netCDF file
!------------------
ierr(1) = NF_CREATE(filename, NF_CLOBBER, ncid)
IF ( SUM(ierr(1:1)) /= 0 ) THEN
  PRINT*, ' *** Error creating file ***'
  PRINT*, filename
  DO n = 1, 1
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF

!Define the dimensions
!---------------------
ierr(1) = NF_DEF_DIM(ncid, 'real_imag', 2,       dimidreim)
ierr(2) = NF_DEF_DIM(ncid, 'total_wn',  L+1,     dimidtotwn)
ierr(3) = NF_DEF_DIM(ncid, 'medid_wn',  L+1,     dimidmeridwn)
ierr(4) = NF_DEF_DIM(ncid, 'vert_mode', ylev,    dimidvertlev)
ierr(5) = NF_DEF_DIM(ncid, 'month',     nmonth,  dimidt)
ierr(6) = NF_DEF_DIM(ncid, 'field',     nfields, dimidf)
IF ( SUM(ierr(1:6)) /= 0 ) THEN
  PRINT*, ' *** Error defining the dimensions ***'
  DO n = 1, 6
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF

!Define the variables (include variables giving the dim. values)
!---------------------------------------------------------------

! Dimension variables
ddA(1)  = dimidreim
ierr(1) = NF_DEF_VAR(ncid, 'real_imag', NF_DOUBLE, 1, ddA, varidreim)
ddA(1)  = dimidtotwn
ierr(2) = NF_DEF_VAR(ncid, 'totwn',     NF_DOUBLE, 1, ddA, varidtotwn)
ddA(1)  = dimidmeridwn
ierr(3) = NF_DEF_VAR(ncid, 'meridwn',   NF_DOUBLE, 1, ddA, varidmeridwn)
ddA(1)  = dimidvertlev
ierr(4) = NF_DEF_VAR(ncid, 'vertlev',   NF_DOUBLE, 1, ddA, varidvertlev)
ddA(1)  = dimidt
ierr(5) = NF_DEF_VAR(ncid, 'month',     NF_DOUBLE, 1, ddA, varidt)
ddA(1)  = dimidf
ierr(6) = NF_DEF_VAR(ncid, 'field',     NF_DOUBLE, 1, ddA, varidf)
IF ( SUM(ierr(1:6)) /= 0 ) THEN
  PRINT*, ' *** Error defining the dimension variables ***'
  DO n = 1, 6
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF

! Main variables
ddE(1)  = dimidreim
ddE(2)  = dimidtotwn
ddE(3)  = dimidmeridwn
ddE(4)  = dimidt
ddE(5)  = dimidf
IF (for_adj == 1) THEN
  ierr(1) = NF_DEF_VAR(ncid, 'flux_con', NF_DOUBLE, 5, ddE, varidflux)
ELSE
  ierr(1) = NF_DEF_VAR(ncid, 'flux_con_adj', NF_DOUBLE, 5, ddE, varidflux)
END IF
ddD(1)  = dimidreim
ddD(2)  = dimidtotwn
ddD(3)  = dimidmeridwn
ddD(4)  = dimidvertlev
IF (for_adj == 1) THEN
  ierr(2) = NF_DEF_VAR(ncid, 'tracer_con', NF_DOUBLE, 4, ddD, varidtracer)
ELSE
  ierr(2) = NF_DEF_VAR(ncid, 'tracer_con_adj', NF_DOUBLE, 4, ddD, varidtracer)
END IF
IF ( SUM(ierr(1:2)) /= 0 ) THEN
  PRINT*, ' *** Error defining the main variables ***'
  DO n = 1, 2
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF


!Change mode of netCDF operation from define to write
!----------------------------------------------------
ierr(1) = NF_ENDDEF(ncid)
IF ( SUM(ierr(1:1)) /= 0 ) THEN
  PRINT*, ' *** Error changing the mode ***'
  DO n = 1, 1
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF

!Output the values of the dimension variables
!--------------------------------------------
startA(1) = 1
countA(1) = 2
ierr(1)   = NF_PUT_VARA_DOUBLE(ncid, varidreim,    startA, countA, reim(1:2))
countA(1) = L+1
ierr(2)   = NF_PUT_VARA_DOUBLE(ncid, varidtotwn,   startA, countA, wn(0:L))
countA(1) = L+1
ierr(3)   = NF_PUT_VARA_DOUBLE(ncid, varidmeridwn, startA, countA, wn(0:L))
countA(1) = ylev
ierr(4)   = NF_PUT_VARA_DOUBLE(ncid, varidvertlev, startA, countA, CVT_VertCors % alts(1:ylev))
countA(1) = nmonth
ierr(5)   = NF_PUT_VARA_DOUBLE(ncid, varidt,       startA, countA, CVT_TemporalCors % months(1:nmonth))
countA(1) = nfields
ierr(6)   = NF_PUT_VARA_DOUBLE(ncid, varidf,       startA, countA, fields(1:nfields))
IF ( SUM(ierr(1:6)) /= 0 ) THEN
  PRINT*, ' *** Error outputting dimension variables ***'
  DO n = 1, 6
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF

!Output the values of the main variables
!---------------------------------------
startE(1) = 1
countE(1) = 2
startE(2) = 1
countE(2) = L+1
startE(3) = 1
countE(3) = L+1
startE(4) = 1
countE(4) = nmonth
startE(5) = 1
countE(5) = nfields
ierr(1)   = NF_PUT_VARA_DOUBLE(ncid, varidflux, startE, countE, flux(1:2, 0:L, 0:L, 1:nmonth, 1:nfields))
startD(1) = 1
countD(1) = 2
startD(2) = 1
countD(2) = L+1
startD(3) = 1
countD(3) = L+1
startD(4) = 1
countD(4) = ylev
ierr(2)   = NF_PUT_VARA_DOUBLE(ncid, varidtracer, startD, countD, tracer(1:2, 0:L, 0:L, 1:ylev))
IF ( SUM(ierr(1:2)) /= 0 ) THEN
  PRINT*, ' *** Error outputting main variables ***'
  DO n = 1, 2
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF

!Close the file
!--------------
ierr(1) = NF_CLOSE(ncid)
IF ( SUM(ierr(1:1)) /= 0 ) THEN
  PRINT*, ' *** Error closing netcdf file ***'
  DO n = 1, 1
    PRINT*, ierr(n), NF_STRERROR(ierr(n))
  END DO
  STOP
ENDIF

END SUBROUTINE Write_state_vector_control




!=================================================================================================
SUBROUTINE CVT_matrices_output (filename,              &
                           ! ***** Longitudes, latitudes, and heights
                                TOMCAT_longs,          & ! Longitudes of TOMCAT grid
                                TOMCAT_lats,           & ! Latitudes of TOMCAT grid
                                TOMCAT_alts,           & ! Altitudes of TOMCAT grid
                                SHtools_longs,         & ! Longitudes of SHtools grid
                                SHtools_lats,          & ! Latitudes of SHtools grid
                                months,                & ! Month indices
                                GaussianCoLats,        & ! Gaussian co-latitudes
                           ! ***** Spatial standard deviations
                                FracOfPrior_tracer,    & ! Fraction of a-priori for stddev tracer
                                std_tracer,            & ! Standard deviation of tracer (use if above is 0.0)
                                FracOfPrior_flux,      & ! Fraction of a-priori for stddev flux
                                std_flux,              & ! Standard deviation of flux (use if above is 0.0)
                           ! ***** Horizontal transform for tracer
                                horiztracer_cor_tpe,   & ! Imposed horizontal cor fn type for tracer
                                lengthscale_tracer,    & ! Imposed horizontal cor lengthscales
                                horizcors_tracer,      & ! Imposed horizontal correlation fn for tracer
                                var_hspec_tracer,      & ! Horizontal variances of tracer
                                std_hspec_tracer,      & ! Horizontal standard deviations of tracer
                           ! ***** Horizontal transform for flux
                                horizflux_cor_tpe,     & ! Imposed horizontal cor fn type for flux
                                lengthscale_flux,      & ! Imposed horizontal cor lengthscales
                                horizcors_flux,        & ! Imposed horizontal correlation fn for flux
                                var_hspec_flux,        & ! Horizontal variances of flux
                                std_hspec_flux,        & ! Horizontal standard deviations of flux
                           ! ***** Vertical transform for tracer
                                ForceVCorr,            & ! Switch to set if vertical correlations are enforced
                                vert_covs,             & ! Option to switch dependency of vert eigenvals
                                glob_av_vert_cov,      & ! Global av vertical covariance matrix
                                vert_eigenvec_tracer,  & ! Global vertical eigenvectors
                                var_vspec_tracer,      & ! Vertical variances
                                std_vspec_tracer,      & ! Vertical standard deviations
                                vert_adjust_tracer,    & ! Adjustment for vertical transform
                           ! ***** Temporal transform for flux
                                temporal_covs,         & ! Option for temporal flux covariances
                                temporal_cor_tpe,      & ! Imposed temporal correlation fn type
                                timescale_flux,        & ! Imposed flux correlation timescales
                                temporalcors_flux,     & ! Imposed temporal correlation fn
                                temp_cor_matrix,       & ! Imposed temporal correlation matrix
                                temp_eigenvec_flux,    & ! Eigenvectors of temporal cor fns
                                var_tspec_flux,        & ! Eigenvalues of temporal cor fn
                                std_tspec_flux )         ! Square-root of eigenvalues

!********************************************************
!* Subroutine to write the CVT matrices                 *
!*                                                      *
!*                                                      *
!* R. Bannister, 21-11-19                               *
!* Tweaks, R. Bannister, July 2021                      *
!*                                                      *
!********************************************************

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields

USE onedvar_data, ONLY: L


IMPLICIT NONE

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
! INCLUDE '/usr/include/netcdf.inc'
INCLUDE 'netcdf.inc'

!Declare parameters
!------------------
! ***** Longitudes, latitudes, and heights
CHARACTER(LEN=*), INTENT(IN)  :: filename
REAL,             INTENT(IN)  :: TOMCAT_longs(1:ylon)
REAL,             INTENT(IN)  :: TOMCAT_lats(1:ylat)
REAL,             INTENT(IN)  :: TOMCAT_alts(1:ylev)
REAL,             INTENT(IN)  :: SHtools_longs(1:2*L+1)
REAL,             INTENT(IN)  :: SHtools_lats(1:L+1)
REAL,             INTENT(IN)  :: months(1:nmonth)
REAL,             INTENT(IN)  :: GaussianCoLats(1:L+1)

! ***** Spatial standard deviations
REAL,             INTENT(IN)  :: FracOfPrior_tracer
REAL,             INTENT(IN)  :: std_tracer(1:ylat, 1:ylon, 1:ylev)
REAL,             INTENT(IN)  :: FracOfPrior_flux
REAL,             INTENT(IN)  :: std_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
! ***** Horizontal transform for tracer
INTEGER,          INTENT(IN)  :: horiztracer_cor_tpe
REAL,             INTENT(IN)  :: lengthscale_tracer(1:ylev)
REAL,             INTENT(IN)  :: horizcors_tracer(1:L+1, 1:ylev)
REAL,             INTENT(IN)  :: var_hspec_tracer(0:L, 1:ylev)
REAL,             INTENT(IN)  :: std_hspec_tracer(0:L, 1:ylev)
! ***** Horizontal transform for flux
INTEGER,          INTENT(IN)  :: horizflux_cor_tpe(1:nfields)
REAL,             INTENT(IN)  :: lengthscale_flux(1:nmonth, 1:nfields)
REAL,             INTENT(IN)  :: horizcors_flux(0:L, 1:nmonth, 1:nfields)
REAL,             INTENT(IN)  :: var_hspec_flux(0:L, 1:nmonth, 1:nfields)
REAL,             INTENT(IN)  :: std_hspec_flux(0:L, 1:nmonth, 1:nfields)

! ***** Vertical transform for tracer
LOGICAL,          INTENT(IN)  :: ForceVCorr
INTEGER,          INTENT(IN)  :: vert_covs
REAL,             INTENT(IN)  :: glob_av_vert_cov(1:ylev, 1:ylev)
REAL,             INTENT(IN)  :: vert_eigenvec_tracer(1:ylev, 1:ylev)
REAL,             INTENT(IN)  :: var_vspec_tracer(1:ylev, 1:ylat, 1:ylon)
REAL,             INTENT(IN)  :: std_vspec_tracer(1:ylev, 1:ylat, 1:ylon)
REAL,             INTENT(IN)  :: vert_adjust_tracer(1:ylev, 1:ylat, 1:ylon)

! ***** Temporal transform for flux
INTEGER,          INTENT(IN)  :: temporal_covs(1:nfields)
INTEGER,          INTENT(IN)  :: temporal_cor_tpe(1:nfields)
REAL,             INTENT(IN)  :: timescale_flux(1:nfields)
REAL,             INTENT(IN)  :: temporalcors_flux(1:nmonth, 1:nfields)
REAL,             INTENT(IN)  :: temp_cor_matrix(1:nmonth, 1:nmonth, 1:nfields)
REAL,             INTENT(IN)  :: temp_eigenvec_flux(1:nmonth, 1:nmonth, 1:nfields)
REAL,             INTENT(IN)  :: var_tspec_flux(1:nmonth, 1:nfields)
REAL,             INTENT(IN)  :: std_tspec_flux(1:nmonth, 1:nfields)


!Declare local variables
!------------------------
INTEGER                       :: ncid, ddA(1), ddB(2), ddC(3), ddD(4)
INTEGER                       :: dimidx_tomcat, dimidy_tomcat, dimidz_tomcat, dimidx_sh, dimidy_sh
INTEGER                       :: dimid_vertmode, dimid_month, dimid_timescale, dimid_field, dimid_totwn
INTEGER                       :: dimid_options, dimid_colats
INTEGER                       :: varidx_tomcat, varidy_tomcat, varidz_tomcat, varidx_sh, varidy_sh, varid_colats
INTEGER                       :: varid_vertmode, varid_month, varid_timescale, varid_field, varid_totwn
INTEGER                       :: varid_FracOfPrior_tracer, varid_FracOfPrior_flux
INTEGER                       :: varid_std_tracer, varid_vert_adjust_tracer, varid_vert_eigenvec_tracer
INTEGER                       :: varid_var_vspec_tracer, varid_std_vspec_tracer
INTEGER                       :: varid_var_hspec_tracer, varid_std_hspec_tracer, varid_std_flux, varid_timescale_flux
INTEGER                       :: varid_temporal_cor_tpe, varid_temp_eigenvec_flux, varid_horizcors_flux
INTEGER                       :: varid_horizflux_cor_tpe, varid_lengthscale_flux, varid_temporalcors_flux
INTEGER                       :: varid_var_tspec_flux, varid_std_tspec_flux, varid_var_hspec_flux, varid_std_hspec_flux
INTEGER                       :: varid_temporal_covs, varid_vert_covs, varid_options, varid_ForceVCorr
INTEGER                       :: varid_horiztracer_cor_tpe, varid_lengthscale_tracer, varid_horizcors_tracer
INTEGER                       :: varid_glob_av_vert_cov, varid_temp_cor_matrix
INTEGER                       :: startA(1), countA(1), startB(2), countB(2)
INTEGER                       :: startC(3), countC(3), startD(4), countD(4)
INTEGER                       :: ierr(1:29), i

INTEGER                       :: x, y, z, t
REAL                          :: vmodes(1:ylev), tscales(1:nmonth), totwn(0:L), field(1:nfields)
REAL                          :: temp(1:1)
INTEGER                       :: options(1:1)


!*****************************************************************************************
PRINT*, 'Writing CVT data'
!*****************************************************************************************

! Create netCDF file
!-------------------------------------
ierr(1) = NF_CREATE(filename, NF_CLOBBER, ncid)
DO i = 1, 1
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error opening netcdf file ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO

!Define the dimensions
!---------------------
ierr(1) = NF_DEF_DIM(ncid, 'x_tomcat',       ylon,    dimidx_tomcat)
ierr(2) = NF_DEF_DIM(ncid, 'y_tomcat',       ylat,    dimidy_tomcat)
ierr(3) = NF_DEF_DIM(ncid, 'z_tomcat',       ylev,    dimidz_tomcat)
ierr(4) = NF_DEF_DIM(ncid, 'x_sh',           2*L+1,   dimidx_sh)
ierr(5) = NF_DEF_DIM(ncid, 'y_sh',           L+1,     dimidy_sh)
ierr(6) = NF_DEF_DIM(ncid, 'vertmode',       ylev,    dimid_vertmode)
ierr(7) = NF_DEF_DIM(ncid, 'month',          nmonth,  dimid_month)
ierr(8) = NF_DEF_DIM(ncid, 'timescale',      nmonth,  dimid_timescale)
ierr(9) = NF_DEF_DIM(ncid, 'fluxfield',      nfields, dimid_field)
ierr(10)= NF_DEF_DIM(ncid, 'tot_wn',         L+1,     dimid_totwn)
ierr(11)= NF_DEF_DIM(ncid, 'options',        1,       dimid_options)
ierr(12)= NF_DEF_DIM(ncid, 'GaussianCoLats', L+1,     dimid_colats)

DO i = 1, 12
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error defining dimension id ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


!Define the variables (include variables giving the dim. values)
!---------------------------------------------------------------

! Dimension variables

ddA(1)   = dimidx_tomcat
ierr(1)  = NF_DEF_VAR(ncid, 'x_tomcat',       NF_DOUBLE, 1, ddA, varidx_tomcat)
ddA(1)   = dimidy_tomcat
ierr(2)  = NF_DEF_VAR(ncid, 'y_tomcat',       NF_DOUBLE, 1, ddA, varidy_tomcat)
ddA(1)   = dimidz_tomcat
ierr(3)  = NF_DEF_VAR(ncid, 'z_tomcat',       NF_DOUBLE, 1, ddA, varidz_tomcat)
ddA(1)   = dimidx_sh
ierr(4)  = NF_DEF_VAR(ncid, 'x_sh',           NF_DOUBLE, 1, ddA, varidx_sh)
ddA(1)   = dimidy_sh
ierr(5)  = NF_DEF_VAR(ncid, 'y_sh',           NF_DOUBLE, 1, ddA, varidy_sh)
ddA(1)   = dimid_vertmode
ierr(6)  = NF_DEF_VAR(ncid, 'vertmode',       NF_DOUBLE, 1, ddA, varid_vertmode)
ddA(1)   = dimid_month
ierr(7)  = NF_DEF_VAR(ncid, 'month',          NF_DOUBLE, 1, ddA, varid_month)
ddA(1)   = dimid_timescale
ierr(8)  = NF_DEF_VAR(ncid, 'timescale',      NF_DOUBLE, 1, ddA, varid_timescale)
ddA(1)   = dimid_field
ierr(9)  = NF_DEF_VAR(ncid, 'fluxfield',      NF_DOUBLE, 1, ddA, varid_field)
ddA(1)   = dimid_totwn
ierr(10) = NF_DEF_VAR(ncid, 'totwn',          NF_DOUBLE, 1, ddA, varid_totwn)
ddA(1)   = dimid_options
ierr(11) = NF_DEF_VAR(ncid, 'options',        NF_INT,    1, ddA, varid_options)
ddA(1)   = dimid_colats
ierr(12) = NF_DEF_VAR(ncid, 'GaussianCoLats', NF_DOUBLE, 1, ddA, varid_colats)


DO i = 1, 12
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error defining dimension variable ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


! Main variables

ddA(1)  = dimid_options
ierr(1) = NF_DEF_VAR(ncid, 'FracOfPrior_tracer', NF_DOUBLE, 1, ddA, varid_FracOfPrior_tracer)

ddA(1)  = dimid_options
ierr(2) = NF_DEF_VAR(ncid, 'FracOfPrior_flux', NF_DOUBLE, 1, ddA, varid_FracOfPrior_flux)

ddC(1)  = dimidy_tomcat
ddC(2)  = dimidx_tomcat
ddC(3)  = dimidz_tomcat
ierr(3) = NF_DEF_VAR(ncid, 'std_tracer', NF_DOUBLE, 3, ddC, varid_std_tracer)

ddD(1)  = dimidy_tomcat
ddD(2)  = dimidx_tomcat
ddD(3)  = dimid_month
ddD(4)  = dimid_field
ierr(4) = NF_DEF_VAR(ncid, 'std_flux', NF_DOUBLE, 4, ddD, varid_std_flux)

ddA(1)  = dimid_options
ierr(5) = NF_DEF_VAR(ncid, 'horiztracer_cor_tpe', NF_INT, 1, ddA, varid_horiztracer_cor_tpe)

ddA(1)  = dimidz_tomcat
ierr(6) = NF_DEF_VAR(ncid, 'lengthscale_tracer', NF_DOUBLE, 1, ddA, varid_lengthscale_tracer)

ddB(1)  = dimidy_sh
ddB(2)  = dimidz_tomcat
ierr(7) = NF_DEF_VAR(ncid, 'horizcors_tracer', NF_DOUBLE, 2, ddB, varid_horizcors_tracer)

ddB(1)  = dimid_totwn
ddB(2)  = dimidz_tomcat
ierr(8) = NF_DEF_VAR(ncid, 'var_hspec_tracer', NF_DOUBLE, 2, ddB, varid_var_hspec_tracer)

ddB(1)  = dimid_totwn
ddB(2)  = dimidz_tomcat
ierr(9) = NF_DEF_VAR(ncid, 'std_hspec_tracer', NF_DOUBLE, 2, ddB, varid_std_hspec_tracer)

ddA(1)  = dimid_field
ierr(10)= NF_DEF_VAR(ncid, 'horizflux_cor_tpe', NF_INT, 1, ddA, varid_horizflux_cor_tpe)

ddB(1)  = dimid_month  ! Actually dimid_timscale for temporal_covs(field) = 1
ddB(2)  = dimid_field
ierr(11)= NF_DEF_VAR(ncid, 'lengthscale_flux', NF_DOUBLE, 2, ddB, varid_lengthscale_flux)

ddC(1)  = dimidy_sh
ddC(2)  = dimid_month  ! Actually dimid_timscale for temporal_covs(field) = 1
ddC(3)  = dimid_field
ierr(12)= NF_DEF_VAR(ncid, 'horizcors_flux', NF_DOUBLE, 3, ddC, varid_horizcors_flux)

ddC(1)  = dimid_totwn
ddC(2)  = dimid_month  ! Actually dimid_timscale for temporal_covs(field) = 1
ddC(3)  = dimid_field
ierr(13)= NF_DEF_VAR(ncid, 'var_hspec_flux', NF_DOUBLE, 3, ddC, varid_var_hspec_flux)

ddC(1)  = dimid_totwn
ddC(2)  = dimid_month  ! Actually dimid_timscale for temporal_covs(field) = 1
ddC(3)  = dimid_field
ierr(14)= NF_DEF_VAR(ncid, 'std_hspec_flux', NF_DOUBLE, 3, ddC, varid_std_hspec_flux)

ddA(1)  = dimid_options
ierr(15)= NF_DEF_VAR(ncid, 'ForceVCorr', NF_INT, 1, ddA, varid_ForceVCorr)

ddA(1)  = dimid_options
ierr(16)= NF_DEF_VAR(ncid, 'vert_covs', NF_INT, 1, ddA, varid_vert_covs)

ddB(1)  = dimidz_tomcat
ddB(2)  = dimidz_tomcat
ierr(17)= NF_DEF_VAR(ncid, 'glob_av_vert_cov', NF_DOUBLE, 2, ddB, varid_glob_av_vert_cov)

ddB(1)  = dimidz_tomcat
ddB(2)  = dimid_vertmode
ierr(18)= NF_DEF_VAR(ncid, 'vert_eigenvec_tracer', NF_DOUBLE, 2, ddB, varid_vert_eigenvec_tracer)

ddC(1)  = dimid_vertmode
ddC(2)  = dimidy_tomcat
ddC(3)  = dimidx_tomcat
ierr(19)= NF_DEF_VAR(ncid, 'var_vspec_tracer', NF_DOUBLE, 3, ddC, varid_var_vspec_tracer)

ddC(1)  = dimid_vertmode
ddC(2)  = dimidy_tomcat
ddC(3)  = dimidx_tomcat
ierr(20)= NF_DEF_VAR(ncid, 'std_vspec_tracer', NF_DOUBLE, 3, ddC, varid_std_vspec_tracer)

ddC(1)  = dimidz_tomcat
ddC(2)  = dimidy_tomcat
ddC(3)  = dimidx_tomcat
ierr(21)= NF_DEF_VAR(ncid, 'vert_adjust_tracer', NF_DOUBLE, 3, ddC, varid_vert_adjust_tracer)

ddA(1)  = dimid_field
ierr(22)= NF_DEF_VAR(ncid, 'temporal_covs', NF_INT, 1, ddA, varid_temporal_covs)

ddA(1)  = dimid_field
ierr(23)= NF_DEF_VAR(ncid, 'temporal_cor_tpe', NF_INT, 1, ddA, varid_temporal_cor_tpe)

ddA(1)  = dimid_field
ierr(24)= NF_DEF_VAR(ncid, 'timescale_flux', NF_DOUBLE, 1, ddA, varid_timescale_flux)

ddB(1)  = dimid_month
ddB(2)  = dimid_field
ierr(25)= NF_DEF_VAR(ncid, 'temporalcors_flux', NF_DOUBLE, 2, ddB, varid_temporalcors_flux)

ddC(1)  = dimid_month
ddC(2)  = dimid_month
ddC(3)  = dimid_field
ierr(26)= NF_DEF_VAR(ncid, 'temp_cor_matrix', NF_DOUBLE, 3, ddC, varid_temp_cor_matrix)

ddC(1)  = dimid_month
ddC(2)  = dimid_timescale
ddC(3)  = dimid_field
ierr(27)= NF_DEF_VAR(ncid, 'temp_eigenvec_flux', NF_DOUBLE, 3, ddC, varid_temp_eigenvec_flux)

ddB(1)  = dimid_month  ! Actually dimid_timscale for temporal_covs(field) = 1
ddB(2)  = dimid_field
ierr(28)= NF_DEF_VAR(ncid, 'var_tspec_flux', NF_DOUBLE, 2, ddB, varid_var_tspec_flux)

ddB(1)  = dimid_month  ! Actually dimid_timscale for temporal_covs(field) = 1
ddB(2)  = dimid_field
ierr(29)= NF_DEF_VAR(ncid, 'std_tspec_flux', NF_DOUBLE, 2, ddB, varid_std_tspec_flux)

DO i = 1, 29
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error defining main variable ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO



!Change mode of netCDF operation from define to write
!------------------------------------------------------
ierr(1) = NF_ENDDEF(ncid)

DO i = 1, 1
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error changing mode ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO

!---------------------------------------------
! Output the values of the dimension variables
! --------------------------------------------

startA(1) = 1
countA(1) = ylon
ierr(1)   = NF_PUT_VARA_DOUBLE(ncid, varidx_tomcat, startA, countA, TOMCAT_longs(1:ylon))
countA(1) = ylat
ierr(2)   = NF_PUT_VARA_DOUBLE(ncid, varidy_tomcat, startA, countA, TOMCAT_lats(1:ylat))
countA(1) = ylev
ierr(3)   = NF_PUT_VARA_DOUBLE(ncid, varidz_tomcat, startA, countA, TOMCAT_alts(1:ylev))

countA(1) = 2*L+1
ierr(4)   = NF_PUT_VARA_DOUBLE(ncid, varidx_sh, startA, countA, SHtools_longs(1:2*L+1))
countA(1) = L+1
ierr(5)   = NF_PUT_VARA_DOUBLE(ncid, varidy_sh, startA, countA, SHtools_lats(1:L+1))

! Set-up other dimensions
DO z = 1, ylev
  vmodes(z) = REAL(z)
END DO
DO x = 0, L
  totwn(x) = REAL(x)
END DO
DO t = 1, nmonth
  tscales(t) = REAL(t)
END DO
DO t = 1, nfields
  field(t) = REAL(t)
END DO
options(1) = 0

countA(1) = ylev
ierr(6)   = NF_PUT_VARA_DOUBLE(ncid, varid_vertmode, startA, countA, vmodes(1:ylev))
countA(1) = nmonth
ierr(7)   = NF_PUT_VARA_DOUBLE(ncid, varid_month, startA, countA, months(1:nmonth))
countA(1) = nmonth
ierr(8)   = NF_PUT_VARA_DOUBLE(ncid, varid_timescale, startA, countA, tscales(1:nmonth))
countA(1) = L+1
ierr(9)   = NF_PUT_VARA_DOUBLE(ncid, varid_totwn, startA, countA, totwn(0:L))
countA(1) = nfields
ierr(10)  = NF_PUT_VARA_DOUBLE(ncid, varid_field, startA, countA, field(1:nfields))
countA(1) = 1
ierr(11)  = NF_PUT_VARA_INT(ncid, varid_options, startA, countA, options(1:1))
countA(1) = L+1
ierr(12)  = NF_PUT_VARA_DOUBLE(ncid, varid_colats, startA, countA, GaussianCoLats(1:L+1))

DO i = 1, 12
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error outputting dimension variables ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


!--------------------------------------------
! Output the values of the main variables
! -------------------------------------------

startA(1) = 1
countA(1) = 1
temp(1)   = FracOfPrior_tracer
ierr(1)   = NF_PUT_VARA_DOUBLE(ncid, varid_FracOfPrior_tracer, startA, countA, temp)

startA(1) = 1
countA(1) = 1
temp(1)   = FracOfPrior_flux
ierr(2)   = NF_PUT_VARA_DOUBLE(ncid, varid_FracOfPrior_flux, startA, countA, temp)

startC(1) = 1
countC(1) = ylat
startC(2) = 1
countC(2) = ylon
startC(3) = 1
countC(3) = ylev
ierr(3)   = NF_PUT_VARA_DOUBLE(ncid, varid_std_tracer, startC, countC, std_tracer(1:ylat, 1:ylon, 1:ylev))

startD(1) = 1
countD(1) = ylat
startD(2) = 1
countD(2) = ylon
startD(3) = 1
countD(3) = nmonth
startD(4) = 1
countD(4) = nfields
ierr(4)   = NF_PUT_VARA_DOUBLE(ncid, varid_std_flux, startD, countD, std_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))

startA(1) = 1
countA(1) = 1
options(1)= horiztracer_cor_tpe
ierr(5)   = NF_PUT_VARA_INT(ncid, varid_horiztracer_cor_tpe, startA, countA, options)

startA(1) = 1
countA(1) = ylev
ierr(6)   = NF_PUT_VARA_DOUBLE(ncid, varid_lengthscale_tracer, startA, countA, &
                               lengthscale_tracer(1:ylev))

startB(1) = 1
countB(1) = L+1
startB(2) = 1
countB(2) = ylev
ierr(7)   = NF_PUT_VARA_DOUBLE(ncid, varid_horizcors_tracer, startB, countB, horizcors_tracer(1:L+1, 1:ylev))

startB(1) = 1
countB(1) = L+1
startB(2) = 1
countB(2) = ylev
ierr(8)   = NF_PUT_VARA_DOUBLE(ncid, varid_var_hspec_tracer, startB, countB, var_hspec_tracer(0:L, 1:ylev))

startB(1) = 1
countB(1) = L+1
startB(2) = 1
countB(2) = ylev
ierr(9)   = NF_PUT_VARA_DOUBLE(ncid, varid_std_hspec_tracer, startB, countB, std_hspec_tracer(0:L, 1:ylev))

startA(1) = 1
countA(1) = nfields
ierr(10)  = NF_PUT_VARA_INT(ncid, varid_horizflux_cor_tpe, startA, countA, &
                            horizflux_cor_tpe(1:nfields))

startB(1) = 1
countB(1) = nmonth
startB(2) = 1
countB(2) = nfields
ierr(11)  = NF_PUT_VARA_DOUBLE(ncid, varid_lengthscale_flux, startB, countB, &
                               lengthscale_flux(1:nmonth, 1:nfields))

startC(1) = 1
countC(1) = L+1
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(12)  = NF_PUT_VARA_DOUBLE(ncid, varid_horizcors_flux, startC, countC, &
                               horizcors_flux(0:L, 1:nmonth, 1:nfields))

startC(1) = 1
countC(1) = L+1
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(13)  = NF_PUT_VARA_DOUBLE(ncid, varid_var_hspec_flux, startC, countC, var_hspec_flux(0:L, 1:nmonth, 1:nfields))

startC(1) = 1
countC(1) = L+1
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(14)  = NF_PUT_VARA_DOUBLE(ncid, varid_std_hspec_flux, startC, countC, std_hspec_flux(0:L, 1:nmonth, 1:nfields))

startA(1)  = 1
countA(1)  = 1
IF (ForceVCorr) THEN
  options(1) = 1
ELSE
  options(1) = 0
END IF
ierr(15)   = NF_PUT_VARA_INT(ncid, varid_ForceVCorr, startA, countA, options(1:1))

startA(1)  = 1
countA(1)  = 1
options(1) = vert_covs
ierr(16)   = NF_PUT_VARA_INT(ncid, varid_vert_covs, startA, countA, options(1:1))

startB(1) = 1
countB(1) = ylev
startB(2) = 1
countB(2) = ylev
ierr(17)  = NF_PUT_VARA_DOUBLE(ncid, varid_glob_av_vert_cov, startB, countB, &
                               glob_av_vert_cov(1:ylev, 1:ylev))

startB(1) = 1
countB(1) = ylev
startB(2) = 1
countB(2) = ylev
ierr(18)  = NF_PUT_VARA_DOUBLE(ncid, varid_vert_eigenvec_tracer, startB, countB, vert_eigenvec_tracer(1:ylev, 1:ylev))

startC(1) = 1
countC(1) = ylev
startC(2) = 1
countC(2) = ylat
startC(3) = 1
countC(3) = ylon
ierr(19)  = NF_PUT_VARA_DOUBLE(ncid, varid_var_vspec_tracer, startC, countC, var_vspec_tracer(1:ylev, 1:ylat, 1:ylon))

startC(1) = 1
countC(1) = ylev
startC(2) = 1
countC(2) = ylat
startC(3) = 1
countC(3) = ylon
ierr(20)  = NF_PUT_VARA_DOUBLE(ncid, varid_std_vspec_tracer, startC, countC, std_vspec_tracer(1:ylev, 1:ylat, 1:ylon))

startC(1) = 1
countC(1) = ylev
startC(2) = 1
countC(2) = ylat
startC(3) = 1
countC(3) = ylon
ierr(21)  = NF_PUT_VARA_DOUBLE(ncid, varid_vert_adjust_tracer, startC, countC, vert_adjust_tracer(1:ylev, 1:ylat, 1:ylon))

startA(1)  = 1
countA(1)  = nfields
ierr(22)   = NF_PUT_VARA_INT(ncid, varid_temporal_covs, startA, countA, temporal_covs(1:nfields))

startA(1) = 1
countA(1) = nfields
ierr(23)  = NF_PUT_VARA_INT(ncid, varid_temporal_cor_tpe, startA, countA, temporal_cor_tpe(1:nfields))

startA(1) = 1
countA(1) = nfields
ierr(24)  = NF_PUT_VARA_DOUBLE(ncid, varid_timescale_flux, startA, countA, timescale_flux(1:nfields))

startB(1) = 1
countB(1) = nmonth
startB(2) = 1
countB(2) = nfields
ierr(25)  = NF_PUT_VARA_DOUBLE(ncid, varid_temporalcors_flux, startB, countB, temporalcors_flux(1:nmonth,1:nfields))

startC(1) = 1
countC(1) = nmonth
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(26)  = NF_PUT_VARA_DOUBLE(ncid, varid_temp_cor_matrix, startC, countC, &
                               temp_cor_matrix(1:nmonth, 1:nmonth, 1:nfields))

startC(1) = 1
countC(1) = nmonth
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(27)  = NF_PUT_VARA_DOUBLE(ncid, varid_temp_eigenvec_flux, startC, countC, &
                               temp_eigenvec_flux(1:nmonth, 1:nmonth, 1:nfields))

startB(1) = 1
countB(1) = nmonth
startB(2) = 1
countB(2) = nfields
ierr(28)  = NF_PUT_VARA_DOUBLE(ncid, varid_var_tspec_flux, startB, countB, var_tspec_flux(1:nmonth, 1:nfields))

startB(1) = 1
countB(1) = nmonth
startB(2) = 1
countB(2) = nfields
ierr(29)  = NF_PUT_VARA_DOUBLE(ncid, varid_std_tspec_flux, startB, countB, std_tspec_flux(1:nmonth, 1:nfields))

DO i = 1, 29
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error outputting main variables ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO

!Close-up the file
!-----------------
ierr(1) = NF_CLOSE(ncid)

DO i = 1, 1
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error closing netcdf file ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


END SUBROUTINE CVT_matrices_output




!=================================================================================================
SUBROUTINE cvt_matrices_input (filename)

!********************************************************
!* Subroutine to write the CVT matrices                 *
!*                                                      *
!*                                                      *
!* R. Bannister, 03-12-19                               *
!* Tweaks, R. Bannister, July 2021                      *
!*                                                      *
!* Note that it is assumed that the arrays are already  *
!* allocated                                            *
!********************************************************

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields

USE onedvar_data, ONLY: L,                         &
                        CVT_std,                   &
                        CVT_HorizCors,             &
                        CVT_VertCors,              &
                        CVT_TemporalCors

IMPLICIT NONE

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE 'netcdf.inc'


! Declare subroutine parameter
! ----------------------------
CHARACTER(LEN=*),    INTENT(IN)  :: filename

!Declare local variables
!------------------------
INTEGER                       :: ncid
INTEGER                       :: dimidx_tomcat, dimidy_tomcat, dimidz_tomcat, dimidx_sh, dimidy_sh
INTEGER                       :: dimid_vertmode, dimid_month, dimid_timescale, dimid_field, dimid_totwn
INTEGER                       :: dimid_options
INTEGER                       :: varidx_tomcat, varidy_tomcat, varidz_tomcat, varidx_sh, varidy_sh
INTEGER                       :: varid_vertmode, varid_month, varid_timescale, varid_field, varid_totwn
INTEGER                       :: varid_FracOfPrior_tracer, varid_FracOfPrior_flux
INTEGER                       :: varid_std_tracer, varid_vert_adjust_tracer, varid_vert_eigenvec_tracer
INTEGER                       :: varid_var_vspec_tracer, varid_std_vspec_tracer
INTEGER                       :: varid_var_hspec_tracer, varid_std_hspec_tracer, varid_std_flux, varid_timescale_flux
INTEGER                       :: varid_temporal_cor_tpe, varid_temp_eigenvec_flux, varid_horizcors_flux
INTEGER                       :: varid_horizflux_cor_tpe, varid_lengthscale_flux, varid_temporalcors_flux
INTEGER                       :: varid_var_tspec_flux, varid_std_tspec_flux, varid_var_hspec_flux, varid_std_hspec_flux
INTEGER                       :: varid_temporal_covs, varid_vert_covs, varid_options, varid_ForceVCorr
INTEGER                       :: options(1:1)
INTEGER                       :: startA(1), countA(1), startB(2), countB(2)
INTEGER                       :: startC(3), countC(3), startD(4), countD(4)
INTEGER                       :: ierr(1:24), i
REAL                          :: temp(1:1)


!*****************************************************************************************
PRINT*, 'Reading CVT data'
!*****************************************************************************************

! Open netCDF file
!-------------------------------------
ierr = NF_OPEN(filename, NF_NOWRITE, ncid)
DO i = 1, 1
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error opening netcdf file ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


!Get the necessary dimension ids
!-------------------------------
ierr(1) = NF_INQ_DIMID(ncid, 'x_tomcat',  dimidx_tomcat)
ierr(2) = NF_INQ_DIMID(ncid, 'y_tomcat',  dimidy_tomcat)
ierr(3) = NF_INQ_DIMID(ncid, 'z_tomcat',  dimidz_tomcat)
ierr(4) = NF_INQ_DIMID(ncid, 'x_sh',      dimidx_sh)
ierr(5) = NF_INQ_DIMID(ncid, 'y_sh',      dimidy_sh)
ierr(6) = NF_INQ_DIMID(ncid, 'vertmode',  dimid_vertmode)
ierr(7) = NF_INQ_DIMID(ncid, 'month',     dimid_month)
ierr(8) = NF_INQ_DIMID(ncid, 'timescale', dimid_timescale)
ierr(9) = NF_INQ_DIMID(ncid, 'fluxfield', dimid_field)
ierr(10)= NF_INQ_DIMID(ncid, 'tot_wn',    dimid_totwn)
ierr(11)= NF_INQ_DIMID(ncid, 'options',   dimid_options)

DO i = 1, 11
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error obtaining dimension ids ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


!Get the variable ids for the dimension variables
!------------------------------------------------
ierr(1)  = NF_INQ_VARID(ncid, 'x_tomcat',  varidx_tomcat)
ierr(2)  = NF_INQ_VARID(ncid, 'y_tomcat',  varidy_tomcat)
ierr(3)  = NF_INQ_VARID(ncid, 'z_tomcat',  varidz_tomcat)
ierr(4)  = NF_INQ_VARID(ncid, 'x_sh',      varidx_sh)
ierr(5)  = NF_INQ_VARID(ncid, 'y_sh',      varidy_sh)
ierr(6)  = NF_INQ_VARID(ncid, 'vertmode',  varid_vertmode)
ierr(7)  = NF_INQ_VARID(ncid, 'month',     varid_month)
ierr(8)  = NF_INQ_VARID(ncid, 'timescale', varid_timescale)
ierr(9)  = NF_INQ_VARID(ncid, 'fluxfield', varid_field)
ierr(10) = NF_INQ_VARID(ncid, 'totwn',     varid_totwn)
ierr(11) = NF_INQ_VARID(ncid, 'options',   varid_options)

DO i = 1, 11
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error obtaining dimension variable ids ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


!Get the variable ids for the main variables
!-------------------------------------------
ierr(1:24) = 0
ierr(1)  = NF_INQ_VARID(ncid, 'temporal_covs',        varid_temporal_covs)
ierr(2)  = NF_INQ_VARID(ncid, 'ForceVCorr',           varid_ForceVCorr)
ierr(3)  = NF_INQ_VARID(ncid, 'vert_covs',            varid_vert_covs)
ierr(4)  = NF_INQ_VARID(ncid, 'FracOfPrior_tracer',   varid_FracOfPrior_tracer)
ierr(5)  = NF_INQ_VARID(ncid, 'std_tracer',           varid_std_tracer)
ierr(6)  = NF_INQ_VARID(ncid, 'vert_adjust_tracer',   varid_vert_adjust_tracer)
ierr(7)  = NF_INQ_VARID(ncid, 'vert_eigenvec_tracer', varid_vert_eigenvec_tracer)
ierr(8)  = NF_INQ_VARID(ncid, 'var_vspec_tracer',     varid_var_vspec_tracer)
ierr(9)  = NF_INQ_VARID(ncid, 'std_vspec_tracer',     varid_std_vspec_tracer)
ierr(10) = NF_INQ_VARID(ncid, 'var_hspec_tracer',     varid_var_hspec_tracer)
ierr(11) = NF_INQ_VARID(ncid, 'std_hspec_tracer',     varid_std_hspec_tracer)
ierr(12) = NF_INQ_VARID(ncid, 'FracOfPrior_flux',     varid_FracOfPrior_flux)
ierr(13) = NF_INQ_VARID(ncid, 'std_flux',             varid_std_flux)
!ierr(14) = NF_INQ_VARID(ncid, 'timescale_flux',       varid_timescale_flux)
!ierr(15) = NF_INQ_VARID(ncid, 'temporalcors_flux',    varid_temporalcors_flux)
!ierr(16) = NF_INQ_VARID(ncid, 'temporal_cor_tpe',     varid_temporal_cor_tpe)
ierr(17) = NF_INQ_VARID(ncid, 'temp_eigenvec_flux',   varid_temp_eigenvec_flux)
!ierr(18) = NF_INQ_VARID(ncid, 'horizcors_flux',       varid_horizcors_flux)
!ierr(19) = NF_INQ_VARID(ncid, 'horizflux_cor_tpe',    varid_horizflux_cor_tpe)
!ierr(20) = NF_INQ_VARID(ncid, 'lengthscale_flux',     varid_lengthscale_flux)
ierr(21) = NF_INQ_VARID(ncid, 'var_tspec_flux',       varid_var_tspec_flux)
ierr(22) = NF_INQ_VARID(ncid, 'std_tspec_flux',       varid_std_tspec_flux)
ierr(23) = NF_INQ_VARID(ncid, 'var_hspec_flux',       varid_var_hspec_flux)
ierr(24) = NF_INQ_VARID(ncid, 'std_hspec_flux',       varid_std_hspec_flux)

DO i = 1, 24
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error obtaining main variable ids ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


!---------------------------------------------
! Get the required dimension variables
! --------------------------------------------

startA(1) = 1
countA(1) = ylon
ierr(1)   = NF_GET_VARA_DOUBLE(ncid, varidx_tomcat, startA, countA, CVT_HorizCors % TOMCAT_longs(1:ylon))
countA(1) = ylat
ierr(2)   = NF_GET_VARA_DOUBLE(ncid, varidy_tomcat, startA, countA, CVT_HorizCors % TOMCAT_lats(1:ylat))

countA(1) = 2*L+1
ierr(3)   = NF_GET_VARA_DOUBLE(ncid, varidx_sh, startA, countA, CVT_HorizCors % SHtools_longs(1:2*L+1))
countA(1) = L+1
ierr(4)   = NF_GET_VARA_DOUBLE(ncid, varidy_sh, startA, countA, CVT_HorizCors % SHtools_lats(1:L+1))

countA(1) = ylev
ierr(5)   = NF_GET_VARA_DOUBLE(ncid, varidz_tomcat, startA, countA, CVT_VertCors % alts(1:ylev))
countA(1) = nmonth
ierr(6)   = NF_GET_VARA_DOUBLE(ncid, varid_month, startA, countA, CVT_TemporalCors % months(1:nmonth))

DO i = 1, 6
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error obtaining dimension variable ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO



!--------------------------------------------
! Get the values of the main variables
! -------------------------------------------
ierr(1:24) = 0

startA(1)  = 1
countA(1)  = nfields
ierr(1)    = NF_GET_VARA_INT(ncid, varid_temporal_covs, startA, countA, &
                             CVT_TemporalCors % Temporal_covs(1:nfields))

countA(1)  = 1
ierr(2)    = NF_GET_VARA_INT(ncid, varid_ForceVCorr, startA, countA, options(1:1))
CVT_VertCors % ForceVCorr = (options(1) == 1)

startA(1)  = 1
ierr(3)    = NF_GET_VARA_INT(ncid, varid_vert_covs, startA, countA, options(1:1))
CVT_VertCors % vert_covs = options(1)

startA(1)  = 1
ierr(4)    = NF_GET_VARA_DOUBLE(ncid, varid_FracOfPrior_tracer, startA, countA, temp)
CVT_std % FracOfPrior_tracer = temp(1)

startC(1) = 1
countC(1) = ylat
startC(2) = 1
countC(2) = ylon
startC(3) = 1
countC(3) = ylev
ierr(5)   = NF_GET_VARA_DOUBLE(ncid, varid_std_tracer, startC, countC, &
                               CVT_std % std_tracer(1:ylat, 1:ylon, 1:ylev))

startC(1) = 1
countC(1) = ylev
startC(2) = 1
countC(2) = ylat
startC(3) = 1
countC(3) = ylon
ierr(6)   = NF_GET_VARA_DOUBLE(ncid, varid_vert_adjust_tracer, startC, countC, &
                               CVT_VertCors % vert_adjust_tracer(1:ylev, 1:ylat, 1:ylon))

startB(1) = 1
countB(1) = ylev
startB(2) = 1
countB(2) = ylev
ierr(7)   = NF_GET_VARA_DOUBLE(ncid, varid_vert_eigenvec_tracer, startB, countB, &
                               CVT_VertCors % vert_eigenvec_tracer(1:ylev, 1:ylev))

startC(1) = 1
countC(1) = ylev
startC(2) = 1
countC(2) = ylat
startC(3) = 1
countC(3) = ylon
ierr(8)   = NF_GET_VARA_DOUBLE(ncid, varid_var_vspec_tracer, startC, countC, &
                               CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1:ylon))

startC(1) = 1
countC(1) = ylev
startC(2) = 1
countC(2) = ylat
startC(3) = 1
countC(3) = ylon
ierr(9)   = NF_GET_VARA_DOUBLE(ncid, varid_std_vspec_tracer, startC, countC, &
                               CVT_VertCors % std_vspec_tracer(1:ylev, 1:ylat, 1:ylon))

startB(1) = 1
countB(1) = L+1
startB(2) = 1
countB(2) = ylev
ierr(10)  = NF_GET_VARA_DOUBLE(ncid, varid_var_hspec_tracer, startB, countB, &
                               CVT_HorizCors % var_hspec_tracer(0:L, 1:ylev))

startB(1) = 1
countB(1) = L+1
startB(2) = 1
countB(2) = ylev
ierr(11)  = NF_GET_VARA_DOUBLE(ncid, varid_std_hspec_tracer, startB, countB, &
                               CVT_HorizCors % std_hspec_tracer(0:L, 1:ylev))

startA(1)  = 1
ierr(12)   = NF_GET_VARA_DOUBLE(ncid, varid_FracOfPrior_flux, startA, countA, temp)
CVT_std % FracOfPrior_flux = temp(1)

startD(1) = 1
countD(1) = ylat
startD(2) = 1
countD(2) = ylon
startD(3) = 1
countD(3) = nmonth
startD(4) = 1
countD(4) = nfields
ierr(13)  = NF_GET_VARA_DOUBLE(ncid, varid_std_flux, startD, countD, &
                               CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))

!startA(1) = 1
!countA(1) = nfields
!ierr(14)  = NF_GET_VARA_DOUBLE(ncid, varid_timescale_flux, startA, countA, &
!                               CVT_TemporalCors % timescale_flux(1:nfields))

!startB(1) = 1
!countB(1) = nmonth
!startB(2) = 1
!countB(2) = nfields
!ierr(15)  = NF_GET_VARA_DOUBLE(ncid, varid_temporalcors_flux, startB, countB, &
!                               CVT_TemporalCors % temporalcors_flux(1:nmonth,1:nfields))

!startA(1) = 1
!countA(1) = nfields
!ierr(16)  = NF_GET_VARA_INT(ncid, varid_temporal_cor_tpe, startA, countA, &
!                            CVT_TemporalCors % temporal_cor_tpe(1:nfields))

startC(1) = 1
countC(1) = nmonth
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(17)  = NF_GET_VARA_DOUBLE(ncid, varid_temp_eigenvec_flux, startC, countC, &
                               CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, 1:nmonth, 1:nfields))

!startC(1) = 1
!countC(1) = L+1
!startC(2) = 1
!countC(2) = nmonth
!startC(3) = 1
!countC(3) = nfields
!ierr(18)  = NF_GET_VARA_DOUBLE(ncid, varid_horizcors_flux, startC, countC, &
!                               CVT_HorizCors % horizcors_flux(1:L+1, 1:nmonth, 1:nfields))

!startA(1) = 1
!countA(1) = nfields
!ierr(19)  = NF_GET_VARA_INT(ncid, varid_horizflux_cor_tpe, startA, countA, &
!                            CVT_HorizCors % horizflux_cor_tpe(1:nfields))

!startB(1) = 1
!countB(1) = nmonth
!startB(2) = 1
!countB(2) = nfields
!ierr(20)  = NF_GET_VARA_DOUBLE(ncid, varid_lengthscale_flux, startB, countB, &
!                               CVT_HorizCors % lengthscale_flux(1:nmonth, 1:nfields))

startB(1) = 1
countB(1) = nmonth
startB(2) = 1
countB(2) = nfields
ierr(21)  = NF_GET_VARA_DOUBLE(ncid, varid_var_tspec_flux, startB, countB, &
                               CVT_TemporalCors % var_tspec_flux(1:nmonth, 1:nfields))

startB(1) = 1
countB(1) = nmonth
startB(2) = 1
countB(2) = nfields
ierr(22)  = NF_GET_VARA_DOUBLE(ncid, varid_std_tspec_flux, startB, countB, &
                               CVT_TemporalCors % std_tspec_flux(1:nmonth, 1:nfields))

startC(1) = 1
countC(1) = L+1
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(23)  = NF_GET_VARA_DOUBLE(ncid, varid_var_hspec_flux, startC, countC, &
                               CVT_HorizCors % var_hspec_flux(0:L, 1:nmonth, 1:nfields))

startC(1) = 1
countC(1) = L+1
startC(2) = 1
countC(2) = nmonth
startC(3) = 1
countC(3) = nfields
ierr(24)  = NF_GET_VARA_DOUBLE(ncid, varid_std_hspec_flux, startC, countC, &
                               CVT_HorizCors % std_hspec_flux(0:L, 1:nmonth, 1:nfields))

DO i = 1, 24
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error obtaining main variable values ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO


!Close-up the file
!-----------------
ierr(1) = NF_CLOSE(ncid)

DO i = 1, 1
  IF (ierr(i) /= 0) THEN
    PRINT*, '***Error closing netcdf file ***'
    PRINT*,'ierr',  i, ierr(i),  NF_STRERROR(ierr(i))
    STOP
  END IF
END DO

END SUBROUTINE cvt_matrices_input

