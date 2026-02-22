!+ Parameters used in RTTOV suite
!
MODULE cparam

  ! Description:
  ! Set up parameters that define the maximum dimension of some arrays
  !  used in the RTTOV suite.
  !
  ! Owner:
  ! SAF NWP
  !
  ! History:
  ! Version      Date        Comment
  ! 1            01/05/2000  Original code. .
  ! 2            21/08/2000  jpnssv increased from 1 to 6. Stephen English.
  !
  ! Code description:
  ! Language:              Fortran 90.
  ! Software Standards:    "European Standards for Writing and Documenting
  !                              Exchangeable Fortran 90 code".
  !
  IMPLICIT NONE

  ! Global parameters:
  INTEGER, PARAMETER :: jpnoaa = 12       ! Max no. of sats in coef file
  INTEGER, PARAMETER :: jptovs = 16       ! Max no. of sats in coef file
  INTEGER, PARAMETER :: jpdmsp =  9       ! Max no. of sats in coef file
  INTEGER, PARAMETER :: jpmet  =  4       ! Max no. of sats in coef file
  INTEGER, PARAMETER :: jpgoes =  5       ! Max no. of sats in coef file
  INTEGER, PARAMETER :: jpgms  =  5       ! Max no. of sats in coef file
  INTEGER, PARAMETER :: jpnsat = 12       ! Total max sats to be used
  INTEGER, PARAMETER :: jplev  = 43       ! No. of pressure levels
  INTEGER, PARAMETER :: jpnav  =  4       ! No. of profile variables
  INTEGER, PARAMETER :: jpnsav =  5       ! No. of surface air variables
  INTEGER, PARAMETER :: jpnssv =  6       ! No. of skin variables
  INTEGER, PARAMETER :: jpncv  =  2       ! No. of cloud variables
  INTEGER, PARAMETER :: jppf   = 100      ! Max no. profiles
  INTEGER, PARAMETER :: jpch   = 47       ! Max. no. of tovs channels
  INTEGER, PARAMETER :: jphir  = 20       ! Max. no. of hirs channels
  INTEGER, PARAMETER :: jpmsu  = 4        ! Max. no. of msu channels
  INTEGER, PARAMETER :: jpamsu = 20       ! Max. no. of amsu channels
  INTEGER, PARAMETER :: jpssu  =  3       ! Max. no. of ssu channels
  INTEGER, PARAMETER :: jpvtpr = 16       ! Max. no. of vtpr channels
  INTEGER, PARAMETER :: jpssmi = 7        ! Max. no. of ssm/i channels
  INTEGER, PARAMETER :: jptmi = 9         ! Max. no. of tmi channels
  INTEGER, PARAMETER :: jpavhrr = 3       ! Max. no. of avhrr channels
  INTEGER, PARAMETER :: jpgoesim = 4      ! Max. no. of goes imager channels
  INTEGER, PARAMETER :: jpgoessnd = 18    ! Max. no. of goes sounder channels
  INTEGER, PARAMETER :: jpchus = 47       ! Max. no. of channels used tovs
  INTEGER, PARAMETER :: jpchpf = jppf*jpchus ! Max no. of profs * chans used
  INTEGER, PARAMETER :: jpcofm = 10       ! Mixed gas coeffs (max
  INTEGER, PARAMETER :: jpcofw = 10       ! Water vapour coeffs (max
  INTEGER, PARAMETER :: jpcofo = 10       ! Ozone coeffs
  INTEGER, PARAMETER :: jpst   = 10       ! Max no. of surface types
  INTEGER, PARAMETER :: iu1    = 10       ! Unit for rt files
  INTEGER, PARAMETER :: nulout = 6        ! Unit for error messages
  INTEGER, PARAMETER :: jmwcldtop = 25    ! Upper level for lwp calcs
  !
  REAL, PARAMETER :: pi      = 3.1415926535
  REAL, PARAMETER :: deg2rad = pi/180.0   ! Degrees to radians conversion factor
  REAL, PARAMETER :: rad2deg = 180.0/pi   ! Radians to degrees conversion factor
  REAL, PARAMETER :: rcnv    = 6.03504E5  ! kg/kg--> ppmv
  REAL, PARAMETER :: gravity = 9.81       ! m/s^2
  REAL, PARAMETER :: Re      = 6371000.0  ! Radius of the Earth (m)
  REAL, PARAMETER :: RCH4    = 518.3      ! Specific gas constant for methane (J/kg/K)
  REAL, PARAMETER :: Rdry    = 287.05     ! Specific gas constant for dry air (J/kg/K)
  !
  ! Module arguments:

  ! Global scalars:F
  INTEGER :: njptovs     ! no. of sats in coef file
  INTEGER :: njpdmsp     ! no. of sats in coef file
  INTEGER :: njpmet      ! no. of sats in coef file
  INTEGER :: njpgoes     ! no. of sats in coef file
  INTEGER :: njpgms      ! no. of sats in coef file
  INTEGER :: njpnsat     ! Total max sats to be used
  INTEGER :: njplev      ! No. of pressure levels
  INTEGER :: njpnav      ! No. of profile variables
  INTEGER :: njpnsav     ! No. of surface air variables
  INTEGER :: njpnssv     ! No. of skin variables
  INTEGER :: njpncv      ! No. of cloud variables
  INTEGER :: njppf       ! Max no. profiles
  INTEGER :: njpch       ! Max. no. of tovs channels 
  INTEGER :: njphir      ! Max. no. of hirs channels
  INTEGER :: njpmsu      ! Max. no. of msu channels
  INTEGER :: njpssu      ! Max. no. of ssu channels
  INTEGER :: njpchus     ! Max. no. of channels used tovs
  INTEGER :: njpchpf     ! Max no. of profs * chans used
  INTEGER :: njpcofm     ! Mixed gas coeffs (max)
  INTEGER :: njpcofw     ! Water vapour coeffs (max)
  INTEGER :: njpcofo     ! Ozone coeffs (max)
  INTEGER :: njpst       ! Max no. of surface types
  INTEGER :: niu1        ! Unit for rt files 
  INTEGER :: nmwcldtop   ! Upper level for lwp calcs

  ! End of module arguments:

  !-----End of header----------------------------------------------------

END MODULE CPARAM
