/* ==================================================================================
   3d source and sink code
   Test to read and write data to check input/output routines

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_ImpliedCov.out

   Modification history
   --------------------
   20/12/21 New Code. Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <source.h>


int main ()
{
  struct HorizTransData_type HorizData;
  struct state_type          field;
  struct Wind_type           wind;
  struct CVTData_type        CVTdata;
  struct metadata_type       MetaData;
  int                        x, y, z, t;
  char                       CVTfilename[256] = "../Calibration/CVT_calib.nc";


//   ===============================================================================
//   Initialisation
//   ===============================================================================

  // Initialise and read in generic data needed for the horizontal transform
  printf ("Initialising horizontal transform\n");
  InitializeHorizTrans (&HorizData,
                        ::ALPfilename);

  // Allocate the cvt
  Allocate_CVT (&CVTdata,
                ::L, ::nlev, ::ntimes_major);

  // Read-in cvt data
  printf ("Reading CVT data\n");
  cvt_matrices_input (&CVTdata,
                      CVTfilename);
  printf ("Read-in completed\n");

  // Set-up the meta data structure
  Allocate_metadata (&MetaData,
                     2*::L+1, ::L+1, ::nlev, ::ntimes_major);

  copy_metadata (CVTdata.nlon,         &(MetaData.nlon),
                 CVTdata.nlat,         &(MetaData.nlat),
                 CVTdata.nlev,         &(MetaData.nlev),
                 CVTdata.nss,          &(MetaData.nss),
                 CVTdata.L,            &(MetaData.L),
                 CVTdata.times,        MetaData.times,
                 CVTdata.longitude,    MetaData.longitude,
                 CVTdata.latitude,     MetaData.latitude,
                 CVTdata.level,        MetaData.level,
                 MetaData.cos_u_lat,   MetaData.cos_v_lat,
                 false);

  // Allocate space for data
  Allocate_state (&field,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_wind (&wind,
                 &MetaData,
                 false);

  // Read-in example state
  ReadStateVector (&field, &MetaData, true, true, "output.nc");

  // Read-in example wind
  Read_winds (&wind, &MetaData, 1.0, "../data/ECMWF_RedResWinds_1995/Winds0009.nc");

  // Write-out copy of state
  WriteStateVector (&field, "output_copy.nc");

  // Write-out copy of winds
  Write_winds (&wind, "Winds0009.nc");


//   ===============================================================================
//   Deallocate
//   ===============================================================================

  DeallocateHorizTrans (&HorizData);
  Deallocate_CVT (&CVTdata);
  Deallocate_metadata (&MetaData);
  Deallocate_state (&field);
  Deallocate_wind (&wind);
}
