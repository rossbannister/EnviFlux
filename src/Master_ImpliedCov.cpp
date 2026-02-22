/* ==================================================================================
   3d source and sink code
   Implied covariances

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make ImpliedCov.out

   Modification history
   --------------------
   20/12/21 New Code. Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <stdio.h>
   #include <math.h>
   #include <stdlib.h>
   #include <string.h>
   #include <source.h>


int main ()
{
  struct HorizTransData_type HorizData;
  struct state_type          input, output, interim;
  struct CVTData_type        CVTdata;
  struct metadata_type       MetaData;
  int                        x, y, z, t;
  char                       CVTfilename[256] = "/media/ross/banny/Enviflux/Bias/Covs/CVT__TracerFac_0.2__FluxFac_1.2.nc";
  char                       ImplCovFilename[256] = "/media/ross/banny/Enviflux/Bias/ImpliedCov/ImpliedCov.nc";


//   ===============================================================================
//   Initialisation
//   ===============================================================================

  // Initialise and read in generic data needed for the horizontal transform
  printf ("Initialising horizontal transform\n");
  InitializeHorizTrans (&HorizData,
                        ::ALPfilename);

  // Read-in and allocate cvt data
  cvt_matrices_input (&CVTdata,
                      CVTfilename);
  printf ("Read-in completed\n");

  // Set-up the meta data structure
  Allocate_metadata (&MetaData,
                     CVTdata.nlon,
                     CVTdata.nlat,
                     CVTdata.nlev,
                     CVTdata.nss);

  copy_metadata (CVTdata.nlon,       &(MetaData.nlon),
                 CVTdata.nlat,       &(MetaData.nlat),
                 CVTdata.nlev,       &(MetaData.nlev),
                 CVTdata.nss,        &(MetaData.nss),
                 CVTdata.L,          &(MetaData.L),
                 CVTdata.times,      MetaData.times,
                 CVTdata.longitude,  MetaData.longitude,
                 CVTdata.latitude,   MetaData.latitude,
                 CVTdata.level,      MetaData.level,
                 MetaData.cos_u_lat, MetaData.cos_v_lat,
                 false);



//   ===============================================================================
//   Do implied covariances
//   ===============================================================================

  printf ("Computing some implied covariances\n");

  Allocate_state (&input,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&output,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&interim,
                  &MetaData,
                  's', 'm', 'm');


  // The input state is automatically initialised to zero.
  // Put some unit spikes into the fields

  // Initial fields
  input.tracer0_rs[5][5][2]    = 1.0;
  input.tracer0_rs[20][15][20] = 1.0;
  input.tracer0_rs[30][25][40] = 1.0;
  // Flux
  input.source_rs[15][5][3]    = 1.0;
  input.source_rs[30][15][6]   = 1.0;
  input.source_rs[45][25][9]   = 1.0;

  // Adjoint operator
  cvt_total_adj (&interim,     // out
                 &input,       // in
                 &MetaData,
                 &HorizData,
                 &CVTdata);

  // Forward operator
  cvt_total (&interim,         // in
             &output,          // out
             &MetaData,
             &HorizData,
             &CVTdata);

  // Write the gradient vector as a check
  printf ("Writing out the implied cov field\n");
  WriteStateVector ( &output,
                     ImplCovFilename );

//   ===============================================================================
//   Deallocate
//   ===============================================================================

  DeallocateHorizTrans (&HorizData);
  Deallocate_CVT (&CVTdata);
  Deallocate_metadata (&MetaData);
  Deallocate_state (&input);
  Deallocate_state (&output);
  Deallocate_state (&interim);

}
