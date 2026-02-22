/* ==================================================================================
   3d source and sink code
   Inverse tests

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make InverseTests_CVT.out

   Modification history
   --------------------
   28/09/25 Code adapted from Master_AdjointTests_CVT. Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <stdio.h>
//   #include <iostream.h>
   #include <math.h>
   #include <stdlib.h>
   #include <string.h>
   #include <source.h>


int main ()
{
  struct HorizTransData_type HorizData;
  struct state_type          state_real_in, state_real_out, state_spec, state_modes;
  struct CVTData_type        CVTdata;
  struct metadata_type       MetaData;
  int                        lon, lat, lev, t;
  double                     freq1, freq2;
  char                       CVTfilename[256]    = "../data/CVT_calib_56levs.nc";
  char                       IntputFile[256]     = "InverseTest_in.nc";
  char                       OutputFile_HT[256]  = "InverseTest_HorizTrans_out.nc";
  char                       OutputFile_VTT[256] = "InverseTest_VertTempTrans_out.nc";


//   ===============================================================================
//   Initialisation
//   ===============================================================================

  // Initialise and read in generic data needed for the horizontal transform
  printf ("Initialising horizontal transform\n");
  InitializeHorizTrans (&HorizData,
                        ::ALPfilename);

  // Read-in cvt data (including allocation of CVT)
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


//   ===============================================================================
//   Set-up arrays
//   ===============================================================================

  Allocate_state (&state_real_in,
                  &MetaData,
                  'r', 'r', 'r');

  Allocate_state (&state_real_out,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&state_spec,
                  &MetaData,
                  's', 'r', 'r');

  Allocate_state (&state_modes,
                  &MetaData,
                  'r', 'm', 'm');

  // Put some numbers into the real space tracer array
  for (lon=1; lon<=MetaData.nlon; lon++)
  { for (lat=1; lat<=MetaData.nlat; lat++)
    { for (lev=1; lev<=MetaData.nlev; lev++)
      { freq1 = ::pi * double(lev) / (10.0 * double(MetaData.nlon));
        freq2 = ::pi * double(lev) / (1.8 * double(MetaData.nlat));
        state_real_in.tracer0_rs[lon][lat][lev] = sin(double(lon) * freq1) +
                                                  sin(double(lat) * freq2);
      }
    }
  }

  // Put some numbers into the real space flux array
  for (lon=1; lon<=MetaData.nlon; lon++)
  { for (lat=1; lat<=MetaData.nlat; lat++)
    { for (t=0; t<MetaData.nss; t++)
      { freq1 = ::pi * double(t+1) / (4.0 * double(MetaData.nlon));
        freq2 = ::pi * double(t+1) / (0.2 * double(MetaData.nlat));
        state_real_in.source_rs[lon][lat][t] = sin(double(lon) * freq1) +
                                               sin(double(lat) * freq2);
      }
    }
  }

  halos (&state_real_in);

  // Output the real space state
  WriteStateVector ( &state_real_in,
                     IntputFile );




//   ===============================================================================
//   Inverse test for horizontal transform
//   ===============================================================================

  printf ("Inverse test of horizontal transform\n");

  // Do the inverse horizontal transform
  cvt_h_inv ( &state_spec,     // out
              &state_real_in,  // in
              &HorizData,
              &CVTdata,
              'p' );           // pseudo inverse

  // Do the forward horizontal transform
  cvt_h ( &state_spec,         // in
          &state_real_out,     // out
          &HorizData,
          &CVTdata );

  // Output the recovered real space state
  WriteStateVector ( &state_real_out,
                     OutputFile_HT );


//   ===============================================================================
//   Inverse test for vertical and temporal transform
//   (these are independent and hence done together)
//   ===============================================================================

  printf ("Inverse test of vertical and temporal transforms\n");

  // Do the inverse vertical transform
  cvt_v_inv ( &state_modes,    // out
              &state_real_in,  // in
              &CVTdata );

  // Do the inverse temporal transform
  cvt_t_inv ( &state_real_in,  // in
              &state_modes,    // out
              &CVTdata );

  // Do the forward vertical transform
  cvt_v ( &state_real_out,     // out
          &state_modes,        // in
          &CVTdata );

  // Do the forward temporal transform
  cvt_t ( &state_real_out,     // out
          &state_modes,        // in
          &CVTdata );

  // Output the recovered real space state
  WriteStateVector ( &state_real_out,
                     OutputFile_VTT );


//   ===============================================================================
//   Deallocate
//   ===============================================================================

  DeallocateHorizTrans (&HorizData);
  Deallocate_CVT (&CVTdata);
  Deallocate_metadata (&MetaData);

  Deallocate_state (&state_real_in);
  Deallocate_state (&state_real_out);
  Deallocate_state (&state_spec);
  Deallocate_state (&state_modes);


}
