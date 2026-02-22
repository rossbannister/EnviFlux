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
  const char                 CVTfilename[196] = "../Calibration/CVT_calib.nc";


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
                  's', 'r', 'm');


  // Initialise input state
  for (z=0; z<CVTdata.nlev; z++)
  { for (l=0; l<=CVTdata.L; l++)
    { state_u.tracer0_ss[0][l][0][z] = 0.0;
      if (l > 0)
      { // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { input.tracer0_ss[cs][l][m][z] = 0.0;
          }
        }
      }
    }
  }
  for (t=0; t<(CVTdata.nss); t++)
  { for (l=0; l<=CVTdata.L; l++)
    { state_u.source_ss[0][l][0][t] = 0.0;
      if (l > 0)
      { // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { input.source_ss[cs][l][m][t] = 0.0;
          }
        }
      }
    }
  }

  // Forward operator
  cvt_h (&state_u,     // in  spectral space (in horiz)
         &state_v,     // out real space (in horiz)
         &HorizData,   // in  info about horiz transform
         &CVTdata);    // in  info about cvt

  // Adjoint operator
  cvt_h_adj (&state_w,     // out spectral space (in horiz)
             &state_v,     // in  real space (in horiz)
             &HorizData,   // in  info about horiz transform
             &CVTdata);    // in  info about cvt


  lhs = innerproduct_general (&state_v,      // in   input state 1
                              &state_v,      // in   input state 2
                              &HorizData);   // in   contains metadata
  rhs = innerproduct_general (&state_u,      // in   input state 1
                              &state_w,      // in   input state 2
                              &HorizData);   // in   contains metadata
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  // Tidy up
  Deallocate_state (&state_u);
  Deallocate_state (&state_v);
  Deallocate_state (&state_w);



//   ===============================================================================
//   Adjoint test for vertical and temporal transform
//   ===============================================================================

  printf ("Combined adjoint test of vertical and temporal transforms\n");

  Allocate_state (&state_u,
                  &MetaData,
                  'r', 'r', 'm');
  Allocate_state (&state_v,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&state_w,
                  &MetaData,
                  'r', 'r', 'm');

  // Put some random data in state_u
  // Populate the initial tracer
  for (z=0; z<CVTdata.nss; z++)
  { for (x=0; x<CVTdata.nlon; x++)
    { for (y=0; y<CVTdata.nlat; y++)
      { state_u.tracer0_rs[x][y][z] = randomno(&rndseq);
      }
    }
  }

  // Populate the source
  for (t=0; t<CVTdata.nss; t++)
  { for (x=0; x<CVTdata.nlon; x++)
    { for (y=0; y<CVTdata.nlat; y++)
      { state_u.source_rs[x][y][t] = randomno(&rndseq);
      }
    }
  }

  // Forward operators (these act on different parts of the states)
  cvt_v (&state_v,     // out spectral space (in horiz)
         &state_u,     // in  real space (in horiz)
         &CVTdata);    // in  info about cvt
  cvt_t (&state_v,     // out spectral space (in horiz)
         &state_u,     // in  real space (in horiz)
         &CVTdata);    // in  info about cvt

  // Adjoint operators (these act on different parts of the states)
  cvt_v (&state_w,     // out spectral space (in horiz)
         &state_v,     // in  real space (in horiz)
         &CVTdata);    // in  info about cvt
  cvt_t_adj (&state_v,     // in  spectral space (in horiz)
             &state_w,     // out real space (in horiz)
             &CVTdata);    // in  info about cvt

  lhs = innerproduct_general (&state_v,      // in   input state 1
                              &state_v,      // in   input state 2
                              &HorizData);   // in   contains metadata
  rhs = innerproduct_general (&state_u,      // in   input state 1
                              &state_w,      // in   input state 2
                              &HorizData);   // in   contains metadata
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  // Tidy up
  Deallocate_state (&state_u);
  Deallocate_state (&state_v);
  Deallocate_state (&state_w);


//   ===============================================================================
//   Adjoint test for total CVT transform
//   ===============================================================================

  printf ("Adjoint test of total transform\n");

  Allocate_state (&state_u,
                  &MetaData,
                  's', 'r', 'm');
  Allocate_state (&state_v,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&state_w,
                  &MetaData,
                  's', 'r', 'm');

  // Put some random data in state_u
  for (z=0; z<CVTdata.nlev; z++)
  { for (l=0; l<=CVTdata.L; l++)
    { state_u.tracer0_ss[0][l][0][z] = randomno(&rndseq);
      if (l > 0)
      { // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { state_u.tracer0_ss[cs][l][m][z] = randomno(&rndseq);
          }
        }
      }
    }
  }
  for (t=0; t<CVTdata.nss; t++)
  { for (l=0; l<=CVTdata.L; l++)
    { state_u.source_ss[0][l][0][t] = randomno(&rndseq);
      if (l > 0)
      { // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { state_u.source_ss[cs][l][m][t] = randomno(&rndseq);
          }
        }
      }
    }
  }

  // Forward operator
  cvt_total (&state_u,     // in
             &state_v,     // out
             &MetaData,
             &HorizData,
             &CVTdata);

  // Adjoint operator
  cvt_total_adj (&state_w,     // out
                 &state_v,     // in
                 &MetaData,
                 &HorizData,
                 &CVTdata);

  lhs = innerproduct_general (&state_v,      // in   input state 1
                              &state_v,      // in   input state 2
                              &HorizData);   // in   contains metadata
  rhs = innerproduct_general (&state_u,      // in   input state 1
                              &state_w,      // in   input state 2
                              &HorizData);   // in   contains metadata
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  // Tidy up
  Deallocate_state (&state_u);
  Deallocate_state (&state_v);
  Deallocate_state (&state_w);

//   ===============================================================================
//   Deallocate
//   ===============================================================================

  DeallocateHorizTrans (&HorizData);
  Deallocate_CVT (&CVTdata);
  Deallocate_metadata (&MetaData);

}
