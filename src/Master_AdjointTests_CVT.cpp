/* ==================================================================================
   3d source and sink code
   Adjoint tests

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make AdjointTests_CVT.out

   Modification history
   --------------------
   18/12/21 New Code. Ross Bannister

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
  struct state_type          state_u, state_v, state_w;
  double                     **RealSpace_v;
  double                     ***SpecSpace_u, ***SpecSpace_w;
  struct CVTData_type        CVTdata;
  struct metadata_type       MetaData;
  int                        x, y, z, t, cs, m, l;
  double                     lhs, rhs;
  int                        rndseq = 7384395;   // Random number seed
  char                       CVTfilename[256] = "../data/CVT_calib_56levs_37months.nc";


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
//   Adjoint test for halo routine
//   ===============================================================================

  printf ("Adjoint test of halo\n");
  printf ("Note: MetaData.nss   = %i\n", MetaData.nss);
  printf ("      ::ntimes_major = %i\n", ::ntimes_major);

  Allocate_state (&state_u,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&state_v,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&state_w,
                  &MetaData,
                  'r', 'r', 'r');

  // Put some random data in state_u (tracer and source)
  for (x=1; x<=MetaData.nlon; x++)
  { for (y=1; y<=MetaData.nlat; y++)
    { for (z=1; z<=MetaData.nlev; z++)
      { state_u.tracer0_rs[x][y][z] = 200.0 * (randomno(&rndseq) - 0.5);
      }
      for (t=0; t<MetaData.nss; t++)
      { state_u.source_rs[x][y][t] = 20.0 * (randomno(&rndseq) - 0.5);
      }
    }
  }

  // Copy to state_v
  copy_general (&state_u,
                &state_v,
                &HorizData,
                true);
  // Forward operator
  halos (&state_v);

  // Copy to state_2
  copy_general (&state_v,
                &state_w,
                &HorizData,
                true);
  // Adjoint operator
  halos_adj (&state_w);


  lhs = innerproduct_general (&state_v,      // in   input state 1
                              &state_v,      // in   input state 2
                              &HorizData,    // in   contains metadata
                              true);
  rhs = innerproduct_general (&state_u,      // in   input state 1
                              &state_w,      // in   input state 2
                              &HorizData,    // in   contains metadata
                              true);
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  // Tidy up
  Deallocate_state (&state_u);
  Deallocate_state (&state_v);
  Deallocate_state (&state_w);




//   ===============================================================================
//   Adjoint test for spherical transform
//   ===============================================================================

  printf ("Adjoint test of spherical operator\n");

  // Initialise data structures
  Array_2d_double_create (&RealSpace_v,
                          HorizData.nlon+1, HorizData.nlat+1);
  Array_3d_double_create (&SpecSpace_u,
                          2, HorizData.L+1, HorizData.L+1);
  Array_3d_double_create (&SpecSpace_w,
                          2, HorizData.L+1, HorizData.L+1);

  // Set the data in specspace_u
  for (cs=0; cs<2; cs++)
  { for (m=0; m<=HorizData.L; m++)
    { if ((m > 0) || (cs == 0))        // no sin term for mm=0
      { for (l=m; l<=HorizData.L; l++)
        { SpecSpace_u[cs][l][m] = randomno (&rndseq);
        }
      }
    }
  }

  // Compute v
  spherical_for ( RealSpace_v,   // output
                  SpecSpace_u,   // input
                  &HorizData );

  // Compute w
  spherical_adj ( RealSpace_v,   // input
                  SpecSpace_w,   // output
                  &HorizData );

  lhs = innerproduct_realsp (RealSpace_v,
                             RealSpace_v,
                             &HorizData,
                             false );

  rhs = innerproduct_spherical (SpecSpace_w,
                                SpecSpace_u,
                                &HorizData );

  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  // Tidy up
  Array_2d_double_destroy (&RealSpace_v,
                           HorizData.nlon+1);
  Array_3d_double_destroy (&SpecSpace_u,
                           2, HorizData.L+1);
  Array_3d_double_destroy (&SpecSpace_w,
                           2, HorizData.L+1);



//   ===============================================================================
//   Adjoint test for horizontal transform
//   ===============================================================================

  printf ("Adjoint test of horizontal transform\n");

  Allocate_state (&state_u,
                  &MetaData,
                  's', 'r', 'r');
  Allocate_state (&state_v,
                  &MetaData,
                  'r', 'r', 'r');
  Allocate_state (&state_w,
                  &MetaData,
                  's', 'r', 'r');

  // Put some random data in state_u
  for (z=1; z<=MetaData.nlev; z++)
  { for (l=0; l<=(::L); l++)
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
  for (t=0; t<(::ntimes_major); t++)
  { for (l=0; l<=(::L); l++)
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
                              &HorizData,    // in   contains metadata
                              false);
  rhs = innerproduct_general (&state_u,      // in   input state 1
                              &state_w,      // in   input state 2
                              &HorizData,    // in   contains metadata
                              false);
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
  for (z=0; z<(::ntimes_major); z++)
  { for (x=0; x<(2*::L+1); x++)
    { for (y=0; y<(::L+1); y++)
      { state_u.tracer0_rs[x][y][z] = randomno(&rndseq);
      }
    }
  }

  // Populate the source
  for (t=0; t<(::ntimes_major); t++)
  { for (x=0; x<(2*::L+1); x++)
    { for (y=0; y<(::L+1); y++)
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
                              &HorizData,    // in   contains metadata
                              false);
  rhs = innerproduct_general (&state_u,      // in   input state 1
                              &state_w,      // in   input state 2
                              &HorizData,    // in   contains metadata
                              false);
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
  for (z=1; z<=MetaData.nlev; z++)
  { for (l=0; l<=(::L); l++)
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
  for (t=0; t<(::ntimes_major); t++)
  { for (l=0; l<=(::L); l++)
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
                              &HorizData,    // in   contains metadata
                              true);
  rhs = innerproduct_general (&state_u,      // in   input state 1
                              &state_w,      // in   input state 2
                              &HorizData,    // in   contains metadata
                              true);
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
