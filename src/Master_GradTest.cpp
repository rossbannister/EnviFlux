/* ==================================================================================
   3d source and sink code
   Gradient test for PenAndGrad routine

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_GradTest.out

   Modification history
   --------------------
   03/01/24 Adapted from Master_AdjointTests_SL+ObsOp and Master_AdjointTests_CVT.
            Ross Bannister

   Documentation
   -------------

   =============================================================================== */

   #include <stdio.h>
   #include <math.h>
   #include <stdlib.h>
   #include <string.h>
   #include <source.h>


int main ()
{ char                       state_ic[256] = "../TestMakeBg/Bg.nc";
  char                       wind_dir[256] = "../data/ECMWF_RedResWinds_1995_56levs";
  char                       obsfile_in[256] = "../TestGenObs/SyntheticObs.dat";
  char                       obs_file_out[256] = "../GradTest/Obs.dat";
  char                       file_dps[256] = "../TestGenObs/DP_file.nc";
  char                       CVTfilename[256] = "../data/CVT_calib_56levs_37months.nc";
  char                       GradFilename[256] = "../GradTest/Gradient.nc";
  char                       GradTestFilename[256] = "../GradTest/GradientTest.dat";
  double                     Dt = 43200.0;
  double                     dt = 3600.0;
  double                     kappa_dt = 86.4;
  double                     kappa_h = 0.0;  //4000000.0;
  double                     kappa_v = 0.0;  //2.0;
  char                       interpolation_lc = 'c';
  double                     factor_w = 1.0;
  bool                       output_dps;
  bool                       inc_adv = true;
  bool                       inc_vert = true;
  struct metadata_type       MetaData;
  struct HorizTransData_type HorizData;
  struct state_type          state, grad, chi0, dchi1, chi1;
  struct state_type          *dummy;
  struct obs_trunk_type      obs_trunk;
  struct CVTData_type        CVTdata;
  double                     Jb0, Jo0, J0, Jb, Jo, J;
  int                        rndseq = 985395;   // Random number seed
  bool                       new_dp;            // Create new departure point file?
  int                        lev, l, m, t, cs;
  double                     ip_bare, alpha0, alpha;
  FILE*                      GradTestStats;
  int                        dp_file_err, ncid;

  // ===============================================================================


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

  // Read-in the state (containing the initial conditions and the surface fluxes)
  printf ("Reading-in the state vector (initial conditions and surface flux time sequence)\n");
  state.horiz_repres = 'r';
  state.vert_repres  = 'r';
  state.temp_repres  = 'r';
  ReadStateVector (&state,
                   &MetaData,
                   true, true,
                   state_ic);
  printf ("Reading-in done\n");

  printf ("Dimensions (nlon, nlat, nlev: %i, %i, %i)\n", MetaData.nlon, MetaData.nlat, MetaData.nlev);

  // Read-in observations
  printf ("Reading observations\n");
  ReadObservations (&MetaData,
                    &obs_trunk,
                    obsfile_in);
  printf ("Observations read-in\n");


  // Find out if the departure point file exists
  dp_file_err = nc_open (file_dps,
                         NC_NOWRITE,
                         &ncid);
  if (dp_file_err == 0)
  { // Yes
    printf ("The departure points file exists\n");
    new_dp      = false;  // Don't need to make a new dp file
    dp_file_err = nc_close (ncid);
  }
  else
  { // No
    printf ("The departure points file does not exist\n");
    new_dp      = true;
  }

//   ===============================================================================
//   Reference state calculations
//   ===============================================================================

  // Allocate space for the gradient
  printf ("Allocating space for the gradient\n");
  Allocate_state (&grad,
                  &MetaData,
                  's', 'm', 'm');

  // Allocate space for the reference state in spectral space
  printf ("Allocating space for the reference state\n");
  Allocate_state (&chi0,
                  &MetaData,
                  's', 'm', 'm');

  // Generate a random state to generate chi0
  printf ("Generating the reference state\n");
  PutWhiteNoiseControlVector (&chi0,
                              &rndseq);

  // Calculate penalty and gradient at chi0
  // The gradient has to be computed only once
  printf ("Computing the gradient\n");
  PenAndGrad (&state,           //in    background state
              &obs_trunk,       //inout observation structure
              &chi0,            //in    control variable
              false,            //in    chi0 is not zero
              &MetaData,        //in    metadata
              &HorizData,       //in    info about horiz transform
              &CVTdata,         //in    info about cvt
              &Jb0,             //out   background part of cost fn
              &Jo0,             //out   observation part of cost fn
              &J0,              //out   total cost fn
              &grad,            //out   gradient in control space
              true,             //in    switch to compute gradient
              Dt,               //in    major timestep
              dt,               //in    minor timestep
              kappa_dt,         //in    diffusion timestep size
              kappa_h,          //in    horiz diffusion coefficient
              kappa_v,          //in    vert diffusion coefficient
              inc_adv,          //in    include advection
              inc_vert,         //in    include vertical transport?
              interpolation_lc, //in    interpolation type
              factor_w,         //in    mult factor of vert winds
              wind_dir,         //in    directory containing driver winds
              file_dps,         //in    in or out file for departure points
              new_dp);          //in    true if dp file needed to be made

  new_dp = false; // No need to make a new departure points file
  printf ("Gradient computation done\n");
  printf ("  J, Jb, Jo = %f, %f, %f\n", J0, Jb0, Jo0);

  // Write the gradient vector as a check
  printf ("Writing out the gradient\n");
  WriteStateVector ( &grad,
                     GradFilename );

  printf ("Writing observations (inlcuding model obs)\n");
  WriteObservations ( &MetaData,
                      &obs_trunk,
                      obs_file_out );
  printf ("Done writing observations\n");



//   ===============================================================================
//   Set the perturbation direction, dchi1
//   ===============================================================================

  // Allocate space for the reference state in spectral space
  printf ("Allocating space for the perturbation direction\n");
  Allocate_state (&dchi1,
                  &MetaData,
                  's', 'm', 'm');

  printf ("Allocating space for the new perturbation\n");
  Allocate_state (&chi1,
                  &MetaData,
                  's', 'm', 'm');

  // Generate a random state to generate chi0
  printf ("Setting the perturbation direction\n");
  PutWhiteNoiseControlVector (&dchi1,
                              &rndseq);

  // Calculate the inner product between dchi1 and grad
  printf ("Computing the bare part of the inner product\n");
  ip_bare = innerproduct_general (&dchi1,
                                  &grad,
                                  &HorizData,
                                  false);     // Do not include halo level
  printf ("  Inner product = %f\n", ip_bare);

//   ===============================================================================
//   Scale the perturbation to different ammounts
//   ===============================================================================

  // Set-up output file
  GradTestStats = fopen (GradTestFilename, "w");

  // Loop over scales (log scale)
  for (alpha0=-10.0; alpha0<=9.0; alpha0++)
  { if (alpha0<1.0)
    { // Treat alpha0 as a log variable
      alpha = pow(10.0,alpha0);
    }
    else
    { // Treat alpha0 as a linear variable (add 1 to avoid doing alpha=1 twice)
      alpha = alpha0 + 1.0;
    }

    printf ("Value for alpha = %f\n", alpha);
    // Construct a new chi vector
    // chi1 = chi0 + alpha * dchi1

    add_general_fac (&chi0,
                     &dchi1,
                     &chi1,
                     &HorizData,
                     true,
                     alpha);

/*
    // Deal with the tracer part
    for (lev=1; lev<=MetaData.nlev; lev++)
    { for (l=0; l<=MetaData.L; l++)
      { // cosine for m=0
        chi1.tracer0_ss[0][l][0][lev] = chi0.tracer0_ss[0][l][0][lev] +
                                        dchi1.tracer0_ss[0][l][0][lev] * alpha;
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { chi1.tracer0_ss[cs][l][m][lev] = chi0.tracer0_ss[cs][l][m][lev] +
                                               dchi1.tracer0_ss[cs][l][m][lev] * alpha;
            }
          }
        }
      }
    }

    // Deal with the flux part
    for (t=0; t<MetaData.nss; t++)
    { for (l=0; l<=MetaData.L; l++)
      { // cosine for m=0
        chi1.source_ss[0][l][0][t] = chi0.source_ss[0][l][0][t] +
                                     dchi1.source_ss[0][l][0][t] * alpha;
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { chi1.source_ss[cs][l][m][t] = chi0.source_ss[cs][l][m][t] +
                                            dchi1.source_ss[cs][l][m][t] * alpha;
            }
          }
        }
      }
    }
*/

    // Compute the cost function for this
    printf ("  Computing the cost function for this alpha\n");
    PenAndGrad (&state,           //in    background state
                &obs_trunk,       //inout observation structure
                &chi1,            //in    control variable
                false,            //in    chi1 is not zero
                &MetaData,        //in    metadata
                &HorizData,       //in    info about horiz transform
                &CVTdata,         //in    info about cvt
                &Jb,              //out   background part of cost fn
                &Jo,              //out   observation part of cost fn
                &J,               //out   total cost fn
                dummy,            //out   gradient in control space
                false,            //in    switch to compute gradient
                Dt,               //in    major timestep
                dt,               //in    minor timestep
                kappa_dt,         //in    diffusion timestep size
                kappa_h,          //in    horiz diffusion coefficient
                kappa_v,          //in    vert diffusion coefficient
                inc_adv,          //in    include advection
                inc_vert,         //in    include vertical transport?
                interpolation_lc, //in    interpolation type
                factor_w,         //in    mult factor of vert winds
                wind_dir,         //in    directory containing driver winds
                file_dps,         //in    in or out file for departure points
                new_dp);          //in    true if dp file needed to be made
    printf ("  J, Jb, Jo = %f, %f, %f\n", J, Jb, Jo);

    fprintf (GradTestStats, "%f  %f\n", alpha,
                                      (J - J0) / (alpha * ip_bare));

  }

  fclose (GradTestStats);

  // Tidy up
  Deallocate_state (&state);
  Deallocate_state (&grad);
  Deallocate_state (&chi0);
  Deallocate_state (&dchi1);
  Deallocate_state (&chi1);
  Deallocate_metadata (&MetaData);
  delete[] obs_trunk.mass_profile;
  Destroy_Obs (&(obs_trunk.first));
  DeallocateHorizTrans (&HorizData);
  Deallocate_CVT (&CVTdata);
}
