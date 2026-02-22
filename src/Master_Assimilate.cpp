/* ==================================================================================
   EnviFlux 3d source and sink code
   Master routine for assimilation run to estimate initial conditions and surface fluxes

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_Assimilate.out

   Modification history
   --------------------
   11/03/24 Adapted from Master_GradTest.cpp
            Ross Bannister


   Documentation
   -------------
   =============================================================================== */

   #include <stdio.h>
   #include <math.h>
   #include <stdlib.h>
   #include <string.h>
   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )

{ char                       bg_filename[256];
  char                       wind_dir[256];
  char                       obsfile_in[256];
  char                       file_dps[256];
  char                       obsfile_out_bg[256];
  char                       CVTfilename[256];
  char                       anal_filename[256];
  char                       analinc_filename[256];
  char                       Assim_diags_filename[256];
  char                       obsfile_out_anal[256];
  double                     Dt, dt, kappa_dt, kappa_h, kappa_v;
  bool                       inc_vert;
  double                     factor_w;
  char                       interpolation_lc, descent_algorithm;
  double                     conv_criterion;
  int                        max_iters;
  struct HorizTransData_type HorizData;
  struct CVTData_type        CVTdata;
  double                     factor_std_tracer, factor_std_flux;
  struct metadata_type       MetaData;
  struct obs_trunk_type      obs_trunk;
  bool                       new_dp;      // Create new departure point file?
  int                        dp_file_err, ncid;
  int                        loop, lev, x, y, t;
  double                     Jb, Jo, J, Jb0, Jo0, J0, Jbp, Jop, Jp, Jbm, Jom, Jm;
  double                     alpha, beta, norm2_res_i, norm2_res_0, norm2_res_ip1;
  double                     norm2_chi, norm2_p_i, mu_modified, minus_mu_modified;
  struct state_type          bg, anal, analinc, pert;
  struct state_type          chiA, chiB, resA, resB, pA, pB;
  struct state_type          *chi_i, *chi_ip1, *res_i, *res_ip1, *p_i, *p_ip1;
  struct state_type          *swap, *dummy=NULL;
  double                     GradRatio;
  bool                       ok, converged;
  FILE*                      AssimDiags;

  // Process command line arguments
  ok = Assimilate_arguments (
         argument_count,      // in
         argument_list,       // in
         bg_filename,         // out Filename of background state
         wind_dir,            // out Directory containing driver winds
         obsfile_in,          // out Input file for observations
         file_dps,            // out Input or output file for departure points
         CVTfilename,         // out CVT filename (specification of B)
         &factor_std_tracer,  // out to multiply tracer err std
         &factor_std_flux,    // out to multiply flux err std
         obsfile_out_bg,      // out Output file of observation information at the background
         anal_filename,       // out Filename of analysis state
         analinc_filename,    // out Filename of analysis increment
         Assim_diags_filename,// out Filename of assimilation diagnostics
         obsfile_out_anal,    // out Output file of observation information at the analysis
         &Dt,                 // out Major timestep (between wind files, seconds)
         &dt,                 // out Minor timestep (for integration scheme, seconds)
         &kappa_dt,           // out Timestep for diffusion
         &kappa_h,            // out Horizontal diffusion coefficient
         &kappa_v,            // out Vertical diffusion coefficient
         &inc_vert,           // out Include vertical transport?
         &factor_w,           // out Multiplication factor of vertical winds
         &interpolation_lc,   // out Linear or cubic interpolation
         &descent_algorithm,  // out Descent algorithm type
         &conv_criterion,     // out Criterion for gradient norm for convergence
         &max_iters );        // out Maximum number of iterations

  if (ok)
  { // ===== Initialise and read in generic data needed for the horizontal transform
    printf ("Initialising horizontal transform\n");
    InitializeHorizTrans (&HorizData,
                          ::ALPfilename);


    // ===== Read-in cvt data (including allocation of CVT)
    printf ("Reading CVT data\n");
    cvt_matrices_input (&CVTdata,
                        CVTfilename);
    printf ("Reading CVT completed\n");

    // Adjust the standard deviations
    printf ("Adjusting the background error standard deviations\n");
    for (lev=0; lev<CVTdata.nlev; lev++)
    { CVTdata.tracer_stddev[lev] *= factor_std_tracer;
    }
    for (x=0; x<CVTdata.nlon; x++)
    { for (y=0; y<CVTdata.nlat; y++)
      { for (t=0; t<CVTdata.nss; t++)
        { CVTdata.source_stddev[x][y][t] *= factor_std_flux;
        }
      }
    }
    printf ("Adjustment completed\n");


    // ===== Set-up the meta data structure
    printf ("Allocating MetaData\n");
    Allocate_metadata (&MetaData,
                       2*::L+1, ::L+1, ::nlev, ::ntimes_major);
    printf ("Allocating MetaData completed\n");

    printf ("Copying MetaData\n");
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
                   true);
    printf ("Copying MetaData completed\n");


    // ===== Read in the background state
    printf ("Reading-in the background state\n");
    bg.horiz_repres = 'r';
    bg.vert_repres  = 'r';
    bg.temp_repres  = 'r';
    ReadStateVector (&bg,
                     &MetaData,
                     true, true,
                     bg_filename);
    printf ("Reading-in the background completed\n");


    // ===== Read-in observations
    printf ("Reading observations\n");
    ReadObservations (&MetaData,
                      &obs_trunk,
                      obsfile_in);
    printf ("Observations read-in\n");


    // ===== Find out if the departure point file exists
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


    // ===== Allocate space for analysis and analysis increment
    printf ("Allocating space for the analysis in real space\n");
    Allocate_state (&anal,
                    &MetaData,
                    'r', 'r', 'r');
    printf ("Allocating space for the analysis increment in real space\n");
    Allocate_state (&analinc,
                    &MetaData,
                    'r', 'r', 'r');


    // ===== Allocate space for control variables, residuals, search directions, and perturbation
    printf ("Allocating space for control variables, residuals, search directions, and perturbation\n");
    Allocate_state (&chiA,
                    &MetaData,
                    's', 'm', 'm');
    Allocate_state (&chiB,
                    &MetaData,
                    's', 'm', 'm');
    printf ("Allocating space for residuals\n");
    Allocate_state (&resA,
                    &MetaData,
                    's', 'm', 'm');
    Allocate_state (&resB,
                    &MetaData,
                    's', 'm', 'm');
    printf ("Allocating space for search directions\n");
    Allocate_state (&pA,
                    &MetaData,
                    's', 'm', 'm');
    printf ("Allocating space for the current search direction\n");
    Allocate_state (&pB,
                    &MetaData,
                    's', 'm', 'm');
    Allocate_state (&pert,
                    &MetaData,
                    's', 'm', 'm');
    printf ("Done allocations\n");

    // ===== Set-up initial pointers
    chi_i   = &chiA;  // Current control variable (initially zero)
    chi_ip1 = &chiB;  // Updated control variable (initially zero)
    res_i   = &resA;  // Current residual (initially zero)
    res_ip1 = &resB;  // Updated residual (initially zero)
    p_i     = &pA;    // Current search direction (initially zero)
    p_ip1   = &pB;    // Updated search direction (initially zero)

    // ===== Set-up diagnostics file
    AssimDiags = fopen (Assim_diags_filename, "w");


    // ===============================================================================
    // Inner loops (code up for conjugate gradients only for now)
    // ===============================================================================

    // Calculate the initial gradient of the cost function
    printf ("Entering PenAndGrad\n");
    PenAndGrad (&bg,              //in    background state
                &obs_trunk,       //inout observation structure
                chi_i,            //in    control variable
                true,             //in    is control variable zero?
                &MetaData,        //in    metadata
                &HorizData,       //in    info about horiz transform
                &CVTdata,         //in    info about cvt
                &Jb0,             //out   background part of cost fn
                &Jo0,             //out   observation part of cost fn
                &J0,              //out   total cost fn
                res_i,            //out   gradient in control space
                true,             //in    switch to compute gradient
                Dt,               //in    major timestep
                dt,               //in    minor timestep
                kappa_dt,         //in    diffusion timestep size
                kappa_h,          //in    horiz diffusion coefficient
                kappa_v,          //in    vert diffusion coefficient
                true,             //in    include advection
                inc_vert,         //in    include vertical transport?
                interpolation_lc, //in    interpolation type
                factor_w,         //in    mult factor of vert winds
                wind_dir,         //in    directory containing driver winds
                file_dps,         //in    in or out file for departure points
                new_dp);          //in    true if dp file needed to be made


    printf ("Writing observations (inlcuding model obs at the background)\n");
    WriteObservations ( &MetaData,
                        &obs_trunk,
                        obsfile_out_bg );
    printf ("Done writing observations\n");


    printf ("Exiting PenAndGrad\n");
    Jb     = Jb0;
    Jo     = Jo0;
    J      = J0;
    new_dp = false; // No need to make another departure points file
                    // (if it didn't exist, PenAndGrad creates it)

    // Residual is -gradient (so multiply output of PenAndGrad (res_i) by -1)
    multiply_general (res_i,
                      &HorizData,
                      true,
                      -1.0, -1.0);

    //printf ("Writing out the first residual\n");
    //WriteStateVector (res_i,
    //                  "../TestAssim/Redidual.nc");


    norm2_res_i = innerproduct_general (res_i,
                                        res_i,
                                        &HorizData,
                                        true);
    norm2_res_0 = norm2_res_i;

    // The initial search direction is the initial residual
    copy_general (res_i,
                  p_i,       // p_i set to residual
                  &HorizData,
                  true);

    loop      = 0;
    converged = false;

    printf ("\n=== INITIAL DIAGNOSTICS: %e %e %e  %e %e\n", Jb, Jo, J, norm2_res_i, sqrt(norm2_res_i));

    // Output initial data to diagnostics file
    fprintf (AssimDiags, "# Iteration, Jb, Jo, J, grad2, grad\n");
    fprintf (AssimDiags, "# ---------------------------------\n");
    fprintf (AssimDiags, "%i %f %f %f %f %f\n",
                          loop, Jb, Jo, J, norm2_res_i, sqrt(norm2_res_i));



    // ============================================================================
    //   VARIATIONAL ASSIMILATION ITERATIONS
    // ============================================================================
    do
    { loop++;
      printf ("===== Starting assimilation loop number %i\n", loop);

      // Do a line search in the direction of p_i

      // Want to make a perturbation of chi_i of a controlled (small) size
      // Perturb in the direction of p_i, but p_i is often relatively large
      // Work out how large the chi_i vector is
      norm2_chi = innerproduct_general (chi_i, chi_i,
                                        &HorizData,
                                        true);
      norm2_p_i = innerproduct_general (p_i, p_i,
                                        &HorizData,
                                        true);
      mu_modified       = sqrt(1.0/norm2_p_i);
      minus_mu_modified = -1.0 * mu_modified;
      printf ("norm2_chi   = %e\n", norm2_chi);
      printf ("norm2_p_i   = %e\n", norm2_p_i);
      printf ("mu_modified = %e\n", mu_modified);

      // Perturb chi_i positively
      add_general_fac (chi_i,
                       p_i,
                       &pert,
                       &HorizData,
                       true,
                       mu_modified);

      // Compute Jp
      PenAndGrad (&bg,              //in    background state
                  &obs_trunk,       //inout observation structure
                  &pert,            //in    control variable
                  false,            //in    is control variable zero?
                  &MetaData,        //in    metadata
                  &HorizData,       //in    info about horiz transform
                  &CVTdata,         //in    info about cvt
                  &Jbp,             //out   background part of cost fn
                  &Jop,             //out   observation part of cost fn
                  &Jp,              //out   total cost fn
                  dummy,            //out   gradient in control space
                  false,            //in    switch to compute gradient
                  Dt,               //in    major timestep
                  dt,               //in    minor timestep
                  kappa_dt,         //in    diffusion timestep size
                  kappa_h,          //in    horiz diffusion coefficient
                  kappa_v,          //in    vert diffusion coefficient
                  true,             //in    include advection
                  inc_vert,         //in    include vertical transport?
                  interpolation_lc, //in    interpolation type
                  factor_w,         //in    mult factor of vert winds
                  wind_dir,         //in    directory containing driver winds
                  file_dps,         //in    in or out file for departure points
                  new_dp);          //in    true if dp file needed to be made
      printf ("\n=== J+ DIAGNOSTICS: %e %e %e\n", Jbp, Jop, Jp);


      // Perturb chi_i negatively
      add_general_fac (chi_i,
                       p_i,
                       &pert,
                       &HorizData,
                       true,
                       minus_mu_modified);

      // Compute Jm
      PenAndGrad (&bg,              //in    background state
                  &obs_trunk,       //inout observation structure
                  &pert,            //in    control variable
                  false,            //in    is control variable zero?
                  &MetaData,        //in    metadata
                  &HorizData,       //in    info about horiz transform
                  &CVTdata,         //in    info about cvt
                  &Jbm,             //out   background part of cost fn
                  &Jom,             //out   observation part of cost fn
                  &Jm,              //out   total cost fn
                  dummy,            //out   gradient in control space
                  false,            //in    switch to compute gradient
                  Dt,               //in    major timestep
                  dt,               //in    minor timestep
                  kappa_dt,         //in    diffusion timestep size
                  kappa_h,          //in    horiz diffusion coefficient
                  kappa_v,          //in    vert diffusion coefficient
                  true,             //in    include advection
                  inc_vert,         //in    include vertical transport?
                  interpolation_lc, //in    interpolation type
                  factor_w,         //in    mult factor of vert winds
                  wind_dir,         //in    directory containing driver winds
                  file_dps,         //in    in or out file for departure points
                  new_dp);          //in    true if dp file needed to be made
      printf ("\n=== J- DIAGNOSTICS: %e %e %e\n", Jbm, Jom, Jm);

      beta = mu_modified * (Jm - Jp) / (2.0 * (Jm + Jp - 2.0*J));
      printf ("\n=== beta: %e\n", beta);

      // Update control variable
      add_general_fac (chi_i,
                       p_i,
                       chi_ip1,
                       &HorizData,
                       true,
                       beta);


      // Calculate a new gradient
      PenAndGrad (&bg,              //in    background state
                  &obs_trunk,       //inout observation structure
                  chi_ip1,          //in    control variable
                  false,            //in    is control variable zero?
                  &MetaData,        //in    metadata
                  &HorizData,       //in    info about horiz transform
                  &CVTdata,         //in    info about cvt
                  &Jb,              //out   background part of cost fn
                  &Jo,              //out   observation part of cost fn
                  &J,               //out   total cost fn
                  res_ip1,          //out   gradient in control space
                  true,             //in    switch to compute gradient
                  Dt,               //in    major timestep
                  dt,               //in    minor timestep
                  kappa_dt,         //in    diffusion timestep size
                  kappa_h,          //in    horiz diffusion coefficient
                  kappa_v,          //in    vert diffusion coefficient
                  true,             //in    include advection
                  inc_vert,         //in    include vertical transport?
                  interpolation_lc, //in    interpolation type
                  factor_w,         //in    mult factor of vert winds
                  wind_dir,         //in    directory containing driver winds
                  file_dps,         //in    in or out file for departure points
                  new_dp);          //in    true if dp file needed to be made
      printf ("\n=== Jnew DIAGNOSTICS: %e %e %e\n", Jb, Jo, J);

      // Residual is -gradient (so multiply output of PenAndGrad (res_ip1) by -1)
      multiply_general (res_ip1,
                        &HorizData,
                        true,
                        -1.0, -1.0);

      norm2_res_ip1 = innerproduct_general (res_ip1,
                                            res_ip1,
                                            &HorizData,
                                            true);
      GradRatio = norm2_res_ip1 / norm2_res_0;
      printf ("\n=== norm2, ratio: %e %e\n", norm2_res_ip1, GradRatio);



      converged = (GradRatio < conv_criterion);

      if ((loop < max_iters) && !converged)
      { // Find new search direction
        alpha = norm2_res_ip1 / norm2_res_i;
        printf ("\n=== alpha: %e\n", alpha);

        add_general_fac (res_ip1,
                         p_i,
                         p_ip1,
                         &HorizData,
                         true,
                         alpha);
        // Swap pointers for next iteration
        swap    = chi_i;
        chi_i   = chi_ip1;
        chi_ip1 = swap;

        swap    = res_i;
        res_i   = res_ip1;
        res_ip1 = swap;

        swap    = p_i;
        p_i     = p_ip1;
        p_ip1   = swap;

        norm2_res_i = norm2_res_ip1;
      }

      // Report on assimilation data
      fprintf (AssimDiags, "%i %f %f %f %f %f\n",
                            loop, Jb, Jo, J, norm2_res_i, sqrt(norm2_res_i));

    }
    while ((loop < max_iters) && !converged);


    // ============================================================================
    //   END OF VARIATIONAL ASSIMILATION ITERATIONS
    // ============================================================================


    printf ("Writing observations (inlcuding model obs at the analysis)\n");
    WriteObservations ( &MetaData,
                        &obs_trunk,
                        obsfile_out_anal );
    printf ("Done writing observations\n");



    fclose (AssimDiags);


    // Find the analysis increment
    cvt_total ( chi_ip1,
                &analinc,
                &MetaData,
                &HorizData,
                &CVTdata );

    // Add the analysis increment onto the background
    add_general (&bg,
                 &analinc,
                 &anal,
                 &HorizData,
                 true);

    // Output the analysis and the analysis increment
    printf ("Writing out the analysis\n");
    WriteStateVector (&anal,
                      anal_filename);
    printf ("Writing out the analysis increment\n");
    WriteStateVector (&analinc,
                      analinc_filename);


    // ===============================================================================
    // Tidy up
    // ===============================================================================

    printf ("Tidying up\n");
    DeallocateHorizTrans (&HorizData);
    Deallocate_CVT (&CVTdata);
    Deallocate_state (&bg);
    Destroy_Obs (&(obs_trunk.first));
    delete[] obs_trunk.mass_profile;
    Deallocate_state (&anal);
    Deallocate_state (&analinc);
    Deallocate_state (&chiA);
    Deallocate_state (&chiB);
    Deallocate_state (&resA);
    Deallocate_state (&resB);
    Deallocate_state (&pA);
    Deallocate_state (&pB);
    Deallocate_state (&pert);
    Deallocate_metadata (&MetaData);
  }
  printf ("Program completed\n");
}
