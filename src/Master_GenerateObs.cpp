/* ==================================================================================
   3d source and sink code
   Generate observations

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
    ./Master_GenerateObs.out \
      <state filename (initial and boundary conds)>\
      <wind directory>                       \
      <obs output filename>                  \
      <major timestep (between winds, s)>    \
      <minor timestep (integration, s)>      \
      <diffusion timestep (s)                \
      <horiz diffusion coefficient>          \
      <vert diffusion coefficient>           \
      <include vertical transport?>          \ 0 or 1
      <w factor>                             \
      <output frequency, seconds>            \ 
      <interpolation type (l or c)>          \ linear or cubic
      <read or write departure points file>  \ r or w
      --- Observation network specifications ---
        <What is observed>                   \ t=tracer, f=flux, x=total column tracer
                                               e=end of spec
        <Network type>                       \ i=individual ob,
                                               g=grid of obs,
                                               s=satellite
        --- If "i" -----------------------------
          <lon (deg)> <lat (deg)>, <ht (m) (tracer only)>  \
          <day> <hour> <min>                 \
          <error standard deviation>         \
        --- If "g" -----------------------------
          <start lon (deg)>                  \ west-most, deg
          <separation in lon (deg)>          \
          <number in lon direction>          \
          <start lat (deg)>                  \ south-most, deg
          <separation in lat (deg)>          \
          <number in lat direction>          \
          <start ht (m) (tracer only)>       \
          <separation in ht (m) (tracer only)>  \
          <number in ht direction (tracer only)>\
          <start time (day)> <hour> <min>    \
          <separation in time (min)>         \
          <number in time>                   \
          <error standard deviation>         \
        ----- If "s" ---------------------------------
          <orbit inclination (deg)>          \
          <orbit period (min)>               \
          <t=0 long (deg)>                   \
          <height (m) (tracer only)>         \
          <separation in time (min)>         \
          <number in time>                   \
          <error standard deviation>


    Alternative locations of winds
                         ../data/ECMWF_RedResWinds_1995 \
                         /media/ross/1297-5336/data/ECMWF_RedResWinds_1995 \


   Modification history
   --------------------
   26/07/23 New Code. Ross Bannister
   16/03/24 Allow more flexible observation network.  Ross Bannister
   04/05/25 Allow observation biases to be specified.  Ross Bannister
   30/10/25 Prevent double consideration of obs on timestep bounds.  Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )

{ struct metadata_type       MetaData;
  int                        rndseq = 7384395;   // Random number seed
  char                       state_ic[256];
  char                       wind_dir[256];
  char                       wind_file[256];
  char                       output_obs_file[256];
  char                       output_file_anim[256];
  char                       output_file_diags[256];
  char                       file_dps[256];
  char                       interpolation_lc, dp_file_rw;
  double                     Dt, dt, kappa_dt, kappa_h, kappa_v, dist_sq, max_val, min_val;
  double                     sDt, fluxfactor;
  double                     maxmax_val, minmin_val;
  double                     runlength, output_freq, factor_w;
  double                     tracer_ob_bias, flux_ob_bias, tc_ob_bias, bias;
  double                     ***dummy;
  struct state_type          state;
  struct Wind_type           windA, windB;
  struct Wind_type           *wind_lowert, *wind_uppert, *wind_temp;
  struct instant_tracer_type field1, field2;
  struct instant_tracer_type *field_lowert, *field_uppert, *field_temp;
  bool                       ok, output_diags, output_anim, mention_dps, inc_vert, inc_adv, need_winds;
  int                        nminor, gamma, epsilon, nmajor, t, windnum, starttime, op_freq, n_op_times;
  int                        x, y, z, blobx, bloby, blobz, xpos, ypos, zpos;
  int                        nss_needed, sfield;
  FILE*                      diags_file = NULL;
  struct obs_trunk_type      obs_trunk;
  struct obs_type            *currentob;
  double                     last_ob_time, obtime, lowert, uppert;
  double                     full_lev_plus1, full_lev_here, density;
  double                     ***dummyd;
  int                        ***dummyi;
  double                     ob_deviation;
  double                     mean_ob_tracer, stddev_ob_tracer;
  double                     mean_ob_xc, stddev_ob_xc;
  double                     mean_ob_flux, stddev_ob_flux;
  int                        nob_tracer, nob_flux, nob_xc;
  bool                       this_window;

//   ===============================================================================
//   Initialisation
//   ===============================================================================

  // Process command line arguments
  ok = GenerateObs_arguments (
         argument_count,    // in
         argument_list,     // in
         state_ic,          // out Filename of state (initial condition and source)
         wind_dir,          // out Directory containing driver winds
         output_obs_file,   // out Output file for synthetic observations
         file_dps,          // out Input or output file for departure points
         &Dt,               // out Major timestep (between wind files, seconds)
         &dt,               // out Minor timestep (for integration scheme, seconds)
         &kappa_dt,         // out Timestep for diffusion
         &kappa_h,          // out Horizontal diffusion coefficient
         &kappa_v,          // out Vertical diffusion coefficient
         &inc_vert,         // out Include vertical transport?
         &factor_w,         // out Multiplication factor of vertical winds
         &output_freq,      // out How often between outputs (s)
         &interpolation_lc, // out Linear or cubic interpolation
         &dp_file_rw,       // out (r)ead or (w)rite departure points file
         &tracer_ob_bias,   // out tracer observation bias
         &flux_ob_bias,     // out flux observation bias
         &tc_ob_bias,       // out total column tracer observation bias
         &(obs_trunk.first) ); // out observation network


  if (ok)
  { // How many minor timesteps per major timestep?
    gamma = int(Dt/dt);
    printf ("There are %i minor timesteps per major timestep\n", gamma);

    // Find out the latest observation time (to determine the run length)
    last_ob_time = 0.0;
    currentob    = obs_trunk.first;
    while (currentob)
    { if ((*currentob).obtime_secs > last_ob_time)
      { last_ob_time = (*currentob).obtime_secs;
      }
      currentob = (*currentob).next;
    }
    printf ("The last obs time: %f s, = %f days\n", last_ob_time, last_ob_time/(3600.0 * 24.0));
    printf ("This may need rounding up for the run length\n");
    runlength = dt * int(last_ob_time / dt + ::nearly1);
    printf ("The run length: %f s, = %f days\n", runlength, runlength/(3600.0 * 24.0));

    // How many minor timesteps in the total integration?
    nminor = int(runlength / dt);
    MetaData.ntimes_minor = nminor;
    printf ("There are %i minor timesteps in the total integration\n", nminor);

    // How many major time steps are involved in this integration?
    nmajor = int(runlength / Dt + ::nearly1);
    MetaData.ntimes_major = nmajor;
    printf ("There are %i major timesteps in this integration\n", nmajor);
    printf ("There are %i wind fields required\n", nmajor+1);

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

    // Determine the time step between source/sink fields
    sDt = state.times[1] - state.times[0];
    printf ("The time between surface flux fields in %f s (%f days)\n", sDt, sDt / (24.0 * 3600.0));

    // How many major timesteps per flux timestep?
    epsilon = int(sDt/Dt);
    printf ("There are %i major timesteps per flux timestep\n", epsilon);

    // sDt must be a multiple of Dt
    if ((sDt/Dt) != double(epsilon))
    { printf ("Error: the value of sDt (%f) must be a multiple of Dt (%f)\n", sDt, Dt);
    }

    // Calculate the flux factor (the factor for the flux term)
    fluxfactor = Dt / (::rho0 * (MetaData.level[2] - MetaData.level[1]));
    printf ("The flux factor is %f\n", fluxfactor);

    // How many surface source/sink files are needed for this integration?
    nss_needed = int(runlength / sDt + ::nearly1);
    printf ("There are %i source/sink fields needed for this integration\n", nss_needed);    
    printf ("There are %i source/sink fields in the specified file\n", state.nss);

    if (nss_needed > state.nss)
    { printf ("Error: the number of source/sink fields needed is larger than those available\n");
      exit(0);
    }

    // What is the output frequency measured in timesteps?
    op_freq = int(output_freq / Dt);
    printf ("There are %i major timesteps between each output time\n", op_freq);
    if (op_freq == 0)
    { printf ("Error: Dt (%f) must be <= output_freq (%f)\n", Dt, output_freq);
      exit(0);
    }

    // How many output times will there be (1 is added because the initial conditions are reproduced in the first)?
    n_op_times = int(runlength / output_freq) + 1;
    printf ("There are %i output times\n", n_op_times);

    inc_adv = true;  // Will run model with advection!

    // Which data do we output?
    output_anim  = false;   // Output animation file?
    output_diags = false;   // Output diagnostics file?
    mention_dps  = (strcmp ("nil", file_dps) != 0);           // Mention departure points file?

    if (output_anim)  { printf ("Will output animation file\n"); }
    if (output_diags) { printf ("Will output diagnostics file\n"); }
    if (mention_dps)  { printf ("Have mentioned departure points file\n"); }

    if ((dp_file_rw == 'r') && !mention_dps)
    { printf ("Departure points file has not been specified (is nil), but should be\n");
      exit(0);
    }

    // Do we need to read from the wind files?
    need_winds = !(mention_dps && (dp_file_rw == 'r'));


    // Allocate space for data
    if (need_winds)
    { printf ("Allocating windA\n");
      Allocate_wind (&windA,
                     &MetaData,
                     false);
      printf ("Allocating windB\n");
      Allocate_wind (&windB,
                     &MetaData,
                     false);
    }

    printf ("Allocating field1\n");
    Allocate_instant (&field1,
                      &MetaData);
    printf ("Allocating field2\n");
    Allocate_instant (&field2,
                      &MetaData);

    printf ("Done allocations\n");

    printf ("Number of levels in MetaData : %i\n", MetaData.nlev);
    printf ("Number of levels in windA    : %i\n", windA.nlev);


    // Set the initial conditions
    for (x=0; x<state.nlon+2; x++)
    { for (y=0; y<state.nlat+2; y++)
      { for (z=0; z<state.nlev+2; z++)
        { field1.tracer[x][y][z] = state.tracer0_rs[x][y][z];
        }
      }
    }
    field1.timestep_major = 0;
    field1.timestep_minor = 0;
    field1.time           = 0.0;

    // Set up the file for the departure points
    if (mention_dps && (dp_file_rw == 'w'))
    { printf ("Output file for departure points: %s\n", file_dps);
      WriteDeparturePoints (&MetaData,
                            dummyd,
                            dummyd,
                            dummyd,
                            dummyi,
                            dummyi,
                            dummyi,
                            nmajor,
                            Dt,
                            0,
                            0,       // 0=create file only,
                            file_dps);
    }


    // Maximum and minimum values
    maxmin_tracer (&MetaData,
                   &field1,
                   true,
                   &max_val,
                   &min_val);
    maxmax_val = max_val;
    minmin_val = min_val;

    // Set up the output file for the fields
    if (output_anim)
    { printf ("Output file for animation: %s\n", output_file_anim);
      WriteTimeSeq (&MetaData,
                    &field1,
                    n_op_times,
                    output_freq,
                    0,
                    0,       // 0=create file only
                    output_file_anim);

      // Write out t=0 field
      WriteTimeSeq (&MetaData,
                    &field1,
                    n_op_times,
                    output_freq,
                    0,
                    1,       // 1=normal (write only)
                    output_file_anim);
    }


    // Open output file for diagnostics
    if (output_diags)
    { printf ("Output file for diagnostics: %s\n", output_file_diags);
      diags_file = fopen (output_file_diags, "w");

      // Output the first maximum and mininmum values of the tracer field at this time
      fprintf (diags_file, "# Maximum and minimum values of the integration\n");
      fprintf (diags_file, "# dt         = %f\n", dt);
      fprintf (diags_file, "# kappa_h    = %f\n", kappa_h);
      fprintf (diags_file, "# kappa_v    = %f\n", kappa_v);
      fprintf (diags_file, "%f  %f\n", min_val, max_val);
    }

    if (need_winds)
    { // Read-in the first two windfields
      windnum = 0;
      sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
      printf ("First wind file : %s\n", wind_file);
      Read_winds (&windA, &MetaData, factor_w, wind_file);
      windnum++;
      sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
      printf ("Second wind file: %s\n", wind_file);
      Read_winds (&windB, &MetaData, factor_w, wind_file);
    }

    // Set-up the initial configuration of the winds
    wind_lowert = &windA;
    wind_uppert = &windB;

    // Configure the pointers to the tracer fields
    field_lowert = &field1;
    field_uppert = &field2;
    starttime    = 1;

//   ===============================================================================
//   Set-up the mass profile (for satellite total column observations)
//   ===============================================================================

    printf ("============================================================\n");
    printf ("Deriving density and mass (per unit area) profile as follows\n");
    obs_trunk.mass_profile = new double [MetaData.nlev+1];
    for (z=1; z<=MetaData.nlev; z++)
    { full_lev_plus1      = 0.5 * (MetaData.level[z+1] + MetaData.level[z]);
      full_lev_here       = 0.5 * (MetaData.level[z] + MetaData.level[z-1]);
      density             = ::rho0 * exp(-1.0 * MetaData.level[z-1] / ::H);
      obs_trunk.mass_profile[z]
                          = (full_lev_plus1 - full_lev_here) * density;
      printf ("%i  %f  %f\n", z, density, obs_trunk.mass_profile[z]);
    }


//   ===============================================================================
//   Start to run the integration of the transport model
//   ===============================================================================

    for (t=starttime; t<=nmajor; t++)
    { sfield = int(((double(t)-0.1) * Dt) / sDt);
      printf ("Major timestep %i of %i, source field %i\n", t, nmajor, sfield);
      // Run a model major timestep
      printf ("Running SemiLagrangian\n");
      SemiLagrangian_1step ( &MetaData,         // Meta data
                             &state,            // State vector
                             field_lowert,      // Tracer field at major t-1
                             field_uppert,      // Tracer field extrapolated to major t
                             t,                 // Major timestep number (new time)
                             sfield,            // Source/sink field number
                             Dt,                // Major timestep size
                             dt,                // Minor timestep size
                             kappa_dt,          // Diffusion timestep size
                             gamma,             // Number of minor timesteps per major timestep
                             kappa_h,           // Horiz diffusion coefficient
                             kappa_v,           // Vert diffusion coefficient
                             inc_adv,           // Include advection
                             fluxfactor,        // Flux factor
                             inc_vert,          // Include vertical transport?
                             wind_lowert,       // Wind fields at lower major timestep
                             wind_uppert,       // Wind fields at next major step
                             interpolation_lc,  // Interpolation type
                             mention_dps,       // Is departure point file mentioned?
                             file_dps,          // Filename of departure points file
                             dp_file_rw );      // Read or write departure points file?
      printf ("Done SemiLagrangian\n");

      // Apply the observation operator for observations between times t*Dt and (t+1)*Dt
      // Notes
      // At this stage, field_lowert corresponds to the tracer field at time t * Dt
      //                field_uppert                                     (t+1) * Dt
      // This state has source field between state.source_rs[.,.,sfield] and 
      //                                     state.source_rs[.,.,sfield+1]

      // Look for the relevant observations
      lowert    = double(t-1) * Dt;
      uppert    = double(t) * Dt;
      currentob = obs_trunk.first;
      printf ("Considering observations\n");
      while (currentob)
      { obtime = (*currentob).obtime_secs;
        //printf ("Checking times %f  ** %f ** %f\n", lowert, obtime, uppert);
        // Do we consider this observation in this time window?
        if (t == starttime)
        { this_window = (obtime >= lowert) && (obtime <= uppert);
        }
        else
        { this_window = (obtime > lowert) && (obtime <= uppert);
        }
        if (this_window)
        { //Generate the model's version of this observation
          //printf ("Generating observation, %c %f %f %f %f\n",
          //         (*currentob).ob_of,
          //         (*currentob).longitude,
          //         (*currentob).latitude,
          //         (*currentob).level,
          //         (*currentob).obtime_secs);
          ObservationOperator (currentob,
                               &obs_trunk,
                               &MetaData,
                               field_lowert,
                               field_uppert,
                               lowert,
                               uppert,
                               &state);

          // Determine the bias
          switch ((*currentob).ob_of)
          { case 't' : // Tracer observation
              bias = tracer_ob_bias;
            break;
            case 'f' : // Flux observation
              bias = flux_ob_bias;
            break;
            case 'x' : // Total column tracer observation
              bias = tc_ob_bias;
            break;
          }

          // Add bias and some random noise to the model's observation
          (*currentob).ob = Normal ((*currentob).model_ob + bias,
                                    (*currentob).stddev,
                                    &rndseq);

          //(*currentob).model_ob = ::notdef;
        }
        currentob = (*currentob).next;
      }


      // Output this timestep
      if (output_anim)
      { if ((t % op_freq) == 0)
        { printf ("Outputting data step %i\n", int(double(t)/double(op_freq)));
          WriteTimeSeq (&MetaData,
                        field_uppert,
                        n_op_times,
                        output_freq,
                        int(double(t)/double(op_freq)),
                        1,            // 1=normal (write only)
                        output_file_anim);
        }
      }

      // Maximum and minimum values
      maxmin_tracer (&MetaData,
                     field_uppert,
                     true,
                     &max_val,
                     &min_val);
      if (max_val > maxmax_val)
      { maxmax_val = max_val;
      }
      if (min_val < minmin_val)
      { minmin_val = min_val;
      }

      if (output_diags)
      { fprintf (diags_file, "%f  %f  %f\n", double(t) * dt, max_val, min_val);
      }


      if (t < nmajor)
      { // Bother with next winds fields only if there is still another step to do
        printf ("Just crossed a major timestep\n");
        windnum++;
        // Swap pointers
        wind_temp   = wind_lowert;
        wind_lowert = wind_uppert;
        wind_uppert = wind_temp;

        if (need_winds)
        { // Read-in the next wind field
          sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
          printf ("Wind file: %s\n", wind_file);
          Read_winds (wind_uppert, &MetaData, factor_w, wind_file);
        }

        // Sort out the configuration of the instantateous fields for the next step
        field_temp   = field_lowert;
        field_lowert = field_uppert;
        field_uppert = field_temp;
      }
    }

    if (output_diags)
    { fclose (diags_file);
    }

    printf ("Smallest value encountered = %f\n", minmin_val);
    printf ("Largest value encountered  = %f\n", maxmax_val);

//   ===============================================================================
//   Output the synthetic observations
//   ===============================================================================
    printf ("Writing observations\n");
    WriteObservations ( &MetaData,
                        &obs_trunk,
                        output_obs_file );
    printf ("Done writing observations\n");

//   ===============================================================================
//   Calculate the variability in the true model observations (stored in model_ob)
//   ===============================================================================
    printf ("Calculating the variability in the observations ...\n");
    // Calculate the means
    currentob      = obs_trunk.first;
    mean_ob_tracer = 0.0;    nob_tracer = 0;
    mean_ob_flux   = 0.0;    nob_flux   = 0;
    mean_ob_xc     = 0.0;    nob_xc     = 0;
    while (currentob)
    { switch ((*currentob).ob_of)
      { case 't' : // Tracer observation
          mean_ob_tracer += (*currentob).model_ob;
          nob_tracer     += 1;
        break;
        case 'f' : // Flux observation
          mean_ob_flux += (*currentob).model_ob;
          nob_flux     += 1;
        break;
        case 'x' : // Total column tracer observation
          mean_ob_xc += (*currentob).model_ob;
          nob_xc     += 1;
        break;
      }
      currentob = (*currentob).next;
    }

    // Normalise the means
    if (nob_tracer > 0)
    { mean_ob_tracer /= double(nob_tracer);
    }
    if (nob_flux > 0)
    { mean_ob_flux /= double(nob_flux);
    }
    if (nob_xc > 0)
    { mean_ob_xc /= double(nob_xc);
    }

    // Calculate the standard deviations
    currentob        = obs_trunk.first;
    stddev_ob_tracer = 0.0;
    stddev_ob_flux   = 0.0;
    stddev_ob_xc     = 0.0;
    while (currentob)
    { switch ((*currentob).ob_of)
      { case 't' : // Tracer observation
          ob_deviation      = (*currentob).model_ob - mean_ob_tracer;
          stddev_ob_tracer += ob_deviation * ob_deviation;
        break;
        case 'f' : // Flux observation
          ob_deviation    = (*currentob).model_ob - mean_ob_flux;
          stddev_ob_flux += ob_deviation * ob_deviation;
        break;
        case 'x' : // Total column tracer observation
          ob_deviation  = (*currentob).model_ob - mean_ob_xc;
          stddev_ob_xc += ob_deviation * ob_deviation;
        break;
      }
      currentob = (*currentob).next;
    }

    // Normalise the standard deviations
    if (nob_tracer > 1)
    { stddev_ob_tracer /= double(nob_tracer);
      stddev_ob_tracer  = sqrt(stddev_ob_tracer);
    }
    if (nob_flux > 1)
    { stddev_ob_flux /= double(nob_flux);
      stddev_ob_flux  = sqrt(stddev_ob_flux);
    }
    if (nob_xc > 1)
    { stddev_ob_xc /= double(nob_xc);
      stddev_ob_xc  = sqrt(stddev_ob_xc);
    }

    printf ("\nTracer observations\n");
    printf ("Mean   : %f\n", mean_ob_tracer);
    printf ("Stddev : %f\n", stddev_ob_tracer);
    printf ("Flux observations\n");
    printf ("Mean   : %f\n", mean_ob_flux);
    printf ("Stddev : %f\n", stddev_ob_flux);
    printf ("Total column observations\n");
    printf ("Mean   : %f\n", mean_ob_xc);
    printf ("Stddev : %f\n\n\n", stddev_ob_xc);

//   ===============================================================================
//   Deallocate
//   ===============================================================================

    Deallocate_state (&state);
    if (need_winds)
    { Deallocate_wind (&windA);
      Deallocate_wind (&windB);
    }
    Deallocate_instant (&field1, &MetaData);
    Deallocate_instant (&field2, &MetaData);
    Deallocate_metadata (&MetaData);
    delete[] obs_trunk.mass_profile;

  }


  printf ("Deallocating arguments\n");
  Destroy_Obs (&(obs_trunk.first));
  printf ("Done deallocation\n");
  printf ("Done program\n");



}
