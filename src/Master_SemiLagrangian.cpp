/* ==================================================================================
   3d source and sink code
   Run semi-Lagrangian forecast model

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_SemiLagrangian.out

   To run (e.g.)

./Master_SemiLagrangian.out <state filename (initial and boundary conds)>\
                            <wind directory>                             \
                            <major timestep (between winds, seconds)>    \
                            <minor timestep (integration, seconds)>      \
                            <diffusion timestep (seconds)                \
                            <horiz diffusion coefficient>                \
                            <vert diffusion coefficient>                 \
                            <include vertical transport?>                \ 0 or 1
                            <w factor>                                   \
                            <run length (days)>                          \ Must be multiple of major timestep
                            <output frequency (seconds)>                 \
                            <interpolation type (l or c)>                \ linear or cubic
                            <animation output filename>                  \ nil to not output anim file
                            <diagnostics output filename>                \ nil to not output diagnostics file
                            <departure points filename>                  \ Dep points filename (e.g. for the adjoint run), nil to not output
                            <create and write dp file, or read existing one (w or r)>
                                                                           

    Alternative locations of winds
                         ../data/ECMWF_RedResWinds_1995 \
                         /media/ross/1297-5336/data/ECMWF_RedResWinds_1995 \


   Modification history
   --------------------
   30/07/22 New Code. Ross Bannister
   31/10/22 Add department points code. Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )

{ struct metadata_type       MetaData;
  char                       state_ic[256];
  char                       wind_dir[256];
  char                       wind_file[256];
  char                       output_file_anim[256];
  char                       output_file_diags[256];
  char                       file_dps[256];
  char                       interpolate_lc, dp_file_rw;
  double                     Dt, dt, kappa_dt, kappa_h, kappa_v, dist_sq, max_val, min_val;
  double                     sDt, fluxfactor;
  double                     maxmax_val, minmin_val;
  double                     runlength, output_freq, factor_w;
  double                     ***dummyd = NULL;
  int                        ***dummyi = NULL;
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


//   ===============================================================================
//   Initialisation
//   ===============================================================================

  // Process command line arguments
  ok = SemiLagrangian_arguments (argument_count,
                                 argument_list,
                                 state_ic,
                                 wind_dir,
                                 &Dt,
                                 &dt,
                                 &kappa_dt,
                                 &kappa_h,
                                 &kappa_v,
                                 &inc_vert,
                                 &factor_w,
                                 &runlength,
                                 &output_freq,
                                 &interpolate_lc,
                                 output_file_anim,
                                 output_file_diags,
                                 file_dps,
                                 &dp_file_rw);


  inc_adv = true;  // Will run model with advection!

  if (ok)
  { // How many minor timesteps per major timestep?
    gamma = int(Dt/dt);
    printf ("There are %i minor timesteps per major timestep\n", gamma);

    // How many minor timesteps in the total integration?
    nminor = int(runlength / dt);
    MetaData.ntimes_minor = nminor;
    printf ("There are %i minor timesteps in the total integration\n", nminor);

    // How many major time steps are involved in this integration?
    nmajor = int(runlength / Dt);
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

    // Which data do we output?
    output_anim  = (strcmp ("nil", output_file_anim) != 0);   // Output animation file?
    output_diags = (strcmp ("nil", output_file_diags) != 0);  // Output diagnostics file?
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
//   Start to run the integration of the transport model
//   ===============================================================================


    for (t=starttime; t<=nmajor; t++)
    { sfield = int(((double(t)-0.1) * Dt) / sDt);
      printf ("Major timestep %i of %i, source field %i\n", t, nmajor, sfield);
      // Run a model major timestep
      // MODIFY for possible reading of departure points file
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
                             interpolate_lc,    // Interpolation type
                             mention_dps,       // Departure points file specified
                             file_dps,          // Filename of departure points file
                             dp_file_rw );      // Whether to (r)ead or (w)rite departure points file


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
    printf ("Done program\n");
  }

}
