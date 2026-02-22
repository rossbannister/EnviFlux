   #include <source.h>


// -------------------------------------------------------------------------------
void SemiLagrangian_1step ( struct metadata_type       *MetaData,       // Meta data
                            struct state_type          *state,          // In:  State vector
                            struct instant_tracer_type *state_prev,     // In:  Tracer field at t-1
                            struct instant_tracer_type *state_next,     // Out: Tracer field extrapolated to t
                            int                        major_timestep,  // Major timestep number
                            int                        sfield,          // source/sink field number
                            double                     Dt,              // Major timestep size
                            double                     dt,              // Minor timestep size
                            double                     kappa_dt,        // Timestep for diffusion
                            int                        gamma,           // Number of minor timesteps per major timestep
                            double                     kappa_h,         // Horiz diffusion coefficient
                            double                     kappa_v,         // Vert diffusion coefficient
                            bool                       adv_switch,      // Advection on/off
                            double                     fluxfactor,      // Flux factor
                            bool                       inc_vert,        // Include vertical transport
                            struct Wind_type           *windA,          // Wind fields at lower major timestep
                            struct Wind_type           *windB,          // Wind fields at next major step
                            char                       interpolate_lc,  // Linear or cubic interpolation
                            bool                       mention_dps,     // Departure points file mentioned?
                            char                       file_dps[256],   // Departure points output file
                            char                       dp_file_rw )     // (r)ead or (w)rite departure points file?

{ // Run the semi-Lagrangian transport model for one major time step (multiple minor timesteps)
  int                        nlon, nlat, nlev;
  int                        x, y, z, tminor, d, xm1, ym1, zm1, xp1, yp1, zp1;
  int                        lower_lev, upper_lev;
  int                        index_lon, index_lat, index_lev;
  int                        t, ndiff;
  double                     particle_lon, particle_lat, particle_lev;
  double                     particle_lon_new, particle_lat_new, particle_lev_new;
  double                     mdt;
  // Variables for linear interpolation
  double                     lon_linear[2], lat_linear[2], lev_linear[2], tracer_linear[2][2][2];
  // Variables for cubic interpolation
  double                     lon_cubic[6], lat_cubic[6], lev_cubic[6], tracer_cubic[6][6][6];
  // Variables to hold the departure points, and their indices
  double                     ***dp_lon = NULL, ***dp_lat = NULL, ***dp_lev = NULL;
  int                        ***dp_index_lon = NULL, ***dp_index_lat = NULL, ***dp_index_lev = NULL;
  // Intermediate states (between SL and diffusion)
  struct instant_tracer_type state_interA, state_interB;
  // For use with diffusion
  struct instant_tracer_type *diff_input_state = NULL, *diff_output_state = NULL, *diff_temp_state = NULL;
  // Variables to do with the diffusion
  double                     Dx_eq, Dz_up, Dz_dn, Dy_N, Dy_S, Dx_recip;
  double                     *Dz = NULL, *Dz_up_recip = NULL, *Dz_dn_recip = NULL, *Dz_recip = NULL;
  double                     *Dy = NULL, *Dy_N_recip = NULL, *Dy_S_recip = NULL, *Dy_recip = NULL, *Dx = NULL, *Dx_recip2 = NULL;
  double                     tracer_mid;
  double                     diff_h = 0.0, diff_v = 0.0;
  bool                       read_ok;


  // State between SL and diffusion
  Allocate_instant (&state_interA,
                    MetaData);

  nlon = (*state).nlon;
  nlat = (*state).nlat;
  nlev = (*state).nlev;

  mdt = -1.0 * dt;

  if (inc_vert)
  { lower_lev = 1;
    upper_lev = nlev;
  }
  else
  { lower_lev = 1;
    upper_lev = 1;
    kappa_v   = 0.0;  // Can't do vertical diffusion if only one level
  }

  if (mention_dps)
  { // Set-up arrays for departure points and their indices
    Array_3d_double_create  (&dp_lon,
                             nlon, nlat, nlev);
    Array_3d_double_create  (&dp_lat,
                             nlon, nlat, nlev);
    Array_3d_double_create  (&dp_lev,
                             nlon, nlat, nlev);
    Array_3d_int_create  (&dp_index_lon,
                          nlon, nlat, nlev);
    Array_3d_int_create  (&dp_index_lat,
                          nlon, nlat, nlev);
    Array_3d_int_create  (&dp_index_lev,
                          nlon, nlat, nlev);
  }


  if (adv_switch)
  { // Do advection

    if (mention_dps && (dp_file_rw == 'r'))
    { // Read the departure points and index fields
      //printf ("Reading departure points for this timestep\n");
      ReadDeparturePoints ( MetaData,
                            dp_lon,
                            dp_lat,
                            dp_lev,
                            dp_index_lon,
                            dp_index_lat,
                            dp_index_lev,
                            major_timestep-1,     // Timestep number
                            file_dps,
                            &read_ok);
    }

    // Loop over arrival points
    for (x=1; x<=nlon; x++)
    { xm1 = x - 1;
      for (y=1; y<=nlat; y++)
      { ym1 = y - 1;
        for (z=lower_lev; z<=upper_lev; z++)
        { zm1 = z - 1;

          if (mention_dps && (dp_file_rw == 'r'))
          { // Take the departure points and their indices from pre-calculated values (read in earlier)
            particle_lon = dp_lon[xm1][ym1][zm1];
            particle_lat = dp_lat[xm1][ym1][zm1];
            particle_lev = dp_lev[xm1][ym1][zm1];
            index_lon    = dp_index_lon[xm1][ym1][zm1];
            index_lat    = dp_index_lat[xm1][ym1][zm1];
            index_lev    = dp_index_lev[xm1][ym1][zm1];
          }
          else
          { // Compute the departure points and their indices
            particle_lon = (*MetaData).longitude[x];
            particle_lat = (*MetaData).latitude[y];
            particle_lev = (*MetaData).level[z];

            // Loop backwards over the minor timesteps
            for (tminor=gamma; tminor>0; tminor--)
            { // Find new particle longitude, latitude and level from back-trajectory
              trajectory ( tminor,              // Minor timestep No., local time (relative to windA)
                           mdt,                 // Minor timestep size
                           Dt,                  // Major timestep size
                           windA,               // Wind fields at start of major timestep
                           windB,               // Wind fields at end of major timestep
                           particle_lon,        // particle longitude at start time
                           particle_lat,        // particle latitude at start time
                           particle_lev,        // particle level at start time
                           &particle_lon_new,   // particle longitude at end time (backwards in time)
                           &particle_lat_new,   // particle latitude at end time (backwards in time)
                           &particle_lev_new,   // particle level at end time (backwards in time)
                           false );             // Don't show debugging diagnostics

              // Set the new particle position
              particle_lon = particle_lon_new;
              particle_lat = particle_lat_new;
              particle_lev = particle_lev_new;
            }

            // Now we have the departure point

            // Find the indices of the longitude, latitude, and level
            index_lon = Find_index_ascend (nlon+2,
                                           (*MetaData).longitude,
                                           particle_lon_new);
            index_lat = Find_index_descend (nlat+2,
                                           (*MetaData).latitude,
                                           particle_lat_new);
            index_lev = Find_index_ascend (nlev+2,
                                           (*MetaData).level,
                                           particle_lev_new);

            // Store the departure points and their indices
            if (mention_dps && (dp_file_rw == 'w'))
            { dp_lon[xm1][ym1][zm1]       = particle_lon;
              dp_lat[xm1][ym1][zm1]       = particle_lat;
              dp_lev[xm1][ym1][zm1]       = particle_lev;
              dp_index_lon[xm1][ym1][zm1] = index_lon;
              dp_index_lat[xm1][ym1][zm1] = index_lat;
              dp_index_lev[xm1][ym1][zm1] = index_lev;
            }
          }
 
          //printf ("====\n%f %i\n", particle_lon_new, index_lon);
          //printf ("%f %i\n", particle_lat_new, index_lat);
          //printf ("%f %i\n", particle_lev_new, index_lev);
          //if (index_lon <= 0) {index_lon = 1;}
          //if (index_lat <= 0) {index_lat = 1;}
          //if (index_lev <= 0) {index_lev = 1;}


          // Interpolate to the departure point - this can be done with either linear or cubic interpolation
          if (interpolate_lc == 'l')
          { // Linear interpolation
            // --------------------
            // Setup short dimensions
            for (d=0; d<2; d++)
            { lon_linear[d] = (*MetaData).longitude[index_lon+d];
              lat_linear[d] = (*MetaData).latitude[index_lat+d];
              lev_linear[d] = (*MetaData).level[index_lev+d];
            }
            // Prepare cube of data
            Fill_short_array ((*state_prev).tracer, tracer_linear,
                              index_lon, index_lat, index_lev);

            // Interpolate
            state_interA.tracer[x][y][z] = Interpolate3D ( tracer_linear,   // The 3D field to be interpolated
                                                           lon_linear,      // The axis x values
                                                           lat_linear,      // The axis y values
                                                           lev_linear,      // The axis z values
                                                           particle_lon,    // The x point to interpolate to
                                                           particle_lat,    // The y point to interpolate to
                                                           particle_lev );  // The z point to interpolate to
          }
          else
          { // Cubic interpolation
            // -------------------
            // Setup short dimensions
            Fill_short_dims_cubic (nlon, nlat, nlev,
                                   (*MetaData).longitude, (*MetaData).latitude, (*MetaData).level,
                                   lon_cubic, lat_cubic, lev_cubic,
                                   index_lon, index_lat, index_lev);
            // Prepare cube of data
            Fill_short_array_cubic (nlon, nlat, nlev,
                                    (*state_prev).tracer, tracer_cubic,
                                    index_lon, index_lat, index_lev);

            // Interpolate
            state_interA.tracer[x][y][z] = interpolate_3d_cubic ( lon_cubic, lat_cubic, lev_cubic, // The axis values
                                                                  tracer_cubic,                    // The 3D field to be interpolated
                                                                  particle_lon,                    // The x point to interpolate to
                                                                  particle_lat,                    // The y point to interpolate to
                                                                  particle_lev,                    // The x point to interpolate to
                                                                  3 );                             // Use three separate formula (as Laczos)
          }

        }
      }
    }
  }
  else
  { // No advection
    for (x=1; x<=nlon; x++)
    { for (y=1; y<=nlat; y++)
      { for (z=lower_lev; z<=upper_lev; z++)
        { state_interA.tracer[x][y][z] = (*state_prev).tracer[x][y][z];
        }
      }
    }
  }


  // Sort out the halos
  halos_tracer (&state_interA,
                MetaData);


  if ( (kappa_v != 0.0) || (kappa_h != 0.0))
  { // Include horizontal and vertical diffusion

    // Allocate and determine variables to do with grid spacings
    Array_1d_double_create (&Dz, nlev+2);
    Array_1d_double_create (&Dz_up_recip, nlev+2);
    Array_1d_double_create (&Dz_dn_recip, nlev+2);
    Array_1d_double_create (&Dz_recip, nlev+2);
    Array_1d_double_create (&Dy, nlat+2);
    Array_1d_double_create (&Dy_N_recip, nlat+2);
    Array_1d_double_create (&Dy_S_recip, nlat+2);
    Array_1d_double_create (&Dy_recip, nlat+2);
    Array_1d_double_create (&Dx, nlat+2);
    Array_1d_double_create (&Dx_recip2, nlat+2);

    for (z=lower_lev-1; z<=upper_lev-1; z++)
    { zp1 = z + 1;
      // Set some parameters to do with vertical diffusion
      Dz_up            = (*MetaData).level[zp1+1] - (*MetaData).level[zp1];
      Dz_dn            = (*MetaData).level[zp1] - (*MetaData).level[z];
      Dz[zp1]          = (Dz_up + Dz_dn) / 2.0;
      Dz_up_recip[zp1] = 1.0 / Dz_up;
      Dz_dn_recip[zp1] = 1.0 / Dz_dn;
      Dz_recip[zp1]    = 1.0 / Dz[zp1];
    }
    Dx_eq = ::Re * ::deg2rad * ((*MetaData).longitude[2] - (*MetaData).longitude[1]);
    for (y=0; y<nlat; y++)
    { yp1 = y + 1;
      // Set some parameters to do with horizontal diffusion
      Dy_N            = ::Re * ::deg2rad * ( (*MetaData).latitude[y] - (*MetaData).latitude[yp1] );
      Dy_S            = ::Re * ::deg2rad * ( (*MetaData).latitude[yp1] - (*MetaData).latitude[yp1+1] );
      Dy[yp1]         = (Dy_N + Dy_S) / 2.0;
      Dy_N_recip[yp1] = 1.0 / Dy_N;
      Dy_S_recip[yp1] = 1.0 / Dy_S;
      Dy_recip[yp1]   = 1.0 / Dy[yp1];
      Dx[yp1]         = Dx_eq * (*MetaData).cos_u_lat[yp1];
      Dx_recip        = 1.0 / Dx[yp1];
      Dx_recip2[yp1]  = Dx_recip * Dx_recip;
    }

    // Report on the stability condition for the diffusion
    if ((major_timestep == 1) && (verbose_output == 'v'))
    { printf ("Diffusion stability conditions (must be <= 1)\n");
      printf ("Longitudinal, equator : %f\n", fabs(2.0 * kappa_h * kappa_dt / (Dx_eq * Dx_eq)));
      printf ("Longitudinal, pole    : %f\n", fabs(2.0 * kappa_h * kappa_dt * Dx_recip2[1]));
      printf ("Latitudinal, equator  : %f\n", fabs(2.0 * kappa_h * kappa_dt * Dy_recip[nlat/2] * Dy_recip[nlat/2]));
      printf ("Latitudinal, pole     : %f\n", fabs(2.0 * kappa_h * kappa_dt * Dy_recip[1] * Dy_recip[1]));
      printf ("Vertical (ground)     : %f\n", fabs(2.0 * kappa_v * kappa_dt * Dz_recip[1]* Dz_recip[1]));
      printf ("Vertical (ceiling)    : %f\n", fabs(2.0 * kappa_v * kappa_dt * Dz_recip[nlev] * Dz_recip[nlev]));
    }

    // Allocate space for another intermediate state
    Allocate_instant (&state_interB,
                      MetaData);

    // Determine the number of diffusion steps
    ndiff = int(Dt / kappa_dt);

    // Set up the state pointers
    diff_input_state  = &state_interA;
    diff_output_state = &state_interB;

    // Time steps for diffusion
    for (t=0; t<ndiff; t++)
    { for (z=lower_lev-1; z<=upper_lev-1; z++)
      { zp1 = z + 1;
        for (y=0; y<nlat; y++)
        { yp1 = y + 1;
          for (x=0; x<nlon; x++)
          { xp1 = x + 1;
            // This is used a lot
            tracer_mid = (*diff_input_state).tracer[xp1][yp1][zp1];

            // Horizontal diffusion
            if (kappa_h != 0.0)
            { diff_h = kappa_h *
                       ( Dx_recip2[yp1] * ((*diff_input_state).tracer[xp1+1][yp1][zp1] +
                                           (*diff_input_state).tracer[x][yp1][zp1] -
                                           2.0 * tracer_mid) +
                         Dy_recip[yp1]  * ((*MetaData).cos_v_lat[yp1] * Dy_N_recip[yp1] *
                                             ((*diff_input_state).tracer[xp1][y][zp1] - tracer_mid) -
                                           (*MetaData).cos_v_lat[yp1+1] * Dy_S_recip[yp1] *
                                             (tracer_mid - (*diff_input_state).tracer[xp1][yp1+1][zp1])) /
                                          (*MetaData).cos_u_lat[yp1] );
            }

            // Vertical diffusion
            if (kappa_v != 0.0)
            { diff_v = kappa_v * Dz_recip[zp1] *
                       ( Dz_up_recip[zp1] * ((*diff_input_state).tracer[xp1][yp1][zp1+1] - tracer_mid) -
                         Dz_dn_recip[zp1] * (tracer_mid - (*diff_input_state).tracer[xp1][yp1][z]) );
            }

            (*diff_output_state).tracer[xp1][yp1][zp1] = (*diff_input_state).tracer[xp1][yp1][zp1] +
                                                         kappa_dt * (diff_h + diff_v);

          }
        }
      }

      // Sort out the halos (all but the last step)
      if (t < (ndiff-1))
      { halos_tracer (diff_output_state,
                      MetaData);
      }

      // Swap states
      diff_temp_state   = diff_input_state;
      diff_input_state  = diff_output_state;
      diff_output_state = diff_temp_state;
    }

    // Copy result to output
    for (x=1; x<=nlon; x++)
    { for (y=1; y<=nlat; y++)
      { for (z=lower_lev; z<=upper_lev; z++)
        { (*state_next).tracer[x][y][z] = (*diff_input_state).tracer[x][y][z];
        }
      }
    }

    Deallocate_instant (&state_interB, MetaData);

    Array_1d_double_destroy (&Dz);
    Array_1d_double_destroy (&Dz_up_recip);
    Array_1d_double_destroy (&Dz_dn_recip);
    Array_1d_double_destroy (&Dz_recip);
    Array_1d_double_destroy (&Dy);
    Array_1d_double_destroy (&Dy_N_recip);
    Array_1d_double_destroy (&Dy_S_recip);
    Array_1d_double_destroy (&Dy_recip);
    Array_1d_double_destroy (&Dx);
    Array_1d_double_destroy (&Dx_recip2);

  }
  else
  { // No diffusion, so just copy output from SL scheme
    // Copy result to output
    for (x=1; x<=nlon; x++)
    { for (y=1; y<=nlat; y++)
      { for (z=lower_lev; z<=upper_lev; z++)
        { (*state_next).tracer[x][y][z] = state_interA.tracer[x][y][z];
        }
      }
    }
  }



  // Add on surface flux
  for (x=1; x<=nlon; x++)
  { for (y=1; y<=nlat; y++)
    { (*state_next).tracer[x][y][1] += fluxfactor * (*state).source_rs[x][y][sfield];
    }
  }


  // Sort out the halos
  halos_tracer (state_next,
                MetaData);



  Deallocate_instant (&state_interA, MetaData);

  if (mention_dps && (dp_file_rw == 'w'))
  { // Write the departure point fields
    if (verbose_output == 'v')
    { printf ("Writing departure points for this timestep\n");
    }
    WriteDeparturePoints (MetaData,
                          dp_lon,
                          dp_lat,
                          dp_lev,
                          dp_index_lon,
                          dp_index_lat,
                          dp_index_lev,
                          1,                // Dummy variable, not used in this call
                          Dt,
                          major_timestep-1, // Minus 1 because we want to output at the time of departure
                          1,                // 1=Output data to pre-created file
                          file_dps);
  }

  if (mention_dps)
  { // Delete arrays for departure points and their indices
    Array_3d_double_destroy (&dp_lon,
                             nlon, nlat);
    Array_3d_double_destroy (&dp_lat,
                             nlon, nlat);
    Array_3d_double_destroy (&dp_lev,
                             nlon, nlat);
    Array_3d_int_destroy (&dp_index_lon,
                          nlon, nlat);
    Array_3d_int_destroy (&dp_index_lat,
                          nlon, nlat);
    Array_3d_int_destroy (&dp_index_lev,
                          nlon, nlat);
  }

}





// -------------------------------------------------------------------------------
void SemiLagrangian_1step_adj ( struct metadata_type       *MetaData,         // Meta data
                                struct state_type          *state_hat,        // In/Out: State vector
                                struct instant_tracer_type *state_prev_hat,   // In/Out: Tracer field at t-1
                                struct instant_tracer_type *state_next_hat,   // In: Tracer field extrapolated to t
                                int                        major_timestep,    // Major timestep number
                                int                        sfield,            // source/sink field number
                                double                     Dt,                // Major timestep size
                                double                     kappa_dt,          // Timestep for diffusion
                                double                     kappa_h,           // Horiz diffusion coefficient
                                double                     kappa_v,           // Vert diffusion coefficient
                                bool                       adv_switch,        // Advection on/off
                                double                     fluxfactor,        // Flux factor
                                bool                       inc_vert,          // Include vertical transport
                                char                       interpolate_lc,    // Linear or cubic interpolation
                                char                       input_file_dps[256]) // Departure points input file

{ // Run the adjoint of the semi-Lagrangian transport model for one major time step (multiple minor timesteps)
  int                        nlon, nlat, nlev;
  int                        x, y, z, d, xm1, ym1, zm1, xp1, yp1, zp1;
  int                        lower_lev, upper_lev;
  int                        index_lon, index_lat, index_lev;
  int                        t, ndiff;
  double                     particle_lon, particle_lat, particle_lev;
  // Variables for linear interpolation
  double                     lon_linear[2], lat_linear[2], lev_linear[2], tracer_linear_hat[2][2][2];
  // Variables for cubic interpolation
  double                     lon_cubic[6], lat_cubic[6], lev_cubic[6], tracer_cubic_hat[6][6][6];
  // Variables to hold the departure points, and their indices
  double                     ***dp_lon = NULL, ***dp_lat = NULL, ***dp_lev = NULL;
  int                        ***dp_index_lon = NULL, ***dp_index_lat = NULL, ***dp_index_lev = NULL;
  // Intermediate states (between SL and diffusion)
  struct instant_tracer_type state_interA_hat, state_interB_hat;
  // For use with diffusion
  struct instant_tracer_type *diff_input_state_hat = NULL, *diff_output_state_hat = NULL, *diff_temp_state_hat = NULL;
  // Variables to do with the diffusion
  double                     Dx_eq, Dz_up, Dz_dn, Dy_N, Dy_S, Dx_recip;
  double                     *Dz = NULL, *Dz_up_recip = NULL, *Dz_dn_recip = NULL, *Dz_recip = NULL;
  double                     *Dy = NULL, *Dy_N_recip = NULL, *Dy_S_recip = NULL, *Dy_recip = NULL, *Dx = NULL, *Dx_recip2 = NULL;
  double                     tracer_mid_hat;
  double                     diff_h_hat, diff_v_hat;
  bool                       read_ok;


  nlon = (*state_hat).nlon;
  nlat = (*state_hat).nlat;
  nlev = (*state_hat).nlev;

  if (inc_vert)
  { lower_lev = 1;
    upper_lev = nlev;
  }
  else
  { lower_lev = 1;
    upper_lev = 1;
    kappa_v   = 0.0;  // Can't do vertical diffusion if only one level
  }

  // Set-up arrays for departure points
  Array_3d_double_create  (&dp_lon,
                           nlon, nlat, nlev);
  Array_3d_double_create  (&dp_lat,
                           nlon, nlat, nlev);
  Array_3d_double_create  (&dp_lev,
                           nlon, nlat, nlev);
  Array_3d_int_create  (&dp_index_lon,
                        nlon, nlat, nlev);
  Array_3d_int_create  (&dp_index_lat,
                        nlon, nlat, nlev);
  Array_3d_int_create  (&dp_index_lev,
                        nlon, nlat, nlev);

  // Read-in the departure points and index fields
  //printf ("Reading departure points for this timestep (adjoint)\n");
  ReadDeparturePoints ( MetaData,
                        dp_lon,
                        dp_lat,
                        dp_lev,
                        dp_index_lon,
                        dp_index_lat,
                        dp_index_lev,
                        major_timestep-1,     // Timestep number
                        input_file_dps,
                        &read_ok);


  // State between SL and diffusion
  Allocate_instant (&state_interA_hat,
                    MetaData);


  // Sort out the halos
  halos_tracer_adj (state_next_hat,
                    MetaData);

  // Add on surface flux
  for (x=1; x<=nlon; x++)
  { for (y=1; y<=nlat; y++)
    { (*state_hat).source_rs[x][y][sfield] += fluxfactor * (*state_next_hat).tracer[x][y][1];
    }
  }



  if ( (kappa_v != 0.0) || (kappa_h != 0.0))
  { // Include horizontal and vertical diffusion

    Allocate_instant (&state_interB_hat,
                      MetaData);

    // Allocate and determine variables to do with grid spacings
    Array_1d_double_create (&Dz, nlev+2);
    Array_1d_double_create (&Dz_up_recip, nlev+2);
    Array_1d_double_create (&Dz_dn_recip, nlev+2);
    Array_1d_double_create (&Dz_recip, nlev+2);
    Array_1d_double_create (&Dy, nlat+2);
    Array_1d_double_create (&Dy_N_recip, nlat+2);
    Array_1d_double_create (&Dy_S_recip, nlat+2);
    Array_1d_double_create (&Dy_recip, nlat+2);
    Array_1d_double_create (&Dx, nlat+2);
    Array_1d_double_create (&Dx_recip2, nlat+2);

    for (z=lower_lev-1; z<=upper_lev-1; z++)
    { zp1 = z + 1;
      // Set some parameters to do with vertical diffusion
      Dz_up            = (*MetaData).level[zp1+1] - (*MetaData).level[zp1];
      Dz_dn            = (*MetaData).level[zp1] - (*MetaData).level[z];
      Dz[zp1]          = (Dz_up + Dz_dn) / 2.0;
      Dz_up_recip[zp1] = 1.0 / Dz_up;
      Dz_dn_recip[zp1] = 1.0 / Dz_dn;
      Dz_recip[zp1]    = 1.0 / Dz[zp1];
    }
    Dx_eq = ::Re * ::deg2rad * ((*MetaData).longitude[2] - (*MetaData).longitude[1]);
    for (y=0; y<nlat; y++)
    { yp1 = y + 1;
      // Set some parameters to do with horizontal diffusion
      Dy_N            = ::Re * ::deg2rad * ( (*MetaData).latitude[y] - (*MetaData).latitude[yp1] );
      Dy_S            = ::Re * ::deg2rad * ( (*MetaData).latitude[yp1] - (*MetaData).latitude[yp1+1] );
      Dy[yp1]         = (Dy_N + Dy_S) / 2.0;
      Dy_N_recip[yp1] = 1.0 / Dy_N;
      Dy_S_recip[yp1] = 1.0 / Dy_S;
      Dy_recip[yp1]   = 1.0 / Dy[yp1];
      Dx[yp1]         = Dx_eq * (*MetaData).cos_u_lat[yp1];
      Dx_recip        = 1.0 / Dx[yp1];
      Dx_recip2[yp1]  = Dx_recip * Dx_recip;
    }

    // Set up the state pointers
    diff_input_state_hat  = &state_interA_hat;
    diff_output_state_hat = &state_interB_hat;

    // Copy result to output
    for (x=1; x<=nlon; x++)
    { for (y=1; y<=nlat; y++)
      { for (z=lower_lev; z<=upper_lev; z++)
        { (*diff_input_state_hat).tracer[x][y][z] = (*state_next_hat).tracer[x][y][z];
        }
      }
    }


    // Determine the number of diffusion steps
    ndiff = int(Dt / kappa_dt);

    // Time steps for diffusion
    for (t=ndiff-1; t>=0; t--)
    { // Swap states
      diff_temp_state_hat   = diff_input_state_hat;
      diff_input_state_hat  = diff_output_state_hat;
      diff_output_state_hat = diff_temp_state_hat;


      // Sort out the halos (not the first step)
      if (t < (ndiff-1))
      { halos_tracer_adj (diff_output_state_hat,
                          MetaData);
      }


      // Initialise *diff_input_state_hat to zero
      for (z=lower_lev-1; z<=upper_lev+1; z++)
      { for (y=0; y<=nlat+1; y++)
        { for (x=0; x<=nlon+1; x++)
          { (*diff_input_state_hat).tracer[x][y][z] = 0.0;
          }
        }
      }


      for (z=lower_lev-1; z<=upper_lev-1; z++)
      { zp1 = z + 1;
        for (y=0; y<nlat; y++)
        { yp1 = y + 1;
          for (x=0; x<nlon; x++)
          { xp1 = x + 1;

            tracer_mid_hat = 0.0;

            // This is a copy of the forward code
            //(*diff_output_state).tracer[xp1][yp1][zp1] = (*diff_input_state).tracer[xp1][yp1][zp1] +
            //                                             kappa_dt * (diff_h + diff_v);
            (*diff_input_state_hat).tracer[xp1][yp1][zp1] += (*diff_output_state_hat).tracer[xp1][yp1][zp1];
            diff_h_hat                                     = (*diff_output_state_hat).tracer[xp1][yp1][zp1] * kappa_dt;
            diff_v_hat                                     = (*diff_output_state_hat).tracer[xp1][yp1][zp1] * kappa_dt;

            // Vertical diffusion
            if (kappa_v != 0.0)
            { // This is a copy of the forward code
              //diff_v = kappa_v * Dz_recip[zp1] *
              //         ( Dz_up_recip[zp1] * ((*diff_input_state).tracer[xp1][yp1][zp1+1] - tracer_mid) -
              //           Dz_dn_recip[zp1] * (tracer_mid - (*diff_input_state).tracer[xp1][yp1][z]) );

             (*diff_input_state_hat).tracer[xp1][yp1][zp1+1] += kappa_v * Dz_recip[zp1] * Dz_up_recip[zp1] * diff_v_hat;
             tracer_mid_hat                                  -= kappa_v * Dz_recip[zp1] * Dz_up_recip[zp1] * diff_v_hat;
             tracer_mid_hat                                  -= kappa_v * Dz_recip[zp1] * Dz_dn_recip[zp1] * diff_v_hat;
             (*diff_input_state_hat).tracer[xp1][yp1][z]     += kappa_v * Dz_recip[zp1] * Dz_dn_recip[zp1] * diff_v_hat;
            }


            // Horizontal diffusion
            if (kappa_h != 0.0)
            { // This is a copy of the forward code
              //diff_h = kappa_h *
              //         ( Dx_recip2[yp1] * ((*diff_input_state).tracer[xp1+1][yp1][zp1] +
              //                             (*diff_input_state).tracer[x][yp1][zp1] -
              //                             2.0 * tracer_mid) +
              //           Dy_recip[yp1]  * ((*MetaData).cos_v_lat[yp1] * Dy_N_recip[yp1] *
              //                               ((*diff_input_state).tracer[xp1][y][zp1] - tracer_mid) -
              //                             (*MetaData).cos_v_lat[yp1+1] * Dy_S_recip[yp1] *
              //                               (tracer_mid - (*diff_input_state).tracer[xp1][yp1+1][zp1])) /
              //                            (*MetaData).cos_u_lat[yp1] );

              (*diff_input_state_hat).tracer[xp1+1][yp1][zp1] += kappa_h * Dx_recip2[yp1] * diff_h_hat;
              (*diff_input_state_hat).tracer[x][yp1][zp1]     += kappa_h * Dx_recip2[yp1] * diff_h_hat;
              tracer_mid_hat                                  -= 2.0 * kappa_h * Dx_recip2[yp1] * diff_h_hat;
              (*diff_input_state_hat).tracer[xp1][y][zp1]     += kappa_h * Dy_recip[yp1] * (*MetaData).cos_v_lat[yp1] *
                                                                 Dy_N_recip[yp1] * diff_h_hat / (*MetaData).cos_u_lat[yp1];
              tracer_mid_hat                                  -= kappa_h * Dy_recip[yp1] * (*MetaData).cos_v_lat[yp1] *
                                                                 Dy_N_recip[yp1] * diff_h_hat / (*MetaData).cos_u_lat[yp1];
              tracer_mid_hat                                  -= kappa_h * Dy_recip[yp1] * (*MetaData).cos_v_lat[yp1+1] *
                                                                 Dy_S_recip[yp1] * diff_h_hat / (*MetaData).cos_u_lat[yp1];
              (*diff_input_state_hat).tracer[xp1][yp1+1][zp1] += kappa_h * Dy_recip[yp1] * (*MetaData).cos_v_lat[yp1+1] *
                                                                 Dy_S_recip[yp1] * diff_h_hat / (*MetaData).cos_u_lat[yp1];

            }

            // This is used a lot
            (*diff_input_state_hat).tracer[xp1][yp1][zp1] += tracer_mid_hat;

          }
        }
      }

    }

    Array_1d_double_destroy (&Dz);
    Array_1d_double_destroy (&Dz_up_recip);
    Array_1d_double_destroy (&Dz_dn_recip);
    Array_1d_double_destroy (&Dz_recip);
    Array_1d_double_destroy (&Dy);
    Array_1d_double_destroy (&Dy_N_recip);
    Array_1d_double_destroy (&Dy_S_recip);
    Array_1d_double_destroy (&Dy_recip);
    Array_1d_double_destroy (&Dx);
    Array_1d_double_destroy (&Dx_recip2);

    // In the adjoint code and when diffusion is performed, the pointer
    // diff_input_state_hat points to the output of the above
    // Don't know whether this points to state_interA_hat or state_interB_hat
    // Make sure the output is in state_interA_hat
    if (diff_input_state_hat != &state_interA_hat)
    { // Transfer field into state_interA_hat
      for (z=lower_lev-1; z<=upper_lev+1; z++)
      { for (y=0; y<=nlat+1; y++)
        { for (x=0; x<=nlon+1; x++)
          { state_interA_hat.tracer[x][y][z] = (*diff_input_state_hat).tracer[x][y][z];
          }
        }
      }
    }

    Deallocate_instant (&state_interB_hat, MetaData);

  }
  else
  { // No diffusion, so just copy output from SL scheme
    // Copy result to output
    for (x=1; x<=nlon; x++)
    { for (y=1; y<=nlat; y++)
      { for (z=lower_lev; z<=upper_lev; z++)
        { state_interA_hat.tracer[x][y][z] = (*state_next_hat).tracer[x][y][z];
        }
      }
    }
  }


  // Sort out the halos
  halos_tracer_adj (&state_interA_hat,
                    MetaData);


  // Initialise the output of the adjoint
  //for (z=lower_lev-1; z<=upper_lev+1; z++)
  //{ for (y=0; y<=nlat+1; y++)
  //  { for (x=0; x<=nlon+1; x++)
  //    { (*state_prev_hat).tracer[x][y][z] = 0.0;
  //    }
  //  }
  //}

  if (adv_switch)
  { // Do advection
    if (read_ok)
    { // Loop over arrival points
      for (x=1; x<=nlon; x++)
      { xm1 = x - 1;
        for (y=1; y<=nlat; y++)
        { ym1 = y - 1;
          for (z=lower_lev; z<=upper_lev; z++)
          { zm1 = z - 1;

            particle_lon = dp_lon[xm1][ym1][zm1];
            particle_lat = dp_lat[xm1][ym1][zm1];
            particle_lev = dp_lev[xm1][ym1][zm1];
            index_lon    = dp_index_lon[xm1][ym1][zm1];
            index_lat    = dp_index_lat[xm1][ym1][zm1];
            index_lev    = dp_index_lev[xm1][ym1][zm1];

            // We have the departure points (read-in, above), interpolate from this

            // Interpolate to the departure point - this can be done with either linear or cubic interpolation
            if (interpolate_lc == 'l')
            { // Linear interpolation
              // --------------------

              // Setup short dimensions
              for (d=0; d<2; d++)
              { lon_linear[d] = (*MetaData).longitude[index_lon+d];
                lat_linear[d] = (*MetaData).latitude[index_lat+d];
                lev_linear[d] = (*MetaData).level[index_lev+d];
              }

              // Interpolate
              Interpolate3D_adj ( tracer_linear_hat,                // Out: The 3D field to be interpolated
                                  state_interA_hat.tracer[x][y][z], // In:
                                  lon_linear,                       // The axis x values
                                  lat_linear,                       // The axis y values
                                  lev_linear,                       // The axis z values
                                  particle_lon,                     // The x point to interpolate to
                                  particle_lat,                     // The y point to interpolate to
                                  particle_lev );                   // The z point to interpolate to


              // Prepare cube of data
              Fill_short_array_adj ((*state_prev_hat).tracer, tracer_linear_hat,
                                    index_lon, index_lat, index_lev);


            }
            else
            { // Cubic interpolation
              // -------------------
              // Setup short dimensions
              Fill_short_dims_cubic (nlon, nlat, nlev,
                                     (*MetaData).longitude, (*MetaData).latitude, (*MetaData).level,
                                     lon_cubic, lat_cubic, lev_cubic,
                                     index_lon, index_lat, index_lev);

              // Interpolate
              interpolate_3d_cubic_adj ( lon_cubic, lat_cubic, lev_cubic,  // The axis values
                                         tracer_cubic_hat,                 // out: The 3D field to be interpolated
                                         state_interA_hat.tracer[x][y][z], // in:
                                         particle_lon,                     // The x point to interpolate to
                                         particle_lat,                     // The y point to interpolate to
                                         particle_lev,                     // The x point to interpolate to
                                         3 );                              // Use three separate formula (as Laczos)
              // Prepare cube of data
              Fill_short_array_cubic_adj (nlon, nlat, nlev,
                                          (*state_prev_hat).tracer, tracer_cubic_hat,
                                          index_lon, index_lat, index_lev);

            }

          }
        }
      }
    }
  }
  else
  { // No advection
    for (x=1; x<=nlon; x++)
    { for (y=1; y<=nlat; y++)
      { for (z=lower_lev; z<=upper_lev; z++)
        { (*state_prev_hat).tracer[x][y][z] += state_interA_hat.tracer[x][y][z];
        }
      }
    }
  }


  Deallocate_instant (&state_interA_hat, MetaData);

  // Delete arrays for departure points
  Array_3d_double_destroy (&dp_lon,
                           nlon, nlat);
  Array_3d_double_destroy (&dp_lat,
                           nlon, nlat);
  Array_3d_double_destroy (&dp_lev,
                           nlon, nlat);
  Array_3d_int_destroy (&dp_index_lon,
                        nlon, nlat);
  Array_3d_int_destroy (&dp_index_lat,
                        nlon, nlat);
  Array_3d_int_destroy (&dp_index_lev,
                        nlon, nlat);


  if (!read_ok)
  { printf ("Error from transport_model_adj - the departure-point file specified is not compatible with the dimensions\n");
    exit(0);
  }


}
