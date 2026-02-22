   #include <source.h>


// -------------------------------------------------------------------------------
void trajectory ( int                   minor_t,           // Minor timestep No., local time (relative to windA)
                  double                dt,                // Minor timestep size
                  double                Dt,                // Major timestep size
                  struct Wind_type      *windA,            // Wind fields at start of major timestep
                  struct Wind_type      *windB,            // Wind fields at end of major timestep
                  double                lon_start,         // particle longitude at start time
                  double                lat_start,         // particle latitude at start time
                  double                lev_start,         // particle level at start time
                  double                *lon_end,          // particle longitude at end time
                  double                *lat_end,          // particle latitude at end time
                  double                *lev_end,          // particle level at end time
                  bool                  show_diags )

{ // Compute a forward trajectory element based on fourth order Runge-Kutta
  double k1_lon, k1_lat, k1_lev;
  double k2_lon, k2_lat, k2_lev;
  double k3_lon, k3_lat, k3_lev;
  double k4_lon, k4_lat, k4_lev;
  double lon_new, lat_new, lev_new, t_new;
  double times[2], t, top_z;

  // Set the time axis
  times[0] = 0.0;
  times[1] = Dt;
  t        = double(minor_t) * fabs(dt);

  if (show_diags)
  { printf ("\n\nMinor timestep %i, time %f   : limits %f, %f\n", minor_t, t, times[0], times[1]);
    printf ("Starting position %f %f %f  %f\n", lon_start, lat_start, lev_start, t);
  }

  // Retrieve the top height of the model
  top_z = (*windA).level_uv[(*windA).nlev];


  // First term in the Runge-Kutta
  // -----------------------------
  //printf ("First term in RK\n");
  RungeKutta_k (dt,                // Minor timestep size
                times,             // Time axis (relative units)
                windA,             // Wind fields at start of major timestep
                windB,             // Wind fields at end of major timestep
                lon_start,         // particle longitude
                lat_start,         // particle latitude
                lev_start,         // particle level
                t,                 // local time (relative to first wind field)
                &k1_lon,           // k value
                &k1_lat,           // k value
                &k1_lev,           // k value
                show_diags);

  // Form new position
  t_new   = t + dt / 2.0;
  lon_new = lon_start + k1_lon / 2.0;
  lat_new = lat_start + k1_lat / 2.0;
  lev_new = lev_start + k1_lev / 2.0;
  if (show_diags)
  { printf ("New position 1:  %f %f %f  %f\n", lon_new, lat_new, lev_new, t_new);
  }
  // Check bounds
  CheckBounds (&lon_new, &lat_new, &lev_new, top_z);


  // Second term in the Runge-Kutta
  // ------------------------------
  //printf ("Second term in RK\n");
  RungeKutta_k (dt,                // Minor timestep size
                times,             // Time axis (relative units)
                windA,             // Wind fields at start of major timestep
                windB,             // Wind fields at end of major timestep
                lon_new,           // particle longitude
                lat_new,           // particle latitude
                lev_new,           // particle level
                t_new,             // local time (relative to first wind field)
                &k2_lon,           // k value
                &k2_lat,           // k value
                &k2_lev,           // k value
                show_diags);

  // Form new position
  lon_new = lon_start + k2_lon / 2.0;
  lat_new = lat_start + k2_lat / 2.0;
  lev_new = lev_start + k2_lev / 2.0;
  if (show_diags)
  { printf ("New position 2:  %f %f %f  %f\n", lon_new, lat_new, lev_new, t_new);
  }
  // Check bounds
  CheckBounds (&lon_new, &lat_new, &lev_new, top_z);


  // Third term in the Runge-Kutta
  // -----------------------------
  //printf ("Third term in RK\n");
  RungeKutta_k (dt,                // Minor timestep size
                times,             // Time axis (relative units)
                windA,             // Wind fields at start of major timestep
                windB,             // Wind fields at end of major timestep
                lon_new,           // particle longitude
                lat_new,           // particle latitude
                lev_new,           // particle level
                t_new,             // local time (relative to first wind field)
                &k3_lon,           // k value
                &k3_lat,           // k value
                &k3_lev,           // k value
                show_diags);

  // Form new position
  t_new   = t + dt;
  lon_new = lon_start + k3_lon;
  lat_new = lat_start + k3_lat;
  lev_new = lev_start + k3_lev;
  if (show_diags)
  { printf ("New position 3:  %f %f %f  %f\n", lon_new, lat_new, lev_new, t_new);
  }
  // Check bounds
  CheckBounds (&lon_new, &lat_new, &lev_new, top_z);


  // Fourth term in the Runge-Kutta
  // ------------------------------
  //printf ("Fourth term in RK\n");
  RungeKutta_k (dt,                // Minor timestep size
                times,             // Time axis (relative units)
                windA,             // Wind fields at start of major timestep
                windB,             // Wind fields at end of major timestep
                lon_new,           // particle longitude
                lat_new,           // particle latitude
                lev_new,           // particle level
                t_new,             // local time (relative to first wind field)
                &k4_lon,           // k value
                &k4_lat,           // k value
                &k4_lev,           // k value
                show_diags);

  // Form new (final) position
  *lon_end = lon_start + k1_lon / 6.0 + k2_lon / 3.0 + k3_lon / 3.0 + k4_lon / 6.0;
  *lat_end = lat_start + k1_lat / 6.0 + k2_lat / 3.0 + k3_lat / 3.0 + k4_lat / 6.0;
  *lev_end = lev_start + k1_lev / 6.0 + k2_lev / 3.0 + k3_lev / 3.0 + k4_lev / 6.0;
  if (show_diags)
  { printf ("New position 4:  %f %f %f  %f\n", *lon_end, *lat_end, *lev_end, t_new);
  }
  // Check bounds
  CheckBounds (lon_end, lat_end, lev_end, top_z);
}




// -------------------------------------------------------------------------------
void RungeKutta_k ( double                dt,                // Minor timestep size
                    double                times[2],          // Time axis (relative units)
                    struct Wind_type      *windA,            // Wind fields at start of major timestep
                    struct Wind_type      *windB,            // Wind fields at end of major timestep
                    double                lon,               // particle longitude
                    double                lat,               // particle latitude
                    double                lev,               // particle level
                    double                t,                 // time
                    double                *kx,               // k value
                    double                *ky,               // k value
                    double                *kz,               // k value
                    bool                  show_diags)

{ // Compute a k-value needed for the Runge-Kutta method
  int    index_lon_u, index_lat_u, index_lev_uv;
  int    index_lon_v, index_lat_v, index_lev_w;
  double u, v, w;
  double uA[2][2][2], uB[2][2][2];
  double vA[2][2][2], vB[2][2][2];
  double wA[2][2][2], wB[2][2][2];
  double lon_u[2], lon_v[2], lat_u[2], lat_v[2], lev_uv[2], lev_w[2];

  // Find the indices of the longitude, latitude, and level
  index_lon_u  = Find_index_ascend ((*windA).nlon+2,
                                    (*windA).longitude_u,
                                    lon);
  index_lon_v  = Find_index_ascend ((*windA).nlon+2,
                                    (*windA).longitude_v,
                                    lon);
  index_lat_u  = Find_index_descend ((*windA).nlat+2,
                                     (*windA).latitude_u,
                                     lat);
  index_lat_v  = Find_index_descend ((*windA).nlat+2,
                                     (*windA).latitude_v,
                                     lat);
  index_lev_uv = Find_index_ascend ((*windA).nlev+2,
                                    (*windA).level_uv,
                                    lev);
  index_lev_w  = Find_index_ascend ((*windA).nlev+2,
                                    (*windA).level_w,
                                    lev);

  Fill_short_dims ((*windA).longitude_u, (*windA).longitude_v,
                   (*windA).latitude_u, (*windA).latitude_v,
                   (*windA).level_uv, (*windA).level_w,
                   lon_u, lon_v, lat_u, lat_v, lev_uv, lev_w,
                   index_lon_u, index_lon_v, index_lat_u, index_lat_v, index_lev_uv, index_lev_w);

  // Prepare cube of data
  Fill_short_array ((*windA).u, uA,
                    index_lon_u, index_lat_u, index_lev_uv);
  Fill_short_array ((*windB).u, uB,
                    index_lon_u, index_lat_u, index_lev_uv);

  // Interpolate
  u = Interpolate3Dt ( uA,          // The 3D field to be interpolated at time 0
                       uB,          // The 3D field to be interpolated at time 1
                       lon_u,       // The axis x values
                       lat_u,       // The axis y values
                       lev_uv,      // The axis z values
                       times,       // The axis time values
                       lon,         // The x point to interpolate to
                       lat,         // The y point to interpolate to
                       lev,         // The z point to interpolate to
                       t);          // The time point to interpolate to

  // Prepare cube of data
  Fill_short_array ((*windA).v, vA,
                    index_lon_v, index_lat_v, index_lev_uv);
  Fill_short_array ((*windB).v, vB,
                    index_lon_v, index_lat_v, index_lev_uv);

  // Interpolate
  v = Interpolate3Dt ( vA,          // The 3D field to be interpolated at time 0
                       vB,          // The 3D field to be interpolated at time 1
                       lon_v,       // The axis x values
                       lat_v,       // The axis y values
                       lev_uv,      // The axis z values
                       times,       // The axis time values
                       lon,         // The x point to interpolate to
                       lat,         // The y point to interpolate to
                       lev,         // The z point to interpolate to
                       t);          // The time point to interpolate to

  // Prepare cube of data
  Fill_short_array ((*windA).w, wA,
                    index_lon_v, index_lat_u, index_lev_w);
  Fill_short_array ((*windB).w, wB,
                    index_lon_v, index_lat_u, index_lev_w);

  // Interpolate
  w = Interpolate3Dt ( wA,          // The 3D field to be interpolated at time 0
                       wB,          // The 3D field to be interpolated at time 1
                       lon_v,       // The axis x values
                       lat_u,       // The axis y values
                       lev_w,       // The axis z values
                       times,       // The axis time values
                       lon,         // The x point to interpolate to
                       lat,         // The y point to interpolate to
                       lev,         // The z point to interpolate to
                       t);          // The time point to interpolate to

  // Compute k value components
  *kx = dt * u * ::m2deg / cos(lat * ::deg2rad);
  *ky = dt * v * ::m2deg;
  *kz = dt * w;

  if (show_diags)
  { printf ("  Runge_kutta  ix iy iz :    %i %i    %i %i    %i %i\n", index_lon_u, index_lon_v, index_lat_u, index_lat_v, index_lev_uv, index_lev_w);
    printf ("               u  v  w  :    %f %f %f\n", u, v, w);
    printf ("               kx ky kz :    %f %f %f\n", *kx, *ky, *kz);
  }

}



// -------------------------------------------------------------------------------
void Fill_short_array ( double ***large,
                        double small[2][2][2],
                        int    index_x,
                        int    index_y,
                        int    index_z )
{ // Fill a small 2,2,2 array from a cube of a larger array (used for linear interpolation)
  small[0][0][0] = large[index_x]  [index_y]  [index_z];
  small[0][0][1] = large[index_x]  [index_y]  [index_z+1];
  small[0][1][0] = large[index_x]  [index_y+1][index_z];
  small[0][1][1] = large[index_x]  [index_y+1][index_z+1];
  small[1][0][0] = large[index_x+1][index_y]  [index_z];
  small[1][0][1] = large[index_x+1][index_y]  [index_z+1];
  small[1][1][0] = large[index_x+1][index_y+1][index_z];
  small[1][1][1] = large[index_x+1][index_y+1][index_z+1];
}


// -------------------------------------------------------------------------------
void Fill_short_array_adj ( double ***large,
                            double small[2][2][2],
                            int    index_x,
                            int    index_y,
                            int    index_z )
{ // Adjoint of: fill a small 2,2,2 array from a cube of a larger array (used for linear interpolation)
  large[index_x]  [index_y]  [index_z]   += small[0][0][0];
  large[index_x]  [index_y]  [index_z+1] += small[0][0][1];
  large[index_x]  [index_y+1][index_z]   += small[0][1][0];
  large[index_x]  [index_y+1][index_z+1] += small[0][1][1];
  large[index_x+1][index_y]  [index_z]   += small[1][0][0];
  large[index_x+1][index_y]  [index_z+1] += small[1][0][1];
  large[index_x+1][index_y+1][index_z]   += small[1][1][0];
  large[index_x+1][index_y+1][index_z+1] += small[1][1][1];
}


// -------------------------------------------------------------------------------
void Fill_short_dims (double *lon_u_large,
                      double *lon_v_large,
                      double *lat_u_large,
                      double *lat_v_large,
                      double *lev_uv_large,
                      double *lev_w_large,
                      double lon_u[2],
                      double lon_v[2],
                      double lat_u[2],
                      double lat_v[2],
                      double lev_uv[2],
                      double lev_w[2],
                      int    index_ux,
                      int    index_vx,
                      int    index_uy,
                      int    index_vy,
                      int    index_uvz,
                      int    index_wz )
{ // Fill a small 2 array from axes (used for linear interpolation)

//if (index_ux<2 || index_vx<2 || index_uy<2 || index_vy<2 || index_uvz<2 || index_wz<2)
//{ printf ("%i  %i  %i  %i  %i  %i\n", index_ux, index_vx, index_uy, index_vy, index_uvz, index_wz);
//}

//if (index_vy == -1)
//{ printf ("%f %f\n", lat_u[0], lat_u[1]);
//  printf ("%f %f\n", lat_v[0], lat_v[1]);
//}

  lon_u[0]  = lon_u_large[index_ux];
  lon_u[1]  = lon_u_large[index_ux+1];
  lon_v[0]  = lon_v_large[index_vx];
  lon_v[1]  = lon_v_large[index_vx+1];
  lat_u[0]  = lat_u_large[index_uy];
  lat_u[1]  = lat_u_large[index_uy+1];
  lat_v[0]  = lat_v_large[index_vy];
  lat_v[1]  = lat_v_large[index_vy+1];
  lev_uv[0] = lev_uv_large[index_uvz];
  lev_uv[1] = lev_uv_large[index_uvz+1];
  lev_w[0]  = lev_w_large[index_wz];
  lev_w[1]  = lev_w_large[index_wz+1];
}


// -------------------------------------------------------------------------------
void Fill_short_dims_cubic (int    nlon,
                            int    nlat,
                            int    nlev,
                            double *lon_large,
                            double *lat_large,
                            double *lev_large,
                            double lon_cubic[6],
                            double lat_cubic[6],
                            double lev_cubic[6],
                            int    index_lon,
                            int    index_lat,
                            int    index_lev)
{ // Fill a small 6 array from axes (used for cubic interpolation)
  int d, data_from;

  for (d=0; d<6; d++)
  { // Set-up the longitude axis values for the short dimensions
    data_from = index_lon + d - 2;
    if (data_from < 0)
    { lon_cubic[d] = lon_large[data_from+nlon] - 360.0;
    }
    else
    { if (data_from > nlon+1)
      { lon_cubic[d] = lon_large[data_from-nlon] + 360.0;
      }
      else
      { lon_cubic[d] = lon_large[data_from];
      }
    }

    // Set-up the latitude axis values for the short dimensions
    data_from = index_lat + d - 2;
    if (data_from < 0)
    { lat_cubic[d] = 2.0 * lat_large[1] - lat_large[2-data_from];
    }
    else
    { if (data_from > nlat+1)
      { lat_cubic[d] = 2.0 * lat_large[nlat] - lat_large[2*nlat-data_from];
      }
      else
      { lat_cubic[d] = lat_large[data_from];
      }
    }

    // Set-up the latitude axis values for the short dimensions
    data_from = index_lev + d - 2;
    if (data_from < 0)
    { lev_cubic[d] = 2.0 * lev_large[1] - lev_large[2-data_from];
    }
    else
    { if (data_from > nlev+1)
      { lev_cubic[d] = 2.0 * lev_large[nlev] - lev_large[2*nlev-data_from];
      }
      else
      { lev_cubic[d] = lev_large[data_from];
      }
    }
  }
}

// -------------------------------------------------------------------------------
void Fill_short_array_cubic ( int    nlon,
                              int    nlat,
                              int    nlev,
                              double ***large,        // in
                              double small[6][6][6],  // out
                              int    index_lon,
                              int    index_lat,
                              int    index_lev )
{ // Fill a small 6,6,6 array from a cube of a larger array (used for linear interpolation)
  int dx, dy, dz;
  int data_from_x, data_from_y, data_from_z;

  for (dx=0; dx<6; dx++)
  { data_from_x = index_lon + dx - 2;
    if (data_from_x < 0)
    { data_from_x += nlon;
    }
    else
    { if (data_from_x > nlon+1)
      { data_from_x -= nlon;
      }
    }

    for (dy=0; dy<6; dy++)
    { data_from_y = index_lat + dy - 2;
      if (data_from_y < 0)
      { data_from_y  = 1 - data_from_y;
        data_from_x += nlon;
        if (data_from_x > nlon+1)
        { data_from_x -= nlon;
        }
      }
      else
      { if (data_from_y > nlat+1)
        { data_from_y  = 2*nlat - data_from_y + 1;
          data_from_x += nlon;
          if (data_from_x > nlon+1)
          { data_from_x -= nlon;
          }
        }
      }

      for (dz=0; dz<6; dz++)
      { data_from_z = index_lev + dz - 2;
        if (data_from_z < 0)
        { data_from_z = 0;
        }
        else
        { if (data_from_z > nlev+1)
          { data_from_z = nlev + 1;
          }
        }

        small[dx][dy][dz] = large[data_from_x][data_from_y][data_from_z];
      }
    }
  }
}


// -------------------------------------------------------------------------------
void Fill_short_array_cubic_adj ( int    nlon,
                                  int    nlat,
                                  int    nlev,
                                  double ***large,       // out
                                  double small[6][6][6], // in
                                  int    index_lon,
                                  int    index_lat,
                                  int    index_lev )
{ // Fill a small 6,6,6 array from a cube of a larger array (used for linear interpolation)
  int dx, dy, dz;
  int data_from_x, data_from_y, data_from_z;

  for (dx=0; dx<6; dx++)
  { data_from_x = index_lon + dx - 2;
    if (data_from_x < 0)
    { data_from_x += nlon;
    }
    else
    { if (data_from_x > nlon+1)
      { data_from_x -= nlon;
      }
    }

    for (dy=0; dy<6; dy++)
    { data_from_y = index_lat + dy - 2;
      if (data_from_y < 0)
      { data_from_y  = 1 - data_from_y;
        data_from_x += nlon;
        if (data_from_x > nlon+1)
        { data_from_x -= nlon;
        }
      }
      else
      { if (data_from_y > nlat+1)
        { data_from_y  = 2*nlat - data_from_y + 1;
          data_from_x += nlon;
          if (data_from_x > nlon+1)
          { data_from_x -= nlon;
          }
        }
      }

      for (dz=0; dz<6; dz++)
      { data_from_z = index_lev + dz - 2;
        if (data_from_z < 0)
        { data_from_z = 0;
        }
        else
        { if (data_from_z > nlev+1)
          { data_from_z = nlev + 1;
          }
        }

        large[data_from_x][data_from_y][data_from_z] += small[dx][dy][dz];
      }
    }
  }
}


// -------------------------------------------------------------------------------
void CheckBounds ( double *lon,
                   double *lat,
                   double *lev,
                   double top )
{ // Check that the lon, lat and lev are all in bounds, adjust otherwise
  double latdiff, lon1;

  // Check if latitude has 'gone over the poles'
  if (*lat > 90.0)
  { latdiff = *lat - 90.0;
    *lat    = 90.0 - latdiff;
    *lon   += 180.0;
  }
  else
  { if (*lat < -90.0)
    { latdiff = -90.0 - *lat;
      *lat    = -90.0 + latdiff;
      *lon   += 180.0;
    }
  }

  // Check if longitude has 'gone over the Greenwich meridian'
  lon1 = *lon;
  if (lon1 < 0.0)
  { lon1 -= int(lon1/360.0) * 360.0 - 360.0;
  }
  if (lon1 >= 360.0)
  { lon1 -= int(lon1/360.0) * 360.0;
  }
  *lon = lon1;

  // Check if height has gone outside of the atmosphere
  if (*lev < 0.0)
  { *lev = 0.0;
  }
  else
  { if (*lev > top)
    { *lev = top;
    }
  }
}
