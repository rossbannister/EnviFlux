/* ==================================================================================
   3d source and sink code
   Make wind fields at required resolution

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_MakeWinds.out

   Modification history
   --------------------
   22/12/21 New Code. Ross Bannister

   Documentation
   -------------
   Run with
   ./Master_MakeWinds.out <CVT file name (for meta data)> \
                          <filename_containing_full-resolution_ECMWF_winds> \
                          <output_directory_for_reduced_resolution_winds>

   =============================================================================== */

   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )
{ struct HorizTransData_type HorizData;
  struct CVTData_type        CVTdata;
  struct metadata_type       MetaData;
  int                        x, y, z, t, yp1, zp1, xpair, nt, totalwindcount, c;
  struct Wind_type           ECMWF_wind, RedRes_wind;
  int                        *low_index_lon_u, *low_index_lat_u;
  int                        *low_index_lon_v, *low_index_lat_v;
  int                        *low_index_lev_uv, *low_index_lev_w;
  int                        xindex, yindex, zindex;
  double                     xx[2], yy[2], zz[2], vals[2][2][2], pos[3];
  double                     dx, dy, dz;
  char                       CVTfilename[256];
  char                       wind_file_list[256], lod[256], infile[64], outfile[64];
  char                       wind_file_dir[256], complete_infile[256], red_res_dir[256];
  char                       complete_outfile[256];
  char                       *temp;
  bool                       ok;
  FILE*                      wind_file;


  // ===============================================================================
  // Initialisation
  // ===============================================================================

  // Find the name of the file containing the full-resolution ECMWF winds
  ok = MakeWinds_arguments (argument_count,
                            argument_list,
                            CVTfilename,
                            wind_file_list,
                            red_res_dir);

  if (ok)
  { // Initialise and read in generic data needed for the Gaussian latitudes
    printf ("Initialising horizontal transform for Gaussian latitudes\n");
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

    // Allocate space for winds
    Allocate_wind (&ECMWF_wind,
                   &MetaData,
                   true);
    Allocate_wind (&RedRes_wind,
                   &MetaData,
                   false);

    // Allocate arrays to store the lower indices (where in the ECMWF fields the low
    // resolution points are)
    low_index_lon_u  = new int[MetaData.nlon+2];
    low_index_lat_u  = new int[MetaData.nlat+2];
    low_index_lon_v  = new int[MetaData.nlon+2];
    low_index_lat_v  = new int[MetaData.nlat+2];
    low_index_lev_uv = new int[MetaData.nlev+2];
    low_index_lev_w  = new int[MetaData.nlev+2];

    // Set the longitudes, latitudes, and levels of the reduced resolution winds
    for (x=1; x<=MetaData.nlon; x++)
    { RedRes_wind.longitude_u[x] = (MetaData.longitude[x] + MetaData.longitude[x+1]) / 2.0;
      RedRes_wind.longitude_v[x] = MetaData.longitude[x];
    }
    dx = RedRes_wind.longitude_u[2] - RedRes_wind.longitude_u[1];
    RedRes_wind.longitude_u[0]               = RedRes_wind.longitude_u[1] - dx;
    RedRes_wind.longitude_u[MetaData.nlon+1] = RedRes_wind.longitude_u[MetaData.nlon] + dx;
    RedRes_wind.longitude_v[0]               = RedRes_wind.longitude_v[1] - dx;
    RedRes_wind.longitude_v[MetaData.nlon+1] = RedRes_wind.longitude_v[MetaData.nlon] + dx;

    for (y=1; y<=MetaData.nlat; y++)
    { RedRes_wind.latitude_u[y] = MetaData.latitude[y];
      RedRes_wind.latitude_v[y] = (MetaData.latitude[y] + MetaData.latitude[y+1]) / 2.0;
    }
    dy = RedRes_wind.latitude_u[1] - RedRes_wind.latitude_u[2];
    RedRes_wind.latitude_u[0]               = RedRes_wind.latitude_u[1] + dy;
    RedRes_wind.latitude_v[0]               = RedRes_wind.latitude_v[1] + dy;
    dy = RedRes_wind.latitude_u[MetaData.nlat-1] - RedRes_wind.latitude_u[MetaData.nlat];
    RedRes_wind.latitude_u[MetaData.nlat+1] = RedRes_wind.latitude_u[MetaData.nlat] - dy;
    RedRes_wind.latitude_v[MetaData.nlat+1] = RedRes_wind.latitude_v[MetaData.nlat] - dy;

    for (z=1; z<=MetaData.nlev; z++)
    { RedRes_wind.level_uv[z] = MetaData.level[z];
    }
    dz = RedRes_wind.level_uv[1] - RedRes_wind.level_uv[2];
    RedRes_wind.level_uv[0]               = RedRes_wind.level_uv[1] + dz;
    dz = RedRes_wind.level_uv[MetaData.nlev-1] - RedRes_wind.level_uv[MetaData.nlev];
    RedRes_wind.level_uv[MetaData.nlev+1] = RedRes_wind.level_uv[MetaData.nlev] - dz;
    for (z=1; z<=MetaData.nlev; z++)
    { RedRes_wind.level_w[z] = (RedRes_wind.level_uv[z] + RedRes_wind.level_uv[z+1]) / 2.0;
    }
    dz = RedRes_wind.level_w[1] - RedRes_wind.level_w[2];
    RedRes_wind.level_w[0]               = RedRes_wind.level_w[1] + dz;
    dz = RedRes_wind.level_w[MetaData.nlev-1] - RedRes_wind.level_w[MetaData.nlev];
    RedRes_wind.level_w[MetaData.nlev+1] = RedRes_wind.level_w[MetaData.nlev] - dz;


    // ===============================================================================
    // Convert the winds
    // ===============================================================================

    // Deal with the file that contains the wind files
    printf ("Wind file list : %s\n", wind_file_list);
    wind_file = fopen (wind_file_list, "r");
    // Comment lines
    for (c=0; c<6; c++)
    { temp = fgets (lod, 256, wind_file);
    }
    // Directory containing the raw ECMWF winds
    temp = fgets (wind_file_dir, 256, wind_file);
    printf ("Directory containing the wind files : %s\n", wind_file_dir);

    totalwindcount = 0;

    while (!feof(wind_file))
    { // Read-in and construct the full filename of the next ECMWF wind file
      fscanf (wind_file, "%s %i\n", infile, &nt);
      //printf ("%s, %i\n", infile, nt);
      construct_total_filename (wind_file_dir,    // in
                                infile,           // in
                                complete_infile); // out

      printf ("ECMWF wind file %s\n\n", complete_infile);

      for (t=0; t<nt; t++)
      { // Read-in ECMWF wind field

        Read_ecmwf_winds (&ECMWF_wind,
                          t,
                          complete_infile);

        if (totalwindcount == 0)
        { // Find the lower indices for each lower-resolution grid point
          // printf ("Check lower index values for longitude_u\n");
          for (x=1; x<=RedRes_wind.nlon; x++)
          { low_index_lon_u[x] = Find_index_ascend (::ecmwf_nlon+2,
                                                    ECMWF_wind.longitude_u,
                                                    RedRes_wind.longitude_u[x]);
            //if (low_index_lon_u[x] > -1)
            //{ printf ("  %f (%i) is between %f and %f\n", RedRes_wind.longitude_u[x],
            //          low_index_lon_u[x],
            //          ECMWF_wind.longitude_u[low_index_lon_u[x]],
            //          ECMWF_wind.longitude_u[low_index_lon_u[x]+1]);
            //}
            if (low_index_lon_u[x] <= -1)
            { printf ("  %f is out of range\n", RedRes_wind.longitude_u[x]);
            }
          }

          //printf ("Check lower index values for longitude_v\n");
          for (x=1; x<=RedRes_wind.nlon; x++)
          { low_index_lon_v[x] = Find_index_ascend (::ecmwf_nlon+2,
                                                    ECMWF_wind.longitude_v,
                                                    RedRes_wind.longitude_v[x]);
            //if (low_index_lon_v[x] > -1)
            //{ printf ("  %f (%i) is between %f and %f\n", RedRes_wind.longitude_v[x],
            //          low_index_lon_v[x],
            //          ECMWF_wind.longitude_v[low_index_lon_v[x]],
            //          ECMWF_wind.longitude_v[low_index_lon_v[x]+1]);
            //}
            if (low_index_lon_v[x] <= -1)
            { printf ("  %f is out of range\n", RedRes_wind.longitude_v[x]);
            }
          }

          //printf ("Check lower index values for latitude_u\n");
          for (y=1; y<=RedRes_wind.nlat; y++)
          { low_index_lat_u[y] = Find_index_descend (::ecmwf_nlat+2,
                                                   ECMWF_wind.latitude_u,
                                                   RedRes_wind.latitude_u[y]);
            //if (low_index_lat_u[y] > -1)
            //{ printf ("  %f (%i) is between %f and %f\n", RedRes_wind.latitude_u[y],
            //          low_index_lat_u[y],
            //          ECMWF_wind.latitude_u[low_index_lat_u[y]],
            //          ECMWF_wind.latitude_u[low_index_lat_u[y]+1]);
            //}
            if (low_index_lat_u[y] <= -1)
            { printf ("  %f is out of range\n", RedRes_wind.latitude_u[y]);
            }
          }

          //printf ("Check lower index values for latitude_v\n");
          for (y=1; y<=RedRes_wind.nlat; y++)
          { low_index_lat_v[y] = Find_index_descend (::ecmwf_nlat+2,
                                                   ECMWF_wind.latitude_v,
                                                   RedRes_wind.latitude_v[y]);
            //if (low_index_lat_v[y] > -1)
            //{ printf ("  %f (%i) is between %f and %f\n", RedRes_wind.latitude_v[y],
            //          low_index_lat_v[y],
            //          ECMWF_wind.latitude_v[low_index_lat_v[y]],
            //          ECMWF_wind.latitude_v[low_index_lat_v[y]+1]);
            //}
            if (low_index_lat_v[y] <= -1)
            { printf ("  %f is out of range\n", RedRes_wind.latitude_v[y]);
            }
          }

          //printf ("Check lower index values for level_uv\n");
          for (z=1; z<=RedRes_wind.nlev; z++)
          { low_index_lev_uv[z] = Find_index_descend (::ecmwf_nlev+2,
                                                      ECMWF_wind.level_uv,
                                                      RedRes_wind.level_uv[z]);
            //if (low_index_lev_uv[z] > -1)
            //{ printf ("  %f (%i) is between %f and %f\n", RedRes_wind.level_uv[z],
            //          low_index_lev_uv[z],
            //          ECMWF_wind.level_uv[low_index_lev_uv[z]],
            //          ECMWF_wind.level_uv[low_index_lev_uv[z]+1]);
            //}
            if (low_index_lev_uv[z] <= -1)
            { printf ("  %f is out of range\n", RedRes_wind.level_uv[z]);
            }
          }

          //printf ("Check lower index values for level_w\n");
          for (z=1; z<=RedRes_wind.nlev; z++)
          { low_index_lev_w[z] = Find_index_descend (::ecmwf_nlev+2,
                                                     ECMWF_wind.level_w,
                                                     RedRes_wind.level_w[z]);
            //if (low_index_lev_w[z] > -1)
            //{ printf ("  %f (%i) is between %f and %f\n", RedRes_wind.level_w[z],
            //          low_index_lev_w[z],
            //          ECMWF_wind.level_w[low_index_lev_w[z]],
            //          ECMWF_wind.level_w[low_index_lev_w[z]+1]);
            //}
            if (low_index_lev_w[z] <= -1)
            { printf ("  %f is out of range\n", RedRes_wind.level_w[z]);
            }
          }
        }


        // Do the interpolation for the zonal wind
        for (x=1; x<=RedRes_wind.nlon; x++)
        { // Set the longitude values of the ECMWF grid about this point
          xindex = low_index_lon_u[x];
          xx[0]  = ECMWF_wind.longitude_u[xindex];
          xx[1]  = ECMWF_wind.longitude_u[xindex+1];
          // Set the x position of the reduced resolution point
          pos[0] = RedRes_wind.longitude_u[x];

          for (y=1; y<=RedRes_wind.nlat; y++)
          { // Set the latitude values of the ECMWF grid about this point
            yindex = low_index_lat_u[y];
            yy[0]  = ECMWF_wind.latitude_u[yindex];
            yy[1]  = ECMWF_wind.latitude_u[yindex+1];
            // Set the y position of the reduced resolution point
            pos[1] = RedRes_wind.latitude_u[y];

            for (z=1; z<=RedRes_wind.nlev; z++)
            { // Set the level values of the ECMWF grid about this point
              zindex = low_index_lev_uv[z];
              zz[0]  = ECMWF_wind.level_uv[zindex];
              zz[1]  = ECMWF_wind.level_uv[zindex+1];
              // Set the z position of the reduced resolution point
              pos[2] = RedRes_wind.level_uv[z];

              // Set all 8 values in the cube around this point for the u wind
              vals[0][0][0] = ECMWF_wind.u[xindex]  [yindex]  [zindex];
              vals[1][0][0] = ECMWF_wind.u[xindex+1][yindex]  [zindex];
              vals[0][1][0] = ECMWF_wind.u[xindex]  [yindex+1][zindex];
              vals[1][1][0] = ECMWF_wind.u[xindex+1][yindex+1][zindex];
              vals[0][0][1] = ECMWF_wind.u[xindex]  [yindex]  [zindex+1];
              vals[1][0][1] = ECMWF_wind.u[xindex+1][yindex]  [zindex+1];
              vals[0][1][1] = ECMWF_wind.u[xindex]  [yindex+1][zindex+1];
              vals[1][1][1] = ECMWF_wind.u[xindex+1][yindex+1][zindex+1];

              // Interpolate
              RedRes_wind.u[x][y][z] = interpolate_3d (xx, yy, zz,
                                                       vals, pos);
            }
          }
        }

        // Do the interpolation for the meridional wind
        for (x=1; x<=RedRes_wind.nlon; x++)
        { // Set the longitude values of the ECMWF grid about this point
          xindex = low_index_lon_v[x];
          xx[0]  = ECMWF_wind.longitude_v[xindex];
          xx[1]  = ECMWF_wind.longitude_v[xindex+1];
          // Set the x position of the reduced resolution point
          pos[0] = RedRes_wind.longitude_v[x];

          for (y=1; y<=RedRes_wind.nlat; y++)
          { // Set the latitude values of the ECMWF grid about this point
            yindex = low_index_lat_v[y];
            yy[0]  = ECMWF_wind.latitude_v[yindex];
            yy[1]  = ECMWF_wind.latitude_v[yindex+1];
            // Set the y position of the reduced resolution point
            pos[1] = RedRes_wind.latitude_v[y];

            for (z=1; z<=RedRes_wind.nlev; z++)
            { // Set the level values of the ECMWF grid about this point
              zindex = low_index_lev_uv[z];
              zz[0]  = ECMWF_wind.level_uv[zindex];
              zz[1]  = ECMWF_wind.level_uv[zindex+1];
              // Set the z position of the reduced resolution point
              pos[2] = RedRes_wind.level_uv[z];

              // Set all 8 values in the cube around this point for the u wind
              vals[0][0][0] = ECMWF_wind.v[xindex]  [yindex]  [zindex];
              vals[1][0][0] = ECMWF_wind.v[xindex+1][yindex]  [zindex];
              vals[0][1][0] = ECMWF_wind.v[xindex]  [yindex+1][zindex];
              vals[1][1][0] = ECMWF_wind.v[xindex+1][yindex+1][zindex];
              vals[0][0][1] = ECMWF_wind.v[xindex]  [yindex]  [zindex+1];
              vals[1][0][1] = ECMWF_wind.v[xindex+1][yindex]  [zindex+1];
              vals[0][1][1] = ECMWF_wind.v[xindex]  [yindex+1][zindex+1];
              vals[1][1][1] = ECMWF_wind.v[xindex+1][yindex+1][zindex+1];

              // Interpolate
              RedRes_wind.v[x][y][z] = interpolate_3d (xx, yy, zz,
                                                       vals, pos);
            }
          }
        }

        // Do the interpolation for the vertical wind
        for (x=1; x<=RedRes_wind.nlon; x++)
        { // Set the longitude values of the ECMWF grid about this point
          xindex = low_index_lon_u[x];
          xx[0]  = ECMWF_wind.longitude_v[xindex];
          xx[1]  = ECMWF_wind.longitude_v[xindex+1];
          // Set the x position of the reduced resolution point
          pos[0] = RedRes_wind.longitude_v[x];

          for (y=1; y<=RedRes_wind.nlat; y++)
          { // Set the latitude values of the ECMWF grid about this point
            yindex = low_index_lat_u[y];
            yy[0]  = ECMWF_wind.latitude_u[yindex];
            yy[1]  = ECMWF_wind.latitude_u[yindex+1];
            // Set the y position of the reduced resolution point
            pos[1] = RedRes_wind.latitude_u[y];

            for (z=1; z<=RedRes_wind.nlev; z++)
            { // Set the level values of the ECMWF grid about this point
              zindex = low_index_lev_w[z];
              zz[0]  = ECMWF_wind.level_w[zindex];
              zz[1]  = ECMWF_wind.level_w[zindex+1];
              // Set the z position of the reduced resolution point
              pos[2] = RedRes_wind.level_w[z];

              // Set all 8 values in the cube around this point for the u wind
              vals[0][0][0] = ECMWF_wind.w[xindex]  [yindex]  [zindex];
              vals[1][0][0] = ECMWF_wind.w[xindex+1][yindex]  [zindex];
              vals[0][1][0] = ECMWF_wind.w[xindex]  [yindex+1][zindex];
              vals[1][1][0] = ECMWF_wind.w[xindex+1][yindex+1][zindex];
              vals[0][0][1] = ECMWF_wind.w[xindex]  [yindex]  [zindex+1];
              vals[1][0][1] = ECMWF_wind.w[xindex+1][yindex]  [zindex+1];
              vals[0][1][1] = ECMWF_wind.w[xindex]  [yindex+1][zindex+1];
              vals[1][1][1] = ECMWF_wind.w[xindex+1][yindex+1][zindex+1];

              // Interpolate
              RedRes_wind.w[x][y][z] = interpolate_3d (xx, yy, zz,
                                                       vals, pos);
            }
          }
        }

        //printf ("u values before halo set\n");
        //for (x=0; x<RedRes_wind.nlon+2; x++)
        //{ printf ("%i  %f  %f\n", x, RedRes_wind.longitude_u[x], RedRes_wind.u[x][10][10]);
        //}


        // ===============================================================================
        // Sort out the halo points for these wind components
        // ===============================================================================

        // (a) in the longitude direction
        for (y=0; y<RedRes_wind.nlat; y++)
        { yp1 = y + 1;
          for (z=0; z<RedRes_wind.nlev; z++)
          { zp1 = z + 1;
            RedRes_wind.u[0][yp1][zp1]                  = RedRes_wind.u[RedRes_wind.nlon][yp1][zp1];
            RedRes_wind.v[0][yp1][zp1]                  = RedRes_wind.v[RedRes_wind.nlon][yp1][zp1];
            RedRes_wind.w[0][yp1][zp1]                  = RedRes_wind.w[RedRes_wind.nlon][yp1][zp1];
            RedRes_wind.u[RedRes_wind.nlon+1][yp1][zp1] = RedRes_wind.u[1][yp1][zp1];
            RedRes_wind.v[RedRes_wind.nlon+1][yp1][zp1] = RedRes_wind.v[1][yp1][zp1];
            RedRes_wind.w[RedRes_wind.nlon+1][yp1][zp1] = RedRes_wind.w[1][yp1][zp1];
          }
        }

        // (b) in the latitude direction
        for(x=0; x<RedRes_wind.nlon+2; x++)
        { xpair = x + RedRes_wind.nlon/2 + 1;
          if (xpair > RedRes_wind.nlon)
          { xpair -= RedRes_wind.nlon;
          }
          for (z=0; z<RedRes_wind.nlev; z++)
          { zp1 = z + 1;
            RedRes_wind.u[x][0][zp1]                  = -1.0 * RedRes_wind.u[xpair][2][zp1];
            RedRes_wind.v[x][0][zp1]                  = -1.0 * RedRes_wind.v[xpair][2][zp1];
            RedRes_wind.w[x][0][zp1]                  = RedRes_wind.w[xpair][2][zp1];
            RedRes_wind.u[x][RedRes_wind.nlat+1][zp1] = -1.0 * RedRes_wind.u[xpair][RedRes_wind.nlat-1][zp1];
            RedRes_wind.v[x][RedRes_wind.nlat+1][zp1] = -1.0 * RedRes_wind.v[xpair][RedRes_wind.nlat-1][zp1];
            RedRes_wind.w[x][RedRes_wind.nlat+1][zp1] = RedRes_wind.w[xpair][RedRes_wind.nlat-1][zp1];
          }
        }

        // (c) in the height direction
        // Use Neuman boundary conditions for u and v, and Dirichlet for w
        for(x=0; x<RedRes_wind.nlon+2; x++)
        { for (y=0; y<RedRes_wind.nlat+2; y++)
          { // Top of model
            RedRes_wind.u[x][y][0] = RedRes_wind.u[x][y][1];
            RedRes_wind.v[x][y][0] = RedRes_wind.v[x][y][1];
            RedRes_wind.w[x][y][0] = -1.0 * RedRes_wind.w[x][y][1];
            // Bottom of model
            RedRes_wind.u[x][y][RedRes_wind.nlev+1] = RedRes_wind.u[x][y][RedRes_wind.nlev];
            RedRes_wind.v[x][y][RedRes_wind.nlev+1] = RedRes_wind.v[x][y][RedRes_wind.nlev];
            RedRes_wind.w[x][y][RedRes_wind.nlev+1] = -1.0 * RedRes_wind.w[x][y][RedRes_wind.nlev];
          }
        }

        //printf ("u values after halo set\n");
        //for (x=0; x<RedRes_wind.nlon+2; x++)
        //{ printf ("%i  %f  %f\n", x, RedRes_wind.longitude_u[x], RedRes_wind.u[x][10][10]);
        //}


        // ===============================================================================
        // Write out the interpolated winds
        // ===============================================================================

        sprintf (outfile, "Winds%04i.nc", totalwindcount);

        construct_total_filename (red_res_dir,       // in
                                  outfile,           // in
                                  complete_outfile); // out

        printf ("Writing %s\n", complete_outfile);

        Write_winds (&RedRes_wind,
                     complete_outfile);

        totalwindcount++;

      }
    }

    // ===============================================================================
    // Deallocate
    // ===============================================================================

    DeallocateHorizTrans (&HorizData);
    Deallocate_CVT (&CVTdata);
    Deallocate_metadata (&MetaData);
    Deallocate_wind (&ECMWF_wind);
    Deallocate_wind (&RedRes_wind);
    delete[] low_index_lon_u;
    delete[] low_index_lat_u;
    delete[] low_index_lon_v;
    delete[] low_index_lat_v;
    delete[] low_index_lev_uv;
    delete[] low_index_lev_w;

  }
}
