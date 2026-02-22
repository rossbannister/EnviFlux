   #include <source.h>


// -------------------------------------------------------------------------------
void WriteSingleRealSpField
// Function to write-out a single real-space field
 ( int    nlon,
   int    nlat,
   double **field,
   double *lons,
   double *lats,
   char   filename[256],
   char   varname[32] )

{ //Declare local variables
  int          x, y;
  // netCDF-related variables
  size_t       nlon1, nlat1;
  double       lineofdata[nlon];
  int          ierr[9], ncid, err;
  int          dimidx, dimidy;
  int          varidx, varidy, varid;
  int          list1[1], list2[2];
  size_t       start[2], count[2];
  int          counter, c;

  printf ("WriteSingleRealSpField: %s\n", filename);

  // Copy over dimension lengths to correct type
  nlon1 = nlon;
  nlat1 = nlat;

  // Create a new netCDF file
  ierr[0] = nc_create (filename,
                       NC_CLOBBER,
                       &ncid);

  // Define the dimensions
  ierr[1] = nc_def_dim (ncid,
                        "longitude",
                        nlon1,
                        &dimidx);
  ierr[2] = nc_def_dim (ncid,
                        "latitude",
                        nlat1,
                        &dimidy);

  // Define the variables
  list1[0] = dimidx;
  ierr[3]  = nc_def_var (ncid,
                         "longitude",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidx);
  list1[0] = dimidy;
  ierr[4]  = nc_def_var (ncid,
                         "latitude",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidy);

  // x, y directions correspond to 1, 0 respectively.
  list2[0] = dimidy;
  list2[1] = dimidx;
  ierr[5] = nc_def_var (ncid,
                        varname,
                        NC_DOUBLE,
                        2,
                        list2,
                        &varid);

  // End define mode
  ierr[6] = nc_enddef (ncid);

  // Output data
  ierr[7] = nc_put_var_double (ncid,
                               varidx,
                               lons);

  ierr[8] = nc_put_var_double (ncid,
                               varidy,
                               lats);

  // Error reporting
  for (c=0; c<=8; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error stage A (%i) %i: %s\n", c, err, nc_strerror(err));
      exit(0);
    }
  }


  // x, y directions correspond to 1, 0 respectively.
  start[1] = 0;
  count[0] = 1;
  count[1] = nlon;
  counter  = 7;
  // Output main field
  for (y=0; y<nlat; y++)
  { start[0] = y;
    // Copy part of data into special array
    for (x=0; x<nlon; x++)
    { lineofdata[x] = field[x][y];
    }
    counter++;
    err = nc_put_vara_double (ncid,
                              varid,
                              start,
                              count,
                              lineofdata);
    if (err != 0)
    { printf ("netcdf error writing data: %i %s\n", err, nc_strerror(err));
      exit(0);
    }
  }

  //Close the netCDF file
  counter++;
  err = nc_close (ncid);

  if (err != 0)
  { printf ("netcdf error closing netcdf file: %i %s\n", err, nc_strerror(err));
    exit(0);
  }

}


// -------------------------------------------------------------------------------
void WriteSingleSpecSpField
// Function to write-out a single spectral-space field
 ( int    L,
   double ***field,
   char   filename[256],
   char   varname[32] )

{ //Declare local variables
  int          l, m, p;
  // netCDF-related variables
  size_t       nLp1;
  double       lineofdata[L+1];
  int          ierr[2*(L+1)+13], ncid, err;
  int          dimidp, dimidl, dimidm;
  int          varidp, varidl, varidm, varid;
  int          list1[1], list3[3];
  size_t       start[3], count[3];
  int          dimensionseq[L+1];
  int          counter, c;

  printf ("WriteSingleRealSpField: %s\n", filename);

  // Copy over dimension lengths to correct type
  nLp1 = L + 1;

  // Set-up the sequence of numbers that will define the dimensions
  for (c=0; c<=L; c++)
  { dimensionseq[c] = c;
  }

  // Create a new netCDF file
  ierr[0] = nc_create (filename,
                       NC_CLOBBER,
                       &ncid);

  // Define the dimensions
  ierr[1] = nc_def_dim (ncid,
                        "realimag",
                        2,
                        &dimidp);
  ierr[2] = nc_def_dim (ncid,
                        "l",
                        nLp1,
                        &dimidl);
  ierr[3] = nc_def_dim (ncid,
                        "m",
                        nLp1,
                        &dimidm);

  // Define the variables
  list1[0] = dimidp;
  ierr[4]  = nc_def_var (ncid,
                         "realimag",
                         NC_INT,
                         1,
                         list1,
                         &varidp);
  list1[0] = dimidl;
  ierr[5]  = nc_def_var (ncid,
                         "l",
                         NC_INT,
                         1,
                         list1,
                         &varidl);
  list1[0] = dimidm;
  ierr[6]  = nc_def_var (ncid,
                         "m",
                         NC_INT,
                         1,
                         list1,
                         &varidm);

  // x, y, z directions correspond to 2, 1, 0 respectively.
  list3[0] = dimidm;
  list3[1] = dimidl;
  list3[2] = dimidp;
  ierr[7]  = nc_def_var (ncid,
                         varname,
                         NC_DOUBLE,
                         3,
                         list3,
                         &varid);

  // End define mode
  ierr[8] = nc_enddef (ncid);

  // Output data
  ierr[9]  = nc_put_var_int (ncid,
                             varidp,
                             dimensionseq);

  ierr[10] = nc_put_var_int (ncid,
                             varidl,
                             dimensionseq);

  ierr[11] = nc_put_var_int (ncid,
                             varidm,
                             dimensionseq);

  // x, y, z directions correspond to 2, 1, 0 respectively.
  start[1] = 0;
  count[0] = 1;
  count[1] = L+1;
  count[2] = 1;
  counter  = 11;
  // Output main field
  for (p=0; p<2; p++)
  { start[2] = p;
    for (m=0; m<=L; m++)
    { start[0] = m;
      // Copy part of data into special array
      for (l=0; l<=L; l++)
      { if ((m == 0) && (p == 1))  // no sin term for mm=0
        { lineofdata[l] = 0.0;
        }
        else
        { if (l >= m)
          { lineofdata[l] = field[p][l][m];
          }
          else
          { lineofdata[l] = 0.0;
          }
        }
      }
      counter++;
      ierr[counter] = nc_put_vara_double (ncid,
                                          varid,
                                          start,
                                          count,
                                          lineofdata);
    }
  }

  //Close the netCDF file
  counter++;
  ierr[counter] = nc_close (ncid);

  // Error reporting
  for (c=0; c<=counter; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error %i: %s\n", err, nc_strerror(err));
    }
  }
}



// -------------------------------------------------------------------------------
void cvt_matrices_input (struct CVTData_type *CVTdata,
                         char                filename[256])
{ // Read-in CVT data from cvt file
  // This includes allocation of the CVT data structure if it is not yet allocated

  int    ierr[14], err, ncid;
  int    dimid_x_sh, dimid_y_sh, dimid_z_tomcat, dimid_month;
  int    varid_x_sh, varid_y_sh, varid_z_tomcat, varid_month;
  int    varid_vert_covs, varid_temporal_covs;
  int    varid_std_tracer, varid_std_hspec_tracer, varid_vert_eigenvec_tracer;
  int    varid_std_vspec_tracer, varid_std_flux, varid_std_hspec_flux;
  int    varid_temp_eigenvec_flux, varid_std_tspec_flux;
  int    nlon, nlat, nlev, nss, L;
  int    x, y, c, z, zp, k, t, counter;
  int    vert_covs, temporal_covs;
  size_t n, start[4], count[4];
  double *lodlev;
  double *lodwn;
  double *lodtime;
  double *lodlon;
  double dnhoriz;

  printf ("Reading cvt file: %s\n", filename);

  // Open the file for reading
  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);
  if (err != 0)
  { printf ("ERROR opening file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }

  // Get the variable ids
  ierr[0]  = nc_inq_varid (ncid, "x_sh",                 &varid_x_sh);
  ierr[1]  = nc_inq_varid (ncid, "y_sh",                 &varid_y_sh);
  ierr[2]  = nc_inq_varid (ncid, "z_tomcat",             &varid_z_tomcat);
  ierr[3]  = nc_inq_varid (ncid, "month",                &varid_month);
  ierr[4]  = nc_inq_varid (ncid, "std_tracer",           &varid_std_tracer);
  ierr[5]  = nc_inq_varid (ncid, "std_hspec_tracer",     &varid_std_hspec_tracer);
  ierr[6]  = nc_inq_varid (ncid, "vert_covs",            &varid_vert_covs);
  ierr[7]  = nc_inq_varid (ncid, "vert_eigenvec_tracer", &varid_vert_eigenvec_tracer);
  ierr[8]  = nc_inq_varid (ncid, "std_vspec_tracer",     &varid_std_vspec_tracer);
  ierr[9]  = nc_inq_varid (ncid, "std_flux",             &varid_std_flux);
  ierr[10] = nc_inq_varid (ncid, "std_hspec_flux",       &varid_std_hspec_flux);
  ierr[11] = nc_inq_varid (ncid, "temporal_covs",        &varid_temporal_covs);
  ierr[12] = nc_inq_varid (ncid, "temp_eigenvec_flux",   &varid_temp_eigenvec_flux);
  ierr[13] = nc_inq_varid (ncid, "std_tspec_flux",       &varid_std_tspec_flux);

  counter = 0;
  for (c=0; c<14; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting variable id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}

  // Get the dimension sizes
  ierr[0] = nc_inq_dimid  (ncid, "x_sh", &dimid_x_sh);
  ierr[1] = nc_inq_dimlen (ncid, dimid_x_sh, &n);
  nlon    = n;
  ierr[2] = nc_inq_dimid  (ncid, "y_sh", &dimid_y_sh);
  ierr[3] = nc_inq_dimlen (ncid, dimid_y_sh, &n);
  nlat    = n;
  ierr[4] = nc_inq_dimid  (ncid, "z_tomcat", &dimid_z_tomcat);
  ierr[5] = nc_inq_dimlen (ncid, dimid_z_tomcat, &n);
  nlev    = n;
  ierr[6] = nc_inq_dimid  (ncid, "month", &dimid_month);
  ierr[7] = nc_inq_dimlen (ncid, dimid_month, &n);
  nss     = n;
  L       = nlat - 1;

  counter = 0;
  for (c=0; c<8; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting data from CVT file to check consistency with toy problem %u: %i (%s)\n",
              c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}


  if ((*CVTdata).times)
  { // The CVT data structure has already been allocated
    // Check that the dimension sizes are correct
    if ( ((*CVTdata).L != L) || ((*CVTdata).nlon != nlon) || ((*CVTdata).nlat != nlat) ||
         ((*CVTdata).nlev != nlev) || ((*CVTdata).nss != nss) )
    { printf ("Error from routine cvt_matrices_input: The allocated CVT data structure dimensions do not match file\n");
      printf ("L:    allocation %i, CVT file %i\n", (*CVTdata).L,    L);
      printf ("nlon: allocation %i, CVT file %i\n", (*CVTdata).nlon, nlon);
      printf ("nlat: allocation %i, CVT file %i\n", (*CVTdata).nlat, nlat);
      printf ("nlev: allocation %i, CVT file %i\n", (*CVTdata).nlev, nlev);
      printf ("nss:  allocation %i, CVT file %i\n", (*CVTdata).nss,  nss);
    }
  }
  else
  { // Allocate the CVT structure
    Allocate_CVT (CVTdata,
                  L,
                  nlev,
                  nss);
  }

  // Allocate structures to hold lines of data
  lodlev  = new double[nlev];
  lodwn   = new double[L+1];
  lodtime = new double[nss];
  lodlon  = new double[nlon];
  dnhoriz = double(((*CVTdata).nlon-1) * ((*CVTdata).nlat-1));

  // Get values of these variables
  // --------------------------------------------------------

  // ===== The dimension variables =====
  start[0] = 0;
  count[0] = (*CVTdata).nlon;
  ierr[0]  = nc_get_vara_double (ncid,
                                 varid_x_sh,
                                 start, count,
                                 (*CVTdata).longitude);
  count[0] = (*CVTdata).nlat;
  ierr[1]  = nc_get_vara_double (ncid,
                                 varid_y_sh,
                                 start, count,
                                 (*CVTdata).latitude);

  count[0] = (*CVTdata).nlev;
  ierr[2]  = nc_get_vara_double (ncid,
                                 varid_z_tomcat,
                                 start, count,
                                 lodlev);
  for (z=0; z<(*CVTdata).nlev; z++)
  { zp = (*CVTdata).nlev - z - 1;                  // Want to make the bottom level lev=0
    (*CVTdata).level[zp] = lodlev[z];              // (input data is opposite)
  }
  count[0] = (*CVTdata).nss;
  ierr[3]  = nc_get_vara_double (ncid,
                                 varid_month,
                                 start, count,
                                 (*CVTdata).times);
  counter = 0;
  for (c=0; c<4; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting dimension variable %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}



  // ===== Stdev of the tracer =====
  // Take the global average vertical profile
  for (z=0; z<(*CVTdata).nlev; z++)
  { (*CVTdata).tracer_stddev[z] = 0.0;
  }
  count[2] = 1;                               // lat
  count[1] = 1;                               // lon
  start[0] = 0;  count[0] = (*CVTdata).nlev;  // lev
  for (y=0; y<(*CVTdata).nlat-1; y++)         // There is one less latitude in the CVT file data
  { start[2] = y;
    for (x=0; x<(*CVTdata).nlon-1; x++)       // There is one less longitude in the CVT file data
    { start[1] = x;
      err = nc_get_vara_double (ncid,
                                varid_std_tracer,
                                start, count,
                                lodlev);
      if (err != 0)
      { printf ("ERROR getting varid_std_tracer %u %u: %i (%s)\n", y, x, err, nc_strerror(err));
        exit(0);
      }
      // Add to global average structure
      for (z=0; z<(*CVTdata).nlev; z++)
      { zp = (*CVTdata).nlev - z - 1;                  // Want to make the bottom level lev=0
        (*CVTdata).tracer_stddev[zp] += lodlev[z] ;    // (input data is opposite)
      }
    }
  }
  // Normalise
  for (z=0; z<(*CVTdata).nlev; z++)
  { (*CVTdata).tracer_stddev[z] /= dnhoriz;
    //printf ("Stdev of the tracer %u = %e\n", z, (*CVTdata).tracer_stddev[z]);
  }


  // ===== Stdev horizontal spectrum for the tracer =====
  // Note 1=wn, 0=lev
                 count[1] = 1;                // wn
  start[0] = 0;  count[0] = (*CVTdata).nlev;  // lev
  for (k=0; k<=(*CVTdata).L; k++)
  { start[1] = k;
    err = nc_get_vara_double (ncid,
                              varid_std_hspec_tracer,
                              start, count,
                              lodlev);
    if (err != 0)
    { printf ("ERROR getting varid_std_hspec_tracer %u: %i (%s)\n", k, err, nc_strerror(err));
      printf ("L = %i\n", (*CVTdata).L);
      exit(0);
    }
    // Copy to correct place
    for (z=0; z<(*CVTdata).nlev; z++)
    { zp = (*CVTdata).nlev - z - 1;                  // Want to make the bottom level lev=0
      (*CVTdata).tracer_hspec[k][zp] = lodlev[z];
      //printf ("Stdev of the hspec tracer %u %u = %e\n", k, zp, (*CVTdata).tracer_hspec[k][zp]);
    }
  }


  // ===== Vertical eigenvectors for the tracer =====

  // First check that vert_covs option is set to 3
  start[0] = 0;  count[0] = 1;
  err = nc_get_vara_int (ncid,
                         varid_vert_covs,
                         start, count,
                         &vert_covs);
  if (err != 0)
  { printf ("ERROR getting varid_vert_covs: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
  if (vert_covs != 3)
  { printf ("Option for vert_covs set incorrectly in the input covariances file\n");
    printf ("This is currently set to %i, but should be set to 3\n", vert_covs);
    exit(0);
  }


  start[1] = 0;  count[1] = (*CVTdata).nlev;  // lev
                 count[0] = 1;                // mode
  for (k=0; k<(*CVTdata).nlev; k++)
  { start[0] = k;
    err = nc_get_vara_double (ncid,
                              varid_vert_eigenvec_tracer,
                              start, count,
                              lodlev);
    if (err != 0)
    { printf ("ERROR getting varid_vert_eigenvec_tracer %u: %i (%s)\n", k, err, nc_strerror(err));
      exit(0);
    }
    // Copy to correct place
    for (z=0; z<(*CVTdata).nlev; z++)
    { zp = (*CVTdata).nlev - z - 1;                  // Want to make the bottom level lev=0
      (*CVTdata).tracer_vertmodes[zp][k] = lodlev[z];
    //printf ("Vert eigenvec tracer %u %u = %e\n", zp, k, (*CVTdata).tracer_vertmodes[zp][k]);
    }
  }


  // ===== Vertical spectrum of standard deviations for the tracer =====
  // Take the first horizontal position
  start[2] = 0; count[2] = (*CVTdata).nlev;  // mode
  start[1] = 0, count[1] = 1;                // lat
  start[0] = 0, count[0] = 1;                // long
  err      = nc_get_vara_double (ncid,
                                 varid_std_vspec_tracer,
                                 start, count,
                                 lodlev);
  if (err != 0)
  { printf ("ERROR getting varid_std_vspec_tracer: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
  // Copy to correct place
  for (k=0; k<(*CVTdata).nlev; k++)
  { (*CVTdata).tracer_vspec[k] = lodlev[k];
    //printf ("Stdev of vspec tracer = %e\n", (*CVTdata).tracer_vspec[k]);
  }


  // ===== Standard deviations for the flux =====
  // Note 3=lat, 2=long, 1=month, 0=field
  // Take the first horizontal postition, month, and field.
                 count[3] = 1;                  // lat
  start[2] = 0;  count[2] = (*CVTdata).nlon-1;  // long, there is one less longitude in the CVT file data
                 count[1] = 1;                  // month
  start[0] = 0;  count[0] = 1;                  // field
  for (t=0; t<(*CVTdata).nss; t++)
  { start[1] = t;
    for (y=0; y<(*CVTdata).nlat-1; y++)         // There is one less latitude in the CVT file data
    { start[3] = y;
      err      = nc_get_vara_double (ncid,
                                     varid_std_flux,
                                     start, count,
                                     lodlon);
      // Provide extra longitude's data (use periodicity)
      lodlon[nlon-1] = lodlon[0];
      for (x=0; x<(*CVTdata).nlon; x++)
      { (*CVTdata).source_stddev[x][y][t] = lodlon[x];
      }
    }
    // Provide extra latitude's data
    for (x=0; x<(*CVTdata).nlon; x++)
    { (*CVTdata).source_stddev[x][(*CVTdata).nlat-1][t] = (*CVTdata).source_stddev[x][(*CVTdata).nlat-2][t];
    }
  }
  if (err != 0)
  { printf ("ERROR getting varid_std_flux: %i (%s)\n", err, nc_strerror(err));
      exit(0);
  }

//  // Check that standard deviation has been input correctly
//  double **field;
//  Array_2d_double_create (&field,
//                          nlon,
//                          nlat);
//  for (x=0; x<(*CVTdata).nlon; x++)
//  { for (y=0; y<(*CVTdata).nlat; y++)
//    { field[x][y] = (*CVTdata).source_stddev[x][y][0];
//    }
//  }
//  WriteSingleRealSpField ( nlon,
//                           nlat,
//                           field,
//                           (*CVTdata).longitude,
//                           (*CVTdata).latitude,
//                           "Source_std.nc",
//                           "Source_std" );
//  delete [] field;

  // ===== Stddev horiz spectrum for the flux =====
  // Note 2=wn, 1=month, 0=field
  // Assume that the values are independent of time, so read month 1
  // Read data for month 1, field 1
  start[2] = 0; count[2] = (*CVTdata).L+1;  // wn
  start[1] = 0, count[1] = 1;               // timescale
  start[0] = 0, count[0] = 1;               // field
  err = nc_get_vara_double (ncid,
                            varid_std_hspec_flux,
                            start, count,
                            lodwn);
  if (err != 0)
  { printf ("ERROR getting varid_std_hspec_flux: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
  // Copy to correct place
  for (k=0; k<=(*CVTdata).L; k++)
  { (*CVTdata).source_hspec[k] = lodwn[k];
    //printf ("Std of the hspec source %u = %e\n", k, (*CVTdata).source_hspec[k]);
  }


  // ===== Temporal eigenvectors for the flux =====

  // First check that temporal_covs option is set to 0 or 1
  start[0] = 0;  count[0] = 1;
  err = nc_get_vara_int (ncid,
                         varid_temporal_covs,
                         start, count,
                         &temporal_covs);
  if (err != 0)
  { printf ("ERROR getting varid_temporal_covs: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
  if ((temporal_covs != 0) && (temporal_covs != 1))
  { printf ("Option for temporal_covs set incorrectly in the input covariances file\n");
    printf ("This is currently set to %i, but should be set to 0 or 1\n", temporal_covs);
    exit(0);
  }

  (*CVTdata).temporal_covs = (temporal_covs == 1);


  // Note 2=month, 1=temporal mode number, 0=field
  // Read data for field 1
  start[2] = 0; count[2] = (*CVTdata).nss;  // month
                count[1] = 1;               // temporal mode
  start[0] = 0, count[0] = 1;               // field
  for (k=0; k<(*CVTdata).nss; k++)
  { start[1] = k;
    err = nc_get_vara_double (ncid,
                              varid_temp_eigenvec_flux,
                              start, count,
                              lodtime);
    if (err != 0)
    { printf ("ERROR getting varid_temp_eigenvec_flux %u: %i (%s)\n", k, err, nc_strerror(err));
      exit(0);
    }
    // Copy to correct place
    for (t=0; t<(*CVTdata).nss; t++)
    { (*CVTdata).source_tempmodes[t][k] = lodtime[t];
      //printf ("Temp eigenvec source %u %u = %e\n", t, k, (*CVTdata).source_tempmodes[t][k]);
    }
  }


  // ===== Temporal mode standard deviations for the flux =====
  // Read data for field 1
  start[1] = 0;  count[1] = (*CVTdata).nss;  // month
  start[0] = 0;  count[0] = 1;               // field
  err = nc_get_vara_double (ncid,
                            varid_std_tspec_flux,
                            start, count,
                            lodtime);
  if (err != 0)
  { printf ("ERROR getting varid_std_tspec_flux: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
  // Copy to correct place
  for (k=0; k<(*CVTdata).nss; k++)
  { (*CVTdata).source_tspec[k] = lodtime[k];
    //printf ("Std of the tspec source %u = %e\n", k, (*CVTdata).source_tspec[k]);
  }


  err = nc_close (ncid);
  if (err != 0)
  { printf ("ERROR closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

  // Deallocate structures used to hold lines of data
  delete [] lodlev;
  delete [] lodwn;
  delete [] lodtime;
  delete [] lodlon;
}



// -------------------------------------------------------------------------------
void WriteStateVector
// Function to write-out the state vector (any representation)
 ( struct state_type *field,
   char              filename[256] )

{ //Declare local variables
  int          x, y, l, m, cs, z, t, xp1, yp1, zp1;
  size_t       nlon, nlat, nlev, nss, L;
  double       lineofdata[(*field).nlon];
  double       Wns_and_cs[(*field).L+1];
  double       vertmodes[(*field).nlev];
  double       tempmodes[(*field).nss];
  int          ierr[7], ncid, err, c;
  int          dimid_lon, dimid_lat, dimid_lev, dimid_l, dimid_m;
  int          dimid_cs, dimid_t, dimid_vertmode, dimid_tempmode;
  int          varid_lon, varid_lat, varid_lev, varid_l, varid_m;
  int          varid_cs, varid_t, varid_vertmode, varid_tempmode;
  int          varid_tracer0, varid_source;
  int          list[4];
  size_t       start[4], count[4];

  printf ("WriteStateVector: %s\n", filename);

  // Copy over dimension lengths to correct type
  nlon = (*field).nlon;
  nlat = (*field).nlat;
  nlev = (*field).nlev;
  nss  = (*field).nss;
  L    = (*field).L;

  // Create a new netCDF file
  err = nc_create (filename,
                   NC_CLOBBER,
                   &ncid);
  if (err != 0)
  { printf ("netcdf error opening output file: %i (%s)\n", err, nc_strerror(err));
  }

  // Define the dimensions
  if ((*field).horiz_repres == 'r')
  { ierr[0] = nc_def_dim (ncid,
                          "longitude",
                          nlon,
                          &dimid_lon);
    ierr[1] = nc_def_dim (ncid,
                          "latitude",
                          nlat,
                          &dimid_lat);
    ierr[2] = 0;
  }
  else
  { ierr[0] = nc_def_dim (ncid,
                          "total_wn",
                          L+1,
                          &dimid_l);
    ierr[1] = nc_def_dim (ncid,
                          "zonal_wn",
                          L+1,
                          &dimid_m);
    ierr[2] = nc_def_dim (ncid,
                          "cossin",
                          2,
                          &dimid_cs);
  }

  if ((*field).vert_repres == 'r')
  { ierr[3] = nc_def_dim (ncid,
                          "level",
                          nlev,
                          &dimid_lev);
  }
  else
  { ierr[3] = nc_def_dim (ncid,
                          "vert_mode",
                          nlev,
                          &dimid_vertmode);
  }
  if ((*field).temp_repres == 'r')
  { ierr[4] = nc_def_dim (ncid,
                          "time",
                          nss,
                          &dimid_t);
  }
  else
  { ierr[4] = nc_def_dim (ncid,
                          "temp_mode",
                          nss,
                          &dimid_tempmode);
  }

  for (c=0; c<=4; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error defining dimensions %u: %i (%s)\n", c, err, nc_strerror(err));
    }
  }


  // Define the variables
  if ((*field).horiz_repres == 'r')
  { // We are in real space in the horizontal
    list[0] = dimid_lon;
    ierr[0] = nc_def_var (ncid,
                          "longitude",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_lon);
    list[0] = dimid_lat;
    ierr[1] = nc_def_var (ncid,
                          "latitude",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_lat);
    // 2, 1, 0 correspond to long, lat, height
    ierr[2] = 0;
    list[2] = dimid_lon;
    list[1] = dimid_lat;
    if ((*field).vert_repres == 'r')
    { list[0] = dimid_lev;
    }
    else
    { list[0] = dimid_vertmode;
    }
    ierr[3] = nc_def_var (ncid,
                          "tracer0",
                          NC_DOUBLE,
                          3,
                          list,
                          &varid_tracer0);
    // 2, 1, 0 correspond to long, lat, time
    list[2] = dimid_lon;
    list[1] = dimid_lat;
    if ((*field).temp_repres == 'r')
    { list[0] = dimid_t;
    }
    else
    { list[0] = dimid_tempmode;
    }
    ierr[4] = nc_def_var (ncid,
                          "source",
                          NC_DOUBLE,
                          3,
                          list,
                          &varid_source);
  }
  else
  { // We are in spectral space in the horizontal
    list[0] = dimid_l;
    ierr[0] = nc_def_var (ncid,
                          "total_wn",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_l);
    list[0] = dimid_m;
    ierr[1] = nc_def_var (ncid,
                          "zonal_wn",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_m);
    list[0] = dimid_cs;
    ierr[2] = nc_def_var (ncid,
                          "cossin",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_cs);
    // 3, 2, 1, 0 correspond to l, m, cs, z
    list[3] = dimid_l;
    list[2] = dimid_m;
    list[1] = dimid_cs;
    if ((*field).vert_repres == 'r')
    { list[0] = dimid_lev;
    }
    else
    { list[0] = dimid_vertmode;
    }
    ierr[3] = nc_def_var (ncid,
                          "tracer0",
                          NC_DOUBLE,
                          4,
                          list,
                          &varid_tracer0);
    // 3, 2, 1, 0 correspond to l, m, cs, t
    list[3] = dimid_l;
    list[2] = dimid_m;
    list[1] = dimid_cs;
    if ((*field).temp_repres == 'r')
    { list[0] = dimid_t;
    }
    else
    { list[0] = dimid_tempmode;
    }
    ierr[4] = nc_def_var (ncid,
                          "source",
                          NC_DOUBLE,
                          4,
                          list,
                          &varid_source);
  }

  if ((*field).vert_repres == 'r')
  { list[0] = dimid_lev;
    ierr[5] = nc_def_var (ncid,
                          "level",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_lev);
  }
  else
  { list[0] = dimid_vertmode;
    ierr[5] = nc_def_var (ncid,
                          "vert_mode",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_vertmode);
  }

  if ((*field).temp_repres == 'r')
  { list[0] = dimid_t;
    ierr[6] = nc_def_var (ncid,
                          "time",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_t);
  }
  else
  { list[0] = dimid_tempmode;
    ierr[6] = nc_def_var (ncid,
                          "temp_mode",
                          NC_DOUBLE,
                          1,
                          list,
                          &varid_tempmode);
  }

  for (c=0; c<=6; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error defining variables %u: %i (%s)\n", c, err, nc_strerror(err));
    }
  }


  // End define mode
  err = nc_enddef (ncid);
  if (err != 0)
  { printf ("netcdf error changing mode: %i (%s)\n", err, nc_strerror(err));
  }


  // Output dimension data
  if ((*field).horiz_repres == 'r')
  { ierr[0]  = nc_put_var_double (ncid,
                                  varid_lon,
                                  &((*field).longitude[1]));
    ierr[1]  = nc_put_var_double (ncid,
                                  varid_lat,
                                  &((*field).latitude[1]));
    ierr[2]  = 0;
  }
  else
  { for (l=0; l<=L; l++)
    { Wns_and_cs[l] = double(l);
    }
    ierr[0]  = nc_put_var_double (ncid,
                                  varid_l,
                                  Wns_and_cs);
    ierr[1]  = nc_put_var_double (ncid,
                                  varid_m,
                                  Wns_and_cs);
    ierr[2]  = nc_put_var_double (ncid,
                                  varid_cs,
                                  Wns_and_cs);
  }

  if ((*field).vert_repres == 'r')
  { ierr[3] = nc_put_var_double (ncid,
                                 varid_lev,
                                 &((*field).level[1]));
  }
  else
  { for (z=0; z<nlev; z++)
    { vertmodes[z] = double(z+1);
    }
    ierr[3] = nc_put_var_double (ncid,
                                 varid_vertmode,
                                 vertmodes);
  }

  if ((*field).temp_repres == 'r')
  { ierr[4] = nc_put_var_double (ncid,
                                 varid_t,
                                 (*field).times);
  }
  else
  { for (t=0; t<nss; t++)
    { tempmodes[t] = double(t+1);
    }
    ierr[4] = nc_put_var_double (ncid,
                                 varid_tempmode,
                                 tempmodes);
  }

  for (c=0; c<=4; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error outputting dimension variables %u: %i (%s)\n", c, err, nc_strerror(err));
    }
  }

  // Output field data
  if ((*field).horiz_repres == 'r')
  { // Real space
    // x, y, z directions correspond to 2, 1, 0 respectively.
    // x, y, t directions correspond to 2, 1, 0 respectively.
    start[2] = 0;
    count[2] = nlon;
    count[1] = 1;
    count[0] = 1;
    for (z=0; z<nlev; z++)
    { zp1      = z + 1;
      start[0] = z;
      for (y=0; y<nlat; y++)
      { yp1 = y + 1;
        start[1] = y;
        // Copy part of data into special array
        for (x=0; x<nlon; x++)
        { lineofdata[x] = (*field).tracer0_rs[x+1][yp1][zp1];
        }
        err = nc_put_vara_double (ncid,
                                  varid_tracer0,
                                  start,
                                  count,
                                  lineofdata);
        if (err != 0)
        { printf ("netcdf error outputting tracer0_rs: %i (%s)\n", err, nc_strerror(err));
        }
      }
    }
    for (t=0; t<nss; t++)
    { start[0] = t;
      for (y=0; y<nlat; y++)
      { yp1      = y + 1;
        start[1] = y;
        // Copy part of data into special array
        for (x=0; x<nlon; x++)
        { lineofdata[x] = (*field).source_rs[x+1][yp1][t];
        }
        err = nc_put_vara_double (ncid,
                                  varid_source,
                                  start,
                                  count,
                                  lineofdata);
        if (err != 0)
        { printf ("netcdf error outputting source_rs: %i (%s)\n", err, nc_strerror(err));
        }
      }
    }
  }
  else
  { // Spectral space
    // l, m, cs, z directions correspond to 3, 2, 1, 0 respectively.
    // l, m, cs, t directions correspond to 3, 2, 1, 0 respectively.
    start[3] = 0;
    count[3] = L+1;
    count[2] = 1;
    count[1] = 1;
    count[0] = 1;
    for (z=0; z<nlev; z++)
    { zp1      = z + 1;
      start[0] = z;
      for (m=0; m<=L; m++)
      { start[2] = m;
        for (cs=0; cs<2; cs++)
        { start[1] = cs;
          // Copy part of data into special array
          for (l=0; l<=L; l++)
          { lineofdata[l] = (*field).tracer0_ss[cs][l][m][zp1];
          }
          err = nc_put_vara_double (ncid,
                                    varid_tracer0,
                                    start,
                                    count,
                                    lineofdata);
          if (err != 0)
          { printf ("netcdf error outputting tracer0_ss: %i (%s)\n", err, nc_strerror(err));
          }
        }
      }
    }
    for (t=0; t<nss; t++)
    { start[0] = t;
      for (m=0; m<=L; m++)
      { start[2] = m;
        for (cs=0; cs<2; cs++)
        { start[1] = cs;
          // Copy part of data into special array
          for (l=0; l<=L; l++)
          { lineofdata[l] = (*field).source_ss[cs][l][m][t];
          }
          err = nc_put_vara_double (ncid,
                                    varid_source,
                                    start,
                                    count,
                                    lineofdata);
          if (err != 0)
          { printf ("netcdf error outputting source_ss: %i (%s)\n", err, nc_strerror(err));
          }
        }
      }
    }
  }

  //Close the netCDF file
  err = nc_close (ncid);

  // Error reporting
  if (err != 0)
  { printf ("netcdf error closing file: %i (%s)\n", err, nc_strerror(err));
  }
}


// -------------------------------------------------------------------------------
void ReadStateVector (struct state_type    *field,
                      struct metadata_type *metadata,
                      bool                 Read_tracer,
                      bool                 Read_source,
                      char                 filename[256])
{ // Read state vector
  // This routine will also allocate the data structures if it has not already been done.
  // The field must specify on input the representations of the fields (real, etc. as below)
  // Different representations are possible, as specified in (*field).horiz_repres,
  //   (*field).vert_repres, and (*field).temp_repres i.e.
  // (*field).horiz_repres = 'r'eal space, 's'pectral space
  // (*field).vert_repres  = 'r'eal space, 'm'odal space
  // (*field).temp_repres  = 'r'eal space, 'm'odal space

  int    ierr[4], err, ncid;
  int    dimid_lon, dimid_lat, dimid_l, dimid_lev, dimid_vertmode, dimid_t, dimid_tempmode;
  int    varid_lon, varid_lat, varid_lev, varid_t, varid_tracer0, varid_source;
  int    x, y, c, z, counter, yp1, zp1, m, cs, l, t;
  int    nlon_read=0, nlat_read=0, nlev_read=0, L_read=0, nss_read=0;
  size_t start[3], count[3], in;
  double *lod;
  bool   dontagree, preexist_field, preexist_meta;

  printf ("Reading state vector file: %s\n", filename);

  // Open the file for reading
  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);

  if (err != 0)
  { printf ("ERROR opening file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }


  // The 'times' array is an indicator if structures have already been allocated
  preexist_field = (*field).times;
  preexist_meta  = (*metadata).times;

  // Get the dimension lengths - put the results into first into metadata
  if ((*field).horiz_repres == 'r')
  { // Get longitude and latitude
    ierr[0] = nc_inq_dimid  (ncid, "longitude", &dimid_lon);
    ierr[1] = nc_inq_dimlen (ncid, dimid_lon, &in);
    nlon_read = in;
    ierr[2] = nc_inq_dimid  (ncid, "latitude", &dimid_lat);
    ierr[3] = nc_inq_dimlen (ncid, dimid_lat, &in);
    nlat_read = in;
    L_read    = nlat_read - 1;
    counter = 0;
    for (c=0; c<4; c++)
    { if (ierr[c] != 0)
      { printf ("ERROR getting horizontal dimension lengths %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
        counter++;
      }
    }
    if (counter > 0) {exit(0);}
  }
  else
  { // Get maximum total wavenumber
    ierr[0] = nc_inq_dimid  (ncid, "total_wn", &dimid_l);
    ierr[1] = nc_inq_dimlen (ncid, dimid_l, &in);
    L_read    = in - 1;
    nlon_read = 2*L_read + 1;
    nlat_read = L_read + 1;
    counter = 0;
    for (c=0; c<2; c++)
    { if (ierr[c] != 0)
      { printf ("ERROR getting horizontal max wavenumber %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
        counter++;
      }
    }
    if (counter > 0) {exit(0);}
  }

  if ((*field).vert_repres == 'r')
  { // Get number of levels
    ierr[0] = nc_inq_dimid  (ncid, "level", &dimid_lev);
    ierr[1] = nc_inq_dimlen (ncid, dimid_lev, &in);
    nlev_read = in;
    counter = 0;
    for (c=0; c<2; c++)
    { if (ierr[c] != 0)
      { printf ("ERROR getting number of vertical levels %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
        counter++;
      }
    }
    if (counter > 0) {exit(0);}
  }
  else
  { // Get number of vertical modes
    ierr[0] = nc_inq_dimid  (ncid, "vert_mode", &dimid_vertmode);
    ierr[1] = nc_inq_dimlen (ncid, dimid_vertmode, &in);
    nlev_read = in;
    counter = 0;
    for (c=0; c<2; c++)
    { if (ierr[c] != 0)
      { printf ("ERROR getting number of vertical modes %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
        counter++;
      }
    }
    if (counter > 0) {exit(0);}
  }

  if ((*field).temp_repres == 'r')
  { // Get number of source/sink fields
    ierr[0] = nc_inq_dimid  (ncid, "time", &dimid_t);
    ierr[1] = nc_inq_dimlen (ncid, dimid_t, &in);
    nss_read = in;
    counter = 0;
    for (c=0; c<2; c++)
    { if (ierr[c] != 0)
      { printf ("ERROR getting number of source/sink fields %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
        counter++;
      }
    }
    if (counter > 0) {exit(0);}
  }
  else
  { // Get number of temporal modes
    ierr[0] = nc_inq_dimid  (ncid, "temp_mode", &dimid_tempmode);
    ierr[1] = nc_inq_dimlen (ncid, dimid_tempmode, &in);
    nss_read = in;
    counter = 0;
    for (c=0; c<2; c++)
    { if (ierr[c] != 0)
      { printf ("ERROR getting number of temporal modes %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
        counter++;
      }
    }
    if (counter > 0) {exit(0);}
  }


  // Allocate arrays if necessary, or check consistency

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (preexist_meta)
  { // Need to check that metadata agrees with file data
    dontagree = !CheckSameLengths ((*field).horiz_repres,
                                   'M', 'F',                    // Meta data vs File
                                   (*metadata).nlon, nlon_read,
                                   (*metadata).nlat, nlat_read,
                                   (*metadata).nlev, nlev_read,
                                   (*metadata).L,    L_read,
                                   (*metadata).nss,  nss_read);
    if (dontagree)
    { printf ("Quitting\n");
      exit (0);
    }
  }
  else
  { // Allocate the metadata structure
    (*metadata).nlon = nlon_read;
    (*metadata).nlat = nlat_read;
    (*metadata).nlev = nlev_read;
    (*metadata).nss  = nss_read;
    Allocate_metadata (metadata,
                       (*metadata).nlon,
                       (*metadata).nlat,
                       (*metadata).nlev,
                       (*metadata).nss);

    // Get MetaData from the file if necessary
    // Get the dimension data from the file
    if ((*field).horiz_repres == 'r')
    { ierr[0]  = nc_inq_varid (ncid, "longitude", &varid_lon);
      start[0] = 0;  count[0] = (*metadata).nlon;
      ierr[1]  = nc_get_vara_double (ncid,
                                     varid_lon,
                                     start, count,
                                     &((*metadata).longitude[1]));
      (*metadata).longitude[0]                  = 2.0 * (*metadata).longitude[1] - (*metadata).longitude[2];
      (*metadata).longitude[(*metadata).nlon+1] = 2.0 * (*metadata).longitude[(*metadata).nlon] - (*metadata).longitude[(*metadata).nlon-1];

      ierr[2]  = nc_inq_varid (ncid, "latitude", &varid_lat);
      count[0] = (*metadata).nlat;
      ierr[3]  = nc_get_vara_double (ncid,
                                    varid_lat,
                                    start, count,
                                    &((*metadata).latitude[1]));
      (*metadata).latitude[0]                  = 2.0 * (*metadata).latitude[1] - (*metadata).latitude[2];
      (*metadata).latitude[(*metadata).nlat+1] = 2.0 * (*metadata).latitude[(*metadata).nlat] - (*metadata).latitude[(*metadata).nlat-1];

      counter = 0;
      for (c=0; c<4; c++)
      { if (ierr[c] != 0)
        { printf ("ERROR getting horiz dimensions id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
          counter++;
        }
      }
      if (counter > 0) {exit(0);}

      // Complete the cosine latitudes
      //printf ("COMPLETING THE COSINE LATS\n");
      Fill_cos_lats ((*metadata).latitude,
                     (*metadata).nlat,
                     (*metadata).cos_v_lat,
                     (*metadata).cos_u_lat,
                     false );
      //printf ("DONE COMPLETING THE COSINE LATS\n");

    }
    if ((*field).vert_repres == 'r')
    { ierr[0]  = nc_inq_varid (ncid, "level", &varid_lev);
      start[0] = 0;  count[0] = (*metadata).nlev;
      ierr[1]  = nc_get_vara_double (ncid,
                                     varid_lev,
                                     start, count,
                                     &((*metadata).level[1]));
      (*metadata).level[0]                  = 2.0 * (*metadata).level[1] - (*metadata).level[2];
      (*metadata).level[(*metadata).nlev+1] = 2.0 * (*metadata).level[(*metadata).nlev] - (*metadata).level[(*metadata).nlev-1];
      counter = 0;
      for (c=0; c<2; c++)
      { if (ierr[c] != 0)
        { printf ("ERROR getting vert dimensions id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
          counter++;
        }
      }
      if (counter > 0) {exit(0);}
    }
    if ((*field).temp_repres == 'r')
    { ierr[0]  = nc_inq_varid (ncid, "time", &varid_t);
      start[0] = 0;  count[0] = (*metadata).nss;
      ierr[1]  = nc_get_vara_double (ncid,
                                     varid_t,
                                     start, count,
                                     (*metadata).times);
      counter = 0;
      for (c=0; c<2; c++)
      { if (ierr[c] != 0)
        { printf ("ERROR getting temporal dimensions id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
          counter++;
        }
      }
      if (counter > 0) {exit(0);}
    }
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (preexist_field)
  { // Need to check that metadata agrees with file data
    dontagree = !CheckSameLengths ((*field).horiz_repres,
                                   'S', 'F',                    // State data vs File
                                   (*field).nlon, nlon_read,
                                   (*field).nlat, nlat_read,
                                   (*field).nlev, nlev_read,
                                   (*field).L,    L_read,
                                   (*field).nss,  nss_read);
    if (dontagree)
    { printf ("Quitting - the new file doesn't have the same specs as previous files.\n");
      exit (0);
    }
  }
  else
  { // Allocate the field structure
    Allocate_state (field,
                    metadata,
                    (*field).horiz_repres,
                    (*field).vert_repres,
                    (*field).temp_repres);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



  // Get the main variable ids of the data
  if (Read_tracer)
  { ierr[0] = nc_inq_varid (ncid, "tracer0", &varid_tracer0);
  }
  else
  { ierr[0] = 0;
  }
  if (Read_source)
  { ierr[1] = nc_inq_varid (ncid, "source",  &varid_source);
  }
  else
  { ierr[1] = 0;
  }

  counter = 0;
  for (c=0; c<2; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting main variable id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}


  // Get values of variables
  // --------------------------------------------------------

  // ===== The initital tracer =====
  if (Read_tracer)
  { if ((*field).horiz_repres == 'r')
    { // Initial tracer is in real space
      lod = new double[(*field).nlon];
                     count[0] = 1;               // z
                     count[1] = 1;               // y
      start[2] = 0;  count[2] = (*field).nlon;   // x
      for (z=0; z<(*field).nlev; z++)
      { zp1      = z + 1;
        start[0] = z;
        for (y=0; y<(*field).nlat; y++)
        { yp1      = y + 1;
          start[1] = y;
          err = nc_get_vara_double (ncid,
                                    varid_tracer0,
                                    start, count,
                                    lod);
          // Error reporting
          if (err != 0)
          { printf ("netcdf error inputting tracer0 (real space) %i: %s\n", err, nc_strerror(err));
            exit(0);
          }
          for (x=0; x<(*field).nlon; x++)
          { (*field).tracer0_rs[x+1][yp1][zp1] = lod[x];
          }
        }
      }
      delete[] lod;
    }
    else
    { // Initial tracer is in spectral space
      lod = new double[(*field).L+1];
                     count[0] = 1;               // z or vert mode
                     count[1] = 1;               // cs
                     count[2] = 1;               // m
      start[3] = 0;  count[3] = (*field).L + 1;  // l

      for (z=0; z<(*field).nlev; z++)
      { zp1      = z + 1;
        start[0] = z;
        for (m=0; m<=(*field).L; m++)
        { start[2] = m;
          for (cs=0; cs<2; cs++)
          { start[1] = cs;
            err = nc_get_vara_double (ncid,
                                      varid_tracer0,
                                      start, count,
                                      lod);
            // Error reporting
            if (err != 0)
            { printf ("netcdf error inputting tracer0 (spectral space) %i: %s\n", err, nc_strerror(err));
              exit(0);
            }
            for (l=0; l<=(*field).L; l++)
            { (*field).tracer0_ss[cs][l][m][zp1] = lod[l];
            }
          }
        }
      }
      delete[] lod;
    }
  }


  // ===== The surface flux =====
  if (Read_source)
  { if ((*field).horiz_repres == 'r')
    { // Surface flux is in real space
      lod = new double[(*field).nlon];
                     count[0] = 1;               // t or temporal mode
                     count[1] = 1;               // y
      start[2] = 0;  count[2] = (*field).nlon;   // x
      for (t=0; t<(*field).nss; t++)
      { start[0] = t;
        for (y=0; y<(*field).nlat; y++)
        { yp1       = y + 1;
          start[1] = y;
          err = nc_get_vara_double (ncid,
                                    varid_source,
                                    start, count,
                                    lod);
          // Error reporting
          if (err != 0)
          { printf ("netcdf error inputting source (real space) %i: %s\n", err, nc_strerror(err));
            exit(0);
          }
          for (x=0; x<(*field).nlon; x++)
          { (*field).source_rs[x+1][yp1][t] = lod[x];
          }
        }
      }
      delete[] lod;
    }
    else
    { // Surface flux is in spectral space
      lod = new double[(*field).L+1];
                     count[0] = 1;               // t or temporal mode
                     count[1] = 1;               // cs
                     count[2] = 1;               // m
      start[3] = 0;  count[3] = (*field).L + 1;  // l

      for (t=0; t<(*field).nss; t++)
      { start[0] = t;
        for (m=0; m<=L; m++)
        { start[2] = m;
          for (cs=0; cs<2; cs++)
          { start[1] = cs;
            err = nc_get_vara_double (ncid,
                                      varid_source,
                                      start, count,
                                      lod);
            // Error reporting
            if (err != 0)
            { printf ("netcdf error inputting source (spectral space) %i: %s\n", err, nc_strerror(err));
              exit(0);
            }
            for (l=0; l<=L; l++)
            { (*field).source_ss[cs][l][m][t] = lod[l];
            }
          }
        }
      }
      delete[] lod;
    }
  }

  err = nc_close (ncid);
  if (err != 0)
  { printf ("ERROR closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

  // Sort-out halos
  halos (field);

}


// -------------------------------------------------------------------------------
bool CheckSameLengths ( char horiz_repres,      // Horizontal representation type
                        char info1, char info2, // Abbreviations of source of data
                        int  nlon1, int  nlon2, // Two versions of No. of longs
                        int  nlat1, int  nlat2, // Two versions of No. of lats
                        int  nlev1, int  nlev2, // Two versions of No. of levs
                        int  L1,    int  L2,    // Two versions of No. of wavenumbers
                        int  nss1,  int  nss2 ) // Two versions of No. of source/sink fields
{ // Check that two versions of dimension lengths are identical
  bool same;

  if (horiz_repres == 'r')
  { // Real space
    same = (nlon1 == nlon2) && (nlat1 == nlat2);
  }
  else
  { // Spectral space
    same = L1 == L2;
  }

  same = same && (nlev1 == nlev2) && (nss1 == nss2);


  if (!same)
  { printf ("Error from routine ReadStateVector: Data structure dimensions do not match file\n");
    printf ("Comparing %c with %c\n", info1, info2);
    if (horiz_repres == 'r')
    { printf ("nlon: %i, %i\n", nlon1, nlon2);
      printf ("nlat: %i, %i\n", nlat1, nlat2);
    }
    else
    { printf ("L   : %i, %i\n", L1, L2);
    }
    printf ("nlev: %i, %i\n", nlev1, nlev2);
    printf ("nss : %i, %i\n", nss1,  nss2);
  }
  return same;
}



// -------------------------------------------------------------------------------
void ReadINVICAT (struct state_type    *field,
                  struct metadata_type *metadata,
                  bool                 Read_tracer,
                  bool                 Read_source,
                  char                 filename[256])
{ // Read state vector for INVICAT
  // This routine will also allocate the data structures if it has not already been done.
  // Horizontal: the nearest grid box will be taken
  // Vertical: the first set of Envi-flux levels will be used (starting at the top as INVICAT is upside-down)
  // Temporal: the first INVICAT flux field will be used for all times

  int    ierr[6], err, ncid;
  bool   preexist_field;
  size_t start3[3], count3[3], start4[4], count4[4], in;
  int    x, y, z, c, counter, xp1, yp1, zp1, t, x_ic, y_ic, z_ic;
  int    dimid_lon, dimid_lat, dimid_lev;
  int    varid_tracer0, varid_source;
  int    nlon_INVICAT, nlat_INVICAT, nlev_INVICAT;
  double **INVICAT_field, *lod;
  double ratiox, ratioy;

  printf ("Reading INVICAT state vector file: %s\n", filename);

  // Open the file for reading
  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);

  if (err != 0)
  { printf ("ERROR opening INVICAT file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }

  // The 'times' array is an indicator if structure have already been allocated
  preexist_field = (*field).times;

  // Allocate the state if necessary
  if (!preexist_field)
  { Allocate_state (field,
                    metadata,
                    (*field).horiz_repres,
                    (*field).vert_repres,
                    (*field).temp_repres);
  }


  // Get the dimension lengths of the INVICAT file
  ierr[0]      = nc_inq_dimid  (ncid, "longitude", &dimid_lon);
  ierr[1]      = nc_inq_dimlen (ncid, dimid_lon, &in);
  nlon_INVICAT = in;
  ierr[2]      = nc_inq_dimid  (ncid, "latitude", &dimid_lat);
  ierr[3]      = nc_inq_dimlen (ncid, dimid_lat, &in);
  nlat_INVICAT = in;
  ierr[4]      = nc_inq_dimid  (ncid, "level", &dimid_lev);
  ierr[5]      = nc_inq_dimlen (ncid, dimid_lev, &in);
  nlev_INVICAT = in;

  counter = 0;
  for (c=0; c<6; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting horizontal dimension lengths %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}

  printf ("Number of longitudes in the INVICAT file : %i\n", nlon_INVICAT);
  printf ("          latitudes                      : %i\n", nlat_INVICAT);
  printf ("          levels                         : %i\n", nlev_INVICAT);


  if (nlev_INVICAT < (*field).nlev)
  { printf ("Error: The number of levels in the INVICAT file is insufficient\n");
    printf ("Number of levels required           : %i\n", (*field).nlev);
    exit (0);
  }

  // Get the main variable ids of the data
  if (Read_tracer)
  { ierr[0] = nc_inq_varid (ncid, "prior_concs", &varid_tracer0);
  }
  else
  { ierr[0] = 0;
  }
  if (Read_source)
  { ierr[1] = nc_inq_varid (ncid, "prior_flux",  &varid_source);
  }
  else
  { ierr[1] = 0;
  }

  counter = 0;
  for (c=0; c<2; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting main variable id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}


  // Create a 2D and 1D arrays for reading in INVICAT data
  Array_2d_double_create (&INVICAT_field,
                          nlon_INVICAT,
                          nlat_INVICAT);
  lod = new double[nlon_INVICAT];


  // ===== Deal with the initial condition field =====
  if (Read_tracer)
  {                 count3[0] = 1;              // z
                    count3[1] = 1;              // y
    start3[2] = 0;  count3[2] = nlon_INVICAT;   // x
    for (z=0; z<(*field).nlev; z++)
    { zp1       = z+1;
      z_ic      = nlev_INVICAT - z - 1; // INVICAT is upside-down
      start3[0] = z_ic;
      for (y=0; y<nlat_INVICAT; y++)
      { start3[1] = y;
        err       = nc_get_vara_double (ncid,
                                        varid_tracer0,
                                        start3, count3,
                                        lod);
        // Error reporting
        if (err != 0)
        { printf ("netcdf error inputting tracer from INVICAT file %i: %s\n", err, nc_strerror(err));
          exit(0);
        }
        for (x=0; x<nlon_INVICAT; x++)
        { INVICAT_field[x][y] = lod[x];
        }
      }

      // Put INVICAT_field on the Envi-flux grid for this level
      ratiox = double(nlon_INVICAT) / double((*field).nlon);
      ratioy = double(nlat_INVICAT) / double((*field).nlat);
      for (y=0; y<(*field).nlat; y++)
      { yp1  = y+1;
        y_ic = int(double(y) * ratioy);
        for (x=0; x<(*field).nlon; x++)
        { xp1  = x+1;
          x_ic = int(double(x) * ratiox);
          (*field).tracer0_rs[xp1][yp1][zp1] = INVICAT_field[x_ic][y_ic];
        }
      }
    }
  }

  // ===== Deal with the flux field =====
  if (Read_source)
  { start4[0] = 0;  count4[0] = 1;              // field
    start4[1] = 0;  count4[1] = 1;              // t
                    count4[2] = 1;              // y
    start4[3] = 0;  count4[3] = nlon_INVICAT;   // x
    for (y=0; y<nlat_INVICAT; y++)
    { start4[2] = y;
      err       = nc_get_vara_double (ncid,
                                      varid_source,
                                      start4, count4,
                                      lod);
      // Error reporting
      if (err != 0)
      { printf ("netcdf error inputting source from INVICAT file %i: %s\n", err, nc_strerror(err));
        exit(0);
      }
      for (x=0; x<nlon_INVICAT; x++)
      { INVICAT_field[x][y] = lod[x];
      }
    }

    // Put INVICAT_field on the Envi-flux grid for all times
    ratiox = double(nlon_INVICAT) / double((*field).nlon);
    ratioy = double(nlat_INVICAT) / double((*field).nlat);
    for (y=0; y<(*field).nlat; y++)
    { yp1  = y+1;
      y_ic = int(double(y) * ratioy);
      for (x=0; x<(*field).nlon; x++)
      { xp1  = x+1;
        x_ic = int(double(x) * ratiox);
        for (t=0; t<(*field).nss; t++)
        { (*field).source_rs[x][y][t] = INVICAT_field[x_ic][y_ic];
        }
      }
    }
  }

  // Sort out the halos
  halos (field);

  // Close and tidy up
  err = nc_close (ncid);
  if (err != 0)
  { printf ("ERROR closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

  Array_2d_double_destroy (&INVICAT_field,
                          nlon_INVICAT);
  delete[] lod;
}





// -------------------------------------------------------------------------------
void Read_ecmwf_winds (struct Wind_type *wind,
                       int              t,
                       char             filename[256])
{ // Read-in ECMWF wind data at a single time
  // Assume data structures have already been allocated

  int    ierr[6], err, ncid;
  int    varid_long, varid_lat, varid_lev, varid_u, varid_v, varid_w;
  int    x, y, c, z, counter, xp1, yp1, zp1, xpair;
  size_t start[4], count[4];
  float  plevs_half[::ecmwf_nlev+2], plevs_full[::ecmwf_nlev+2], lod[::ecmwf_nlon];
  float  u_scale_factor[1], u_add_offset[1];
  float  v_scale_factor[1], v_add_offset[1];
  float  w_scale_factor[1], w_add_offset[1];
  short  lodi[::ecmwf_nlon];
  double conversionfactor, dz, dy, dx, dpdz;


  printf ("Reading ecmwf file: %s\n", filename);
  printf ("Time %i\n", t);


  // Open the file for reading
  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);
  if (err != 0)
  { printf ("ERROR opening file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }

  // Get the variable ids of the dimensions
  ierr[0]  = nc_inq_varid (ncid, "longitude",    &varid_long);
  ierr[1]  = nc_inq_varid (ncid, "latitude",     &varid_lat);
  ierr[2]  = nc_inq_varid (ncid, "level",        &varid_lev);
  // Get the variable ids of the data
  ierr[3]  = nc_inq_varid (ncid, "u",            &varid_u);
  ierr[4]  = nc_inq_varid (ncid, "v",            &varid_v);
  ierr[5]  = nc_inq_varid (ncid, "w",            &varid_w);

  counter = 0;
  for (c=0; c<6; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting variable id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}


  // Get attributes
  // --------------------------------------------------------
  ierr[0] = nc_get_att_float (ncid, varid_u, "scale_factor", u_scale_factor);
  ierr[1] = nc_get_att_float (ncid, varid_u, "add_offset",   u_add_offset);
  ierr[2] = nc_get_att_float (ncid, varid_v, "scale_factor", v_scale_factor);
  ierr[3] = nc_get_att_float (ncid, varid_v, "add_offset",   v_add_offset);
  ierr[4] = nc_get_att_float (ncid, varid_w, "scale_factor", w_scale_factor);
  ierr[5] = nc_get_att_float (ncid, varid_w, "add_offset",   w_add_offset);

  //printf ("Scale factors and offsets\n");
  //printf ("u: %f  %f\n", u_scale_factor[0], u_add_offset[0]);
  //printf ("v: %f  %f\n", v_scale_factor[0], v_add_offset[0]);
  //printf ("w: %f  %f\n", w_scale_factor[0], w_add_offset[0]);

  counter = 0;
  for (c=0; c<6; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting attribute id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}


  // Get values of these variables
  // --------------------------------------------------------

  // ===== The dimension variables =====
  start[0] = 0;
  count[0] = ::ecmwf_nlon;
  ierr[0]   = nc_get_vara_float (ncid,
                                 varid_long,
                                 start, count,
                                 lod);
  for (x=0; x<(::ecmwf_nlon); x++)
  { (*wind).longitude_u[x+1] = lod[x];
    (*wind).longitude_v[x+1] = lod[x];
  }
  dx = (*wind).longitude_u[2] - (*wind).longitude_u[1];
  (*wind).longitude_u[0]              = (*wind).longitude_u[1] - dx;
  (*wind).longitude_u[::ecmwf_nlon+1] = (*wind).longitude_u[::ecmwf_nlon] + dx;
  (*wind).longitude_v[0]              = (*wind).longitude_u[0];
  (*wind).longitude_v[::ecmwf_nlon+1] = (*wind).longitude_u[::ecmwf_nlon+1];
  //printf ("Longitudes\n");
  //for (x=0; x<5; x++)
  //{ printf ("%f  ", (*wind).longitude_u[x]);
  //}
  //printf ("... ");
  //for (x=(::ecmwf_nlon-4); x<(::ecmwf_nlon+2); x++)
  //{ printf ("%f  ", (*wind).longitude_u[x]);
  //}
  //printf ("\n\n");

  count[0] = ::ecmwf_nlat;
  ierr[1]  = nc_get_vara_float (ncid,
                                varid_lat,
                                start, count,
                                lod);
  for (y=0; y<(::ecmwf_nlat); y++)
  { (*wind).latitude_u[y+1] = lod[y];
    (*wind).latitude_v[y+1] = lod[y];
  }
  dy = (*wind).latitude_u[1] - (*wind).latitude_u[2];
  (*wind).latitude_u[0]              = (*wind).latitude_u[1] + dy;
  (*wind).latitude_v[0]              = (*wind).latitude_u[0];
  dy = (*wind).latitude_u[::ecmwf_nlat-1] - (*wind).latitude_u[::ecmwf_nlat];
  (*wind).latitude_u[::ecmwf_nlat+1] = (*wind).latitude_u[::ecmwf_nlat] - dy;
  (*wind).latitude_v[::ecmwf_nlat+1] = (*wind).latitude_u[::ecmwf_nlat+1];
  //printf ("Latitudes\n");
  //for (y=0; y<5; y++)
  //{ printf ("%f  ", (*wind).latitude_u[y]);
  //}
  //printf ("... ");
  //for (y=(::ecmwf_nlat-4); y<(::ecmwf_nlat+2); y++)
  //{ printf ("%f  ", (*wind).latitude_u[y]);
  //}
  //printf ("\n\n");


  count[0] = ::ecmwf_nlev;
  ierr[2]  = nc_get_vara_float (ncid,
                                varid_lev,
                                start, count,
                                &(plevs_full[1]));
  // Extrapolate to the halo levels
  dz                         = plevs_full[2] - plevs_full[1];
  plevs_full[0]              = plevs_full[1] - dz;                                       
  dz                         = plevs_full[::ecmwf_nlev] - plevs_full[::ecmwf_nlev-1];
  plevs_full[::ecmwf_nlev+1] = plevs_full[::ecmwf_nlev] + dz;

  // Compute the pressures on half levels
  for (z=0; z<=(::ecmwf_nlev); z++)
  { plevs_half[z] = (plevs_full[z] + plevs_full[z+1]) / 2.0;
  }
  // Extrapolate to the missing halo level
  dz                         = plevs_half[::ecmwf_nlev] - plevs_half[::ecmwf_nlev-1];
  plevs_half[::ecmwf_nlev+1] = plevs_half[1] + dz;


  //printf ("Here are the pressure levels (half on left, full on right)\n");
  //for (z=0; z<=(::ecmwf_nlev+1); z++)
  //{ printf ("%i f               %f\n", z, plevs_full[z]);
  //  printf ("%i h %f\n"              , z, plevs_half[z]);
  //}

  counter = 0;
  for (c=0; c<3; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting dimension variable %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}

  // Convert the pressures (mb) to heights (remember index 0 is for the top of the model)
  // First compute the heights on the full levels
  // (we know the bottom boundary condition for the full levels, zero vertical wind)
  (*wind).level_w[::ecmwf_nlev]  = 0.0;
  // printf ("%i: %f hPa => %f m\n", ::ecmwf_nlev, plevs_full[::ecmwf_nlev], (*wind).level_w[::ecmwf_nlev]);
  // printf ("Full level heights\n");
  for (z=::ecmwf_nlev-1; z>=0; z--)
  { (*wind).level_w[z] = (*wind).level_w[z+1] - 100.0 * (plevs_full[z] - plevs_full[z+1]) *
                                                exp((*wind).level_w[z+1] / ::H) / (::g * ::rho0);
    //printf ("%i: %f hPa => %f m\n", z, plevs_full[z], (*wind).level_w[z]);
  }
  // Complete the bottom halo level
  dz = (*wind).level_w[::ecmwf_nlev] - (*wind).level_w[::ecmwf_nlev-1];
  (*wind).level_w[::ecmwf_nlev+1] = (*wind).level_w[::ecmwf_nlev] + dz;

  // Compute the heights on the half-levels
  // printf ("Half level heights\n");
  for (z=0; z<(::ecmwf_nlev+1); z++)
  { (*wind).level_uv[z] = ((*wind).level_w[z] + (*wind).level_w[z+1]) / 2.0;
    //printf ("%i: %f hPa => %f m\n", z, plevs_half[z], (*wind).level_w[z]);
  }
  // Complete the bottom halo level
  dz = (*wind).level_uv[::ecmwf_nlev] - (*wind).level_uv[::ecmwf_nlev-1];
  (*wind).level_uv[::ecmwf_nlev+1] = (*wind).level_uv[::ecmwf_nlev] + dz;

  //printf ("Here are the vertical levels (half on left, full on right)\n");
  //for (z=0; z<=(::ecmwf_nlev+1); z++)
  //{ printf ("%i f               %f\n", z, (*wind).level_w[z]);
  //  printf ("%i h %f\n"              , z, (*wind).level_uv[z]);
  //}


  // Read-in the winds themselves
  start[3] = 0;  count[3] = ::ecmwf_nlon;
                 count[2] = 1;
                 count[1] = 1;
  start[0] = t;  count[0] = 1;

  // Read-in the zonal winds
  printf ("Reading in zonal winds\n");
  for (z=0; z<(::ecmwf_nlev); z++)
  { start[1] = z;
    for (y=0; y<(::ecmwf_nlat); y++)
    { start[2] = y;
      err      = nc_get_vara_short (ncid,
                                    varid_u,
                                    start, count,
                                    lodi);
      if (err != 0)
      { printf ("ERROR getting u (z=%u, y=%u): %i (%s)\n", z, y, err, nc_strerror(err));
        exit(0);
      }
      // Copy the data to the appropriate place
      for (x=0; x<(::ecmwf_nlon); x++)
      { (*wind).u[x+1][y+1][z+1] = double(lodi[x]) * double(u_scale_factor[0]) +
                                   double(u_add_offset[0]);
      }
    }
  }
  //printf ("Some example zonal winds\n");
  //for (x=10; x<30; x++)
  //{ printf ("%f ", (*wind).u[x][x][x]);
  //}
  //printf ("\n");


  // Read-in the meridional winds
  printf ("Reading in the meridional winds\n");
  for (z=0; z<(::ecmwf_nlev); z++)
  { start[1] = z;
    for (y=0; y<(::ecmwf_nlat); y++)
    { start[2] = y;
      err      = nc_get_vara_short (ncid,
                                    varid_v,
                                    start, count,
                                    lodi);
      if (err != 0)
      { printf ("ERROR getting v (z=%u, y=%u): %i (%s)\n", z, y, err, nc_strerror(err));
        exit(0);
      }
      // Copy the data to the appropriate place
      for (x=0; x<(::ecmwf_nlon); x++)
      { (*wind).v[x+1][y+1][z+1] = double(lodi[x]) * double(v_scale_factor[0]) +
                                   double(v_add_offset[0]);
      }
    }
  }
  //printf ("Some example meridional winds\n");
  //for (x=10; x<30; x++)
  //{ printf ("%f ", (*wind).v[x][x][x]);
  //}
  //printf ("\n");

  // Read-in the vertical winds
  printf ("Reading in the vertical winds\n");
  for (z=0; z<(::ecmwf_nlev); z++)
  { start[1] = z;
    conversionfactor = -1.0 * exp((*wind).level_uv[z] / ::H) / (::rho0 * ::g);
    for (y=0; y<(::ecmwf_nlat); y++)
    { start[2] = y;
      err      = nc_get_vara_short (ncid,
                                    varid_w,
                                    start, count,
                                    lodi);
      if (err != 0)
      { printf ("ERROR getting w (z=%u, y=%u): %i (%s)\n", z, y, err, nc_strerror(err));
        exit(0);
      }
      // Copy the data to the appropriate place
      for (x=0; x<(::ecmwf_nlon); x++)
      { (*wind).w[x+1][y+1][z+1] = (double(lodi[x]) * double(w_scale_factor[0]) +
                                    double(w_add_offset[0])) * conversionfactor;
      }
    }
  }
  //printf ("Some example vertical winds\n");
  //for (x=10; x<30; x++)
  //{ printf ("%f ", (*wind).w[x][x][x]);
  //}
  //printf ("\n");

  // Sort out halo winds
  // (a) in the longitude direction
  for (y=0; y<(::ecmwf_nlat); y++)
  { yp1 = y + 1;
    for (z=0; z<(::ecmwf_nlev); z++)
    { zp1 = z + 1;
      (*wind).u[0][yp1][zp1]             = (*wind).u[::ecmwf_nlon][yp1][zp1];
      (*wind).v[0][yp1][zp1]             = (*wind).v[::ecmwf_nlon][yp1][zp1];
      (*wind).w[0][yp1][zp1]             = (*wind).w[::ecmwf_nlon][yp1][zp1];
      (*wind).u[::ecmwf_nlon+1][yp1][zp1] = (*wind).u[1][yp1][zp1];
      (*wind).v[::ecmwf_nlon+1][yp1][zp1] = (*wind).v[1][yp1][zp1];
      (*wind).w[::ecmwf_nlon+1][yp1][zp1] = (*wind).w[1][yp1][zp1];
    }
  }

  // (b) in the latitude direction
  for (x=0; x<(::ecmwf_nlon+2); x++)
  { xpair = x + ::ecmwf_nlon/2 + 1;
    if (xpair > ::ecmwf_nlon)
    { xpair -= ::ecmwf_nlon;
    }
    for (z=0; z<(::ecmwf_nlev); z++)
    { zp1 = z + 1;
      (*wind).u[x][0][zp1]              = -1.0 * (*wind).u[xpair][2][zp1];
      (*wind).v[x][0][zp1]              = -1.0 * (*wind).v[xpair][2][zp1];
      (*wind).w[x][0][zp1]              = (*wind).w[xpair][2][zp1];
      (*wind).u[x][::ecmwf_nlat+1][zp1] = -1.0 * (*wind).u[xpair][::ecmwf_nlat-1][zp1];
      (*wind).v[x][::ecmwf_nlat+1][zp1] = -1.0 * (*wind).v[xpair][::ecmwf_nlat-1][zp1];
      (*wind).w[x][::ecmwf_nlat+1][zp1] = (*wind).w[xpair][::ecmwf_nlat-1][zp1];
    }
  }

  // (c) in the height direction
  // Use Neuman boundary conditions for u and v, and Dirichlet for w
  for (x=0; x<(::ecmwf_nlon+2); x++)
  { for (y=0; y<(::ecmwf_nlat+2); y++)
    { // Top of model
      (*wind).u[x][y][0] = (*wind).u[x][y][1];
      (*wind).v[x][y][0] = (*wind).v[x][y][1];
      (*wind).w[x][y][0] = -1.0 * (*wind).w[x][y][1];
      // Bottom of model
      (*wind).u[x][y][::ecmwf_nlev+1] = (*wind).u[x][y][::ecmwf_nlev];
      (*wind).v[x][y][::ecmwf_nlev+1] = (*wind).v[x][y][::ecmwf_nlev];
      (*wind).w[x][y][::ecmwf_nlev+1] = -1.0 * (*wind).w[x][y][::ecmwf_nlev];
    }
  }


  err = nc_close (ncid);
  if (err != 0)
  { printf ("ERROR closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
}



// -------------------------------------------------------------------------------
void Write_winds ( struct Wind_type *wind,
                   char             filename[256])

{ // Write u, v, w, wind fields
  //Declare local variables
  int          x, y, z, c;
  // netCDF-related variables
  size_t       nlon1, nlat1, nlev1;
  int          ierr[6], ncid, err;
  int          dimidxu, dimidxv, dimidyu, dimidyv, dimidzuv, dimidzw;
  int          varidxu, varidxv, varidyu, varidyv, varidzuv, varidzw;
  int          varidu, varidv, varidw;
  int          list1[1], list3[3];
  size_t       start[3], count[3];
  double       lineofdata[(*wind).nlon+2];


  printf ("Write_winds: %s\n", filename);

  // Copy over dimension lengths to correct type
  nlon1 = (*wind).nlon+2;
  nlat1 = (*wind).nlat+2;
  nlev1 = (*wind).nlev+2;

  // Create a new netCDF file
  err = nc_create (filename,
                   NC_CLOBBER,
                   &ncid);
  // Error reporting
  if (err != 0)
  { printf ("netcdf error opening output wind file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }

  // Define the dimensions
  ierr[0] = nc_def_dim (ncid,
                        "longitude_u",
                        nlon1,
                        &dimidxu);
  ierr[1] = nc_def_dim (ncid,
                        "longitude_v",
                        nlon1,
                        &dimidxv);
  ierr[2] = nc_def_dim (ncid,
                        "latitude_u",
                        nlat1,
                        &dimidyu);
  ierr[3] = nc_def_dim (ncid,
                        "latitude_v",
                        nlat1,
                        &dimidyv);
  ierr[4] = nc_def_dim (ncid,
                        "level_uv",
                        nlev1,
                        &dimidzuv);
  ierr[5] = nc_def_dim (ncid,
                        "level_w",
                        nlev1,
                        &dimidzw);
  // Error reporting
  for (c=0; c<6; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error defining the dimensions (variable %i) %i: %s\n", c, err, nc_strerror(err));
      exit(0);
    }
  }

  // Define the dimension variables
  list1[0] = dimidxu;
  ierr[0]  = nc_def_var (ncid,
                         "longitude_u",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidxu);
  list1[0] = dimidxv;
  ierr[1]  = nc_def_var (ncid,
                         "longitude_v",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidxv);
  list1[0] = dimidyu;
  ierr[2]  = nc_def_var (ncid,
                         "latitude_u",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidyu);
  list1[0] = dimidyv;
  ierr[3]  = nc_def_var (ncid,
                         "latitude_v",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidyv);
  list1[0] = dimidzuv;
  ierr[4]  = nc_def_var (ncid,
                         "level_uv",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidzuv);
  list1[0] = dimidzw;
  ierr[5]  = nc_def_var (ncid,
                         "level_w",
                         NC_DOUBLE,
                         1,
                         list1,
                         &varidzw);
  // Error reporting
  for (c=0; c<6; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error defining the dimension variables (variable %i) %i: %s\n", c, err, nc_strerror(err));
      exit(0);
    }
  }


  // Define the main variables
  // x, y, z directions correspond to 2, 1, 0 respectively.
  list3[0] = dimidzuv;
  list3[1] = dimidyu;
  list3[2] = dimidxu;
  ierr[0] = nc_def_var (ncid,
                        "u",
                        NC_DOUBLE,
                        3,
                        list3,
                        &varidu);
  list3[1] = dimidyv;
  list3[2] = dimidxv;
  ierr[1] = nc_def_var (ncid,
                        "v",
                        NC_DOUBLE,
                        3,
                        list3,
                        &varidv);
  list3[0] = dimidzw;
  list3[1] = dimidyu;
  list3[2] = dimidxv;
  ierr[2] = nc_def_var (ncid,
                        "w",
                        NC_DOUBLE,
                        3,
                        list3,
                        &varidw);
  // Error reporting
  for (c=0; c<3; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error defining the main variables (variable %i) %i: %s\n", c, err, nc_strerror(err));
      exit(0);
    }
  }


  // End define mode
  err = nc_enddef (ncid);
  // Error reporting
  if (err != 0)
  { printf ("netcdf error changing mode %i: %s\n", err, nc_strerror(err));
    exit(0);
  }


  // Output dimension data
  ierr[0] = nc_put_var_double (ncid,
                               varidxu,
                               (*wind).longitude_u);
  ierr[1] = nc_put_var_double (ncid,
                               varidxv,
                               (*wind).longitude_v);
  ierr[2] = nc_put_var_double (ncid,
                               varidyu,
                               (*wind).latitude_u);
  ierr[3] = nc_put_var_double (ncid,
                               varidyv,
                               (*wind).latitude_v);
  ierr[4] = nc_put_var_double (ncid,
                               varidzuv,
                               (*wind).level_uv);
  ierr[5] = nc_put_var_double (ncid,
                               varidzw,
                               (*wind).level_w);
  // Error reporting
  for (c=0; c<3; c++)
  { err = ierr[c];
    if (err != 0)
    { printf ("netcdf error outputting the dimension variables (variable %i) %i: %s\n", c, err, nc_strerror(err));
      exit(0);
    }
  }

  // x, y, z directions correspond to 2, 1, 0 respectively.
  start[2] = 0;  count[2] = (*wind).nlon+2;
  count[1] = 1;  count[0] = 1;
  // Output u
  for (z=0; z<(*wind).nlev+2; z++)
  { start[0] = z;
    for (y=0; y<(*wind).nlat+2; y++)
    { start[1] = y;
      // Copy part of data into special array
      for (x=0; x<(*wind).nlon+2; x++)
      { lineofdata[x] = (*wind).u[x][y][z];
        err = nc_put_vara_double (ncid,
                                  varidu,
                                  start,
                                  count,
                                  lineofdata);
        // Error reporting
        if (err != 0)
        { printf ("netcdf error outputting u %i: %s\n", err, nc_strerror(err));
          exit(0);
        }
      }
    }
  }

  // Output v
  for (z=0; z<(*wind).nlev+2; z++)
  { start[0] = z;
    for (y=0; y<(*wind).nlat+2; y++)
    { start[1] = y;
      // Copy part of data into special array
      for (x=0; x<(*wind).nlon+2; x++)
      { lineofdata[x] = (*wind).v[x][y][z];
        err = nc_put_vara_double (ncid,
                                  varidv,
                                  start,
                                  count,
                                  lineofdata);
        // Error reporting
        if (err != 0)
        { printf ("netcdf error outputting v %i: %s\n", err, nc_strerror(err));
          exit(0);
        }
      }
    }
  }

  // Output w
  for (z=0; z<(*wind).nlev+2; z++)
  { start[0] = z;
    for (y=0; y<(*wind).nlat+2; y++)
    { start[1] = y;
      // Copy part of data into special array
      for (x=0; x<(*wind).nlon+2; x++)
      { lineofdata[x] = (*wind).w[x][y][z];
        err = nc_put_vara_double (ncid,
                                  varidw,
                                  start,
                                  count,
                                  lineofdata);
        // Error reporting
        if (err != 0)
        { printf ("netcdf error outputting w %i: %s\n", err, nc_strerror(err));
          exit(0);
        }
      }
    }
  }

  //Close the netCDF file
  err = nc_close (ncid);
  // Error reporting
  if (err != 0)
  { printf ("netcdf closing file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }
}


// -------------------------------------------------------------------------------
void Read_winds (struct Wind_type     *wind,
                 struct metadata_type *MetaData,
                 double               factor_w,
                 char                 filename[256])
{ // Read u, v, w, wind fields
  // Assume data structures have already been allocated

  int    ierr[9], err, ncid;
  int    varid_long_u, varid_lat_u, varid_lev_uv;
  int    varid_long_v, varid_lat_v, varid_lev_w;
  int    varid_u, varid_v, varid_w;
  int    x, y, c, z, counter;
  size_t start[3], count[3];
  double lod[(*MetaData).nlon+2];

  printf ("Reading file: %s\n", filename);

  // Open the file for reading
  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);
  if (err != 0)
  { printf ("ERROR opening file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }

  (*wind).nlon = (*MetaData).nlon;
  (*wind).nlat = (*MetaData).nlat;
  (*wind).nlev = (*MetaData).nlev;

  // Get the variable ids of the dimensions
  ierr[0]  = nc_inq_varid (ncid, "longitude_u", &varid_long_u);
  ierr[1]  = nc_inq_varid (ncid, "longitude_v", &varid_long_v);
  ierr[2]  = nc_inq_varid (ncid, "latitude_u",  &varid_lat_u);
  ierr[3]  = nc_inq_varid (ncid, "latitude_v",  &varid_lat_v);
  ierr[4]  = nc_inq_varid (ncid, "level_uv",    &varid_lev_uv);
  ierr[5]  = nc_inq_varid (ncid, "level_w",     &varid_lev_w);
  // Get the variable ids of the data
  ierr[6]  = nc_inq_varid (ncid, "u",           &varid_u);
  ierr[7]  = nc_inq_varid (ncid, "v",           &varid_v);
  ierr[8]  = nc_inq_varid (ncid, "w",           &varid_w);

  counter = 0;
  for (c=0; c<9; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting variable id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}


  // Get values of variables
  // --------------------------------------------------------

  // ===== The dimension variables =====
  start[0] = 0;
  count[0] = (*wind).nlon + 2;
  ierr[0]  = nc_get_vara_double (ncid,
                                 varid_long_u,
                                 start, count,
                                 (*wind).longitude_u);
  ierr[1]  = nc_get_vara_double (ncid,
                                 varid_long_v,
                                 start, count,
                                 (*wind).longitude_v);
  count[0] = (*wind).nlat + 2;
  ierr[2]  = nc_get_vara_double (ncid,
                                 varid_lat_u,
                                 start, count,
                                 (*wind).latitude_u);
  ierr[3]  = nc_get_vara_double (ncid,
                                 varid_lat_v,
                                 start, count,
                                 (*wind).latitude_v);
  count[0] = (*wind).nlev + 2;
  ierr[4]  = nc_get_vara_double (ncid,
                                 varid_lev_uv,
                                 start, count,
                                 (*wind).level_uv);

  ierr[5]  = nc_get_vara_double (ncid,
                                 varid_lev_w,
                                 start, count,
                                 (*wind).level_w);
  counter = 0;
  for (c=0; c<6; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting dimension variable %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}

  // ===== The wind variables =====
                 count[0] = 1;                  // z
                 count[1] = 1;                  // y
  start[2] = 0;  count[2] = (*wind).nlon + 2;   // x
  for (z=0; z<(*wind).nlev+2; z++)
  { start[0] = z;
    for (y=0; y<(*wind).nlat+2; y++)
    { start[1] = y;
      err = nc_get_vara_double (ncid,
                                varid_u,
                                start, count,
                                lod);
      // Error reporting
      if (err != 0)
      { printf ("netcdf error inputting u %i: %s\n", err, nc_strerror(err));
        exit(0);
      }
      for (x=0; x<(*wind).nlon+2; x++)
      { (*wind).u[x][y][z] = lod[x];
      }
    }
  }

  for (z=0; z<(*wind).nlev+2; z++)
  { start[0] = z;
    for (y=0; y<(*wind).nlat+2; y++)
    { start[1] = y;
      err = nc_get_vara_double (ncid,
                                varid_v,
                                start, count,
                                lod);
      // Error reporting
      if (err != 0)
      { printf ("netcdf error inputting v %i: %s\n", err, nc_strerror(err));
        exit(0);
      }
      for (x=0; x<(*wind).nlon+2; x++)
      { (*wind).v[x][y][z] = lod[x];
      }
    }
  }

  for (z=0; z<(*wind).nlev+2; z++)
  { start[0] = z;
    for (y=0; y<(*wind).nlat+2; y++)
    { start[1] = y;
      err = nc_get_vara_double (ncid,
                                varid_w,
                                start, count,
                                lod);
      // Error reporting
      if (err != 0)
      { printf ("netcdf error inputting w %i: %s\n", err, nc_strerror(err));
        exit(0);
      }
      for (x=0; x<(*wind).nlon+2; x++)
      { (*wind).w[x][y][z] = factor_w * lod[x];
      }
    }
  }

  err = nc_close (ncid);
  if (err != 0)
  { printf ("ERROR closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
}


// -------------------------------------------------------------------------------
void WriteTimeSeq
// Write-out a time sequence of a tracer field, one time per call
 ( struct metadata_type       *MetaData,
   struct instant_tracer_type *field,
   int                        ntimes_total, // No of times to be output (inc. 0)
   double                     deltat,       // Separation between times (seconds, but output will be in days)
   int                        timestep,     // Timestep number
   int                        status,       // 0=create file only,
                                            // 1=normal (write only)
   char                       filename[256])


{ //Declare local variables
  int          x, y, z, t, yp1, zp1;
  // netCDF-related variables
  size_t       nlon1, nlat1, nlev1, ntimes_total1;
  double       lineofdata[(*MetaData).nlon], time[ntimes_total];
  int          ierr[5], ncid, err;
  int          dimidx, dimidy, dimidz, dimidt;
  int          varidx, varidy, varidz, varidt, varidtracer;
  int          list[4];
  size_t       start[4], count[4];
  int          ecount, c;



  if (status == 0)
  { // This is the first time this routine is called -- need to create and set-up file
    printf ("WriteTimeSeq: %s\n", filename);
    err = nc_create (filename,
                     NC_CLOBBER,
                     &ncid);
    if (err != 0)
    { printf ("netcdf error opening output file for creation: %i (%s)\n", err, nc_strerror(err));
      exit(0);
    }

    // Copy over dimension lengths to correct type
    nlon1         = (*MetaData).nlon;
    nlat1         = (*MetaData).nlat;
    nlev1         = (*MetaData).nlev;
    ntimes_total1 = ntimes_total;


    // Define the dimensions
    ierr[0] = nc_def_dim (ncid,
                          "longitude",
                          nlon1,
                          &dimidx);
    ierr[1] = nc_def_dim (ncid,
                          "latitude",
                          nlat1,
                          &dimidy);
    ierr[2] = nc_def_dim (ncid,
                          "level",
                          nlev1,
                          &dimidz);
    ierr[3] = nc_def_dim (ncid,
                          "time_days",
                          ntimes_total1,
                          &dimidt);

    ecount = 0;
    for (c=0; c<=3; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error defining dimensions %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }

    // Define the variables
    list[0] = dimidx;
    ierr[0]  = nc_def_var (ncid,
                           "longitude",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidx);
    list[0] = dimidy;
    ierr[1]  = nc_def_var (ncid,
                           "latitude",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidy);
    list[0] = dimidz;
    ierr[2]  = nc_def_var (ncid,
                           "level",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidz);
    list[0] = dimidt;
    ierr[3]  = nc_def_var (ncid,
                           "time_days",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidt);
    list[3] = dimidx;
    list[2] = dimidy;
    list[1] = dimidz;
    list[0] = dimidt;
    ierr[4]  = nc_def_var (ncid,
                           "tracer",
                           NC_DOUBLE,
                           4,
                           list,
                           &varidtracer);

    ecount = 0;
    for (c=0; c<=4; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error defining variable %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }


    // End define mode
    err = nc_enddef (ncid);
    if (err != 0)
    { printf ("netcdf error changing mode: %i (%s)\n", err, nc_strerror(err));
      exit(0);
    }


    // Output dimension data
    ierr[0] = nc_put_var_double (ncid,
                                 varidx,
                                 &((*MetaData).longitude[1]));

    ierr[1] = nc_put_var_double (ncid,
                                 varidy,
                                 &((*MetaData).latitude[1]));

    ierr[2] = nc_put_var_double (ncid,
                                 varidz,
                                 &((*MetaData).level[1]));

    for (t=0; t<ntimes_total; t++)
    { time[t] = double(t) * deltat / 86400.0;  // New units are days
    }

    ierr[3] = nc_put_var_double (ncid,
                                 varidt,
                                 time);

    ecount = 0;
    for (c=0; c<=3; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error outputting dimension data %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }

  }


  else

  { // Output a time slice only -- the file will have already been created
    err = nc_open (filename,
                   NC_WRITE,
                   &ncid);
    if (err != 0)
    { printf ("netcdf error reopening output file: %i (%s)\n", err, nc_strerror(err));
      exit(0);
    }

    // Get the variable id
    err = nc_inq_varid(ncid,
                       "tracer",
                       &varidtracer);
    if (err != 0)
    { printf ("netcdf error enquiring about tracer variable: %i (%s)\n", err, nc_strerror(err));
      exit(0);
    }

    // Output the time slice
    // x, y, z, t directions correspond to 3, 2, 1, 0 respectively.
    start[3] = 0;        count[3] = (*MetaData).nlon;  // x
                         count[2] = 1;                 // y
                         count[1] = 1;                 // z
    start[0] = timestep; count[0] = 1;                 // t

    for (z=0; z<(*MetaData).nlev; z++)
    { zp1      = z + 1;
      start[1] = z;
      for (y=0; y<(*MetaData).nlat; y++)
      { yp1      = y + 1;
        start[2] = y;
        // Copy part of data into special array
        for (x=0; x<(*MetaData).nlon; x++)
        { lineofdata[x] = (*field).tracer[x+1][yp1][zp1];
        }
        err = nc_put_vara_double (ncid,
                                  varidtracer,
                                  start,
                                  count,
                                  lineofdata);
        if (err != 0)
        { printf ("netcdf error outputting tracer: %i (%s)\n", err, nc_strerror(err));
          printf ("  y = %i, z = %i, t = %i\n", y, z, timestep);
          exit(0);
        }
      }
    }
  }

  //Close the netCDF file
  err = nc_close (ncid);
  if (err != 0)
  { printf ("netcdf error closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

}


// -------------------------------------------------------------------------------
void Read1Time
// Read-in a tracer field, at a specified time, from a time sequence
 ( struct metadata_type       *MetaData,       // In  (to check against)
   struct instant_tracer_type *field,          // Out (assumed already allocated)
   int                        time_request,    // In  Time index requested (1 is first)
   char                       filename[256],   // In  Filename
   bool                       *ok )            // Out Successful execution flag


{ //Declare local variables
  int     x, y, z, yp1, zp1;
  // netCDF-related variables
  int     ierr[8], ncid, err;
  int     dimid_lon, dimid_lat, dimid_lev, dimid_t;
  size_t  start[4], count[4], in;
  int     ecount, c;
  int     varid_tracer;
  double  *lod;

  // Open the file for reading
  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);
  if (err != 0)
  { printf ("ERROR opening file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }

  *ok = true;

  // Get, and check the dimension lengths
  ierr[0] = nc_inq_dimid  (ncid, "longitude", &dimid_lon);
  ierr[1] = nc_inq_dimlen (ncid, dimid_lon, &in);
  if ((*MetaData).nlon != in)
  { *ok = false;
    printf ("Error: the number of longitudes is different to that expected\n");
  }

  ierr[2] = nc_inq_dimid  (ncid, "latitude", &dimid_lat);
  ierr[3] = nc_inq_dimlen (ncid, dimid_lat, &in);
  if ((*MetaData).nlat != in)
  { *ok = false;
    printf ("Error: the number of latitudes is different to that expected\n");
  }

  ierr[4] = nc_inq_dimid  (ncid, "level", &dimid_lev);
  ierr[5] = nc_inq_dimlen (ncid, dimid_lev, &in);
  if ((*MetaData).nlev != in)
  { *ok = false;
    printf ("Error: the number of levels is different to that expected\n");
  }

  ierr[6] = nc_inq_dimid  (ncid, "time_days", &dimid_t);
  ierr[7] = nc_inq_dimlen (ncid, dimid_t, &in);
  if (time_request > in)
  { *ok = false;
    printf ("Error: specified time slice is not in the input file\n");
  }


  ecount = 0;
  for (c=0; c<7; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting dimensions id/len %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      ecount++;
    }
  }
  if (ecount > 0) {exit(0);}



  if (*ok)
  { // The dimension lengths are as expected

    // Extract the field
    err = nc_inq_varid (ncid, "tracer", &varid_tracer);
    if (err != 0)
    { printf ("ERROR getting tracer var id %i: %s\n", err, nc_strerror(err));
      exit(0);
    }

    lod = new double[(*MetaData).nlon];
    start[0] = time_request-1;  count[0] = 1;                  // t
                                count[1] = 1;                  // z
                                count[2] = 1;                  // y
    start[3] = 0;               count[3] = (*MetaData).nlon;   // x
    for (z=0; z<(*MetaData).nlev; z++)
    { zp1      = z + 1;
      start[1] = z;
      for (y=0; y<(*MetaData).nlat; y++)
      { yp1      = y + 1;
        start[2] = y;
        err = nc_get_vara_double (ncid,
                                  varid_tracer,
                                  start, count,
                                  lod);
        // Error reporting
        if (err != 0)
        { printf ("netcdf error inputting tracer %i: %s\n", err, nc_strerror(err));
          exit(0);
        }
        for (x=0; x<(*MetaData).nlon; x++)
        { (*field).tracer[x+1][yp1][zp1] = lod[x];
        }
      }
    }
    delete[] lod;

    // Sort halos
    halos_tracer (field,
                  MetaData);

  }

  //Close the netCDF file
  err = nc_close (ncid);
  if (err != 0)
  { printf ("netcdf error closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

}


// -------------------------------------------------------------------------------
void WriteDeparturePoints
// Write-out semi-Lagrangian departure points for a specific time
 ( struct metadata_type *MetaData,
   double               ***dp_lon,
   double               ***dp_lat,
   double               ***dp_lev,
   int                  ***dp_index_lon,
   int                  ***dp_index_lat,
   int                  ***dp_index_lev,
   int                  nmajor,       // No of times to be output
   double               deltat,       // Separation between times (seconds, but output will be in days)
   int                  timestep,     // Timestep number
   int                  status,       // 0=create file only,
                                      // 1=normal (write only)
   char                 filename[256])

{ //Declare local variables
  int          x, y, z, t;
  // netCDF-related variables
  size_t       nlon1, nlat1, nlev1, nmajor1;
  double       lod[(*MetaData).nlon], time[nmajor];
  int          lodi[(*MetaData).nlon];
  int          ierr[10], ncid, err;
  int          dimidx, dimidy, dimidz, dimidt;
  int          varidx, varidy, varidz, varidt, varid_depx, varid_depy, varid_depz;
  int          varid_index_depx, varid_index_depy, varid_index_depz;
  int          list[4];
  size_t       start[4], count[4];
  int          ecount, c;

  if (status == 0)
  { // This is the first time this routine is called -- need to create and set-up file

    printf ("WriteDeparturePoints: %s\n", filename);
    err = nc_create (filename,
                     NC_CLOBBER,
                     &ncid);
    if (err != 0)
    { printf ("netcdf error opening output dp file for creation: %i (%s)\n", err, nc_strerror(err));
      exit(0);
    }


    // Copy over dimension lengths to correct type
    nlon1   = (*MetaData).nlon;
    nlat1   = (*MetaData).nlat;
    nlev1   = (*MetaData).nlev;
    nmajor1 = nmajor;

    // Define the dimensions
    ierr[0] = nc_def_dim (ncid,
                          "longitude",
                          nlon1,
                          &dimidx);
    ierr[1] = nc_def_dim (ncid,
                          "latitude",
                          nlat1,
                          &dimidy);
    ierr[2] = nc_def_dim (ncid,
                          "level",
                          nlev1,
                          &dimidz);
    ierr[3] = nc_def_dim (ncid,
                          "time_days",
                          nmajor1,
                          &dimidt);

    ecount = 0;
    for (c=0; c<=3; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error defining dimensions %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }

    // Define the variables
    list[0] = dimidx;
    ierr[0]  = nc_def_var (ncid,
                           "longitude",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidx);
    list[0] = dimidy;
    ierr[1]  = nc_def_var (ncid,
                           "latitude",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidy);
    list[0] = dimidz;
    ierr[2]  = nc_def_var (ncid,
                           "level",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidz);
    list[0] = dimidt;
    ierr[3]  = nc_def_var (ncid,
                           "time_days",
                           NC_DOUBLE,
                           1,
                           list,
                           &varidt);
    list[3] = dimidx;
    list[2] = dimidy;
    list[1] = dimidz;
    list[0] = dimidt;
    ierr[4]  = nc_def_var (ncid,
                           "dep_lon",
                           NC_DOUBLE,
                           4,
                           list,
                           &varid_depx);
    ierr[5]  = nc_def_var (ncid,
                           "dep_lat",
                           NC_DOUBLE,
                           4,
                           list,
                           &varid_depy);
    ierr[6]  = nc_def_var (ncid,
                           "dep_lev",
                           NC_DOUBLE,
                           4,
                           list,
                           &varid_depz);

    ierr[7]  = nc_def_var (ncid,
                           "dep_index_lon",
                           NC_INT,
                           4,
                           list,
                           &varid_index_depx);
    ierr[8]  = nc_def_var (ncid,
                           "dep_index_lat",
                           NC_INT,
                           4,
                           list,
                           &varid_index_depy);
    ierr[9]  = nc_def_var (ncid,
                           "dep_index_lev",
                           NC_INT,
                           4,
                           list,
                           &varid_index_depz);

    ecount = 0;
    for (c=0; c<=9; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error defining variable %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }


    // End define mode
    err = nc_enddef (ncid);
    if (err != 0)
    { printf ("netcdf error changing mode: %i (%s)\n", err, nc_strerror(err));
      exit(0);
    }


    // Output dimension data
    ierr[0] = nc_put_var_double (ncid,
                                 varidx,
                                 &((*MetaData).longitude[1]));

    ierr[1] = nc_put_var_double (ncid,
                                 varidy,
                                 &((*MetaData).latitude[1]));

    ierr[2] = nc_put_var_double (ncid,
                                 varidz,
                                 &((*MetaData).level[1]));

    for (t=0; t<nmajor; t++)
    { time[t] = double(t) * deltat / 86400.0;  // New units are days
    }

    ierr[3] = nc_put_var_double (ncid,
                                 varidt,
                                 time);

    ecount = 0;
    for (c=0; c<=3; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error outputting dimension data %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }

  }


  else

  { // Output a time slice only -- the file will have already been created
    err = nc_open (filename,
                   NC_WRITE,
                   &ncid);
    if (err != 0)
    { printf ("netcdf error reopening dp output file: %i (%s)\n", err, nc_strerror(err));
      exit(0);
    }

    // Get the variable ids
    ierr[0] = nc_inq_varid(ncid,
                           "dep_lon",
                           &varid_depx);
    ierr[1] = nc_inq_varid(ncid,
                           "dep_lat",
                           &varid_depy);
    ierr[2] = nc_inq_varid(ncid,
                           "dep_lev",
                           &varid_depz);

    ierr[3] = nc_inq_varid(ncid,
                           "dep_index_lon",
                           &varid_index_depx);
    ierr[4] = nc_inq_varid(ncid,
                           "dep_index_lat",
                           &varid_index_depy);
    ierr[5] = nc_inq_varid(ncid,
                           "dep_index_lev",
                           &varid_index_depz);

    ecount = 0;
    for (c=0; c<=5; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error enquiring about dp variable %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }

    // Output the time slice
    // x, y, z, t directions correspond to 3, 2, 1, 0 respectively.
    start[3] = 0;        count[3] = (*MetaData).nlon;  // x
                         count[2] = 1;                 // y
                         count[1] = 1;                 // z
    start[0] = timestep; count[0] = 1;                 // t

    for (z=0; z<(*MetaData).nlev; z++)
    { start[1] = z;
      for (y=0; y<(*MetaData).nlat; y++)
      { start[2] = y;

        // x departure points
        for (x=0; x<(*MetaData).nlon; x++)
        { lod[x] = dp_lon[x][y][z];
        }
        ierr[0] = nc_put_vara_double (ncid,
                                      varid_depx,
                                      start,
                                      count,
                                      lod);

        // y departure points
        for (x=0; x<(*MetaData).nlon; x++)
        { lod[x] = dp_lat[x][y][z];
        }
        ierr[1] = nc_put_vara_double (ncid,
                                      varid_depy,
                                      start,
                                      count,
                                      lod);

        // z departure points
        for (x=0; x<(*MetaData).nlon; x++)
        { lod[x] = dp_lev[x][y][z];
        }
        ierr[2] = nc_put_vara_double (ncid,
                                      varid_depz,
                                      start,
                                      count,
                                      lod);

        // x departure point indices
        for (x=0; x<(*MetaData).nlon; x++)
        { lodi[x] = dp_index_lon[x][y][z];
        }
        ierr[3] = nc_put_vara_int (ncid,
                                   varid_index_depx,
                                   start,
                                   count,
                                   lodi);

        // y departure point indices
        for (x=0; x<(*MetaData).nlon; x++)
        { lodi[x] = dp_index_lat[x][y][z];
        }
        ierr[4] = nc_put_vara_int (ncid,
                                   varid_index_depy,
                                   start,
                                   count,
                                   lodi);

        // z departure point indices
        for (x=0; x<(*MetaData).nlon; x++)
        { lodi[x] = dp_index_lev[x][y][z];
        }
        ierr[5] = nc_put_vara_int (ncid,
                                   varid_index_depz,
                                   start,
                                   count,
                                   lodi);


        ecount = 0;
        for (c=0; c<=5; c++)
        { err = ierr[c];
          if (err != 0)
          { printf ("netcdf error outputting dp variable %u: %i (%s)\n", c, err, nc_strerror(err));
            ecount++;
          }
        }
        if (ecount > 0)
        { exit(0);
        }
      }
    }

  }

  //Close the netCDF file
  err = nc_close (ncid);
  if (err != 0)
  { printf ("netcdf error closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

}


// -------------------------------------------------------------------------------
void ReadDeparturePoints
// Read departure points for a given timestep
 ( struct metadata_type *MetaData,
   double               ***dp_lon,
   double               ***dp_lat,
   double               ***dp_lev,
   int                  ***dp_index_lon,
   int                  ***dp_index_lat,
   int                  ***dp_index_lev,
   int                  timestep,     // Timestep number
   char                 filename[256],
   bool                 *ok)

{ //Declare local variables
  int          x, y, z, t;
  // netCDF-related variables
  double       lod[(*MetaData).nlon];
  int          lodi[(*MetaData).nlon];
  int          ierr[11], ncid, err;
  int          dimidx, dimidy, dimidz, dimidt;
  int          varid_depx, varid_depy, varid_depz, varid_index_depx, varid_index_depy, varid_index_depz;
  size_t       start[4], count[4], in;
  int          ecount, c;

  *ok = true;

  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);
  if (err != 0)
  { printf ("netcdf error opening dp output file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

  // Get, and check the dimension lengths
  ierr[0] = nc_inq_dimid  (ncid, "longitude", &dimidx);
  ierr[1] = nc_inq_dimlen (ncid, dimidx, &in);
  if ((*MetaData).nlon != in)
  { *ok = false;
    printf ("Error: the number of longitudes is different to that expected\n");
  }

  ierr[2] = nc_inq_dimid  (ncid, "latitude", &dimidy);
  ierr[3] = nc_inq_dimlen (ncid, dimidy, &in);
  if ((*MetaData).nlat != in)
  { *ok = false;
    printf ("Error: the number of latitudes is different to that expected\n");
  }

  ierr[4] = nc_inq_dimid  (ncid, "level", &dimidz);
  ierr[5] = nc_inq_dimlen (ncid, dimidz, &in);
  if ((*MetaData).nlev != in)
  { *ok = false;
    printf ("Error: the number of levels is different to that expected\n");
  }

  ierr[6] = nc_inq_dimid  (ncid, "time_days", &dimidt);
  ierr[7] = nc_inq_dimlen (ncid, dimidt, &in);
  if (timestep >= in)
  { *ok = false;
    printf ("Error: specified time slice is not in the input file\n");
  }

  ecount = 0;
  for (c=0; c<8; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting dimensions id/len %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      ecount++;
    }
  }
  if (ecount > 0) {exit(0);}


  if (*ok)
  { // The dimension lengths are as expected


    // Get the variable ids
    ierr[0] = nc_inq_varid(ncid,
                           "dep_lon",
                           &varid_depx);
    ierr[1] = nc_inq_varid(ncid,
                           "dep_lat",
                           &varid_depy);
    ierr[2] = nc_inq_varid(ncid,
                           "dep_lev",
                           &varid_depz);

    ierr[3] = nc_inq_varid(ncid,
                           "dep_index_lon",
                           &varid_index_depx);
    ierr[4] = nc_inq_varid(ncid,
                           "dep_index_lat",
                           &varid_index_depy);
    ierr[5] = nc_inq_varid(ncid,
                           "dep_index_lev",
                           &varid_index_depz);

    ecount = 0;
    for (c=0; c<=5; c++)
    { err = ierr[c];
      if (err != 0)
      { printf ("netcdf error enquiring about dp variable %u: %i (%s)\n", c, err, nc_strerror(err));
        ecount++;
      }
    }
    if (ecount > 0)
    { exit(0);
    }

    // Input the time slice
    // x, y, z, t directions correspond to 3, 2, 1, 0 respectively.
    start[3] = 0;        count[3] = (*MetaData).nlon;  // x
                         count[2] = 1;                 // y
                         count[1] = 1;                 // z
    start[0] = timestep; count[0] = 1;                 // t

    for (z=0; z<(*MetaData).nlev; z++)
    { start[1] = z;
      for (y=0; y<(*MetaData).nlat; y++)
      { start[2] = y;

        // x departure points
        ierr[0] = nc_get_vara_double (ncid,
                                      varid_depx,
                                      start,
                                      count,
                                      lod);
        for (x=0; x<(*MetaData).nlon; x++)
        { dp_lon[x][y][z] = lod[x];
        }

        // y departure points
        ierr[1] = nc_get_vara_double (ncid,
                                      varid_depy,
                                      start,
                                      count,
                                      lod);
        for (x=0; x<(*MetaData).nlon; x++)
        { dp_lat[x][y][z] = lod[x];
        }

        // z departure points
        ierr[2] = nc_get_vara_double (ncid,
                                      varid_depz,
                                      start,
                                      count,
                                      lod);
        for (x=0; x<(*MetaData).nlon; x++)
        { dp_lev[x][y][z] = lod[x];
        }


        // x departure point indices
        ierr[3] = nc_get_vara_int (ncid,
                                   varid_index_depx,
                                   start,
                                   count,
                                   lodi);
        for (x=0; x<(*MetaData).nlon; x++)
        { dp_index_lon[x][y][z] = lodi[x];
        }

        // y departure point indices
        ierr[4] = nc_get_vara_int (ncid,
                                   varid_index_depy,
                                   start,
                                   count,
                                   lodi);
        for (x=0; x<(*MetaData).nlon; x++)
        { dp_index_lat[x][y][z] = lodi[x];
        }

        // z departure point indices
        ierr[5] = nc_get_vara_int (ncid,
                                   varid_index_depz,
                                   start,
                                   count,
                                   lodi);
        for (x=0; x<(*MetaData).nlon; x++)
        { dp_index_lev[x][y][z] = lodi[x];
        }
        ecount = 0;
        for (c=0; c<=5; c++)
        { err = ierr[c];
          if (err != 0)
          { printf ("netcdf error outputting dp variable %u: %i (%s)\n", c, err, nc_strerror(err));
            ecount++;
          }
        }
        if (ecount > 0)
        { exit(0);
        }
      }
    }

  }

  //Close the netCDF file
  err = nc_close (ncid);
  if (err != 0)
  { printf ("netcdf error closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

}


// -------------------------------------------------------------------------------
void WriteObservations ( struct metadata_type  *MetaData,
                         struct obs_trunk_type *obs_trunk,
                         char                  output_obs_file[256] )
// Write-out observations
{ //Declare local variables
  int             lev, c;
  struct obs_type *ob;
  FILE*           obs_file = NULL;


  if (strcmp ("nil", output_obs_file) != 0)
  { // Proceed to output observations file (the reason to output is to include model observation data)
    obs_file = fopen (output_obs_file, "w");
    fprintf (obs_file, "===== Observation file =====\n");
    fprintf (obs_file, "nlevs: %i\n", (*MetaData).nlev);
    fprintf (obs_file, "===== Mass profile =========\n");
    for (lev=1; lev<=(*MetaData).nlev; lev++)
    { fprintf (obs_file, "%03i %lf\n", lev, (*obs_trunk).mass_profile[lev]);
    }

    ob = (*obs_trunk).first;
    c  = 0;
    while (ob)
    { c++;
      fprintf (obs_file, "===== Observation number %u =====\n", c);
      fprintf (obs_file, "ob_of: %c\n", (*ob).ob_of);
      fprintf (obs_file, "ob_type: %c\n", (*ob).obtpe);
      fprintf (obs_file, "time: %i %i %i %i %lf %lf %lf\n", (*ob).time_index, (*ob).day, (*ob).hour, (*ob).min, (*ob).obtime_secs,
                                                         (*ob).time_alpha, (*ob).time_beta);
      fprintf (obs_file, "lon: %i %lf %lf %lf\n", (*ob).lon_index, (*ob).longitude, (*ob).lon_alpha, (*ob).lon_beta);
      fprintf (obs_file, "lat: %i %lf %lf %lf\n", (*ob).lat_index, (*ob).latitude, (*ob).lat_alpha, (*ob).lat_beta);
      if ((*ob).ob_of == 't')
      { fprintf (obs_file, "lev: %i %lf %lf %lf\n", (*ob).lev_index, (*ob).level, (*ob).lev_alpha, (*ob).lev_beta);
      }
      fprintf (obs_file, "ob: % 28.12f % 28.12f\n", (*ob).ob, (*ob).stddev);
      fprintf (obs_file, "model_ob: % 28.12f\n", (*ob).model_ob);
      fprintf (obs_file, "innov: % 28.12f\n", (*ob).innov);
      fprintf (obs_file, "grad: % 28.12f\n", (*ob).dJo_dmodel_ob);
      ob = (*ob).next;
    }

    fclose (obs_file);
  }
  else
  { // Don't output the observations file by request
    printf ("** Do not output observation file by nil request\n");
  }
}


// -------------------------------------------------------------------------------
void ReadObservations ( struct metadata_type  *MetaData,
                        struct obs_trunk_type *obs_trunk,
                        char                  input_obs_file[256] )
// Read-in observations
{ //Declare local variables
  int             lev, c, i1, i2, i3, i4;
  double          d1, d2, d3;
  struct obs_type *ob;
  char            c1;
  char            lod[256];
  FILE*           obs_file = NULL;

  obs_file = fopen (input_obs_file, "r");
  // fprintf (obs_file, "===== Observation file =====\n");
  fgets (lod, 256, obs_file);
  //printf ("LOD1 : %s", lod);
  // fprintf (obs_file, "nlevs: %i\n", (*MetaData).nlev);
  fscanf (obs_file, "%6s %i\n", lod, &i1);
  //printf ("LOD2 : %s\n", lod);
  //printf ("Nlevs : %i\n", i1);
  if (i1 != (*MetaData).nlev)
  { printf ("Attempt to read in incompatible observation file\n");
    printf ("Number of levels specified %i\n", i1);
    printf ("Expecting %i\n", (*MetaData).nlev);
    exit(0);
  }
  else
  { //printf ("Nlevs in obs file has expected value\n");
  }
  // fprintf (obs_file, "===== Mass profile =========\n");
  fgets (lod, 256, obs_file);
  //printf ("LOD3 : %s\n", lod);

  if (!((*obs_trunk).mass_profile))
  { //printf ("Alocating obs_trunk.mass_profile array (%i)\n", (*MetaData).nlev);
    (*obs_trunk).mass_profile = new double[(*MetaData).nlev+1];
  }
  for (lev=1; lev<=(*MetaData).nlev; lev++)
  { // fprintf (obs_file, "%03i %lf\n", lev, (*obs_trunk).mass_profile[lev]);
    fscanf (obs_file, "%03s %lf\n", lod, &d1);
    //printf ("%03i %s >%lf\n", lev, lod, d1);
    (*obs_trunk).mass_profile[lev] = d1;
  }

  c = 0;
  while (!feof(obs_file))
  { if (!((*obs_trunk).first))
    { (*obs_trunk).first = new struct obs_type;
      ob                 = (*obs_trunk).first;
      //printf ("Created first pointer\n");
    }
    else
    { (*ob).next         = new struct obs_type;
      ob                 = (*ob).next;
      //printf ("Created new pointer\n");
    }
    c++;
    // fprintf (obs_file, "===== Observation number %u =====\n", c);
    fgets (lod, 256, obs_file);
    //printf ("LOD: %s\n", lod);

    // fprintf (obs_file, "ob_of: %c\n", (*ob).ob_of);
    fscanf (obs_file, "%6s %c\n", lod, &c1);
    //printf ("LOD (ob of) %s\n", lod);
    //printf ("Ob of:%c\n", c1);
    (*ob).ob_of = c1;

    // fprintf (obs_file, "ob_type: %c\n", (*ob).obtpe);
    fscanf (obs_file, "%8s %c\n", lod, &c1);
    //printf ("Ob type:%c\n", c1);
    (*ob).obtpe = c1;

    // fprintf (obs_file, "time: %i %i %i %i %lf %lf %lf\n", (*ob).time_index, (*ob).day, (*ob).hour, (*ob).min, (*ob).obtime_secs,
    //                                                    (*ob).time_alpha, (*ob).time_beta);
    fscanf (obs_file, "%6s %i %i %i %i %lf %lf %lf\n", lod, &i1, &i2, &i3, &i4, &d1, &d2, &d3);
    //printf ("%s :%i %i %i %i %lf %lf %lf\n", lod, i1, i2, i3, i4, d1, d2, d3);
    (*ob).time_index  = i1;
    (*ob).day         = i2;
    (*ob).hour        = i3;
    (*ob).min         = i4;
    (*ob).obtime_secs = d1;
    (*ob).time_alpha  = d2;
    (*ob).time_beta   = d3;

    // fprintf (obs_file, "lon: %i %lf %lf %lf\n", (*ob).lon_index, (*ob).longitude, (*ob).lon_alpha, (*ob).lon_beta);
    fscanf (obs_file, "%5s %i %lf %lf %lf\n", lod, &i1, &d1, &d2, &d3);
    //printf ("%s :%i %lf %lf %lf\n", lod, i1, d1, d2, d3);
    (*ob).lon_index = i1;
    (*ob).longitude = d1;
    (*ob).lon_alpha = d2;
    (*ob).lon_beta  = d3;

    // fprintf (obs_file, "lat: %i %lf %lf %lf\n", (*ob).lat_index, (*ob).latitude, (*ob).lat_alpha, (*ob).lat_beta);
    fscanf (obs_file, "%5s %i %lf %lf %lf\n", lod, &i1, &d1, &d2, &d3);
    //printf ("%s :%i %lf %lf %lf\n", lod, i1, d1, d2, d3);
    (*ob).lat_index = i1;
    (*ob).latitude  = d1;
    (*ob).lat_alpha = d2;
    (*ob).lat_beta  = d3;

    if ((*ob).ob_of == 't')
    { // fprintf (obs_file, "lev: %i %lf %lf %lf\n", (*ob).lev_index, (*ob).level, (*ob).lev_alpha, (*ob).lev_beta);
      fscanf (obs_file, "%5s %i %lf %lf %lf\n", lod, &i1, &d1, &d2, &d3);
      //printf ("%s :%i %lf %lf %lf\n", lod, i1, d1, d2, d3);
      (*ob).lev_index = i1,
      (*ob).level     = d1;
      (*ob).lev_alpha = d2;
      (*ob).lev_beta  = d3;
    }
    // fprintf (obs_file, "ob: % 28.12f % 28.12f\n", (*ob).ob, (*ob).stddev);
    fscanf (obs_file, "%4s %lf %lf\n", lod, &d1, &d2);
    //printf ("%s : %lf %lf\n", lod, d1, d2);
    (*ob).ob       = d1;
    (*ob).stddev   = d2;
    (*ob).variance = d2 * d2;

    // fprintf (obs_file, "model_ob: % 28.12f\n", (*ob).model_ob);
    fscanf (obs_file, "%10s %lf\n", lod, &d1);
    //printf ("%s : %lf\n", lod, d1);
    (*ob).model_ob = d1;

    // fprintf (obs_file, "innov: % 28.12f\n", (*ob).innov);
    fscanf (obs_file, "%7s %lf\n", lod, &d1);
    //printf ("%s : %lf\n", lod, d1);
    (*ob).innov = d1;

    //fprintf (obs_file, "grad: % 28.12f\n", (*ob).dJo_dmodel_ob);
    fscanf (obs_file, "%6s %lf\n", lod, &d1);
    //printf ("%s : %lf\n", lod, d1);
    (*ob).dJo_dmodel_ob = d1;
  }

  fclose (obs_file);

  printf ("Read-in %i observations\n", c);
}



// -------------------------------------------------------------------------------
void Read_INVICAT_all (struct state_type *state,
                       char              filename[256])
{ // Read-in INVICAT data

  int                  ierr[4], err, ncid;
  int                  varid_long, varid_lat, varid_field, varid_month, varid_flux, varid_tracer;
  int                  dimid_long, dimid_lat, dimid_month, dimid_lev;
  int                  c, counter, x, y, z, t;
  size_t               n, start[4], count[4];
  struct metadata_type MetaData_INVICAT;
  float                *lod;
  double               dx, dy;

  printf ("Reading INCIVAT flux file: %s\n", filename);


  // Open the file for reading
  err = nc_open (filename,
                 NC_NOWRITE,
                 &ncid);
  if (err != 0)
  { printf ("ERROR opening file %i: %s\n", err, nc_strerror(err));
    exit(0);
  }

  // Get the variable ids of the dimensions
  ierr[0]  = nc_inq_varid (ncid, "longitude",    &varid_long);
  ierr[1]  = nc_inq_varid (ncid, "latitude",     &varid_lat);
  // Get the variable id of the data
  ierr[2]  = nc_inq_varid (ncid, "prior_flux",   &varid_flux);
  ierr[3]  = nc_inq_varid (ncid, "prior_concs",  &varid_tracer);

  counter = 0;
  for (c=0; c<4; c++)
  { if (ierr[c] != 0)
    { printf ("ERROR getting variable id %u: %i (%s)\n", c, ierr[c], nc_strerror(ierr[c]));
      counter++;
    }
  }
  if (counter > 0) {exit(0);}

  // Get the dimension sizes
  ierr[0] = nc_inq_dimid  (ncid, "longitude", &dimid_long);
  ierr[1] = nc_inq_dimlen (ncid, dimid_long, &n);
  MetaData_INVICAT.nlon = n;
  ierr[0] = nc_inq_dimid  (ncid, "latitude", &dimid_lat);
  ierr[1] = nc_inq_dimlen (ncid, dimid_lat, &n);
  MetaData_INVICAT.nlat = n;
  ierr[0] = nc_inq_dimid  (ncid, "level", &dimid_lev);
  ierr[1] = nc_inq_dimlen (ncid, dimid_lev, &n);
  MetaData_INVICAT.nlev = n;
  ierr[0] = nc_inq_dimid  (ncid, "month", &dimid_month);
  ierr[1] = nc_inq_dimlen (ncid, dimid_month, &n);
  MetaData_INVICAT.nss          = n;
  MetaData_INVICAT.ntimes_major = n;

  printf ("No. of longs  %i\n", MetaData_INVICAT.nlon);
  printf ("No. of lats   %i\n", MetaData_INVICAT.nlat);
  printf ("No. of levs   %i\n", MetaData_INVICAT.nlev);
  printf ("No. of months %i\n", MetaData_INVICAT.nss);

  // Allocate metadata structure of the INVICAT state
  Allocate_metadata (&MetaData_INVICAT,
                     MetaData_INVICAT.nlon,
                     MetaData_INVICAT.nlat,
                     MetaData_INVICAT.nlev,
                     MetaData_INVICAT.nss);

  // Make a temporary array
  lod = new float [MetaData_INVICAT.nlon];


  // Get values of dimension variables
  start[0] = 0;
  count[0] = MetaData_INVICAT.nlon;
  err      = nc_get_vara_float (ncid,
                                varid_long,
                                start, count,
                                lod);
  if (err != 0)
  { printf ("ERROR getting longitudes from INVICAT file %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
  for (x=0; x<MetaData_INVICAT.nlon; x++)
  { MetaData_INVICAT.longitude[x+1] = lod[x]; 
  }
  dx = MetaData_INVICAT.longitude[2] - MetaData_INVICAT.longitude[1];
  MetaData_INVICAT.longitude[0]                       = MetaData_INVICAT.longitude[1] - dx;
  MetaData_INVICAT.longitude[MetaData_INVICAT.nlon+1] = MetaData_INVICAT.longitude[MetaData_INVICAT.nlon] + dx;


  count[0] = MetaData_INVICAT.nlat;
  err      = nc_get_vara_float (ncid,
                                varid_lat,
                                start, count,
                                lod);
  if (err != 0)
  { printf ("ERROR getting latitudes from INVICAT file %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }
  for (y=0; y<MetaData_INVICAT.nlat; y++)
  { MetaData_INVICAT.latitude[y+1] = lod[y];
  }
  dy = MetaData_INVICAT.latitude[2] - MetaData_INVICAT.latitude[1];
  MetaData_INVICAT.latitude[0]                       = MetaData_INVICAT.latitude[1] - dy;
  MetaData_INVICAT.latitude[MetaData_INVICAT.nlat+1] = MetaData_INVICAT.latitude[MetaData_INVICAT.nlat] + dy;


  for (t=0; t<MetaData_INVICAT.nss; t++)
  { MetaData_INVICAT.times[y] = double(t);
  }


  // Declare the state
  Allocate_state (state,
                  &MetaData_INVICAT,
                  'r', 'r', 'r');


  // Read-in the flux fields
  start[3] = 0;  count[3] = MetaData_INVICAT.nlon;  // Longitude
                 count[2] = 1;                      // Latitude
                 count[1] = 1;                      // Month
  start[0] = 0;  count[0] = 1;                      // Field

  printf ("Reading in INVICAT flux\n");
  for (t=0; t<MetaData_INVICAT.nss; t++)
  { start[1] = t;
    for (y=0; y<MetaData_INVICAT.nlat; y++)
    { start[2] = y;
      err      = nc_get_vara_float (ncid,
                                    varid_flux,
                                    start, count,
                                    lod);
      if (err != 0)
      { printf ("ERROR getting flux (t=%u, y=%u): %i (%s)\n", t, y, err, nc_strerror(err));
        exit(0);
      }
      // Copy the data to the appropriate place
      for (x=0; x<MetaData_INVICAT.nlon; x++)
      { (*state).source_rs[x+1][y+1][t] = double(lod[x]);
      }
    }
  }



  // Read-in the tracer field
  start[2] = 0;  count[2] = MetaData_INVICAT.nlon;  // Longitude
                 count[1] = 1;                      // Latitude
                 count[0] = 1;                      // Level

  printf ("Reading in INVICAT initial tracer\n");
  for (z=0; z<MetaData_INVICAT.nlev; z++)
  { start[0] = z;
    for (y=0; y<MetaData_INVICAT.nlat; y++)
    { start[1] = y;
      err      = nc_get_vara_float (ncid,
                                    varid_tracer,
                                    start, count,
                                    lod);
      if (err != 0)
      { printf ("ERROR getting tracer (z=%u, y=%u): %i (%s)\n", z, y, err, nc_strerror(err));
        exit(0);
      }
      // Copy the data to the appropriate place (remember that the levels need to be reversed)
      for (x=0; x<MetaData_INVICAT.nlon; x++)
      { (*state).tracer0_rs[x+1][y+1][MetaData_INVICAT.nlev-z] = double(lod[x]);
      }
    }
  }


  err = nc_close (ncid);
  if (err != 0)
  { printf ("ERROR closing file: %i (%s)\n", err, nc_strerror(err));
    exit(0);
  }

  // Sort out halos
  halos (state);


  // Tidy up
  delete[] lod;
  Deallocate_metadata (&MetaData_INVICAT);

}
