#include <source.h>

// -------------------------------------------------------------------------------
int MakeWinds_arguments (
       int           argument_count,       // in
       char          **argument_list,      // in
       char          CVTfilename[256],     // out CVT filename for meta data
       char          wind_file_list[256],  // out
       char          red_res_dir[256] )    // out

{ // Declare local variables
  int    success, count;
  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;

  if (argument_count > count)
  { strcpy (CVTfilename, argument_list[count]);
    printf ("Filename of CVT for meta data: %s\n", CVTfilename);
  }
  else
  { success = 0;
    printf ("Please include the CVT filename (for meta data)\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (wind_file_list, argument_list[count]);
    printf ("Filename containing full-resolution ECMWF winds: %s\n", wind_file_list);
  }
  else
  { success = 0;
    printf ("Please include the filename containing full-resolution ECMWF winds\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (red_res_dir, argument_list[count]);
    printf ("Output directory for reduced resolution winds: %s\n", red_res_dir);
  }
  else
  { success = 0;
    printf ("Please include output directory for the reduced resolution winds\n");
  }

  return success;
}



// -------------------------------------------------------------------------------
int MakeFields_arguments (
       int              argument_count,       // in
       char             **argument_list,      // in
       char             CVTfilename[256],     // out CVT filename for meta data
       char             init_cond_file[256],  // out init cond filename
       char             *init_cond_file_type, // out 'e'nvi-flux or 'i'nvicat
       double           *factor_ic,           // out Multiplication factor for initial condition
       double           *factor_flux,         // out Multiplication factor for flux
       int              *num_init_cond_blobs, // out number of init cond blobs
       struct blob_type **ic_blob1,           // out linked list of blobs for init conds
       int              *num_ss,              // out number of surface flux fields
       double           *sDt,                 // out Source timestep (seconds, converted from days)
       int              **num_ss_blobs,       // out number of surface flux blobs in each field
       struct blob_type ***ss_field1,         // out array of linked lists of blobs for init conds
       double           *min_flux_set,        // out Sets min fluxes to this value (preserve sign)
       char             output_file[256] )    // out output filename

{ // Declare local variables
  int              success, count, blob, ss;
  struct blob_type *current_blob, *current_ss;
  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;

  if (argument_count > count)
  { strcpy (CVTfilename, argument_list[count]);
    printf ("Filename of CVT for meta data: %s\n", CVTfilename);
  }
  else
  { success = 0;
    printf ("Please include the CVT filename (for meta data)\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (init_cond_file, argument_list[count]);
    printf ("Filename containing intitial condition and fluxes: %s\n", init_cond_file);
  }
  else
  { success = 0;
    printf ("Please include the filename containing the initial conditions and fluxes (or nil)\n");
  }

  if (strcmp ("nil", init_cond_file) != 0)
  { count++;
    if (argument_count > count)
    { *init_cond_file_type = argument_list[count][0];
      printf ("Initial condition file type: %c\n", *init_cond_file_type);
      if (*init_cond_file_type != 'i' && *init_cond_file_type != 'e')
      { printf ("Must enter i or e (invicat or envi-flux)\n");
        success = 0;
      }
    }
    else
    { success = 0;
      printf ("Please include initial condition file type (i)nvicat or (e)nvi-flux\n");
    }

    count++;
    if (argument_count > count)
    { *factor_ic = atof(argument_list[count]);
      printf ("Factor for the initial condition read in: %f\n", *factor_ic);
    }
    else
    { success = 0;
      printf ("Please include the multiplication factor for the initial condition that is read in\n");
    }

    count++;
    if (argument_count > count)
    { *factor_flux = atof(argument_list[count]);
      printf ("Factor for the flux that is read in: %f\n", *factor_flux);
    }
    else
    { success = 0;
      printf ("Please include the multiplication factor for the flux that is read in\n");
    }
  }

  count++;
  if (argument_count > count)
  { *num_init_cond_blobs = atoi(argument_list[count]);
    printf ("Number of blobs to be added to initial condition file: %i\n", *num_init_cond_blobs);
  }
  else
  { success = 0;
    printf ("Please include number of blobs to be added to initial condition file\n");
  }

  if ((success > 0) && (*num_init_cond_blobs > 0))
  { for (blob=0; blob<*num_init_cond_blobs; blob++)
    { if (blob == 0)
      { // Set the root of the linked list
        *ic_blob1    = new struct blob_type;
        current_blob = *ic_blob1;
      }
      else
      { (*current_blob).next = new struct blob_type;
        current_blob = (*current_blob).next;
      }
      // Set characteristics of this initial condition blob

      count++;
      if (argument_count > count)
      { (*current_blob).longitude = atof(argument_list[count]);
        printf ("IC blob %i longitude: %f deg\n", blob+1, (*current_blob).longitude);
      }
      else
      { success = 0;
        printf ("Please include the longitude for IC blob %i (deg)\n", blob+1);
      }

      count++;
      if (argument_count > count)
      { (*current_blob).latitude = atof(argument_list[count]);
        printf ("IC blob %i latitude: %f deg\n", blob+1, (*current_blob).latitude);
      }
      else
      { success = 0;
        printf ("Please include the latitude for IC blob %i (deg)\n", blob+1);
      }

      count++;
      if (argument_count > count)
      { (*current_blob).level = atof(argument_list[count]);
        printf ("IC blob %i level: %f m\n", blob+1, (*current_blob).level);
      }
      else
      { success = 0;
        printf ("Please include the level for IC blob %i (m)\n", blob+1);
      }

      count++;
      if (argument_count > count)
      { (*current_blob).amplitude = atof(argument_list[count]);
        printf ("IC blob %i amplitude: %f tracer units\n", blob+1, (*current_blob).amplitude);
      }
      else
      { success = 0;
        printf ("Please include the amplitude for IC blob %i (tracer units)\n", blob+1);
      }

      count++;
      if (argument_count > count)
      { (*current_blob).size_h = atof(argument_list[count]);
        printf ("IC blob %i horizontal size: %f deg\n", blob+1, (*current_blob).size_h);
      }
      else
      { success = 0;
        printf ("Please include the horizontal size for IC blob %i (deg)\n", blob+1);
      }

      count++;
      if (argument_count > count)
      { (*current_blob).size_v = atof(argument_list[count]);
        printf ("IC blob %i vertical size: %f m\n", blob+1, (*current_blob).size_v);
      }
      else
      { success = 0;
        printf ("Please include the vertical size for IC blob %i (m)\n", blob+1);
      }

      (*current_blob).next = NULL;
    }
  }

  count++;
  if (argument_count > count)
  { *num_ss = atoi(argument_list[count]);
    printf ("Number of source/sink fields: %i\n", *num_ss);
  }
  else
  { success    = 0;
    *num_ss    = 0;
    printf ("Please include number of source/sink fields\n");
  }

  if ((success > 0) && (*num_ss > 0))
  { count++;
    if (argument_count > count)
    { *sDt = atof(argument_list[count]);
      printf ("Time step between between source/sink files (days): %f\n", *sDt);
      *sDt *= 3600.0 * 24.0;
      printf ("                                            (secs): %f\n", *sDt);
    }
    else
    { success = 0;
      printf ("Please include time step between source/sink files (days)\n");
    }


    // Make array of pointers to linked lists of blobs
    *ss_field1 = new struct blob_type*[*num_ss];
    // Set each to null
    for (ss=0; ss<*num_ss; ss++)
    { (*ss_field1)[ss] = NULL;
    }

    // Make array of numbers of blobs in each source/sink field
    *num_ss_blobs = new int[*num_ss];

    // Loop over source/sink fields
    for (ss=0; ss<*num_ss; ss++)
    { count++;
      if (argument_count > count)
      { (*num_ss_blobs)[ss] = atoi(argument_list[count]);
        printf ("Number of blobs for SS %i: %i\n", ss, (*num_ss_blobs)[ss]);
      }
      else
      { success = 0;
       printf ("Please include number of blobs for source/sink field %i\n", ss);
      }

      if ((success > 0) && ((*num_ss_blobs)[ss] > 0))
      { for (blob=0; blob<(*num_ss_blobs)[ss]; blob++)
        { if (blob == 0)
          { // Set the root of the linked list
            (*ss_field1)[ss] = new struct blob_type;
            current_blob     = (*ss_field1)[ss];
          }
          else
          { (*current_blob).next = new struct blob_type;
            current_blob = (*current_blob).next;
          }
          // Set characteristics of this source/sink blob

          count++;
          if (argument_count > count)
          { (*current_blob).longitude = atof(argument_list[count]);
            printf ("SS %i blob %i longitude: %f deg\n", ss+1, blob+1, (*current_blob).longitude);
          }
          else
          { success = 0;
            printf ("Please include the longitude for SS %i IC blob %i (deg)\n", ss+1, blob+1);
          }

          count++;
          if (argument_count > count)
          { (*current_blob).latitude = atof(argument_list[count]);
            printf ("SS %i blob %i latitude: %f deg\n", ss+1, blob+1, (*current_blob).latitude);
          }
          else
          { success = 0;
            printf ("Please include the latitude for SS %i IC blob %i (deg)\n", ss+1, blob+1);
          }

          (*current_blob).level = 0.0;

          count++;
          if (argument_count > count)
          { (*current_blob).amplitude = atof(argument_list[count]);
            printf ("SS %i blob %i amplitude: %f tracer units/s\n", ss+1, blob+1, (*current_blob).amplitude);
          }
          else
          { success = 0;
            printf ("Please include the amplitude for SS %i IC blob %i (tracer units/s)\n", ss+1, blob+1);
          }

          count++;
          if (argument_count > count)
          { (*current_blob).size_h = atof(argument_list[count]);
            printf ("SS %i blob %i horizontal size: %f deg\n", ss+1, blob+1, (*current_blob).size_h);
          }
          else
          { success = 0;
            printf ("Please include the horizontal size for SS %i IC blob %i (deg)\n", ss+1, blob+1);
          }

          (*current_blob).size_v = 0.0;
          (*current_blob).next = NULL;
        }
      }
    }
  }

  count++;
  if (argument_count > count)
  { *min_flux_set = atof(argument_list[count]);
    printf ("Min flux set: %f\n", *min_flux_set);
  }
  else
  { success = 0;
    printf ("Please include the min flux allowed in final file (smaller values will be set to this, sign preserving)\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (output_file, argument_list[count]);
    printf ("Output file: %s\n", output_file);
  }
  else
  { success = 0;
    printf ("Please include the output filename\n");
  }

  return success;
}


// -------------------------------------------------------------------------------
void Destroy_MakeFields_arguments (
       struct blob_type **ic_blob1,           // out linked list of blobs for init conds
       int              num_ss,               // in  number of surface flux fields
       int              **num_ss_blobs,       // out number of surface flux blobs in each field
       struct blob_type ***ss_field1 )        // out array of linked lists of blobs for init conds
{ // Delete the structure created inside MakeFields_arguments
  int              blob, ss;
  struct blob_type *current_blob, *next_blob, *current_ss;

  if (*ic_blob1)
  { // Delete the ic blobs
    current_blob = *ic_blob1;
    while (current_blob)
    { next_blob = (*current_blob).next;
      delete current_blob;
      current_blob = next_blob;
    }
  }

  // Delete the ss infos
  for (ss=0; ss<num_ss; ss++)
  { if ((*ss_field1)[ss])
    { current_blob = (*ss_field1)[ss];
      // Delete the ss blobs
      while (current_blob)
      { next_blob = (*current_blob).next;
        delete current_blob;
        current_blob = next_blob;
      }
    }
  }

  // Delete the ss_fields
  if (*ss_field1)
  { delete [] (*ss_field1);
  }

  // Delete the num_ss_blobs
  if (*num_ss_blobs)
  { delete [] (*num_ss_blobs);
  }

}


// -------------------------------------------------------------------------------
int RunTraj_arguments (
       int     argument_count,         // in
       char    **argument_list,        // in
       double  *Dt,                    // out Major timestep (between wind files, seconds)
       double  *dt,                    // out Minor timestep (for integration scheme, seconds)
       double  *w_factor,              // out Multiplication factor of vertical winds
       double  *runlength,             // out Run length (days)
       char    wind_dir[256],          // out Wind directory
       char    output_file[256])       // out Output filename

{ // Declare local variables
  int    success, count, truefalse;
  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 0;

  count++;
  if (argument_count > count)
  { strcpy (wind_dir, argument_list[count]);
    printf ("Directory containing the winds: %s\n", wind_dir);
  }
  else
  { success = 0;
    printf ("Please include directory containing the winds\n");
  }

  count++;
  if (argument_count > count)
  { *Dt = atof(argument_list[count]);
    printf ("Time step between wind files (major time step, seconds): %f\n", *Dt);
  }
  else
  { success = 0;
    printf ("Please include the time step between wind files (major time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *dt = atof(argument_list[count]);
    printf ("Time step of integration scheme (minor time step, seconds): %f\n", *dt);
  }
  else
  { success = 0;
    printf ("Please include the time step of integration scheme (minor time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *w_factor = atof(argument_list[count]);
    printf ("Multiplication factor of vertical winds: %f\n", *w_factor);
  }
  else
  { success = 0;
    printf ("Please include the multiplication factor of vertical winds\n");
  }

  count++;
  if (argument_count > count)
  { *runlength = atof(argument_list[count]);
    printf ("Run length (in days): %f\n", *runlength);
  }
  else
  { success = 0;
    printf ("Please include the run length (in days)\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (output_file, argument_list[count]);
    printf ("Output file (animation): %s\n", output_file);
  }
  else
  { success = 0;
    printf ("Please include output filename\n");
  }
  return success;
}


// -------------------------------------------------------------------------------
int SemiLagrangian_arguments (
       int     argument_count,         // in
       char    **argument_list,        // in
       char    state_ic[256],          // out Filename of state (initial condition and source)
       char    wind_dir[256],          // out Directory containing driver winds
       double  *Dt,                    // out Major timestep (between wind files, seconds)
       double  *dt,                    // out Minor timestep (for integration scheme, seconds)
       double  *kappa_dt,              // out Timestep for diffusion
       double  *kappa_h,               // out Horizontal diffusion coefficient
       double  *kappa_v,               // out Vertical diffusion coefficient
       bool    *inc_vert,              // out Include vertical transport?
       double  *factor_w,              // out Multiplication factor of vertical winds
       double  *runlength,             // out Run length (seconds, converted from days)
       double  *output_freq,           // out Time separation between outputs (seconds)
       char    *interpolation_lc,      // out (l)inear or (c)ubic interpolation
       char    output_file_anim[256],  // out Output filename
       char    output_file_diags[256], // out Output filename (diagnostics)
       char    file_dps[256],          // out File for departure points
       char    *dp_file_rw )           // out (r)ead or (w) departure points for forward model

{ // Declare local variables
  int success, count, truefalse;
  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 0;

  count++;
  if (argument_count > count)
  { strcpy (state_ic, argument_list[count]);
    printf ("Filename of the initial state and source: %s\n", state_ic);
  }
  else
  { success = 0;
    printf ("Please include the filename of the initial state and source\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (wind_dir, argument_list[count]);
    printf ("Directory containing the winds: %s\n", wind_dir);
  }
  else
  { success = 0;
    printf ("Please include directory containing the winds\n");
  }

  count++;
  if (argument_count > count)
  { *Dt = atof(argument_list[count]);
    printf ("Time step between wind files (major time step, seconds): %f\n", *Dt);
  }
  else
  { success = 0;
    printf ("Please include the time step between wind files (major time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *dt = atof(argument_list[count]);
    printf ("Time step of integration scheme (minor time step, seconds): %f\n", *dt);
  }
  else
  { success = 0;
    printf ("Please include the time step of integration scheme (minor time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_dt = atof(argument_list[count]);
    printf ("Time step of diffusion scheme (seconds): %f\n", *kappa_dt);
  }
  else
  { success = 0;
    printf ("Please include the time step of diffusion scheme (seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_h = atof(argument_list[count]);
    printf ("Horizontal diffusion coefficient: %f\n", *kappa_h);
  }
  else
  { success = 0;
    printf ("Please include the horizontal diffusion coefficient\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_v = atof(argument_list[count]);
    printf ("Vertical diffusion coefficient: %f\n", *kappa_v);
  }
  else
  { success = 0;
    printf ("Please include the vertical diffusion coefficient\n");
  }

  count++;
  if (argument_count > count)
  { truefalse = atoi(argument_list[count]);
    *inc_vert = truefalse == 1;
    if (*inc_vert)
    { printf ("Include vertical transport\n");
    }
    else
    { printf ("Do not include vertical transport\n");
      printf ("***** NOTE: THIS OPTION SHOULD BE ACCOMPANIED BY factor_w = 0.0  *****\n");
    }
  }
  else
  { success = 0;
    printf ("Please include whether to include vertical transport (0/1)\n");
  }

  count++;
  if (argument_count > count)
  { *factor_w = atof(argument_list[count]);
    printf ("Multiplication factor of vertical winds: %f\n", *factor_w);
  }
  else
  { success = 0;
    printf ("Please include the multiplication factor of vertical winds\n");
  }

  count++;
  if (argument_count > count)
  { *runlength = atof(argument_list[count]);
    printf ("Run length (days): %f\n", *runlength);
    *runlength *= 3600.0 * 24.0;
    printf ("           (secs): %f\n", *runlength);
  }
  else
  { success = 0;
    printf ("Please include the run length (in days)\n");
  }

  count++;
  if (argument_count > count)
  { *output_freq = atof(argument_list[count]);
    printf ("Time separation between outputs (seconds): %f\n", *output_freq);
  }
  else
  { success = 0;
    printf ("Please include the time separation between outputs (seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *interpolation_lc = argument_list[count][0];
    printf ("Linear or cubic interpolation: %c\n", *interpolation_lc);
    if (*interpolation_lc != 'l' && *interpolation_lc != 'c')
    { printf ("Must enter l or c\n");
      success = 0;
    }
  }
  else
  { success = 0;
    printf ("Please include interpolation type (l=linear, c=cubic)\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (output_file_anim, argument_list[count]);
    printf ("Output file (animation): %s\n", output_file_anim);
  }
  else
  { success = 0;
    printf ("Please include output filename for animation\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (output_file_diags, argument_list[count]);
    printf ("Output file (diagnostics): %s\n", output_file_diags);
  }
  else
  { success = 0;
    printf ("Please include output filename for diagnostics\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (file_dps, argument_list[count]);
    printf ("File (semi-Lagrangian departure points): %s\n", file_dps);
  }
  else
  { success = 0;
    printf ("Please include filename for semi-Lagrangian departure points (or nil)\n");
  }

  count++;
  if (argument_count > count)
  { *dp_file_rw = argument_list[count][0];
    printf ("Read or write departure points file: %c\n", *dp_file_rw);
    if (*dp_file_rw != 'r' && *dp_file_rw != 'w')
    { printf ("Must enter r or w\n");
      success = 0;
    }
  }
  else
  { success = 0;
    printf ("Please include whether to (r)ead or (w)rite departure points file\n");
  }


  return success;
}


// -------------------------------------------------------------------------------
int ReplaceFirstField_arguments (
       int              argument_count,            // in
       char             **argument_list,           // in
       char             State_in_filename[256],    // out Input state filename
       char             Forecast_in_filename[256], // out Input forecast filename
       int              *FcTime,                   // out Chosen time index to be extracted from fc (1 is first)
       char             State_out_filename[256] )  // out Output state filename

{ // Declare local variables
  int              success, count;

  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;

  if (argument_count > count)
  { strcpy (State_in_filename, argument_list[count]);
    printf ("Filename of input state: %s\n", State_in_filename);
  }
  else
  { success = 0;
    printf ("Please include the input state filename (ics + fluxes)\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (Forecast_in_filename, argument_list[count]);
    printf ("Filename of input forecast: %s\n", Forecast_in_filename);
  }
  else
  { success = 0;
    printf ("Please include the input forecast filename\n");
  }

  count++;
  if (argument_count > count)
  { *FcTime = atoi(argument_list[count]);
    printf ("Index of chosen time to be extracted from forcast (1 is first): %i\n", *FcTime);
  }
  else
  { success = 0;
    printf ("Please include index of chosen time to be extracted from forcast (1 is first)\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (State_out_filename, argument_list[count]);
    printf ("Filename of output state: %s\n", State_out_filename);
  }
  else
  { success = 0;
    printf ("Please include the output state filename\n");
  }

  return success;
}


// -------------------------------------------------------------------------------
int GenerateObs_arguments (
       int             argument_count,    // in
       char            **argument_list,   // in
       char            state_ic[256],     // out Filename of state (initial condition and source)
       char            wind_dir[256],     // out Directory containing driver winds
       char            obs_file[256],     // out Output file of observations
       char            file_dps[256],     // out File for departure points
       double          *Dt,               // out Major timestep (between wind files, seconds)
       double          *dt,               // out Minor timestep (for integration scheme, seconds)
       double          *kappa_dt,         // out Timestep for diffusion
       double          *kappa_h,          // out Horizontal diffusion coefficient
       double          *kappa_v,          // out Vertical diffusion coefficient
       bool            *inc_vert,         // out Include vertical transport?
       double          *factor_w,         // out Multiplication factor of vertical winds
       double          *output_freq,      // out Time separation between outputs (seconds)
       char            *interpolation_lc, // out Linear or cubic interpolation
       char            *dp_file_rw,       // out (r)ead or (w) departure points for forward model
       double          *tracer_ob_bias,   // out bias added to observations of the tracer
       double          *flux_ob_bias,     // out bias added to observations of the flux
       double          *tc_ob_bias,       // out bias added to observation of total column tracer
       struct obs_type **obs )            // out observation network

{ // Declare local variables
  int              success, count, truefalse;
  bool             go;
  char             networktpe, ob_of;
  double           lon_start, lon_sep, lat_start, lat_sep, lev_start, lev_sep;
  double           t_start_s, t_sep_s, stddev;
  int              t_start_dy, t_start_hr, t_start_min;
  int              lon_n, lat_n, lev_n, t_n;
  int              lon, lat, lev, t;
  double           longitude, latitude, level, time_s, time_m, t_sep_min, earth_rot_fac;
  double           inc, period, omega, longitude0, sin_lat, cos_lat, sin_lon;
  struct obs_type  *current = NULL;
  double           v0[3], v1[3], v2[3];

  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;
  go      = true;

  if (argument_count > count)
  { strcpy (state_ic, argument_list[count]);
    printf ("Filename of the initial state and source: %s\n", state_ic);
  }
  else
  { success = 0;
    printf ("Please include the filename of the initial state and source\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (wind_dir, argument_list[count]);
    printf ("Directory containing the winds: %s\n", wind_dir);
  }
  else
  { success = 0;
    printf ("Please include directory containing the winds\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (obs_file, argument_list[count]);
    printf ("Output file for the synthetic observations: %s\n", obs_file);
  }
  else
  { success = 0;
    printf ("Please include output file for the synthetic observations\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (file_dps, argument_list[count]);
    printf ("File (semi-Lagrangian departure points): %s\n", file_dps);
  }
  else
  { success = 0;
    printf ("Please include filename for semi-Lagrangian departure points (or nil)\n");
  }

  count++;
  if (argument_count > count)
  { *Dt = atof(argument_list[count]);
    printf ("Time step between wind files (major time step, seconds): %f\n", *Dt);
  }
  else
  { success = 0;
    printf ("Please include the time step between wind files (major time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *dt = atof(argument_list[count]);
    printf ("Time step of integration scheme (minor time step, seconds): %f\n", *dt);
  }
  else
  { success = 0;
    printf ("Please include the time step of integration scheme (minor time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_dt = atof(argument_list[count]);
    printf ("Time step of diffusion scheme (seconds): %f\n", *kappa_dt);
  }
  else
  { success = 0;
    printf ("Please include the time step of diffusion scheme (seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_h = atof(argument_list[count]);
    printf ("Horizontal diffusion coefficient: %f\n", *kappa_h);
  }
  else
  { success = 0;
    printf ("Please include the horizontal diffusion coefficient\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_v = atof(argument_list[count]);
    printf ("Vertical diffusion coefficient: %f\n", *kappa_v);
  }
  else
  { success = 0;
    printf ("Please include the vertical diffusion coefficient\n");
  }

  count++;
  if (argument_count > count)
  { truefalse = atoi(argument_list[count]);
    *inc_vert = truefalse == 1;
    if (*inc_vert)
    { printf ("Include vertical transport\n");
    }
    else
    { printf ("Do not include vertical transport\n");
      printf ("***** NOTE: THIS OPTION SHOULD BE ACCOMPANIED BY factor_w = 0.0  *****\n");
    }
  }
  else
  { success = 0;
    printf ("Please include whether to include vertical transport (0/1)\n");
  }

  count++;
  if (argument_count > count)
  { *factor_w = atof(argument_list[count]);
    printf ("Multiplication factor of vertical winds: %f\n", *factor_w);
  }
  else
  { success = 0;
    printf ("Please include the multiplication factor of vertical winds\n");
  }

  count++;
  if (argument_count > count)
  { *output_freq = atof(argument_list[count]);
    printf ("Time separation between outputs (seconds): %f\n", *output_freq);
  }
  else
  { success = 0;
    printf ("Please include the time separation between outputs (seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *interpolation_lc = argument_list[count][0];
    printf ("Linear or cubic interpolation: %c\n", *interpolation_lc);
    if (*interpolation_lc != 'l' && *interpolation_lc != 'c')
    { printf ("Must enter l or c\n");
      success = 0;
    }
  }
  else
  { success = 0;
    printf ("Please include interpolation type (l=linear, c=cubic)\n");
  }

  count++;
  if (argument_count > count)
  { *dp_file_rw = argument_list[count][0];
    printf ("Read or write departure points file: %c\n", *dp_file_rw);
    if (*dp_file_rw != 'r' && *dp_file_rw != 'w')
    { printf ("Must enter r or w\n");
      success = 0;
    }
  }
  else
  { success = 0;
    printf ("Please include whether to (r)ead or (w)rite departure points file\n");
  }

  count++;
  if (argument_count > count)
  { *tracer_ob_bias = atof(argument_list[count]);
    printf ("Bias added to tracer observations: %f\n", *tracer_ob_bias);
  }
  else
  { success = 0;
    printf ("Please include the bias added to tracer observations\n");
  }

  count++;
  if (argument_count > count)
  { *flux_ob_bias = atof(argument_list[count]);
    printf ("Bias added to flux observations: %f\n", *flux_ob_bias);
  }
  else
  { success = 0;
    printf ("Please include the bias added to flux observations\n");
  }


  count++;
  if (argument_count > count)
  { *tc_ob_bias = atof(argument_list[count]);
    printf ("Bias added to total column tracer observations: %f\n", *tc_ob_bias);
  }
  else
  { success = 0;
    printf ("Please include the bias added to total column tracer observations\n");
  }


  // Initialise the list of observations to null
  *obs = NULL;

  while (go && success == 1)
  { // Find what the observation is of
    count++;
    if (argument_count > count)
    { ob_of = argument_list[count][0];
      printf ("Observation of: %c\n", ob_of);
      if (ob_of != 't' && ob_of != 'f' && ob_of != 'x' && ob_of !='e')
      { printf ("Must enter t (tracer), f (flux), x (total column tracer), or e (end)\n");
        success = 0;
      }
    }
    else
    { printf ("Must enter t (tracer), f (flux), x (total column tracer), or e (end)\n");
      success = 0;
    }

    if (ob_of != 'e')
    { // Find the ob network type
      count++;
      if (argument_count > count)
      { networktpe = argument_list[count][0];
        printf ("Ob network type: %c\n", networktpe);
        if (networktpe != 'i' && networktpe != 'g' && networktpe != 's')
        { printf ("Must enter i (individual observation), g (grid), or s (satellite)\n");
          success = 0;
        }
      }
      else
      { success = 0;
        printf ("Please include ob network type i (individual ob), g (grid), s (satellite), or e (end)\n");
      }


      switch (networktpe)
      { case 'i':  // ----- Individual observation ------------------------------
          // Create a new item
          if (current)
          { (*current).next = new struct obs_type;
            current         = (*current).next;
          }
          else
          { *obs    = new struct obs_type;
            current = (*obs);
          }

          count++;
          if (argument_count > count)
          { (*current).longitude = atof(argument_list[count]);
            if ((*current).longitude < 0.0)
            { (*current).longitude += 360.0;
            }
            printf ("Longitude : %f\n", (*current).longitude);
          }
          else
          { success = 0;
            printf ("Please include the observation's longitude\n");
          }

          count++;
          if (argument_count > count)
          { (*current).latitude = atof(argument_list[count]);
            printf ("Latitude : %f\n", (*current).latitude);
          }
          else
          { success = 0;
            printf ("Please include the observation's latitude\n");
          }

          if (ob_of == 't')
          { count++;
            if (argument_count > count)
            { (*current).level = atof(argument_list[count]);
              printf ("Height : %f\n", (*current).level);
            }
            else
            { success = 0;
              printf ("Please include the observation's height (m)\n");
            }
          }
          else
          { (*current).level = ::notdef;
          }

          count++;
          if (argument_count > count)
          { (*current).day = atoi(argument_list[count]);
            printf ("Day : %i\n", (*current).day);
          }
          else
          { success = 0;
            printf ("Please include the observation's time (day)\n");
          }

          count++;
          if (argument_count > count)
          { (*current).hour = atoi(argument_list[count]);
            printf ("Hour : %i\n", (*current).hour);
          }
          else
          { success = 0;
            printf ("Please include the observation's time (hour)\n");
          }

          count++;
          if (argument_count > count)
          { (*current).min = atoi(argument_list[count]);
            printf ("Minute : %i\n", (*current).min);
          }
          else
          { success = 0;
            printf ("Please include the observation's time (minute)\n");
          }

          // Work out the ob time in seconds
          (*current).obtime_secs = ( ( double((*current).day) * 24.0 +
                                       double((*current).hour)) * 60.0 +
                                       double((*current).min)) * 60.0;

          count++;
          if (argument_count > count)
          { (*current).stddev = atof(argument_list[count]);
            printf ("Err stddev : %f\n", (*current).stddev);
          }
          else
          { success = 0;
            printf ("Please include the observation's error stddev\n");
          }

          // Work out the ob's variance
          (*current).variance = (*current).stddev * (*current).stddev;
          (*current).ob_of    = ob_of;
          (*current).obtpe    = 'i';

          // Initialise the remaining variables
          (*current).lon_index     = -1;
          (*current).lat_index     = -1;
          (*current).lev_index     = -1;
          (*current).time_index    = -1;
          (*current).lon_alpha     = 0.0;
          (*current).lon_beta      = 0.0;
          (*current).lat_alpha     = 0.0;
          (*current).lat_beta      = 0.0;
          (*current).lev_alpha     = 0.0;
          (*current).lev_beta      = 0.0;
          (*current).time_alpha    = 0.0;
          (*current).time_beta     = 0.0;
          (*current).ob            = ::notdef;
          (*current).model_ob      = ::notdef;
          (*current).innov         = ::notdef;
          (*current).dJo_dmodel_ob = ::notdef;
          (*current).next = NULL;

        break;

        case 'g':  // ----- A grid of observations ------------------------------

          count++;
          if (argument_count > count)
          { lon_start = atof(argument_list[count]);
            printf ("Grid starting longitude (degrees): %f\n", lon_start);
          }
          else
          { success = 0;
            printf ("Please include the grid starting longitude (degrees)\n");
          }

          count++;
          if (argument_count > count)
          { lon_sep = atof(argument_list[count]);
            printf ("Grid longitude separation (degrees): %f\n", lon_sep);
          }
          else
          { success = 0;
            printf ("Please include the longitude separation (degrees)\n");
          }

          count++;
          if (argument_count > count)
          { lon_n = atoi(argument_list[count]);
            printf ("No. of obs in longitude direction: %i\n", lon_n);
          }
          else
          { success = 0;
            printf ("Please include the number of obs in longitude\n");
          }

          count++;
          if (argument_count > count)
          { lat_start = atof(argument_list[count]);
            printf ("Grid starting latitude (degrees): %f\n", lat_start);
          }
          else
          { success = 0;
            printf ("Please include the grid starting latitude (degrees)\n");
          }

          count++;
          if (argument_count > count)
          { lat_sep = atof(argument_list[count]);
            printf ("Grid latitude separation (degrees): %f\n", lat_sep);
          }
          else
          { success = 0;
            printf ("Please include the latitude separation (degrees)\n");
          }

          count++;
          if (argument_count > count)
          { lat_n = atoi(argument_list[count]);
            printf ("No. of obs in latitude direction: %i\n", lat_n);
          }
          else
          { success = 0;
            printf ("Please include the number of obs in latitude\n");
          }

          if (ob_of == 't')
          { count++;
            if (argument_count > count)
            { lev_start = atof(argument_list[count]);
              printf ("Grid starting level (m): %f\n", lev_start);
            }
            else
            { success = 0;
              printf ("Please include the grid starting level (m)\n");
            }

            count++;
            if (argument_count > count)
            { lev_sep = atof(argument_list[count]);
              printf ("Grid level separation (m): %f\n", lev_sep);
            }
            else
            { success = 0;
              printf ("Please include the level separation (m)\n");
            }

            count++;
            if (argument_count > count)
            { lev_n = atoi(argument_list[count]);
              printf ("No. of obs in vertical direction: %i\n", lev_n);
            }
            else
            { success = 0;
              printf ("Please include the number of obs in the vertical\n");
            }
          }
          else
          { lev_start = ::notdef;
            lev_sep   = 0.0;
            lev_n     = 1;
          }

          count++;
          if (argument_count > count)
          { t_start_dy = atoi(argument_list[count]);
            printf ("Grid starting time (days): %i\n", t_start_dy);
          }
          else
          { success = 0;
            printf ("Please include the grid starting time (days)\n");
          }

          count++;
          if (argument_count > count)
          { t_start_hr = atoi(argument_list[count]);
            printf ("Grid starting time (hours): %i\n", t_start_hr);
          }
          else
          { success = 0;
            printf ("Please include the grid starting time (hours)\n");
          }

          count++;
          if (argument_count > count)
          { t_start_min = atoi(argument_list[count]);
            printf ("Grid starting time (mins): %i\n", t_start_min);
          }
          else
          { success = 0;
            printf ("Please include the grid starting time (mins)\n");
          }

          t_start_s = (( double(t_start_dy) * 24.0 +
                         double(t_start_hr) ) * 60.0 +
                         double(t_start_min) ) * 60.0;

          count++;
          if (argument_count > count)
          { t_sep_min = atof(argument_list[count]);
            printf ("Grid time separation (mins): %f\n", t_sep_min);
          }
          else
          { success = 0;
            printf ("Please include the time separation (mins)\n");
          }

          t_sep_s = t_sep_min * 60.0;

          count++;
          if (argument_count > count)
          { t_n = atoi(argument_list[count]);
            printf ("No. of obs in time: %i\n", t_n);
          }
          else
          { success = 0;
            printf ("Please include the number of obs in time\n");
          }

          count++;
          if (argument_count > count)
          { stddev = atof(argument_list[count]);
            printf ("Standard deviation of the observation error: %f\n", stddev);
          }
          else
          { success = 0;
            printf ("Please include the standard deviation of the observation error\n");
          }

          if (success == 1)
          { // Loop round the observations to create the grid
            for (lon=0; lon<lon_n; lon++)
            { longitude  = (lon_start + double(lon) * lon_sep);
              longitude -= 360.0 * double(int(longitude/360.0));
              for (lat=0; lat<lat_n; lat++)
              { latitude = lat_start + double(lat) * lat_sep;
                if ((latitude >= -90.0) && (latitude <= 90.0))
                { for (lev=0; lev<lev_n; lev++)
                  { level = lev_start + double(lev) * lev_sep;
                    for (t=0; t<t_n; t++)
                    { time_s = t_start_s + double(t) * t_sep_s;

                      // Create a new item
                      if (current)
                      { (*current).next = new struct obs_type;
                        current         = (*current).next;
                      }
                      else
                      { *obs    = new struct obs_type;
                        current = (*obs);
                      }
                      (*current).longitude     = longitude;
                      (*current).latitude      = latitude;
                      (*current).level         = level;
                      Secs_to_ddhhmm (time_s,
                                      &((*current).day),
                                      &((*current).hour),
                                      &((*current).min));
                      (*current).obtime_secs   = time_s;
                      (*current).stddev        = stddev;
                      (*current).variance      = stddev * stddev;
                      (*current).ob_of         = ob_of;
                      (*current).obtpe         = 'g';

                      // Initialise the remaining variables
                      (*current).lon_index     = -1;
                      (*current).lat_index     = -1;
                      (*current).lev_index     = -1;
                      (*current).time_index    = -1;
                      (*current).lon_alpha     = 0.0;
                      (*current).lon_beta      = 0.0;
                      (*current).lat_alpha     = 0.0;
                      (*current).lat_beta      = 0.0;
                      (*current).lev_alpha     = 0.0;
                      (*current).lev_beta      = 0.0;
                      (*current).time_alpha    = 0.0;
                      (*current).time_beta     = 0.0;
                      (*current).ob            = ::notdef;
                      (*current).model_ob      = ::notdef;
                      (*current).innov         = ::notdef;
                      (*current).dJo_dmodel_ob = ::notdef;
                      (*current).next          = NULL;
                    }
                  }
                }
              }
            }
          }
        break;

        case 's':  // ----- Satellite observations ------------------------------

          count++;
          if (argument_count > count)
          { inc = atof(argument_list[count]);
            printf ("Orbit inclination (degrees): %f\n", inc);
          }
          else
          { success = 0;
            printf ("Please include the orbit inclination (degrees)\n");
          }

          count++;
          if (argument_count > count)
          { period = atof(argument_list[count]);
            printf ("Orbit period (minutes): %f\n", period);
            omega = 2.0 * ::pi / (60.0 * period);
          }
          else
          { success = 0;
            printf ("Please include the orbit period (minutes)\n");
          }

          count++;
          if (argument_count > count)
          { longitude0 = atof(argument_list[count]);
            printf ("Longitude of orbit at t=0 (degrees): %f\n", longitude0);
          }
          else
          { success = 0;
            printf ("Please include the longitude of orbit at t=0 (degrees)\n");
          }

          if (ob_of == 't')
          { count++;
            if (argument_count > count)
            { level = atof(argument_list[count]);
              printf ("Height : %f\n", level);
            }
            else
            { success = 0;
              printf ("Please include the observation's height (m)\n");
            }
          }
          else
          { level = ::notdef;
          }

          count++;
          if (argument_count > count)
          { t_sep_min = atof(argument_list[count]);
            printf ("Separation in time (minutes): %f\n", t_sep_min);
          }
          else
          { success = 0;
            printf ("Please include the separation in time (minutes)\n");
          }

          count++;
          if (argument_count > count)
          { t_n = atoi(argument_list[count]);
            printf ("Number in time: %i\n", t_n);
          }
          else
          { success = 0;
            printf ("Please include the number in time\n");
          }

          count++;
          if (argument_count > count)
          { stddev = atof(argument_list[count]);
            printf ("Standard deviation of the observation error: %f\n", stddev);
          }
          else
          { success = 0;
            printf ("Please include the standard deviation of the observation error\n");
          }

          // Convert inclination and longitude0 to radians
          inc        *= ::deg2rad;
          longitude0 *= ::deg2rad;


          if (success == 1)
          { // Factor allowing for Earth's rotation
            earth_rot_fac = 2.0 * ::pi / ::siderial_day;
            for (t=0; t<t_n; t++)
            { time_m = double(t) * t_sep_min;
              time_s = double(t) * t_sep_min * 60.0;
              // Set up the intial vector
              v0[0]  = sin(omega*time_s);
              v0[1]  = -1.0 * cos(omega*time_s);
              v0[2]  = 0.0;

              // Rotate in xz
              xz_rotate ( v0, v1, inc);
              // Rotate in xy (longitude0 is incremented due to Earth's rotation)
              xy_rotate ( v1, v2, longitude0 - time_m * earth_rot_fac);

              sin_lat    = v2[2];
              cos_lat    = sqrt(v2[0] * v2[0] + v2[1] * v2[1]);
              sin_lon    = v2[0] / cos_lat;
              latitude   = asin(sin_lat) * ::rad2deg;
              longitude  = asin(sin_lon) * ::rad2deg;
              // Adjustments
              if (v2[1] > 0.0) {longitude = copysign(180.0, longitude) - longitude;}
              if (longitude < 0.0) {longitude += 360.0;}

              // Loop in time to create the observation locations

              // Create a new item
              if (current)
              { (*current).next = new struct obs_type;
                current         = (*current).next;
              }
              else
              { *obs    = new struct obs_type;
                current = (*obs);
              }
              (*current).longitude     = longitude;
              (*current).latitude      = latitude;
              (*current).level         = level;
              Secs_to_ddhhmm (time_s,
                              &((*current).day),
                              &((*current).hour),
                              &((*current).min));
              (*current).obtime_secs   = time_s;
              (*current).stddev        = stddev;
              (*current).variance      = stddev * stddev;
              (*current).ob_of         = ob_of;
              (*current).obtpe         = 's';

              // Initialise the remaining variables
              (*current).lon_index     = -1;
              (*current).lat_index     = -1;
              (*current).lev_index     = -1;
              (*current).time_index    = -1;
              (*current).lon_alpha     = 0.0;
              (*current).lon_beta      = 0.0;
              (*current).lat_alpha     = 0.0;
              (*current).lat_beta      = 0.0;
              (*current).lev_alpha     = 0.0;
              (*current).lev_beta      = 0.0;
              (*current).time_alpha    = 0.0;
              (*current).time_beta     = 0.0;
              (*current).ob            = ::notdef;
              (*current).model_ob      = ::notdef;
              (*current).innov         = ::notdef;
              (*current).dJo_dmodel_ob = ::notdef;
              (*current).next          = NULL;
            }
          }
        break;
      }
    }
    else
    { go = false;
    }

  }
  printf ("===== End of program argument interpretation =====\n\n");
  return success;
}

// -------------------------------------------------------------------------------
void Secs_to_ddhhmm (double time_s,
                     int    *days,
                     int    *hrs,
                     int    *mins)
{ // Convert from total number of seconds to days, hours, and minutes
  // How many days?
  *days   = int(time_s / 86400.0);
  time_s -= 86400.0 * double(*days);
  *hrs    = int(time_s / 3600.0);
  time_s -= 3600.0 * double(*hrs);
  *mins   = int(time_s / 60.0);
}

// -------------------------------------------------------------------------------
// Routines to do rotation to generate positions for satellite obs
void xy_rotate( double v_in[3],
                double v_out[3],
                double angle_rad )
{ // Rotate in the xy plane
  double s, c;

  s = sin(angle_rad);
  c = cos(angle_rad);
  v_out[0] = c * v_in[0] - s * v_in[1];
  v_out[1] = s * v_in[0] + c * v_in[1];
  v_out[2] = v_in[2];
}

void xz_rotate( double v_in[3],
                double v_out[3],
                double angle_rad )
{ // Rotate in the xy plane
  double s, c;

  s = sin(angle_rad);
  c = cos(angle_rad);
  v_out[0] = c * v_in[0] - s * v_in[2];
  v_out[1] = v_in[1];
  v_out[2] = s * v_in[0] + c * v_in[2];
}


// -------------------------------------------------------------------------------
int MakeBg_arguments (
       int    argument_count,     // in
       char   **argument_list,    // in
       char   truth_in[256],      // out Filename of truth state (initial condition and source)
       char   bg_out[256],        // out Filename of background state (initial condition and source)
       char   bgpert_out[256],    // out Filename of background perturbation (initial condition and source)
       char   CVTfilename[256],   // out CVT filename
       double *factor_std_tracer, // out To multiply tracer err std
       double *factor_std_flux,   // out To multiply flux err std
       double *factor_pert_tracer,// out To multiply tracer bg pert
       double *factor_pert_flux)  // out To multiply flux bg pert

{ // Declare local variables
  int              success, count;

  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;

  if (argument_count > count)
  { strcpy (truth_in, argument_list[count]);
    printf ("Input filename of the truth state: %s\n", truth_in);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the truth state\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (bg_out, argument_list[count]);
    printf ("Output filename of the background state: %s\n", bg_out);
  }
  else
  { success = 0;
    printf ("Please include the output filename of the background state\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (bgpert_out, argument_list[count]);
    printf ("Output filename of the background perturbation: %s\n", bgpert_out);
  }
  else
  { success = 0;
    printf ("Please include the output filename of the background perturbation\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (CVTfilename, argument_list[count]);
    printf ("Input filename of the CVT file: %s\n", CVTfilename);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the CVT file\n");
  }

  count++;
  if (argument_count > count)
  { *factor_std_tracer = atof(argument_list[count]);
    printf ("Multipication factor for tracer err std: %f\n", *factor_std_tracer);
  }
  else
  { success = 0;
    printf ("Please include the multipication factor for tracer err std\n");
  }

  count++;
  if (argument_count > count)
  { *factor_std_flux = atof(argument_list[count]);
    printf ("Multipication factor for flux err std: %f\n", *factor_std_flux);
  }
  else
  { success = 0;
    printf ("Please include the multipication factor for flux err std\n");
  }

  count++;
  if (argument_count > count)
  { *factor_pert_tracer = atof(argument_list[count]);
    printf ("Multipication factor for tracer bg pert: %f\n", *factor_pert_tracer);
  }
  else
  { success = 0;
    printf ("Please include the multipication factor for tracer bg pert\n");
  }

  count++;
  if (argument_count > count)
  { *factor_pert_flux = atof(argument_list[count]);
    printf ("Multipication factor for flux bg pert: %f\n", *factor_pert_flux);
  }
  else
  { success = 0;
    printf ("Please include the multipication factor for flux bg pert\n");
  }

  printf ("===== End of program argument interpretation =====\n\n");
  return success;
}





// -------------------------------------------------------------------------------
int Assimilate_arguments (
       int    argument_count,           // in
       char   **argument_list,          // in
       char   bg_filename[256],         // out Filename of background state
       char   wind_dir[256],            // out Directory containing driver winds
       char   obsfile_in[256],          // out Input file for observations
       char   file_dps[256],            // out Input or output file for departure points
       char   CVTfilename[256],         // out CVT filename (specification of B)
       double *factor_std_tracer,       // out To multiply tracer err std
       double *factor_std_flux,         // out To multiply flux err std
       char   obsfile_out_bg[256],      // out Output file of observation information at the background
       char   anal_filename[256],       // out Filename of analysis state
       char   analinc_filename[256],    // out Filename of analysis increment
       char   Assim_diags_filename[256],// out Filename of assimilation diagnostics
       char   obsfile_out_anal[256],    // out Output file of observation information at the analysis
       double *Dt,                      // out Major timestep (between wind files, seconds)
       double *dt,                      // out Minor timestep (for integration scheme, seconds)
       double *kappa_dt,                // out Timestep for diffusion
       double *kappa_h,                 // out Horizontal diffusion coefficient
       double *kappa_v,                 // out Vertical diffusion coefficient
       bool   *inc_vert,                // out Include vertical transport?
       double *factor_w,                // out Multiplication factor of vertical winds
       char   *interpolation_lc,        // out Linear or cubic interpolation
       char   *descent_algorithm,       // out Descent algorithm type
       double *conv_criterion,          // out Criterion for gradient norm for convergence
       int    *max_iters )              // out Maximum number of iterations

{ // Declare local variables
  int              success, count, truefalse;

  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;

  if (argument_count > count)
  { strcpy (bg_filename, argument_list[count]);
    printf ("Background state: %s\n", bg_filename);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the background state\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (wind_dir, argument_list[count]);
    printf ("Wind directory: %s\n", wind_dir);
  }
  else
  { success = 0;
    printf ("Please include the wind directory\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (obsfile_in, argument_list[count]);
    printf ("Input filename of the observations: %s\n", obsfile_in);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the observations\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (file_dps, argument_list[count]);
    printf ("Filename of the departure points: %s\n", file_dps);
  }
  else
  { success = 0;
    printf ("Please include the filename of the departure points\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (CVTfilename, argument_list[count]);
    printf ("Input filename of the CVT file: %s\n", CVTfilename);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the CVT file\n");
  }

  count++;
  if (argument_count > count)
  { *factor_std_tracer = atof(argument_list[count]);
    printf ("Multipication factor for tracer err std: %f\n", *factor_std_tracer);
  }
  else
  { success = 0;
    printf ("Please include the multipication factor for tracer err std\n");
  }

  count++;
  if (argument_count > count)
  { *factor_std_flux = atof(argument_list[count]);
    printf ("Multipication factor for flux err std: %f\n", *factor_std_flux);
  }
  else
  { success = 0;
    printf ("Please include the multipication factor for flux err std\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (obsfile_out_bg, argument_list[count]);
    printf ("Output filename of the observations at the background: %s\n", obsfile_out_bg);
  }
  else
  { success = 0;
    printf ("Please include the output filename of the observations at the background\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (anal_filename, argument_list[count]);
    printf ("Output filename of the analysis: %s\n", anal_filename);
  }
  else
  { success = 0;
    printf ("Please include the output filename of the analysis\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (analinc_filename, argument_list[count]);
    printf ("Output filename of the analysis perturbation: %s\n", analinc_filename);
  }
  else
  { success = 0;
    printf ("Please include the output filename of analysis perturbation\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (Assim_diags_filename, argument_list[count]);
    printf ("Output filename of the assimilation diagnostics: %s\n", Assim_diags_filename);
  }
  else
  { success = 0;
    printf ("Please include the output filename of assimilation diagnostics\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (obsfile_out_anal, argument_list[count]);
    printf ("Output filename of the observations at the analysis: %s\n", obsfile_out_anal);
  }
  else
  { success = 0;
    printf ("Please include the output filename of the observations at the analysis\n");
  }

  count++;
  if (argument_count > count)
  { *Dt = atof(argument_list[count]);
    printf ("Time step between wind files (major time step, seconds): %f\n", *Dt);
  }
  else
  { success = 0;
    printf ("Please include the time step between wind files (major time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *dt = atof(argument_list[count]);
    printf ("Time step of integration scheme (minor time step, seconds): %f\n", *dt);
  }
  else
  { success = 0;
    printf ("Please include the time step of integration scheme (minor time step, seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_dt = atof(argument_list[count]);
    printf ("Time step of diffusion scheme (seconds): %f\n", *kappa_dt);
  }
  else
  { success = 0;
    printf ("Please include the time step of diffusion scheme (seconds)\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_h = atof(argument_list[count]);
    printf ("Horizontal diffusion coefficient: %f\n", *kappa_h);
  }
  else
  { success = 0;
    printf ("Please include the horizontal diffusion coefficient\n");
  }

  count++;
  if (argument_count > count)
  { *kappa_v = atof(argument_list[count]);
    printf ("Vertical diffusion coefficient: %f\n", *kappa_v);
  }
  else
  { success = 0;
    printf ("Please include the vertical diffusion coefficient\n");
  }

  count++;
  if (argument_count > count)
  { truefalse = atoi(argument_list[count]);
    *inc_vert = truefalse == 1;
    if (*inc_vert)
    { printf ("Include vertical transport\n");
    }
    else
    { printf ("Do not include vertical transport\n");
      printf ("***** NOTE: THIS OPTION SHOULD BE ACCOMPANIED BY factor_w = 0.0  *****\n");
    }
  }
  else
  { success = 0;
    printf ("Please include whether to include vertical transport (0/1)\n");
  }

  count++;
  if (argument_count > count)
  { *factor_w = atof(argument_list[count]);
    printf ("Multiplication factor of vertical winds: %f\n", *factor_w);
  }
  else
  { success = 0;
    printf ("Please include the multiplication factor of vertical winds\n");
  }

  count++;
  if (argument_count > count)
  { *interpolation_lc = argument_list[count][0];
    printf ("Linear or cubic interpolation: %c\n", *interpolation_lc);
    if (*interpolation_lc != 'l' && *interpolation_lc != 'c')
    { printf ("Must enter l or c\n");
      success = 0;
    }
  }
  else
  { success = 0;
    printf ("Please include interpolation type (l=linear, c=cubic)\n");
  }

  count++;
  if (argument_count > count)
  { *descent_algorithm = argument_list[count][0];
    printf ("Conjugate gradient or quasi-Newton descent algorithm: %c\n", *descent_algorithm);
    if (*descent_algorithm != 'c' && *descent_algorithm != 'q')
    { printf ("Must enter c or q\n");
      success = 0;
    }
  }
  else
  { success = 0;
    printf ("Please include descent algorithm type (c=conjugate gradient, q=quasi-Newton)\n");
  }

  count++;
  if (argument_count > count)
  { *conv_criterion = atof(argument_list[count]);
    printf ("Convergence criterion: %f\n", *conv_criterion);
  }
  else
  { success = 0;
    printf ("Please include the convergence criterion\n");
  }

  count++;
  if (argument_count > count)
  { *max_iters = atoi(argument_list[count]);
    printf ("Maximum number of variational iterations: %i\n", *max_iters);
  }
  else
  { success = 0;
    printf ("Please include the maximum number of variational iterations\n");
  }

  printf ("===== End of program argument interpretation =====\n\n");
  return success;
}


// -------------------------------------------------------------------------------
int CalcError_arguments (
       int    argument_count,           // in
       char   **argument_list,          // in
       char   State_file[256],          // out Filename of state
       char   Truth_file[256],          // out Filename of truth
       char   Background_file[256],     // out Filename of background
       char   CVTfilename[256],         // out CVT filename
       char   Output_file[256] )        // out Filename of data file containing error diagnostics

{ // Declare local variables
  int              success, count;

  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;

  if (argument_count > count)
  { strcpy (State_file, argument_list[count]);
    printf ("State filename: %s\n", State_file);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the state\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (Truth_file, argument_list[count]);
    printf ("Truth filename: %s\n", Truth_file);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the truth\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (Background_file, argument_list[count]);
    printf ("Background filename: %s\n", Background_file);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the background\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (CVTfilename, argument_list[count]);
    printf ("CVT filename: %s\n", CVTfilename);
  }
  else
  { success = 0;
    printf ("Please include the input CVT filename\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (Output_file, argument_list[count]);
    printf ("Error diagnostics output file: %s\n", Output_file);
  }
  else
  { success = 0;
    printf ("Please include the output filename for the diagnostics\n");
  }

  printf ("===== End of program argument interpretation =====\n\n");

  return success;
}


int Invicat2Enviflux_arguments (
       int    argument_count,           // in
       char   **argument_list,          // in
       char   CVTfilename[256],         // out CVT filename
       char   INVICAT_file[256],        // out INVICAT filename
       char   ENVIFLUX_file[256],       // out ENVIFLUX filename for output
       double *maxflux,                 // out Fluxes will be scaled to have this max value
       double *maxtracer )              // out Tracer will be scaled to have this max value

{ // Declare local variables
  int              success, count, truefalse;

  // Return 1 from subroutine if OK, 0 otherwise
  success = 1;
  count   = 1;

  if (argument_count > count)
  { strcpy (CVTfilename, argument_list[count]);
    printf ("Input filename of the CVT file: %s\n", CVTfilename);
  }
  else
  { success = 0;
    printf ("Please include the input filename of the CVT file\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (INVICAT_file, argument_list[count]);
    printf ("INVICAT input file: %s\n", INVICAT_file);
  }
  else
  { success = 0;
    printf ("Please include the INVICAT input filename\n");
  }

  count++;
  if (argument_count > count)
  { strcpy (ENVIFLUX_file, argument_list[count]);
    printf ("ENVIFLUX output filename: %s\n", ENVIFLUX_file);
  }
  else
  { success = 0;
    printf ("Please include the ENVIFLUX output filename\n");
  }

  count++;
  if (argument_count > count)
  { *maxflux = atof(argument_list[count]);
    printf ("Maximum flux set to : %f\n", *maxflux);
    printf ("Value of 999.9 means no scaling is done for the fluxes\n");
  }
  else
  { success = 0;
    printf ("Please include the maximum flux value (999.0 for no scaling)\n");
  }

  count++;
  if (argument_count > count)
  { *maxtracer = atof(argument_list[count]);
    printf ("Maximum tracer set to : %f\n", *maxtracer);
    printf ("Value of 999.9 means no scaling is done for the tracer\n");
  }
  else
  { success = 0;
    printf ("Please include the maximum tracer value (999.0 for no scaling)\n");
  }

  printf ("===== End of program argument interpretation =====\n\n");
  return success;
}
