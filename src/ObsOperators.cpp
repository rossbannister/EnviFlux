   #include <source.h>


// -------------------------------------------------------------------------------
void ObservationOperator
 ( struct obs_type            *ob,            // inout Observation and data
   struct obs_trunk_type      *ob_trunk,      // in    Containing the mass profile (for satellites)
   struct metadata_type       *MetaData,      // in    Lats, lons, etc.
   struct instant_tracer_type *field_lowert,  // in    Tracer field at lower time
   struct instant_tracer_type *field_uppert,  // in    Tracer field at upper time
   double                     lowert_tracer,  // in    Lower time of tracer field in s
   double                     uppert_tracer,  // in    Upper time of tracer field in s
   struct state_type          *state )        // in    Contains flux field


{ // Function to produce model observation (may be of tracer or flux)
  double delta_recip, time1, time2, full_lev_plus1, full_lev_here, density;
  int    index, lev;

  // ===== Find the array indices and interpolation weights =====
  // Find the longitude index if necessary
  if ((*ob).lon_index == -1)
  { index           = Find_index_ascend ( (*MetaData).nlon+2,     // Number of elements
                                          (*MetaData).longitude,  // Ascending array
                                          (*ob).longitude );      // Longitude of ob
    (*ob).lon_index = index;
    //printf ("Observation longitude %f is at index %i between %f and %f\n",
    //        (*ob).longitude, index,
    //        (*MetaData).longitude[index], (*MetaData).longitude[index+1]);
    delta_recip     = 1.0 / ((*MetaData).longitude[index+1] - (*MetaData).longitude[index]);
    (*ob).lon_alpha = delta_recip * ((*MetaData).longitude[index+1] - (*ob).longitude);
    (*ob).lon_beta  = delta_recip * ((*ob).longitude - (*MetaData).longitude[index]);
    //printf ("  Longitude %i %fl %fl\n", index, (*ob).lon_alpha, (*ob).lon_beta);
  }


  // Find the latitude index if necessary
  if ((*ob).lat_index == -1)
  { index           = Find_index_descend_mod ( (*MetaData).nlat+2,    // Number of elements
                                               (*MetaData).latitude,  // Descending array
                                               (*ob).latitude );      // Latitude of ob
    (*ob).lat_index = index;
    //printf ("Observation latitude %f is at index %i between %f and %f\n",
    //        (*ob).latitude, index,
    //        (*MetaData).latitude[index], (*MetaData).latitude[index+1]);
    delta_recip     = 1.0 / ((*MetaData).latitude[index+1] - (*MetaData).latitude[index]);
    (*ob).lat_alpha = delta_recip * ((*MetaData).latitude[index+1] - (*ob).latitude);
    (*ob).lat_beta  = delta_recip * ((*ob).latitude - (*MetaData).latitude[index]);
    //printf ("  Latitude %i %fl %fl\n", index, (*ob).lat_alpha, (*ob).lat_beta);
  }


  // Find the level index if necessary
  if (((*ob).ob_of == 't') && ((*ob).lev_index == -1))
  { index           = Find_index_ascend ( (*MetaData).nlev+2,     // Number of elements
                                          (*MetaData).level,      // Ascending array
                                          (*ob).level );          // Longitude of ob
    (*ob).lev_index = index;
    //printf ("Observation level %f is at index %i between %f and %f\n",
    //        (*ob).level, index,
    //        (*MetaData).level[index], (*MetaData).level[index+1]);
    delta_recip     = 1.0 / ((*MetaData).level[index+1] - (*MetaData).level[index]);
    (*ob).lev_alpha = delta_recip * ((*MetaData).level[index+1] - (*ob).level);
    (*ob).lev_beta  = delta_recip * ((*ob).level - (*MetaData).level[index]);
    //printf ("  Lev %i %fl %fl\n", index, (*ob).lev_alpha, (*ob).lev_beta);
  }


  // Find the time index if necessary
  if ((*ob).time_index == -1)
  { if ((*ob).ob_of == 't' || (*ob).ob_of == 'x')
    { // Observation is of the tracer or total column, time index refers to model timestep
      (*ob).time_index = 0;    // Not needed for tracer or total column observation
      delta_recip      = 1.0 / (uppert_tracer - lowert_tracer);
      (*ob).time_alpha = delta_recip * (uppert_tracer - (*ob).obtime_secs);
      (*ob).time_beta  = delta_recip * ((*ob).obtime_secs - lowert_tracer);
      //printf ("  t (tracer) %i %fl %fl\n", index, (*ob).time_alpha, (*ob).time_beta);
    }
    else
    { // Observation is of the flux, time index refers to flux timestep
      index            = Find_index_ascend ( (*MetaData).nss,     // Number of elements
                                             (*MetaData).times,   // Ascending array
                                             (*ob).obtime_secs);  // Time of ob
      (*ob).time_index = index;
      if (index == (*MetaData).ntimes_major-1)
      { // The last time
        (*ob).time_alpha = 1.0;
        (*ob).time_beta  = 0.0;
      }
      else
      { delta_recip     = 1.0 / ((*MetaData).times[index+1] - (*MetaData).times[index]);
        (*ob).time_alpha = delta_recip * ((*MetaData).times[index+1] - (*ob).obtime_secs);
        (*ob).time_beta  = delta_recip * ((*ob).obtime_secs - (*MetaData).times[index]);
      }
      //printf ("  t (flux) %i %fl %fl\n", index, (*ob).time_alpha, (*ob).time_beta);
    }
  }



  // ===== Do the interpolation =====
  if ((*ob).ob_of == 't')
  { // Observation is of the tracer (point observation)
    // Do the spatial interpolation at the lower time
    time1 = ob_op_interp_tracer_3D (ob,
                                    field_lowert);
    if ((*ob).time_alpha < 1.0)
    { // Do the spatial interpolation at the upper time
      time2 = ob_op_interp_tracer_3D (ob,
                                      field_uppert);
    }
    else
    { time2 = 0.0;
    }
    (*ob).model_ob = (*ob).time_alpha * time1 + (*ob).time_beta * time2;
  }
  else
  { if ((*ob).ob_of == 'x')
    { // Observation is a 'total column' with an averaging kernel for total amount
      // Perform operation for the lower time
      time1 = 0.0;
      for (lev=1; lev<=(*MetaData).nlev; lev++)
      { time1 += ob_op_interp_tracer_2D (ob,
                                         field_lowert,
                                         lev) *
                 (*ob_trunk).mass_profile[lev];
      }
      time2 = 0.0;
      if ((*ob).time_alpha < 1.0)
      { // Perform operation for the upper time
        for (lev=1; lev<=(*MetaData).nlev; lev++)
        { time2 += ob_op_interp_tracer_2D (ob,
                                           field_uppert,
                                           lev) *
                   (*ob_trunk).mass_profile[lev];
        }
      }
      (*ob).model_ob = (*ob).time_alpha * time1 + (*ob).time_beta * time2;
    }
    else
    { // Observation is of flux
      time1 = ob_op_interp_flux_2D (ob,
                                    (*state).source_rs,
                                    (*ob).time_index);
      if ((*ob).time_alpha < 1.0)
      { // Do the spatial interpolation at the upper time
        time2 = ob_op_interp_flux_2D (ob,
                                      (*state).source_rs,
                                      (*ob).time_index+1);
      }
      else
      { time2 = 0.0;
      }
      (*ob).model_ob = (*ob).time_alpha * time1 + (*ob).time_beta * time2;
    }
  }
}



// -------------------------------------------------------------------------------
void ObservationOperator_adj
 ( struct obs_type            *ob_hat,            // in    Observation and data
   struct obs_trunk_type      *ob_trunk,          // in    Containing the mass profile (for satellites)
   struct metadata_type       *MetaData,          // in    Lats, lons, etc.
   struct instant_tracer_type *field_lowert_hat,  // inout Tracer field at lower time
   struct instant_tracer_type *field_uppert_hat,  // inout Tracer field at upper time
   struct state_type          *state_hat )        // inout Contains flux field


{ // Adjoint of: Function to produce model observation (may be of tracer or flux)
  double time1_hat, time2_hat, inter_hat;
  int    lev;

  // ===== Do the interpolation =====
  if ((*ob_hat).ob_of == 't')
  { // Observation is of the tracer (point observation)
    // (*ob).model_ob = (*ob).time_alpha * time1 + (*ob).time_beta * time2;
    time1_hat = (*ob_hat).dJo_dmodel_ob * (*ob_hat).time_alpha;
    time2_hat = (*ob_hat).dJo_dmodel_ob * (*ob_hat).time_beta;

    // Do the spatial interpolation at the lower time
    // time1 = ob_op_interp_tracer_3D (ob,
    //                                 field_lowert);
    ob_op_interp_tracer_3D_adj (ob_hat,
                                field_lowert_hat,
                                time1_hat);

    if ((*ob_hat).time_alpha < 1.0)
    { // Do the spatial interpolation at the upper time
      // time2 = ob_op_interp_tracer_3D (ob,
      //                                 field_uppert);
      ob_op_interp_tracer_3D_adj (ob_hat,
                                  field_uppert_hat,
                                  time2_hat);
    }
  }
  else
  { if ((*ob_hat).ob_of == 'x')
    { // Observation is a 'total column' with an averaging kernel for total amount

      // (*ob).model_ob = (*ob).time_alpha * time1 + (*ob).time_beta * time2;
      time1_hat = (*ob_hat).dJo_dmodel_ob * (*ob_hat).time_alpha;
      time2_hat = (*ob_hat).dJo_dmodel_ob * (*ob_hat).time_beta;

      // Perform operation for the lower time
      // time1 = 0.0;
      for (lev=1; lev<=(*MetaData).nlev; lev++)
      { // inter  = ob_op_interp_tracer_2D (ob,
        //                                  field_lowert,
        //                                  lev)
        // time1 += inter * (*ob_trunk).mass_profile[lev];
        inter_hat = time1_hat * (*ob_trunk).mass_profile[lev];
        ob_op_interp_tracer_2D_adj (ob_hat,
                                    field_lowert_hat,
                                    lev,
                                    inter_hat);
      }

      // time2 = 0.0;
      if ((*ob_hat).time_alpha < 1.0)
      { // Perform operation for the upper time
        for (lev=1; lev<=(*MetaData).nlev; lev++)
        { // inter  = ob_op_interp_tracer_2D (ob,
          //                                  field_uppert,
          //                                  lev)
          // time2 += inter * (*ob_trunk).mass_profile[lev];
          inter_hat = time2_hat * (*ob_trunk).mass_profile[lev];
          ob_op_interp_tracer_2D_adj (ob_hat,
                                      field_uppert_hat,
                                      lev,
                                      inter_hat);
        }
      }
    }
    else
    { // Observation is of flux
      // (*ob).model_ob = (*ob).time_alpha * time1 + (*ob).time_beta * time2;
      time1_hat = (*ob_hat).dJo_dmodel_ob * (*ob_hat).time_alpha;
      time2_hat = (*ob_hat).dJo_dmodel_ob * (*ob_hat).time_beta;

      // time1 = ob_op_interp_flux_2D (ob,
      //                               (*state).source_rs,
      //                               (*ob).time_index);
      ob_op_interp_flux_2D_adj (ob_hat,
                                (*state_hat).source_rs,
                              (*ob_hat).time_index,
                              time1_hat);

      if ((*ob_hat).time_alpha < 1.0)
      { // Do the spatial interpolation at the upper time
        // time2 = ob_op_interp_flux_2D (ob,
        //                               (*state).source_rs,
        //                               (*ob).time_index+1);
        ob_op_interp_flux_2D_adj (ob_hat,
                                  (*state_hat).source_rs,
                                  (*ob_hat).time_index+1,
                                  time2_hat);
      }
    }
  }
}



// -------------------------------------------------------------------------------
double ob_op_interp_tracer_3D ( struct obs_type            *ob,
                                struct instant_tracer_type *field ) // in
{ // Interpolate in 3D as part of the observation operator (this is for the tracer)
  double lon11, lon12, lon21, lon22, lat1, lat2, result;

  // Do interpolation in longitude
  lon11  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index]  [(*ob).lev_index] +
           (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index]  [(*ob).lev_index];
  lon12  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index+1][(*ob).lev_index] +
           (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index+1][(*ob).lev_index];
  lon21  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index]  [(*ob).lev_index+1] +
           (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index]  [(*ob).lev_index+1];
  lon22  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index+1][(*ob).lev_index+1] +
           (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index+1][(*ob).lev_index+1];

  // Do interpolation in latitude
  lat1   = (*ob).lat_alpha * lon11 + (*ob).lat_beta * lon12;
  lat2   = (*ob).lat_alpha * lon21 + (*ob).lat_beta * lon22;

  // Do interpolation in height
  result = (*ob).lev_alpha * lat1 + (*ob).lev_beta * lat2;

  return result;
}


// -------------------------------------------------------------------------------
void ob_op_interp_tracer_3D_adj ( struct obs_type            *ob,
                                  struct instant_tracer_type *field_hat,  // inout
                                  double                     result_hat ) // in
{ // Adjoint of: Interpolate in 3D as part of the observation operator (this is for the tracer)
  double lon11_hat, lon12_hat, lon21_hat, lon22_hat, lat1_hat, lat2_hat;

  // Do interpolation in height
  // result = (*ob).lev_alpha * lat1 + (*ob).lev_beta * lat2;
  lat1_hat = result_hat * (*ob).lev_alpha;
  lat2_hat = result_hat * (*ob).lev_beta;

  // Do interpolation in latitude
  // lat1   = (*ob).lat_alpha * lon11 + (*ob).lat_beta * lon12;
  // lat2   = (*ob).lat_alpha * lon21 + (*ob).lat_beta * lon22;
  lon11_hat = lat1_hat * (*ob).lat_alpha;
  lon12_hat = lat1_hat * (*ob).lat_beta;
  lon21_hat = lat2_hat * (*ob).lat_alpha;
  lon22_hat = lat2_hat * (*ob).lat_beta;

  // Do interpolation in longitude
  // lon11  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index]  [(*ob).lev_index] +
  //          (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index]  [(*ob).lev_index];
  // lon12  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index+1][(*ob).lev_index] +
  //          (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index+1][(*ob).lev_index];
  // lon21  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index]  [(*ob).lev_index+1] +
  //          (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index]  [(*ob).lev_index+1];
  // lon22  = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index+1][(*ob).lev_index+1] +
  //          (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index+1][(*ob).lev_index+1];

  (*field_hat).tracer[(*ob).lon_index]  [(*ob).lat_index]  [(*ob).lev_index]   += lon11_hat * (*ob).lon_alpha;
  (*field_hat).tracer[(*ob).lon_index+1][(*ob).lat_index]  [(*ob).lev_index]   += lon11_hat * (*ob).lon_beta;
  (*field_hat).tracer[(*ob).lon_index]  [(*ob).lat_index+1][(*ob).lev_index]   += lon12_hat * (*ob).lon_alpha;
  (*field_hat).tracer[(*ob).lon_index+1][(*ob).lat_index+1][(*ob).lev_index]   += lon12_hat * (*ob).lon_beta;
  (*field_hat).tracer[(*ob).lon_index]  [(*ob).lat_index]  [(*ob).lev_index+1] += lon21_hat * (*ob).lon_alpha;
  (*field_hat).tracer[(*ob).lon_index+1][(*ob).lat_index]  [(*ob).lev_index+1] += lon21_hat * (*ob).lon_beta;
  (*field_hat).tracer[(*ob).lon_index]  [(*ob).lat_index+1][(*ob).lev_index+1] += lon22_hat * (*ob).lon_alpha;
  (*field_hat).tracer[(*ob).lon_index+1][(*ob).lat_index+1][(*ob).lev_index+1] += lon22_hat * (*ob).lon_beta;
}

// -------------------------------------------------------------------------------
double ob_op_interp_flux_2D ( struct obs_type *ob,
                              double          ***field,
                              int             tindex )
{ // Interpolate in 2D as part of the observation operator (this is for the flux)
  double lon1, lon2, result;

  // Do interpolation in longitude
  lon1 = (*ob).lon_alpha * field[(*ob).lon_index]  [(*ob).lat_index]  [tindex] +
         (*ob).lon_beta  * field[(*ob).lon_index+1][(*ob).lat_index]  [tindex];

  lon2 = (*ob).lon_alpha * field[(*ob).lon_index]  [(*ob).lat_index+1][tindex] +
         (*ob).lon_beta  * field[(*ob).lon_index+1][(*ob).lat_index+1][tindex];

  // Do interpolation in latitude
  result = (*ob).lat_alpha * lon1 + (*ob).lat_beta * lon2;

  return result;
}


// -------------------------------------------------------------------------------
void ob_op_interp_flux_2D_adj ( struct obs_type *ob,
                                double          ***field_hat, // inout
                                int             tindex,
                                double          result_hat )  // in
{ // Adjoint of: Interpolate in 2D as part of the observation operator (this is for the flux)
  double lon1_hat, lon2_hat;

  // Do interpolation in latitude
  // result = (*ob).lat_alpha * lon1 + (*ob).lat_beta * lon2;
  lon1_hat = result_hat * (*ob).lat_alpha;
  lon2_hat = result_hat * (*ob).lat_beta;

  // Do interpolation in longitude
  // lon1 = (*ob).lon_alpha * field[(*ob).lon_index]  [(*ob).lat_index]  [tindex] +
  //        (*ob).lon_beta  * field[(*ob).lon_index+1][(*ob).lat_index]  [tindex];
  field_hat[(*ob).lon_index]  [(*ob).lat_index]  [tindex] += lon1_hat * (*ob).lon_alpha;
  field_hat[(*ob).lon_index+1][(*ob).lat_index]  [tindex] += lon1_hat * (*ob).lon_beta;

  // lon2 = (*ob).lon_alpha * field[(*ob).lon_index]  [(*ob).lat_index+1][tindex] +
  //        (*ob).lon_beta  * field[(*ob).lon_index+1][(*ob).lat_index+1][tindex];
  field_hat[(*ob).lon_index]  [(*ob).lat_index+1][tindex] += lon2_hat * (*ob).lon_alpha;
  field_hat[(*ob).lon_index+1][(*ob).lat_index+1][tindex] += lon2_hat * (*ob).lon_beta;
}



// -------------------------------------------------------------------------------
double ob_op_interp_tracer_2D ( struct obs_type            *ob,
                                struct instant_tracer_type *field,
                                int                        lev )
{ // Interpolate in 2D as part of the observation operator (this is for the tracer)
  double lon1, lon2, result;

  // Do interpolation in longitude
  lon1 = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index]  [lev] +
         (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index]  [lev];

  lon2 = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index+1][lev] +
         (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index+1][lev];

  // Do interpolation in latitude
  result = (*ob).lat_alpha * lon1 + (*ob).lat_beta * lon2;

  return result;
}

// -------------------------------------------------------------------------------
void ob_op_interp_tracer_2D_adj ( struct obs_type            *ob,
                                  struct instant_tracer_type *field_hat,  // inout
                                  int                        lev,
                                  double                     result_hat ) // in
{ // Adjoint of: Interpolate in 2D as part of the observation operator (this is for the tracer)
  double lon1_hat, lon2_hat;

  // Do interpolation in latitude
  // result = (*ob).lat_alpha * lon1 + (*ob).lat_beta * lon2;
  lon1_hat = result_hat * (*ob).lat_alpha;
  lon2_hat = result_hat * (*ob).lat_beta;

  // Do interpolation in longitude
  // lon1 = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index]  [lev] +
  //        (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index]  [lev];
  (*field_hat).tracer[(*ob).lon_index]  [(*ob).lat_index]  [lev] += lon1_hat * (*ob).lon_alpha;
  (*field_hat).tracer[(*ob).lon_index+1][(*ob).lat_index]  [lev] += lon1_hat * (*ob).lon_beta;

  // lon2 = (*ob).lon_alpha * (*field).tracer[(*ob).lon_index]  [(*ob).lat_index+1][lev] +
  //        (*ob).lon_beta  * (*field).tracer[(*ob).lon_index+1][(*ob).lat_index+1][lev];
  (*field_hat).tracer[(*ob).lon_index]  [(*ob).lat_index+1][lev] += lon2_hat * (*ob).lon_alpha;
  (*field_hat).tracer[(*ob).lon_index+1][(*ob).lat_index+1][lev] += lon2_hat * (*ob).lon_beta;
}
