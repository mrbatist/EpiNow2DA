// calculate Rt directly from inferred infections
vector calculate_Rt(vector infections, int seeding_time,
                    vector gt_rev_pmf) {
  int t = num_elements(infections);
  int ot = t - seeding_time;
  vector[ot] R;
  vector[ot] sR;
  vector[ot] infectiousness = rep_vector(1e-5, ot);
  // calculate Rt using Cori et al. approach
  for (s in 1:ot) {
    infectiousness[s] += update_infectiousness(
      infections, gt_rev_pmf, seeding_time, s
    );
    R[s] = infections[s + seeding_time] / infectiousness[s];
  }
  return(R);
}
// Convert an estimate of Rt to growth
real[] R_to_growth(vector R, real gt_mean, real gt_sd) {
  int t = num_elements(R);
  real r[t];
  if (gt_sd > 0) {
    real k = pow(gt_sd / gt_mean, 2);
    for (s in 1:t) {
      r[s] = (pow(R[s], k) - 1) / (k * gt_mean);
    }
  } else {
    // limit as gt_sd -> 0
    for (s in 1:t) {
      r[s] = log(R[s]) / gt_mean;
    }
  }
  return(r);
}
// Calculate growth rate
real[] calculate_growth(vector infections, int seeding_time) {
  int t = num_elements(infections);
  int ot = t - seeding_time;
  vector[t] log_inf = log(infections); 
  vector[ot] growth = log_inf[(seeding_time + 1):t] - log_inf[seeding_time:(t - 1)];
  return(to_array_1d(growth));
}
