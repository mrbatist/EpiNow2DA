functions {
#include functions/convolve.stan
#include functions/pmfs.stan
#include functions/gaussian_process.stan
#include functions/infections.stan
#include functions/observation_model.stan
#include functions/generated_quantities.stan
}

data {
  // dimensions
  int n; // number of samples
  int t; // unobserved time
  int seeding_time; // time period used for seeding and not observed
  int future_time; // fixed future time
  // Rt
#include data/simulation_rt.stan
  // delay from infection to report
#include data/simulation_delays.stan
  // observation model
#include data/simulation_observation_model.stan
}

transformed data {
  int delay_max_total = sum(delay_max) - num_elements(delay_max) + 1;
}

generated quantities {
  // generated quantities
  vector[n] infections[t]; //latent infections
  vector[n] reports[t - seeding_time]; // observed cases
  int imputed_reports[n, t - seeding_time];
  vector[n] r[t - seeding_time];
  vector[seeding_time] uobs_inf;
  for (i in 1:n) {
    // generate infections from Rt trace
    vector[gt_max[1]] gt_rev_pmf;
    vector[delay_max_total] delay_rev_pmf;

    gt_rev_pmf = reverse_mf(discretised_pmf(
      gt_mean[i, 1], gt_sd[i, 1], gt_max[1], gt_dist[1], 1
    ));
    delay_rev_pmf = combine_pmfs(
      to_vector([ 1 ]), delay_mean[i], delay_sd[i], delay_max, delay_dist,
      delay_max_total, 0, 1
    );

    uobs_inf = generate_seed(initial_infections[i], initial_growth[i], seeding_time);
     // generate infections from Rt trace
    infections[i] = renewal_model(R[i], uobs_inf, gt_rev_pmf, pop, future_time);
    // convolve from latent infections to mean of observations
    reports[i] = convolve_to_report(infections[i], delay_rev_pmf, seeding_time);
    // weekly reporting effect
    if (week_effect > 1) {
      reports[i] = day_of_week_effect(
        reports[i], day_of_week, to_vector(day_of_week_simplex[i])
      );
    }
    // scale observations
    if (obs_scale) {
      reports[i] = scale_obs(reports[i], frac_obs[i, 1]);
    }
   // simulate reported cases
   imputed_reports[i] = report_rng(reports[i], rep_phi[i], obs_dist);
   r[i] = calculate_growth(infections[i], seeding_time);
  }
}
