library(EpiNow2)
library(tictoc)

# get example case counts
reported_cases <- example_confirmed[1:60]
snapshot_cases <- example_confirmed[71:130]
bp_cases <- data.table::copy(reported_cases)
bp_cases <- bp_cases[, breakpoint := ifelse(date == as.Date("2020-03-16"), 1, 0)]

# set generation time
gt <- list(
  fixed = generation_time_opts(
    disease = "SARS-CoV-2", source = "ganyani", fixed = TRUE
  ),
  var = generation_time_opts(
    disease = "SARS-CoV-2", source = "ganyani"
  ),
  ar1 = generation_time_opts()
)

# set delays between infection and case report
incubation_period <- get_incubation_period(
 disease = "SARS-CoV-2", source = "lauer"
)
reporting_delay <- list(
  mean = convert_to_logmean(2, 1), mean_sd = 0.1,
  sd = convert_to_logsd(2, 1), sd_sd = 0.1, max = 10
)

delays <- list(
  var = delay_opts(incubation_period, reporting_delay),
  fixed = delay_opts(incubation_period, reporting_delay, fixed = TRUE),
  none = delay_opts()
)

trunc_dist <- trunc_opts(
  mean = convert_to_logmean(0.5, 0.5), mean_sd = 0.1,
  sd = convert_to_logsd(0.5, 0.5), sd_sd = 0.1,
  max = 3
)

methods <- list(
  nuts = list(
    stan = stan_opts(
      chains = 1, warmup = 100, samples = 400,
      control = list(adapt_delta = 0.95)
    )
  ),
  vb = list(
    stan = stan_opts(
      method = "vb"
    )
  )
)

test_scenarios <- list(
  default = list(
    reported_cases,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1))
  ),
  approximate_gp = list(
    reported_cases,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1)),
    gp = gp_opts(ls_min = 10, basis_prop = 0.1)
  ),
  truncation = list(
    reported_cases,
    truncation = trunc_dist,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1))
  ),
  backcalc = list(
    reported_cases,
    rt = NULL,
    backcalc = backcalc_opts(),
    obs = obs_opts(scale = list(mean = 0.4, sd = 0.05)),
    horizon = 0
  ),
  later_snapshot = list(
    snapshot_cases,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1))
  ),
  stationary_rt = list(
    reported_cases,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1), gp_on = "R0")
  ),
  fixed_rt = list(
    reported_cases,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1)),
    gp = NULL
  ),
  breakpoints_only = list(
    bp_cases,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1)),
    gp = NULL
  ),
  random_walk = list(
    reported_cases,
    rt = rt_opts(prior = list(mean = 2, sd = 0.1), rw = 7),
    gp = NULL
  )
)

rep <- 5

scenarios <- list()
timings <- list()
for (method in names(methods)) {
  scenarios[[method]] <- list()
  method_options <- methods[[method]]
  for (gt_scenario in names(gt)) {
    scenarios[[method]][[gt_scenario]] <- list()
    gt_options <- gt[[gt_scenario]]
    for (delay_scenario in names(delays)) {
      scenarios[[method]][[gt_scenario]][[delay_scenario]] <- list()
      delay_options <- delays[[delay_scenario]]
      for (test_scenario in names(test_scenarios)) {
        message(method, gt_scenario, delay_scenario, test_scenario)
        test_options <- test_scenarios[[test_scenario]]
        exec_times <- c()
        for (run in seq_len(rep)) {
          tic()
          scenarios[[method]][[gt_scenario]][[delay_scenario]][[test_scenario]] <-
            do.call(estimate_infections,
                    c(list(generation_time = gt_options),
                      list(delays = delay_options),
                      test_options,
                      method_options))
          run_time <- toc()
          exec_times[run] <- run_time$toc - run_time$tic
        }
        scenarios[[method]][[gt_scenario]][[delay_scenario]][[test_scenario]]$exec_times <-
          exec_times
      }
    }
  }
}

saveRDS(scenarios, "scenarios.rds")
