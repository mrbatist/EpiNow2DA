old_opts <- options()
options(mc.cores = ifelse(interactive(), 4, 1))

# get example case counts
reported_cases <- example_confirmed[1:60]

# set up example generation time
generation_time_fixed <- generation_time_opts(
 disease = "SARS-CoV-2", source = "ganyani", fixed = TRUE
)
generation_time_var <- generation_time_opts(
 disease = "SARS-CoV-2", source = "ganyani"
)
# set delays between infection and case report
incubation_period <- get_incubation_period(
 disease = "SARS-CoV-2", source = "lauer"
)
reporting_delay <- list(
  mean = convert_to_logmean(2, 1), mean_sd = 0,
  sd = convert_to_logsd(2, 1), sd_sd = 0, max = 10
)

def_fixed <- epinow(reported_cases,
  generation_time = generation_time_fixed,
  delays = delay_opts(incubation_period, reporting_delay, fixed = TRUE),
  rt = rt_opts(prior = list(mean = 2, sd = 0.1)),
  stan = stan_opts(control = list(adapt_delta = 0.95)),
  output = "timing"
)

def_one <- epinow(reported_cases,
  generation_time = generation_time_opts(),
  delays = delay_opts(incubation_period, reporting_delay, fixed = TRUE),
  rt = rt_opts(prior = list(mean = 2, sd = 0.1)),
  stan = stan_opts(control = list(adapt_delta = 0.95)),
  output = "timing"
)

def_one_nd <- epinow(reported_cases,
  generation_time = generation_time_opts(),
  rt = rt_opts(prior = list(mean = 2, sd = 0.1)),
  output = "timing"
)
