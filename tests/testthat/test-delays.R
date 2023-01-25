test_stan_data <- function(generation_time = dist_spec(mean = 1),
                           delays = dist_spec(),
                           truncation = dist_spec(),
                           params = c()) {
  data <- create_stan_data(
    reported_cases = example_confirmed,
    generation_time = generation_time,
    delays = delays,
    truncation = truncation,
    rt = rt_opts(),
    gp = gp_opts(),
    obs = obs_opts(),
    backcalc = backcalc_opts(),
    shifted_cases = NULL,
    horizon = 7,
    seeding_time = 0
  )
  return(unlist(unname(data[params])))
}

test_that("generation times can be specified in different ways", {
  gt_params <-
    c("gt_mean_mean", "gt_mean_sd", "gt_sd_mean", "gt_sd_sd", "gt_max")
  expect_equal(
    test_stan_data(params = gt_params),
    c(1, 0, 0, 0, 1)
  )
  expect_equal(
    test_stan_data(
      generation_time = dist_spec(mean = 3),
      params = gt_params
    ),
    c(3, 0, 0, 0, 3)
  )
  expect_equal(
    test_stan_data(
      generation_time = dist_spec(mean = 3, sd = 1, max = 5),
      params = gt_params
    ),
    c(3, 0, 1, 0, 5)
  )
  expect_equal(
    round(test_stan_data(
      generation_time = get_generation_time(
        disease = "SARS-CoV-2", source = "ganyani",
        max = 10, fixed = TRUE
      ),
      params = gt_params
    ), digits = 2),
    c(3.64, 0, 3.08, 0, 10)
  )
})

test_that("delay parameters can be specified in different ways", {
  delay_params <-
    c(
      "delay_mean_mean", "delay_mean_sd", "delay_sd_mean", "delay_sd_sd",
      "delay_max"
    )
  expect_equal(
    test_stan_data(
      delays = dist_spec(mean = 3),
      params = delay_params
    ),
    c(3, 0, 0, 0, 3)
  )
  expect_equal(
    test_stan_data(
      delays = dist_spec(mean = 3, sd = 1, max = 5),
      params = delay_params
    ),
    c(3, 0, 1, 0, 5)
  )
})

test_that("truncation parameters can be specified in different ways", {
  trunc_params <-
    c(
      "trunc_mean_mean", "trunc_mean_sd", "trunc_sd_mean", "trunc_sd_sd",
      "trunc_max"
    )
  expect_equal(
    test_stan_data(
      truncation = dist_spec(mean = 3, sd = 1, max = 5),
      params = trunc_params
    ),
    c(3, 0, 1, 0, 5)
  )
})

test_that("contradictory generation times are caught", {
  expect_error(dist_spec(mean = 3.5), "must be an integer")
  expect_error(
    dist_spec(mean = 3, mean_sd = 1),
    "must be 0"
  )
})

test_that("contradictory delays are caught", {
  expect_error(
    test_stan_data(delays = dist_spec(mean = 3.5)),
    "must be an integer"
  )
  expect_error(
    test_stan_data(delays = dist_spec(mean = 3, mean_sd = 1)),
    "must be 0"
  )
})
