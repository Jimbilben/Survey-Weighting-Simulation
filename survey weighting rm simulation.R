library(tidyverse)
library(brms)
library(survey)
library(weights)
library(anesrake)

#### True population proportions ####
# set the true proportions as targets for the raking from anesrake
target_rm <- list(
  sex = wpct(c("female", "male"),
             c(0.513, 0.487)), 
  student = wpct(c("no", "yes"),
                   c(.67, .33))
)

#### Set up our true effect ####
# starting score/intercept is 100
# male has no benefit
# females get +20

# low ed is no difference
# high ed is +10

# being young is -10
# middle age is no difference
# old age is -10
true_treatment <- .3 + (.513 * .1) + (.33 * .1)

#### Generate some samples with different levels of bias/skew for certain features ####
sex_vars <- c("female", "male")
student_vars <- c("no", "yes")

data_maker_rm <- function(sim = 1, bias = "none") {
  
  # sample matches population
  if(bias == "none") {
    sex_prob <- c(0.513, 0.487)
    student_prob <- c(.67, .33)
  }
  # a through c are increasingly large amounts of bias/skew in the sampled data
  else if(bias == "a") {
    sex_prob <- c(0.4, 0.6)
    student_prob <- c(.5, .5)
  }
  else if(bias == "b") {
    sex_prob <- c(.3, .6)
    student_prob <- c(.3, .7)
  }
  else if(bias == "c") {
    sex_prob <- c(0.2, .8)
    student_prob <- c(.15, .85)
  }
  
  generated_data <- tibble(sex = sample(sex_vars,
                                        1000,
                                        replace = TRUE,
                                        prob = sex_prob),
                           student = sample(student_vars,
                                            1000,
                                            replace = TRUE,
                                            prob = student_prob),
                           rm_error = rnorm(1000, 0, .5),
                           ppn = 1:1000) %>% 
    mutate(sex = factor(sex,
                        levels = c("female", "male")),
           student = factor(student,
                            levels = c("no", "yes")))
  
  generated_data_long <- bind_rows(generated_data,
                                   generated_data) %>% 
    mutate(error = rnorm(2000, 0, .5),
           time = c(rep("time_a", 1000), rep("time_b", 1000)),
           sex_num = case_when(time == "time_a" ~ 0,
                               sex == "male" ~ 0,
                               sex == "female" ~ .1),
           student_num = case_when(time == "time_a" ~ 0,
                                   student == "no" ~ 0,
                                   student == "yes" ~ .1),
           treatment_num = case_when(time == "time_a" ~ 0,
                                     time == "time_b" ~ .3),
           score = (sex_num + student_num + treatment_num + rm_error + error), # calculate their final score based on pop features
           bias = bias,
           sim = sim) %>% 
    # we have to ensure the factor levels exactly match the weighting targets
    mutate(sex = factor(sex,
                        levels = c("female", "male")),
           student = factor(student,
                            levels = c("no", "yes")))
  
  # generate the weights when bias is not none
  generated_data$id <- 1:1000
  
  if(bias != "none") {
    
    # in reality we'd probably want to ensure smaller caps
    # but this is for demonstrative purposes
    # to just make sure the raking generally works
    cap <- case_when(bias == "a" ~ 20,
                     bias == "b" ~ 25,
                     bias == "c" ~ 30)
    
    total_raking <- anesrake(inputter = target_rm,
                             data = as.data.frame(generated_data),
                             caseid = generated_data$id,
                             cap = cap, # Maximum allowed weight per iteration.
                             choosemethod = "total", # How are parameters compared for selection?
                             type = "pctlim", # What selection criterion is used?
                             pctlim = 0.03, # Threshold for selection
                             nlim = 10,
                             maxit = 2000,
                             force1 = TRUE)
    
    generated_data_long$weight <- c(total_raking$weightvec, total_raking$weightvec)
    generated_data_long$deff <- generaldesigneffect(generated_data_long$weight)
    generated_data_long$weight_mod <- generated_data_long$weight / generated_data_long$deff
    generated_data_long$weight_modsqrt <- generated_data_long$weight / sqrt(generated_data_long$deff)
    
  }
  # when the data matches the population we will avoid weighting
  else if(bias == "none") {
    generated_data_long$weight <- 1
    generated_data_long$deff <- 1
    generated_data_long$weight_mod <- 1
    generated_data_long$weight_modsqrt <- 1
  }

  return(generated_data_long)
  
}

# now we map the function to create lots of simulated data - 250 per type of skew/bias
sim_data_rm <- map2(.x = 1:4000,
                 .y = c(rep("none", 1000), rep("a", 1000), rep("b", 1000), rep("c", 1000)),
                 .f = data_maker_rm)

# there are a few of the sims where the raking did not fully converge,
# but we will proceed for now

# the design effects and weight changes are clearly getting larger with increasingly unrepresentative data:
sim_data_rm[251]
sim_data_rm[501]
sim_data_rm[751]

#### Regression modeling ####
# set up the basic brm models for updating (to avoid recompilation)

# imagining we have some vague information
# suggesting an expected score of around 105 in the population
reg_prior_rm <- c(set_prior("normal(0, 1)", class = "Intercept"),
                  set_prior("normal(0, .5)", class = "b"),
                  set_prior("normal(0, .5)", class = "sd"))

test_data_rm <- sim_data_rm[[500]]

# this model will be the standard brm weights
base_reg_weight_rm <- brm(formula = score | weights(weight) ~ 1 + time + (1 | ppn),
                       family = gaussian(),
                       data = test_data_rm,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       prior = reg_prior_rm,
                       chains = 2,
                       cores = 2,
                       iter = 1200,
                       warmup = 200,
                       backend = "cmdstanr",
                       threads = threading(4),
                       seed = 1010)

# this model will be our weights, corrected according to the overall design effect
base_reg_weight_mod_rm <- brm(formula = score | weights(weight_mod) ~ 1 + time + (1 | ppn),
                           family = gaussian(),
                           data = test_data_rm,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           prior = reg_prior_rm,
                           chains = 2,
                           cores = 2,
                           iter = 1200,
                           warmup = 200,
                           backend = "cmdstanr",
                           threads = threading(4),
                           seed = 1010)

base_reg_weight_modsqrt_rm <- brm(formula = score | weights(weight_modsqrt) ~ 1 + time + (1 | ppn),
                              family = gaussian(),
                              data = test_data_rm,
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              prior = reg_prior_rm,
                              chains = 2,
                              cores = 2,
                              iter = 1200,
                              warmup = 200,
                              backend = "cmdstanr",
                              threads = threading(4),
                              seed = 1010)

# now we set up a function that will run our regressions
# and pull out some basic features for assessment
run_regs_rm <- function(data) {
  
  #run model with standard brm weighting
  standard_weight <- update(base_reg_weight,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior,
                            chains = 2,
                            cores = 2,
                            iter = 900,
                            warmup = 250,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)

  standard_weight <- summary(standard_weight)
  
  # run model with 'corrected' i.e., reduced/penalised weights
  modified_weight <- update(base_reg_weight_mod_rm,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior_rm,
                            chains = 1,
                            # cores = 2,
                            iter = 900,
                            warmup = 250,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)

  modified_weight <- summary(modified_weight)
  
  # modified_weight <- update(base_reg_weight_modsqrt_rm,
  #                           newdata = data,
  #                           control = list(adapt_delta = 0.99, max_treedepth = 15),
  #                           prior = reg_prior_rm,
  #                           chains = 1,
  #                           # cores = 2,
  #                           iter = 800,
  #                           warmup = 200,
  #                           backend = "cmdstanr",
  #                           threads = threading(4),
  #                           seed = 1010)
  
  # run as svyglm
  # svy_weight_design <- svydesign(ids = ~1, data = data, weights = data$weight)
  # svy_weight_output <- svyglm(formula = score ~ 1,
  #                             design = svy_weight_design)
  
  # extract desired information:
  output <- tibble(estimate = c(standard_weight$fixed[[1]], modified_weight$fixed[[1]], summary(svy_weight_output)$coefficients[[1]]),
                   error = c(standard_weight$fixed[[2]], modified_weight$fixed[[2]], summary(svy_weight_output)$coefficients[[2]]),
                   lower = c(standard_weight$fixed[[3]], modified_weight$fixed[[3]], confint(svy_weight_output)[[1]]),
                   upper = c(standard_weight$fixed[[4]], modified_weight$fixed[[4]], confint(svy_weight_output)[[2]]),
                   bias = data[[1, 'bias']],
                   deff = data[[1, 'deff']],
                   type = c("standard brm", "modified brm", "svyglm"),
                   width = abs(lower - upper),
                   true_within = true_value < upper & true_value > lower, # does the true value fall within our 95% bounds?
                   sim = data[[1, 'sim']])
  
  # output <- tibble(estimate = c(modified_weight$fixed[[1]]),
  #                  error = c(modified_weight$fixed[[2]]),
  #                  lower = c(modified_weight$fixed[[3]]),
  #                  upper = c(modified_weight$fixed[[4]]),
  #                  bias = data[[1, 'bias']],
  #                  deff = data[[1, 'deff']],
  #                  type = c("modified brm"),
  #                  width = abs(lower - upper),
  #                  true_within = true_treatment < upper & true_treatment > lower, # does the true value fall within our 95% bounds?
  #                  sim = data[[1, 'sim']])
  
  print(data[[1, 'sim']])
  
  return(output)
  
}

# map the regression function over all of our simulated data sets
run_regs_rm(test_data_rm)
regression_output_rm <- map_dfr(.x = sim_data_rm, .f = run_regs_rm)

regression_output_rm_coef <- regression_output_rm[seq(2, 2400, 2) , ]

regression_output_rm_coef <-
regression_output_rm_coef %>% mutate(est_type = case_when(upper < true_treatment ~ "Under",
                                                          lower > true_treatment ~ "Over",
                                                          TRUE ~ "correct"))

output_summary_rm2 <- regression_output_rm_coef %>% 
  group_by(bias, est_type) %>% 
  summarise(prop = n()/ 3)

output_summary_rm <- regression_output_rm %>% 
  group_by(bias) %>% 
  summarise(falls_within = sum(true_within) / 3,
            point_est = mean(estimate),
            avg_error = mean(error),
            avg_width = mean(width),
            avg_lower = mean(lower),
            avg_upper = mean(upper)) %>% 
  mutate(bias = factor(bias,
                       levels = c("none", "a", "b", "c"),
                       labels = c("None", "Moderate", "High", "Very high")))

ggplot(output_summary_rm) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = falls_within), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")


regression_output_rm_sqrt <- map_dfr(.x = sim_data_rm, .f = run_regs_rm)

regression_output_rm_sqrt_coef <- regression_output_rm_sqrt[seq(2, 2400, 2) , ]

regression_output_rm_sqrt_coef <-
  regression_output_rm_sqrt_coef %>% mutate(est_type = case_when(upper < true_treatment ~ "Under",
                                                            lower > true_treatment ~ "Over",
                                                            TRUE ~ "correct"))

output_summary_rm_sqrt2 <- regression_output_rm_sqrt_coef %>% 
  group_by(bias, est_type) %>% 
  summarise(prop = n()/ 3)

output_summary_rm_sqrt <- regression_output_rm_sqrt %>% 
  group_by(bias) %>% 
  summarise(falls_within = sum(true_within) / 3,
            point_est = mean(estimate),
            avg_error = mean(error),
            avg_width = mean(width),
            avg_lower = mean(lower),
            avg_upper = mean(upper)) %>% 
  mutate(bias = factor(bias,
                       levels = c("none", "a", "b", "c"),
                       labels = c("None", "Moderate", "High", "Very high")))

ggplot(output_summary_rm_sqrt) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = falls_within), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")
