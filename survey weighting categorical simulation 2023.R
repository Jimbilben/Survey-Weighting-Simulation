library(tidyverse)
library(brms)
library(tidybayes)
library(survey)
library(weights)
library(anesrake)
library(svyVGAM)

#### True population proportions ####
# set the true proportions as targets for the raking from anesrake
target <- list(
  sex = wpct(c("female", "male"),
             c(0.513, 0.487)), 
  age = wpct(c("young", "middle", "old"),
             c(.15, .7, .15)),
  education = wpct(c("low", "high"),
                   c(.6, .4))
)

#### Set up our true effect ####
# A = 50
# B = 35
# C = 15

# male has no benefit
# females get A - 5, b - 5, c + 10

# being young is A - 5, B + 5
# middle age is no difference
# old age is A - 5 B + 5

# low ed is no difference
# high ed is A - 10, b same, c +10

true_a <- .5 + (.513 * -.05) + (.4 * -.1) + (.3 * -.05)
true_b <- .35 + (.513 * -.05) + (.4 * 0) + (.3 * .05)
true_c <- .15 + (.513 * .1) + (.4 * .1) + (.3 * 0) 
true_a + true_b + true_c

b_intercept <- log(true_b / true_a)
c_intercept <- log(true_c / true_a)

# exp(-.292) / sum(exp(-.292))

#### Generate some samples with different levels of bias/skew for certain features ####
sex_vars <- c("female", "male")
age_vars <- c("young", "middle", "old")
education_vars <- c("low", "high")

data_maker_cat <- function(sim = 1, bias = "none") {
  
  # sample matches population
  if(bias == "none") {
    sex_prob <- c(0.513, 0.487)
    age_prob <- c(.15, .7, .15)
    education_prob <- c(.6, .4)
  }
  # a through c are increasingly large amounts of bias/skew in the sampled data
  else if(bias == "a") {
    sex_prob <- c(.7, .3)
    age_prob <- c(.3, .4, .3)
    education_prob <- c(.4, .6)
  }
  else if(bias == "b") {
    sex_prob <- c(.8, .2)
    age_prob <- c(.35, .3, .35)
    education_prob <- c(.3, .7)
  }
  else if(bias == "c") {
    sex_prob <- c(.85, .15)
    age_prob <- c(.5, .25, .25)
    education_prob <- c(.225, .775)
  }
  
  generated_data <- tibble(sex = sample(sex_vars,
                                        1000,
                                        replace = TRUE,
                                        prob = sex_prob),
                           age = sample(age_vars,
                                        1000,
                                        replace = TRUE,
                                        prob = age_prob),
                           education = sample(education_vars,
                                              1000,
                                              replace = TRUE,
                                              prob = education_prob),
                           sex_prob_a = case_when(sex == "male" ~  0,
                                                  sex == "female" ~ -.05),
                           age_prob_a = case_when(age == "young" ~ -.05,
                                                  age == "middle" ~ 0,
                                                  age == "old" ~ -.05),
                           edu_prob_a = case_when(education == "low" ~ 0,
                                                  education == "high" ~ -.1),
                           sex_prob_b = case_when(sex == "male" ~ 0,
                                                  sex == "female" ~ -.05),
                           age_prob_b = case_when(age == "young" ~ .05,
                                                  age == "middle" ~ 0,
                                                  age == "old" ~ .05),
                           edu_prob_b = case_when(education == "low" ~ 0,
                                                  education == "high" ~ 0),
                           sex_prob_c = case_when(sex == "male" ~ 0,
                                                  sex == "female" ~ .1),
                           age_prob_c = case_when(age == "young" ~ 0,
                                                  age == "middle" ~ 0,
                                                  age == "old" ~ 0),
                           edu_prob_c = case_when(education == "low" ~ 0,
                                                  education == "high" ~ .1),
                           a_prob = .5 + sex_prob_a + age_prob_a + edu_prob_a,
                           b_prob = .35 + sex_prob_b + age_prob_b + edu_prob_b,
                           c_prob = .15 + sex_prob_c + age_prob_c + edu_prob_c,
                           bias = bias,
                           sim = sim) %>% 
    # we have to ensure the factor levels exactly match the weighting targets
    mutate(sex = factor(sex,
                        levels = c("female", "male")),
           education = factor(education,
                              levels = c("low", "high")),
           age = factor(age,
                        levels = c("young", "middle", "old")))
  
  sample_me <- function(a, b, c) {
    
    category <-
      sample(c("A", "B", "C"), 1, prob = c(a, b, c))
    
    return(category)
    
  }
  
  
  generated_data$outcome <- pmap_chr(.l = list(a = generated_data$a_prob,
                                               b = generated_data$b_prob,
                                               c = generated_data$c_prob),
                                     .f = sample_me)
  
  # generate the weights when bias is not none
  generated_data$id <- 1:1000
  
  if(bias != "none") {
    
    # in reality we'd probably want to ensure smaller caps
    # but this is for demonstrative purposes
    # to just make sure the raking generally works
    cap <- case_when(bias == "a" ~ 20,
                     bias == "b" ~ 25,
                     bias == "c" ~ 30)
    
    total_raking <- anesrake(inputter = target,
                             data = as.data.frame(generated_data),
                             caseid = generated_data$id,
                             cap = cap, # Maximum allowed weight per iteration.
                             choosemethod = "total", # How are parameters compared for selection?
                             type = "pctlim", # What selection criterion is used?
                             pctlim = 0.03, # Threshold for selection
                             nlim = 10,
                             maxit = 2000,
                             force1 = TRUE)
    
    generated_data$weight <- total_raking$weightvec
    generated_data$deff <- generaldesigneffect(generated_data$weight)
    generated_data$weight_mod <- generated_data$weight / generated_data$deff
    
  }
  # when the data matches the population we will avoid weighting
  else if(bias == "none") {
    generated_data$weight <- 1
    generated_data$deff <- 1
    generated_data$weight_mod <- 1
  }
  
  return(generated_data)
  
}

# now we map the function to create lots of simulated data - 250 per type of skew/bias
sim_data_cat <- map2(.x = 1:4000,
                     .y = c(rep("none", 1000), rep("a", 1000), rep("b", 1000), rep("c", 1000)),
                     .f = data_maker_cat)

# there are a few of the sims where the raking did not fully converge,
# but we will proceed for now

# the design effects and weight changes are clearly getting larger with increasingly unrepresentative data:
sim_data[251]
sim_data[501]
sim_data[751]

#### Regression modeling ####
# set up the basic brm models for updating (to avoid recompilation)

# get_prior(formula = outcome | weights(weight) ~ 1,
#           data = test_data_cat,
#           family = "categorical")
reg_prior_cat <- c(set_prior("normal(0, 1)", class = "Intercept", dpar = "muB"),
                   set_prior("normal(0, 1)", class = "Intercept", dpar = "muC"))

test_data_cat <- sim_data_cat[[1]]

# this model will be the standard brm weights
base_reg_weight_cat <- brm(formula = outcome | weights(weight) ~ 1,
                             family = categorical(),
                             data = test_data_cat,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             prior = reg_prior_cat,
                             chains = 2,
                             cores = 2,
                             iter = 1200,
                             warmup = 200,
                             backend = "cmdstanr",
                             threads = threading(4),
                             seed = 1010)

svy_weight_design_test <- svydesign(ids = ~1, data = test_data_cat, weights = test_data_cat$weight)
svy_weight_output_test <- svy_vglm(formula = outcome ~ 1,
                            family = multinomial(refLevel = "A"),
                            design = svy_weight_design_test)

svy_weight_output_test_summary <-
tidy.svyVGAM(svy_weight_output_test, exponentiate = FALSE, conf.int = TRUE)

output_cat <- tibble(b_estimate = c(mean(standard_weight_draws[,1]), mean(modified_weight_draws[,1]), svy_weight_outputcat_summary$estimate[[1]]),
                      c_estimate = c(mean(standard_weight_draws[,2]), mean(modified_weight_draws[,2]), svy_weight_outputcat_summary$estimate[[2]]),
                      b_lower = c(hdi(standard_weight_draws[,1])[1], hdi(modified_weight_draws[,1])[1], svy_weight_outputcat_summary$conf.low[[1]]),
                      c_lower = c(hdi(standard_weight_draws[,1])[1], hdi(modified_weight_draws[,1])[1], svy_weight_outputcat_summary$conf.low[[2]]),
                      b_upper = c(hdi(standard_weight_draws[,1])[2], hdi(modified_weight_draws[,1])[2], svy_weight_outputcat_summary$conf.high[[1]]),
                      c_upper = c(hdi(standard_weight_draws[,1])[2], hdi(modified_weight_draws[,1])[2], svy_weight_outputcat_summary$conf.high[[2]]),
                      bias = data[[1, 'bias']],
                      deff = data[[1, 'deff']],
                      type = c("standard brm", "modified brm", "svyglm"),
                      a_width = abs(b_lower - b_upper),
                      b_width = abs(c_lower - c_upper),
                      b_true_within = b_intercept < b_upper & b_intercept > b_lower, # does the true value fall within our 95% bounds?
                      c_true_within = c_intercept < c_upper & c_intercept > c_lower, # does the true value fall within our 95% bounds?
                      sim = data[[1, 'sim']])

# this model will be our weights, corrected according to the overall design effect
base_reg_weight_mod_cat <- brm(formula = outcome | weights(weight_mod) ~ 1,
                               family = categorical(),
                               data = test_data_cat,
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               prior = reg_prior_cat,
                               chains = 2,
                               cores = 2,
                               iter = 1200,
                               warmup = 200,
                               backend = "cmdstanr",
                               threads = threading(4),
                               seed = 1010)

# now we set up a function that will run our regressions
# and pull out some basic features for assessment

tidy.svyVGAM <- function(
    x, 
    conf.int = FALSE, 
    conf.level = 0.95,
    exponentiate = FALSE, 
    ...
){
  # Replace `summary(x)$coefficients` with `summary(x)$coeftable`
  ret <- as_tibble(summary(x)$coeftable, rownames = "term")
  
  # All of this stays the same:
  colnames(ret) <- c("term", "estimate", "std.error", "statistic", "p.value")
  coefs <- tibble::enframe(stats::coef(x), name = "term", value = "estimate")
  ret <- left_join(coefs, ret, by = c("term", "estimate"))
  if (conf.int){
    ci <- broom:::broom_confint_terms(x, level = conf.level, ...)
    ret <- dplyr::left_join(ret, ci, by = "term")
  }
  if (exponentiate){ret <- broom:::exponentiate(ret)}
  
  # This part only works for the multinomial case, and only if your covariates
  # have no ":" in their names - NOT FOR GENERAL USE
  ret %>% 
    separate(term, into = c("term", "y.level"), sep = ":") %>% 
    arrange(y.level) %>% 
    relocate(y.level, .before = term)
}

run_regs_cat <- function(data) {
  
  # run model with standard brm weighting
  standard_weight <- update(base_reg_weight_cat,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior_cat,
                            chains = 2,
                            cores = 2,
                            iter = 1200,
                            warmup = 200,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)
  
  standard_weight <- as_draws_df(standard_weight) %>% as_tibble()
  
  # run model with 'corrected' i.e., reduced/penalised weights
  modified_weight <- update(base_reg_weight_mod_cat,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior_cat,
                            chains = 2,
                            cores = 2,
                            iter = 1200,
                            warmup = 200,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)
  
  modified_weight <- as_draws_df(modified_weight) %>% as_tibble()
  
  # run as svyglm
  svy_weight_design <- svydesign(ids = ~1, data = data, weights = data$weight)
  svy_weight_output <- svy_vglm(formula = outcome ~ 1,
                                family = multinomial(refLevel = "A"),
                                design = svy_weight_design)
  
  svy_weight_output_summary <-
    tidy.svyVGAM(svy_weight_output, exponentiate = FALSE, conf.int = TRUE)
  
  # extract desired information:
  output <- tibble(b_estimate = c(mean(pull(standard_weight[,1])), mean(pull(modified_weight[,1])), svy_weight_output_summary$estimate[[1]]),
                   c_estimate = c(mean(pull(standard_weight[,2])), mean(pull(modified_weight[,2])), svy_weight_output_summary$estimate[[2]]),
                   b_lower = c(hdi(pull(standard_weight[,1]))[1], hdi(pull(modified_weight[,1]))[1], svy_weight_output_summary$conf.low[[1]]),
                   c_lower = c(hdi(pull(standard_weight[,2]))[1], hdi(pull(modified_weight[,2]))[1], svy_weight_output_summary$conf.low[[2]]),
                   b_upper = c(hdi(pull(standard_weight[,1]))[2], hdi(pull(modified_weight[,1]))[2], svy_weight_output_summary$conf.high[[1]]),
                   c_upper = c(hdi(pull(standard_weight[,2]))[2], hdi(pull(modified_weight[,2]))[2], svy_weight_output_summary$conf.high[[2]]),
                   bias = data[[1, 'bias']],
                   deff = data[[1, 'deff']],
                   type = c("standard brm", "modified brm", "svyglm"),
                   b_width = abs(b_lower - b_upper),
                   c_width = abs(c_lower - c_upper),
                   b_true_within = b_intercept < b_upper & b_intercept > b_lower, # does the true value fall within our 95% bounds?
                   c_true_within = c_intercept < c_upper & c_intercept > c_lower, # does the true value fall within our 95% bounds?
                   sim = data[[1, 'sim']])
  
  print(data[[1, 'sim']])
  
  return(output)
  
}

# map the regression function over all of our simulated data sets
t1 <- Sys.time()
regression_output_cat <- map_dfr(.x = sim_data_cat, .f = run_regs_cat)
t2 <- Sys.time()

regression_output_cat <-
  regression_output_cat %>% 
  mutate(both_within = case_when(b_true_within == TRUE & c_true_within == TRUE ~ TRUE,
                                 TRUE ~ FALSE))
save(regression_output_cat,
     file = "simulated 1k categorical.RData")

output_summary_cat <- 
  regression_output_cat %>% 
  group_by(type, bias) %>% 
  summarise(b_falls_within = sum(b_true_within) / 10,
            b_point_est = mean(b_estimate),
            b_avg_width = mean(b_width),
            b_avg_lower = mean(b_lower),
            b_avg_upper = mean(b_upper),
            c_falls_within = sum(c_true_within) / 10,
            c_point_est = mean(c_estimate),
            c_avg_width = mean(c_width),
            c_avg_lower = mean(c_lower),
            c_avg_upper = mean(c_upper),
            both_fall_within = sum(both_within) / 10) %>% 
  mutate(bias = factor(bias,
                       levels = c("none", "a", "b", "c"),
                       labels = c("None", "Moderate", "High", "Very high")),
         type = factor(type,
                       levels = c("svyglm", "modified brm", "standard brm")))

ggplot(output_summary_cat) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = b_falls_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

ggplot(output_summary_cat) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = c_falls_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

ggplot(output_summary_cat) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = both_fall_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

regression_output_wide_cat <-
  regression_output_cat %>% 
  pivot_wider(names_from = c(type), id_cols = sim, values_from = c(perc_estimate, perc_lower, perc_upper, perc_width))

ggplot(regression_output_wide_cat) +
  geom_point(aes(x = `perc_width_svyglm` * 100, y = `perc_width_standard brm` * 100), alpha = .2, color = "#619ecd") +
  geom_point(aes(x = `perc_width_svyglm` * 100, y = `perc_width_modified brm` * 100), alpha = .2, color = "#61cd7d") +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Width of svyglm 95% CI in %age points",
       y = "Width of brm 95% ETI in %age points",
       subtitle = "Green dots show modified brms weights, blue dots standard brms weights")

save(sim_data_cat, regression_output_cat, file = "design effect calibration catial.RData")