library(tidyverse)
library(brms)
library(survey)
library(weights)
library(anesrake)

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
# starting score/intercept is 100
# male has no benefit
# females get +20

# low ed is no difference
# high ed is +10

# being young is -10
# middle age is no difference
# old age is -10

true_value <- 100 + (.513 * 20) + (.4 * 10) + (.3 * -10)

#### Generate some samples with different levels of bias/skew for certain features ####
sex_vars <- c("female", "male")
age_vars <- c("young", "middle", "old")
education_vars <- c("low", "high")

data_maker <- function(sim = 1, bias = "none") {
  
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
                           sex_num = case_when(sex == "male" ~ 0,
                                               sex == "female" ~ 20),
                           age_num = case_when(age == "young" ~ -10,
                                               age == "middle" ~ 0,
                                               age == "old" ~ -10),
                           education_num = case_when(education == "low" ~ 0,
                                               education == "high" ~ 10),
                           error = rnorm(1000, 0, 12.5),
                           score = (100 + sex_num + age_num + education_num + error), # calculate their final score based on pop features
                           bias = bias,
                           sim = sim) %>% 
    # we have to ensure the factor levels exactly match the weighting targets
    mutate(sex = factor(sex,
                        levels = c("female", "male")),
           education = factor(education,
                              levels = c("low", "high")),
           age = factor(age,
                        levels = c("young", "middle", "old")))
  
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

# now we map the function to create lots of simulated data - 1000 per type of skew/bias
sim_data <- map2(.x = 1:4000,
                 .y = c(rep("none", 1000), rep("a", 1000), rep("b", 1000), rep("c", 1000)),
                 .f = data_maker)

# there are a few of the sims where the raking did not fully converge,
# but we will proceed for now

# the design effects and weight changes are clearly getting larger with increasingly unrepresentative data:
sim_data[251]
sim_data[501]
sim_data[751]

#### Regression modeling ####
# set up the basic brm models for updating (to avoid recompilation)

# imagining we have some vague information
# suggesting an expected score of around 105 in the population
reg_prior <- set_prior("normal(105, 10)", class = "Intercept")

test_data <- sim_data[[500]]

# this model will be the standard brm weights
base_reg_weight <- brm(formula = score | weights(weight) ~ 1,
                       family = gaussian(),
                       data = test_data,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       prior = reg_prior,
                       chains = 2,
                       cores = 2,
                       iter = 1300,
                       warmup = 300,
                       backend = "cmdstanr",
                       threads = threading(4),
                       seed = 1010)

# this model will be our weights, corrected according to the overall design effect
base_reg_weight_mod <- brm(formula = score | weights(weight_mod) ~ 1,
                       family = gaussian(),
                       data = test_data,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       prior = reg_prior,
                       chains = 2,
                       cores = 2,
                       iter = 1300,
                       warmup = 300,
                       backend = "cmdstanr",
                       threads = threading(4),
                       seed = 1010)

# now we set up a function that will run our regressions
# and pull out some basic features for assessment
run_regs <- function(data) {
  
  # run model with standard brm weighting
  standard_weight <- update(base_reg_weight,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior,
                            chains = 2,
                            cores = 2,
                            iter = 1000,
                            warmup = 250,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)
  
  standard_weight <- summary(standard_weight)
  
  # run model with 'corrected' i.e., reduced/penalised weights
  modified_weight <- update(base_reg_weight_mod,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior,
                            chains = 2,
                            cores = 2,
                            iter = 1000,
                            warmup = 250,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)
  
  modified_weight <- summary(modified_weight)
  
  # run as svyglm
  svy_weight_design <- svydesign(ids = ~1, data = data, weights = data$weight)
  svy_weight_output <- svyglm(formula = score ~ 1,
                              design = svy_weight_design)
  
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
  
  print(data[[1, 'sim']])
  
  return(output)
  
}

# map the regression function over all of our simulated data sets
t1_norm <- Sys.time()
regression_output <- map_dfr(.x = sim_data, .f = run_regs)
t2_norm <- Sys.time()

save(regression_output, file = "simulated 1k normal.RData")
# regression_output_2 <- regression_output %>% mutate(true_within_2 = true_value < upper & true_value > lower)

output_summary <- 
  regression_output %>% 
  group_by(type, bias) %>% 
  summarise(falls_within = sum(true_within) / 10,
            point_est = mean(estimate),
            avg_error = mean(error),
            avg_width = mean(width),
            avg_lower = mean(lower),
            avg_upper = mean(upper)) %>% 
  mutate(bias = factor(bias,
                       levels = c("none", "a", "b", "c"),
                       labels = c("None", "Moderate", "High", "Very high")),
         type = factor(type,
                       levels = c("svyglm", "modified brm", "standard brm")))

ggplot(output_summary) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = falls_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

regression_output_wide <- 
  regression_output %>% 
  pivot_wider(names_from = c(type), id_cols = sim, values_from = c(estimate, error, lower, upper, width))

ggplot(regression_output_wide) +
  geom_point(aes(x = `width_svyglm`, y = `width_standard brm`), alpha = .2, color = "#619ecd") +
  geom_point(aes(x = `width_svyglm`, y = `width_modified brm`), alpha = .2, color = "#61cd7d") +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Width of svyglm 95% CI",
       y = "Width of brm 95% ETI",
       subtitle = "Green dots show modified brms weights, blue dots standard brms weights")

save(sim_data, regression_output, file = "design effect calibration.RData")