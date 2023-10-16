remotes::install_github("carlganz/svrepmisc")
library(svrepmisc)
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
# 1 = 7.5
# 2 = 10
# 3 = 20
# 4 = 35
# 5 = 20

# cutpoints <- -1.75, -.75, 0, 1, 1.5
test <-
tibble(normal = rnorm(100000000, -.2, 1),
       sex = sample(c("female", "male"), 100000000, replace = TRUE, prob = c(0.513, 0.487)),
       age = sample(c("young", "middle", "old"), 100000000, replace = TRUE, prob = c(.15, .7, .15)),
       education = sample(c("low", "high"), 100000000, replace = TRUE, prob = c(.6, .4)),
       sex_effect = case_when(sex == "female" ~ .3,
                              TRUE ~ 0),
       age_effect = case_when(age == "young" ~ .3,
                              age == "old" ~ .15,
                              TRUE ~ 0),
       edu_effect = case_when(education == "high" ~ .3,
                              TRUE ~ 0),
       score = normal + sex_effect + age_effect + edu_effect,
       ordinal = case_when(score < -1.25 ~ 1,
                           score < -.7 ~ 2,
                           score < 0 ~ 3,
                           score < 1 ~ 4,
                           TRUE ~ 5))

count_100mil <-
  test %>% 
  count_data(ordinal)

a_perc <- count_100mil[[1, 3]]/100
b_perc <- count_100mil[[2, 3]]/100
c_perc <- count_100mil[[3, 3]]/100
d_perc <- count_100mil[[4, 3]]/100
e_perc <- count_100mil[[5, 3]]/100

test_small <-
  test %>% 
  sample_n(1000)

test_small %>% 
  count_data(ordinal)

rm(test)

ggplot(test_small) +
  geom_histogram(aes(x = ordinal), bins = 5)



#### Generate some samples with different levels of bias/skew for certain features ####
sex_vars <- c("female", "male")
age_vars <- c("young", "middle", "old")
education_vars <- c("low", "high")

data_maker_ord <- function(sim = 1, bias = "none") {
  
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
  
  generated_data <- tibble(id = as.factor(1:1000),
                           sex = sample(sex_vars,
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
                           normal = rnorm(1000, -.2, 1),
                           sex_effect = case_when(sex == "female" ~ .3,
                                                  TRUE ~ 0),
                           age_effect = case_when(age == "young" ~ .3,
                                                  age == "old" ~ .15,
                                                  TRUE ~ 0),
                           edu_effect = case_when(education == "high" ~ .3,
                                                  TRUE ~ 0),
                           score = normal + sex_effect + age_effect + edu_effect,
                           outcome = case_when(score < -1.25 ~ 1,
                                               score < -.7 ~ 2,
                                               score < 0 ~ 3,
                                               score < 1 ~ 4,
                                               TRUE ~ 5),
                           bias = bias,
                           sim = sim) %>% 
    # we have to ensure the factor levels exactly match the weighting targets
    mutate(sex = factor(sex,
                        levels = c("female", "male")),
           education = factor(education,
                              levels = c("low", "high")),
           age = factor(age,
                        levels = c("young", "middle", "old")))
  
  if(bias != "none") {
    
    # in reality we'd probably want to ensure smaller caps
    # but this is for demonstrative purposes
    # to just make sure the raking generally works
    cap <- case_when(bias == "a" ~ 25,
                     bias == "b" ~ 30,
                     bias == "c" ~ 35)
    
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
sim_data_ord <- map2(.x = 1:4000,
                     .y = c(rep("none", 1000), rep("a", 1000), rep("b", 1000), rep("c", 1000)),
                     .f = data_maker_ord)

# there are a few of the sims where the raking did not fully converge,
# but we will proceed for now

# the design effects and weight changes are clearly getting larger with increasingly unrepresentative data:
sim_data_ord[251]
sim_data_ord[501]
sim_data_ord[751]

#### Regression modeling ####
# set up the basic brm models for updating (to avoid recompilation)

# get_prior(formula = outcome | weights(weight) ~ 1,
#           data = test_data_cat,
#           family = "categorical")
reg_prior_ord <- c(set_prior("normal(0, 1.5)", class = "Intercept"))

test_data_ord <- sim_data_ord[[1]]

# this model will be the standard brm weights
base_reg_weight_ord <- brm(formula = outcome | weights(weight) ~ 1,
                           family = brms::cumulative("probit"),
                           data = test_data_ord,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           prior = reg_prior_ord,
                           chains = 2,
                           cores = 2,
                           iter = 1200,
                           warmup = 200,
                           init = 0,
                           backend = "cmdstanr",
                           threads = threading(4),
                           seed = 1010)

# this model will be our weights, corrected according to the overall design effect
base_reg_weight_mod_ord <- brm(formula = outcome | weights(weight_mod) ~ 1,
                               family = categorical(),
                               data = test_data_ord,
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               prior = reg_prior_ord,
                               chains = 2,
                               cores = 2,
                               iter = 1200,
                               warmup = 200,
                               init = 0,
                               backend = "cmdstanr",
                               threads = threading(4),
                               seed = 1010)

# now we set up a function that will run our regressions
# and pull out some basic features for assessment

run_regs_ord <- function(data) {
  
  # run model with standard brm weighting
  standard_weight <- update(base_reg_weight_ord,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior_ord,
                            chains = 2,
                            cores = 2,
                            iter = 1200,
                            warmup = 200,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)
  
  standard_weight <- 
    posterior_epred(standard_weight,
                    ndraws = 2000,
                    newdata = data[1, ])
  
  # run model with 'corrected' i.e., reduced/penalised weights
  modified_weight <- update(base_reg_weight_mod_ord,
                            newdata = data,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            prior = reg_prior_ord,
                            chains = 2,
                            cores = 2,
                            iter = 1200,
                            warmup = 200,
                            backend = "cmdstanr",
                            threads = threading(4),
                            seed = 1010)
  
  modified_weight <- 
    posterior_epred(modified_weight,
                    ndraws = 2000,
                    newdata = data[1, ])

  output <- tibble(a_estimate = c(mean(standard_weight[,,1]), mean(modified_weight[,,1])),
                       b_estimate = c(mean(standard_weight[,,2]), mean(modified_weight[,,2])),
                       c_estimate = c(mean(standard_weight[,,3]), mean(modified_weight[,,3])),
                       d_estimate = c(mean(standard_weight[,,4]), mean(modified_weight[,,4])),
                       e_estimate = c(mean(standard_weight[,,5]), mean(modified_weight[,,5])),
                       a_lower = c(hdi(standard_weight[,,1])[1], hdi(modified_weight[,,1])[1]),
                       b_lower = c(hdi(standard_weight[,,2])[1], hdi(modified_weight[,,2])[1]),
                       c_lower = c(hdi(standard_weight[,,3])[1], hdi(modified_weight[,,3])[1]),
                       d_lower = c(hdi(standard_weight[,,4])[1], hdi(modified_weight[,,4])[1]),
                       e_lower = c(hdi(standard_weight[,,5])[1], hdi(modified_weight[,,5])[1]),
                       a_upper = c(hdi(standard_weight[,,1])[2], hdi(modified_weight[,,1])[2]),
                       b_upper = c(hdi(standard_weight[,,2])[2], hdi(modified_weight[,,2])[2]),
                       c_upper = c(hdi(standard_weight[,,3])[2], hdi(modified_weight[,,3])[2]),
                       d_upper = c(hdi(standard_weight[,,4])[2], hdi(modified_weight[,,4])[2]),
                       e_upper = c(hdi(standard_weight[,,5])[2], hdi(modified_weight[,,5])[2]),
                       bias = data[[1, 'bias']],
                       deff = data[[1, 'deff']],
                       type = c("standard brm", "modified brm"),
                       a_width = abs(a_lower - a_upper),
                       b_width = abs(b_lower - b_upper),
                       c_width = abs(c_lower - c_upper),
                       d_width = abs(d_lower - d_upper),
                       e_width = abs(e_lower - e_upper),
                       
                       a_true_within = a_perc < a_upper & a_perc > a_lower, # does the true value fall within our 95% bounds?
                       b_true_within = b_perc < b_upper & b_perc > b_lower, # does the true value fall within our 95% bounds?
                       c_true_within = c_perc < c_upper & c_perc > c_lower, # does the true value fall within our 95% bounds?
                       d_true_within = d_perc < d_upper & d_perc > d_lower, # does the true value fall within our 95% bounds?
                       e_true_within = e_perc < e_upper & e_perc > e_lower, # does the true value fall within our 95% bounds?
                       
                       sim = data[[1, 'sim']])
  
  print(data[[1, 'sim']])
  
  return(output)
  
}

# map the regression function over all of our simulated data sets
t1 <- Sys.time()
regression_output_ord <- map_dfr(.x = sim_data_ord, .f = run_regs_ord)
t2 <- Sys.time()

save(regression_output_ord,
     file = "simulated 1k ordinal.RData")

regression_output_ord <-
  regression_output_ord %>% 
  mutate(all_within = case_when(a_true_within == TRUE & b_true_within == TRUE & c_true_within == TRUE & d_true_within == TRUE & e_true_within == TRUE ~ TRUE,
                                 TRUE ~ FALSE))

output_summary_ord <- 
  regression_output_ord %>% 
  group_by(type, bias) %>% 
  summarise(a_falls_within = sum(a_true_within) / 10,
            a_point_est = mean(a_estimate),
            a_avg_width = mean(a_width),
            a_avg_lower = mean(a_lower),
            a_avg_upper = mean(a_upper),
            b_falls_within = sum(b_true_within) / 10,
            b_point_est = mean(b_estimate),
            b_avg_width = mean(b_width),
            b_avg_lower = mean(b_lower),
            b_avg_upper = mean(b_upper),
            c_falls_within = sum(c_true_within) / 10,
            c_point_est = mean(c_estimate),
            c_avg_width = mean(c_width),
            c_avg_lower = mean(c_lower),
            c_avg_upper = mean(c_upper),
            all_fall_within = sum(all_within) / 10,
            d_falls_within = sum(d_true_within) / 10,
            d_point_est = mean(d_estimate),
            d_avg_width = mean(d_width),
            d_avg_lower = mean(d_lower),
            d_avg_upper = mean(d_upper),
            e_falls_within = sum(e_true_within) / 10,
            e_point_est = mean(e_estimate),
            e_avg_width = mean(e_width),
            e_avg_lower = mean(e_lower),
            e_avg_upper = mean(e_upper)) %>% 
  mutate(bias = factor(bias,
                       levels = c("none", "a", "b", "c"),
                       labels = c("None", "Moderate", "High", "Very high")),
         type = factor(type,
                       levels = c("modified brm", "standard brm")))

ggplot(output_summary_ord) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = all_fall_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

ggplot(output_summary_ord) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = b_falls_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

ggplot(output_summary_ord) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = c_falls_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

ggplot(output_summary_ord) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = d_falls_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")

ggplot(output_summary_ord) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_col(aes(x = bias, y = e_falls_within, fill = type), position = position_dodge(.7), width = .7, alpha = .75) +
  labs(x = "Level of bias/unrepresentativeness of data",
       y = "Percent of simulations where true value\nfalls within 95% estimated range",
       fill = "Analysis:")