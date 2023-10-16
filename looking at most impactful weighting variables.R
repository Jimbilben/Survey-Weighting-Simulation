library(tidystats)

generaldesigneffect(data_nofake_weighted$weight)


race <- tibble(race = c("White or Caucasian", "Hispanic or Latino", "Other",
                           "Black or African American", "Asian or Asian American"),
              pct = c(.578, .187, .055, .121, .059) * 100,
              type = "target")
gender <- tibble(gender = c('Female', 'Male'),
                pct = c(0.513, 0.487) * 100,
                type = "target")
age <- tibble(age = c("18-24", "25-44", "45-64", "65+"),
             pct = c(.1304, .3505, .3478, .1713) * 100,
             type = "target")
education <- tibble(education = c('Completed graduate school', 'Graduated from college',
                                'Some college, no degree', 'Graduated from high school',
                                'Less than high school'),
                   pct = c(0.11376, 0.19265, 0.30353, 0.27579, 0.11427) * 100,
                   type = "target")
income.ces <- tibble(income.ces = c('Under $20,000',
                                 'Between $20,000 and $49,999',
                                 'Between $50,000 and $79,999',
                                 'Between $80,000 and $99,999', 
                                 'Between $100,000 and $150,000',
                                 'Over $150,000'),
                    pct = c(.093,
                            0.21,
                            0.205,
                            0.111,
                            0.185,
                            0.196) * 100,
                    type = "target")
region <- tibble(region = c("Midwest", "Mountains", "Northeast",
                             "Pacific", "South", "Southwest", "Southeast"),
                pct = c(.20949, .052430, .15822, .16080, .21637, .055450, .14725) * 100,
                type = "target")
party <- tibble(party = c("Republican", "Independent", "Democrat"),
               pct = c(.276, .436, .288) * 100,
               type = "target")

libcon <- tibble(libcon = c("Liberal", "Moderate", "Conservative"),
                pct = c(.333, .347, .320) * 100,
                type = "target")

race_data <-
data_nofake_weighted %>%
  count_data(race) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(race) %>% 
  mutate(y = race,
         variable = "race")

age_data <-
  data_nofake_weighted %>%
  count_data(age) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(age) %>% 
  mutate(y = age,
         variable = "age")

education_data <-
  data_nofake_weighted %>%
  count_data(education) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(education) %>% 
  mutate(y = education,
         variable = "education")

income.ces_data <-
  data_nofake_weighted %>%
  count_data(income.ces) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(income.ces) %>% 
  mutate(y = income.ces,
         variable = "income")

party_data <-
  data_nofake_weighted %>%
  count_data(party) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(party) %>% 
  mutate(y = party,
         variable = "party")

libcon_data <-
  data_nofake_weighted %>%
  count_data(libcon) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(libcon) %>% 
  mutate(y = libcon,
         variable = "libcon")

region_data <-
  data_nofake_weighted %>%
  count_data(region) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(region) %>% 
  mutate(y = region,
         variable = "region")

gender_data <-
  data_nofake_weighted %>%
  count_data(gender) %>% 
  mutate(type = "observed") %>% 
  select(-n) %>% 
  bind_rows(gender) %>% 
  mutate(y = gender,
         variable = "gender")

all_data <- bind_rows(race_data[ , -1],
                      age_data[ , -1],
                      income.ces_data[ , -1],
                      education_data[ , -1],
                      gender_data[ , -1],
                      region_data[ , -1],
                      libcon_data[ , -1],
                      party_data[ , -1])

ggplot(all_data, aes(x = pct, y = y, fill = type)) +
  geom_col(position = position_dodge(.8), width = .8) +
  geom_text(aes(label = jimbilben::nice_num(pct, 1, FALSE)), size = 2, position = position_dodge(.8), hjust = -.15) +
  facet_wrap(~variable, scales = "free_y") +
  theme(
    axis.title = element_blank()
  )

weight_formula <- log_weight ~ 0 + age + education + income.ces + race + party + libcon + gender + region

ggplot(data_nofake_weighted) +
  geom_density(aes(x = weight))
ggplot(data_nofake_weighted) +
  geom_density(aes(x = log(weight)))
ggplot(data_nofake_weighted) +
  geom_density(aes(x = weight_change))


data_nofake_weighted$log_weight <- log(data_nofake_weighted$weight)

1/5
1 / (1/5)
data_nofake_weighted <- data_nofake_weighted %>% 
  mutate(weight_change = log(1 / (1 / weight)))
max(data_nofake_weighted$weight)
summary(lm(weight_formula,
   data_nofake_weighted))

data_nofake_weighted %>% 
  group_by(race) %>% 
  summarise(mean_weight = mean(weight))

weight_mapper <- function(variable, var_string, data = data_nofake_weighted) {
  
  output <- data %>% 
    group_by({{variable}}) %>% 
    summarise(mean_weight = mean(weight)) %>% 
    mutate(y = {{variable}},
           variable = var_string)
  
  output <- output[ , -1]
  
  return(output)
  
}

weights_df <- map2_df(.x = vars(age, income.ces, education, race, gender, region, party, libcon),
                   .y = c("age", "income.ces", "education", "race", "gender", "region", "party", "libcon"),
                   .f = weight_mapper)

ggplot(weights_df, aes(x = mean_weight, y = y)) +
  geom_col(fill = "#619ecd", alpha = .6) +
  geom_vline(aes(xintercept = 1), linetype = "dashed", alpha = .5) +
  geom_text(aes(label = jimbilben::nice_num(mean_weight, 2, FALSE)), size = 3.5, hjust = -.1) +
  facet_wrap(~variable, scales = "free_y") +
  theme(
    axis.title.y = element_blank()
  )

library(survey)
generaldesigneffect(data_weighting$total_weight)

library(jimbilben)

brm_code()

data_nofake_weighted$weight_4 <-  data_nofake_weighted$weight / 4
data_nofake_weighted$weight_2 <-  data_nofake_weighted$weight / 2


data_nofake_weighted$weight_43k <-  data_nofake_weighted$weight / 8
data_nofake_weighted$weight_2 <-  data_nofake_weighted$weight / 2

names(data_nofake_weighted)
count_data(data_nofake_weighted,
           terms_ea)
formula_4weight <-  terms_ea | weights(weight_4) ~ 1
formula_2weight <-  terms_ea | weights(weight_2) ~ 1
  
formula_4weight3k <-  terms_ea | weights(weight_43k) ~ 1
formula_2weight3k <-  terms_ea | weights(weight_23k) ~ 1

weight_4_6k <- 
  brm(formula = formula_4weight,
      family = bernoulli(),
      data = data_nofake_weighted,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      #prior = my_priors,
      chains = 4,
      cores = 4,
      iter = 1250,
      warmup = 250,
      #init = 0,
      backend = 'cmdstanr',
      threads = threading(4),
      seed = 1010,
      stan_model_args=list(stanc_options = list('O1')))

weight_2_6k <- 
  brm(formula = formula_2weight,
      family = bernoulli(),
      data = data_nofake_weighted,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      #prior = my_priors,
      chains = 4,
      cores = 4,
      iter = 1250,
      warmup = 250,
      #init = 0,
      backend = 'cmdstanr',
      threads = threading(4),
      seed = 1010,
      stan_model_args=list(stanc_options = list('O1')))
  


weight_4_3k <- 
  brm(formula = formula_4weight3k,
      family = bernoulli(),
      data = data_nofake_weighted,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      #prior = my_priors,
      chains = 4,
      cores = 4,
      iter = 1250,
      warmup = 250,
      #init = 0,
      backend = 'cmdstanr',
      threads = threading(4),
      seed = 1010,
      stan_model_args=list(stanc_options = list('O1')))

weight_2_3k <- 
  brm(formula = formula_2weight,
      family = bernoulli(),
      data = data_nofake_weighted,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      #prior = my_priors,
      chains = 4,
      cores = 4,
      iter = 1250,
      warmup = 250,
      #init = 0,
      backend = 'cmdstanr',
      threads = threading(4),
      seed = 1010,
      stan_model_args=list(stanc_options = list('O1')))


# 6k weighted at /4
# 6k weighted at /2
# 3k weighted at /4
# 3k weights at /2
# 1.5k weighted at /4
# 1.5k weights at /2

weighter <- function(weight_div, name_weighting) {
  
  data_nofake_weighted$weight <-  data_nofake_weighted$weight / weight_div
  
  formula_weight <-  terms_ea | weights(weight) ~ 1
  
  model <- brm(formula = formula_weight,
               family = bernoulli(),
               data = data_nofake_weighted,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               #prior = my_priors,
               chains = 4,
               cores = 4,
               iter = 1250,
               warmup = 250,
               #init = 0,
               backend = 'cmdstanr',
               threads = threading(4),
               seed = 1010,
               stan_model_args=list(stanc_options = list('O1')))
  
  output <-
    as_tibble(posterior_epred(model,
                              newdata = data_nofake_weighted[1, ])) %>% 
    jimbilben::nice_post(V1,
                         percentage = TRUE) %>% 
    mutate(name_weight = name_weighting)
  
  return(output)
  
}

weight_divs <- c(4, 2, 8, 4, 16, 8, 24, 12)
weight_names <- c("6k, Deff 4", "6k, Deff 2", "3k, Deff 4", "3k, Deff 2", "1.5k, Deff 4", "1.5k, Deff 2", "1k, Deff 4", "1k, Deff 2")

widths <- 
map2_dfr(.x = weight_divs,
         .y = weight_names,
         .f = weighter)


weighter_2 <- function(weight_div, name_weighting) {
  
  data_nofake_weighted$weight <-  data_nofake_weighted$weight / weight_div
  
  formula_weight <-  terms_ebm | weights(weight) ~ 1
  
  model <- brm(formula = formula_weight,
               family = bernoulli(),
               data = data_nofake_weighted,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               #prior = my_priors,
               chains = 4,
               cores = 4,
               iter = 1250,
               warmup = 250,
               #init = 0,
               backend = 'cmdstanr',
               threads = threading(4),
               seed = 1010,
               stan_model_args=list(stanc_options = list('O1')))
  
  output <-
    as_tibble(posterior_epred(model,
                              newdata = data_nofake_weighted[1, ])) %>% 
    jimbilben::nice_post(V1,
                         percentage = TRUE) %>% 
    mutate(name_weight = name_weighting)
  
  return(output)
  
}

widths_ebm <- 
  map2_dfr(.x = weight_divs,
           .y = weight_names,
           .f = weighter_2)

names(data_nofake_weighted)
unique(data_nofake_weighted$urban_rural)


weighter_3 <- function(weight_div, name_weighting) {
  
  data_nofake_weighted$weight <-  data_nofake_weighted$weight / weight_div
  
  formula_weight <-  urban_rural | weights(weight) ~ 1
  
  model <- brm(formula = formula_weight,
               family = categorical(),
               data = data_nofake_weighted,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               #prior = my_priors,
               chains = 4,
               cores = 4,
               iter = 1250,
               warmup = 250,
               #init = 0,
               backend = 'cmdstanr',
               threads = threading(4),
               seed = 1010,
               stan_model_args=list(stanc_options = list('O1')))
  
  output <-
    as_tibble(posterior_epred(model,
                              newdata = data_nofake_weighted[1, ])) %>%
    pivot_longer(cols = everything(),
                 names_to = "outcome",
                 values_to = "perc") %>% 
    group_by(outcome) %>% 
    jimbilben::nice_post(perc,
                         percentage = TRUE) %>% 
    mutate(name_weight = name_weighting)
  
  return(output)
  
}

weighter_3(1, "test")

widths_urb <- 
  map2_dfr(.x = weight_divs,
           .y = weight_names,
           .f = weighter_3)

print(widths_urb,
      n = 24)

widths <- 
  map2_dfr(.x = weight_divs,
           .y = weight_names,
           .f = weighter)

widths_ebm <- 
  map2_dfr(.x = weight_divs,
           .y = weight_names,
           .f = weighter_2)



ggplot(widths) +
  geom_point(aes(x = hdi_width, y = name_weight))

widths_ebm$name_weight <- factor(widths_ebm$name_weight,
                                 levels = c("6k, Deff 4", "3k, Deff 4", "1.5k, Deff 4","1k, Deff 4", "6k, Deff 2", "3k, Deff 2", "1.5k, Deff 2",  "1k, Deff 2"))
ggplot(widths_ebm) +
  geom_point(aes(x = hdi_width, y = name_weight)) +
  geom_point(data = widths, aes(x = hdi_width, y = name_weight), color = "red")

widths_urb$name_weight <- factor(widths_urb$name_weight,
                                 levels = c("6k, Deff 4", "3k, Deff 4", "1.5k, Deff 4","1k, Deff 4", "6k, Deff 2", "3k, Deff 2", "1.5k, Deff 2",  "1k, Deff 2"))

ggplot(widths_urb) +
  geom_point(aes(x = hdi_width, y = name_weight, color = outcome))
