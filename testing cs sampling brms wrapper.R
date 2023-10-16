devtools::install_github("RyanHornby/csSampling")
library(csSampling)
library(rstan)
library(brms)
library(survey)
rstan_options(auto_write = TRUE)

load("~/Dropbox/R dbx/Survey Weighting Simulation/example data.RData")

svy_weight_des <- svydesign(ids = ~1, data = test_data, weights = test_data$weight)

# just looking to make a simple model estimating the weighted average
svy_brms <- cs_sampling_brms(svydes = svy_weight_des, 
                             brmsmod = brmsformula(score | weights(weight) ~ 1, center = FALSE), 
                             data = test_data,
                             family = gaussian())

mean(test_data$weight)
mean(test_data$weight) == 1
mean(test_data$weight) %% 1 # there is a tiny remainder that means it != 1

test_data$weight <- test_data$weight - (sum(test_data$weight) - length(test_data$weight))/length(test_data$weight)

mean(test_data$weights_2)
mean(test_data$weights_2) == 1
mean(test_data$weights_2) %% 1 # there is a tiny remainder that means it != 1

svy_brms <- cs_sampling_brms(svydes = svy_weight_des, 
                             brmsmod = brmsformula(score | weights(weights_2) ~ 1, center = FALSE), 
                             data = test_data,
                             family = gaussian())

# make sure the weights average to 1!
test_data$weights_check <- 1

svy_weight_des2 <- svydesign(ids = ~1, data = test_data, weights = test_data$weights_check)

svy_brms2 <- cs_sampling_brms(svydes = svy_weight_des2, 
                             brmsmod = brmsformula(score | weights(weights_check) ~ 1, center = FALSE), 
                             data = test_data,
                             family = gaussian())

# Throws a different error:
# Error in object@.MISC$stan_fit_instance$unconstrain_pars(pars) : 
# mismatch in number dimensions declared and found in context;
# processing stage=parameter initialization; variable name=b; dims declared=(1); dims found=()





cs_sampling <- function(svydes, mod_stan, par_stan = NA, data_stan,
                        ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                        rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
                        sampling_args = list()){
  require(rstan)
  require(survey)
  require(plyr)
  require(pkgcond)
  
  #Check weights
  #Check that the weights exist in both the survey object and the stan data
  #weights() returns full replicate weights set if svrepdesign
  if(rep_design){svyweights <- svydes$pweights}else{svyweights <-weights(svydes)}
  
  if (is.null(svyweights)) {
    if (!is.null(weights(data_stan))) {
      stop("No survey weights")
    }
  }
  if (is.null(weights(data_stan))) {
    if (!is.null(svyweights)) {
      warning("No stan data weights, using survey weights instead")
      data_stan$weights = weights(svydes)
    }
  }
  #Check that the weights are the same
  if (!isTRUE(all.equal(as.numeric(weights(data_stan)), as.numeric(svyweights)))) {
    stop("Survey weights and stan data weights do not match")
  }
  #Check that the mean is 1
  if (mean(weights(data_stan)) != 1) {
    stop("Mean of the weights is not 1")
  }
  
  
  print("stan fitting")
  out_stan  <- do.call(sampling, c(list(object = mod_stan, data = data_stan,
                                        pars = par_stan,
                                        chains = ctrl_stan$chains,
                                        iter = ctrl_stan$iter, warmup = ctrl_stan$warmup, thin = ctrl_stan$thin), sampling_args)
  )
  
  #Extract parameter draws and convert to unconstrained parameters
  
  #Get posterior mean (across all chains)
  par_samps_list <- rstan::extract(out_stan, permuted = TRUE)
  
  #If par_stan is not provided (NA) use all parameters (except "lp__", which is last)
  if(anyNA(par_stan)){
    par_stan <- names(par_samps_list)[-length(names(par_samps_list))]
  }
  
  #concatenate across multiple chains - save for later for export
  par_samps <- as.matrix(out_stan, pars = par_stan)
  
  #convert to list type input > convert to unconstrained parameterization > back to matrix/array
  for(i in 1:dim(par_samps)[1]){#just need the length here
    if(i == 1){upar_samps <- unconstrain_pars(out_stan, list_2D_row_subset(par_samps_list, i))
    }else{upar_samps <- rbind(upar_samps, unconstrain_pars(out_stan, list_2D_row_subset(par_samps_list, i)))}
  }
  row.names(upar_samps) <- 1:dim(par_samps)[1]
  
  upar_hat <- colMeans(upar_samps)
  
  #Estimate Hessian
  Hhat  <- -1*optimHess(upar_hat, gr = function(x){grad_log_prob(out_stan, x)})
  
  #create svrepdesign
  if(rep_design == TRUE){svyrep <- svydes
  }else{
    svyrep <- as.svrepdesign(design = svydes, type = ctrl_rep$type, replicates = ctrl_rep$replicates)
  }
  
  #Estimate Jhat = Var(gradient)
  print("gradient evaluation")
  rep_tmp <- withReplicates(design = svyrep, theta = grad_par, stanmod = mod_stan,
                            standata = data_stan, par_hat = upar_hat)#note upar_hat
  Jhat <- vcov(rep_tmp)
  
  #compute adjustment
  #use pivot for numerical stability - close to positive semi-definite if some parameters are highly correlated
  #(Q <- chol(m, pivot = TRUE))
  ## we can use this by
  #pivot <- attr(Q, "pivot")
  #Q[, order(pivot)]
  Hi <- solve(Hhat)
  V1 <- Hi%*%Jhat%*%Hi
  R1 <- chol(V1,pivot = TRUE)
  pivot <- attr(R1, "pivot")
  R1 <- R1[, order(pivot)]
  
  R2 <- chol(Hi, pivot = TRUE)
  pivot2 <- attr(R2, "pivot")
  R2 <- R2[, order(pivot2)]
  R2i <- solve(R2)
  R2iR1 <- R2i%*%R1
  
  #adjust samples
  upar_adj <- aaply(upar_samps, 1, DEadj, par_hat = upar_hat, R2R1 = R2iR1, .drop = TRUE)
  
  #back transform to constrained parameter space
  for(i in 1:dim(upar_adj)[1]){
    if(i == 1){par_adj <- unlist(constrain_pars(out_stan, upar_adj[i,])[par_stan])#drop derived quantities
    }else{par_adj <- rbind(par_adj, unlist(constrain_pars(out_stan, upar_adj[i,])[par_stan]))}
  }
  
  #make sure names are the same for sampled and adjusted parms
  row.names(par_adj) <- 1:dim(par_samps)[1]
  colnames(par_samps) <- colnames(par_adj)
  
  rtn = list(stan_fit = out_stan, sampled_parms = par_samps, adjusted_parms = par_adj)
  class(rtn) = c("cs_sampling", class(rtn))
  
  return(rtn)
  
}#end of cs_sampling

cs_sampling_brms <- function(svydes, brmsmod, data, family, par_brms = NA,prior = NULL, stanvars = NULL, knots = NULL,
                             ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                             rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
                             stancode_args = list(), standata_args = list(), sampling_args = list()) {
  
  
  stancode <- do.call(make_stancode, c(list(brmsmod, data = data, family = family, prior = prior, stanvars = stanvars, knots = knots), stancode_args))
  print("compiling stan model")
  mod_brms  <- stan_model(model_code = stancode)
  data_brms <- do.call(make_standata, c(list(brmsmod, data = data, family = family, prior = prior, stanvars = stanvars, knots = knots), standata_args))
  
  return(cs_sampling(svydes = svydes, mod_stan = mod_brms, par_stan = par_brms, data_stan = data_brms,
                     rep_design = rep_design, ctrl_rep = ctrl_rep, ctrl_stan = ctrl_stan, sampling_args))
  
}

plot.cs_sampling <- function(x, varnames = NULL) {
  
  datpl <- data.frame(rbind(as.matrix(x$sampled_parms), as.matrix(x$adjusted_parms))
                      , as.factor(c(rep("NO", dim(x$sampled_parms)[1]), rep("YES", dim(x$adjusted_parms)[1]))))
  names(datpl)[dim(x$sampled_parms)[2]+1] <- c("Adjust")
  rownames(datpl) <- NULL
  
  #subset to varnames
  if(!is.null(varnames)){datpl <- datpl[, c(varnames, "Adjust")]}
  
  require(GGally)
  
  my_ellipse <- function(data, mapping){
    ggplot(data = data, mapping = mapping) +
      geom_point()+
      stat_ellipse(level = 0.90, type = "norm", size = 2)
  }
  
  my_violin <- function(data, mapping){
    ggplot(data = data, mapping = mapping) +
      geom_violin(trim=TRUE,draw_quantiles = c(0.05, 0.5, 0.95),alpha=0.5, size = 1.5)
  }
  
  p1 <- ggpairs(datpl, mapping = aes(color = Adjust, alpha = 0.5), columns = c(1:(dim(datpl)[2]-1)),
                lower = list(continuous = my_ellipse))
  return(p1)
}

grad_par <- function(pwts, svydata, stanmod, standata,par_hat){
  #ignore svydata argument - it allows access to svy object data
  standata$weights <- pwts
  
  
  suppress_messages(out_stan  <- sampling(object = stanmod, data = standata,
                                          chains = 0, warmup = 0,), "the number of chains is less than 1")
  
  gradpar <- grad_log_prob(out_stan,par_hat)
  return(gradpar)
}#end of grad theta

DEadj <- function(par, par_hat, R2R1){
  par_adj <- (par - par_hat)%*%R2R1 + par_hat
  return(par_adj)
}

list_2D_row_subset <- function (nmlist, rindex)  {
  temp_list <- list()
  for (k in 1:length(nmlist)) {
    tmpdim <- dim(nmlist[[k]])
    ldim <- length(tmpdim)
    lcommas <- paste(rep(",", ldim - 1), collapse = " ")
    #copy over to new list - drop = FALSE retains ALL dimensions
    eval(parse(text = paste("temp_list$", names(nmlist)[k],
                            " <- ", "(nmlist$", names(nmlist)[k],
                            ")[rindex", lcommas, ",drop = FALSE]", sep = "")))
    #drop only first dimension of array - not others of size 1
    if(ldim > 1){
      eval(parse(text = paste("temp_list$", names(nmlist)[k],
                              " <- ", "array(temp_list$", names(nmlist)[k], ", dim = tmpdim[-1])", sep = "")))
    }
  }
  return(temp_list)
}

svy_brms <- cs_sampling_brms(svydes = svy_weight_des, 
                             brmsmod = brmsformula(score | weights(weights_2) ~ 1, center = FALSE), 
                             data = test_data,
                             family = gaussian())










##### here ####
list_2D_row_subset_DROP <- function (nmlist, rindex) 
{
  temp_list <- list()
  for (k in 1:length(nmlist)) {
    tmpdim <- dim(nmlist[[k]])
    ldim <- length(tmpdim)
    lcommas <- paste(rep(",", ldim - 1), collapse = " ")
    #copy over to new list - drop = FALSE retains ALL dimensions
    eval(parse(text = paste("temp_list$", names(nmlist)[k], 
                            " <- ", "(nmlist$", names(nmlist)[k], 
                            ")[rindex", lcommas, ",drop = FALSE]", sep = "")))
    #drop only first dimension of array - not others of size 1
    if(ldim > 1){
      eval(parse(text = paste("temp_list$", names(nmlist)[k], 
                              " <- ", "array(temp_list$", names(nmlist)[k], ", dim = tmpdim[-1])", sep = "")))
    }
    #if only had 1 dim which is the MCMC draw, make a double (no dim), rather than an array of dim 1 or 0
    if(ldim == 1){
      eval(parse(text = paste("temp_list$", names(nmlist)[k], 
                              " <- ", "as.double(temp_list$", names(nmlist)[k], ")", sep = "")))
    }	  	
  }
  return(temp_list)
}


cs_sampling_DROP <- function (svydes, mod_stan, par_stan = NA, data_stan, ctrl_stan = list(chains = 1, 
                                                                                           iter = 2000, warmup = 1000, thin = 1), rep_design = FALSE, 
                              ctrl_rep = list(replicates = 100, type = "mrbbootstrap"), 
                              sampling_args = list()) 
{
  require(rstan)
  require(survey)
  require(plyr)
  require(pkgcond)
  if (rep_design) {
    svyweights <- svydes$pweights
  }else {
    svyweights <- weights(svydes)
  }
  if (is.null(svyweights)) {
    if (!is.null(weights(data_stan))) {
      stop("No survey weights")
    }
  }
  if (is.null(weights(data_stan))) {
    if (!is.null(svyweights)) {
      warning("No stan data weights, using survey weights instead")
      data_stan$weights = weights(svydes)
    }
  }
  if (!isTRUE(all.equal(as.numeric(weights(data_stan)), as.numeric(svyweights)))) {
    stop("Survey weights and stan data weights do not match")
  }
  if (mean(weights(data_stan)) != 1) {
    stop("Mean of the weights is not 1")
  }
  print("stan fitting")
  out_stan <- do.call(sampling, c(list(object = mod_stan, data = data_stan, 
                                       pars = par_stan, chains = ctrl_stan$chains, iter = ctrl_stan$iter, 
                                       warmup = ctrl_stan$warmup, thin = ctrl_stan$thin), sampling_args))
  par_samps_list <- rstan::extract(out_stan, permuted = TRUE)
  if (anyNA(par_stan)) {
    par_stan <- names(par_samps_list)[-length(names(par_samps_list))]
  }
  par_samps <- as.matrix(out_stan, pars = par_stan)
  for (i in 1:dim(par_samps)[1]) {
    if (i == 1) {
      upar_samps <- unconstrain_pars(out_stan, list_2D_row_subset_DROP(par_samps_list, 
                                                                       i))
    }
    else {
      upar_samps <- rbind(upar_samps, unconstrain_pars(out_stan, 
                                                       list_2D_row_subset_DROP(par_samps_list, i)))
    }
  }
  row.names(upar_samps) <- 1:dim(par_samps)[1]
  upar_hat <- colMeans(upar_samps)
  Hhat <- -1 * optimHess(upar_hat, gr = function(x) {
    grad_log_prob(out_stan, x)
  })
  if (rep_design == TRUE) {
    svyrep <- svydes
  }
  else {
    svyrep <- as.svrepdesign(design = svydes, type = ctrl_rep$type, 
                             replicates = ctrl_rep$replicates)
  }
  print("gradient evaluation")
  rep_tmp <- withReplicates(design = svyrep, theta = grad_par, 
                            stanmod = mod_stan, standata = data_stan, par_hat = upar_hat)
  Jhat <- vcov(rep_tmp)
  Hi <- solve(Hhat)
  V1 <- Hi %*% Jhat %*% Hi
  R1 <- chol(V1, pivot = TRUE)
  pivot <- attr(R1, "pivot")
  R1 <- R1[, order(pivot)]
  R2 <- chol(Hi, pivot = TRUE)
  pivot2 <- attr(R2, "pivot")
  R2 <- R2[, order(pivot2)]
  R2i <- solve(R2)
  R2iR1 <- R2i %*% R1
  upar_adj <- aaply(upar_samps, 1, DEadj, par_hat = upar_hat, 
                    R2R1 = R2iR1, .drop = TRUE)
  for (i in 1:dim(upar_adj)[1]) {
    if (i == 1) {
      par_adj <- unlist(constrain_pars(out_stan, upar_adj[i, 
      ])[par_stan])
    }
    else {
      par_adj <- rbind(par_adj, unlist(constrain_pars(out_stan, 
                                                      upar_adj[i, ])[par_stan]))
    }
  }
  row.names(par_adj) <- 1:dim(par_samps)[1]
  colnames(par_samps) <- colnames(par_adj)
  rtn = list(stan_fit = out_stan, sampled_parms = par_samps, 
             adjusted_parms = par_adj)
  class(rtn) = c("cs_sampling", class(rtn))
  return(rtn)
}

cs_sampling_brms_DROP <- function (svydes, brmsmod, data, family, par_brms = NA, prior = NULL, 
                                   stanvars = NULL, knots = NULL, ctrl_stan = list(chains = 1, 
                                                                                   iter = 2000, warmup = 1000, thin = 1), rep_design = FALSE, 
                                   ctrl_rep = list(replicates = 100, type = "mrbbootstrap"), 
                                   stancode_args = list(), standata_args = list(), sampling_args = list()) 
{
  stancode <- do.call(make_stancode, c(list(brmsmod, data = data, 
                                            family = family, prior = prior, stanvars = stanvars, 
                                            knots = knots), stancode_args))
  print("compiling stan model")
  mod_brms <- stan_model(model_code = stancode)
  data_brms <- do.call(make_standata, c(list(brmsmod, data = data, 
                                             family = family, prior = prior, stanvars = stanvars, 
                                             knots = knots), standata_args))
  return(cs_sampling_DROP(svydes = svydes, mod_stan = mod_brms, 
                          par_stan = par_brms, data_stan = data_brms, rep_design = rep_design, 
                          ctrl_rep = ctrl_rep, ctrl_stan = ctrl_stan, sampling_args))
}

svy_brms <- cs_sampling_brms_DROP(svydes = svy_weight_des, 
                             brmsmod = brmsformula(score | weights(weight) ~ 1, center = FALSE), 
                             data = test_data,
                             family = gaussian())

summary(svy_brms)
dim(svy_brms)
names(svy_brms$svy_brms)
svy_brms$sampled_parms[1:5, 1:3]
plot(svy_brms)

test_data

quantile(svy_brms$sampled_parms[ , 1], .025)
quantile(svy_brms$sampled_parms[ , 1], .975)

brm_version <-
brm(formula = score | weights(weight_mod) ~ 1,
    data = test_data)
