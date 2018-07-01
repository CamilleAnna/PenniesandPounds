# Pound&Pennies review - Simulations

# LAST EDIT: 2018/06/19: changing parameter sampling to get more stuff in the top right corner: increasing range of f sampling to 4 instead of 2
# 

# EDIT: 2018/06/12: changed sampling of parameter a to  sample a bit more in the low proportion of F usage and high inapropriate usage for high F usage
# decided to use a gamma distribution, rgamma(1, 1, 0.3)

# EDIT: 2018/06/08: changed sampling of parameter a to a<-exp(runif(1,-2.5,2.5)) to sample more in the f usage > 0.5

# EDIT: 2018/06/06
# putting run time to 600 generations to make sure we see end equilibrium in new drug model 
# (was checking that all three strains end up at the same equilibrium)
# changed initial conditions before burn in phase:
# all resistant at zero (they will emerge from conversion rates)
# and sensitive to 0.1 each so that 0.2 of the ill compartment is filled

# EDIT: 2018/06/05
# edit by luke: change in bb and e parameters
# re implement middle line drug intervention for the new drug model
# modifying measured recorded: remove morbidity and number of untreatable infections, add average length of infection and proba od death (+ number of transmission, required tp compute the two previous ones)
# required adding initial conditions for dTrans: zero before burn in phase, and then number of infections before the run
# modified initial conditions before burn in to have it to match equilibrium conditions in the absence of resistance
# so that we see evolution of resistance hapenning during the burn in phase --- > ACTUALLY STILL HAVE TO FIGURE IT OUT, PUT TEMPORARY VALUES FOR NOW
# EDIT: 2018/06/04
# making burn-in phase longer
# add intervention in the carriage class for both adjuvant and new drug models
# value to compute is L - F / N
# re make simulation part longer but thin intervals
# EDIT: 2018/05/14
# re-adding no intervention model
# reduce time of simulation to 150 because equilibrium reached much before 500 generations in most cases
# EDIT: 2018/04/20
# Adding another measure of front line drug proportion usgae, within ill class only and not within both carrier and ill class
# EDIT: 2018/04/19
# Adding new measure of morbidity
# EDIT: 2018/04/18
# increasing alpha range to 15 and f range from 0 to 2 (f now is not a proportion but a variable size of compartment)
# EDIT: 2018/04/15
# - changing ways to calculate F/L usage and inapropriate usage: take proportions, and take in during burn in phase (last time point of burn in phase = initial conditions)
# - remove Cxxx,M class
# - change how initial conditions are set for the new_drug model accordingly
# - hashtag all the swap models thing cause we will not run it anymore
# - hashtag all the introduction as middle line as well


#!/usr/bin/env Rscript

library("optparse")


option_list = list(
  make_option(c("-s", "--simulations"), type="integer", default=NULL,
              help="number of simulations to run", metavar="integer"),
  make_option(c("-d", "--seed"), type="integer", default=NULL,
              help="seed for numbers generation", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$simulations)){
  print_help(opt_parser)
  stop("provide a number of simulations to run", call.=FALSE)
}
if (is.null(opt$seed)){
  print_help(opt_parser)
  stop("provide a number for the seed. This might be the array number if this job is launched as an array job.", call.=FALSE)
}




# run simulations with new model version (with mis-usage)


library(ggplot2)
library(deSolve)
library(gridExtra)


#setwd('/Volumes/csce/biology/users/s1687811/phD/PoundsAndPennies/PPmodel_review_June/')
#source('Model_functions.R')



# Draw parameters combinations ----
# draw a number of parameter combinations, stored in a dataframe
# then run the models going through this table of parameters values
# each draw contains all parameters, but not all of them is necessarily used by all models, only the releavant ones
# so all models are run on same sample of parameters +/- the ones specific to each model
# we do not do any check on the combination drawn except that that no parameter is equal to zero
# We necessarily have bb < bw, because we set bb = 0.5 * bw 

nb_sim<- opt$simulations
my_seed<- opt$seed


# We draw a parameter combination and run a burn in phase of 100 generation
# the state at the end of this burn in phase is the initial state for the simulation with that parameter combination
# some initial states are irrelevant though: absence of resistance or absence of disease at all
# hence, simulations go as follow
# we draw a parameter combination
# check that no parameter equal zero
# run the burn in phase (model of no intervention, 100 generations, time steop 0.05)
# record: the out of burn in run + the parameter combination used + the last line or out which will be our vector of initial conditions for simuls
# check those initial conditions: want presence of disease, and presence of resistance
# if conditions not met, re-draw new parameter, re-run, replace in the records
# if conditions are met, draw the next parameter combination
# do this until we reach nb_simuls


# In order to be able to re-get the exact same parameter combinations in future, we need to know the seed but also the number of the draw.
# Hence we also record the seed and the draw number in the parameters table



# initialise pars ----
pars<- data.frame(bw = numeric(), bb = numeric(), c = numeric(), d = numeric(), g = numeric(),
                  z = numeric(), a = numeric(), e = numeric(),  p = numeric(), mu = numeric(),
                  f = numeric(), my_seed = numeric(), draw = numeric())

# initialise some vectors
initial_states<- vector('list')                  # to store initial states
outs_burn_in<- vector('list', length=nrow(pars)) # to store outputs of burn in phase
prop_f_usage<- NULL                              # Measure to record, proportion of F usage (total)
prop_f_usage_ill<- NULL                          # Measure to record, proportion of F usage (in ill class)
prop_inapropriate_usage<- NULL                   # Measure to record, proportion of inapropriate usage (by standard selection)



# Burin-in run parameters
burn_in_time <- seq(0, 300, by = 0.05)

state_burn_in = c(C_00_u = 0.1, C_00_f = 0.1, I_00_f = 0.1, I_00_l = 0.1, 
                  C_10_u = 0, C_10_f = 0, I_10_f = 0, I_10_l = 0,  
                  C_01_u = 0, C_01_f = 0, I_01_f = 0, I_01_l = 0, 
                  C_11_u = 0, C_11_f = 0, I_11_f = 0, I_11_l = 0,
                  D = 0, dTrans = 0)

index = 0
i=1
set.seed(my_seed)

while(nrow(pars) < nb_sim){
  
  # Re-set GOON as FALSE at each parameter draw
  GOON=FALSE
  
  # while GOON condition is FALSE, re-draw parameters
  while(GOON == FALSE){
    
    #while(index<40){
    index = index + 1 # To record how many draw we do to get x parameter combinations that lead to conditions we want
    
    print(paste('burn in', i, index, sep=' : '))
    
    # draw pars + check not equal to zero. If zero, re-draw
    cond = TRUE
    while(cond == TRUE){
      bw<-runif(1,2,5)
      bb<-runif(1,0,1) * bw # edit by Luke
      c<-runif(1,0,1)
      d<-runif(1,0,2)
      g<-runif(1,0,2)
      z<-runif(1,0,1)
      a<-rgamma(1, 1, 0.3) 
      e<-runif(1,0,1) * g # edit by Luke
      p<-runif(1,0,1)
      mu<-runif(1,0,0.1)#0.01
      f<-runif(1,0,4)
      
      pars_burn<- c(bw=bw, bb=bb, c=c, d=d, g=g, z=z, a=a, e=e, p=p, mu=mu, f=f, my_seed = my_seed, draw = index)
      
      cond = (sum(pars_burn == 0) > 0)
      if(cond == TRUE){print("sampling new parameters")}
    }
    
    #}
    
    # if not zero, run model for burn in phase
    # and record pars, out, initial state
    
    
    
    outs_burn_in[[i]]<- ode(y = state_burn_in, times = burn_in_time, func = mod1, parms = pars_burn)     # run model
    outs_burn_in[[i]]<- outs_burn_in[[i]][outs_burn_in[[i]][,'time'] %in% seq(0, tail(burn_in_time,1)),] # record only 300 points
    
    pars[i,]<- pars_burn                                  # record parameter combination
    initial_states[[i]]<- tail(outs_burn_in[[i]], 1)[,-1] # record initial states
    
    # modify dTrans initial state before the run --> number of infections at the end of burn in phase
    initial_states[[i]]['dTrans'] <- sum(initial_states[[i]][c('I_00_f', 'I_00_l', 'I_10_f', 'I_10_l', 'I_01_f', 'I_01_l', 'I_11_f', 'I_11_l')])
    
    # check outcome of burn in: We don't want extinction and we want presence of resistance
    # If presence of resistance, means it is not extinct
    # So simply check for presence of resistance
    # i.e. sum of all other strains than 00 must be > 0
    # If this conditions is NOT met, GOON is FALSE
    # we go back at beginning of this loop and re-draw a new parameter combination
    # we replace the pars, outs and initial_states in the record
    # the 'index' values increases by one, so that we can keep track of how many combination we draw in total
    # (and record it to know what iteration of the seed produced a given outcome, for reproducibility)
    
    GOON<- round(sum(initial_states[[i]][c('C_10_u', 'C_10_f','I_10_f', 'I_10_l', 'C_01_u', 'C_01_f','I_01_f', 'I_01_l', 'C_11_u','C_11_f','I_11_f', 'I_11_l')]), 3) > 0 # (the Deads and dTrans are not in there obviously...)
    
    
  }
  
  # if GOON is TRUE, we exit of this while loop and i increases by one
  # we draw the next parameter combination
  
  # Relative usage of front and last drugs: 
  # total front line usage (both in Carriers and Ill class)/ total usage of drug (in both carriers and Ill classes)
  # F/F+L
  prop_f_usage<- c(prop_f_usage, 
                   sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f','I_00_f', 'I_10_f', 'I_01_f', 'I_11_f')])/sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f','I_00_f', 'I_10_f', 'I_01_f', 'I_11_f','I_00_l', 'I_10_l', 'I_01_l', 'I_11_l')]))
  
  
  prop_f_usage_ill<- c(prop_f_usage_ill, sum(initial_states[[i]][c('I_00_f', 'I_10_f', 'I_01_f', 'I_11_f')])/sum(initial_states[[i]][c('I_00_f', 'I_10_f', 'I_01_f', 'I_11_f','I_00_l', 'I_10_l', 'I_01_l', 'I_11_l')]))
  
  # Proportion of inapropriate usage:
  # front line drug used in carriers / total use of front line drug
  # Cxx,F / (Cxx,F + Ixx,F)
  
  
  prop_inapropriate_usage<- c(prop_inapropriate_usage,
                              sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f')])/sum(initial_states[[i]][c('C_00_f', 'C_10_f', 'C_01_f', 'C_11_f','I_00_f', 'I_10_f', 'I_01_f', 'I_11_f')]))
  
  
  
  i = i + 1
  
  
}

# ----


# Now that we have all the initial conditions, we must adjust them for each model

# DIAGNOSTIC MODEL
# we have no ill individuals resistant to F being front line treated (I_10_f) or multi resistant being front line treated (I_11_f).
# But in burn-in phase (before deployment of the diagnostic) they do
# so it IS releavant to keep those classe in the burn in phase
# but at the time of the diagnostic starting, all I_10_f and I_10_f must be converted. All I_10_f  are transferred to I_10_l, and I_11_f are transferred to I_11_l
# So we modify the initial conditions to have I_10_l = I_10_l + I_10_f , I_11_l = I_11_l + I_11_f 
# then, for code convenience, I keep the I_10_f  and I_11_f  classes but set their rate of change to zero (, I_10_f  = 0, I_11_f  = 0), so that they don't matter
# Similarly, for last line diagnostic, I set I_01_f = I_01_f + I_01_l, I_11_f = I_11_f + I_11_l and then I_01_l = 0, I_11_l = 0



# NEW DRUG MODEL (ADDING IT)
# For the new drug model, there is only one function, describing the dynamic with three drugs F, M, L
# The only difference to model is if we introduce as front, middle or last is the change in initial conditions
# - For introduction as front:
# all resistance to previously F drug become resistant to M
# resistance to F is now zero because it is a new drug
# L remain what there were
# - For introduction as last:
# resistance to F remain F
# All resistance to L now are resistant to M
# resistance to L is zero because it is a new drug

# - For introduction as middle:
# it actually does not change anything, those that were resistant to former front line druf remain resistant to F because it is still the same
# and same of those resistant for last line drug

# In addition, we transition from a system of 4*4+1 equation to 5*8+1 equations because of the introduction of the new drug
# So the table of initial conditions must be extended

# NO INTERVENTION MODEL
# No need of adjusting the initial conditions to run this model

# ADJUVANT MODEL
# No need of adjusting the initial conditions to run this model


# Copying initial states into new vectors that will be modified according to each model needs

# DIAGNOSTIC MODEL
initial_states_front_diag<- initial_states
initial_states_last_diag<- initial_states

# NEW DRUG MODEL
initial_states_front_new<- list()
initial_states_last_new<- list()
initial_states_middle_new<- list()





# Now modify those initial states
for(i in 1:length(initial_states))
{
  
  # FOR DIAGNOSTIC MODEL
  # I_10_l = I_10_l + I_10_f , I_11_l = I_11_l + I_11_f , I_10_f  = 0, I_11_f  = 0
  initial_states_front_diag[[i]][c('I_10_l', 'I_11_l')] = 
    initial_states_front_diag[[i]][c('I_10_l', 'I_11_l')] + 
    initial_states_front_diag[[i]][c('I_10_f', 'I_11_f')] 
  initial_states_front_diag[[i]][c('I_10_f', 'I_11_f')] = 0
  
  # LAST LINE DIAGNOSTIC
  # I_01_f = I_01_f + I_01_l, I_11_f = I_11_f + I_11_l, I_01_l = 0, I_11_l = 0
  initial_states_last_diag[[i]][c('I_01_f', 'I_11_f')] = 
    initial_states_last_diag[[i]][c('I_01_f', 'I_11_f')] + 
    initial_states_last_diag[[i]][c('I_01_l', 'I_11_l')] 
  initial_states_last_diag[[i]][c('I_01_l', 'I_11_l')] = 0
  
  
  
  # FOR NEW DRUG MODEL
  # Treatments change with the intervention, those that were treated with front line drug, switch to the new front line drug when this one is added to the market
  # and from now on, follow normal procedure with escalation to middle drug (fomer front line), and then last line (which has not changed)
  # For the resistance, when introducing as front line, 1xx = 0 , 010 = 10, 011 = 11, 001 = 01
  # basically we add a zero in front of resistance profiles ...
  # and all Xxxx,m = 0, because everyone either switch to the new front line drug, or remain with the former last line one
  
  initial_states_front_new[[i]]<-     c(C_000_u = initial_states[[i]][['C_00_u']],
                                        C_000_f = initial_states[[i]][['C_00_f']],
                                        
                                        I_000_f = initial_states[[i]][['I_00_f']],
                                        I_000_m = 0,
                                        I_000_l = initial_states[[i]][['I_00_l']],
                                        
                                        C_100_u = 0,
                                        C_100_f = 0,
                                        
                                        I_100_f = 0,
                                        I_100_m = 0,
                                        I_100_l = 0,
                                        
                                        C_010_u = initial_states[[i]][['C_10_u']],
                                        C_010_f = initial_states[[i]][['C_10_f']],
                                        
                                        I_010_f = initial_states[[i]][['I_10_f']],
                                        I_010_m = 0,
                                        I_010_l = initial_states[[i]][['I_10_l']],
                                        
                                        C_001_u = initial_states[[i]][['C_01_u']],
                                        C_001_f = initial_states[[i]][['C_01_f']],
                                        
                                        I_001_f = initial_states[[i]][['I_01_f']],
                                        I_001_m = 0,
                                        I_001_l = initial_states[[i]][['I_01_l']],
                                        
                                        C_110_u = 0,
                                        C_110_f = 0,
                                        
                                        I_110_f = 0,
                                        I_110_m = 0,
                                        I_110_l = 0,
                                        
                                        C_101_u = 0,
                                        C_101_f = 0,
                                        
                                        I_101_f = 0,
                                        I_101_m = 0,
                                        I_101_l = 0,
                                        
                                        C_011_u = initial_states[[i]][['C_11_u']],
                                        C_011_f = initial_states[[i]][['C_11_f']],
                                        
                                        I_011_f = initial_states[[i]][['I_11_f']],
                                        I_011_m = 0,
                                        I_011_l = initial_states[[i]][['I_11_l']],
                                        
                                        C_111_u = 0,
                                        C_111_f = 0,
                                        
                                        I_111_f = 0,
                                        I_111_m = 0,
                                        I_111_l = 0,
                                        
                                        D = initial_states[[i]][['D']],
                                        dTrans = initial_states[[i]][['dTrans']])
  
  
  
  initial_states_middle_new[[i]]<-    c(C_000_u = initial_states[[i]][['C_00_u']],
                                        C_000_f = initial_states[[i]][['C_00_f']],
                                        
                                        I_000_f = initial_states[[i]][['I_00_f']],
                                        I_000_m = 0,
                                        I_000_l = initial_states[[i]][['I_00_l']],
                                        
                                        C_100_u = initial_states[[i]][['C_10_u']],
                                        C_100_f = initial_states[[i]][['C_10_f']],
                                        
                                        I_100_f = initial_states[[i]][['I_10_f']],
                                        I_100_m = 0,
                                        I_100_l = initial_states[[i]][['I_10_l']],
                                        
                                        
                                        
                                        C_010_u = 0,
                                        C_010_f = 0,
                                        
                                        I_010_f = 0,
                                        I_010_m = 0,
                                        I_010_l = 0,
                                        
                                        C_001_u = initial_states[[i]][['C_01_u']],
                                        C_001_f = initial_states[[i]][['C_01_f']],
                                        
                                        I_001_f = initial_states[[i]][['I_01_f']],
                                        I_001_m = 0,
                                        I_001_l = initial_states[[i]][['I_01_l']],
                                        
                                        C_110_u = 0,
                                        C_110_f = 0,
                                        
                                        I_110_f = 0,
                                        I_110_m = 0,
                                        I_110_l = 0,
                                        
                                        C_101_u = initial_states[[i]][['C_11_u']],
                                        C_101_f = initial_states[[i]][['C_11_f']],
                                        
                                        I_101_f = initial_states[[i]][['I_11_f']],
                                        I_101_m = 0,
                                        I_101_l = initial_states[[i]][['I_11_l']],
                                        
                                        C_011_u = 0,
                                        C_011_f = 0,
                                        
                                        I_011_f = 0,
                                        I_011_m = 0,
                                        I_011_l = 0,
                                        
                                        C_111_u = 0,
                                        C_111_f = 0,
                                        
                                        I_111_f = 0,
                                        I_111_m = 0,
                                        I_111_l = 0,
                                        
                                        D = initial_states[[i]][['D']],
                                        dTrans = initial_states[[i]][['dTrans']])
  
  
  
  # When introducing drug as last line
  # All treated with former last line swithc to the new last line
  # so no one being treated with middle drug (former last line one)
  # But keep this former last line as an intermediate before escalating to the new last line one
  # 00 = 000, 10 = 100, 01 = 010, 11 = 110 (add a zero to the right of all resistance profiles)
  # 001, 011, 101 and 111 = zero because no one is resistant to this new last line drug
  
  initial_states_last_new[[i]]<-   c(C_000_u = initial_states[[i]][['C_00_u']],
                                     C_000_f = initial_states[[i]][['C_00_f']],
                                     
                                     I_000_f = initial_states[[i]][['I_00_f']],
                                     I_000_m = 0,
                                     I_000_l = initial_states[[i]][['I_00_l']],
                                     
                                     C_100_u = initial_states[[i]][['C_10_u']],
                                     C_100_f = initial_states[[i]][['C_10_f']],
                                     
                                     I_100_f = initial_states[[i]][['I_10_f']],
                                     I_100_m = 0,
                                     I_100_l = initial_states[[i]][['I_10_l']],
                                     
                                     
                                     
                                     C_010_u = initial_states[[i]][['C_01_u']],
                                     C_010_f = initial_states[[i]][['C_01_f']],
                                     
                                     I_010_f = initial_states[[i]][['I_01_f']],
                                     I_010_m = 0,
                                     I_010_l = initial_states[[i]][['I_01_l']],
                                     
                                     C_001_u = 0,
                                     C_001_f = 0,
                                     
                                     I_001_f = 0,
                                     I_001_m = 0,
                                     I_001_l = 0,
                                     
                                     C_110_u = initial_states[[i]][['C_11_u']],
                                     C_110_f = initial_states[[i]][['C_11_f']],
                                     
                                     I_110_f = initial_states[[i]][['I_11_f']],
                                     I_110_m = 0,
                                     I_110_l = initial_states[[i]][['I_11_l']],
                                     
                                     C_101_u = 0,
                                     C_101_f = 0,
                                     
                                     I_101_f = 0,
                                     I_101_m = 0,
                                     I_101_l = 0,
                                     
                                     C_011_u = 0,
                                     C_011_f = 0,
                                     
                                     I_011_f = 0,
                                     I_011_m = 0,
                                     I_011_l = 0,
                                     
                                     C_111_u = 0,
                                     C_111_f = 0,
                                     
                                     I_111_f = 0,
                                     I_111_m = 0,
                                     I_111_l = 0,
                                     
                                     D = initial_states[[i]][['D']],
                                     dTrans = initial_states[[i]][['dTrans']])
  
  
}



# Run the model, record the output ----
# now run all the models on this table of parameters, using the initial conditions after the burn-in phase
# model requiring adapted initial conditions use those
# keep time step of 0.05 to avoid calculation issues with small values

# results are stored in a list. Each element of the list correspond to a simulation (i simulations)
# Within each element of the list, there is one entry per model (no_int, adj_f, adj_l, diag_f, diag_l, new_f, new_l, swap_f, swap_l)

times<- seq(0, 600, 0.05)   # time with 0.05 step sizes
records<- seq(0, tail(times,1), 1)  # to record only 300 points

# To store simulations
runs<- vector('list', length = nrow(pars))


# To store values of interest
#relative_usage<- data.frame (no_int = numeric(), swap_f = numeric(), swap_l = numeric(), adj_f = numeric(), adj_l = numeric(), diag_f = numeric(), diag_l = numeric(), new_f = numeric(), new_m = numeric(), new_l = numeric())
#total_innapropriate<- data.frame (adj_f = numeric(), adj_l = numeric(), diag_f = numeric(), diag_l = numeric(), new_f = numeric(), new_l = numeric())

total_prevalence<- data.frame (no_int = numeric(), adj_f = numeric(), adj_l = numeric(), diag_f = numeric(), diag_l = numeric(), new_f = numeric(), new_l = numeric(), new_m = numeric())
total_death<- data.frame (no_int = numeric(), adj_f = numeric(), adj_l = numeric(), diag_f = numeric(), diag_l = numeric(), new_f = numeric(), new_l = numeric(), new_m = numeric())

total_trans<- data.frame (no_int = numeric(), adj_f = numeric(), adj_l = numeric(), diag_f = numeric(), diag_l = numeric(), new_f = numeric(), new_l = numeric(), new_m = numeric())

infection_lenght<- data.frame (no_int = numeric(), adj_f = numeric(), adj_l = numeric(), diag_f = numeric(), diag_l = numeric(), new_f = numeric(), new_l = numeric(), new_m = numeric())
proba_death<- data.frame (no_int = numeric(), adj_f = numeric(), adj_l = numeric(), diag_f = numeric(), diag_l = numeric(), new_f = numeric(), new_l = numeric(), new_m = numeric())


for(i in 1:nrow(pars))
{
  
  # NO INTERVENTION: mod1, ini = initial_states
  runs[[i]]$no_int<- ode(y = initial_states[[i]], times = times, func = mod1, parms = pars[i,])
  runs[[i]]$no_int<- runs[[i]]$no_int[runs[[i]]$no_int[,'time'] %in% records,]
  
  
  
  # ADJUVANT: mod2 and mod3, ini = normal initial_states
  runs[[i]]$adj_f<-  ode(y = initial_states[[i]], times = times, func = mod2, parms = pars[i,])
  runs[[i]]$adj_f<- runs[[i]]$adj_f[runs[[i]]$adj_f[,'time'] %in% records,]
  
  
  runs[[i]]$adj_l<-  ode(y = initial_states[[i]], times = times, func = mod3, parms = pars[i,])
  runs[[i]]$adj_l<- runs[[i]]$adj_l[runs[[i]]$adj_l[,'time'] %in% records,]
  
  # DIAGNOSTIC: mod4 and mod5, ini = initial_states_front_diag and initial_states_last_diag
  runs[[i]]$diag_f<- ode(y = initial_states_front_diag[[i]], times = times, func = mod4, parms = pars[i,])
  runs[[i]]$diag_f<- runs[[i]]$diag_f[runs[[i]]$diag_f[,'time'] %in% records,]
  
  
  runs[[i]]$diag_l<- ode(y = initial_states_last_diag[[i]], times = times, func = mod5, parms = pars[i,])
  runs[[i]]$diag_l<- runs[[i]]$diag_l[runs[[i]]$diag_l[,'time'] %in% records,]
  
  
  # NEW DRUG: mod6, initial_states_front_new, initial_states_middle_new, initial_states_last_new
  runs[[i]]$new_f<-  ode(y = initial_states_front_new[[i]], times = times, func = mod6, parms = pars[i,])
  runs[[i]]$new_f<- runs[[i]]$new_f[runs[[i]]$new_f[,'time'] %in% records,]
  
  
  runs[[i]]$new_l<-  ode(y = initial_states_last_new[[i]], times = times, func = mod6, parms = pars[i,])
  runs[[i]]$new_l<- runs[[i]]$new_l[runs[[i]]$new_l[,'time'] %in% records,]
  
  runs[[i]]$new_m<-  ode(y = initial_states_middle_new[[i]], times = times, func = mod6, parms = pars[i,])
  runs[[i]]$new_m<- runs[[i]]$new_m[runs[[i]]$new_m[,'time'] %in% records,]
  
  
  
  # PART2 2: Compute values of interest
  
  
  
  # Total PREVALENCE over the course of the simulations
  # integral of the total number of infected people (whatever resistance profile), in the ILL (treated) class
  treated<- c('I_00_l', 'I_10_l', 'I_01_l', 'I_11_l', 'I_00_f', 'I_10_f', 'I_01_f', 'I_11_f')
  treated_new<- c('I_000_f', 'I_100_f', 'I_010_f', 'I_001_f', 'I_110_f', 'I_101_f', 'I_011_f', 'I_111_f',
                  'I_000_m', 'I_100_m', 'I_010_m', 'I_001_m', 'I_110_m', 'I_101_m', 'I_011_m', 'I_111_m',
                  'I_000_l', 'I_100_l', 'I_010_l', 'I_001_l', 'I_110_l', 'I_101_l', 'I_011_l', 'I_111_l')
  
  
  total_prevalence[i,]<- c(integral(runs[[i]]$no_int, treated),
                           integral(runs[[i]]$adj_f, treated),
                           integral(runs[[i]]$adj_l,  treated),
                           integral(runs[[i]]$diag_f, treated),
                           integral(runs[[i]]$diag_l, treated),
                           integral(runs[[i]]$new_f, treated_new),
                           integral(runs[[i]]$new_l, treated_new),
                           integral(runs[[i]]$new_m, treated_new))
  
  
  # Total number of DEATH
  # Integral of the number of deads of the course of the simulation
  deads<- 'D'       
  total_death[i,]<- c(integral(runs[[i]]$no_int, deads),
                      integral(runs[[i]]$adj_f, deads),
                      integral(runs[[i]]$adj_l,  deads),
                      integral(runs[[i]]$diag_f, deads),
                      integral(runs[[i]]$diag_l, deads),
                      integral(runs[[i]]$new_f, deads),
                      integral(runs[[i]]$new_l, deads),
                      integral(runs[[i]]$new_m, deads))
  
  
  
  # No. TRANSMISSONS
  trans<- 'dTrans'       
  total_trans[i,]<- c(integral(runs[[i]]$no_int, trans),
                      integral(runs[[i]]$adj_f, trans),
                      integral(runs[[i]]$adj_l,  trans),
                      integral(runs[[i]]$diag_f, trans),
                      integral(runs[[i]]$diag_l, trans),
                      integral(runs[[i]]$new_f, trans),
                      integral(runs[[i]]$new_l, trans),
                      integral(runs[[i]]$new_m, trans))
  
  
  # INFECTION LENGHT
  infection_lenght[i,]<- total_prevalence[i,] / total_trans[i,]
  
  
  # PROBA DEATH
  proba_death[i,]<- total_death[i,] / total_trans[i,]
  
  
  save.image(paste('sims_nb_', nb_sim, '_seed_', my_seed, '.RData', sep=''))
  print(paste('run ',i, ' done!', sep=''))
  
  
}