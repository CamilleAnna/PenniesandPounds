# Pound&Pennies review - Functions

# LAST EDIT: 2018/06/08
# fixed bug:  dI_001_m = a*I_001_f - I_001_m*(a+z+g+d+c) + mu*I_000_m --> should be - mu*I_001_m

# EDIT: 2018/06/05: 
# edit by Luke changing mu*I_00_f to mu*I_00_l
# adding dTrans variable which counts the total number of transmission event that occured (! it is actually different for diag and new model cause transmissons are not the same)

# EDIT: 2018/06/04
# having the intervention in carrier class in adjuvant and new drug model + not in diagnostic cause it makes no sense

#EDIT: 2018/04/18
# modifying definition of Ill and Carrier compartments. Instead of a constant whole population size with proportions f and 1-f, 
# the Ill class is now of constant size 1 and the carrier class of variable size f (which is not a proportion anymore then)

# Implementing model functions, new version with misusage

# No intervention, two drugs, stratified, with mis-usage
mod1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dC_00_u  = (1 - p) * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_u * d
    dC_00_f  = p * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_f * (g+d) - mu * C_00_f
    dI_00_f  = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_00_f * (a+z+g+d) - mu * I_00_f
    dI_00_l  = a * I_00_f - I_00_l * (z+g+d) - mu*I_00_l ##### edit by Luke changing mu*I_00_f to mu*I_00_l
    
    dC_10_u  = (1 - p) * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_u * (d+c)
    dC_10_f  = p * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_f * (d+c) + mu * C_00_f
    dI_10_f  = (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_10_f * (a+z+d+c) + mu * I_00_f
    dI_10_l  = a * I_10_f - I_10_l * (z+g+d+c) - mu*I_10_l
    
    dC_01_u  = (1 - p) * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_u * (d+c)
    dC_01_f  = p * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_f * (g+d+c) - mu * C_01_f
    dI_01_f  = (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_01_f * (a+z+g+d+c) - mu * I_01_f
    dI_01_l  = a * I_01_f - I_01_l * (z+d+c) + mu*I_00_l
    
    dC_11_u  = (1 - p) * (bw * (C_11_u + C_11_f) + bb * (I_11_f + I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_u * (d+2*c)
    dC_11_f  = p * (bw * (C_11_u + C_11_f) + bb * (I_11_f + I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_f * (d+2*c) + mu * C_01_f
    dI_11_f  = (bw * (I_11_f + I_11_l) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_11_f * (a+z+d+2*c) + mu * I_01_f
    dI_11_l  = a * I_11_f - I_11_l * (z+d+2*c) + mu*I_10_l
    
    dD = z * (I_00_f + I_10_f + I_01_f + I_11_f + I_00_l + I_10_l + I_01_l + I_11_l)
    
    dTrans = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) + 
             (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_11_f + I_11_l) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l)
    
    
    list(c(dC_00_u, dC_00_f, dI_00_f, dI_00_l, 
           dC_10_u, dC_10_f, dI_10_f, dI_10_l,
           dC_01_u, dC_01_f, dI_01_f, dI_01_l,
           dC_11_u, dC_11_f, dI_11_f, dI_11_l,
           dD, dTrans))
  })
}

# front line adjuvant, two drugs, stratified, with mis-usage
mod2 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dC_00_u  = (1 - p) * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_u * d
    dC_00_f  = p * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_f * (g+d) - mu * C_00_f
    dI_00_f  = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_00_f * (a+z+g+d) - mu * I_00_f
    dI_00_l  = a * I_00_f - I_00_l * (z+g+d) - mu*I_00_l
    
    dC_10_u  = (1 - p) * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_u * (d+c)
    dC_10_f  = p * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_f * (d+c+e) + mu * C_00_f # ** intervention also occurs in carrier class
    dI_10_f  = (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_10_f * (a+z+d+c+e) + mu * I_00_f
    dI_10_l  = a * I_10_f - I_10_l * (z+g+d+c) - mu*I_10_l
    
    dC_01_u  = (1 - p) * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_u * (d+c)
    dC_01_f  = p * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_f * (g+d+c) - mu * C_01_f
    dI_01_f  = (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_01_f * (a+z+g+d+c) - mu * I_01_f
    dI_01_l  = a * I_01_f - I_01_l * (z+d+c) + mu*I_00_l
    
    dC_11_u  = (1 - p) * (bw * (C_11_u + C_11_f) + bb * (I_11_f + I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_u * (d+2*c)
    dC_11_f  = p * (bw * (C_11_u + C_11_f) + bb * (I_11_f + I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_f * (d+2*c+e) + mu * C_01_f # **
    dI_11_f  = (bw * (I_11_f + I_11_l) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_11_f * (a+z+d+2*c+e) + mu * I_01_f
    dI_11_l  = a * I_11_f - I_11_l * (z+d+2*c) + mu*I_10_l
    
    dD = z * (I_00_f + I_10_f + I_01_f + I_11_f + I_00_l + I_10_l + I_01_l + I_11_l)
    
    dTrans = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) + 
             (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_11_f + I_11_l) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l)
    
    
    list(c(dC_00_u, dC_00_f, dI_00_f, dI_00_l, 
           dC_10_u, dC_10_f, dI_10_f, dI_10_l,
           dC_01_u, dC_01_f, dI_01_f, dI_01_l,
           dC_11_u, dC_11_f, dI_11_f, dI_11_l,
           dD, dTrans))
  })
}

# last line adjuvant, two drugs, stratified, with mis-usage
mod3 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dC_00_u  = (1 - p) * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_u * d
    dC_00_f  = p * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_f * (g+d) - mu * C_00_f
    dI_00_f  = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_00_f * (a+z+g+d) - mu * I_00_f
    dI_00_l  = a * I_00_f - I_00_l * (z+g+d) - mu*I_00_l
    
    dC_10_u  = (1 - p) * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_u * (d+c)
    dC_10_f  = p * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_f * (d+c) + mu * C_00_f
    dI_10_f  = (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_10_f * (a+z+d+c) + mu * I_00_f
    dI_10_l  = a * I_10_f - I_10_l * (z+g+d+c) - mu*I_10_l
    
    dC_01_u  = (1 - p) * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_u * (d+c)
    dC_01_f  = p * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_f * (g+d+c) - mu * C_01_f
    dI_01_f  = (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_01_f * (a+z+g+d+c) - mu * I_01_f
    dI_01_l  = a * I_01_f - I_01_l * (z+d+c+e) + mu*I_00_l
    
    dC_11_u  = (1 - p) * (bw * (C_11_u + C_11_f) + bb * (I_11_f + I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_u * (d+2*c)
    dC_11_f  = p * (bw * (C_11_u + C_11_f) + bb * (I_11_f + I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_f * (d+2*c) + mu * C_01_f
    dI_11_f  = (bw * (I_11_f + I_11_l) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_11_f * (a+z+d+2*c) + mu * I_01_f
    dI_11_l  = a * I_11_f - I_11_l * (z+d+2*c+e) + mu*I_10_l
    
    dD = z * (I_00_f + I_10_f + I_01_f + I_11_f + I_00_l + I_10_l + I_01_l + I_11_l)
    
    
    dTrans = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) + 
             (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_11_f + I_11_l) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l - I_01_l - I_11_l)
    
    
    list(c(dC_00_u, dC_00_f, dI_00_f, dI_00_l, 
           dC_10_u, dC_10_f, dI_10_f, dI_10_l,
           dC_01_u, dC_01_f, dI_01_f, dI_01_l,
           dC_11_u, dC_11_f, dI_11_f, dI_11_l,
           dD, dTrans))
  })
}


# front line diagnostic, two drugs, stratified, with mis-usage
mod4 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dC_00_u  = (1 - p) * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_u * d
    dC_00_f  = p * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_f * (g+d) - mu * C_00_f
    dI_00_f  = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_01_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_00_f * (a+z+g+d) - mu * I_00_f
    dI_00_l  = a * I_00_f - I_00_l * (z+g+d) - mu*I_00_l
    
    dC_10_u  = (1 - p) * (bw * (C_10_u + C_10_f) + bb * ( I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_u * (d+c)
    dC_10_f  = p * (bw * (C_10_u + C_10_f) + bb * (I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_f * (d+c) + mu * C_00_f
    dI_10_f  = 0
    dI_10_l  = (bw * (I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_01_f  - I_00_l - I_10_l - I_01_l - I_11_l)- I_10_l * (z+g+d+c) - mu*I_10_l + mu*I_00_f
    
    dC_01_u  = (1 - p) * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_u * (d+c)
    dC_01_f  = p * (bw * (C_01_u + C_01_f) + bb * (I_01_f + I_01_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_f * (g+d+c) - mu * C_01_f
    dI_01_f  = (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_01_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_01_f * (a+z+g+d+c) - mu * I_01_f
    dI_01_l  = a * I_01_f - I_01_l * (z+d+c) + mu*I_00_l
    
    dC_11_u  = (1 - p) * (bw * (C_11_u + C_11_f) + bb * (I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_u * (d+2*c)
    dC_11_f  = p * (bw * (C_11_u + C_11_f) + bb * ( I_11_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_f * (d+2*c) + mu * C_01_f
    dI_11_f  =  0
    dI_11_l  = (bw * (I_11_l) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_01_f - I_00_l - I_10_l - I_01_l - I_11_l) - I_11_l * (z+d+2*c) + mu*I_10_l + mu*I_01_f
    
    dD = z * (I_00_f + I_01_f + I_00_l + I_10_l + I_01_l + I_11_l)
    
    dTrans = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_01_f - I_00_l - I_10_l - I_01_l - I_11_l) + 
             (bw * (I_10_l) + bb * (C_10_u + C_10_f))          * (1 - I_00_f - I_01_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_01_f + I_01_l) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_01_f - I_00_l - I_10_l - I_01_l - I_11_l) +
             (bw * (I_11_l) + bb * (C_11_u + C_11_f))          * (1 - I_00_f - I_01_f - I_00_l - I_10_l - I_01_l - I_11_l)
    
    
    list(c(dC_00_u, dC_00_f, dI_00_f, dI_00_l, 
           dC_10_u, dC_10_f, dI_10_f, dI_10_l,
           dC_01_u, dC_01_f, dI_01_f, dI_01_l,
           dC_11_u, dC_11_f, dI_11_f, dI_11_l,
           dD, dTrans))
  })
}

# Last line diagnostic, two drugs, stratified, with mis-usage
mod5 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dC_00_u  = (1 - p) * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_u * d
    dC_00_f  = p * (bw * (C_00_u + C_00_f) + bb * (I_00_f + I_00_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_00_f * (g+d) - mu * C_00_f
    dI_00_f  = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l) - I_00_f * (a+z+g+d) - mu * I_00_f
    dI_00_l  = a * I_00_f - I_00_l * (z+g+d) - mu*I_00_l
    
    dC_10_u  = (1 - p) * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_u * (d+c)
    dC_10_f  = p * (bw * (C_10_u + C_10_f) + bb * (I_10_f + I_10_l)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_10_f * (d+c) + mu * C_00_f
    dI_10_f  = (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l) - I_10_f * (a+z+d+c) + mu * I_00_f
    dI_10_l  = a * I_10_f - I_10_l * (z+g+d+c) - mu*I_10_l
    
    dC_01_u  = (1 - p) * (bw * (C_01_u + C_01_f) + bb * (I_01_f)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_u * (d+c)
    dC_01_f  = p * (bw * (C_01_u + C_01_f) + bb * (I_01_f)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_01_f * (g+d+c) - mu * C_01_f
    dI_01_f  = (bw * (I_01_f) + bb * (C_01_u + C_01_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l) - I_01_f * (z+g+d+c) - mu * I_01_f + mu*I_00_l
    dI_01_l  = 0
    
    dC_11_u  = (1 - p) * (bw * (C_11_u + C_11_f) + bb * (I_11_f )) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_u * (d+2*c)
    dC_11_f  = p * (bw * (C_11_u + C_11_f) + bb * (I_11_f)) * (f - C_00_u - C_10_u - C_01_u - C_11_u - C_00_f - C_10_f - C_01_f - C_11_f) - C_11_f * (d+2*c) + mu * C_01_f
    dI_11_f  = (bw * (I_11_f) + bb * (C_11_u + C_11_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l) - I_11_f * (z+d+2*c) + mu * I_01_f + mu*I_10_l
    dI_11_l  = 0
    
    dD = z * (I_00_f + I_10_f + I_01_f + I_11_f + I_00_l + I_10_l)
    
    dTrans = (bw * (I_00_f + I_00_l) + bb * (C_00_u + C_00_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l) + 
             (bw * (I_10_f + I_10_l) + bb * (C_10_u + C_10_f)) * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l) +
             (bw * (I_01_f) + bb * (C_01_u + C_01_f))          * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l) +
             (bw * (I_11_f) + bb * (C_11_u + C_11_f))          * (1 - I_00_f - I_10_f - I_01_f - I_11_f - I_00_l - I_10_l)
    
    
    list(c(dC_00_u, dC_00_f, dI_00_f, dI_00_l, 
           dC_10_u, dC_10_f, dI_10_f, dI_10_l,
           dC_01_u, dC_01_f, dI_01_f, dI_01_l,
           dC_11_u, dC_11_f, dI_11_f, dI_11_l,
           dD, dTrans))
  })
}


# New drug on market, three drugs, stratified, with mis-usage
mod6 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dC_000_u = (1-p) * ( ( bw *(C_000_u + C_000_f)) + ( bb *(I_000_f + I_000_m + I_000_l)) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_000_u*d
    dC_000_f =  p * ( ( bw *(C_000_u + C_000_f)) + ( bb *(I_000_f + I_000_m + I_000_l)) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_000_f*(g+d) - mu*C_000_f
    dI_000_f =  ( ( bw *(I_000_f + I_000_m + I_000_l)) + ( bb *(C_000_u + C_000_f)) ) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_000_f*(a+z+g+d) - mu*I_000_f 
    dI_000_m = a*I_000_f - I_000_m*(a+z+g+d) - mu*I_000_m
    dI_000_l = a*I_000_m - I_000_l*(z+g+d) - mu*I_000_l
    
    dC_100_u = (1-p) * ( ( bw *(C_100_u + C_100_f)) + ( bb *(I_100_f + I_100_m + I_100_l)) ) * (f -(C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_100_u*(d+c)
    dC_100_f =  p * ( ( bw *(C_100_u + C_100_f)) + ( bb *(I_100_f + I_100_m + I_100_l)) ) * (f -(C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_100_f*(d+c) + mu*C_100_f
    dI_100_f =  ( ( bw *(I_100_f + I_100_m + I_100_l)) + ( bb *(C_100_u + C_100_f)) ) * ( 1 -(I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_100_f*(a+z+d+c) + mu*I_000_f 
    dI_100_m = a*I_100_f - I_100_m*(a+z+g+d+c) - mu*I_100_m
    dI_100_l = a*I_100_m - I_100_l*(z+g+d+c) - mu*I_100_l
    
    dC_010_u = (1-p) * ( ( bw *(C_010_u + C_010_f)) + ( bb *(I_010_f + I_010_m + I_010_l)) ) * (f -(C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_010_u*(d+c)
    dC_010_f =  p * ( ( bw *(C_010_u + C_010_f)) + ( bb *(I_010_f + I_010_m + I_010_l)) ) * (f -(C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_010_f*(g+d+c) - mu*C_010_f
    dI_010_f =  ( ( bw *(I_010_f + I_010_m + I_010_l)) + ( bb *(C_010_u + C_010_f)) ) * ( 1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_010_f*(a+z+g+d+c) - mu*I_010_f 
    dI_010_m = a*I_010_f - I_010_m*(a+z+d+c) + mu*I_000_m
    dI_010_l = a*I_010_m - I_010_l*(z+g+d+c) - mu*I_010_l + mu*I_000_l
    
    dC_001_u = (1-p) * ( ( bw *(C_001_u + C_001_f) ) + ( bb *(I_001_f + I_001_m + I_001_l) ) ) * (f -(C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_001_u*(d+c)
    dC_001_f =  p * ( ( bw *(C_001_u + C_001_f) ) + ( bb *(I_001_f + I_001_m + I_001_l) ) ) * (f -(C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_001_f*(g+d+c) - mu*C_001_f
    dI_001_f =  ( ( bw *(I_001_f + I_001_m + I_001_l) ) + ( bb *(C_001_u + C_001_f) ) ) * ( 1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_001_f*(a+z+g+d+c) - mu*I_001_f 
    dI_001_m = a*I_001_f - I_001_m*(a+z+g+d+c) - mu*I_001_m
    dI_001_l = a*I_001_m - I_001_l*(z+d+c) + mu*I_000_l
    
    
    

    dC_110_u = (1-p) * ( ( bw *(C_110_u + C_110_f) ) + ( bb *(I_110_f + I_110_m + I_110_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_110_u*(d+2*c)
    dC_110_f =  p * ( ( bw *(C_110_u + C_110_f) ) + ( bb *(I_110_f + I_110_m + I_110_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_110_f*(d+2*c) + mu*C_010_f
    dI_110_f =  ( ( bw *(I_110_f + I_110_m + I_110_l) ) + ( bb *(C_110_u + C_110_f) ) ) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_110_f*(a+z+d+2*c) + mu*I_010_f 
    dI_110_m = a*I_110_f - I_110_m*(a+z+d+2*c) + mu*I_100_m
    dI_110_l = a*I_110_m - I_110_l*(z+g+d+2*c) - mu*I_110_l
    

    
    dC_101_u = (1-p) * ( ( bw *(C_101_u + C_101_f) ) + ( bb *(I_101_f + I_101_m + I_101_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_101_u*(d+2*c)
    dC_101_f =  p * ( ( bw *(C_101_u + C_101_f) ) + ( bb *(I_101_f + I_101_m + I_101_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_101_f*(d+2*c) + mu*C_001_f
    dI_101_f =  ( ( bw *(I_101_f + I_101_m + I_101_l) ) + ( bb *(C_101_u + C_101_f) ) ) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_101_f*(a+z+d+2*c) + mu*I_001_f 
    dI_101_m = a*I_101_f - I_101_m*(a+z+g+d+2*c) - mu*I_101_m
    dI_101_l = a*I_101_m - I_101_l*(z+d+2*c) + mu*I_100_l
    
    # here
    dC_011_u = (1-p) * ( ( bw *(C_011_u + C_011_f) ) + ( bb *(I_011_f + I_011_m + I_011_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_011_u*(d+2*c)
    dC_011_f =  p * ( ( bw *(C_011_u + C_011_f) ) + ( bb *(I_011_f + I_011_m + I_011_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_011_f*(g+d+2*c) - mu*C_011_f
    dI_011_f =  ( ( bw *(I_011_f + I_011_m + I_011_l) ) + ( bb *(C_011_u + C_011_f) ) ) * ( 1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_011_f*(a+z+g+d+2*c) - mu*I_011_f 
    dI_011_m = a*I_011_f - I_011_m*(a+z+d+2*c) + mu*I_001_m
    dI_011_l = a*I_011_m - I_011_l*(z+d+2*c) + mu*I_010_l
    
    dC_111_u = (1-p) * ( ( bw *(C_111_u + C_111_f) ) + ( bb *(I_111_f + I_111_m + I_111_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_111_u*(d+3*c)
    dC_111_f =  p * ( ( bw *(C_111_u + C_111_f) ) + ( bb *(I_111_f + I_111_m + I_111_l) ) ) * (f - (C_000_u + C_000_f + C_100_u + C_100_f + C_010_u + C_010_f + C_001_u + C_001_f + C_110_u + C_110_f + C_101_u + C_101_f + C_011_u + C_011_f + C_111_u + C_111_f)) - C_111_f*(d+3*c) + mu*C_011_f
    dI_111_f =  ( ( bw *(I_111_f + I_111_m + I_111_l) ) + ( bb *(C_111_u + C_111_f) ) ) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) - I_111_f*(a+z+d+3*c) + mu*I_011_f 
    dI_111_m = a*I_111_f - I_111_m*(a+z+d+3*c) + mu*I_101_m
    dI_111_l = a*I_111_m - I_111_l*(z+d+3*c) + mu*I_110_l
    
    dD = z * (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)
    
    
    
    dTrans =  ((bw *(I_000_f + I_000_m + I_000_l)) + (bb *(C_000_u + C_000_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) +
              ((bw *(I_100_f + I_100_m + I_100_l)) + (bb *(C_100_u + C_100_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) +
              ((bw *(I_010_f + I_010_m + I_010_l)) + (bb *(C_010_u + C_010_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) +
              ((bw *(I_001_f + I_001_m + I_001_l)) + (bb *(C_001_u + C_001_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) +
              ((bw *(I_110_f + I_110_m + I_110_l)) + (bb *(C_110_u + C_110_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) +
              ((bw *(I_101_f + I_101_m + I_101_l)) + (bb *(C_101_u + C_101_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) +
              ((bw *(I_011_f + I_011_m + I_011_l)) + (bb *(C_011_u + C_011_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l)) +
              ((bw *(I_111_f + I_111_m + I_111_l)) + (bb *(C_111_u + C_111_f))) * (1 - (I_000_f + I_100_f + I_010_f + I_001_f  + I_110_f  + I_101_f  + I_011_f  + I_111_f + I_000_m + I_100_m + I_010_m + I_001_m  + I_110_m  + I_101_m  + I_011_m  + I_111_m + I_000_l + I_100_l + I_010_l + I_001_l  + I_110_l  + I_101_l  + I_011_l  + I_111_l))
    
    
    
    list(c(dC_000_u, dC_000_f, dI_000_f, dI_000_m, dI_000_l,
           dC_100_u, dC_100_f, dI_100_f, dI_100_m, dI_100_l,
           dC_010_u, dC_010_f, dI_010_f, dI_010_m, dI_010_l,
           dC_001_u, dC_001_f, dI_001_f, dI_001_m, dI_001_l,
           dC_110_u, dC_110_f, dI_110_f, dI_110_m, dI_110_l,
           dC_101_u, dC_101_f, dI_101_f, dI_101_m, dI_101_l,
           dC_011_u, dC_011_f, dI_011_f, dI_011_m, dI_011_l,
           dC_111_u, dC_111_f, dI_111_f, dI_111_m, dI_111_l,
           dD, dTrans))
  })
}







integral<- function(x, strains){ 
  out<- as.data.frame(x)
  out_sum<- data.frame(time = out$time, sum = rowSums(as.data.frame(out[,strains])))
  sum(diff(out_sum$time) * (head(out_sum$sum,-1)+tail(out_sum$sum,-1)))/2
}



