library("R.utils")

res<- list()
iniz<- list()
parz<- list()
files<- list.files('/exports/csce/eddie/biology/groups/mcnally/camille/PoundsandPennies/NEW_STUFF/sims_output')



nb<- do.call('rbind', strsplit(files, '_'))[,5]
nb2<- do.call('rbind', strsplit(nb, '.R'))[,1]
t<- data.frame(files, as.numeric(nb2))
files<- as.character(t[order(t[,2]),1])

setwd('/exports/csce/eddie/biology/groups/mcnally/camille/PoundsandPennies/NEW_STUFF/sims_output')

counter = 1

for(ind in 1:length(files)){
  
load(paste(files[ind]))

SIM_ISSUE<- NULL


#for(bla in 1:nrow(total_prevalence)){
for(bla in 1:length(runs)){
SIM_ISSUE<- c(SIM_ISSUE, sum(tail(outs_burn_in[[bla]],1) < 0) > 0)
}
  
  res[[ind]]<- data.frame(
   

    tot_prevalence_no_int =  total_prevalence$no_int,
    tot_prevalence_adj_f =  total_prevalence$adj_f,
    tot_prevalence_adj_l =  total_prevalence$adj_l,
    tot_prevalence_diag_f =  total_prevalence$diag_f,
    tot_prevalence_diag_l =  total_prevalence$diag_l,
    tot_prevalence_new_f =  total_prevalence$new_f,
    tot_prevalence_new_l =  total_prevalence$new_l,
    tot_prevalence_new_m =  total_prevalence$new_m,
    
    nb_death_no_int =  total_death$no_int,
    nb_death_adj_f =  total_death$adj_f,
    nb_death_adj_l =  total_death$adj_l,
    nb_death_diag_f =  total_death$diag_f,
    nb_death_diag_l =  total_death$diag_l,
    nb_death_new_f =  total_death$new_f,
    nb_death_new_l =  total_death$new_l,
    nb_death_new_m =  total_death$new_m,
    
    
    total_trans_no_int = total_trans$no_int,
    total_trans_adj_f = total_trans$adj_f,
    total_trans_adj_l = total_trans$adj_l,
    total_trans_diag_f = total_trans$diag_f,
    total_trans_diag_l = total_trans$diag_l,
    total_trans_new_f = total_trans$new_f,
    total_trans_new_l = total_trans$new_l,
    total_trans_new_m = total_trans$new_m,
    
    
    infection_lenght_no_int = infection_lenght$no_int,
    infection_lenght_adj_f = infection_lenght$adj_f,
    infection_lenght_adj_l = infection_lenght$adj_l,
    infection_lenght_diag_f = infection_lenght$diag_f,
    infection_lenght_diag_l = infection_lenght$diag_l,
    infection_lenght_new_f = infection_lenght$new_f,
    infection_lenght_new_l = infection_lenght$new_l,
    infection_lenght_new_m = infection_lenght$new_m,
    
    
    proba_death_no_int = proba_death$no_int,
    proba_death_adj_f = proba_death$adj_f,
    proba_death_adj_l = proba_death$adj_l,
    proba_death_diag_f = proba_death$diag_f,
    proba_death_diag_l = proba_death$diag_l,
    proba_death_new_f = proba_death$new_f,
    proba_death_new_l = proba_death$new_l,
    proba_death_new_m = proba_death$new_m,

    
    prop_f_usage =  prop_f_usage[1:nrow(total_prevalence)],
    prop_f_usage_ill = prop_f_usage_ill[1:nrow(total_prevalence)],
    prop_inapropriate_usage = prop_inapropriate_usage[1:nrow(total_prevalence)],
    
    seed = pars[1:nrow(total_prevalence),'my_seed'],
    draw = pars[1:nrow(total_prevalence),'draw'],
    
    sim_issue = SIM_ISSUE[1:nrow(total_prevalence)]
    
    )

  
  
    
 
 iniz[[ind]]<- as.data.frame(do.call('rbind', initial_states))[1:nrow(total_prevalence),]
 iniz[[ind]]$sim_issue <- SIM_ISSUE[1:nrow(total_prevalence)]

 
 parz[[ind]]<- pars[1:nrow(total_prevalence),]
 parz[[ind]]$sim_issue <- SIM_ISSUE[1:nrow(total_prevalence)]


benefit<- do.call('rbind', res)
write.table(benefit, 
'/exports/csce/eddie/biology/groups/mcnally/camille/PoundsandPennies/NEW_STUFF/benefits_assembled.txt', col.names = 
TRUE, row.names = FALSE, quote = FALSE, sep='\t')

 initial_cond<- do.call('rbind', iniz)
write.table(initial_cond, 
'/exports/csce/eddie/biology/groups/mcnally/camille/PoundsandPennies/NEW_STUFF/initial_states_assembled.txt', 
col.names = TRUE, row.names = FALSE, quote = FALSE, sep='\t')

paramz<- do.call('rbind', parz)
write.table(paramz, 
'/exports/csce/eddie/biology/groups/mcnally/camille/PoundsandPennies/NEW_STUFF/pars_assembled.txt', col.names 
= TRUE, row.names = FALSE, quote = FALSE, sep='\t')


 print(ind) 
 
}

