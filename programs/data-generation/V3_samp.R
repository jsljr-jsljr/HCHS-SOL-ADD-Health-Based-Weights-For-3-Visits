### Generate Samples for Simulation Study (based on plan from 10/14/16) ###
library(tidyverse)

####### function to generate sample #######
samp.gen = function(pop,prob.bg,num.bg,prob.hh.hisp,prob.hh.other,prob.age){
  
  pop.unq=pop[pop$v_num==1,]   #this is the population only including visit 2 records (to simplify sampling)  
  pop.v2=pop[pop$v_num==2,]
  pop.v3=pop[pop$v_num==3,]
  
  ### re-create hh.size ###
  hh.list=unique(pop.unq[,c("BGid","hhid")])  # list of BGid & hhid for each unique hhid
  bg.size=table(hh.list$BGid)   # number of HHs per BG, 
  hh.size=table(pop.unq$hhid)   # number of subjects per HH
  hisp_strat.unq=aggregate(hisp_strat~hhid,pop.unq,mean)$hisp_strat   # hisp_strat with 1 record per HH
  
  ### generate raw weights (these do not depend on the sample (only depend on sampling probabilities)) ###
  pop.unq$W_bg=1/prob.bg[pop.unq$strat]   # BG stage raw weight (based on BG sampling fraction; 1/sampling fraction for that stratum)
  #pop.unq$W_hh=ifelse(pop.unq$hhid %in% hh.list$hhid[hisp_strat.unq],1/prob.hh.hisp[pop.unq$strat],1/prob.hh.other[pop.unq$strat])   #HH stage raw weight
  pop.unq$W_hh=ifelse(pop.unq$hisp_strat,1/prob.hh.hisp[pop.unq$strat],1/prob.hh.other[pop.unq$strat])   #HH stage raw weight
  pop.unq$W_sub=ifelse(pop.unq$age_strat,1/prob.age[2],1/prob.age[1])    # subject stage raw weight
  pop.unq$W_bghhsub=pop.unq$W_bg*pop.unq$W_hh*pop.unq$W_sub     # raw combined weights
  pop.unq$W_bghh=pop.unq$W_bg*pop.unq$W_hh     # raw combined weights, bghh
  pop.unq$W_hhsub=pop.unq$W_hh*pop.unq$W_sub     # raw combined weights, hhsub
  
  
  ### select random sample from population ###
  #select stratified random sample of BGs from pop & save list of BGs
  # bg.select=c(sample(seq(num.bg[1]),round(prob.bg[1]*num.bg[1])),
  #             sample(seq(num.bg[1] + 1, sum(num.bg[1:2])),round(prob.bg[2]*num.bg[2])),
  #             sample(seq(sum(num.bg[1:2]) + 1, sum(num.bg[1:3])),round(prob.bg[3]*num.bg[3])),
  #             sample(seq(sum(num.bg[1:3]) + 1, sum(num.bg[1:4])),round(prob.bg[4]*num.bg[4])),
  #             sample(seq(sum(num.bg[1:4]) + 1, sum(num.bg[1:5])),round(prob.bg[5]*num.bg[5])),
  #             sample(seq(sum(num.bg[1:5]) + 1, sum(num.bg[1:6])),round(prob.bg[6]*num.bg[6])),
  #             sample(seq(sum(num.bg[1:6]) + 1, sum(num.bg[1:7])),round(prob.bg[7]*num.bg[7])),
  #             sample(seq(sum(num.bg[1:7]) + 1, sum(num.bg[1:8])),round(prob.bg[8]*num.bg[8])))
  # bg.select.s=ifelse(bg.select<=num.bg[1], 1, 
  #                    ifelse(bg.select<=sum(num.bg[1:2]), 2,
  #                           ifelse(bg.select<=sum(num.bg[1:3]), 3,
  #                                  ifelse(bg.select<=sum(num.bg[1:4]), 4,
  #                                         ifelse(bg.select<=sum(num.bg[1:5]), 5,
  #                                                ifelse(bg.select<=sum(num.bg[1:6]), 6,
  #                                                       ifelse(bg.select<=sum(num.bg[1:7]), 7, 8)))))))
    
  bg.select=c(sample(1:58,round(prob.bg[1]*num.bg[1])),sample(59:79,round(prob.bg[2]*num.bg[2])),sample(80:209,round(prob.bg[3]*num.bg[3])),sample(210:376,round(prob.bg[4]*num.bg[4])),sample(377:434,round(prob.bg[5]*num.bg[5])),sample(435:455,round(prob.bg[6]*num.bg[6])),sample(456:585,round(prob.bg[7]*num.bg[7])),sample(586:752,round(prob.bg[8]*num.bg[8])))
  bg.select.s=1*(bg.select<=58)+2*(bg.select>=59 & bg.select<=79)+3*(bg.select>=80 & bg.select<=209)+4*(bg.select>=210 & bg.select<=376)+5*(bg.select>=377 & bg.select<=434)+6*(bg.select>=435 & bg.select<=455)+7*(bg.select>=456 & bg.select<=585)+8*(bg.select>=586 & bg.select<=752)  #stratum for each BG
  
  samp.hh.bystrat=function(s){
    hh.select.s=rep(FALSE,dim(pop.unq[pop.unq$strat==s,])[1])         #initially set hh.select.s=FALSE for all subjects in stratum s, so that unselected BG's will be set to FALSE
    for (j in bg.select[bg.select.s==s]){
      hh.list.hisp=hh.list$hhid[hh.list$BGid==j & hisp_strat.unq==TRUE]   #list of unique HHs in BGid j, w/ Hisp surname
      hh.list.other=hh.list$hhid[hh.list$BGid==j & hisp_strat.unq==FALSE] #list of unique HHs in BGid j, w/ other surname
      hh.select.hisp=sample(hh.list.hisp,round(prob.hh.hisp[s]*length(hh.list.hisp)))      #list of randomly sampled HHs in BGid w/ Hispanic surname
      hh.select.other=sample(hh.list.other,round(prob.hh.other[s]*length(hh.list.other)))  #list of randomly sampled HHs in BGid w/ other surname
      hh.select.s[pop.unq$hhid[pop.unq$strat==s] %in% c(hh.select.hisp,hh.select.other)]=TRUE  #select subjects from randomly sampled HHs in BGid j (one entry per subject in BGid j)
    }
    return(hh.select.s)
  }
  hh.select=c(samp.hh.bystrat(1),samp.hh.bystrat(2),samp.hh.bystrat(3),samp.hh.bystrat(4),samp.hh.bystrat(5),samp.hh.bystrat(6),samp.hh.bystrat(7),samp.hh.bystrat(8))  #indicator of HH selection
  
  
  #select random sample of subjects from each selected HH & save indicator of subject selection
  sub.select=rep(FALSE,dim(pop.unq)[1])           # initially set sub.select=FALSE for all subjects, so that unselected HH's will be set to FALSE
  sub.select[pop.unq$subid %in% sample(pop.unq$subid[!pop.unq$age_strat & hh.select],round(prob.age[1]*dim(pop.unq[!pop.unq$age_strat & hh.select,])[1]))]=TRUE   # randomly sample younger subjects among sampled HH's
  sub.select[pop.unq$subid %in% sample(pop.unq$subid[pop.unq$age_strat & hh.select],round(prob.age[2]*dim(pop.unq[pop.unq$age_strat & hh.select,])[1]))]=TRUE     # randomly sample older subjects among sampled HH's
  
  samp=pop.unq[sub.select,]# create sample by restricting pop to selected subjects
  
  # generate normalized weight
  samp$bghhsub_s2=samp$W_bghhsub/mean(samp$W_bghhsub)    # normalized weight (raw combined weight/mean combined weight)
  
  # scaled subject weights
  sumw.sub=aggregate(W_sub ~ BGid + hhid, dat=samp, FUN=sum); sumw.sub=sumw.sub[order(sumw.sub$BGid,sumw.sub$hhid),]
  sumw.sub.2=aggregate(W_sub^2 ~ BGid + hhid, dat=samp, FUN=sum); sumw.sub.2=sumw.sub.2[order(sumw.sub.2$BGid,sumw.sub.2$hhid),]
  nw.sub=aggregate(W_sub ~ BGid + hhid, dat=samp, FUN=length); nw.sub=nw.sub[order(nw.sub$BGid,nw.sub$hhid),]
  
  samp$sub_s1=samp$W_sub*rep(sumw.sub$W_sub, times=nw.sub$W_sub)/rep(sumw.sub.2$"W_sub^2", times=nw.sub$W_sub)
  samp$sub_s2=samp$W_sub*rep(nw.sub$W_sub, times=nw.sub$W_sub)/rep(sumw.sub$W_sub, times=nw.sub$W_sub)
  
  # scaled HH weights
  samp.unqhh=samp[!duplicated(samp$hhid),]    #keep one row per HH
  sumw.hh=aggregate(W_hh ~ BGid, dat=samp.unqhh, FUN=sum)
  sumw.hh.2=aggregate(W_hh^2 ~ BGid, dat=samp.unqhh, FUN=sum)
  nw.hh=aggregate(W_hh ~ BGid, dat=samp.unqhh, FUN=length)
  
  samp$hh_s1=samp$W_hh*merge(samp[,c('BGid','hhid','subid')],sumw.hh,by='BGid',all=TRUE)$W_hh/merge(samp[,c('BGid','hhid','subid')],sumw.hh.2,by='BGid',all=TRUE)$'W_hh^2'
  samp$hh_s2=samp$W_hh*merge(samp[,c('BGid','hhid','subid')],nw.hh,by='BGid',all=TRUE)$W_hh/merge(samp[,c('BGid','hhid','subid')],sumw.hh,by='BGid',all=TRUE)$W_hh
  
  
  # scaled hh-subject weights
  sumw.hhsub=aggregate(W_hhsub ~ BGid, dat=samp, FUN=sum); sumw.hhsub=sumw.hhsub[order(sumw.hhsub$BGid),]
  sumw.hhsub.2=aggregate(W_hhsub^2 ~ BGid, dat=samp, FUN=sum); sumw.hhsub.2=sumw.hhsub.2[order(sumw.hhsub.2$BGid),]
  nw.hhsub=aggregate(W_hhsub ~ BGid, dat=samp, FUN=length); nw.hhsub=nw.hhsub[order(nw.hhsub$BGid),]
  
  samp$hhsub_s1=samp$W_hhsub*rep(sumw.hhsub$W_hhsub, times=nw.hhsub$W_hhsub)/rep(sumw.hhsub.2$"W_hhsub^2", times=nw.hhsub$W_hhsub)
  samp$hhsub_s2=samp$W_hhsub*rep(nw.hhsub$W_hhsub, times=nw.hhsub$W_hhsub)/rep(sumw.hhsub$W_hhsub, times=nw.hhsub$W_hhsub)
  
  samp_final_ <- rbind(samp, samp, samp)
  samp_final_$x6 <- pop$x6[rep(sub.select, 3)]
  samp_final_$y_bmi <- pop$y_bmi[rep(sub.select, 3)]
  samp_final_$y_gfr <- pop$y_gfr[rep(sub.select, 3)]
  
  samp_final_$y_bin_gfr_low <- pop$y_bin_gfr_low[rep(sub.select, 3)]
  samp_final_$y_bin_gfr_med <- pop$y_bin_gfr_med[rep(sub.select, 3)]
  samp_final_$y_bin_gfr_hi <- pop$y_bin_gfr_hi[rep(sub.select, 3)]
  samp_final_$y_bin_gfr <- pop$y_bin_gfr[rep(sub.select, 3)]
  
  samp_final_$miss_ind <- pop$miss_ind[rep(sub.select, 3)]
  samp_final_$miss_ind_mar <- pop$miss_ind_mar[rep(sub.select, 3)]
  
  samp_final_$v_num <- pop$v_num[rep(sub.select, 3)]
  
  samp_final <- samp_final_ %>% arrange(subid, v_num)
  
  return(samp_final)
}

####### read-in population datasets #######

pop = read.csv("/pine/scr/j/s/jsljr/HCHSSOL_AddHealth_WCs/population_3visits_June2021_missing.csv")

pop <- pop %>% arrange(v_num, subid)

####### generate S samples from same population for each scenario #######

allsamp=function(k){
  samp=samp.gen(pop,
                #prob.bg=rep(c(.25,.25,.6,.6)*1.66,times=2),
                prob.bg=rep(c(.25,.25,.6,.6),times=2),
                num.bg=rep(c(58,21,130,167),times=2),
                prob.hh.hisp=rep(c(.18,.225,.14,.14)*2.7,times=2),
                prob.hh.other=rep(c(.025,.0175,.035,.04)*2.7,times=2),
                #prob.age=c(.3575,.55),
                prob.age=c(.3575,.55))
  
  # samp$dat_num=rep(k,times=dim(samp)[1])
  samp$dat_num=k
  
  write.csv(samp,paste0("/pine/scr/j/s/jsljr/HCHSSOL_AddHealth_WCs/Original_Samples/sample_",k,".csv"),row.names=FALSE)

}

sim_batch <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(sim_batch)
n_batch <- 25

set.seed(sim_batch)

sim_start <- (sim_batch-1)*n_batch+1
sim_end <- sim_batch*n_batch

for(k in seq(from=sim_start, to=sim_end, by=1)){
  allsamp(k)
}
