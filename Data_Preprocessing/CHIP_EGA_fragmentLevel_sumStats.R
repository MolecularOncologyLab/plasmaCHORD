#CHIP fragment level summary statistics for Razavi EGA data set



library(tidyr) 
library(dplyr) 
library(magrittr)
library(readxl)
library(writexl)
library(stringr)
library(pracma)
library(ggplot2)


EGAdir<-"I:/Path/to/Data"


serialcohort<-readRDS(file.path(EGAdir, "CHIP_EGAcohort_8_20_22.rds"))

beforeP<-function(x,p){
  str_split_fixed(x, p, 2)[,1]
}

processedfiles<-tibble(filename= list.files(file.path(EGAdir,"rds") ,
                                                pattern = "nofamily_galp_processed.rds"))%>%
  mutate(file=str_replace_all(filename, c("\uf03a"=":", "\uf03e"=">")))%>%
  mutate(index=beforeP(file, "_nofamily"))%>%
  mutate(direct=file.path(EGAdir, "rds"))


filelist2<-inner_join(serialcohort, processedfiles)


Marass<-function(d=t1b, tumor=c(127:141, 272:292), CH=c(173:191, 346-361)){
  m<-d%>%
    filter(mutation=="mutant")%>%
    select(width)
  mfrac<-sum(m$width %in% tumor)/(sum(m$width %in% tumor)+sum(m$width %in% CH))
  w<-d%>%
    filter(mutation=="wild")%>%
    select(width)
  wfrac<-sum(w$width %in% tumor)/(sum(w$width %in% tumor)+sum(w$width %in% CH))
  return(c(mfrac, wfrac))
}

deltaS<-function(d=t1b, p){
  m<-d%>%
    filter(mutation=="mutant")%>%
    select(width)
  w<-d%>%
    filter(mutation=="wild")%>%
    select(width)
  ds<-vector()
  for (pi in p){
    dsi<-(sum(m$width<=pi)/dim(m)[1])-(sum(w$width<=pi)/dim(w)[1])
    ds<-c(ds, dsi)
  }
  return(ds)
}


endmotifs<-crossing(nt1=c("A", "C", "G", "T"), nt2=c("A", "C", "G", "T"),
                    nt3=c("A", "C", "G", "T"), nt4=c("A", "C", "G", "T"))%>%
  mutate(motif=paste0(nt1, nt2, nt3, nt4))%>%
  select(motif)

### LOOP through files ##########


allSumStats<-NULL

for (j in 1:dim(filelist2)[1]){
  file1<-filelist2$filename[j]
  t1<-readRDS(file.path(filelist2$direct[j], file1))
  fname<-filelist2$index[j]
  
  ## problem with MSK-VB-0006__chr12:111856097_C>T as only 1 mutant fragment
  
  t1b<-t1%>%
    filter(mutation=="mutant" | mutation=="wild")%>%
    mutate(rel_start=mut_loc - start)%>%
    mutate(rel_end=end-mut_loc-1)%>%
    filter(width<=500)
  
  if (length(unique(t1b$mutation))!=2){
    next
  }
  
  actual_nmut<-sum(t1b$mutation=="mutant")
  actual_nwild<-sum(t1b$mutation=="wild")
  
  if (actual_nmut<=2){
    next
  }
  
  # Stats looking at fragment length
  
  mutwidth<-t1b%>%
    filter(mutation=="mutant")%>%
    select(width)%>%
    arrange(width)%>%
    distinct()%>%
    mutate(deltaSvalue=deltaS(d=t1b, p=width))
  
  maxDS<-max(mutwidth$deltaSvalue)
  absAveDS<-sum(abs(mutwidth$deltaSvalue))/dim(mutwidth)[1]
  mutwidth2<-mutwidth%>%
    filter(width <=500)
  AUCdeltaS<-trapz(mutwidth2$width, mutwidth2$deltaSvalue)
  
  Ftest<-var.test(width~mutation, t1b, alternative="two.sided")
  Ttest<-t.test(width~mutation, t1b, alternative="two.sided")
  
  
  Md<-Marass(t1b)
  
  
  
  mt1b<-t1b%>%filter(mutation=="mutant")%>%
    mutate(motif1=factor(motif1, levels=endmotifs$motif))%>%
    mutate(motif2=factor(motif2, levels=endmotifs$motif))
  wt1b<-t1b%>%filter(mutation=="wild")%>%
    mutate(motif1=factor(motif1, levels=endmotifs$motif))%>%
    mutate(motif2=factor(motif2, levels=endmotifs$motif))
  
  
  KStest<-ks.test(wt1b$width, mt1b$width)
  MWUwidth<-wilcox.test(wt1b$width, mt1b$width)
  
  
  # Stats looking at endpoint motif  
  mend1<-mt1b%>%dplyr::count(motif1, .drop=FALSE)%>%
    filter(str_detect(motif1, "N", negate=TRUE))
  
  mend2<-mt1b%>%dplyr::count(motif2, .drop=FALSE)%>%
    filter(str_detect(motif2, "N", negate=TRUE))
  
  wend1<-wt1b%>%dplyr::count(motif1, .drop=FALSE)%>%
    filter(str_detect(motif1, "N", negate=TRUE))
  
  wend2<-wt1b%>%dplyr::count(motif2, .drop=FALSE)%>%
    filter(str_detect(motif2, "N", negate=TRUE))
  
  end1<-bind_rows(bind_cols(mend1, mutation="mutant"),
                  bind_cols(wend1, mutation="wild"))
  
  end2<-bind_rows(bind_cols(mend2, mutation="mutant"),
                  bind_cols(wend2, mutation="wild"))
  
  wilcox_end1<-wilcox.test(mend1$n, wend1$n)
  wilcox_end2<-wilcox.test(mend2$n, wend2$n)
  
  
  
  # Relative Endpoint analysis
  mutends<-data.frame(ends=c(-mt1b$rel_start, mt1b$rel_end)) %>%
    mutate(mutation="mutant")
  
  wtends<-data.frame(ends=c(-wt1b$rel_start, wt1b$rel_end)) %>%
    mutate(mutation="wild")
  
  allends<-rbind(mutends, wtends)
  
  allends_count<-allends%>%
    dplyr::count(ends, mutation)%>%
    add_count(mutation, wt=n)%>%
    mutate(freq=n/nn)
  
  
  cutpoint<-data.frame(cp=seq(-250,250,1))%>%
    left_join(dplyr::count(mutends, ends), by=c("cp"="ends"))%>%
    mutate(mutant=ifelse(is.na(n), 0, n))%>%
    select(-n)%>%
    left_join(dplyr::count(wtends,ends), by=c("cp"="ends"))%>%
    mutate(wild=ifelse(is.na(n), 0, n))%>%
    select(-n)
  
  wilcox_cut<-wilcox.test(cutpoint$mutant, cutpoint$wild)
  
  
  
  # Width statistics  
  allstat<-c(fname,
             actual_nmut,
             actual_nwild,
             Ttest$statistic,
             Ttest$p.value,
             Ftest$statistic,
             Ftest$p.value,
             Md,
             KStest$statistic,
             KStest$p.value,
             MWUwidth$statistic,
             MWUwidth$p.value,
             absAveDS,
             AUCdeltaS,
             maxDS,
             
             # Endpoint motif statistics
             wilcox_end1$statistic,
             wilcox_end1$p.value,
             wilcox_end2$statistic,
             wilcox_end2$p.value,
             
             # Cutpoint statistics
             wilcox_cut$statistic,
             wilcox_cut$p.value)
  
  names(allstat)<-c("index","actual_nmut", "actual_nwild",
                    "width_T_stat","width_T_p", "width_F_stat",
                    "width_F_p", "width_Md_mut", "width_Md_wild", 
                    "width_KS_stat", "width_KS_p","width_MannWhitneyU_stat",
                    "width_MannWhitneyU_p",
                    "absAveDS" ,"AUCdeltaS","maxDS",
                    "motif_Wilcox1_stat", "motif_Wilcox1_p",
                    "motif_Wilcox2_stat", "motif_Wilcox2_p",
                    "cut_Wilcox_stat", "cut_Wilcox_p")
  allSumStats<-bind_rows(allSumStats, allstat)
  
  
}



allSumStats_EGA<- filelist2%>%
  left_join(allSumStats)

ggplot(data=allSumStats_EGA, aes(x=distinct.mutant.reads, 
                                 y=as.numeric(actual_nmut)))+
  geom_point()


saveRDS(allSumStats_EGA, file.path(EGAdir, "allSumStats_EGA.rds"))

