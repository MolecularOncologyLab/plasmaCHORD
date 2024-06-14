###########  Testing CHIP_model_list_allTraining.rds with 
# Razavi EGA data

library(tidyr) 
library(dplyr) 
library(magrittr)
library(stringr)
library(caret)
library(naivebayes)
library(xgboost)
library(readxl)
library(writexl)


EGAdir<-file.path("Path/to/Data/", "EGAdata")

model_list<- readRDS(file=file.path(labdir, "CHIP_model_list_allTraining.rds"))


### Pulling in fragment-level features
allSumStats_EGA<-readRDS(file.path(EGAdir, "allSumStats_EGA.rds"))

allSumStats_EGA2<-filter(allSumStats_EGA, !is.na(actual_nmut))%>%
  select(colnames(allSumStats_EGA)[!str_detect(colnames(allSumStats_EGA), "_p")])%>%
  mutate(across(c(actual_nmut, actual_nwild, width_T_stat, width_F_stat,
                  width_Md_mut, width_Md_wild, width_KS_stat,
                  width_MannWhitneyU_stat, absAveDS, AUCdeltaS, maxDS,
                  motif_Wilcox1_stat, motif_Wilcox2_stat, cut_Wilcox_stat),
                as.numeric))%>%
  mutate(plasma.vaf=MAF)%>%
  group_by(patient.id)%>%
  mutate(mean_pt_vaf=mean(plasma.vaf))%>%
  ungroup()%>%
  mutate(plasma.scaled.vaf=plasma.vaf/mean_pt_vaf)


## Pulling in patient-level features (i.e. age)
agedata<-read_xlsx(file.path(EGAdir, "Razavi Fig 4 source data.xlsx"),
                   sheet="Fig_4c")

agedata2<-agedata%>%select(patient.id=patient_id, tissue, age)%>%
  distinct()



## Pulling in variant-level features from open cravat data
OCdata<-read_xlsx(file.path(EGAdir, "EGA_OCoutput.xlsx"),
                  sheet="addl_annot", skip=6)

OCdata2<-OCdata%>%
  select(Chr=Chrom...49, Start_hg19=Pos, Base_from=Reference_allele, 
         Base_to=Alternate_allele, `CHASMplus.MSK-IMPACT.Score`=Score...19,
         `COSMIC.Variant.Count`=Variant_Count, 
         `COSMIC.Variant.Count.Tissue`=`Variant_Count (Tissue)`, 
         `COSMIC.Gene.Occurrences`=Occurrences,
         `COSMIC.Gene.Tissue.count`=Tissue_count)

# Put together fragment, cosmic, pt level data

afterP<-function(x,p){
  str_split_fixed(x, p, 2)[,2]
}

extractHemeCount<-function(x){
  if (!str_detect(x, "haematopoietic") | is.na(x)){
    hemecount<-0
  }else {
    end<-afterP(x, "haematopoietic and lymphoid tissue\",")
    hemecount<-as.numeric(str_extract(end, "[0-9]+"))
  }
  return (hemecount)
}

combined_EGA<-inner_join(OCdata2, allSumStats_EGA2)%>%
  mutate(Variant.Heme.Count=sapply(COSMIC.Variant.Count.Tissue, extractHemeCount,
                                   USE.NAMES=FALSE))%>%
  mutate(Gene.Heme.Count=sapply(COSMIC.Gene.Tissue.count, extractHemeCount,
                                USE.NAMES=FALSE))%>%
  mutate(COSMIC.Variant.Count=as.numeric(COSMIC.Variant.Count))%>%
  mutate(Gene.Heme.Frac=Gene.Heme.Count/COSMIC.Gene.Occurrences)%>%
  mutate(Variant.Heme.Frac=Variant.Heme.Count/COSMIC.Variant.Count)%>%
  replace_na(list("Variant.Heme.Frac"=-1, "COSMIC.Variant.Count"=0,
                  "Variant.Heme.Count"=0, "Gene.Heme.Frac"=-1,
                  "width_Md_mut"=-1))%>%
  left_join(agedata2)%>%
  filter(!is.na(age))%>%
mutate(sig=paste0(Base_from, "->", Base_to))%>%
  mutate(indel=(Base_from=="NA" | Base_to=="NA"))%>%
  mutate(driver=if_else(`CHASMplus.MSK-IMPACT.Score`>=0.75 | indel,1,0))%>%
  mutate(signature=if_else((indel | str_length(sig)>4), "indel_complex", sig))%>%
  mutate(signature=recode(signature, 
                          "G->T"="C->A",
                          "G->C"="C->G",
                          "G->A"="C->T",
                          "A->T"="T->A",
                          "A->G"="T->C",
                          "T->G"="A->C"))%>%
  select(-indel, -sig)




### Change into correct format for inputting into model
features<-colnames(model_list$xgb$call$x)
base_features<-c(features[!str_detect(features, "signature")], "signature")


EGA_x<-select(combined_EGA, all_of(base_features))

dmy<-dummyVars("~.", data=EGA_x)
EGA_x2<-data.frame(predict(dmy, newdata=EGA_x))%>%
  mutate(`signatureindel_complex`=0)%>%
  select(all_of(features))

EGA_y<-recode(combined_EGA$Mutation_type,
              "WBC_matched"="WBC",
              "IMPACT-BAM_matched"="Tumor",
              "biopsy_matched"="Tumor")
  

## Test model performance OVERALL
confusionMatrix(predict(model_list$xgb,EGA_x2), as.factor(EGA_y))

EGApred_call<-predict(model_list$xgb,EGA_x2)
EGApred<-predict(model_list$xgb, EGA_x2, type="prob")

EGA_data_preds<-cbind(combined_EGA, EGApred, EGApred_call, Actual=EGA_y)

write_xlsx(EGA_data_preds, file.path(EGAdir, "EGA_data_FinalPreds.xlsx"))


## Examine validation performance by cancer type
breast<-combined_EGA$tissue=="Breast"
NSCLC<-combined_EGA$tissue=="Lung"
prostate<-combined_EGA$tissue=="Prostate"

confusionMatrix(predict(model_list$xgb, EGA_x2[breast, ]), as.factor(EGA_y[breast]))
confusionMatrix(predict(model_list$xgb, EGA_x2[NSCLC, ]), as.factor(EGA_y[NSCLC]))
confusionMatrix(predict(model_list$xgb, EGA_x2[prostate, ]), as.factor(EGA_y[prostate]))

## Examine validation performance by gene

EGAgenes<-combined_EGA %>%
  select(Gene, `Mutation_type`)%>%
  mutate(`Variant Origin`=recode(combined_EGA$Mutation_type,
          "WBC_matched"="WBC",
          "IMPACT-BAM_matched"="Tumor",
          "biopsy_matched"="Tumor"))%>%
  select(-`Mutation_type`)%>%
  count(Gene, `Variant Origin`)%>%
  pivot_wider(names_from=`Variant Origin`, values_from = n)%>%
  replace_na(list(Tumor=0, WBC=0))%>%
  mutate(allcount=Tumor+WBC)%>%
  arrange(desc(allcount))

TP53<-combined_EGA$Gene=="TP53"
ATM<-combined_EGA$Gene=="ATM"
KMT2D<-combined_EGA$Gene=="KMT2D"
PPM1D<-combined_EGA$Gene=="PPM1D"
GNAS<-combined_EGA$Gene=="GNAS"

trickygenes<-TP53 | ATM | KMT2D | PPM1D | GNAS

confusionMatrix(predict(model_list$xgb, EGA_x2[TP53, ]), as.factor(EGA_y[TP53]))
confusionMatrix(predict(model_list$xgb, EGA_x2[ATM, ]), as.factor(EGA_y[ATM]))
confusionMatrix(predict(model_list$xgb, EGA_x2[KMT2D, ]), as.factor(EGA_y[KMT2D]))
confusionMatrix(predict(model_list$xgb, EGA_x2[PPM1D, ]), as.factor(EGA_y[PPM1D]))
confusionMatrix(predict(model_list$xgb, EGA_x2[GNAS, ]), as.factor(EGA_y[GNAS]))
confusionMatrix(predict(model_list$xgb, EGA_x2[trickygenes, ]), as.factor(EGA_y[trickygenes]))

chipblacklist<-c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "JAK2",
                 "RUNX1", "SF3B1", "SRSF2", "IDH1", "IDH2", "U2AF1",
                 "CBL", "ATM", "CHEK2")

chipBL<-combined_EGA$Gene %in% chipblacklist

confusionMatrix(predict(model_list$xgb, EGA_x2[chipBL, ]), as.factor(EGA_y[chipBL]))

notDAT<- ! combined_EGA$Gene %in% c("DNMT3A", "TET2", "ASXL1")

confusionMatrix(predict(model_list$xgb, EGA_x2[notDAT, ]), as.factor(EGA_y[notDAT]))

notDNMT3a<- combined_EGA$Gene != "DNMT3A"

confusionMatrix(predict(model_list$xgb, EGA_x2[notDNMT3a, ]), as.factor(EGA_y[notDNMT3a]))


###########  Examining sens, spec, accuracy for different cutoffs ##########

accsenspec<-function( cutoff){
  d<-EGApred%>%
    mutate(pred=if_else(Tumor>=cutoff, "Tumor", "WBC"))
  cm<-confusionMatrix(as.factor(d$pred), as.factor(EGA_y))
  acc<-cm$overall["Accuracy"]
  NIR<-cm$overall["AccuracyNull"]
  sens<-cm$byClass["Sensitivity"]
  spec<-cm$byClass["Specificity"]
  return(data.frame(cutoff=cutoff, acc, NIR, sens, spec))
  
}

cutoffdata<-sapply(seq(0.01, 0.99, 0.01), accsenspec)%>%
  apply(2, unlist)%>%
  as.data.frame(row.names = NULL)%>%
  t()%>%
  as.data.frame(row.names = NULL)

library(ggpubr)

statcolors<-c("Sensitivity"="red", "Specificity"="blue", "Accuracy"="black")
Acc_cutoff_EGAplot<-ggplot(data=cutoffdata, aes(x=cutoff))+
  geom_point(aes(y=acc, color="Accuracy"))+
  geom_point(aes(y=sens, color="Sensitivity"))+
  geom_point(aes(y=spec, color="Specificity"))+
  theme_classic2()+
  labs(x="Cutoff", y="", color="Model performance", title="Validation cohort")+
  scale_color_manual(values=statcolors)

#pdf(file.path(labdir, "accuracy_cutoff_EGAplot.pdf"))  
Acc_cutoff_EGAplot
#dev.off()





