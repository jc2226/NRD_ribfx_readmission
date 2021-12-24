####STUDY AIM #####
  #1) Understand prevalence of 3-month readmissions among geriatric patients with rib fx without polytrauma
  #2) Understand top causes of those readmissions
  #3) Understand associations between index admission charactersitics (injury pattern ,decisions) & readmission 


#### PACKAGES ####
library(data.table)
library(devtools)
devtools::install_github("ablack3/icdpicr")
devtools::install_github("vcastro/CCS")
library(icdpicr)
library(tidyr)
library(dplyr)
library(tableone)
library(glmnet)
library(stringr)
library(survival)
library(survminer)
library(survey)
library(srvyr)
library(tableone)
library(questionr)

####DATA EXPLANATION####
#"data" contains cleaned NRD 2017 with the following inclusion criteria:
# adults aged ≥18 years (will subset to older adults later on)
# index hospitalizations of patients with principal diagnosis (ICD_DX1) of multiple rib fractures, and their 3-month readmission encounters


#### DATA CLEANING####
#import dataframe with index rib fx admit encounters + subsequent 3-month readmissions 
data<-fread("Data/dHTX_NRD2017_trauma_only_formatted_meet_criteria.csv") #fread to rapidly import large files
  colnames(data) 
  
data_iss<-icdpicr::cat_trauma(data,"I10_DX", TRUE, "roc_max_NIS") # this calculates ISS, AIS scores using ICDpic
  write.csv(data_iss,"/Data/dHTX_NRD2017_data_iss.csv")


#identify patients with isolated thoracic injuries (AIS non-chest ≤2, AIS-chest >2). Make sure to get ALL encounters 
data_iss<-data_iss %>%
  mutate(index_iso_thoracic_injury=ifelse((RibFx_indexadmission==1 & mxaisbr_General<=2 & mxaisbr_HeadNeck<=2 & mxaisbr_Face<=2 & 
                                       mxaisbr_Extremities<=2 & mxaisbr_Abdomen <=2 & mxaisbr_Chest>2),1,0)) %>%
  mutate(isolated_pleural_effusion=ifelse((Pleural_effusion_new==1 & Hemothorax_new==0 & Hemopneumothorax_new==0 & Hemopneumothorax_new==0),1,0)) %>%
  #consider only isolated pleural effusion! no concomitant hemothorax or hemopneumothorax dx 
  mutate(isolated_htx=ifelse((Hemothorax_new==1 &Pleural_effusion_new==0 &  Hemopneumothorax_new==0&  Pneumothorax_new==0),1,0)) %>%
  mutate(isolated_ptx=ifelse((Hemothorax_new==0 &Pleural_effusion_new==0 &  Hemopneumothorax_new==0&  Pneumothorax_new==1),1,0)) %>%
  mutate(hptx=ifelse((Pleural_effusion_new==0 &  ((Hemothorax_new==1 &Pneumothorax_new==1)| Hemopneumothorax_new==1)),1,0)) %>%
  mutate(pleural_dz=ifelse(isolated_htx==1,"isolated_htx",
                           ifelse(isolated_ptx==1,"isolated_ptx",
                                  ifelse(isolated_pleural_effusion==1,"isolated_pleural_effusion",
                                         ifelse(hptx==1,"hptx",
                                                ifelse((isolated_htx==0 & isolated_ptx==0 & isolated_pleural_effusion==0 & hptx==0),"none","other"))))))
  xtabs(~pleural_dz,data_iss) #no "others
  
#data_iso_thoracic contains index + 3mo readmit encounters of pts who had multiple rib fx and isolated thoracic injuries during index admission 
data_iso_thoracic<-data_iss %>%
  group_by(NRD_VISITLINK) %>%
  filter(any(index_iso_thoracic_injury==1))

#categorize groups by how many 3 month readmissions they had
data_iso_thoracic<-data_iso_thoracic %>%
  group_by(NRD_VISITLINK) %>%
  mutate(readmit_cat=ifelse(total_3mo_readmit==0,0,
                            ifelse(total_3mo_readmit==1,1,
                                   ifelse(total_3mo_readmit==2,2,3))))

  #write.csv(data_iso_thoracic,"Data/data_iso_thoracic_index.csv") 


### FOCUS ON OLDER ADULT SUBSET (Age ≥65)
data_iso_thoracic<-data_iso_thoracic %>% mutate(geriatric=ifelse(AGE_atindex>=65,1,0))

data_iso_thoracic_geriatric<-data_iso_thoracic %>%
  filter(geriatric==1)

#tabulations
  length(unique(data_iso_thoracic_geriatric$NRD_VISITLINK)) #13,978 unique pts (non-survey weighted)
  length(unique(data_iso_thoracic_geriatric[data_iso_thoracic_geriatric$readmitted==1,]$NRD_VISITLINK)) #2763 (19.8%) had at least one 3-month readmission 

#link to CCS codes: FOR FIRST ANALYSIS, KEEP ALL. FOR SENSITIVITY ANALYSIS, LIMIT 
#https://github.com/vcastro/CCS
#https://www.hcup-us.ahrq.gov/toolssoftware/ccs10/CCSCategoryNames(FullLabels).pdf
  CCS_dict<-CCS::CCS_DX_mapping
    CCS_dict<-CCS_dict %>% filter(vocabulary_id=="ICD10CM")
    CCS_define<-CCS::CCS_DX_categories
    CCS_dict<-merge(CCS_dict,CCS_define)

  data_iso_thoracic_geriatric<-data_iso_thoracic_geriatric %>%
    mutate_at(vars(starts_with("I10_DX")),list(~ifelse(.%in% CCS_dict$code, .,""))) %>%
    relocate(c(NRD_VISITLINK,I10_DX1,three_month_readmit_number), .after = last_col()) %>%
    mutate(ribfx_related_readmit=ifelse((three_month_readmit_number>=1 & 
                                          I10_DX1 %in% CCS_dict$code),
                                        1,0))   #should be same here 
    #ribfx_related_readmit will be NA if if three_month_readmit_number is NA and I10_DX1(primary code) is relevant
  
  data_iso_thoracic_geriatric<-data_iso_thoracic_geriatric %>%
    group_by(NRD_VISITLINK) %>%
    mutate(any_ribfx_related_readmit=ifelse(any(ribfx_related_readmit==1),1,0))
  data_iso_thoracic_geriatric$any_ribfx_related_readmit[is.na(data_iso_thoracic_geriatric$any_ribfx_related_readmit)] <- 0
    #turns NA (no readmits at all) to 0
  
  #convert all ICD10 codes to first 4 digits (DID NOT DO)
  #first4<-function(x) {substr(x,start=1,stop=4)} #use only first 4 letters of icd10 codes
  #data_iso_thoracic_geriatric<-data_iso_thoracic_geriatric %>%
    #mutate_at(vars(starts_with("I10_DX")),first4)
  #write.csv(data_iso_thoracic_geriatric,"Data_Aug31/data_iso_thoracic_geriatric_aug31.csv")

  length(unique(data_iso_thoracic_geriatric[data_iso_thoracic_geriatric$any_ribfx_related_readmit==1,]$NRD_VISITLINK)) #2763= all pts with rib-fx related readmissions
  length(unique(data_iso_thoracic_geriatric$NRD_VISITLINK)) #13978 unique patients 

  data_iso_thoracic_geriatric<-data_iso_thoracic_geriatric %>% select(-V1)

  
####DATA ANALYSIS####
#Specify the sampling design with sampling weights DISCWT,
# hospital clusters HOSP_NRD, and stratification NRD_STRATUM
#consider analysis using srvyr: https://cran.r-project.org/web/packages/srvyr/vignettes/srvyr-vs-survey.html
svytable <-data_iso_thoracic_geriatric %>% as_survey_design(HOSP_NRD, weights = DISCWT,
                      strata = NRD_STRATUM)

print(summary(svytable))
  #need to deal with lonely PSU: https://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
options(survey.lonely.psu="adjust")

#OUTCOME 1: compare index hospitalization vs non-index hospitalization LOS and charges  
  test<- data_iso_thoracic_geriatric #%>% replace(is.na(.),0)
  svytable2<-test %>% as_survey_design(HOSP_NRD, weights = DISCWT,
                                                              strata = NRD_STRATUM)
  test<-svyCreateContTable(c("LOS","TOTCHG"),
                     strata="RibFx_indexadmission",
                     data=svytable2)
  print(test,smd=T)

#TABLE 1: readmitted vs non-readmitted patients, index hospitalizations
  CAT_VAR<-c("FEMALE","HOSP_UR_TEACH","Chest_tube","Thoracotomy","VATS","Ventilation",
              "Multiple_RibFx_new","Sternal_Fx_new","Flail_chest_new","Pulmonary_contusion_new",
             "Hemothorax_new","Pleural_effusion_new","Pneumothorax_new","Hemopneumothorax_new",
             "isolated_pleural_effusion","isolated_htx","isolated_ptx","pleural_dz",
             "Pneumonia_new","Respiratory_Failure_new","PYOTHORAX_new","DIED","DISPUNIFORM")
  CONT_VAR<-c("AGE_atindex","niss","LOS","LOS_atindex","TOTCHG")
  
  #to look at individual encounters! 
  tableone_CAT<-svyCreateCatTable(CAT_VAR,
                    strata="any_ribfx_related_readmit",
                    data=svytable%>%filter(RibFx_indexadmission==1))
  
    tableone_CAT_print<-print(tableone_CAT,smd=T,catDigits=1)
    #write.csv(tableone_CAT_print,"Data_Aug31/tableone_CAT_print.csv")
  
  tableone_CONT<-svyCreateContTable(CONT_VAR,
                                  strata="any_ribfx_related_readmit",
                                  data=svytable%>%filter(RibFx_indexadmission==1))
  
    tableone_CONT_print<-print(tableone_CONT,smd=T,contDigits=1)
    #write.csv(tableone_CONT_print,"Data_Aug31/tableone_CONT_print.csv")
    
    #total hospital charge due to readmissions
    svytable%>%filter(three_month_readmit_number==1) %>%
      summarise(total_charge=survey_total(TOTCHG))
    
    #proporiton of patients by number of readmissions 
    svyCreateCatTable("total_3mo_readmit",
                               data=svytable%>%filter(RibFx_indexadmission==1))
        #3818/4984 had one 3-month readmission
    
    
###LEADING PRINCIPAL DIAGNOSES (DON'T CONDENSE ICD10 CODES INTO 4 DIGITS)####
      data_iso_thoracic_geriatric_fullicd<-data_iso_thoracic %>%
        filter(geriatric==1)
      
      data_iso_thoracic_geriatric_fullicd<-data_iso_thoracic_geriatric_fullicd %>%
        mutate_at(vars(starts_with("I10_DX")),list(~ifelse(.%in% CCS_dict$code, .,""))) %>%
        relocate(c(NRD_VISITLINK,I10_DX1,three_month_readmit_number), .after = last_col()) %>%
        mutate(ribfx_related_readmit=ifelse((three_month_readmit_number>=1 & 
                                               I10_DX1 %in% CCS_dict$code),
                                            1,0)) 
      
      data_iso_thoracic_geriatric_fullicd<-data_iso_thoracic_geriatric_fullicd %>%
        group_by(NRD_VISITLINK) %>%
        mutate(any_ribfx_related_readmit=ifelse(any(ribfx_related_readmit==1),1,0))
      
      data_iso_thoracic_geriatric_fullicd$any_ribfx_related_readmit[is.na(data_iso_thoracic_geriatric_fullicd$any_ribfx_related_readmit)] <- 0

      
      svytable_fullicd <-data_iso_thoracic_geriatric_fullicd %>% as_survey_design(HOSP_NRD, weights = DISCWT,
                                                                                  strata = NRD_STRATUM)
      
      readmit_dx_freq_fullicd<-svytable_fullicd %>%
        filter(three_month_readmit_number==1) %>% 
        survey_count(I10_DX1) %>%
        arrange(desc(n)) %>%
        rename(code=I10_DX1)
        
        write.csv(readmit_dx_freq_fullicd,"Data_Aug31/readmit_dx_freq_fullicd.csv")
      
      #Recode potential leading diagnoses (from earlier analytic version)
      respfail_icd10<-readmit_dx_freq_fullicd %>% 
        filter(category_desc=="Respiratory failure; insufficiency; arrest (adult)")
      respfail_icd10<-respfail_icd10$code
      
      hhkd_icd10<-c("I110","I119","I120","I129","I130","I1310","I1311","I132")
      
      htx_icd10<-c("J942","S271XXA","S271XXD","S271XXS","S272XXA","S272XXD","S272XXS")
      
      pleur_eff_icd10<-"J90"
      
      sepsis_icd10<-readmit_dx_freq_fullicd %>% 
        filter(category_desc=="Septicemia (except in labor)")
      sepsis_icd10<-sepsis_icd10$code
      
      pna_icd10<-readmit_dx_freq_fullicd %>% 
        filter(category_desc=="Pneumonia (except that caused by tuberculosis or sexually transmitted disease)")
      pna_icd10<-pna_icd10$code
      
      uti_icd10<-readmit_dx_freq_fullicd %>% 
        filter(category_desc=="Urinary tract infections")
      uti_icd10<-uti_icd10$code
      
      
      #Recode ICD10 codes with common principal diagnoses (know from previous analysis version)
      data_iso_thoracic_geriatric_fullicd<-data_iso_thoracic_geriatric_fullicd %>% 
        select(-V1) %>%
        mutate(across(starts_with("I10_DX"),~replace(.,. %in% sepsis_icd10,"Sepsis"))) %>%
        mutate(across(starts_with("I10_DX"),~replace(.,. %in% hhkd_icd10,"Hhkd"))) %>%
        mutate(across(starts_with("I10_DX"),~replace(.,. %in% htx_icd10,"Hemothorax"))) %>%
        mutate(across(starts_with("I10_DX"),~replace(.,. %in% pna_icd10,"Pneumonia"))) %>%
        mutate(across(starts_with("I10_DX"),~replace(.,. %in% pleur_eff_icd10,"Pleural_effusion"))) %>%
        mutate(across(starts_with("I10_DX"),~replace(.,. %in% respfail_icd10,"Respiratory_failure")))
        

      
#####RE_ADMISSION REASONS among home vs facility dispo####
        #home dispo
         home_readmit_fullicd<-svytable_fullicd %>%
          filter(three_month_readmit_number==1 & (NRD_VISITLINK %in% id_home_index_dispo$NRD_VISITLINK)) %>% #??can use is.na(RibFx_indexadmission) instead of ribfx_related_readmit==1 to get total N!
          survey_count(I10_DX1) %>%
          arrange(desc(n))%>%
           rename(code=I10_DX1)
      
         write.csv(home_readmit_fullicd,"Data_Aug31/home_readmit_fullicd.csv")
        
        #facility dispo
        facility_readmit_fullicd<-svytable_fullicd %>%
          filter(three_month_readmit_number==1 & (NRD_VISITLINK %in% id_facility_index_dispo$NRD_VISITLINK)) %>%
          survey_count(I10_DX1) %>%
          arrange(desc(n))%>%
          rename(code=I10_DX1)

         write.csv(facility_readmit_fullicd,"Data_Aug31/facility_readmit_fullicd.csv")
        

####AMONG PATIENTS READMITTED WITH SEPSIS, WHAT ARE THE MOST COMMON CONCURRENT DIAGNOSES?####
          #Narrow down to infectious diagnosis CCS categories
          CCS_infection<-c("2","3","4","5","7","8","9","76","77", "78","90","122","123","124","125",
                           "126", "135","146","147","148","149","151", "159","197") #included all infections
          CCS_dict_infection<-CCS_dict %>%filter(category_code %in% CCS_infection)
          
          ##from prior analytic file know likely top two are pna and uti
              #combine and evaluate secondary diagnosis: I10_DX2 to I10_DX40
              sepsis_cohort_fullicd<-data_iso_thoracic_geriatric_fullicd %>%
                filter(three_month_readmit_number==1 & I10_DX1 %in% sepsis_icd10) %>%
                select(starts_with("I10_DX"),NRD_VISITLINK, HOSP_NRD,DISCWT,NRD_STRATUM) %>%
                select(-I10_DX1) %>% #get rid of principle diagnosis column 
                mutate(secondary_pna=ifelse((I10_DX2 %in% pna_icd10|I10_DX3 %in% pna_icd10|I10_DX4 %in% pna_icd10|
                                                I10_DX5 %in% pna_icd10|I10_DX6 %in% pna_icd10|I10_DX7 %in% pna_icd10|
                                                I10_DX8 %in% pna_icd10|I10_DX9 %in% pna_icd10|I10_DX10 %in% pna_icd10|
                                                I10_DX11 %in% pna_icd10|I10_DX12 %in% pna_icd10|I10_DX13 %in% pna_icd10|
                                                I10_DX14 %in% pna_icd10|I10_DX15 %in% pna_icd10|I10_DX16 %in% pna_icd10|
                                                I10_DX17 %in% pna_icd10|I10_DX18 %in% pna_icd10|I10_DX19 %in% pna_icd10|
                                                I10_DX20 %in% pna_icd10|I10_DX21 %in% pna_icd10|I10_DX22 %in% pna_icd10|
                                                I10_DX23 %in% pna_icd10|I10_DX24 %in% pna_icd10|I10_DX25 %in% pna_icd10|
                                                I10_DX26 %in% pna_icd10|I10_DX27 %in% pna_icd10|I10_DX28 %in% pna_icd10|
                                                I10_DX29 %in% pna_icd10|I10_DX30 %in% pna_icd10|I10_DX31 %in% pna_icd10|
                                                I10_DX32 %in% pna_icd10|I10_DX33 %in% pna_icd10|I10_DX34 %in% pna_icd10|
                                                I10_DX35 %in% pna_icd10|I10_DX36 %in% pna_icd10|I10_DX37 %in% pna_icd10|
                                                I10_DX38 %in% pna_icd10|I10_DX39 %in% pna_icd10|I10_DX40 %in% pna_icd10),
                                             1,0),
                          secondary_uti=ifelse((I10_DX2 %in% uti_icd10|I10_DX3 %in% uti_icd10|I10_DX4 %in% uti_icd10|
                                                I10_DX5 %in% uti_icd10|I10_DX6 %in% uti_icd10|I10_DX7 %in% uti_icd10|
                                                I10_DX8 %in% uti_icd10|I10_DX9 %in% uti_icd10|I10_DX10 %in% uti_icd10|
                                                I10_DX11 %in% uti_icd10|I10_DX12 %in% uti_icd10|I10_DX13 %in% uti_icd10|
                                                I10_DX14 %in% uti_icd10|I10_DX15 %in% uti_icd10|I10_DX16 %in% uti_icd10|
                                                I10_DX17 %in% uti_icd10|I10_DX18 %in% uti_icd10|I10_DX19 %in% uti_icd10|
                                                I10_DX20 %in% uti_icd10|I10_DX21 %in% uti_icd10|I10_DX22 %in% uti_icd10|
                                                I10_DX23 %in% uti_icd10|I10_DX24 %in% uti_icd10|I10_DX25 %in% uti_icd10|
                                                I10_DX26 %in% uti_icd10|I10_DX27 %in% uti_icd10|I10_DX28 %in% uti_icd10|
                                                I10_DX29 %in% uti_icd10|I10_DX30 %in% uti_icd10|I10_DX31 %in% uti_icd10|
                                                I10_DX32 %in% uti_icd10|I10_DX33 %in% uti_icd10|I10_DX34 %in% uti_icd10|
                                                I10_DX35 %in% uti_icd10|I10_DX36 %in% uti_icd10|I10_DX37 %in% uti_icd10|
                                                I10_DX38 %in% uti_icd10|I10_DX39 %in% uti_icd10|I10_DX40 %in% uti_icd10),
                                             1,0))
                
          sepsis_cohort_svy_fullicd <-sepsis_cohort_fullicd %>%
            as_survey_design(HOSP_NRD, weights = DISCWT,strata = NRD_STRATUM)
          
          #Rate of secondary PNA: 235.1/579
          svyCreateCatTable("secondary_uti", #placeholder variable 
                                  strata="secondary_pna",
                                  data=sepsis_cohort_svy_fullicd)
          
          #Rate of secondary UTI: 235.1/579
          svyCreateCatTable("secondary_pna",#placeholder variable 
                            strata="secondary_uti",
                            data=sepsis_cohort_svy_fullicd)
          
          
####PER READMISSION REASON,"days2index", "LOS_atindex", "LOS", "TOTCHG", DISPUNIFORM ####
          readmit_specifics_dec<-svytable_fullicd %>%
            filter(three_month_readmit_number==1) %>%
            filter(I10_DX1=="Sepsis"|I10_DX1=="Hhkd"|I10_DX1=="Hemothorax"|I10_DX1=="Pneumonia"|I10_DX1=="Respiratory_failure")%>%
            group_by(I10_DX1) %>% 
            summarise(
              index_los=survey_quantile(LOS_atindex,c(0.25, 0.5, 0.75)),
              days_to_readmit=survey_quantile(days2index,c(0.25, 0.5, 0.75)),
              readmit_los=survey_quantile(LOS,c(0.25, 0.5, 0.75)),
              readmit_charge=survey_quantile(TOTCHG,c(0.25, 0.5, 0.75)))
            
            write.csv(readmit_specifics_dec,"Data_Aug31/readmit_specifics_dec.csv")
      
          #Optional figure showing boxplot
          questionr::ggsurvey(svytable_fullicd %>%
                                filter(three_month_readmit_number==1 & (I10_DX1=="Sepsis"|I10_DX1=="Hhkd"|I10_DX1=="Hemothorax"|I10_DX1=="Pneumonia"|I10_DX1=="Respiratory_failure")),
                              aes(x=days2index,fill=I10_DX1)) +
            geom_boxplot()+
            scale_fill_brewer(palette="Set3")+
            facet_wrap(~I10_DX1,ncol=1)+
            theme_classic2()+
            theme(strip.text=element_blank())+
            scale_x_continuous(name="Days to Readmission",limits=c(0,90),breaks=seq(0,90,10))+
            labs(fill="Principal readmission diagnosis")
          
####MANAGEMENT OF PLEURAL DISEASES####
  #need to compare those who were not readmitted vs those readmitted without htx/Pleureff as primary diagnosis vs those reamditted WITH htx/pleureff as primary diagnosis 
  pleural_dz_mgmt<-data_iso_thoracic_geriatric %>%
    mutate(hptx_mgmt=ifelse(Hemopneumothorax_new==1 & (VATS==1 | Chest_tube==1),"hptx_intervene",
                            ifelse(Hemopneumothorax_new==1 & (VATS==0 & Chest_tube==0),"hptx_no_intervene",NA)),
           pleureff_mgmt=ifelse(Pleural_effusion_new==1 & (VATS==1 | Chest_tube==1),"pleureff_intervene",
                                            ifelse(Pleural_effusion_new==1 & (VATS==0 & Chest_tube==0),"pleureff_no_intervene",NA)),
           hthx_mgmt=ifelse(Hemothorax_new==1 & (VATS==1 | Chest_tube==1),"htx_intervene",
                                ifelse(Hemothorax_new==1 & (VATS==0 & Chest_tube==0),"htx_no_intervene",NA))) %>%
        filter(three_month_readmit_number<=1)
 
    #ID of patients re-admitted with htx or pleural effusion as primary diagnosis 
    pleural_dz_primary_readmit<-data_iso_thoracic_geriatric %>% 
      filter(three_month_readmit_number==1 & (I10_DX1=="S271"|I10_DX1=="S272"|I10_DX1=="J90"|I10_DX1=="J942")) %>%
      select(NRD_VISITLINK)
    
    pleural_eff_primary_readmit<-data_iso_thoracic_geriatric %>% 
      filter(three_month_readmit_number==1 & I10_DX1=="J90") %>%
      select(NRD_VISITLINK)
    
    htx_primary_readmit<-data_iso_thoracic_geriatric %>% 
      filter(three_month_readmit_number==1 & (I10_DX1=="S271"|I10_DX1=="S272"|I10_DX1=="J942")) %>%
      select(NRD_VISITLINK)
    
    pleural_dz_mgmt<-pleural_dz_mgmt %>%
      relocate(RibFx_indexadmission, .after = last_col())
    
    pleural_dz_mgmt$RibFx_indexadmission[is.na(pleural_dz_mgmt$RibFx_indexadmission)]<-0 #really important to do! 
   
    #pleural_dz_class: the 3 categories 
    pleural_dz_mgmt<-pleural_dz_mgmt %>%
        mutate(pleural_dz_class=ifelse(RibFx_indexadmission==1,"index_admission",
                                       ifelse(three_month_readmit_number==1 & (I10_DX1=="S271"|I10_DX1=="S272"|I10_DX1=="J90"|I10_DX1=="J942"),"pleural_dz_readmit",NA)))
         
    #make survey style format 
    svytable_pleural_dz<-pleural_dz_mgmt %>% as_survey_design(HOSP_NRD, weights = DISCWT,
                                       strata = NRD_STRATUM)
  
    #QUESTION 1: among those with primary pleural dz readmit, how many underwent intervention at index vs not index
      #pleural_dz_primary_readmit= ID of patients readmitted with pleural dz
      #this compares index_admission vs pleural_dz_readmit. index_admission may exclude pts who didn't have these diagnoses
    pleural_interventions<-svyCreateCatTable(c("Chest_tube","VATS"),
                      strata="pleural_dz_class",
                      includeNA=FALSE,
                      data=svytable_pleural_dz %>% filter(NRD_VISITLINK %in% pleural_dz_primary_readmit$NRD_VISITLINK)) 
    pleural_interventions<-print(pleural_interventions,smd=T)
    
      write.csv(pleural_interventions,"Data_Aug31/pleural_interventions_aug31.csv") #more interventions upon readmission 

    #Question 2: among those readmitted with hemothorax, what diagnoses dit they at index and what interventions
    svytable_htx<-svyCreateCatTable(c("isolated_pleural_effusion","isolated_htx", "Hemothorax_new","Pleural_effusion_new","hthx_mgmt","pleureff_mgmt"),
                                    strata="pleural_dz_class",
                                    includeNA=FALSE,
                                    data=svytable_pleural_dz %>% filter(NRD_VISITLINK %in% htx_primary_readmit$NRD_VISITLINK)) #those who ahd pleaural dz readmissions
    svytable_htx<-print(svytable_htx,smd=T)
      
      write.csv(svytable_htx,"Data_Aug31/svytable_htx_aug31.csv")
    
    #Question 3: among those readmitted with pleural effusion, what diagnoses did they have at index and what interventions?
    svytable_pleuraleff<-svyCreateCatTable(c("isolated_pleural_effusion","isolated_htx", "Hemothorax_new","Pleural_effusion_new","hthx_mgmt","pleureff_mgmt"),
                                          strata="pleural_dz_class",
                                          includeNA=FALSE,
                                          data=svytable_pleural_dz %>% filter(NRD_VISITLINK %in% pleural_eff_primary_readmit$NRD_VISITLINK)) #those who ahd pleaural dz readmissions
    svytable_pleuraleff<-print(svytable_pleuraleff,smd=T)
      
      write.csv(svytable_pleuraleff,"Data_Aug31/svytable_pleuraleff_aug31.csv")
     

####SENSITIVITY ANALYSIS#### 
    #NUMBERS LIMITING READMIT CODES TO CERTAIN CCS CATEGORIES 
    #identify readmissions related to rib fracture-specific sequelae (refer to NRD_CCS_categories.doc)
    #CCS_categories to keep according to https://www-ncbi-nlm-nih-gov.laneproxy.stanford.edu/pmc/articles/PMC6345420/#pone.0209896.ref019
    CCS_keep_sensitivity<-c("2","3","60","62","122","126","129","130","131","133","234")
    CCS_dict_sensitivity<-CCS_dict %>%filter(category_code %in% CCS_keep_sensitivity)
    
    sn_data_iso_thoracic_geriatric<-data_iso_thoracic %>%
      filter(geriatric==1)
    
    sn_data_iso_thoracic_geriatric<-sn_data_iso_thoracic_geriatric %>%
      mutate_at(vars(starts_with("I10_DX")),list(~ifelse(.%in% CCS_dict_sensitivity$code, .,""))) %>%
      mutate(ribfx_related_readmit=ifelse((three_month_readmit_number>=1 & #this dataframe has different definition of "rib fracture related admit"
                                             I10_DX1 %in% CCS_dict_sensitivity$code),
                                          1,0)) %>%
      filter(three_month_readmit_number<=1)
    
    sn_data_iso_thoracic_geriatric<-sn_data_iso_thoracic_geriatric %>%
      group_by(NRD_VISITLINK) %>%
      mutate(any_ribfx_related_readmit=ifelse(any(ribfx_related_readmit==1),1,0))
    sn_data_iso_thoracic_geriatric$any_ribfx_related_readmit[is.na(sn_data_iso_thoracic_geriatric$any_ribfx_related_readmit)] <- 0
    
    #Recode variables. 
    sn_data_iso_thoracic_geriatric<-sn_data_iso_thoracic_geriatric %>%
      mutate(across(starts_with("I10_DX"),~replace(.,. %in% sepsis_icd10,"Sepsis"))) %>%
      mutate(across(starts_with("I10_DX"),~replace(.,. %in% hhkd_icd10,"Hhkd"))) %>%
      mutate(across(starts_with("I10_DX"),~replace(.,. %in% htx_icd10,"Hemothorax"))) %>%
      mutate(across(starts_with("I10_DX"),~replace(.,. %in% pna_icd10,"Pneumonia"))) %>%
      mutate(across(starts_with("I10_DX"),~replace(.,. %in% pleur_eff_icd10,"Pleural_effusion"))) %>%
      mutate(across(starts_with("I10_DX"),~replace(.,. %in% respfail_icd10,"Respiratory_failure"))) 
    
    #weighted analysis starts here 
    sn_svytable <-sn_data_iso_thoracic_geriatric %>% as_survey_design(HOSP_NRD, weights = DISCWT,
                                                                strata = NRD_STRATUM)
    
    sn_tableone_CAT<-svyCreateCatTable(CAT_VAR,
                                    strata="any_ribfx_related_readmit",
                                    data=sn_svytable%>%filter(RibFx_indexadmission==1))
    sn_tableone_CAT_print<-print(sn_tableone_CAT,smd=T,catDigits=1)
      #1453 of 25,092 had readmission with a principal dx really attributable to rib fx 
    
    #write.csv(sn_tableone_CAT_print,"Data_Aug31/sn_tableone_CAT_print.csv")
    
    sn_tableone_CONT<-svyCreateContTable(CONT_VAR,
                                      strata="any_ribfx_related_readmit",
                                      data=sn_svytable%>%filter(RibFx_indexadmission==1))
    sn_tableone_CONT_print<-print(sn_tableone_CONT,smd=T,contDigits=1)
    
    #write.csv(sn_tableone_CONT_print,"Data_Aug31/sn_tableone_CONT_print.csv")
  
    #frequency table 
    sn_readmit_dx_freq<-sn_svytable %>%
      filter(three_month_readmit_number==1) %>% 
      survey_count(I10_DX1) %>%
      arrange(desc(n))
    #write.csv(sn_readmit_dx_freq,"Data_Aug31/sn_readmit_dx_freq.csv")
    #leading= sepsis, htx, pna, pleural effusion, respiratory failure
  
  #time to readmission: optional figure
  questionr::ggsurvey(sn_svytable %>%
                        filter(three_month_readmit_number==1 & I10_DX1 %in% c("A419","J189","S271","J90")) %>%
                        mutate(I10_DX1=recode(I10_DX1,A419="Sepsis",J189="Pneumonia",S271="Hemothorax",J90="Pleural effusion")),
                      aes(x=days2index,fill=I10_DX1)) +
    geom_boxplot()+
    scale_fill_brewer(palette="Set3")+
    facet_wrap(~I10_DX1,ncol=1)+
    theme_classic2()+
    theme(strip.text=element_blank())+
    scale_x_continuous(name="Days to Readmission",limits=c(0,90),breaks=seq(0,90,10))+
    labs(fill="Principal readmission diagnosis")
  
  