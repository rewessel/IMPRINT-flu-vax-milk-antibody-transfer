setwd('Z:/users/wesselr/data/IMPRINT')
obcrf = read_csv('./unedited_datasets/FW_ Nasal Swab Data/Remziye.Obstetrical.Form_2025-09-05.csv')
nbcrf = read_csv('./unedited_datasets/FW_ Nasal Swab Data/Remziye.Newborn.Form_2025-09-05.csv')
bl = read_csv('./unedited_datasets/FW_ Nasal Swab Data/Remziye.Baseline.Form_2025-09-05.csv')


# Extract continuous variables from each data set
obcrf.cols.cont = c('Subject.ID','mat_deliv_wt_obcrf','gravida_obcrf')
nbcrf.cols.cont = c('Subject.ID','birthwt_nbcrf','gestage_delivery_nbcrf')
bl.cols.cont = c('Subject.ID','income_bl','prepreg_weight_bl','mat_totalheight_bl','second_hand_smoke_bl','feed_trybf_bl')

# Convert to numeric and calculate BMI and GWG columns
cont.var = cbind(obcrf[,obcrf.cols.cont],nbcrf[,nbcrf.cols.cont],bl[,bl.cols.cont]) 
cont.var = cont.var[,-c(1,4,7)] %>% mutate_all(as.numeric)

cont.var$gest_wt_gain = cont.var$mat_deliv_wt_obcrf - cont.var$prepreg_weight_bl
cont.var$prepreg_bmi = cont.var$prepreg_weight_bl*0.4536/(cont.var$mat_totalheight_bl*0.0254)^2

# Convert gestational age into weeks (numeric)
gestage.weeks = as.numeric(substr(cont.var$gestage_delivery_nbcrf, 1, 2))
gestage.days = as.numeric(substr(cont.var$gestage_delivery_nbcrf, 3, 4))
cont.var$gestage_delivery_nbcrf = gestage.weeks+gestage.days/7

# Remove height and weight columns, clean up missing data
cont.var = cont.var[,-c(1,6,7)]
cont.var$income_bl[cont.var$income_bl==77] = 7
cont.var[cont.var==88 | cont.var==99] = NA

colnames(cont.var) = c('gravida','birthweight_lb','gestational_age','income','second_hand_smoke','breastfeeding',
                       'gestational_weight_gain','prepreg_BMI')
rownames(cont.var) = obcrf$Subject.ID

# Extract categorical variables from each data set
obcrf.cols.cat = c('Subject.ID','deliverymode_obcrf','gbs_urine_obcrf','gbs_rectovaginal_obcrf','htn_chronic_obcrf','htn_gestational_obcrf',
                   'preeclampsia_obcrf','dm_obcrf','hellp_obcrf','ptl_obcrf','iugr_obcrf','dm_gestational_obcrf','prom_obcrf',
                   'antenatal_hosp_obcrf')
nbcrf.cols.cat = c('Subject.ID','inf_gender_nbcrf','nicu_nbcrf','hosp_feeding_nbcrf')
bl.cols.cat = c('Subject.ID','insurance_bl','insurance_type_bl','ethnicity_bl','race_bl___1','marital_status_bl','education_bl','flu_vax_preg_bl',
                'mat_ari_bl___0','swab_test_pos_bl','vitamins_bl___1','travel_bl')

cat.var = cbind(obcrf[,obcrf.cols.cat],nbcrf[,nbcrf.cols.cat],bl[,bl.cols.cat])
cat.var[cat.var==88 | cat.var==99] = NA
cat.var = cat.var[,-c(1,15,19)] %>% mutate_all(as.factor)

cat.var$GBS[cat.var$gbs_rectovaginal_obcrf==1 | cat.var$gbs_urine_obcrf==1] = 1 
cat.var$GBS[is.na(cat.var$GBS)] = 0
cat.var = cat.var[,-c(2,3)]

colnames(cat.var) = c('delivery_mode','chronic_hypertension','gestational_hypertension','preeclampsia','diabetes','HELLP',
                      'preterm_labor','IUGR','gestational_diabetes','PROM','antenatal_hosp','child_sex','NICU','hospital_breastfeed',
                      'insurance','insurance_type','ethnicity','race_white','marital_status','education_level','flu_vax_in_preg',
                      'mat_symptoms','positive_covid_test','prenatal_vitamins','international_travel','GBS')
rownames(cat.var) = obcrf$Subject.ID


# Set up birth month timing information
birthdate = matrix(NA,nrow = nrow(nbcrf),ncol=0) %>% as.data.frame()
birthdate$birth.month = sub("\\/.*", "",nbcrf$date_nbcrf) %>% as.numeric() # extract months
birthdate$birth.year = sub(".*/20", "", nbcrf$date_nbcrf) %>% as.numeric() # extract years
birthdate = data.frame(nbcrf$Subject.ID,birthdate)
birthdate$birth.season[birthdate$birth.month==12|birthdate$birth.month==1|birthdate$birth.month==2]='winter'
birthdate$birth.season[birthdate$birth.month==3|birthdate$birth.month==4|birthdate$birth.month==5]='spring'
birthdate$birth.season[birthdate$birth.month==6|birthdate$birth.month==7|birthdate$birth.month==8]='summer'
birthdate$birth.season[birthdate$birth.month==9|birthdate$birth.month==10|birthdate$birth.month==11]='fall'
rownames(birthdate) = birthdate$nbcrf.Subject.ID
birthdate = birthdate[,-1]

cat.var = cbind(cat.var,birthdate)

