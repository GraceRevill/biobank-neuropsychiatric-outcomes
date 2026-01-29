################################################################
#
# Analysis of Biobank dataset
#
# Revill et al
#
################################################################

# Load libraries
library(readr)
library(foreign)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(patchwork)
library(lme4)
library(stringr)
library(reshape2)
library(visdat)
library(jtools)
library(gt)
library(gtsummary)
library(table1)
library(afex)
library(sjPlot)
library(lubridate)



# Remove everything from memory if needed to start afresh
#rm(list = ls())


# Set data and output directories
data_dir = "S:/Track_TBI/Grace/Test/"
output_dir = "S:/Track_TBI/Grace/Test/Output/"

setwd(data_dir)

###########################################################
#
# Load data
#
###########################################################

# Whole cohort
full_cohort <- read.csv("S:/Track_TBI/Grace/Test/whole_cohort2.csv")


# TBI ICD
tbi_icd <- read.csv("S:/Track_TBI/Grace/Test/tbi_icd.csv")


# TBI self-report
tbi_self <- read.csv("S:/Track_TBI/Grace/Test/tbi_self.csv")



# For ICD dates
icd <- read.csv("S:/Track_TBI/Grace/Test/tbi_multiple_icd.csv")
icd_date <- read.csv("S:/Track_TBI/Grace/Test/Participant_table.csv")


###########################################################
#
# Rename and recode variables
#
###########################################################

# Rename and select demographic variables 
full_cohort <- full_cohort %>%
  rename(subjectkey = Participant.ID) %>%
  rename(ethnicity = Ethnic.background...Instance.0) %>%
  rename(age = Age.at.recruitment) %>%
  rename(sex = Sex) %>%
  rename(deprivation_index = Townsend.deprivation.index.at.recruitment) %>%
  rename(alcohol_addict = Ever.addicted.to.alcohol) %>%
  rename(drug_addict = Ever.addicted.to.illicit.or.recreational.drugs) %>%
  rename(alcohol_freq = Alcohol.intake.frequency....Instance.0) %>%
  rename(dementia_source = Source.of.all.cause.dementia.report) %>%
  rename(parkinsons_source = Source.of.all.cause.parkinsonism.report) %>%
  rename(stroke_source = Source.of.stroke.report) %>%
  rename(dementia_date = Date.of.all.cause.dementia.report) %>%
  rename(parkinsons_date = Date.of.all.cause.parkinsonism.report) %>%
  rename(stroke_date = Date.of.stroke) %>%
  rename(dep_ep_source = Source.of.report.of.F32..depressive.episode.) %>%
  rename(dep_disorder_source = Source.of.report.of.F33..recurrent.depressive.disorder.) %>%
  rename(affect_disorder_source = Source.of.report.of.F34..persistent.mood..affective..disorders.) %>%  
  rename(other_mood_affect_source = Source.of.report.of.F38..other.mood..affective..disorders.) %>%  
  rename(phobic_disorder_source = Source.of.report.of.F40..phobic.anxiety.disorders.) %>%    
  rename(anxiety_disorder_source = Source.of.report.of.F41..other.anxiety.disorders.) %>%  
  rename(schizophrenia_source = Source.of.report.of.F20..schizophrenia.) %>%  
  rename(delusion_disorder_source = Source.of.report.of.F22..persistent.delusional.disorders.) %>%  
  rename(acute_psychotic_source = Source.of.report.of.F23..acute.and.transient.psychotic.disorders.) %>%  
  rename(schizoaffective_source = Source.of.report.of.F25..schizoaffective.disorders.) %>%  
  rename(bipolar_source = Source.of.report.of.F31..bipolar.affective.disorder.) %>%  
  rename(organic_mental_source = Source.of.report.of.F06..other.mental.disorders.due.to.brain.damage.and.dysfunction.and.to.physical.disease.) %>%  
  rename(organic_personality_source = Source.of.report.of.F07..personality.and.behavioural.disorders.due.to.brain.disease..damage.and.dysfunction.) %>% 
  rename(ocd_source = Source.of.report.of.F42..obsessive.compulsive.disorder.) %>% 
  rename(dep_ep_date = Date.F32.first.reported..depressive.episode.) %>% 
  rename(dep_disorder_date = Date.F33.first.reported..recurrent.depressive.disorder.) %>% 
  rename(affect_disorder_date = Date.F34.first.reported..persistent.mood..affective..disorders.) %>%  
  rename(other_mood_affect_date = Date.F38.first.reported..other.mood..affective..disorders.) %>% 
  rename(phobic_disorder_date = Date.F40.first.reported..phobic.anxiety.disorders.) %>%   
  rename(anxiety_disorder_date = Date.F41.first.reported..other.anxiety.disorders.) %>%  
  rename(schizophrenia_date = Date.F20.first.reported..schizophrenia.) %>%  
  rename(delusion_disorder_date = Date.F22.first.reported..persistent.delusional.disorders.) %>%  
  rename(acute_psychotic_date = Date.F23.first.reported..acute.and.transient.psychotic.disorders.) %>%   
  rename(schizoaffective_date = Date.F25.first.reported..schizoaffective.disorders.) %>%   
  rename(bipolar_date = Date.F31.first.reported..bipolar.affective.disorder.) %>% 
  rename(organic_mental_date = Date.F06.first.reported..other.mental.disorders.due.to.brain.damage.and.dysfunction.and.to.physical.disease.) %>% 
  rename(organic_personality_date = Date.F07.first.reported..personality.and.behavioural.disorders.due.to.brain.disease..damage.and.dysfunction.) %>% 
  rename(ocd_date = Date.F42.first.reported..obsessive.compulsive.disorder.) %>%
  rename(dep_sad_ever_1 = Ever.had.prolonged.feelings.of.sadness.or.depression.participant...p20446.) %>%
  rename(dep_sad_ever_2 = Ever.had.prolonged.feelings.of.sadness.or.depression.participant...p29011.) %>%
  rename(mod_dep_prob = Probable.recurrent.major.depression..moderate....Instance.0) %>%
  rename(sev_dep_prob = Probable.recurrent.major.depression..severe....Instance.0) %>%
  rename(site = UK.Biobank.assessment.centre...Instance.0)


# Make sure all empty rows are labelled NA 
full_cohort[full_cohort == ""] <- NA
full_cohort[full_cohort == -818] <- NA  # Columns coded as 'Prefer not to say'
full_cohort[full_cohort == -3] <- NA # Columns coded as 'Prefer not to say'

# Recode ethnicity
full_cohort$ethnicity[full_cohort$ethnicity == 1] <- 'White'
full_cohort$ethnicity[full_cohort$ethnicity == 1001] <- 'White'
full_cohort$ethnicity[full_cohort$ethnicity == 2001] <- 'Mixed or multiple ethnic groups'
full_cohort$ethnicity[full_cohort$ethnicity == 3001] <- 'Asian or Asian British'
full_cohort$ethnicity[full_cohort$ethnicity == 4001] <- 'Black, Black British, Caribbean or African'
full_cohort$ethnicity[full_cohort$ethnicity == 2] <- 'Mixed or multiple ethnic groups'
full_cohort$ethnicity[full_cohort$ethnicity == 1002] <- 'White'
full_cohort$ethnicity[full_cohort$ethnicity == 2002] <- 'Mixed or multiple ethnic groups'
full_cohort$ethnicity[full_cohort$ethnicity == 3002] <- 'Asian or Asian British'
full_cohort$ethnicity[full_cohort$ethnicity == 4002] <- 'Black, Black British, Caribbean or African'
full_cohort$ethnicity[full_cohort$ethnicity == 3] <- 'Asian or Asian British'
full_cohort$ethnicity[full_cohort$ethnicity == 1003] <- 'White'
full_cohort$ethnicity[full_cohort$ethnicity == 2003] <- 'Mixed or multiple ethnic groups'
full_cohort$ethnicity[full_cohort$ethnicity == 3003] <- 'Asian or Asian British'
full_cohort$ethnicity[full_cohort$ethnicity == 4003] <- 'Black, Black British, Caribbean or African'
full_cohort$ethnicity[full_cohort$ethnicity == 4] <- 'Black, Black British, Caribbean or African'
full_cohort$ethnicity[full_cohort$ethnicity == 2004] <- 'Mixed or multiple ethnic groups'
full_cohort$ethnicity[full_cohort$ethnicity == 3004] <- 'Asian or Asian British'
full_cohort$ethnicity[full_cohort$ethnicity == 5] <- 'Asian or Asian British'
full_cohort$ethnicity[full_cohort$ethnicity == 6] <- 'Other'
full_cohort$ethnicity[full_cohort$ethnicity == -1] <- NA

# Set data type and level labels
full_cohort$ethnicity <- factor(full_cohort$ethnicity, 
                                levels = c("White", "Black, Black British, Caribbean or African", "Asian or Asian British", "Mixed or multiple ethnic groups", 
                                           "Other"))

# Recode sex 
full_cohort$sex[full_cohort$sex == 0] <- 'Female'
full_cohort$sex[full_cohort$sex == 1] <- 'Male'


# Set data type and level labels
full_cohort$sex <- factor(full_cohort$sex, 
                          levels = c("Female", "Male"))




# Ensure study site ID is stored as a factor  
class(full_cohort$site)
full_cohort$site <- as.factor(full_cohort$site)

# Recode addiction to alcohol
full_cohort$alcohol_addict[full_cohort$alcohol_addict == -121] <- NA
full_cohort$alcohol_addict[full_cohort$alcohol_addict == 0] <- 'No'
full_cohort$alcohol_addict[full_cohort$alcohol_addict == 1] <- 'Yes'



# Set data type and level labels
full_cohort$alcohol_addict <- factor(full_cohort$alcohol_addict, 
                                     levels = c("Yes", "No"))

# Recode addiction to drugs
full_cohort$drug_addict[full_cohort$drug_addict == -121] <- NA
full_cohort$drug_addict[full_cohort$drug_addict == 0] <- 'No'
full_cohort$drug_addict[full_cohort$drug_addict == 1] <- 'Yes'


# Set data type and level labels
full_cohort$drug_addict <- factor(full_cohort$drug_addict, 
                                  levels = c("Yes", "No"))


# Recode alcohol intake frequency 
full_cohort$alcohol_freq[full_cohort$alcohol_freq == 1] <- 'Daily/almost daily'
full_cohort$alcohol_freq[full_cohort$alcohol_freq == 2] <- 'Three/four times a week'
full_cohort$alcohol_freq[full_cohort$alcohol_freq == 3] <- 'Once/twice a week'
full_cohort$alcohol_freq[full_cohort$alcohol_freq == 4] <- 'One to three times a month'
full_cohort$alcohol_freq[full_cohort$alcohol_freq == 5] <- 'Special occasions only'
full_cohort$alcohol_freq[full_cohort$alcohol_freq == 6] <- 'Never'


# Set data type and level labels
full_cohort$alcohol_freq <- factor(full_cohort$alcohol_freq, 
                                   levels = c("Daily/almost daily", "Three/four times a week", "Once/twice a week", "One to three times a month", 
                                              "Special occasions only" ,"Never"))


# Recode source of dementia 
full_cohort$dementia_source[full_cohort$dementia_source == 0] <- 'Self-report only'
full_cohort$dementia_source[full_cohort$dementia_source == 1] <- 'Hospital admission'
full_cohort$dementia_source[full_cohort$dementia_source == 2] <- 'Death only'
full_cohort$dementia_source[full_cohort$dementia_source == 11] <- 'Hospital primary'
full_cohort$dementia_source[full_cohort$dementia_source == 12] <- 'Death primary'
full_cohort$dementia_source[full_cohort$dementia_source == 21] <- 'Hospital secondary'
full_cohort$dementia_source[full_cohort$dementia_source == 22] <- 'Death contributory'


# Recode source of parkinsons
full_cohort$parkinsons_source[full_cohort$parkinsons_source == 0] <- 'Self-report only'
full_cohort$parkinsons_source[full_cohort$parkinsons_source == 1] <- 'Hospital admission'
full_cohort$parkinsons_source[full_cohort$parkinsons_source == 2] <- 'Death only'
full_cohort$parkinsons_source[full_cohort$parkinsons_source == 11] <- 'Hospital primary'
full_cohort$parkinsons_source[full_cohort$parkinsons_source == 12] <- 'Death primary'
full_cohort$parkinsons_source[full_cohort$parkinsons_source == 21] <- 'Hospital secondary'
full_cohort$parkinsons_source[full_cohort$parkinsons_source == 22] <- 'Death contributory'


# Recode source of stroke
full_cohort$stroke_source[full_cohort$stroke_source == 0] <- 'Self-report only'
full_cohort$stroke_source[full_cohort$stroke_source == 1] <- 'Hospital admission'
full_cohort$stroke_source[full_cohort$stroke_source == 2] <- 'Death only'
full_cohort$stroke_source[full_cohort$stroke_source == 11] <- 'Hospital primary'
full_cohort$stroke_source[full_cohort$stroke_source == 12] <- 'Death primary'
full_cohort$stroke_source[full_cohort$stroke_source == 21] <- 'Hospital secondary'
full_cohort$stroke_source[full_cohort$stroke_source == 22] <- 'Death contributory'


# Create yes/no variable for neurological disorders 
full_cohort$dementia <- ifelse(is.na(full_cohort$dementia_date), "No", "Yes")
full_cohort$parkinsons <- ifelse(is.na(full_cohort$parkinsons_date), "No", "Yes")
full_cohort$stroke <- ifelse(is.na(full_cohort$stroke_date), "No", "Yes")

# Check this has been coded correctly 
xtabs(~ dementia, data = full_cohort)
xtabs(~ parkinsons, data = full_cohort)
xtabs(~ stroke, data = full_cohort)



# Create yes/no variable for mental health problems 
full_cohort$dep_ep <- ifelse(is.na(full_cohort$dep_ep_date), "No", "Yes")
full_cohort$dep_disorder <- ifelse(is.na(full_cohort$dep_disorder_date), "No", "Yes")
full_cohort$affect_disorder <- ifelse(is.na(full_cohort$affect_disorder_date), "No", "Yes")
full_cohort$other_affect <- ifelse(is.na(full_cohort$other_mood_affect_date), "No", "Yes")
full_cohort$phobic_disorder <- ifelse(is.na(full_cohort$phobic_disorder_date), "No", "Yes")
full_cohort$anxiety_disorder <- ifelse(is.na(full_cohort$anxiety_disorder_date), "No", "Yes")
full_cohort$schizophrenia <- ifelse(is.na(full_cohort$schizophrenia_date), "No", "Yes")
full_cohort$delusion_disorder <- ifelse(is.na(full_cohort$delusion_disorder_date), "No", "Yes")
full_cohort$acute_psychotic_disorder <- ifelse(is.na(full_cohort$acute_psychotic_date), "No", "Yes")
full_cohort$schizoaffective_disorder <- ifelse(is.na(full_cohort$schizoaffective_date), "No", "Yes")
full_cohort$bipolar_disorder <- ifelse(is.na(full_cohort$bipolar_date), "No", "Yes")
full_cohort$organic_mental_disorder <- ifelse(is.na(full_cohort$organic_mental_date), "No", "Yes")
full_cohort$organic_personality_disorder <- ifelse(is.na(full_cohort$organic_personality_date), "No", "Yes")
full_cohort$ocd <- ifelse(is.na(full_cohort$ocd_date), "No", "Yes")


# Create grouped variables for mood disorders
full_cohort$mood_disorder <- ifelse(rowSums(full_cohort[,c("dep_ep", "dep_disorder", "affect_disorder", "other_affect", "bipolar_disorder")] == "Yes", na.rm = TRUE) > 0, "Yes", "No")


# Create grouped variables for anxiety disorders 
full_cohort$anx_disorder <- ifelse(rowSums(full_cohort[,c("phobic_disorder", "anxiety_disorder", "ocd")] == "Yes", na.rm = TRUE) > 0, "Yes", "No")


# Create grouped variables for psychotic disorders 
full_cohort$psyc_disorder <- ifelse(rowSums(full_cohort[,c("schizophrenia", "delusion_disorder", "acute_psychotic_disorder", "schizoaffective_disorder")] == "Yes", na.rm = TRUE) > 0, "Yes", "No")



# Create grouped variables for any mental health disorder 
full_cohort <- full_cohort %>% 
  mutate(any_mh_disorder = ifelse(
    if_any(c(mood_disorder, anx_disorder, psyc_disorder), ~.x == "Yes"), 
    "Yes", "No"
  ))



# Create grouped variables for any neurological disorders 
full_cohort$any_neuro <- ifelse(rowSums(full_cohort[,c("dementia", "stroke", "parkinsons")] == "Yes", na.rm = TRUE) > 0, "Yes", "No")



### Create new column for mood disorder date, taking the first date if multiple occur

# First check how many rows have dates on them 
cols <- c("dep_ep_date", "dep_disorder_date", "affect_disorder_date", "other_mood_affect_date", "bipolar_date")

full_cohort %>% 
  filter(if_any(all_of(cols), ~ !is.na(.))) %>%
  nrow


# Now create new column for mood disorder date
full_cohort <- full_cohort %>%
  mutate(
    across(c(dep_ep_date, dep_disorder_date, affect_disorder_date, other_mood_affect_date, bipolar_date), ymd), 
    mood_disorder_date = pmin(dep_ep_date, dep_disorder_date, affect_disorder_date, other_mood_affect_date, bipolar_date, na.rm = TRUE)
  )

# Check this has coded properly
sum(!is.na(full_cohort$mood_disorder_date))    


### Create new column for anxiety disorder date, taking the first date if multiple occur

# First check how many rows have dates on them 
cols <- c("phobic_disorder_date", "anxiety_disorder_date", "ocd_date")

full_cohort %>% 
  filter(if_any(all_of(cols), ~ !is.na(.))) %>%
  nrow


# Now create new column for anxiety disorder date
full_cohort <- full_cohort %>%
  mutate(
    across(c(phobic_disorder_date, anxiety_disorder_date, ocd_date), ymd), 
    anx_disorder_date = pmin(phobic_disorder_date, anxiety_disorder_date, ocd_date, na.rm = TRUE)
  )



# Check this has coded properly
sum(!is.na(full_cohort$anx_disorder_date))   



### Create new column for psychotic disorder date, taking the first date if multiple occur

# First check how many rows have dates on them 
cols <- c("schizophrenia_date", "delusion_disorder_date", "acute_psychotic_date", "schizoaffective_date")

full_cohort %>% 
  filter(if_any(all_of(cols), ~ !is.na(.))) %>%
  nrow


# Now create new column for psychotic disorder date
full_cohort <- full_cohort %>%
  mutate(
    across(c(schizophrenia_date, delusion_disorder_date, acute_psychotic_date, schizoaffective_date), ymd), 
    psyc_disorder_date = pmin(schizophrenia_date, delusion_disorder_date, acute_psychotic_date, schizoaffective_date, na.rm = TRUE)
  )



# Check this has coded properly
sum(!is.na(full_cohort$psyc_disorder_date))   



# Create new column for any mental health disorder date, taking the first date if multiple occur

# First check how many rows have dates on them 
cols <- c("dep_ep_date", "dep_disorder_date", "affect_disorder_date", "other_mood_affect_date", "bipolar_date", 
          "phobic_disorder_date", "anxiety_disorder_date", "ocd_date", 
          "schizophrenia_date", "delusion_disorder_date", "acute_psychotic_date", "schizoaffective_date")

full_cohort %>% 
  filter(if_any(all_of(cols), ~ !is.na(.))) %>%
  nrow


# Now create new column for any mental health disorder date
full_cohort <- full_cohort %>%
  mutate(
    across(c(dep_ep_date, dep_disorder_date, affect_disorder_date, other_mood_affect_date, bipolar_date, 
             phobic_disorder_date, anxiety_disorder_date, ocd_date, 
             schizophrenia_date, delusion_disorder_date, acute_psychotic_date, schizoaffective_date), ymd), 
    any_mh_date = pmin(dep_ep_date, dep_disorder_date, affect_disorder_date, other_mood_affect_date, bipolar_date, 
                       phobic_disorder_date, anxiety_disorder_date, ocd_date, 
                       schizophrenia_date, delusion_disorder_date, acute_psychotic_date, schizoaffective_date, na.rm = TRUE)
  )


# Check this has coded properly
sum(!is.na(full_cohort$any_mh_date))   




# Create a column only showing neuro cases before mental health diagnoses dates

#### Mood disorders

# Dementia 
full_cohort$dementia_before_mood <- full_cohort$dementia_date <= full_cohort$mood_disorder_date

# Change from logical to yes/no format 
full_cohort$dementia_before_mood <- ifelse(full_cohort$dementia_before_mood, "Yes", "No")

# Stroke 
full_cohort$stroke_before_mood <- full_cohort$stroke_date <= full_cohort$mood_disorder_date

# Change from logical to yes/no format 
full_cohort$stroke_before_mood <- ifelse(full_cohort$stroke_before_mood, "Yes", "No")

# Parkinsons
full_cohort$parkinsons_before_mood <- full_cohort$parkinsons_date <= full_cohort$mood_disorder_date

# Change from logical to yes/no format 
full_cohort$parkinsons_before_mood <- ifelse(full_cohort$parkinsons_before_mood, "Yes", "No")


###### Anxiety disorders

# Dementia
full_cohort$dementia_before_anx <- full_cohort$dementia_date <= full_cohort$anx_disorder_date

# Change from logical to yes/no format 
full_cohort$dementia_before_anx <- ifelse(full_cohort$dementia_before_anx, "Yes", "No")


# Stroke
full_cohort$stroke_before_anx <- full_cohort$stroke_date <= full_cohort$anx_disorder_date

# Change from logical to yes/no format 
full_cohort$stroke_before_anx <- ifelse(full_cohort$stroke_before_anx, "Yes", "No")

# Parkinsons 
full_cohort$parkinsons_before_anx <- full_cohort$parkinsons_date <= full_cohort$anx_disorder_date

# Change from logical to yes/no format 
full_cohort$parkinsons_before_anx <- ifelse(full_cohort$parkinsons_before_anx, "Yes", "No")


### Psychotic disorder 

# Dementia 
full_cohort$dementia_before_psyc <- full_cohort$dementia_date <= full_cohort$psyc_disorder_date

# Change from logical to yes/no format 
full_cohort$dementia_before_psyc <- ifelse(full_cohort$dementia_before_psyc, "Yes", "No")

# Stroke
full_cohort$stroke_before_psyc <- full_cohort$stroke_date <= full_cohort$psyc_disorder_date

# Change from logical to yes/no format 
full_cohort$stroke_before_psyc <- ifelse(full_cohort$stroke_before_psyc, "Yes", "No")

# Parkinsons
full_cohort$parkinsons_before_psyc <- full_cohort$parkinsons_date <= full_cohort$psyc_disorder_date

# Change from logical to yes/no format 
full_cohort$parkinsons_before_psyc <- ifelse(full_cohort$parkinsons_before_psyc, "Yes", "No")

##### Any mental health disorder 

# Dementia
full_cohort$dementia_before_any_mh <- full_cohort$dementia_date <= full_cohort$any_mh_date

# Change from logical to yes/no format 
full_cohort$dementia_before_any_mh <- ifelse(full_cohort$dementia_before_any_mh, "Yes", "No")

# Stroke
full_cohort$stroke_before_any_mh <- full_cohort$stroke_date <= full_cohort$any_mh_date

# Change from logical to yes/no format 
full_cohort$stroke_before_any_mh <- ifelse(full_cohort$stroke_before_any_mh, "Yes", "No")

# Parkinsons
full_cohort$parkinsons_before_any_mh <- full_cohort$parkinsons_date <= full_cohort$any_mh_date


# Change from logical to yes/no format 
full_cohort$parkinsons_before_any_mh <- ifelse(full_cohort$parkinsons_before_any_mh, "Yes", "No")



# Create a column showing mood disorder diagnoses the number of years after dementia diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    dementia_before_mood_years = if_else(
      !is.na(mood_disorder_date) & !is.na(dementia_date) & 
        mood_disorder_date >= dementia_date, 
      round(as.numeric(difftime(mood_disorder_date, dementia_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ dementia_before_mood_years, data = full_cohort)

# Create a column showing anxiety disorder diagnoses the number of years after dementia diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    dementia_before_anx_years = if_else(
      !is.na(anx_disorder_date) & !is.na(dementia_date) & 
        anx_disorder_date >= dementia_date, 
      round(as.numeric(difftime(anx_disorder_date, dementia_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ dementia_before_anx_years, data = full_cohort)



# Create a column showing psychotic disorder diagnoses the number of years after dementia diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    dementia_before_psyc_years = if_else(
      !is.na(psyc_disorder_date) & !is.na(dementia_date) & 
        psyc_disorder_date >= dementia_date, 
      round(as.numeric(difftime(psyc_disorder_date, dementia_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ dementia_before_psyc_years, data = full_cohort)



# Create a column showing any mental health diagnoses the number of years after dementia diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    dementia_before_any_mh_years = if_else(
      !is.na(any_mh_date) & !is.na(dementia_date) & 
        any_mh_date >= dementia_date, 
      round(as.numeric(difftime(any_mh_date, dementia_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ dementia_before_any_mh_years, data = full_cohort)




# Create a column showing mood disorder diagnoses the number of years after stroke diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    stroke_before_mood_years = if_else(
      !is.na(mood_disorder_date) & !is.na(stroke_date) & 
        mood_disorder_date >= stroke_date, 
      round(as.numeric(difftime(mood_disorder_date, stroke_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ stroke_before_mood_years, data = full_cohort)



# Create a column showing anxiety disorder diagnoses the number of years after stroke diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    stroke_before_anx_years = if_else(
      !is.na(anx_disorder_date) & !is.na(stroke_date) & 
        anx_disorder_date >= stroke_date, 
      round(as.numeric(difftime(anx_disorder_date, stroke_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ stroke_before_anx_years, data = full_cohort)



# Create a column showing psychotic disorder diagnoses the number of years after stroke diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    stroke_before_psyc_years = if_else(
      !is.na(psyc_disorder_date) & !is.na(stroke_date) & 
        psyc_disorder_date >= stroke_date, 
      round(as.numeric(difftime(psyc_disorder_date, stroke_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ stroke_before_psyc_years, data = full_cohort)



# Create a column showing any mental health diagnoses the number of years after stroke diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    stroke_before_any_mh_years = if_else(
      !is.na(any_mh_date) & !is.na(stroke_date) & 
        any_mh_date >= stroke_date, 
      round(as.numeric(difftime(any_mh_date, stroke_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ stroke_before_any_mh_years , data = full_cohort)



# Create a column showing mood disorder diagnoses the number of years after parkinsons diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    parkinsons_before_mood_years  = if_else(
      !is.na(mood_disorder_date) & !is.na(parkinsons_date) & 
        mood_disorder_date >= parkinsons_date, 
      round(as.numeric(difftime(mood_disorder_date, parkinsons_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ parkinsons_before_mood_years, data = full_cohort)



# Create a column showing anxiety disorder diagnoses the number of years after parkinsons diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    parkinsons_before_anx_years = if_else(
      !is.na(anx_disorder_date) & !is.na(parkinsons_date) & 
        anx_disorder_date >= parkinsons_date, 
      round(as.numeric(difftime(anx_disorder_date, parkinsons_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ parkinsons_before_anx_years, data = full_cohort)



# Create a column showing psychotic disorder diagnoses the number of years after parkinsons diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    parkinsons_before_psyc_years = if_else(
      !is.na(psyc_disorder_date) & !is.na(parkinsons_date) & 
        psyc_disorder_date >= parkinsons_date, 
      round(as.numeric(difftime(psyc_disorder_date, parkinsons_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ parkinsons_before_psyc_years , data = full_cohort)



# Create a column showing any mental health diagnoses the number of years after parkinsons diagnosis 
full_cohort <- full_cohort %>%
  mutate(
    parkinsons_before_any_mh_years  = if_else(
      !is.na(any_mh_date) & !is.na(parkinsons_date) & 
        any_mh_date >= parkinsons_date, 
      round(as.numeric(difftime(any_mh_date, parkinsons_date, units = "days")) / 365.25), 
      NA_real_
    )
  )

xtabs( ~ parkinsons_before_any_mh_years, data = full_cohort)



# Rename and select demographic variables from TBI ICD dataframe
tbi_icd <- tbi_icd %>%
  rename(subjectkey = Participant.ID) 


# Create new TBI ICD column 
tbi_icd <- tbi_icd %>% 
  mutate(tbi = 'TBI') %>%
  mutate(tbi_icd = 'ICD')


# Rename and select demographic variables from TBI self-report dataframe 
tbi_self <- tbi_self %>%
  rename(subjectkey = Participant.ID) 

# Create new TBI self-report column 
tbi_self <- tbi_self %>% 
  mutate(tbi = 'TBI') %>%
  mutate(tbi_self = 'Self')


# Merge datasets
data <- full_cohort %>%
  left_join(tbi_icd, by = "subjectkey") %>%
  left_join(tbi_self, by = "subjectkey")


# Combine TBI columns 
data$tbi <- paste(data$tbi.x, data$tbi.y) 


# Recode TBI column 
data$tbi[data$tbi == 'NA NA'] <- 0
data$tbi[data$tbi == 'NA TBI'] <- 1                            
data$tbi[data$tbi == 'TBI TBI'] <- 1   
data$tbi[data$tbi == 'TBI NA'] <- 1     


# Set data type and level labels
data$tbi <- factor(data$tbi, levels = c(1,0), labels = c("Yes", "No" ))

# Remove columns no longer needed 
data <- subset(data, select = -c(tbi.x, tbi.y))


# Recode ICD TBI column 
data$tbi_icd <- ifelse(is.na(data$tbi_icd), "No", "Yes")


# Recode self-report TBI column  
data$tbi_self <- ifelse(is.na(data$tbi_self), "No", "Yes")


# Combine TBI columns 
data$tbi_report <- paste(data$tbi_icd, data$tbi_self) 

xtabs(~ tbi_report, data = data)



# Recode TBI column 
data$tbi_report[data$tbi_report == 'No No'] <- 'None'
data$tbi_report[data$tbi_report == 'Yes No'] <- 'ICD-coded TBI'                            
data$tbi_report[data$tbi_report == 'No Yes'] <- 'Self-reported TBI'   
data$tbi_report[data$tbi_report == 'Yes Yes'] <- 'Both'     


# Set data type and level labels
data$tbi_report <- factor(data$tbi_report, 
                          levels = c("ICD-coded TBI", "Self-reported TBI", "Both", "None" ))




# Create a new column where only dementia cases are reported for TBI/no TBI 
data$dementia_tbi <- ifelse(data$dementia == 'Yes', 
                            ifelse(data$tbi == "Yes", "TBI", "No TBI"), NA) 


# Set data type and level labels
data$dementia_tbi <- factor(data$dementia_tbi, levels = c("TBI", "No TBI" ))


# Create a new column where only stroke cases are reported for TBI/no TBI 
data$stroke_tbi <- ifelse(data$stroke == 'Yes', 
                          ifelse(data$tbi == "Yes", "TBI", "No TBI"), NA) 


# Set data type and level labels
data$stroke_tbi <- factor(data$stroke_tbi, levels = c("TBI", "No TBI" ))


# Create a new column where only parkinsons cases are reported for TBI/no TBI 
data$parkinsons_tbi <- ifelse(data$parkinsons == 'Yes', 
                              ifelse(data$tbi == "Yes", "TBI", "No TBI"), NA) 


# Set data type and level labels
data$parkinsons_tbi <- factor(data$parkinsons_tbi, levels = c("TBI", "No TBI" ))


# Rename variables for better readability 
label(data$age) <- "Age"
label(data$sex) <- "Sex" 
label(data$ethnicity) <- "Ethnicity" 
label(data$tbi) <- "TBI" 
label(data$alcohol_freq) <- "Alcohol intake frequency" 
label(data$deprivation_index) <- "Townsend deprivation index" 
label(data$dementia) <- "Number of dementia cases" 
label(data$stroke) <- "Number of stroke cases" 
label(data$parkinsons) <- "Number of Parkinson's cases" 
label(data$tbi_report) <- "Source of TBI" 
label(data$mood_disorder) <- "Mood disorders" 
label(data$anx_disorder) <- "Anxiety disorders" 
label(data$psyc_disorder) <- "Psychotic disorders" 
label(data$any_mh_disorder) <- "Any mental health disorder" 


# Save new file for analyses 
write.csv(data, paste(output_dir, "data.csv", sep = ""))

# Create new dataframe with just neurological condition cases 
neuro <- data[data$any_neuro == "Yes",]


# Save new file for analyses 
write.csv(neuro, paste(output_dir, "neuro.csv", sep = ""))


###### Create ICD code dates 

# Merge datasets
icd_tbi <- icd %>%
  left_join(icd_date, by = "subjectkey")


##  1) Identify the date array columns
date_cols <- grep("Date.of.first.in.patient.diagnosis", names(icd_tbi), value = TRUE)

## 2) Expand ICD10 codes into long format 
icd_long_all <- icd_tbi %>%
  mutate(icd_list = strsplit(Diagnoses...ICD10, "\\|")) %>%
  unnest_longer(icd_list, indices_to = "icd_index", values_to = "diagnosis_icd") %>%
  mutate(
    icd_index    = as.integer(icd_index),
    first_letter = substr(diagnosis_icd, 1, 1)
  )

## 3) Expand Date arrays into long format 
date_long <- icd_tbi %>%
  select(subjectkey, all_of(date_cols)) %>%
  pivot_longer(
    cols      = all_of(date_cols),
    names_to  = "date_col",
    values_to = "icd_date"
  ) %>%
  mutate(icd_index = as.integer(str_extract(date_col, "\\d+"))) %>%
  select(subjectkey, icd_index, icd_date)

## 4) Join codes + dates by array index & parse dates 
long_joined_all <- icd_long_all %>%
  left_join(date_long, by = c("subjectkey", "icd_index")) %>%
  mutate(icd_date_parsed = suppressWarnings(dmy(icd_date)))

##  5) S/T ONLY summary (true injury dates) 
st_summary <- long_joined_all %>%
  filter(first_letter %in% c("S", "T")) %>% 
  group_by(subjectkey) %>%
  summarise(
    n_ST_codes       = n(),    
    earliest_ST_date = if (all(is.na(icd_date_parsed))) {
      as.Date(NA)
    } else {
      min(icd_date_parsed, na.rm = TRUE)
    },
    .groups = "drop"
  )

##  6) ANY-date summary: earliest date from ANY ICD code 
any_date_summary <- long_joined_all %>%
  group_by(subjectkey) %>%
  summarise(
    earliest_any_date = if (all(is.na(icd_date_parsed))) {
      as.Date(NA)
    } else {
      min(icd_date_parsed, na.rm = TRUE)
    },
    .groups = "drop"
  )

##  7) Create final dataset: use S/T date if available, otherwise any ICD date 
subject_ids <- icd_tbi %>%
  distinct(subjectkey)     # ensures one row per participant

final_optionB <- subject_ids %>%
  left_join(st_summary, by = "subjectkey") %>%
  left_join(any_date_summary, by = "subjectkey") %>%
  mutate(
    n_ST_codes = ifelse(is.na(n_ST_codes), 0L, n_ST_codes),
    injury_date_B = if_else(
      !is.na(earliest_ST_date), earliest_ST_date, earliest_any_date
    ) 
  )

## Output:
final_optionB


# Save new file for analyses 
write.csv(final_optionB, paste(output_dir, "icd_date.csv", sep = ""))


# Load dataset
neuro <- read.csv("S:/Track_TBI/Grace/Test/Output/neuro.csv")

# Load ICD date table
icd_date <- read.csv("S:/Track_TBI/Grace/Test/Output/icd_date.csv")


# Merge datasets
neuro <- neuro %>%
  left_join(icd_date, by = "subjectkey")

# Save new file for analyses 
write.csv(neuro, paste(output_dir, "neuro.csv", sep = ""))

####### DESCRIPTIVE STATISTICS ############

# Full cohort descriptive table 

# Create demographic and descriptive table for full sample 
descriptive_table_full <- data %>%
  select(age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary() %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_caption("**Table 1. Demographics of sample**")

# Create demographic and descriptive table for TBI sample 
descriptive_table_tbi <- data %>%
  select(tbi, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**TBI**") %>%
  modify_caption("**Table 1. Demographics of sample**")



# Create demographic and descriptive table for dementia sample 
descriptive_table_dementia <- data %>%
  select(dementia, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = dementia) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Dementia**") %>%
  modify_caption("**Table 1. Demographics of sample**")



# Create demographic and descriptive table for stroke sample 
descriptive_table_stroke <- data %>%
  select(stroke, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = stroke) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Stroke**") %>%
  modify_caption("**Table 1. Demographics of sample**")


# Create demographic and descriptive table for parkinsons sample 
descriptive_table_parkinsons <- data %>%
  select(parkinsons, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = parkinsons) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Parkinsons**") %>%
  modify_caption("**Table 1. Demographics of sample**")


# Merge tables
tbl <- tbl_merge(
  tbls = list(descriptive_table_full, descriptive_table_tbi, descriptive_table_dementia, descriptive_table_stroke, descriptive_table_parkinsons), 
  tab_spanner = c("**Overall**","**TBI**", "**Dementia**" , "**Stroke**", "**Parkinsons**")
)

tbl 

# Save table
table_1_filename = paste(output_dir, "descriptive_table.html", sep="")
gt::gtsave(as_gt(tbl), file = table_1_filename)
head(table_1_filename)


# Create demographic and descriptive table for TBI sample - how many have neurological disorders 
descriptive_table <- data %>%
  select(tbi, dementia, stroke, parkinsons) %>%
  tbl_summary(by = tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**TBI**") %>%
  modify_caption("**Table 1. Neurological diseases in people with/without history of TBI**")

descriptive_table      

# Save table
table_1_filename = paste(output_dir, "neuro_tbi_table.html", sep="")
gt::gtsave(as_gt(descriptive_table), file = table_1_filename)
head(table_1_filename)



# Create demographic and descriptive table for sample - TBI source and neurological disorder 
descriptive_table <- data %>%
  select(tbi_report, dementia, stroke, parkinsons) %>%
  tbl_summary(by = tbi_report) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**TBI**") %>%
  modify_caption("**Table 1. Neurological diseases in people with/without history of TBI**")

descriptive_table      



### Create psychiatric disorders table with neurological condition stratified by TBI 

# Create demographic and descriptive table for full sample 
descriptive_psyc_table_full <- data %>%
  select(mood_disorder, anx_disorder, psyc_disorder, any_mh_disorder) %>%
  tbl_summary() %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_caption("**Table 1. Demographics of sample**")



# Create TBI table
descriptive_psyc_table_tbi <- data %>%
  select(tbi, mood_disorder, anx_disorder, psyc_disorder, any_mh_disorder) %>%
  tbl_summary(by = tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**TBI**") %>%
  modify_caption("**Table 1. Psychological problems in people with TBI**")



# Create dementia table
descriptive_psyc_table_dementia <- data %>%
  select(dementia_tbi, mood_disorder, anx_disorder, psyc_disorder, any_mh_disorder) %>%
  tbl_summary(by = dementia_tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Dementia**") %>%
  modify_caption("**Table 1. Psychological problems in people with dementia with/without history of TBI**")


# Create stroke table 
descriptive_psyc_table_stroke <- data %>%
  select(stroke_tbi, mood_disorder, anx_disorder, psyc_disorder, any_mh_disorder) %>%
  tbl_summary(by = stroke_tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Stroke**") %>%
  modify_caption("**Table 1. Psychological problems in people with stroke with/without history of TBI**")


# Create parkinson's table 
descriptive_psyc_table_parkinsons <- data %>%
  select(parkinsons_tbi, mood_disorder, anx_disorder, psyc_disorder, any_mh_disorder) %>%
  tbl_summary(by = parkinsons_tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Parkinsons**") %>%
  modify_caption("**Table 1. Psychological problems in people with Parkinsons with/without history of TBI**")


# Merge tables
tbl <- tbl_merge(
  tbls = list(descriptive_psyc_table_full, descriptive_psyc_table_tbi, descriptive_psyc_table_dementia, descriptive_psyc_table_stroke, descriptive_psyc_table_parkinsons), 
  tab_spanner = c("**Overall**", "**TBI**", "**Dementia**" , "**Stroke**", "**Parkinsons**")
)

tbl 

# Save table
table_1_filename = paste(output_dir, "neuro_tbi_psyc_descriptive_table.html", sep="")
gt::gtsave(as_gt(tbl), file = table_1_filename)
head(table_1_filename)

