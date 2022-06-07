#**********************************************************************************************************************#
# Acid_Base Script Testing
#
# Authors: Sierra
# Reviewers:
#**********************************************************************************************************************#

########################################################################################################################*
#load libraries

library(coronaenv)
library(tidyverse)
library(readxl)

########################################################################################################################*
# Compare test set results to current functions

#Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Acid-Base Test Set ----
# Format csv
testset <- read.csv("Acid-Base Test Set.csv") %>%
  mutate(Temp0 = Temp) %>%
  # Change integers to numbers
  mutate(across(where(is.integer), as.numeric))

# Apply functions
testing <- testset %>%
  # Apply chem_calc_ph function
  pmap_dfr(chem_calc_ph) %>%
  # Apply hardness calc
  mutate(Hardness = as.numeric(pmap(list(Ca.Conc, Mg.Conc), calculate_hardness))) %>%
  # Apply chem_calc with corrosivity, but no chemical because it's already been added
  rename(Alk0 = Alkalinity, pH0 = pH, Temp0 = Temp) %>%
  pmap_dfr(chem_calc_ph, corrosivity = TRUE)

# Check delta between test set and data
deltacheck <- testing %>%
  pivot_longer(ends_with(".x"), names_to = "Parameter", values_to = "TestSet") %>%
  mutate(Parameter = gsub("[.]x", "", Parameter)) %>%
  pivot_longer(-c(Condition, Parameter, TestSet), names_to = "ParameterC", values_to = "CurrentRun") %>%
  subset(Parameter == ParameterC) %>%
  select(-c(ParameterC)) %>%
  mutate(Delta = abs(CurrentRun - TestSet))

# Blending Test Set ----
testset_b <- read.csv("Blend Test Set.csv") %>%
  rename(ID = 1) %>%
  group_by(ID) %>%
  nest() %>%
  mutate(fin = map(.x = data, blend_calc_ph, df.return = TRUE)) %>%
  unnest() %>%
  select(!ends_with("1")) %>%
  mutate(DeltapH = abs(pH.x - pH.fin),
         DeltaAlk = abs(Alk.x - Alk.fin))

# Chemical Dose for Target pH Test Set ----
testset_c <- read.csv("Chemical Dose Test Set.csv") %>%
  rename(Condition = 1)

testing_c <- testset_c %>%
  pmap_dfr(.f = solve_chem_dose) %>%
  mutate(delta = abs(solve.x - chem.solve))

########################################################################################################################*
########################################################################################################################*
########################################################################################################################*

# This is only used for rewriting the test set CSVs and should rarely be run

# # Acid base test set
# newfile1 <- testing %>%
#   select(-c(contains(".x"))) %>%
#   rename_with(~paste(.x, ".x", sep = "")) %>%
#   rename(Condition = Condition.x)
#
# newfile <- testset %>%
#   select(1:27) %>%
#   full_join(newfile1)
#
# write.csv(newfile, "Acid-Base Test Set.csv", row.names = FALSE)

# # Blending test set
# newfile <- testset_b %>%
#   select(-c(pH.x, Alk.x, Temp.fin, IS, DeltapH, DeltaAlk)) %>%
#   rename(pH.x = pH.fin, Alk.x = Alk.fin)
#
# write.csv(newfile, "Blend Test Set.csv", row.names = FALSE)

# Chemical dose test set
newfile <- testset_c %>%
  select(-c(solve.x)) %>%
  full_join(subset(testing_c, select = c(Condition, chem.solve))) %>%
  rename(solve.x = chem.solve)


