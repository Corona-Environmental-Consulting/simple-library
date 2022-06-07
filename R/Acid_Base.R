#**********************************************************************************************************************#
# Find pH after chemical addition
# Acid-base equilibrium only, no solubility equilibrium
#
# Authors: Chris & Sierra
# Reviewers: Sheldon 2021/03/09, Libby 2021/04/08, Juliette 2021/06/30
#**********************************************************************************************************************#


# INPUTS ----
# Required inputs
#   pH - SU before chemical addition
#   Alkalinity - mg/L as CaCO3
#   Temperature (deg C)

# Optional inputs
#   Conductivity (umho/cm) - default of 500
#   Final temperature (if different during chemical addition) - default equal to initial temp
#   Ion Concentrations (Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc)
#   Chemical doses listed below (mg/L) - default of 0


########################################################################################################################*
#----------------------------------------------------------------------------------------------------------------------#*
# HELPER FUNCTIONS ----
#----------------------------------------------------------------------------------------------------------------------#*


# * CONVERSIONS/CORRECTIONS ----

# Temperature based acidity constants corrected for ionic strength
acidity_const <- function(TempK, IS, ID) {

  Dielectric <- 78.54 * (1 - (0.004579 * (TempK - 298)) + 11.9 * 10 ^ -6 * (TempK - 298) ^2 + 28 * 10 ^ -9 * (TempK - 298) ^ 3)

  #monovalent activity correction
  Mono <- 10 ^ -(1824335 * (Dielectric * (TempK)) ^ (-3 / 2) * (sqrt(IS) / (1 + sqrt(IS))) - 0.3 * IS)
  #divalent activity correction
  Di <- 10 ^ -(1824335 * (Dielectric * (TempK)) ^ (-3 / 2) * 4 * (sqrt(IS) / (1 + sqrt(IS))) - 0.3 * IS)
  #trivalent activity correction
  Tri <- 10 ^ -(1824335 * (Dielectric * (TempK)) ^ (-3 / 2) * 9 * (sqrt(IS) / (1 + sqrt(IS))) - 0.3 * IS)

  # Acidity constants
  Kw <- 4787.3 / TempK + 7.1321 * log10(TempK) + 0.010365 * TempK - 22.801
  K1 <- 17052 / TempK + 215.21 * log10(TempK) - 0.12675 * TempK - 545.56
  K2 <- 2902.39 / TempK + 0.02379 * TempK - 6.498
  KOCl <- 0.0253 * (TempK) + 3000 / (TempK) - 10.0686
  KSO4 <- 2 + (-3.59) / 0.008314 * (1 / 298.15 - 1 / TempK)
  K1PO4 <- 2.148 + (-8) / 0.008314 * (1 / 298.15 - 1 / TempK)
  K2PO4 <- 7.2 + (3.6) / 0.008314 * (1 / 298.15 - 1 / TempK)
  K3PO4 <- 12.35 + (16) / 0.008314 * (1 / 298.15 - 1 / TempK)
  KCaCO3 <- 171.9065 + 0.077993 * TempK - 2839.319 / TempK - 71.595 * log10(TempK)

  if (ID == "cKw") {
    10 ^ -Kw / (Mono ^ 2)
  } else if (ID == "cK1") {
    10 ^ -K1 / (Mono ^ 2)
  } else if (ID == "cK2"){
    10 ^ -K2 / Di
  } else if (ID == "cKOCl") {
    10 ^ -KOCl / (Mono ^ 2)
  } else if (ID == "cKSO4") {
    10 ^ -KSO4 / Di
  } else if (ID == "cK1PO4") {
    10 ^ -K1PO4 / (Mono ^ 2)
  } else if (ID == "cK2PO4") {
    10 ^ -K2PO4 / Di
  } else if (ID == "cK3PO4") {
    10 ^ -K3PO4 * Di / (Mono * Tri)
  } else if (ID == "cKCaCO3") {
    10 ^ -KCaCO3 / (Di ^ 2)
  } else if (ID == "K2") {
    K2
  } else if (ID == "KCaCO3") {
    KCaCO3
  } else if (ID == "Mono") {
    Mono
  } else if (ID == "Di") {
    Di
  }
}

#############################################################*
#' Unit conversions
#'
#' This function does unit conversions of common ions.
#'
#' @param Conc concentration of the compound in the current units
#' @param Compound chemical formula for the compound to be converted in quotes. Currently accepts:
#' "Na", "K", "Cl", "HCO3", "Ca", "Mg", "SO4", "CO3"
#' @param end desired units/L of conversion in quotes. Currently accepts: "mg", "mol", "eq", "mgCaCO3"
#' @param start starting units/L of conversion in quotes. Currently accepts: "mg", "mol", "eq", "mgCaCO3"
#'
#' @examples convert(50, "Ca", end = "mgCaCO3", start = "mg")
#'
#' @export
#'
convert <- function(Conc, Compound, end = "mol", start = "mg") {

  # Check that start and end are entered correctly
  if (!(start == "mg" | start == "mol" | start == "eq" | start == "mgCaCO3")) {
    stop("Unrecognized start units")
  }
  if (!(end == "mg" | end == "mol" | end == "eq" | end == "mgCaCO3")) {
    stop("Unrecognized end units")
  }


  # Determine number of ions for equivalents calculation (if relevant)
  if ((end == "eq" | end == "mgCaCO3") &
      (Compound == "Na" | Compound == "K" | Compound == "Cl" | Compound == "HCO3")) {
    eq <- 1
  } else if ((end == "eq" | end == "mgCaCO3") &
             (Compound == "Ca" | Compound == "Mg" | Compound == "SO4" | Compound == "CO3")) {
    eq <- 2
  } else if ((start == "eq" | start == "mgCaCO3") &
             (Compound == "Ca" | Compound == "Mg" | Compound == "SO4" | Compound == "CO3")) {
    eq <- 1/2
  } else {
    eq <- 1
  }

  # If converting to or from CaCO3, need molar mass of CaCO3, if not, use 1
  if (end == "mgCaCO3") {
    caco3 <- 100.0869 / 2 # g/eq
  } else if (start == "mgCaCO3") {
    caco3 <- 2 / 100.0869 # eq/g
  } else {
    caco3 <- 1
  }

  # If converting from mg/L to mol/L or eq/L, divide by 1000 and divide by molar mass (raise to -1)
  # mg/L to mgCaCO3/L doesn't need the 1000 conversion, but still divide by molar mass (raise to -1)
  # from mol/L to mg/L, multiply by molar mass (raise to 1) and multiply by 1000 to get mg
  # from mol to eq or CaCO3, molar mass not used (raise to 0) and no 1000 conversion needed
  if (start == "mg" & (end == "mol" | end == "eq")) {
    mm <- -1
    mg <- 1/1000
  } else if (start == "mg" & end == "mgCaCO3"){
    mm <- -1
    mg <- 1
  } else if (start == "mol" & end == "mg") {
    mm <- 1
    mg <- 1000
  } else if (start == "mol" & end == "mgCaCO3") {
    mm <- 0
    mg <- 1000
  } else if (start == "mol" & end == "eq") {
    mm <- 0
    mg <- 1
  } else if (start == "mgCaCO3" & end == "mg"){
    mm <- 1
    mg <- 1
  }

  # Conversions depend on molar mass
  if (start == end) {
    Conc
  } else if (Compound == "Na") {
    Conc * mg * 22.98977 ^ mm * eq * caco3
  } else if (Compound == "Ca") {
    Conc * mg * 40.078 ^ mm * eq * caco3
  } else if (Compound == "Mg") {
    Conc * mg * 24.305 ^ mm * eq * caco3
  } else if (Compound == "K") {
    Conc * mg * 39.0983 ^ mm * eq * caco3
  } else if (Compound == "Cl") {
    Conc * mg * 35.453 ^ mm * eq * caco3
  } else if (Compound == "SO4") {
    Conc * mg * (32.065 + 4*15.9994) ^ mm * eq * caco3
  } else if (Compound == "HCO3") {
    Conc * mg * 61.0168 ^ mm * eq * caco3
  } else if (Compound == "CO3") {
    Conc * mg * 60.008 ^ mm * eq * caco3
  } else if (Compound == "CaCO3") {
    Conc * mg * 100.0869 ^ mm * eq * caco3
  } else if (Compound == "Fe" & (end == "mg" | end == "mol") & (start == "mg" | start == "mol")) {
    Conc * mg * 55.845 ^ mm * eq * caco3
  } else if (Compound == "Mn" & (end == "mg" | end == "mol") & (start == "mg" | start == "mol")) {
    Conc * mg * 54.938 ^ mm * eq * caco3
  } else{
    warning("Compound not included in convert function, NA produced.")
    NA
  }

}

########################################################################################################################*
# * WQ PARAMETER CALCULATIONS ----

#' Hardness Calculation
#'
#' This function takes Calcium and Magnesium concentrations and calculates total hardness
#'
#' @param Ca.Conc calcium concentration mg/L
#' @param Mg.Conc magnesium concentration mg/L
#'
#' @export
#'
calculate_hardness <- function(Ca.Conc, Mg.Conc) {
  Ca.Hardness <- convert(Ca.Conc, "Ca", "mgCaCO3")
  Mg.Hardness <- convert(Mg.Conc, "Mg", "mgCaCO3")
  Ca.Hardness + Mg.Hardness
}

#############################################################*
#' TDS Calculation
#'
#' This function takes alkalinity and ions and calculates TDS
#'
#' @param Alk Alkalinity in mg as CacO3/L
#' @param Na.Conc sodium concentration in mg/L
#' @param Ca.Conc calcium concentration mg/L
#' @param Mg.Conc magnesium concentration mg/L
#' @param K.Conc potassium concentration mg/L
#' @param Cl.Conc chloride concentration mg/L
#' @param SO4.Conc sulfate concentration mg/L
#'
#' @export
#'
calculate_tds <- function(Alk, Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc) {
  (Na.Conc + Ca.Conc + Mg.Conc + K.Conc + Cl.Conc + SO4.Conc) + 0.6 * Alk
}

#############################################################*
# Ionic Strength Calculations
calculate_is <- function(Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc, AlkMol){
  0.5 * (convert(Na.Conc, "Na") + convert(Ca.Conc, "Ca") *4 + convert(Mg.Conc, "Mg")*4 +
           convert(K.Conc, "K") + convert(Cl.Conc, "Cl") + convert(SO4.Conc, "SO4")*4) +
    0.5 * AlkMol
}

correlate_is_cond <- function(Cond) {
  1.6 * 10 ^ (-5) * Cond
}

correlate_is_tds <- function(TDS) {
  2.5 * 10 ^ (-5) * TDS
}

#############################################################*
#' Ion balance
#'
#' This function takes ions, converts them to mg as CaCO3/L and determines the balance.
#' It can also calculate the difference between the measured TDS and the calculated TDS from the ion concentrations.
#' Best used in a dataframe.
#'
#' @param Temp0 temperature in degrees C (optional, defaults to 20)
#' @param Alk0 alkalinity in mg CaCO3/L
#' @param pH0 pH in SU
#' @param TDS total dissolved solids (TDS) measured (optional, but necessary if you want the difference between calculated and measured)
#' @param Na.Conc sodium concentration in mg/L
#' @param Ca.Conc calcium concentration mg/L
#' @param Mg.Conc magnesium concentration mg/L
#' @param K.Conc potassium concentration mg/L (optional, defaults to 0)
#' @param Cl.Conc chloride concentration mg/L
#' @param SO4.Conc sulfate concentration mg/L
#'
#' @examples data <- data %>% pmap_drf(balance_ions)
#'
#' @export
#'
balance_ions <- function(Temp0 = 20, Alk0, pH0, TDS = NA, Na.Conc, Ca.Conc, Mg.Conc, K.Conc = 0, Cl.Conc, SO4.Conc, ...){

  outputs <- calculate_carbonate(Temp0, Alk0, pH0, Cond = 500, TDS, IS = NA,
                                 Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc)
  HCO3.mol <- outputs$iniHCO3
  CO3.mol <- outputs$iniCO3

  Na <- convert(Na.Conc, "Na", end = "mgCaCO3")
  K <- convert(K.Conc, "K", end = "mgCaCO3")
  Ca <- convert(Ca.Conc, "Ca", end = "mgCaCO3")
  Mg <- convert(Mg.Conc, "Mg", end = "mgCaCO3")
  Cl <- convert(Cl.Conc, "Cl", "mgCaCO3")
  SO4 <- convert(SO4.Conc, "SO4", "mgCaCO3")
  HCO3 <- convert(HCO3.mol, "HCO3", end = "mgCaCO3", start = "mol")
  CO3 <- convert(CO3.mol, "CO3", end = "mgCaCO3", start = "mol")

  TDS_calc <- calculate_tds(Alk = Alk0, Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc)

  cations <- Na + K + Ca + Mg
  anions <- Cl + SO4 + HCO3 + CO3
  Ion.Balance <- cations - anions
  Ion.Balance.pct <- (cations - anions) / mean(c(cations, anions))
  TDS.diff <- TDS - TDS_calc

  list(Temp0 = Temp0, Alk0 = Alk0, pH0 = pH0,
       Na.Conc = Na.Conc, Ca.Conc = Ca.Conc, Mg.Conc = Mg.Conc, K.Conc = K.Conc, Cl.Conc = Cl.Conc, SO4.Conc = SO4.Conc,
       ..., Ion.Balance = Ion.Balance, Ion.Balance.pct = Ion.Balance.pct, TDS.diff = TDS.diff,
       Na.IB = Na, K.IB = K, Ca.IB = Ca, Mg.IB = Mg, Cl.IB = Cl, SO4.IB = SO4, HCO3.IB = HCO3, CO3.IB = CO3)
}

########################################################################################################################*
# * PH/ALK CALCULATIONS ----


# pH or H+ Calculation
calculate_ph <- function(CTCO3, CBA, cK1, cK2, cKw,
                         CTSO4 = 0, CTPO4 = 0, CTOCl = 0,
                         cKSO4 = .1, cK1PO4 = .1, cK2PO4 = .1, cK3PO4 = .1, cKOCl = .1,
                         output.ph = TRUE) {

  #Solve for pH after chemical addition
  h_solve_fn <- function(H.Final) {
    cKw / H.Final +
      (2 + H.Final / cKSO4) * (CTSO4 / (H.Final / cKSO4 + 1)) +
      (H.Final ^ 2 / cK2PO4 / cK3PO4 + 2 * H.Final / cK3PO4 + 3) *
      (CTPO4 / (H.Final ^ 3 / cK1PO4 / cK2PO4 / cK3PO4 + H.Final ^ 2 / cK2PO4 / cK3PO4 + H.Final / cK3PO4 + 1)) +
      (H.Final / cK2 + 2) * (CTCO3 / (H.Final ^ 2 / cK1 / cK2 + H.Final / cK2 + 1)) +
      CTOCl / (H.Final / cKOCl + 1) - H.Final - CBA
  }

  hsolved <- uniroot(h_solve_fn, c(10 ^ -1, 10 ^ -14), tol = 1e-14)

  HConc.Final <- as.numeric(hsolved[1])

  pH.Final <- -log10(HConc.Final)

  if (output.ph) {
    return(pH = pH.Final)
  } else {
    return(H.Conc = HConc.Final)
  }

}

# pH or H+ Calculation
# pass post blend alkalinity (meq/L) & CTCO3 (mmol/L), also pass cKw, cK1, cK2

blend_ph <- function(CTCO3.Blend, Alk.Blend, cK1, cK2, cKw,
                      output.ph = TRUE) {

  # Convert alkalinity to equivalents/L
  AlkEq <- Alk.Blend * 2 / 100.0869 / 1000

  #Solve for pH after blend
  h_blend_fn <- function(H.Blend) {
    AlkEq -
      CTCO3.Blend * (H.Blend * cK1 + 2 * cK1 * cK2) / (H.Blend^2 + H.Blend * cK1 + cK1 * cK2) -
      (cKw / H.Blend-H.Blend)
  }

  hsolved <- uniroot(h_blend_fn, c(10 ^ -1, 10 ^ -14), tol = 1e-14)

  HConc.Final <- as.numeric(hsolved[1])

  pH.Final <- -log10(HConc.Final)

  if (output.ph) {
    return(pH = pH.Final)
  } else {
    return(H.Conc = HConc.Final)
  }
}

#############################################################*
# pH and alkalinity (from chem_addition output)

calculate_alkph <- function(Ion.Output) {

  H.Conc <- Ion.Output$HConc.Final
  pH.Final <- -log10(H.Conc)

  CTCO3 <- Ion.Output$CTCO3
  cK1 <- Ion.Output$cK1
  cK2 <- Ion.Output$cK2
  cKw <- Ion.Output$cKw
  HCO3 <- CTCO3 * H.Conc * cK1 / (H.Conc ^ 2 + H.Conc * cK1 + cK1 * cK2)
  CO3 <- CTCO3 * cK1 * cK2 / (H.Conc ^ 2 + H.Conc * cK1 + cK1 * cK2)
  OH <- cKw / H.Conc
  Alk.Final <- (HCO3 + 2 * CO3 + OH - H.Conc) * 50.04345 * 1000

  return(list(pH = pH.Final, Alk = Alk.Final, HCO3 = HCO3, CO3 = CO3))
}


########################################################################################################################*
#----------------------------------------------------------------------------------------------------------------------#*
# MAIN FUNCTIONS ----
#----------------------------------------------------------------------------------------------------------------------#*


# * CALCULATE CARBONATE IONS ----

# Starting parameters (no chem addition)
calculate_carbonate <- function(Temp0, Alk0, pH0, Cond = 500, TDS = NA, IS = NA,
                                Na.Conc = NA, Ca.Conc = NA, Mg.Conc = NA, K.Conc = 0, Cl.Conc = NA, SO4.Conc = NA, ...) {
  # Unit conversions
  TempK <- Temp0 + 273.15 # C to K
  AlkMol <- Alk0 * 2 / 100.0869 / 1000 #eq/L
  HConc <- 10 ^ -pH0

  # Determine ionic strength
  if(!is.na(IS)){
    IS <- IS
  } else if (!is.na(Na.Conc) & !is.na(Ca.Conc) & !is.na(Mg.Conc) & !is.na(Cl.Conc) & !is.na(SO4.Conc)){
    IS <- calculate_is(Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc, AlkMol)
  } else if (!is.na(TDS)) {
    IS <- correlate_is_tds(TDS)
  } else {
    IS <- correlate_is_cond(Cond)
  }

  # Acidity constants corrected for ionic strength
  cKw <- acidity_const(TempK, IS, "cKw")
  cK1 <- acidity_const(TempK, IS, "cK1")
  cK2 <- acidity_const(TempK, IS, "cK2")

  # Fraction of carbonate species
  alpha1 <- cK1 * HConc / (HConc ^ 2 + HConc * cK1 + cK1 * cK2)
  alpha2 <- cK1 * cK2 / (HConc ^ 2 + HConc * cK1 + cK1 * cK2)

  # Calculate OH concentration
  OHConc <- cKw / HConc

  # Calculate total carbonate species concentration
  iniCTCO3 <- (AlkMol + HConc - OHConc) / (alpha1 + 2 * alpha2)
  # Calculate bicarb and carb concentrations
  iniHCO3 <- alpha1 * iniCTCO3
  iniCO3 <- alpha2 * iniCTCO3
  # Calculate remaining ion concentration
  iniCBA <- iniHCO3 + 2 * iniCO3 + OHConc - HConc

  return(list(..., iniCTCO3 = iniCTCO3, iniCBA = iniCBA, IS = IS, Temp0 = Temp0, HConc.Final = HConc,
              Na.Conc = Na.Conc, Ca.Conc = Ca.Conc, Mg.Conc = Mg.Conc, K.Conc = K.Conc, Cl.Conc = Cl.Conc, SO4.Conc = SO4.Conc,
              iniHCO3 = iniHCO3, iniCO3 = iniCO3))
}

########################################################################################################################*
# * CHEMICAL ADDITION ----

chem_addition <- function(iniCTCO3, iniCBA, IS, Temp,
                          HCl.Conc = 0, H2SO4.Conc = 0, H3PO4.Conc = 0, CO2.Gas = 0,
                          NaOH.Conc = 0, Na2CO3.Conc = 0, NaHCO3.Conc = 0, Lime.Conc = 0, MgOH2.Conc = 0,
                          Cl2.Gas = 0, Cl2.Conc = 0, CaOCl2.Conc = 0,
                          Alum.Conc = 0, FerricCl.Conc = 0, FerricSO4.Conc = 0, CaCO3.Conc = 0,
                          Na.Conc = NA, Ca.Conc = NA, Mg.Conc = NA, K.Conc = 0, Cl.Conc = NA, SO4.Conc = NA, ...) {

  # Temp final allows for different temp when chems are added. Default matches Temp0
  TempK <- Temp + 273.15 # C to K

  # Total carbonate added (mol)
  CTCO3 <- iniCTCO3 + (CO2.Gas/44.0095 + Na2CO3.Conc/105.98844 + NaHCO3.Conc/84.00661 + CaCO3.Conc/100.0869) / 1000
  # Total chlorine added
  CTOCl <- (Cl2.Gas + Cl2.Conc + CaOCl2.Conc) / 70.906 / 1000
  #Sodium
  CTNa <- (NaOH.Conc/39.9971 + 2 * Na2CO3.Conc/105.98844 + NaHCO3.Conc/84.0061 + Cl2.Conc/70.906) / 1000
  #Chloride
  CTCl <- (HCl.Conc/36.461 + Cl2.Gas/70.906 + 3 * FerricCl.Conc/162.2) / 1000
  #Calcium
  CTCa <- (Lime.Conc/74.09268 + 0.5 * CaOCl2.Conc/70.906 + CaCO3.Conc/100.0869) / 1000
  #Magnesium
  CTMg <- MgOH2.Conc/58.31968 / 1000
  #Sulfate
  CTSO4 <- (H2SO4.Conc/98.079 + 3 * Alum.Conc/594.4 + 3 * FerricSO4.Conc/400) / 1000
  #Phosphate
  CTPO4 <- H3PO4.Conc/97.995181 / 1000
  # Remaining ions with chemical addition
  CBA <- iniCBA + CTNa + 2 * CTCa + 2 * CTMg - CTCl

  # Recalculate values with chemical addition
  IS <- IS + 0.5 * (CTNa + CTCl + CTOCl) +
    0.5 * (CTSO4 + CTCa + CTMg + CTCO3 - iniCTCO3) * 4 + 0.5 * CTPO4 * 9

  cKw <- acidity_const(TempK, IS, "cKw")
  cK1 <- acidity_const(TempK, IS, "cK1")
  cK2 <- acidity_const(TempK, IS, "cK2")

  cKOCl <- acidity_const(TempK, IS, "cKOCl")
  cKSO4 <- acidity_const(TempK, IS, "cKSO4")
  cK1PO4 <- acidity_const(TempK, IS, "cK1PO4")
  cK2PO4 <- acidity_const(TempK, IS, "cK2PO4")
  cK3PO4 <- acidity_const(TempK, IS, "cK3PO4")

  # Calculate final concentration of ions used for corrosivity calcs
  finCTNa <- ifelse(!is.na(Na.Conc), CTNa + convert(Na.Conc, "Na"), CTNa)
  finCTCl <- ifelse(!is.na(Cl.Conc), CTCl + convert(Cl.Conc, "Cl"), CTCl)
  finCTCa <- ifelse(!is.na(Ca.Conc), CTCa + convert(Ca.Conc, "Ca"), CTCa)
  finCTMg <- ifelse(!is.na(Mg.Conc), CTMg + convert(Mg.Conc, "Mg"), CTMg)
  finCTSO4 <- ifelse(!is.na(SO4.Conc), CTSO4 + convert(SO4.Conc, "SO4"), CTSO4)

  # Determine pH
  HConc.Final <- calculate_ph(CTCO3 = CTCO3, CBA = CBA, CTSO4 = CTSO4, CTPO4 = CTPO4, CTOCl = CTOCl,
               cK1 = cK1, cK2 = cK2, cKw = cKw, cKSO4 = cKSO4, cK1PO4 = cK1PO4, cK2PO4 = cK2PO4, cK3PO4 = cK3PO4,
               cKOCl = cKOCl, output.ph = FALSE)

  return(list(..., TempK = TempK, CTCO3 = CTCO3, CTOCl = CTOCl, finCTNa = finCTNa, finCTCl = finCTCl,
              finCTCa = finCTCa, finCTMg = finCTMg, finCTSO4 = finCTSO4, CTPO4 = CTPO4,
              CBA = CBA, IS = IS, HConc.Final = HConc.Final, cK1 = cK1, cK2 = cK2, cKw = cKw))

}

########################################################################################################################*

# * BLENDING FUNCTION ----

# Inputs: Dataframe with the following columns
#   - Ratio
#   - pH0
#   - Alk0
#   - Temp0
#   - Cond (optional)
#   - TDS (optional)
#   - Ionic Strength (optional)
#   - Na, Ca, Mg, K, Cl, SO4 concentrations (optional)

calculate_blend_raw <- function(blend) {

  # Determine finished water parameters for IS determination
  if (!is.null(blend$Cond)){
    Cond <- sum(blend$Cond * blend$Ratio)
  } else if (!is.null(blend$TDS)) {
    TDS <- sum(blend$TDS * blend$Ratio)
  } else if (!is.null(blend$IS)) {
    ISblend <- sum(blend$IS * blend$Ratio)
  }

  blend <- blend %>%
    mutate(FW_Alk = Ratio * Alk0) %>%
    pmap_dfr(calculate_carbonate) %>%
    mutate(FW_CTCO3 = Ratio * iniCTCO3,
           FW_Temp = Ratio * Temp0,
           FW_Na = Ratio * Na.Conc,
           FW_Ca = Ratio * Ca.Conc,
           FW_Mg = Ratio * Mg.Conc,
           FW_K = Ratio * K.Conc,
           FW_Cl = Ratio * Cl.Conc,
           FW_SO4 = Ratio * SO4.Conc)

  # Unit conversions
  TempK <- sum(blend$FW_Temp) + 273.15 # C to K
  Alk <- sum(blend$FW_Alk)
  AlkMol <- convert(Alk, "CaCO3")

  FW_Na <- sum(blend$FW_Na)
  FW_Ca <- sum(blend$FW_Ca)
  FW_Mg <- sum(blend$FW_Mg)
  FW_Cl <- sum(blend$FW_Cl)
  FW_SO4 <- sum(blend$FW_SO4)
  FW_K <- sum(blend$FW_K)

  # Determine ionic strength (depends on which info is given)
  # First, determine from ions if possible
  if(!is.na(FW_Na) & !is.na(FW_Ca) & !is.na(FW_Mg) & !is.na(FW_Cl) & !is.na(FW_SO4)) {
    IS <- calculate_is(FW_Na, FW_Ca, FW_Mg, FW_K, FW_Cl, FW_SO4, AlkMol)
    # Next, determine from blended cond or TDS if possible
  } else if (exists("Cond")) {
    IS <- correlate_is_cond(Cond)
  } else if (exists("TDS")) {
    IS <- correlate_is_tds(TDS)
    # Finally, determine from blended IS
  } else if (exists("ISblend")) {
    IS <- ISblend
    # If no other params are provided, assume conductivity of 500 and calculate
  } else {
    IS <- correlate_is_cond(500)
  }

  # Acidity constants corrected for ionic strength
  cKw <- acidity_const(TempK, IS, "cKw")
  cK1 <- acidity_const(TempK, IS, "cK1")
  cK2 <- acidity_const(TempK, IS, "cK2")

  list(Alk0 = sum(blend$FW_Alk), CTCO3 = sum(blend$FW_CTCO3), IS = IS, TempK = TempK,
       Na.Conc = FW_Na, Ca.Conc = FW_Ca, Mg.Conc = FW_Mg, K.Conc = FW_K, Cl.Conc = FW_Cl, SO4.Conc = FW_SO4,
       cKw = cKw, cK1 = cK1, cK2 = cK2)

}
########################################################################################################################*

# * CORROSIVITY PARAMETERS ----

calculate_corrosivity <- function(pH, Temp, IS, HCO3, CO3, CTCO3, Alk,
                                  Na.Conc, Cl.Conc, Ca.Conc, Mg.Conc, SO4.Conc, K.Conc, ...){
  #Unit conversions
  TempK <- Temp + 273.15 # C to K
  CTCa <- convert(Ca.Conc, "Ca")

  LSI <- pH -
    (acidity_const(TempK, IS, "K2") - acidity_const(TempK, IS, "KCaCO3") -
       log10(CTCa) - log10(HCO3) -
       log10(acidity_const(TempK, IS, "Mono")) - log10(acidity_const(TempK, IS, "Di")))

  DeltaG <- -2.3 * 1.99E-3 * TempK * LSI

  DIC <- (CTCO3) * 12.0107 * 1000

  CSMR <- Cl.Conc / SO4.Conc

  LarsonSkold <- (convert(Cl.Conc, "Cl", end = "eq", start = "mg") *1000 + convert(SO4.Conc, "SO4", end = "eq", start = "mg") *1000) /
    (convert(HCO3, "HCO3", end = "eq", start = "mol") *1000 + convert(CO3, "CO3", end = "eq", start = "mol") *1000)

  AI <- pH + log10(convert(Ca.Conc, "Ca", end = "mgCaCO3", start = "mg")) + log10(Alk)

  RI <- 2 * (acidity_const(TempK, IS, "K2") - acidity_const(TempK, IS, "KCaCO3") -
               log10(CTCa) - log10(HCO3) -
               log10(acidity_const(TempK, IS, "Mono")) - log10(acidity_const(TempK, IS, "Di"))) - pH

  return(list(..., LSI = LSI, DeltaG = DeltaG, DIC = DIC, CSMR = CSMR, LarsonSkold = LarsonSkold, AI = AI, RI = RI))

}

########################################################################################################################*
########################################################################################################################*
#----------------------------------------------------------------------------------------------------------------------#*
# CHAINED FUNCTIONS ----
#----------------------------------------------------------------------------------------------------------------------#*

#' Calculate pH, Alk, and optional corrosivity parameters given raw water quality and optional chemical addition
#'
#' Given raw pH, Alk, Temp and chem doses, find pH and alk (corrosivity params optional). Works well in a dataframe,
#' but this will require/create a lot of columns.
#'
#' @param Temp0 temperature in degrees C
#' @param Alk0 alkalinity in mg CaCO3/L
#' @param pH0 pH in SU
#' @param Cond conductivity (optional, use when ionic strength and ion concentrations are unknown)
#' @param TDS total dissolved solids (optional, use when IS and ion concentrations are unknown)
#' @param IS ionic strength (optional, enter if known)
#' @param Temp.Final final temperature (optional, enter if different from initial)
#' @param HCl.Conc Hydrochloric Acid | HCl	| HCl -> H+ + Cl- (all chemical doses default to 0, units are mg/L)
#' @param H2SO4.Conc Sulfuric Acid | H2SO4	| H2SO4  -> 2H+ + SO42-
#' @param H3PO4.Conc Phosphoric Acid | H3PO4	| H3PO4  -> 3H+ + PO43-
#' @param CO2.Gas Carbon Dioxide | CO2 (gas) | CO2 + H2O -> H2CO3*
#' @param NaOH.Conc sodium hydroxide/Caustic Soda  | NaOH	| NaOH -> Na+ + OH-
#' @param Na2CO3.Conc Soda ash  | Na2CO3	| Na2CO3  -> 2Na+ + CO32-
#' @param NaHCO3.Conc Sodium Bicarbonate | NaHCO3	| NaHCO3  -> Na+ + H+ + CO32-
#' @param Lime.Conc Lime | Ca(OH)2	| Ca(OH)2 -> Ca2+ + 2OH-
#' @param MgOH2.Conc Magnesium Hydroxide | Mg(OH)2	| Mg(OH)2 -> Mg2+ + 2OH-
#' @param Cl2.Gas Chlorine Gas  | Cl2 (gas)	| Cl2 + H2O -> HOCl + H+ + Cl-
#' @param Cl2.Conc Sodium Hypochlorite | Cl2 | NaOCl	NaOCl -> Na+ + OCl-
#' @param CaOCl2.Conc Calcium Hypochlorite (HTH) | Cl2  | Ca(OCl)2	Ca(OCl)2 -> Ca2+ + 2OCl-
#' @param Alum.Conc Alum  | Al2(SO4)3*14H2O	| Al2(SO4)3*14H2O + 6HCO3- -> 2Al(OH)3(am) +3SO42- + 14H2O + 6CO2
#' @param FerricSO4.Conc Ferric Sulfate | Fe2(SO4)3	| Fe2(SO4)3 + 6HCO3- -> 2Fe(OH)3(am) +3SO42- + 6CO2
#' @param FerricCl.Conc Ferric Chloride | FeCl3	| FeCl3 +3HCO3- -> Fe(OH)3(am) + 3Cl- + 3CO2
#' @param CaCO3.Conc Calcium Carbonate | CaCO3 | CaCO3 -> Ca2+ + CO32-
#'
#' @param Na.Conc sodium concentration (optional, enter if known)
#' @param Ca.Conc calcium concentration (optional, enter if known)
#' @param Mg.Conc magnesium concentration (optional, enter if known)
#' @param K.Conc potassium concentration (optional, defaults to 0)
#' @param Cl.Conc chloride concentration (optional, enter if known)
#' @param SO4.Conc sulfate concentration (optional, enter if known)
#'
#' @param ionunits if Na, Ca, Mg, K, Cl, and SO4 are not entered in mg/L, provide the units here: "mol", "eq", or "mgCaCO3" (per L)
#' @param corrosivity if TRUE, corrosivity parameters (LSI, Delta G, DIC, CSMR, Larson SKold index, AI, and RI) will be included in output.
#'
#' @examples chem_calc_ph(20, 120, 8.1, Alum.Conc = 40)
#' chem_calc_ph(20, 120, 8.1, Cond = 350, H2SO4.Conc = 8, FerricSO4.Conc = 20, corrosivity = TRUE)
#'
#' @export
#'
chem_calc_ph <- function(Temp0, Alk0, pH0, Cond = 500, TDS = NA, IS = NA, Temp.Final = Temp0,
                         HCl.Conc = 0, H2SO4.Conc = 0, H3PO4.Conc = 0, CO2.Gas = 0,
                         NaOH.Conc = 0, Na2CO3.Conc = 0, NaHCO3.Conc = 0, Lime.Conc = 0, MgOH2.Conc = 0,
                         Cl2.Gas = 0, Cl2.Conc = 0, CaOCl2.Conc = 0,
                         Alum.Conc = 0, FerricCl.Conc = 0, FerricSO4.Conc = 0, CaCO3.Conc = 0,
                         Na.Conc = NA, Ca.Conc = NA, Mg.Conc = NA, K.Conc = 0, Cl.Conc = NA, SO4.Conc = NA,
                         ionunits = "mg", corrosivity = FALSE, ...) {

  # Convert Units if necessary
  Na.Conc <- convert(Na.Conc, "Na", end = "mg", start = ionunits)
  Ca.Conc <- convert(Ca.Conc, "Ca", end = "mg", start = ionunits)
  Mg.Conc <- convert(Mg.Conc, "Mg", end = "mg", start = ionunits)
  K.Conc <- convert(K.Conc, "K", end = "mg", start = ionunits)
  Cl.Conc <- convert(Cl.Conc, "Cl", end = "mg", start = ionunits)
  SO4.Conc <- convert(SO4.Conc, "SO4", end = "mg", start = ionunits)

  Raw.Carb <- calculate_carbonate(Temp0, Alk0, pH0, Cond, TDS, IS, Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc)


  Chem.Added <- chem_addition(Raw.Carb$iniCTCO3, Raw.Carb$iniCBA, Raw.Carb$IS, Temp.Final,
                              HCl.Conc, H2SO4.Conc, H3PO4.Conc, CO2.Gas,
                              NaOH.Conc, Na2CO3.Conc, NaHCO3.Conc, Lime.Conc, MgOH2.Conc ,
                              Cl2.Gas, Cl2.Conc, CaOCl2.Conc,
                              Alum.Conc, FerricCl.Conc, FerricSO4.Conc, CaCO3.Conc,
                              Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc)

  pH <- calculate_alkph(Chem.Added)$pH
  Alk <- calculate_alkph(Chem.Added)$Alk
  IS <- Chem.Added$IS
  HCO3 <- calculate_alkph(Chem.Added)$HCO3
  CO3 <- calculate_alkph(Chem.Added)$CO3

  Na.Conc <- convert(Chem.Added$finCTNa, "Na", end = "mg", start = "mol")
  Cl.Conc <- convert(Chem.Added$finCTCl, "Cl", end = "mg", start = "mol")
  Ca.Conc <- convert(Chem.Added$finCTCa, "Ca", end = "mg", start = "mol")
  Mg.Conc <- convert(Chem.Added$finCTMg, "Mg", end = "mg", start = "mol")
  SO4.Conc <- convert(Chem.Added$finCTSO4, "SO4", end = "mg", start = "mol")

  if (corrosivity) {
    Corros.Params <- calculate_corrosivity(pH, Temp.Final, Chem.Added$IS, HCO3, CO3, Chem.Added$CTCO3, Alk,
                                           Na.Conc, Cl.Conc, Ca.Conc, Mg.Conc, SO4.Conc, K.Conc)
    return(list(..., pH = pH, Alkalinity = Alk, IS = IS,
                LSI = Corros.Params$LSI, DeltaG = Corros.Params$DeltaG, DIC = Corros.Params$DIC, CSMR = Corros.Params$CSMR,
                LarsonSkold = Corros.Params$LarsonSkold, AI = Corros.Params$AI, RI = Corros.Params$RI))
  } else {
    return(list(..., pH = pH, Alkalinity = Alk, IS = IS,
                Na.Conc = Na.Conc, Cl.Conc = Cl.Conc,
                Ca.Conc = Ca.Conc, Mg.Conc = Mg.Conc,
                K.Conc = K.Conc, SO4.Conc = SO4.Conc))
  }
}

#' Calculate blended pH, Alk, temp, and optional corrosivity parameters given  and optional chemical addition
#'
#' Given raw pH, Alk, Temp and chem doses for each water and blending ratios, find pH and alk (corrosivity params optional).
#'
#' @param blend Dataframe where rows are the different waters with the following columns:
#'   Ratio (fraction of blend for each water, should be between 0 and 1),
#'   pH0 (pH of each water), Alk0 (alkalinity in mg CaCO3/L), Temp0 (temperature in deg C), Cond (conductivity, optional),
#'   Ionic Strength (optional),
#'   Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc (ion concentrations, optional)
#' @param HCl.Conc Hydrochloric Acid | HCl	| HCl -> H+ + Cl- (all chemical doses default to 0, units are mg/L)
#' @param H2SO4.Conc Sulfuric Acid | H2SO4	| H2SO4  -> 2H+ + SO42-
#' @param H3PO4.Conc Phosphoric Acid | H3PO4	| H3PO4  -> 3H+ + PO43-
#' @param CO2.Gas Carbon Dioxide | CO2 (gas) | CO2 + H2O -> H2CO3*
#' @param NaOH.Conc sodium hydroxide/Caustic Soda  | NaOH	| NaOH -> Na+ + OH-
#' @param Na2CO3.Conc Soda ash  | Na2CO3	| Na2CO3  -> 2Na+ + CO32-
#' @param NaHCO3.Conc Sodium Bicarbonate | NaHCO3	| NaHCO3  -> Na+ + H+ + CO32-
#' @param Lime.Conc Lime | Ca(OH)2	| Ca(OH)2 -> Ca2+ + 2OH-
#' @param MgOH2.Conc Magnesium Hydroxide | Mg(OH)2	| Mg(OH)2 -> Mg2+ + 2OH-
#' @param Cl2.Gas Chlorine Gas  | Cl2 (gas)	| Cl2 + H2O -> HOCl + H+ + Cl-
#' @param Cl2.Conc Sodium Hypochlorite | Cl2 | NaOCl	NaOCl -> Na+ + OCl-
#' @param CaOCl2.Conc Calcium Hypochlorite (HTH) | Cl2  | Ca(OCl)2	Ca(OCl)2 -> Ca2+ + 2OCl-
#' @param Alum.Conc Alum  | Al2(SO4)3x14H2O | Al2(SO4)3x14H2O + 6HCO3- -> 2Al(OH)3(am) +3SO42- + 14H2O + 6CO2
#' @param FerricSO4.Conc Ferric Sulfate | Fe2(SO4)3	| Fe2(SO4)3 + 6HCO3- -> 2Fe(OH)3(am) +3SO42- + 6CO2
#' @param FerricCl.Conc Ferric Chloride | FeCl3	| FeCl3 +3HCO3- -> Fe(OH)3(am) + 3Cl- + 3CO2
#' @param CaCO3.Conc Calcium Carbonate | CaCO3 | CaCO3 -> Ca2+ + CO32-
#'
#' @param ionunits if Na, Ca, Mg, K, Cl, and SO4 are not entered in mg/L, provide the units here: "mol", "eq", or "mgCaCO3" (per L)
#' @param corrosivity if TRUE, corrosivity parameters (LSI, Delta G, DIC, CSMR, Larson SKold index, AI, and RI) will be included in output.
#' @param df.return argument still in development. Returns outputs as a dataframe if TRUE
#'
#' @examples blend_calc_ph(blend_df, Alum.Conc = 40, corrosivity = TRUE)
#'
#' @export
#'
blend_calc_ph <- function(blend, HCl.Conc = 0, H2SO4.Conc = 0, H3PO4.Conc = 0, CO2.Gas = 0,
                          NaOH.Conc = 0, Na2CO3.Conc = 0, NaHCO3.Conc = 0, Lime.Conc = 0, MgOH2.Conc = 0,
                          Cl2.Gas = 0, Cl2.Conc = 0, CaOCl2.Conc = 0,
                          Alum.Conc = 0, FerricCl.Conc = 0, FerricSO4.Conc = 0, CaCO3.Conc = 0,
                          ionunits = "mg", corrosivity = FALSE, df.return = FALSE, ...) {

  Blended <- calculate_blend_raw(blend)
  Temp.Blend <- Blended$TempK - 273.15 # K to C

  pH.blend <- blend_ph(Blended$CTCO3, Blended$Alk0, Blended$cK1, Blended$cK2, Blended$cKw)

  # Convert Units if necessary
  Na.Conc <- convert(Blended$Na.Conc, "Na", end = "mg", start = ionunits)
  Ca.Conc <- convert(Blended$Ca.Conc, "Ca", end = "mg", start = ionunits)
  Mg.Conc <- convert(Blended$Mg.Conc, "Mg", end = "mg", start = ionunits)
  K.Conc <- convert(Blended$K.Conc, "K", end = "mg", start = ionunits)
  Cl.Conc <- convert(Blended$Cl.Conc, "Cl", end = "mg", start = ionunits)
  SO4.Conc <- convert(Blended$SO4.Conc, "SO4", end = "mg", start = ionunits)

  Blend.Carb <- calculate_carbonate(Temp.Blend, Blended$Alk0, pH.blend, IS = Blended$IS,
                                 Na.Conc = Blended$Na.Conc, Ca.Conc = Blended$Ca.Conc, Mg.Conc = Blended$Mg.Conc,
                                 K.Conc = Blended$K.Conc, Cl.Conc = Blended$Cl.Conc, SO4.Conc = Blended$SO4.Conc)

  Blend.Chem.Added <- chem_addition(Blend.Carb$iniCTCO3, Blend.Carb$iniCBA, Blend.Carb$IS, Temp.Blend,
                                    HCl.Conc, H2SO4.Conc, H3PO4.Conc, CO2.Gas,
                                    NaOH.Conc, Na2CO3.Conc, NaHCO3.Conc, Lime.Conc, MgOH2.Conc ,
                                    Cl2.Gas, Cl2.Conc, CaOCl2.Conc,
                                    Alum.Conc, FerricCl.Conc, FerricSO4.Conc, CaCO3.Conc,
                                    Na.Conc = Blend.Carb$Na.Conc, Ca.Conc = Blend.Carb$Ca.Conc, Mg.Conc = Blend.Carb$Mg.Conc,
                                    K.Conc = Blend.Carb$K.Conc, Cl.Conc = Blend.Carb$Cl.Conc, SO4.Conc = Blend.Carb$SO4.Conc)

  pH <- calculate_alkph(Blend.Chem.Added)$pH
  Alk <- calculate_alkph(Blend.Chem.Added)$Alk
  HCO3 <- calculate_alkph(Blend.Chem.Added)$HCO3
  CO3 <- calculate_alkph(Blend.Chem.Added)$CO3

  Na.Conc <- convert(Blend.Chem.Added$finCTNa, "Na", end = "mg", start = "mol")
  Cl.Conc <- convert(Blend.Chem.Added$finCTCl, "Cl", end = "mg", start = "mol")
  Ca.Conc <- convert(Blend.Chem.Added$finCTCa, "Ca", end = "mg", start = "mol")
  Mg.Conc <- convert(Blend.Chem.Added$finCTMg, "Mg", end = "mg", start = "mol")
  SO4.Conc <- convert(Blend.Chem.Added$finCTSO4, "SO4", end = "mg", start = "mol")

  if (corrosivity) {
    Corros.Params <- calculate_corrosivity(pH, Temp.Blend, IS = Blend.Chem.Added$IS,
                                           HCO3, CO3, CTCO3 = Blend.Chem.Added$CTCO3,Alk,
                                           Na.Conc, Cl.Conc, Ca.Conc, Mg.Conc, SO4.Conc, K.Conc)
    if (df.return) {
      return(data.frame(pH.fin = pH, Alk.fin = Alk, Temp.fin = Temp.Blend,
                 LSI = Corros.Params$LSI, DeltaG = Corros.Params$DeltaG, DIC = Corros.Params$DIC, CSMR = Corros.Params$CSMR,
                 LarsonSkold = Corros.Params$LarsonSkold, AI = Corros.Params$AI, RI = Corros.Params$RI))
    } else {
      return(list(..., pH = pH, Alkalinity = Alk, Temperature = Temp.Blend,
                  LSI = Corros.Params$LSI, DeltaG = Corros.Params$DeltaG, DIC = Corros.Params$DIC, CSMR = Corros.Params$CSMR,
                  LarsonSkold = Corros.Params$LarsonSkold, AI = Corros.Params$AI, RI = Corros.Params$RI))
    }
  } else {
    if (df.return) {
      return(data.frame(pH.fin = pH, Alk.fin = Alk, Temp.fin = Temp.Blend,
                        IS = Blend.Chem.Added$IS, Na.Conc = Na.Conc, Cl.Conc = Cl.Conc,
                        Ca.Conc = Ca.Conc, Mg.Conc = Blend.Carb$Mg.Conc,
                        K.Conc = Blend.Carb$K.Conc, SO4.Conc = SO4.Conc))
    } else {
      return(list(..., pH = pH, Alkalinity = Alk, Temperature = Temp.Blend))
    }
  }
}


#' Calculate chemical dose based on target pH
#'
#' Given pH, alk, a target pH, and chemical name, solve for target dose
#'
#' @param Temp0 temperature in degrees C
#' @param Alk0 alkalinity in mg CaCO3/L
#' @param pH0 starting pH in SU
#' @param pH.Final desired final pH in SU
#' @param Cond conductivity (optional, use when ionic strength and ion concentrations are unknown)
#' @param TDS total dissolved solids (optional, use when ionic strength and ion concentrations are unknown)
#' @param IS ionic strength (optional, enter if known)
#'
#' @param Na.Conc sodium concentration (optional, enter if known)
#' @param Ca.Conc calcium concentration (optional, enter if known)
#' @param Mg.Conc magnesium concentration (optional, enter if known)
#' @param K.Conc potassium concentration (optional, defaults to 0)
#' @param Cl.Conc chloride concentration (optional, enter if known)
#' @param SO4.Conc sulfate concentration (optional, enter if known)
#'
#' @param chemical enter name of chemical to be solved for in quotes, caustic, lime, hydrochloric, CO2, or sulfuric
#' @param ionunits if Na, Ca, Mg, K, Cl, and SO4 are not entered in mg/L, provide the units here: "mol", "eq", or "mgCaCO3" (per L)
#'
#' @examples solve_chem_dose(20, 120, 6.8, 8.0)
#'
#' @export
#'
solve_chem_dose <- function(Temp0, Alk0, pH0, pH.Final, Cond = 500, TDS = NA, IS = NA,
                          Na.Conc = NA, Ca.Conc = NA, Mg.Conc = NA, K.Conc = NA, Cl.Conc = NA, SO4.Conc = NA,
                          chemical, ionunits = "mg", ...) {

  # Convert Units if necessary
  Na.Conc <- convert(Na.Conc, "Na", end = "mg", start = ionunits)
  Ca.Conc <- convert(Ca.Conc, "Ca", end = "mg", start = ionunits)
  Mg.Conc <- convert(Mg.Conc, "Mg", end = "mg", start = ionunits)
  K.Conc <- convert(K.Conc, "K", end = "mg", start = ionunits)
  Cl.Conc <- convert(Cl.Conc, "Cl", end = "mg", start = ionunits)
  SO4.Conc <- convert(SO4.Conc, "SO4", end = "mg", start = ionunits)

  # Check that you are using a correct chemical.
  if (!(chemical == "caustic" | chemical == "lime" |
        chemical == "sulfuric" | chemical == "CO2" | chemical == "hydrochloric")) {
    stop("Chemical name not recognized. Check documentation.")
  }

  delta_pH <- function(targetdose, chemical, Temp0, Alk0, pH0, pH.Final, Cond, TDS, IS,
                       Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc){

    # Assign chemical doses
    NaOH.Conc <- ifelse(chemical == "caustic", targetdose, 0)
    Lime.Conc <- ifelse(chemical == "lime", targetdose, 0)
    H2SO4.Conc <- ifelse(chemical == "sulfuric", targetdose, 0)
    CO2.Gas <- ifelse(chemical == "CO2", targetdose, 0)
    HCl.Conc <- ifelse(chemical == "hydrochloric", targetdose, 0)

    (pH.Final - chem_calc_ph(Temp0, Alk0, pH0, Cond, TDS, IS,
                             NaOH.Conc = NaOH.Conc, Lime.Conc = Lime.Conc,
                             H2SO4.Conc = H2SO4.Conc, CO2.Gas = CO2.Gas, HCl.Conc = HCl.Conc,
                             Na.Conc = Na.Conc, Ca.Conc = Ca.Conc, Mg.Conc = Mg.Conc,
                             K.Conc = K.Conc, Cl.Conc = Cl.Conc, SO4.Conc = SO4.Conc)$pH)^2
  }

  if (((chemical == "caustic" | chemical == "lime") & pH.Final > pH0) |
      ((chemical == "sulfuric" | chemical == "hydrochloric") & pH.Final < pH0) |
      (chemical == "CO2" & pH.Final < pH0 & pH.Final >= 6.5)) {
    chemicaldose <- optimize(delta_pH, interval = c(0,200),
                             chemical = chemical, Temp0 = Temp0, Alk0 = Alk0, pH0 = pH0, pH.Final = pH.Final,
                             Cond = Cond, IS = IS,
                             Na.Conc = Na.Conc, Ca.Conc = Ca.Conc, Mg.Conc = Mg.Conc,
                             K.Conc = K.Conc, Cl.Conc = Cl.Conc, SO4.Conc = SO4.Conc)$minimum
  } else {
    warning("Target pH cannot be met with selected chemical, NA produced.")
    chemicaldose <- NA
  }

  return(list(..., chem.solve = chemicaldose, chemical = chemical, Temp0 = Temp0, Alk0 = Alk0, pH0 = pH0, Cond = Cond, IS = IS,
              Ca.Conc = Ca.Conc, Cl.Conc = Cl.Conc, SO4.Conc = SO4.Conc))
}

#END ###################################################################################################################

