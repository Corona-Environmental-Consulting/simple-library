% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Acid_Base.R
\name{blend_calc_ph}
\alias{blend_calc_ph}
\title{Calculate blended pH, Alk, temp, and optional corrosivity parameters given  and optional chemical addition}
\usage{
blend_calc_ph(
  blend,
  HCl.Conc = 0,
  H2SO4.Conc = 0,
  H3PO4.Conc = 0,
  CO2.Gas = 0,
  NaOH.Conc = 0,
  Na2CO3.Conc = 0,
  NaHCO3.Conc = 0,
  Lime.Conc = 0,
  MgOH2.Conc = 0,
  Cl2.Gas = 0,
  Cl2.Conc = 0,
  CaOCl2.Conc = 0,
  Alum.Conc = 0,
  FerricCl.Conc = 0,
  FerricSO4.Conc = 0,
  CaCO3.Conc = 0,
  ionunits = "mg",
  corrosivity = FALSE,
  df.return = FALSE,
  ...
)
}
\arguments{
\item{blend}{Dataframe where rows are the different waters with the following columns:
Ratio (fraction of blend for each water, should be between 0 and 1),
pH0 (pH of each water), Alk0 (alkalinity in mg CaCO3/L), Temp0 (temperature in deg C), Cond (conductivity, optional),
Ionic Strength (optional),
Na.Conc, Ca.Conc, Mg.Conc, K.Conc, Cl.Conc, SO4.Conc (ion concentrations, optional)}

\item{HCl.Conc}{Hydrochloric Acid | HCl    | HCl -> H+ + Cl- (all chemical doses default to 0, units are mg/L)}

\item{H2SO4.Conc}{Sulfuric Acid | H2SO4    | H2SO4  -> 2H+ + SO42-}

\item{H3PO4.Conc}{Phosphoric Acid | H3PO4    | H3PO4  -> 3H+ + PO43-}

\item{CO2.Gas}{Carbon Dioxide | CO2 (gas) | CO2 + H2O -> H2CO3*}

\item{NaOH.Conc}{sodium hydroxide/Caustic Soda  | NaOH    | NaOH -> Na+ + OH-}

\item{Na2CO3.Conc}{Soda ash  | Na2CO3    | Na2CO3  -> 2Na+ + CO32-}

\item{NaHCO3.Conc}{Sodium Bicarbonate | NaHCO3    | NaHCO3  -> Na+ + H+ + CO32-}

\item{Lime.Conc}{Lime | Ca(OH)2    | Ca(OH)2 -> Ca2+ + 2OH-}

\item{MgOH2.Conc}{Magnesium Hydroxide | Mg(OH)2    | Mg(OH)2 -> Mg2+ + 2OH-}

\item{Cl2.Gas}{Chlorine Gas  | Cl2 (gas)    | Cl2 + H2O -> HOCl + H+ + Cl-}

\item{Cl2.Conc}{Sodium Hypochlorite | Cl2 | NaOCl    NaOCl -> Na+ + OCl-}

\item{CaOCl2.Conc}{Calcium Hypochlorite (HTH) | Cl2  | Ca(OCl)2    Ca(OCl)2 -> Ca2+ + 2OCl-}

\item{Alum.Conc}{Alum  | Al2(SO4)3x14H2O | Al2(SO4)3x14H2O + 6HCO3- -> 2Al(OH)3(am) +3SO42- + 14H2O + 6CO2}

\item{FerricCl.Conc}{Ferric Chloride | FeCl3    | FeCl3 +3HCO3- -> Fe(OH)3(am) + 3Cl- + 3CO2}

\item{FerricSO4.Conc}{Ferric Sulfate | Fe2(SO4)3    | Fe2(SO4)3 + 6HCO3- -> 2Fe(OH)3(am) +3SO42- + 6CO2}

\item{CaCO3.Conc}{Calcium Carbonate | CaCO3 | CaCO3 -> Ca2+ + CO32-}

\item{ionunits}{if Na, Ca, Mg, K, Cl, and SO4 are not entered in mg/L, provide the units here: "mol", "eq", or "mgCaCO3" (per L)}

\item{corrosivity}{if TRUE, corrosivity parameters (LSI, Delta G, DIC, CSMR, Larson SKold index, AI, and RI) will be included in output.}

\item{df.return}{argument still in development. Returns outputs as a dataframe if TRUE}
}
\description{
Given raw pH, Alk, Temp and chem doses for each water and blending ratios, find pH and alk (corrosivity params optional).
}
\examples{
blend_calc_ph(blend_df, Alum.Conc = 40, corrosivity = TRUE)

}
