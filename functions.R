# Part 1: Data cleaning

# Converts OR and HR into RR, and leaves RR unconverted
# unit: unit of effect size measure, measured as odds ratio (OR), hazard ratio (HR), or relative risk (RR).
# ee: effect size, as column in table
# i: incidence rate (person per year)
rrConvert <- function(unit, ee, i) {case_when(
  unit == "OR" ~ ee/((1-i)+(i*ee)),
  unit == "RR" ~ ee
)
}

# Adjusts NO2, O3, PM2.5, PM10-2.5, PM10, SO2  (per 10 ug/m3)
# lee: log effect size 
# i: original increment of exposure (ug/m3)
adjEE_ug <- function(lee, i){
  lee * (10/i)
}

# Adjusts BC (per 1 ug/m3)
# lee: log effect size 
# i: original increment of exposure (ug/m3)
adjEE_BC_ug <- function(lee, i){
  lee * (0.1/i)
}

# Adjusts UFP (per 10,000 particles/cm3)
# lee: log effect size 
# i: original increment of exposure (particles/cm3)
adjEE_UFP_part <- function(lee, i){
  lee * (10000/i)
}

# Adjusts CO (per 0.1 mg/m3)
# lee: log effect size 
# i: original increment of exposure (mg/m3)
adjEE_CO_mg <- function(lee, i){
  lee * (0.1/i)
}

# Adjusts CO (per 1 mg/m3)
# lee: log effect size 
# i: original increment of exposure (ppm)
# avg_t: mean temperature
adjEE_CO_ppm <- function(lee, i, avg_t){
  i_standard <- (i * 28 * 273.15)/(0.08206 * (273.15 + avg_t)^2)
  lee * (0.1/i_standard)
}

# Adjusts NO2 (per 10 ug/m3)
# lee: log effect size 
# i: original increment of exposure (ppb)
adjEE_NO2_ppb <- function(lee, i, avg_t){
  i_standard <- (i * 46 * 273.15)/(0.08206 * (273.15 + avg_t)^2)
  lee * (10/i_standard)
}

# Adjusts SO2 (per 10 ug/m3)
# lee: log effect size 
# i: original increment of exposure (ppb)
adjEE_SO2_ppb <- function(lee, i, avg_t){
  i_standard <- (i * 64 * 273.15)/(0.08206 * (273.15 + avg_t)^2)
  lee * (10/i_standard)
}

# Adjusts O3 (per 10 ug/m3)
# lee: log effect size 
# i: original increment of exposure (ppb)
adjEE_O3_ppb <- function(lee, i, avg_t){
  i_standard <- (i * 48 * 273.15)/(0.08206 * (273.15 + avg_t)^2)
  lee * (10/i_standard)
}

# Adjusts non-optimal temperature (per 5 deg C increase/decrease)
# lee: log effect size 
# i: original increment of exposure (deg C)
adjEE_temp <- function(lee, i){
  lee * (5/i)
}

# Meta analysis conversions, adds 8 columns for 
# df: meta-analysis data source
# ee: effect size (column ee_original)
# ll: effect size lower 95% CI (column ll_original)
# ul: effect size upper 95% CI (column ll_original)
# oi: study increment (column increase_original)
# ti: target increment

metaPrep <- function(df, ee, ll, ul, oi, ti){
    mutate(log.EE = log(ee)) %>%
    mutate(log.SE = (log(ul)-log(ll))/3.92) %>%
    mutate(alog.EE = adjInc(log.EE, ic, ti)) %>%
    mutate(alog.SE = adjInc(log.SE, ic, ti))
  }

# Part 2: Meta-analysis functions
# Fixed-effect meta-analysis for overlapping lags
# df: data
# ti: title for meta-analysis, input as string

meta_analyze_lag <- function(df){ 
  metagen(
    data = df,
    TE = alog.EE,
    seTE = alog.SE,
    sm = "OR",
    common = TRUE,
    random = FALSE,
    method.tau = "REML",
    method.common.ci = "classic",
    studlab = fa_year,
    prediction = TRUE,
    digits = 4
  )
}

# Fixed-effect meta-analysis to obtain overall effects
# df: data
# ti: title for meta-analysis, input as string

meta_analyze_overall <- function(df){ 
  metagen(
    data = df,
    TE = alog.EE,
    seTE = alog.SE,
    sm = "OR",
    common = TRUE,
    random = FALSE,
    method.tau = "REML",
    method.common.ci = "classic",
    studlab = fa_year,
    digits = 4
  )
}

# Pooling of effect sizes (fully adjusted), standardized to exposure increment
# df: data
# ti: title for meta-analysis, input as string

meta_analyze <- function(df, ti){ 
  metagen(
    data = df,
    TE = alog.EE,
    seTE = alog.SE,
    sm = "OR",
    common = TRUE,
    random = TRUE,
    method.tau = "REML",
    method.random.ci = "HK",
    method.common.ci = "classic",
    studlab = fa_year,
    title = ti,
    prediction = TRUE,
    digits = 4
  )
}

# Pooling of effect sizes (fully adjusted), standardized to exposure increment
# df: data
# ti: title for meta-analysis, input as string

meta_analyze_DL <- function(df, ti){ 
  metagen(
    data = df,
    TE = alog.EE,
    seTE = alog.SE,
    sm = "OR",
    common = FALSE,
    random = TRUE,
    method.tau = "REML",
    studlab = fa_year,
    title = ti,
    prediction = TRUE,
    digits = 4
  )
}

# Funnel plot
# ma: input meta-analysis object (from metagen)
# sl: include study label (true/false)
fun <- function(ma, ti){
  funnel(ma, 
         main = ti,
         pch=19, 
         contour = c(0.9, 0.95, 0.99),
         col.contour = c("gray75", "gray85", "gray95"),
         cex = 0.8,
         ref = 1,
         xlab = "Relative Risk")
}

# Egger's test for funnel plot bias
# ma: input meta-analysis object (from metagen)
funBias <- function(ma){
      metabias(ma, method.bias = "Egger", k.min = 10)
      }


# Print forest plot for manuscript
# ma: input meta-analysis object (from metagen)
# fix: include pooled fixed effect size estimate (true/false)
# rand: include pooled random effect size estimate (true/false)
# forest_basic <- function(ma, fix, rand) {
#   forest(ma,
#           prediction = F,
#           print.tau2 = T,
#           print.chi2 = T,
#           common = fix,
#           random = rand,
#           leftcols = c("fa_year"),
#           leftlabs = c ("Study"),
#           label.left = "No effect",
#           label.right = "Effect",
#           col.square = "grey20"
#          )
# }

# Print forest plot with RevMan5 format for manuscript
# ma: input meta-analysis object (from metagen)
# xl: lower bound for scale bar on forest plot
# xu: upper bound for scale bar on forest plot
forest_rev5 <- function(ma, xl, xu) {
  forest(ma,
         prediction = ma$prediction,
         print.tau2 = T,
         xlim = c(xl, xu),
         common = FALSE,
         random = TRUE,
         digits = 4,
         fontsize = 10,
         layout = "RevMan5",
         label.left = "Decreased risk",
         label.right = "Increased risk"
  )
}

# Forest lancet
forest_lancet <- function(ma, xl, xu) {
  
  forest(ma,
         prediction = ma$prediction,
         xlim = c(xl, xu),
         layout = "RevMan5",
         label.left = "Decreased risk",
         label.right = "Increased risk",
         col.square = "grey50",
         col.square.lines = "black",
         col.lines = "black",
         col.inside = "black",
         digits = 4,
         fontsize = 10,
         family = "serif",
         leftcols = c("studlab", "w.random","effect.ci"),
         weight = "random"
         )
}

draw_forest <- function(ma, ti) {
  
  forest(ma,
         prediction = ma$prediction,
         main = ti,
         layout = "RevMan5",
         label.left = "Decreased risk",
         label.right = "Increased risk",
         col.square = "grey50",
         col.square.lines = "black",
         col.lines = "black",
         col.inside = "black",
         digits = 4,
         fontsize = 10,
         family = "serif",
         leftcols = c("studlab", "w.random","effect.ci")
  )
}


# Exploratory meta-regression for individual study co-variates
# ma: input meta-analysis object (from metagen)
# Average exposure
mregExp <- function(ma){
  metareg(ma, ~mean_adjusted, rm.na = T)
}

# Continent
mregCon <- function(ma){
  metareg(ma, ~continent, rm.na = T)
}
# Sex/gender proportion
mregSex <- function(ma){
  metareg(ma, ~sex.pf, rm.na = T)
}
# Case definition of asthma
mregCase <- function(ma){
  metareg(ma, ~case.def, rm.na = T)
}
# Adjustment for co-pollutants
mregCopol <- function(ma){
  metareg(ma, ~adj.pol, rm.na = T)
}

# Meta-regression for all other variables
# ma: input meta-analysis object (from metagen)
mregAll <- function(ma){
  metareg(ma, ~avg.exp+ region + sex.pf + case.def + adj.pol, rm.na = T)
}

