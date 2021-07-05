ExcessDeathCOVID19
================
DHY
6/22/2021

## Excess death by COVID19

``` r
setwd("/Users/DHY 1/Documents/GitHub/COVID19impactModeling/")
load(file = "analysis_mod.Rda")

# Aggregating from quarters to the sum for the entire year

df.lin      <- analysis_mod %>% 
  group_by(Country, iso3, WHO_region, Nx, ncdr, sdi) %>%
  summarize(covid = sum(cov_deaths), 
            expected = sum(expected), 
            observed = sum(observed), 
            .groups = "drop") %>%
  ungroup() %>% 
  mutate(excess    = observed - expected, 
         expectedc = expected + covid, 
         ratio     = observed/expectedc, 
         covidr    = covid/Nx, 
         covsdi    = covidr * sdi) %>%
  arrange(iso3)

# Indices to identify rows
df.lin      <- df.lin %>% mutate(row = 1:nrow(df.lin))
ind.obsl    <- df.lin %>% filter(!is.na(excess)) %>% pull(row)


reg.mod     <- lm(log(ratio) ~ . - 1, data = df.lin %>% 
                    dplyr::select(ratio, sdi, covidr, covsdi, ncdr))
pe          <- reg.mod$coefficients  # point estimates
vc          <- vcov(reg.mod)         # variance-covariance

set.seed(1234)
draws       <- 4000

# Distribution of coefficients
simbetas    <- MASS::mvrnorm(draws, pe, vc)

# Input data to scale the coefficients
df.lin.res  <- df.lin %>% dplyr::select(dimnames(simbetas)[[2]]) %>% as.matrix()

# Distribution of fitted values (these are scaling factors from the ratio of observed to covid plus expected)
replic      <- exp(df.lin.res %*% t(simbetas))

# Use the scaling factors to obtain distribution of excess deaths
fitdist1e   <- replic*replicate(draws, df.lin$expectedc) - replicate(draws, df.lin$expected)
obsvals     <- df.lin$excess # vector of observed values

# For the non-missing, replace predictions with the actual observed
for (i in 1:draws){
  fitdist1e[ind.obsl,] <- obsvals[ind.obsl]
}
# Basic functions to summarise 80% bootstrap intervals or sample from gaussian

getsum <- function(x){
  tibble(mean = quantile(x, 0.5, na.rm = T), 
         lwr  = quantile(x, 0.025, na.rm = T), 
         uppr = quantile(x, 0.975, na.rm = T)) %>% round()
}

getsum2 <- function(x){
  c(quantile(x, 0.5, na.rm = T), 
    quantile(x, 0.025, na.rm = T), 
    quantile(x, 0.975, na.rm = T))
}

get.sample <- function(x){
  rnorm(draws, x[1], x[2])
}

# Table of summarized values joined to observed and explanator data

df.lin.out   <- data.table::data.table(t(apply(fitdist1e, 1, getsum2))) %>%
  rename(y = "50%", low = "2.5%", high = "97.5%") %>%
  cbind(df.lin)  %>% 
  mutate(y    = ifelse(!is.na(excess), excess, y),
         low  = ifelse(!is.na(excess), NA, low),
         high = ifelse(!is.na(excess), NA, high), 
         source = ifelse(is.na(excess), "Predicted", "Observed")) %>% arrange(-y) %>%
  data.frame()


# Indices by region ==> INDEX by country
df.lintemp<-df.lin%>% filter(WHO_region == "AMRO")
df.lintemp<-df.lintemp%>% mutate(row = 1:nrow(df.lintemp))
#df.lin      <- df.lin %>% mutate(row = 1:nrow(df.lin))
#==================== for missing (20) ===========================================
Argentina.indl   <- df.lintemp %>% filter(Country == "Argentina") %>% pull(row)
AntiguaandBarbudaa.indl   <- df.lintemp %>% filter(Country == "Antigua and Barbuda") %>% pull(row)
Bahamas.indl   <- df.lintemp %>% filter(Country == "Bahamas") %>% pull(row)
Belize.indl   <- df.lintemp %>% filter(Country == "Belize") %>% pull(row)
Barbados.indl   <- df.lintemp %>% filter(Country == "Barbados") %>% pull(row)
Dominica.indl   <- df.lintemp %>% filter(Country == "Dominica") %>% pull(row)
DominicanRepublic.indl   <- df.lintemp %>% filter(Country == "Dominican Republic") %>% pull(row)
Grenada.indl   <- df.lintemp %>% filter(Country == "Grenada") %>% pull(row)
Guyana.indl   <- df.lintemp %>% filter(Country == "Guyana") %>% pull(row)
Honduras.indl   <- df.lintemp %>% filter(Country == "Honduras") %>% pull(row)
Haiti.indl   <- df.lintemp %>% filter(Country == "Haiti") %>% pull(row)
SaintKittsandNevis.indl   <- df.lintemp %>% filter(Country == "Saint Kitts and Nevis") %>% pull(row)
SaintLucia.indl   <- df.lintemp %>% filter(Country == "Saint Lucia") %>% pull(row)
Nicaragua.indl   <- df.lintemp %>% filter(Country == "Nicaragua") %>% pull(row)
ElSalvador.indl   <- df.lintemp %>% filter(Country == "El Salvador ") %>% pull(row)
Suriname.indl   <- df.lintemp %>% filter(Country == "Suriname") %>% pull(row)
TrinidadandTobago.indl   <- df.lintemp %>% filter(Country == "Trinidad and Tobago") %>% pull(row)
Uruguay.indl   <- df.lintemp %>% filter(Country == "Uruguay") %>% pull(row)
SaintVincentandtheGrenadines.indl   <- df.lintemp %>% filter(Country == "Saint Vincent and the Grenadines") %>% pull(row)
VenezuelaBolivarianRepublicof.indl   <- df.lintemp %>% filter(Country == "Venezuela (Bolivarian Republic of)") %>% pull(row)
#==================== for non missing (15) ===========================================
BoliviaPlurinationalStateof.indl   <- df.lintemp %>% filter(Country == "Bolivia (Plurinational State of)") %>% pull(row)
Brazil.indl   <- df.lintemp %>% filter(Country == "Brazil") %>% pull(row)
Canada.indl   <- df.lintemp %>% filter(Country == "Canada") %>% pull(row)
Chile.indl   <- df.lintemp %>% filter(Country == "Chile") %>% pull(row)
Colombia.indl   <- df.lintemp %>% filter(Country == "Colombia") %>% pull(row)
CostaRica.indl   <- df.lintemp %>% filter(Country == "Costa Rica") %>% pull(row)
Cuba.indl   <- df.lintemp %>% filter(Country == "Cuba") %>% pull(row)
Ecuador.indl   <- df.lintemp %>% filter(Country == "Ecuador") %>% pull(row)
Guatemala.indl   <- df.lintemp %>% filter(Country == "Guatemala") %>% pull(row)
Jamaica.indl   <- df.lintemp %>% filter(Country == "Jamaica") %>% pull(row)
Mexico.indl   <- df.lintemp %>% filter(Country == "Mexico") %>% pull(row)
Panama.indl   <- df.lintemp %>% filter(Country == "Panama") %>% pull(row)
Peru.indl   <- df.lintemp %>% filter(Country == "Peru") %>% pull(row)
Paraguay.indl   <- df.lintemp %>% filter(Country == "Paraguay") %>% pull(row)
UnitedStatesofAmerica.indl   <- df.lintemp %>% filter(Country == "United States of America") %>% pull(row)

# ===========================================================================
Argentina.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Argentina.indl,]))), 2, sum)) %>%
                    mutate(Country = "Argentina",
                           covid = sum(df.lintemp[Argentina.indl,]$covid), 
                            Nx = sum(df.lintemp[Argentina.indl,]$Nx)) 

AntiguaandBarbudaa.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[AntiguaandBarbudaa.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Antigua and Barbuda",   covid = sum(df.lintemp[AntiguaandBarbudaa.indl ,]$covid), 
                     Nx = sum(df.lintemp[AntiguaandBarbudaa.indl ,]$Nx))

Bahamas.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Bahamas.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Bahamas",  covid = sum(df.lintemp[Bahamas.indl ,]$covid), 
                     Nx = sum(df.lintemp[Bahamas.indl ,]$Nx))

Belize.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Belize.indl ,]))), 2, sum)) %>% 
           mutate(Country = "Belize", covid = sum(df.lintemp[Belize.indl ,]$covid), 
                     Nx = sum(df.lintemp[Belize.indl ,]$Nx))

Barbados.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Barbados.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Barbados" ,covid = sum(df.lintemp[Barbados.indl ,]$covid), 
                     Nx = sum(df.lintemp[Barbados.indl ,]$Nx))

Dominica.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Dominica.indl ,]))), 2, sum)) %>% 
           mutate(Country = "Dominica", 
                  covid = sum(df.lintemp[Dominica.indl ,]$covid), 
                  Nx = sum(df.lintemp[Dominica.indl ,]$Nx))

DominicanRepublic.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[DominicanRepublic.indl  ,]))), 2, sum)) %>% 
           mutate( Country = "Dominican Republic"  ,  covid = sum(df.lintemp[DominicanRepublic.indl  ,]$covid), 
                     Nx = sum(df.lintemp[DominicanRepublic.indl  ,]$Nx))

Grenada.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Grenada.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Grenada" , covid = sum(df.lintemp[Grenada.indl ,]$covid), 
                     Nx = sum(df.lintemp[Grenada.indl ,]$Nx))

Guyana.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Guyana.indl  ,]))), 2, sum)) %>% 
           mutate( Country = "Guyana"   , covid = sum(df.lintemp[Guyana.indl  ,]$covid), 
                     Nx = sum(df.lintemp[Guyana.indl  ,]$Nx))

Honduras.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Honduras.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Honduras"  , covid = sum(df.lintemp[Honduras.indl ,]$covid), 
                     Nx = sum(df.lintemp[Honduras.indl ,]$Nx))

Haiti.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Haiti.indl,]))), 2, sum)) %>% 
           mutate( Country = "Haiti"   ,  covid = sum(df.lintemp[Haiti.indl,]$covid), 
                     Nx = sum(df.lintemp[Haiti.indl,]$Nx))

SaintKittsandNevis.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[SaintKittsandNevis.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Saint Kitts and Nevis"  ,   covid = sum(df.lintemp[SaintKittsandNevis.indl ,]$covid), 
                     Nx = sum(df.lintemp[SaintKittsandNevis.indl ,]$Nx))

SaintLucia.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[SaintLucia.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Saint Lucia"  , covid = sum(df.lintemp[SaintLucia.indl ,]$covid), 
                     Nx = sum(df.lintemp[SaintLucia.indl ,]$Nx))

Nicaragua.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Nicaragua.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Nicaragua"  ,  covid = sum(df.lintemp[Nicaragua.indl ,]$covid), 
                     Nx = sum(df.lintemp[Nicaragua.indl ,]$Nx))

ElSalvador.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[ElSalvador.indl ,]))), 2, sum)) %>% 
           mutate( Country = "El Salvador"   , covid = sum(df.lintemp[ElSalvador.indl ,]$covid), 
                     Nx = sum(df.lintemp[ElSalvador.indl ,]$Nx))

Suriname.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Suriname.indl,]))), 2, sum)) %>% 
           mutate( Country = "Suriname" ,    covid = sum(df.lintemp[Suriname.indl,]$covid), 
                     Nx = sum(df.lintemp[Suriname.indl,]$Nx))

TrinidadandTobago.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[TrinidadandTobago.indl  ,]))), 2, sum)) %>% 
           mutate( Country = "Trinidad and Tobago"  , covid = sum(df.lintemp[TrinidadandTobago.indl  ,]$covid), 
                     Nx = sum(df.lintemp[TrinidadandTobago.indl  ,]$Nx))

Uruguay.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Uruguay.indl,]))), 2, sum)) %>% 
           mutate( Country = "Uruguay"   ,  covid = sum(df.lintemp[Uruguay.indl,]$covid), 
                     Nx = sum(df.lintemp[Uruguay.indl,]$Nx))

SaintVincentandtheGrenadines.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[SaintVincentandtheGrenadines.indl  ,]))), 2, sum)) %>% 
           mutate( Country = "Saint Vincent and the Grenadines"  ,covid = sum(df.lintemp[SaintVincentandtheGrenadines.indl  ,]$covid), 
                     Nx = sum(df.lintemp[SaintVincentandtheGrenadines.indl  ,]$Nx))


VenezuelaBolivarianRepublicof.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[VenezuelaBolivarianRepublicof.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Venezuela (Bolivarian Republic of)" ,covid = sum(df.lintemp[VenezuelaBolivarianRepublicof.indl ,]$covid), 
                     Nx = sum(df.lintemp[VenezuelaBolivarianRepublicof.indl ,]$Nx))


# non missing 
BoliviaPlurinationalStateof.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[BoliviaPlurinationalStateof.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Bolivia (Plurinational State of)",  covid = sum(df.lintemp[BoliviaPlurinationalStateof.indl ,]$covid), 
                     Nx = sum(df.lintemp[BoliviaPlurinationalStateof.indl ,]$Nx))

Brazil.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Brazil.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Brazil" ,covid = sum(df.lintemp[Brazil.indl ,]$covid), 
                     Nx = sum(df.lintemp[Brazil.indl ,]$Nx))

Canada.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Canada.indl ,]))), 2, sum)) %>% 
           mutate(Country = "Canada" , covid = sum(df.lintemp[Canada.indl ,]$covid), 
                     Nx = sum(df.lintemp[Canada.indl ,]$Nx))

Chile.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Chile.indl ,]))), 2, sum)) %>% 
           mutate(Country = "Chile" , covid = sum(df.lintemp[Chile.indl ,]$covid), 
                     Nx = sum(df.lintemp[Chile.indl ,]$Nx))

Colombia.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Colombia.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Colombia"  , covid = sum(df.lintemp[Colombia.indl ,]$covid), 
                     Nx = sum(df.lintemp[Colombia.indl ,]$Nx))

CostaRica.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[CostaRica.indl ,]))), 2, sum)) %>% 
           mutate( Country = "Costa Rica"  , covid = sum(df.lintemp[CostaRica.indl ,]$covid), 
                     Nx = sum(df.lintemp[CostaRica.indl ,]$Nx))

Cuba.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Cuba.indl ,]))), 2, sum)) %>% 
           mutate(Country = "Cuba"    ,  covid = sum(df.lintemp[Cuba.indl ,]$covid), 
                     Nx = sum(df.lintemp[Cuba.indl ,]$Nx))


Ecuador.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Ecuador.indl  ,]))), 2, sum)) %>% 
           mutate( Country = "Ecuador"    , covid = sum(df.lintemp[Ecuador.indl  ,]$covid), 
                     Nx = sum(df.lintemp[Ecuador.indl  ,]$Nx))

Guatemala.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Guatemala.indl  ,]))), 2, sum)) %>% 
           mutate( Country = "Guatemala"  , covid = sum(df.lintemp[Guatemala.indl  ,]$covid), 
                     Nx = sum(df.lintemp[Guatemala.indl  ,]$Nx))

Jamaica.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Jamaica.indl   ,]))), 2, sum)) %>% 
           mutate( Country = "Jamaica"    ,  covid = sum(df.lintemp[Jamaica.indl   ,]$covid), 
                     Nx = sum(df.lintemp[Jamaica.indl   ,]$Nx))

Mexico.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Mexico.indl   ,]))), 2, sum)) %>% 
           mutate( Country = "Mexico"    ,covid = sum(df.lintemp[Mexico.indl   ,]$covid), 
                     Nx = sum(df.lintemp[Mexico.indl   ,]$Nx))

Panama.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Panama.indl   ,]))), 2, sum)) %>% 
           mutate( Country = "Panama"  ,  covid = sum(df.lintemp[Panama.indl   ,]$covid), 
                     Nx = sum(df.lintemp[Panama.indl   ,]$Nx))

Peru.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Peru.indl   ,]))), 2, sum)) %>% 
           mutate(Country = "Peru"    ,covid = sum(df.lintemp[Peru.indl   ,]$covid), 
                     Nx = sum(df.lintemp[Peru.indl   ,]$Nx))

Paraguay.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[Paraguay.indl   ,]))), 2, sum)) %>% 
           mutate(Country = "Paraguay"    ,  covid = sum(df.lintemp[Paraguay.indl   ,]$covid), 
                     Nx = sum(df.lintemp[Paraguay.indl   ,]$Nx))

UnitedStatesofAmerica.sum   <- getsum(apply(as.data.frame(t(as.matrix(fitdist1e[UnitedStatesofAmerica.indl   ,]))), 2, sum)) %>% 
           mutate(Country = "United States of America"  ,  covid = sum(df.lintemp[UnitedStatesofAmerica.indl   ,]$covid), 
                     Nx = sum(df.lintemp[UnitedStatesofAmerica.indl   ,]$Nx))

# ===========================================================================

sum.exc.df <- rbind(Argentina.sum  ,AntiguaandBarbudaa.sum   ,Bahamas.sum  ,Belize.sum  ,Barbados.sum ,Dominica.sum  ,DominicanRepublic.sum  ,Grenada.sum   ,Guyana.sum  ,Honduras.sum  ,Haiti.sum  ,SaintKittsandNevis.sum  ,SaintLucia.sum  ,Nicaragua.sum  ,ElSalvador.sum   ,Suriname.sum  ,TrinidadandTobago.sum   ,Uruguay.sum  ,SaintVincentandtheGrenadines.sum  ,VenezuelaBolivarianRepublicof.sum  ,BoliviaPlurinationalStateof.sum   ,Brazil.sum   ,Canada.sum  ,Chile.sum  ,Colombia.sum   ,CostaRica.sum  ,Cuba.sum  ,Ecuador.sum   ,Guatemala.sum   ,Jamaica.sum   ,Mexico.sum  ,Panama.sum  ,Peru.sum   ,Paraguay.sum  ,UnitedStatesofAmerica.sum  
) %>% 
  rename(low = lwr, high = uppr) %>%
  mutate(excess = round(mean - covid), 
         mean = round(mean), low = round(low), high = round(high),
         region = factor(Country, levels = c( "Argentina" ,"Antigua and Barbuda","Bahamas" ,"Belize","Bolivia (Plurinational State of)","Brazil","Barbados","Canada","Chile","Colombia","Costa Rica","Cuba","Dominica","Dominican Republic","Ecuador","Grenada","Guatemala","Guyana","Honduras","Haiti","Jamaica","Saint Kitts and Nevis","Saint Lucia","Mexico","Nicaragua","Panama","Peru","Paraguay","El Salvador","Suriname","Trinidad and Tobago","Uruguay","United States of America","Saint Vincent and the Grenadines","Venezuela (Bolivarian Republic of)")), 
         model = "Basic Linear") %>%
  dplyr::select(model, Country, Nx, covid, excess, low, mean, high)
```

## Mixed effects model in INLA, similar syntax to GLM

This model does approximate Bayesian inference for Latent Gaussian Models

``` r
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) # Run this first if you don't have INLA library CAUTION : HEAVY LOADING
library(INLA)

df.mixed  <-INLA::inla(excessr ~ covidr + sdi*cardiovasc_death_rate + quarter + 
                         diabetes_prevalence + aged_65_older + positive_rate - 1 + 
                         f(WHO_region, model = "iid") + f(group, model = "iid"), 
                  data = analysis_mod, 
                  control.predictor= list(compute=TRUE))

pred.in   <- cbind(df.mixed$summary.fitted.values$mean, df.mixed$summary.fitted.values$sd) 

set.seed(12345)
fit.dist3 <- t(apply(pred.in, 1, get.sample)) # Assume distribution of fit is Normal with mean and sd

ind.obs <- analysis_mod %>% filter(!is.na(excess)) %>% pull(row) # indices for no data


# For the observed replace predictions with actual excess rate estimates

for (i in 1:draws){
  fit.dist3[ind.obs,] <- analysis_mod$excessr[ind.obs]
}

fit.dist3e   <- analysis_mod$Nx     * fit.dist3 


# Results are by quarter, need to derive draws by country and then use those to derive draws by region etc
fit.dist3ef1 <- array(dim = c(194, draws))
for (j in 1:194){
  fit.dist3ef1[j,] <- apply(fit.dist3e[(j*4 - 3):(j*4),], 2, sum)
}


# Summarize with bootstrap confidence intervals based on the distribution of the draws
df.mixed.out   <- data.table::data.table(t(apply(fit.dist3ef1, 1, getsum2))) %>%
  rename(y = "50%", low = "2.5%", high = "97.5%") %>%
  cbind(df.lin) %>% filter(!is.na(y))  %>% 
  mutate(y    = ifelse(!is.na(excess), excess, y),
         low  = ifelse(!is.na(excess), NA, low),
         high = ifelse(!is.na(excess), NA, high), 
         source = ifelse(is.na(excess), "Predicted", "Observed")) %>% arrange(-y)


# Identify the rows for each region to aggregate the matrices accordingly
analysis_mod_paho<-analysis_mod %>% filter(WHO_region == "AMRO")
analysis_mod_paho<-analysis_mod_paho%>% mutate(row = 1:nrow(analysis_mod_paho))


Argentina.ind   <- analysis_mod_paho %>% filter(Country == "Argentina") %>% pull(row)
AntiguaandBarbudaa.ind   <- analysis_mod_paho %>% filter(Country == "Antigua and Barbuda") %>% pull(row)
Bahamas.ind   <- analysis_mod_paho %>% filter(Country == "Bahamas") %>% pull(row)
Belize.ind   <- analysis_mod_paho %>% filter(Country == "Belize") %>% pull(row)
Barbados.ind   <- analysis_mod_paho %>% filter(Country == "Barbados") %>% pull(row)
Dominica.ind   <- analysis_mod_paho %>% filter(Country == "Dominica") %>% pull(row)
DominicanRepublic.ind   <- analysis_mod_paho %>% filter(Country == "Dominican Republic") %>% pull(row)
Grenada.ind   <- analysis_mod_paho %>% filter(Country == "Grenada") %>% pull(row)
Guyana.ind   <- analysis_mod_paho %>% filter(Country == "Guyana") %>% pull(row)
Honduras.ind   <- analysis_mod_paho %>% filter(Country == "Honduras") %>% pull(row)
Haiti.ind   <- analysis_mod_paho %>% filter(Country == "Haiti") %>% pull(row)
SaintKittsandNevis.ind   <- analysis_mod_paho %>% filter(Country == "Saint Kitts and Nevis") %>% pull(row)
SaintLucia.ind   <- analysis_mod_paho %>% filter(Country == "Saint Lucia") %>% pull(row)
Nicaragua.ind   <- analysis_mod_paho %>% filter(Country == "Nicaragua") %>% pull(row)
ElSalvador.ind   <- analysis_mod_paho %>% filter(Country == "El Salvador ") %>% pull(row)
Suriname.ind   <- analysis_mod_paho %>% filter(Country == "Suriname") %>% pull(row)
TrinidadandTobago.ind   <- analysis_mod_paho %>% filter(Country == "Trinidad and Tobago") %>% pull(row)
Uruguay.ind   <- analysis_mod_paho %>% filter(Country == "Uruguay") %>% pull(row)
SaintVincentandtheGrenadines.ind   <- analysis_mod_paho %>% filter(Country == "Saint Vincent and the Grenadines") %>% pull(row)
VenezuelaBolivarianRepublicof.ind   <- analysis_mod_paho %>% filter(Country == "Venezuela (Bolivarian Republic of)") %>% pull(row)
BoliviaPlurinationalStateof.ind   <- analysis_mod_paho %>% filter(Country == "Bolivia (Plurinational State of)") %>% pull(row)
Brazil.ind   <- analysis_mod_paho %>% filter(Country == "Brazil") %>% pull(row)
Canada.ind   <- analysis_mod_paho %>% filter(Country == "Canada") %>% pull(row)
Chile.ind   <- analysis_mod_paho %>% filter(Country == "Chile") %>% pull(row)
Colombia.ind   <- analysis_mod_paho %>% filter(Country == "Colombia") %>% pull(row)
CostaRica.ind   <- analysis_mod_paho %>% filter(Country == "Costa Rica") %>% pull(row)
Cuba.ind   <- analysis_mod_paho %>% filter(Country == "Cuba") %>% pull(row)
Ecuador.ind   <- analysis_mod_paho %>% filter(Country == "Ecuador") %>% pull(row)
Guatemala.ind   <- analysis_mod_paho %>% filter(Country == "Guatemala") %>% pull(row)
Jamaica.ind   <- analysis_mod_paho %>% filter(Country == "Jamaica") %>% pull(row)
Mexico.ind   <- analysis_mod_paho %>% filter(Country == "Mexico") %>% pull(row)
Panama.ind   <- analysis_mod_paho %>% filter(Country == "Panama") %>% pull(row)
Peru.ind   <- analysis_mod_paho %>% filter(Country == "Peru") %>% pull(row)
Paraguay.ind   <- analysis_mod_paho %>% filter(Country == "Paraguay") %>% pull(row)
UnitedStatesofAmerica.ind   <- analysis_mod_paho %>% filter(Country == "United States of America") %>% pull(row)


# Summarise within each aggregate, first by column (draw) and then across draws by aggregate


Argentina.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Argentina.ind,]))), 2, sum)) %>%
                    mutate(Country = "Argentina",
                           covid = sum(analysis_mod_paho[Argentina.ind,]$cov_deaths), 
                            Nx = .25*sum(analysis_mod_paho[Argentina.ind,]$Nx)) 


AntiguaandBarbudaa.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[AntiguaandBarbudaa.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Antigua and Barbuda",   covid = sum(analysis_mod_paho[AntiguaandBarbudaa.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[AntiguaandBarbudaa.ind ,]$Nx))

Bahamas.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Bahamas.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Bahamas",  covid = sum(analysis_mod_paho[Bahamas.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Bahamas.ind ,]$Nx))

Belize.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Belize.ind ,]))), 2, sum)) %>% 
           mutate(Country = "Belize", covid = sum(analysis_mod_paho[Belize.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Belize.ind ,]$Nx))

Barbados.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Barbados.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Barbados" ,covid = sum(analysis_mod_paho[Barbados.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Barbados.ind ,]$Nx))

Dominica.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Dominica.ind ,]))), 2, sum)) %>% 
           mutate(Country = "Dominica", 
                  covid = sum(analysis_mod_paho[Dominica.ind ,]$cov_deaths), 
                  Nx = sum(analysis_mod_paho[Dominica.ind ,]$Nx))

DominicanRepublic.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[DominicanRepublic.ind  ,]))), 2, sum)) %>% 
           mutate( Country = "Dominican Republic"  ,  covid = sum(analysis_mod_paho[DominicanRepublic.ind  ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[DominicanRepublic.ind  ,]$Nx))

Grenada.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Grenada.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Grenada" , covid = sum(analysis_mod_paho[Grenada.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Grenada.ind ,]$Nx))

Guyana.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Guyana.ind  ,]))), 2, sum)) %>% 
           mutate( Country = "Guyana"   , covid = sum(analysis_mod_paho[Guyana.ind  ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Guyana.ind  ,]$Nx))

Honduras.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Honduras.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Honduras"  , covid = sum(analysis_mod_paho[Honduras.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Honduras.ind ,]$Nx))

Haiti.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Haiti.ind,]))), 2, sum)) %>% 
           mutate( Country = "Haiti"   ,  covid = sum(analysis_mod_paho[Haiti.ind,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Haiti.ind,]$Nx))

SaintKittsandNevis.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[SaintKittsandNevis.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Saint Kitts and Nevis"  ,   covid = sum(analysis_mod_paho[SaintKittsandNevis.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[SaintKittsandNevis.ind ,]$Nx))

SaintLucia.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[SaintLucia.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Saint Lucia"  , covid = sum(analysis_mod_paho[SaintLucia.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[SaintLucia.ind ,]$Nx))

Nicaragua.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Nicaragua.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Nicaragua"  ,  covid = sum(analysis_mod_paho[Nicaragua.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Nicaragua.ind ,]$Nx))

ElSalvador.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[ElSalvador.ind ,]))), 2, sum)) %>% 
           mutate( Country = "El Salvador"   , covid = sum(analysis_mod_paho[ElSalvador.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[ElSalvador.ind ,]$Nx))

Suriname.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Suriname.ind,]))), 2, sum)) %>% 
           mutate( Country = "Suriname" ,    covid = sum(analysis_mod_paho[Suriname.ind,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Suriname.ind,]$Nx))

TrinidadandTobago.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[TrinidadandTobago.ind  ,]))), 2, sum)) %>% 
           mutate( Country = "Trinidad and Tobago"  , covid = sum(analysis_mod_paho[TrinidadandTobago.ind  ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[TrinidadandTobago.ind  ,]$Nx))

Uruguay.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Uruguay.ind,]))), 2, sum)) %>% 
           mutate( Country = "Uruguay"   ,  covid = sum(analysis_mod_paho[Uruguay.ind,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Uruguay.ind,]$Nx))

SaintVincentandtheGrenadines.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[SaintVincentandtheGrenadines.ind  ,]))), 2, sum)) %>% 
           mutate( Country = "Saint Vincent and the Grenadines"  ,covid = sum(analysis_mod_paho[SaintVincentandtheGrenadines.ind  ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[SaintVincentandtheGrenadines.ind  ,]$Nx))


VenezuelaBolivarianRepublicof.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[VenezuelaBolivarianRepublicof.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Venezuela (Bolivarian Republic of)" ,covid = sum(analysis_mod_paho[VenezuelaBolivarianRepublicof.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[VenezuelaBolivarianRepublicof.ind ,]$Nx))


# non missing 
BoliviaPlurinationalStateof.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[BoliviaPlurinationalStateof.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Bolivia (Plurinational State of)",  covid = sum(analysis_mod_paho[BoliviaPlurinationalStateof.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[BoliviaPlurinationalStateof.ind ,]$Nx))

Brazil.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Brazil.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Brazil" ,covid = sum(analysis_mod_paho[Brazil.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Brazil.ind ,]$Nx))

Canada.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Canada.ind ,]))), 2, sum)) %>% 
           mutate(Country = "Canada" , covid = sum(analysis_mod_paho[Canada.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Canada.ind ,]$Nx))

Chile.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Chile.ind ,]))), 2, sum)) %>% 
           mutate(Country = "Chile" , covid = sum(analysis_mod_paho[Chile.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Chile.ind ,]$Nx))

Colombia.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Colombia.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Colombia"  , covid = sum(analysis_mod_paho[Colombia.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Colombia.ind ,]$Nx))

CostaRica.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[CostaRica.ind ,]))), 2, sum)) %>% 
           mutate( Country = "Costa Rica"  , covid = sum(analysis_mod_paho[CostaRica.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[CostaRica.ind ,]$Nx))

Cuba.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Cuba.ind ,]))), 2, sum)) %>% 
           mutate(Country = "Cuba"    ,  covid = sum(analysis_mod_paho[Cuba.ind ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Cuba.ind ,]$Nx))


Ecuador.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Ecuador.ind  ,]))), 2, sum)) %>% 
           mutate( Country = "Ecuador"    , covid = sum(analysis_mod_paho[Ecuador.ind  ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Ecuador.ind  ,]$Nx))

Guatemala.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Guatemala.ind  ,]))), 2, sum)) %>% 
           mutate( Country = "Guatemala"  , covid = sum(analysis_mod_paho[Guatemala.ind  ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Guatemala.ind  ,]$Nx))

Jamaica.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Jamaica.ind   ,]))), 2, sum)) %>% 
           mutate( Country = "Jamaica"    ,  covid = sum(analysis_mod_paho[Jamaica.ind   ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Jamaica.ind   ,]$Nx))

Mexico.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Mexico.ind   ,]))), 2, sum)) %>% 
           mutate( Country = "Mexico"    ,covid = sum(analysis_mod_paho[Mexico.ind   ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Mexico.ind   ,]$Nx))

Panama.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Panama.ind   ,]))), 2, sum)) %>% 
           mutate( Country = "Panama"  ,  covid = sum(analysis_mod_paho[Panama.ind   ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Panama.ind   ,]$Nx))

Peru.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Peru.ind   ,]))), 2, sum)) %>% 
           mutate(Country = "Peru"    ,covid = sum(analysis_mod_paho[Peru.ind   ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Peru.ind   ,]$Nx))

Paraguay.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[Paraguay.ind   ,]))), 2, sum)) %>% 
           mutate(Country = "Paraguay"    ,  covid = sum(analysis_mod_paho[Paraguay.ind   ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[Paraguay.ind   ,]$Nx))

UnitedStatesofAmerica.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[UnitedStatesofAmerica.ind   ,]))), 2, sum)) %>% 
           mutate(Country = "United States of America"  ,  covid = sum(analysis_mod_paho[UnitedStatesofAmerica.ind   ,]$cov_deaths), 
                     Nx = sum(analysis_mod_paho[UnitedStatesofAmerica.ind   ,]$Nx))







sum.exc.exc <- rbind(Argentina.sum.exc  ,AntiguaandBarbudaa.sum.exc   ,Bahamas.sum.exc  ,Belize.sum.exc  ,Barbados.sum.exc ,Dominica.sum.exc  ,DominicanRepublic.sum.exc  ,Grenada.sum.exc   ,Guyana.sum.exc  ,Honduras.sum.exc  ,Haiti.sum.exc  ,SaintKittsandNevis.sum.exc  ,SaintLucia.sum.exc  ,Nicaragua.sum.exc  ,ElSalvador.sum.exc   ,Suriname.sum.exc  ,TrinidadandTobago.sum.exc   ,Uruguay.sum.exc  ,SaintVincentandtheGrenadines.sum.exc  ,VenezuelaBolivarianRepublicof.sum.exc  ,BoliviaPlurinationalStateof.sum.exc   ,Brazil.sum.exc   ,Canada.sum.exc  ,Chile.sum.exc  ,Colombia.sum.exc   ,CostaRica.sum.exc  ,Cuba.sum.exc  ,Ecuador.sum.exc   ,Guatemala.sum.exc   ,Jamaica.sum.exc   ,Mexico.sum.exc  ,Panama.sum.exc  ,Peru.sum.exc   ,Paraguay.sum.exc  ,UnitedStatesofAmerica.sum.exc  
) %>% 
  rename(low = lwr, high = uppr) %>%
  mutate(excess = round(mean - covid), 
         mean = round(mean), low = round(low), high = round(high),
         region = factor(Country, levels = c( "Argentina" ,"Antigua and Barbuda","Bahamas" ,"Belize","Bolivia (Plurinational State of)","Brazil","Barbados","Canada","Chile","Colombia","Costa Rica","Cuba","Dominica","Dominican Republic","Ecuador","Grenada","Guatemala","Guyana","Honduras","Haiti","Jamaica","Saint Kitts and Nevis","Saint Lucia","Mexico","Nicaragua","Panama","Peru","Paraguay","El Salvador","Suriname","Trinidad and Tobago","Uruguay","United States of America","Saint Vincent and the Grenadines","Venezuela (Bolivarian Republic of)")), 
         model = "Mixed Effects") %>%
  dplyr::select(model, region, Nx, covid, excess, low, mean, high)
```

``` r
# make output into same format as other datasets for final dataset 

#excess     = round(observed - expected),                      # excess deaths
#excess     = round(mean - covid)

df_excess<-sum.exc.df 
df_excess<-sum.exc.df%>%
            rename(country=Country)%>%
            mutate(value= (excess/Nx)*-1, # (observed - covid)/Nx ratio)
                              year=2020,
                              indicator="ExcessDeathRatio")%>%
            relocate(value, .after = indicator)%>%
            arrange(year, country, indicator)%>%
            select("year" , "country", "indicator", "value")

head(df_excess,5)
```

    ## # A tibble: 5 x 4
    ##    year country             indicator          value
    ##   <dbl> <chr>               <chr>              <dbl>
    ## 1  2020 Antigua and Barbuda ExcessDeathRatio  -703. 
    ## 2  2020 Argentina           ExcessDeathRatio    78.1
    ## 3  2020 Bahamas             ExcessDeathRatio -1435. 
    ## 4  2020 Barbados            ExcessDeathRatio -3486. 
    ## 5  2020 Belize              ExcessDeathRatio    38.2
