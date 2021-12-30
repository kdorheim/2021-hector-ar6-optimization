# Emission driven Hector was underestimating CO2 RF while CO2 constrained Hector was spot on when it came to 
# to CO2 RF. I think that this related to carbon cycle feedbacks, paritcuarly, beta. 

# Set up the working environment
library(assertthat)
library(data.table)
library(dplyr)
library(ggplot2)
library(hector) # Requires devlopmental Hector v3  

# Lets see how much that this can change as beta is adjusted. 
DIR <- here::here()
INPUT_DIR <- file.path(DIR, "input")

# Support functions & hector cores -------------------------------------------------------
optim_beta_hector <- function(comp_data, par){
  
  
  assert_that(has_name(x = par, which = BETA()))
  
  # Set up the Hector core with the new natural emissions values
  lapply(core_list, setvar, dates = NA, var = BETA(), value = par[[BETA()]], unit = getunits(BETA()))
  
  # Calculate the MSE between the comparison data & the Hector outputs. 
  multi_mean_mse <- tryCatch({
    lapply(core_list, hector::reset)
    lapply(core_list, hector::run, runtodate = 2100)
    hector_data <- bind_rows(lapply(core_list, fetchvars, dates = 1740:2100, vars = unique(comp_data$variable)))
    
    # Get 
    hector_data %>% 
      inner_join(comp_data) %>% 
      mutate(SE = (value - ar6_val)^2) %>% 
      group_by(scenario, variable) %>% 
      summarise(MSE = mean(SE)) %>% 
      ungroup() %>% 
      pull(MSE) %>% 
      mean()}, error = function(error){99999})
  
  return(multi_mean_mse)
  
}

ini_files <- list.files(INPUT_DIR, "ini", full.names = TRUE)
lapply(ini_files, function(f){
  nm <- gsub(pattern = "hector_|.ini", replacement = "", basename(f))
  newcore(inifile = f, name = nm)
}) -> 
  core_list

# Fit to CO2 concentrations ---------------------------------------------------------------
# What is the best fit of beta if only CO2 concentrations are used as comparison data? 
# Read in the CO2 concentraion data from the input tables. 
file.path(INPUT_DIR, "tables") %>% 
  list.files(full.names = TRUE, pattern = "csv") %>% 
  lapply(function(x){
    
    scn <- substr(x = basename(x), start = 1, stop = 6)
    data <- as.data.table(read.csv(x, stringsAsFactors = FALSE, comment.char = ";"))
    melt(data, id.vars = 'Date') %>% 
      select(year = Date, variable, value) %>%  
      mutate(scenario = scn) %>% 
      filter(variable == CO2_CONSTRAIN()) %>% 
      mutate(variable = ATMOSPHERIC_CO2())
    
  }) %>% 
  bind_rows() %>%  
  rename(ar6_val = value) -> 
  ar6_co2


# Initalize a starting parameter. 
par <- c(0.35)
names(par) <- BETA()

# Optmize Hector's beta parameters to get the best fit of CO2 concentrations 
result_co2conc <- optim(par = par, fn = optim_beta_hector, comp_data = ar6_co2)


# Whaat about when fitting to CO2 RF -----------------------------------------------------------------------
# Define some helper functions are used to format and extract the AR6 CMIP data. 
format_data <- function(path){
  
  x <- gsub(pattern = ".csv", replacement = "", x = basename(path))
  xx <- unlist(strsplit(x, split = "_"))
  
  if (length(xx) == 3 ){
    out <- as.data.frame(matrix(xx, nrow = 1, dimnames = list(NULL, c("output", "scn", "time"))))
  }
  
  if (length(xx) == 4 ){
    out <- as.data.frame(matrix(xx, nrow = 1, dimnames = list(NULL, c("output", "scn", "time", "percent"))))
  }
  
  return(out)
  
}

get_cmip_data <- function(path){
  
  dat <- read.csv(path, stringsAsFactors = FALSE)
  info <- format_data(path)
  
  dat %>% 
    melt(id.vars = c( "year"), 
         variable.name = "variable", value.name = "value") %>% 
    cbind(info) -> 
    out
  
  return(out)
}

# Set up the ar6 comparison data, note that this comes from the repo https://github.com/IPCC-WG1/Chapter-7.git 
files   <- list.files(file.path(DIR, "AR6","SSPs"), pattern = "ERF", full.names = TRUE)  

files %>% 
  lapply(get_cmip_data) %>% 
  bind_rows() %>% 
  select(year, variable, value, scenario = scn, percent) %>% 
  as.data.table() %>% 
  filter(is.na(percent)) %>% 
  select(year, variable, value, scenario) %>% 
  filter(variable %in% c("co2", "ch4", "n2o", "total")) %>%  
  mutate(variable = paste0("F", toupper(variable))) -> 
  all_ar6_rf

all_ar6_rf %>% 
  filter(variable ==  RF_CO2()) %>%  
  rename(ar6_val = value) -> 
  ar6_rf


# Initalize a starting parameter. 
par <- c(0.35)
names(par) <- BETA()

# Optmize Hector's beta parameters to get the best fit of CO2 concentrations 
result_co2rf <- optim(par = par, fn = optim_beta_hector, comp_data = ar6_rf)

# Awesome both fits of beta = 0.1332! -------------------------------------------------------------------
# Run Hector with the new beta value and compare CO2 RF & the total RF 
lapply(core_list, hector::reset)
lapply(core_list, hector::setvar, dates = NA, var = hector::BETA(), value = 0.13, unit = getunits(hector::BETA()))
lapply(core_list, hector::reset)
lapply(core_list, hector::run, runtodate = 2100)
hector_data <- bind_rows(lapply(core_list, hector::fetchvars, dates = 1850:2100, vars = c(RF_CO2(), RF_TOTAL(), NPP())))
hector_data$source <- "Hector"


lapply(core_list, hector::reset)
lapply(core_list, hector::setvar, dates = NA, var = hector::BETA(), value = 0.36, unit = getunits(hector::BETA()))
lapply(core_list, hector::reset)
lapply(core_list, hector::run, runtodate = 2100)
hector_data1 <- bind_rows(lapply(core_list, hector::fetchvars, dates = 1850:2100, vars = c(RF_CO2(), NPP(), RF_TOTAL())))
hector_data1$source <- "default"

all_ar6_rf %>%
  filter(variable %in% RF_CO2()) %>% 
  bind_rows(all_ar6_rf %>% 
              filter(variable == "FTOTAL") %>% 
              mutate(variable = RF_TOTAL())) -> 
  comp_rf
comp_rf$source <- "AR6"

bind_rows(hector_data, 
          comp_rf) %>% 
  filter(year %in% 1900:2100) %>% 
  filter(variable == RF_CO2()) %>% 
  filter(scenario %in% hector_data$scenario) %>%  
  ggplot(aes(year, value, color = source, groupby = scenario)) + 
  geom_line() + 
  NULL


hector_data %>%  
  bind_rows(hector_data1) %>%  
  filter(variable %in% c(RF_CO2(), NPP())) %>% 
  ggplot(aes(year, value, color = source, groupby = scenario)) + 
  geom_line() + 
  facet_wrap("variable", scales = "free")
