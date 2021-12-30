# Emission driven Hector was underestimating CH4 & N2O RF related to natural CH4 & N2O emissions. In this script 
# update RF optmize parameter values. 

# Set up the working environment
library(assertthat)
library(data.table)
library(dplyr)
library(ggplot2)
library(hector) # Requires devlopmental Hector v3  

# Lets see how much that this can change as beta is adjusted. 
DIR <- here::here()
INPUT_DIR <- file.path(DIR, "input")

# Comparison Data -------------------------------------------------------------------------------------------
ini_files <- list.files(INPUT_DIR, "ini", full.names = TRUE)
lapply(ini_files, function(f){
  nm <- gsub(pattern = "hector_|.ini", replacement = "", basename(f))
  newcore(inifile = f, name = nm)
}) -> 
  core_list


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

# Format the data & match up the AR6 output variables with the Hector variable names. 
files %>% 
  lapply(get_cmip_data) %>% 
  bind_rows() %>% 
  select(year, variable, value, scenario = scn, percent) %>% 
  as.data.table() %>% 
  filter(is.na(percent))  -> 
  ar6_data_og

ar6_names <- c("ch4","n2o")
hector_names <- c(RF_CH4(), RF_N2O()) 

ar6_data_og %>% 
  filter(variable %in% ar6_names) %>% 
  left_join(data.frame(variable = ar6_names, 
                       hvar = hector_names)) %>% 
  select(year, variable = hvar, ar6_val = value, scenario) -> 
  ar6_results 

# 1. N2O & CH4 natural emissions optimzation ----------------------------------------------------------------
# The N2O & CH4 Natural emission optmization, becacuse the N2O & CH4 RF equations depend on both N2O, CH4, & CO2 
# concentrations the CO2 concenstrations should be constrained and the natural emissions need to be optmized 
# together. 
# Set up the Hector cores  -------------------------------------------------------------------------------------------
ini_files <- list.files(INPUT_DIR, "ini", full.names = TRUE)
lapply(ini_files, function(f){
  nm <- gsub(pattern = "hector_|.ini", replacement = "", basename(f))
  newcore(inifile = f, name = nm)
}) -> 
  core_list

# Optimize Natural N2O and CH4  -------------------------------------------------------------------------------------------
optim_nat_hector <- function(comp_data, par){
  
  assert_that(has_name(x = par, which = c("nat_n2o")))
  
  # Set up the Hector core with the new natural emissions values
  lapply(core_list, setvar, dates = NA, var =  NATURAL_CH4(), value = par[["nat_ch4"]], unit =  "Tg CH4")
  dates <- 1750:2100
  lapply(core_list, setvar, dates = dates, var =  NAT_EMISSIONS_N2O(), value = rep(par[["nat_n2o"]], length(dates)), unit =  "Tg N")
  lapply(core_list, reset)
  
  # Run the Hector cores 
  multi_mean_mse <- tryCatch({
    lapply(core_list, hector::run, runtodate = 2100)
    hector_data <- bind_rows(lapply(core_list, hector::fetchvars, dates = 1740:2100, vars = c(RF_CH4(), RF_N2O())))
    
    # Get 
    hector_data %>% 
      inner_join(comp_data) %>% 
      mutate(SE = (value - ar6_val)^2) %>% 
      group_by(scenario, variable) %>% 
      summarise(MSE = mean(SE, na.rm = TRUE)) %>% 
      ungroup() %>% 
      pull(MSE) %>% 
      mean(na.rm = TRUE)}, error = function(error){99999})
  
  
  return(multi_mean_mse)
  
}

par <- c(5, 300)
names(par) <- c("nat_n2o", "nat_ch4")

result <- optim(par = par, fn = optim_nat_hector, comp_data = ar6_results)

# Plot Results  -------------------------------------------------------------------------------------------

# Run Hector with the new natural emission values
lapply(core_list, setvar, dates = NA, var =  NATURAL_CH4(), value = result$par[['nat_ch4']] , unit =  "Tg CH4")
dates <- 1750:2100
lapply(core_list, setvar, dates = dates, var =  NAT_EMISSIONS_N2O(), value = rep(result$par[['nat_n2o']] , length(dates)), unit =  "Tg N")
lapply(core_list, reset)
lapply(core_list, run)
hector <- bind_rows(lapply(core_list, fetchvars, dates = 1740:2100, vars = c(RF_TOTAL(), RF_N2O(), RF_CH4())))
hector$name <- "nat. ch4 & n2o"

ar6_results %>%  
  rename(value = ar6_val) %>% 
  mutate(name = "ar6") -> 
  ar6_results_new

bind_rows(hector_data_new) %>% 
  filter(scenario %in% c(hector_data$scenario)) %>% 
  filter(variable %in% c(RF_CH4(), RF_N2O())) %>% 
  filter(year <= 2100) %>%  
  ggplot(aes(year, value, color = name, groupby = scenario, linetype = name)) + 
  geom_line() + 
  facet_wrap("variable")


# Compare default & new vals  -------------------------------------------------------------------------------------------
# Run Hector with the new natural emission values
lapply(core_list, hector::setvar, dates = NA, var =  NATURAL_CH4(), value = result$par[['nat_ch4']] , unit =  "Tg CH4")
dates <- 1750:2100
lapply(core_list, hector::setvar, dates = dates, var =  NAT_EMISSIONS_N2O(), value = rep(result$par[['nat_n2o']] , length(dates)), unit =  "Tg N")
lapply(core_list, hector::reset)
lapply(core_list, hector::run, runtodate=2100)
hector_data_new <- bind_rows(lapply(core_list, hector::fetchvars, dates = 1740:2100, vars = c(RF_TOTAL(), RF_N2O(), RF_CH4())))
hector_data_new$name <- "new natural ch4 & n2o"

ar6_results$name <- "ar6"
names(ar6_results) <- c("year", "variable", "value", "scenario", "name")

bind_rows(hector_data_new, ar6_results) %>% 
  filter(scenario %in% c(hector_data$scenario)) %>% 
  filter(year <= 2100) %>%  
  filter(variable %in% c(RF_CH4(), RF_N2O())) %>% 
  ggplot(aes(year, value, color = name, groupby = scenario, linetype = name)) + 
  geom_line() + 
  facet_grid(variable~scenario, scales = "free")


# # Results from the fit 
# nat_n2o    nat_ch4 
# 9.713487 340.950746 





