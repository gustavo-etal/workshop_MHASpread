#=========================================================#
#       Control actions disease workshop code #4          #
#=========================================================#
#Note:
library(MHASpread)                         # Load the package with all SEIR model functions


#==========================================================================#
# initial scenarios without control ---------------------------------------
#==========================================================================#

##  Scenario 4 ----
# Select the farm to be infected

population <- MHASpread::population                                                # Get the population data example
population$I_bov_pop[population$node== 196734] <- 40                                  # Infected 10 bovine in fasrm in livramento
# select the in and out farm dynamics (movements, births, deaths etc.. )
events <- MHASpread::events                                                        # load the events database


## run the model in shiny
# load_App()
## run model in R
model_4     <- SEIR_model(population = population,                                                 #  Population database
                          events = events,                                                         #  Events database
                          simulation_name = "scenario_2_test",                                     #  Simulation tag name
                          days_of_simulation = 20,                                                 #  Population database
                          initial_day_simulation=1,                                                #  Initial day of simulation
                          max_distance_in_km= 40,                                                  #  Max distance kernel by local disease spread
                          num_threads=1,                                                           #  Number of CPU to parallel tasks; set 1 to not overload your computer
                          a = 0.012,                                                               #  To set kernel curve max infection rate (S*I)/N when animals are in the same area
                          b =  0.6 ,                                                               #  Shape of the kernel curve
                          beta_bov_to_bov= c(min = 0.01833333, mode = 0.025, max = 0.05666667),    #  Transmission coefficient of bovine infects bovine
                          beta_bov_to_swi= c(min = 0.01833333, mode = 0.025, max = 0.05666667),    #  Transmission coefficient of bovine infects swine
                          beta_bov_to_SR=c(min = 0.012, mode = 0.031, max = 0.065),                #  Transmission coefficient of bovine infects small ruminants
                          lambda1_bov=c(min = 3, mode = 5.9, max = 16),                            #  Rate from exposed (E) to infectious (I) in bovine
                          lambda2_bov=c(min = 6, mode = 15, max = 20),                             #  Rate from infectious (I) to recovered (R) in bovine
                          beta_swi_to_swi=c(min = 0.044, mode = 0.14, max = 0.33) ,                #  Transmission coefficient of swine infects swine
                          beta_swi_to_bov=c(min = 0.014, mode = 0.044, max = 0.033) ,              #  Transmission coefficient of swine infects bovine
                          beta_swi_to_SR= c(min = 0.014, mode = 0.044, max = 0.033),               #  Transmission coefficient of swine infects small ruminants
                          lambda1_swi=c(min = 3, mode = 5.9, max = 16),                            #  Rate from exposed (E) to infectious (I) in swine
                          lambda2_swi=c(min = 5, mode = 6.44, max = 14) ,                          #  Rate from infectious (I) to recovered (R) in swine
                          beta_SR_to_SR=c(min = 0.16, mode = 0.24, max = 0.5) ,                    #  Transmission coefficient of small ruminants infects small ruminants
                          beta_SR_to_bov=c(min=0.012,mode=0.031,max= 0.033) ,                      #  Transmission coefficient of small ruminants infects bovine
                          beta_SR_to_swi=c(min = 0.006, mode = 0.024, max = 0.09),                 #  Transmission coefficient of small ruminants infects swine
                          lambda1_SR=c(min = 4, mode = 5, max = 14) ,                              #  Rate from exposed (E) to infectious (I) in small ruminants
                          lambda2_SR=c(min = 6, mode = 15, max = 20) )                             #  Rate from infectious (I) to recovered (R) in bovine


#==========================================================================#
# Control action modelling          ---------------------------------------
#==========================================================================#
# here use the previous model output to set the initial condition for the control actions
population <- update_population_in_model(original_population = MHASpread::population,
                                         model_4$populationdb )                                        # Complete  and update the  population database
initial_day_control <- 21                                                                              # Seed the initial day of infection in the day 21
days_of_control_action <- 90                                                                           # Select 30 days working on control actions
banco_summary <- c()                                                                                   # Create an empty object to save the model output
days_to_test <- seq(initial_day_control, (initial_day_control+days_of_control_action), by= 7)          # Create the daya that will be used to test and update the buffer control areas



# Here open a loop to simulate the control action adn FDM infection by day
t1 <- Sys.time()
for (day in initial_day_control:(initial_day_control+days_of_control_action)) {
  print(day)


  #====================================#
  #    check point of the loop      ----
  #====================================#
  ## check if there are no more infectious or exposed animals, if so break the loop
  mystatement <- id_of_infectious_farms(population = population, only_infected_comp = F)
  if (identical(mystatement, numeric(0))) {
    print("No more infectious animals")
    break
  }





  #=======================================#
  ##   establish the control areas zones ----
  #=======================================#
  # note: Here, we select the days to test the animal population and then update control zones

  if (day %in% days_to_test) {
    control_zones_areas <- assign_control_zones(population = population, infected_size = 3,buffer_size = 7, surveillance_size = 15)
  }

  #=======================================#
  ##      farm standstill          ----
  #=======================================#
  events <- movement_control_ca(population = population,                #  Population database
                                events = events,                        #  Events database
                                day = day,                              #  Day of start the control movements
                                control_zones_db = control_zones_areas, #  Object with the farms in control areas
                                ban_length = 30,                        #  30 days of movements ban
                                infected_zone = T,                      #  Animal ban will be applied to infected zone
                                buffer_zone = T,                        #  Animal ban will be applied to buffer zone
                                surveillance_zone = T,                  #  Animal ban will be applied to surveillance zone
                                direct_contacts = T,                    #  Ban farm outside of control zones with contact with positive farms
                                traceback_length = 1)                   #  Traceback in-going animals movements of infected farms
  print("movement ban done!")

  #=======================================#
  ##      depopulation per farm       ----
  #=======================================#

  population <- depop_farms_ca(population = population,                 #  Population database
                               limit_per_day_farms = 4,                #  10 farm will be depopulated by day
                               infected_zone = T,                       #  Depopulation will be applied to infected zone
                               buffer_zone = F,                         #  Depopulation will be applied to buffer zone
                               surveillance_zone = F,                   #  Surveillance will be applied to infected zone
                               farms_in_zones = control_zones_areas)    #  Object with the farms in control areas
  print("Depopulation done!")

  #=======================================#
  ##      vaccination per farm        ----
  #=======================================#
  #select the animals vaccinated in the previous round
  farms_prev_vaccinated <- population%>%
    mutate(animals_vaccinated = rowSums(select(., starts_with("V_")))) %>%
    filter(animals_vaccinated>0) %>%
    pull(node)


  # in case of no animals vaccinated in previous rounds set NA
  if (identical(farms_prev_vaccinated, numeric(0))) {
    farms_prev_vaccinated <- NA
  }


  # Run the vaccination function bases on the selected paramenters
  output_vacc_pop <- vaccinate_farms_ca(population = population,
                                        limit_per_day_farms = 100,
                                        infected_zone = T,
                                        vacc_infectious_farms = T,
                                        buffer_zone = T,
                                        surveillance_zone = F,
                                        control_zones_db = control_zones_areas,
                                        vacc_bovine = T,
                                        vacc_swine = F,
                                        vacc_small = T,
                                        vaccine_efficacy =0.9,
                                        vaccinated_farms = farms_prev_vaccinated)
  print("vaccination done!")
  #update databases

  farms_prev_vaccinated <-  output_vacc_pop$vaccinated_farms
  population <- output_vacc_pop$population


  #=======================================#
  ##             SEIR dynamics        ----
  #=======================================#
  # if there are no more exposed or infected farms skip the SEIR model dynamics
  if (identical(id_of_infectious_farms(population, only_infected_comp = F),numeric(0))) {
    next
  }



  model_in_loop     <- SEIR_model(population = population,                                                 #  Population database
                                  events = events,                                                         #  Events database
                                  simulation_name = "scenario_4_control",                                          #  Simulation tag name
                                  days_of_simulation = 1,                                                  #  Population database
                                  initial_day_simulation=day,                                              #  Initial day of simulation
                                  max_distance_in_km= 40,                                                  #  Max distance kernel by local disease spread
                                  num_threads=1,                                                           #  Number of CPU to parallel tasks; set 1 to not overload your computer
                                  a = 0.012,                                                               #  To set kernel curve max infection rate (S*I)/N when animals are in the same area
                                  b =  0.6 ,                                                               #  Shape of the kernel curve
                                  beta_bov_to_bov= c(min = 0.01833333, mode = 0.025, max = 0.05666667),    #  Transmission coefficient of bovine infects bovine
                                  beta_bov_to_swi= c(min = 0.01833333, mode = 0.025, max = 0.05666667),    #  Transmission coefficient of bovine infects swine
                                  beta_bov_to_SR=c(min = 0.012, mode = 0.031, max = 0.065),                #  Transmission coefficient of bovine infects small ruminants
                                  lambda1_bov=c(min = 6, mode = 6, max = 16),                              #  Rate from exposed (E) to infectious (I) in bovine
                                  lambda2_bov=c(min = 6, mode = 15, max = 20),                             #  Rate from infectious (I) to recovered (R) in bovine
                                  beta_swi_to_swi=c(min = 0.044, mode = 0.14, max = 0.33) ,                #  Transmission coefficient of swine infects swine
                                  beta_swi_to_bov=c(min = 0.014, mode = 0.044, max = 0.033) ,              #  Transmission coefficient of swine infects bovine
                                  beta_swi_to_SR= c(min = 0.014, mode = 0.044, max = 0.033),               #  Transmission coefficient of swine infects small ruminants
                                  lambda1_swi=c(min = 3, mode = 5.9, max = 16),                            #  Rate from exposed (E) to infectious (I) in swine
                                  lambda2_swi=c(min = 5, mode = 6.44, max = 14) ,                          #  Rate from infectious (I) to recovered (R) in swine
                                  beta_SR_to_SR=c(min = 0.16, mode = 0.24, max = 0.5) ,                    #  Transmission coefficient of small ruminants infects small ruminants
                                  beta_SR_to_bov=c(min=0.012,mode=0.031,max= 0.033) ,                      #  Transmission coefficient of small ruminants infects bovine
                                  beta_SR_to_swi=c(min = 0.006, mode = 0.024, max = 0.09),                 #  Transmission coefficient of small ruminants infects swine
                                  lambda1_SR=c(min = 4, mode = 5, max = 14) ,                              #  Rate from exposed (E) to infectious (I) in small ruminants
                                  lambda2_SR=c(min = 6, mode = 15, max = 20) )                             #  Rate from infectious (I) to recovered (R) in bovine

  # update the population and model summary

  # Update population database ----
  population <- update_population_in_model(original_population = population, model_in_loop$populationdb)


  banco_summary <- rbind(banco_summary, model_in_loop)


}

t2 <- Sys.time()

t2-t1

#1.433649 hours
#===========================================================================#
# plot output model --------------------------------------------------------
#===========================================================================#
# plot epidemic curve animal first 20 days (initial spread )
plot_SEIR_animals(model_4, plot_suceptible_compartment = F, by_host = F)
## plot epidemic curve animal level by host
plot_SEIR_animals(model_4, plot_suceptible_compartment = F, by_host = T)
#plot the epidemic curve considering
plot_SEIR_farms(model_4)

# here will plot the control actions
myplots  <- plot_SEIR_vacc_depop_ani(modeloutput = model_4, summary_control = banco_summary)

#  1) animals vaccinated
myplots[1]
# 2) animals depopulated
myplots[2]
# 3) both vaccinated and depopulated animals
myplots[3]
