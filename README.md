# Workshop: Use of transmission models to simulate the spread of livestock diseases [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)



<a href="url"><img src="https://github.com/ncespedesc/logos_nc_state/blob/main/MHASpread_logo.png?raw=true" align="left" height="150" width="150" ></a>

## :mortar_board: About this workshop

> In this four-day workshop, you will have an introduction to a range of mathematical models used to simulate the spread of livestock diseases. We will focus on the application of such epidemiological models and demonstrate with real data, how you can use mathematical transmission models to make informed decisions before, during, and after an animal health emergency.  
We use the MHASpread: A multi-host animal spread stochastic multilevel model (version 0.1.0) which is an R package to be used throughout the training. The MHASpread allows for explicit specification of species-specific disease transmission probability, among other important transmission dynamics of disease infecting multiple species, such as FMD. This model considers the entry and exit of animals given between-farm animal movements, movements into slaughterhouses, births, and, deaths, for each species. 
You will learn how to use MHASpread, including the simulation of the introduction, and dissemination of FMD in the state of Rio Grande do Sul, Brazil. You will have access to highly specialized computational and epidemiological tools within an easy-to-use workflow. 
For the second half of the workshop, you will learn how similar models are used in the preparation of ASF in Rio Grande do Sul, Brazil, and in the United States. We will demonstrate how the “PigSpread” model works. PigSpread is also a mathematical model specially developed to be used in the dynamics of the disease of swine.

## :bomb: Aims of the workshop
* Learn how to use the MHASpread v.0.1.0 package, introduction, and control of FMD.
    - [ ]  Overview of the model outputs and their interpretation.
    - [ ]  MHASpread to simulate FMD countermeasure actions (depopulation, vaccination, and traceability). 
    - [ ]  To be exposed to additional transmission models.

## Pre workshop activities 

| **Topic**                          | **Activity and assignment**                                                                                                       | **Date** |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|----------|
| [R & RStudio](https://github.com/machado-lab/workshop_MHASpread/tree/main/install_RStudio_and_packages)| Video call to set Rstudio in your computer [_link to call here:_](https://ncsu.zoom.us/j/93057927377?pwd=V2diU0ZaQmJjdWlTaUJGNXZFUG0zZz09)                                                                                         | Sep/12/22 |
| Introduction to disease spread     | [A Practical Introduction to Mechanistic Modeling of Disease Transmission in Veterinary Science](https://doi.org/10.3389/fvets.2020.546651)                                    | Sep/12/22 |
| Why use models                     | [Three questions to ask before using model outputs for decision support](https://doi.org/10.1038/s41467-020-17785-2)               | Sep/26/22 |
| Use of data to prepare against FMD | [Challenges and opportunities for using national animal datasets to support foot-and-mouth disease control]( https://doi.org/10.1111/tbed.13858)                         | Oct/01/22 |
| Main transmission routes           | [Understanding the transmission of foot-and-mouth disease virus at different scales](https://doi.org/10.1016/j.coviro.2017.11.013) | Oct/10/22 |


## :calendar: Calendar 

| **Day-1 (October 17)**                                                                                                                          | Topics                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|---------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Review of FMD modeling outbreaks, the application of the MHASpread R package. We will use actual population and movement data *Morning      | Welcome & introductions Importance of the response plan in outbreak events. How simulation models can help the design and update of control strategies. How to integrate the model outputs with real case scenarios. The whole game: why R and Rstudio environment.                                                                                                                                                                                                                                                                                                                |
| Introduction to Rstudio installing R/Rstudio and MHSpread  *Afternoon                                                                       | Installation of R, Rstudio, packages. Setting up your computer system.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| **Day-2 (October 18)**                                                                                                                          |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Hands-on MHASpread  & Morning                                                                                                               | MHASpread architecture. Load the package and set the environment. Define initial conditions. Select the originally infected farm. Export output numeric summaries and plots.                                                                                                                                                                                                                                                                                                                                                                                                       |
| Simulating FMD epidemics within the state of Rio Grande do Sul.  *Afternoon                                                                 | Select the initial condition  Seeding infections and reconstruct spread along with epidemic sizes. Simulate control and eradication actions biweek. Export output numeric summaries and plots. MHASpread is integrated with the Rapid Access Biosecurity (RAB) app to expedite decision-making.                                                                                                                                                                                                                                                                                    |
| **Day-3 (October 19)**                                                                                                                          |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| The simulation of FMD epidemics with index cases in swine, cattle, and small ruminants. & Morning and *Afternoon                            | Simulate FMD introduction and spread from swine, bovine, and small ruminants populations  Estimate epidemic trajectories. Implementation of control and eradication actions and proposed next steps for elimination. Extract model outputs and data into RABapp™ Simulation of alternative control scenarios Enlarge control zone size. Traceback in the contact networks.  Increase the duration of control zones and surveillance activities. Demonstrate the real model repeats (from the code stack). Simulate scenarios that will take place next week in the field exercise. |
| **Day-4 (October 20)** Afternoon                                                                                                                |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Demonstration of other transmission models The Rio Grande do Sul, ASF transmission model PigSpread ASF model U.S. PigSpread PRRS model U.S. | We will demonstrate each model with some opportunity for hands-on interaction with these models.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Question round other applications   & Morning                                                                                                        | Questions and future remarks                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |


### Usefull links
- [ ]  [How to install R and Rstudio](https://github.com/machado-lab/workshop_MHASpread/blob/main/install_RStudio_and_packages/README.md) 
- [ ]  [Necessary packages to run the model](https://github.com/machado-lab/workshop_MHASpread/blob/main/install_RStudio_and_packages/package_to_install.R) 

## Authors
Nicolas Cespedes Cardenas [![ORCIDiD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-7884-2353)

Gustavo Machado [![ORCIDiD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-7552-6144)

## :computer: Website
[MachadoLAb](https://machado-lab.github.io/) 

## :muscle: Sponsors
<a href="url"><img src="https://github.com/ncespedesc/logos_nc_state/blob/main/fundesalogo.jpg?raw=true" align="left" width="200" ></a>

<a href="url"><img src="https://github.com/ncespedesc/logos_nc_state/blob/main/ncstate-type-4x1-red-min.png?raw=true" align="left" width="200" ></a>


<a href="url"><img src="https://github.com/ncespedesc/logos_nc_state/blob/main/seapilogo.png?raw=true" align="left" width="300" ></a>


