ABM simulating vaccination strategies against COVID-19. Presented on Work "Projecting the impact of Targeted Vaccination on COVID-19 Outbreaks in the United States". Authors:
Seyed M. Moghadas, Thomas N. Vilches, Kevin Zhang,3 Chad R Wells, Affan Shoukat, Meagan C. Fitzpatrick, Joanne M. Langley, Alison P. Galvani

This ABM is based on the one presented by Affan S. in the link
https://github.com/affans/covid19abm.jl

It was modified to include vaccination, comorbidity and a age-based herd immunity previously calibrated.
Vaccination was modified in Dec 2020 being applied while the epidemics happens

In other to run the simulations enter in Julia and use either
include("simulations.jl"), if you wont use a cluster, or
include("simulations_cluster.jl"), if you will use a cluster (note that it is done using "Cluster manager" package and slurm).
Warning: specify the path in the "create_folder" function

After including the above, run simulations using
run_param(beta,herd_im_v = [0],fs=0.0,init=3,vaccinate = false,ves1=0.0,vese1=0.0,vei1=0.0,ves2=0.0,vese2=0.0,vei2=0.0,vp = 28,sc = false,cov = 0.0,nsims=1000)

in which the arguments are
b = the probability transmission (the default is 0.08 which generates 60% of Attack Rate without any kind of control)
herd_im_v  = integer vector with the percentage of herd immunity one wants to run. The possible values are 0, 5, 10, 20.
fs = proportion of severe cases that are isolated
init = number of initial infections
vaccinate = Boolean saying if one wants to apply vac or not. (default is false)
ves* = efficacy against symptoms for 1- first dose, and 2- second dose
vese* = efficacy against SEVERE symptoms for 1- first dose, and 2- second dose
vei* = efficacy against infection for 1- first dose, and 2- second dose
vp = period between 1st and 2nd doses
sc = set specific coverage? Boolean. The default is false, which generates 40% of coverage (when vaccinating).
cov = the specific coverage
nsims = integer. Number of monte carlo simulations to be performed

The age distribution is set to USA reports. The Vaccination is implemented in order to vaccinate HCW first, then comorbid and elderly people, and last general population (all of them over 18 years old)

Special attention to parameter fsevere (fs in the function run_param). It is the probability that a severe individual isolate themselve, limiting their contacts to 3 per day. 

The output is a folder that is named "results_" and the combination of vac_ef, herd_immunity and days_b.

"simlevel" files have the temporal dynamics of individual class they indicates, e.g. "simlevel_ded_inc_ag3.dat" stands for incidence of deads in age group 3.
hos = hospitalization, icu = icu hospitalization, lat = latent. It has a "heading"

File "general_vac_info.dat" gives us, for each simulation the total number of
vaccines given to susceptibles, vaccines given to recovered (herd immunity), infected individuals that were vaccinated, infected individuals that were not vaccinated, dead individuals that were vaccinated, dead individuals that were not vaccinated, non-icu hospitalized individuals that were vaccinated, non-icu hospitalized individuals that were not vaccinated, icu hospitalized individuals that were vaccinated, icu hospitalized individuals that were not vaccinated

The file run_scen_review.jl runs all the needed scenarios.