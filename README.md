ABM simulating vaccination strategies against COVID-19. Presented on Work "Projecting the impact of Targeted Vaccination on COVID-19 Outbreaks in the United States". Authors:
Seyed M. Moghadas, Thomas N. Vilches, Kevin Zhang,3 Chad R Wells, Affan Shoukat, Meagan C. Fitzpatrick, Joanne M. Langley, Alison P. Galvani

This ABM is based on the one presented by Affan S. in the link
https://github.com/affans/covid19abm.jl

It was modified to include vaccination, comorbidity and a age-based herd immunity previously calibrated.

In other to run the simulations enter in Julia and use either
include("simulations.jl"), if you wont use a cluster, or
include("simulations_cluster.jl"), if you will use a cluster (note that it is done using "Cluster manager" package and slurm).

After including the above, run simulations using
run_param(beta,ap_vac,vac_ef_v,vac_com_v,herd_im_v)

in which the arguments are
beta = the probability transmission (the default is 0.08 which generates 60% of Attack Rate without any kind of control)
ap_vac = Boolean saying if one wants to apply vac or not. (default is false)
vac_ef_v = float vector with vaccines efficacy one wants to run. (default is 0.0)
vac_com_v = Boolean vector saying if you want to focus vaccination on comorbidity of not.
herd_im_v  = integer vector with the percentage of herd immunity one wants to run. The possible values are 0, 5, 10, 20.
nsims = integer. Number of monte carlo simulations to be performed

The age distribution and vaccine coverage is set to USA reports.
In file covid19abm.jl, one can change the parameter set_g_cov to true, which ignores the coverage reported in USA (48% of population) and use a general coverage given in parameter cov_val (default = 0.8), distributing vaccines randomly (uniform distribution) to population.

Special attention to parameter fsevere. It is the probability that a severe individual isolate themselve, limiting their contacts to 3 per day. The default value is 1.0, but calibration must be performed to fsevere = 0.0.

The output is a folder that is named "results_" and the combination of vac_ef, herd_immunity and com_vac. when vaccination focuses on comorbidity, the folder will have the S1, and S2 otherwise.

"simlevel" files have the temporal dynamics of individual class they indicates, e.g. "simlevel_ded_inc_ag3.dat" stands for incidence of deads in age group 3. "prev" means prevalence.
hos = hospitalization, icu = icu hospitalization, lat = latent. It has a "heading"

File "general_vac_info.dat" gives us, for each simulation the total number of
vaccines given to susceptibles, vaccines given to recovered (herd immunity), infected individuals that were vaccinated, infected individuals that were not vaccinated, dead individuals that were vaccinated, dead individuals that were not vaccinated, non-icu hospitalized individuals that were vaccinated, non-icu hospitalized individuals that were not vaccinated, icu hospitalized individuals that were vaccinated, icu hospitalized individuals that were not vaccinated