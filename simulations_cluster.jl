using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles

## load the packages by covid19abm

#using covid19abm

#addprocs(2, exeflags="--project=.")


#@everywhere using covid19abm

addprocs(SlurmManager(500), N=17, topology=:master_worker, exeflags="--project=.")
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("covid19abm.jl")
@everywhere const cv=covid19abm


function run(myp::cv.ModelParameters, nsims=1000, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    dump(myp)
    myp.calibration && error("can not run simulation, calibration is on.")
    # will return 6 dataframes. 1 total, 4 age-specific 
    cdr = pmap(1:nsims) do x                 
            cv.runsim(x, myp)
    end      

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")
    ## write the infectors 
    DelimitedFiles.writedlm("$(folderprefix)/infectors.dat", [cdr[i].infectors for i = 1:nsims])    

    ## write contact numbers
    #writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    ag6 = vcat([cdr[i].g6 for i = 1:nsims]...)
    mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6)
    #mydfs = Dict("all" => allag)
    
    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild 
    #c1 = Symbol.((:LAT, :ASYMP, :INF, :IISO, :HOS, :ICU, :DED), :_INC)
    #c2 = Symbol.((:LAT, :ASYMP, :INF, :IISO, :HOS, :ICU, :DED), :_PREV)
    c1 = Symbol.((:LAT, :HOS, :ICU, :DED), :_INC)
    c2 = Symbol.((:LAT, :HOS, :ICU, :DED), :_PREV)
   for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        for c in vcat(c1..., c2...)
        #for c in vcat(c1...)
        #for c in vcat(c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        #yaf = compute_yearly_average(df)       
        #fn = string("$(folderprefix)/timelevel_", k, ".dat")   
        #CSV.write(fn, yaf)       
    end

    ##########Saving vaccine status
    #=println("saving vac status")
    vac_file = string(folderprefix,"/vac_status.dat")
    fv = open(vac_file,"w")
   
    for i = 1:nsims
        for j = 1:length(cdr[1].vi)
            print(fv,"$(cdr[i].vi[j]) ")
        end
        println(fv,"")
    end

    close(fv)=#
    #############saving vac ef
    println("saving vac effcicacy")
    vac_ef_file = string(folderprefix,"/vac_ef.dat")
    fve = open(vac_ef_file,"w")
   
    for i = 1:nsims
        for j = 1:length(cdr[1].ve)
            print(fve,"$(cdr[i].ve[j]) ")
        end
        println(fve,"")
    end

    close(fve)

    #######saving commorbidity status
    #=println("saving commorbidity status")
    com_file = string(folderprefix,"/com_idx.dat")
    fc = open(com_file,"w")
    
    for i = 1:nsims
        for j = 1:length(cdr[1].com)
            print(fc,"$(cdr[i].com[j]) ")
        end
        println(fc,"")
    end
    close(fc)=#

    ########## save general info about vaccine
    n_vac_sus1 = [cdr[i].n_vac_sus1 for i=1:nsims]
    n_vac_rec1 = [cdr[i].n_vac_rec1 for i=1:nsims]
    n_inf_vac1 = [cdr[i].n_inf_vac1 for i=1:nsims]
    n_dead_vac1 = [cdr[i].n_dead_vac1 for i=1:nsims]
    n_hosp_vac1 = [cdr[i].n_hosp_vac1 for i=1:nsims]
    n_icu_vac1 = [cdr[i].n_icu_vac1 for i=1:nsims]
    
    n_vac_sus2 = [cdr[i].n_vac_sus2 for i=1:nsims]
    n_vac_rec2 = [cdr[i].n_vac_rec2 for i=1:nsims]
    n_inf_vac2 = [cdr[i].n_inf_vac2 for i=1:nsims]
    n_dead_vac2 = [cdr[i].n_dead_vac2 for i=1:nsims]
    n_hosp_vac2 = [cdr[i].n_hosp_vac2 for i=1:nsims]
    n_icu_vac2 = [cdr[i].n_icu_vac2 for i=1:nsims]

    n_dead_nvac = [cdr[i].n_dead_nvac for i=1:nsims]
    n_inf_nvac = [cdr[i].n_inf_nvac for i=1:nsims]
    n_hosp_nvac = [cdr[i].n_hosp_nvac for i=1:nsims]
    n_icu_nvac = [cdr[i].n_icu_nvac for i=1:nsims]
    R0 = [cdr[i].R0 for i=1:nsims]
    #=
    writedlm(string(folderprefix,"/general_vac_info_cov_$(replace(string(myp.cov_val), "." => "_"))_vac_ef_$(replace(string(myp.vaccine_ef), "." => "_")).dat"),[n_vac_sus1 n_vac_rec1 n_inf_vac1 n_dead_vac1 n_hosp_vac1 n_icu_vac1 n_vac_sus2 n_vac_rec2 n_inf_vac2 n_dead_vac2 n_hosp_vac2 n_icu_vac2 n_inf_nvac n_dead_nvac n_hosp_nvac n_icu_nvac])
    writedlm(string(folderprefix,"/com_vac_cov_$(replace(string(myp.cov_val), "." => "_"))_vac_ef_$(replace(string(myp.vaccine_ef), "." => "_")).dat"),[cdr[i].com_v for i=1:nsims])
    writedlm(string(folderprefix,"/com_total_cov_$(replace(string(myp.cov_val), "." => "_"))_vac_ef_$(replace(string(myp.vaccine_ef), "." => "_")).dat"),[cdr[i].com_t for i=1:nsims])
    writedlm(string(folderprefix,"/ncom_vac_cov_$(replace(string(myp.cov_val), "." => "_"))_vac_ef_$(replace(string(myp.vaccine_ef), "." => "_")).dat"),[cdr[i].ncom_v for i=1:nsims])
    writedlm(string(folderprefix,"/ncom_total_cov_$(replace(string(myp.cov_val), "." => "_"))_vac_ef_$(replace(string(myp.vaccine_ef), "." => "_")).dat"),[cdr[i].ncom_t for i=1:nsims])
    =#
    data = DataFrame(vac_sus_dose1 = n_vac_sus1,vac_herd_dose_1 = n_vac_rec1,inf_dose_1 = n_inf_vac1, dead_dose_1 = n_dead_vac1, hosp_dose_1 = n_hosp_vac1,icu_dose_1 = n_icu_vac1, vac_sus_dose_2 = n_vac_sus2, vac_herd_dose_2 = n_vac_rec2, inf_dose_2 = n_inf_vac2, dead_dose_2 = n_dead_vac2, hosp_dose_2 = n_hosp_vac2, icu_dose_2 = n_icu_vac2, inf_n_vac = n_inf_nvac,dead_n_vac = n_dead_nvac,hosp_n_vac = n_hosp_nvac,icu_n_vac = n_icu_nvac)
    
    #writedlm(string(folderprefix,"/general_vac_info.dat"),data)
    CSV.write("$folderprefix/general_vac_info.csv",data)
    writedlm(string(folderprefix,"/com_vac1.dat"),[cdr[i].com_v1 for i=1:nsims])
    writedlm(string(folderprefix,"/ncom_vac1.dat"),[cdr[i].ncom_v1 for i=1:nsims])
    writedlm(string(folderprefix,"/com_vac2.dat"),[cdr[i].com_v2 for i=1:nsims])
    writedlm(string(folderprefix,"/ncom_vac2.dat"),[cdr[i].ncom_v2 for i=1:nsims])
    writedlm(string(folderprefix,"/com_total.dat"),[cdr[i].com_t for i=1:nsims])
    writedlm(string(folderprefix,"/ncom_total.dat"),[cdr[i].ncom_t for i=1:nsims])
    writedlm(string(folderprefix,"/R0.dat"),R0)
    
    writedlm(string(folderprefix,"/init_iso.dat"),[cdr[i].iniiso for i=1:nsims])

    return mydfs
end


function compute_yearly_average(df)
    ya = df |> @groupby(_.time) |> @map({time=key(_), cnt=length(_),
              sus_prev=mean(_.SUS_PREV), 
              lat_prev=mean(_.LAT_PREV), 
              pre_prev=mean(_.PRE_PREV), 
              asymp_prev=mean(_.ASYMP_PREV), 
              mild_prev=mean(_.MILD_PREV), 
              miso_prev=mean(_.MISO_PREV), 
              inf_prev=mean(_.INF_PREV), 
              iiso_prev=mean(_.IISO_PREV), 
              hos_prev=mean(_.HOS_PREV), 
              icu_prev=mean(_.ICU_PREV), 
              rec_prev=mean(_.REC_PREV), 
              ded_prev=mean(_.DED_PREV), 
              sus_inc=mean(_.SUS_INC),
              lat_inc=mean(_.LAT_INC), 
              pre_inc=mean(_.PRE_INC), 
              asymp_inc=mean(_.ASYMP_INC), 
              mild_inc=mean(_.MILD_INC), 
              miso_inc=mean(_.MISO_INC), 
              inf_inc=mean(_.INF_INC),
              iiso_inc=mean(_.IISO_INC),
              hos_inc=mean(_.HOS_INC),
              icu_inc=mean(_.ICU_INC),
              rec_inc=mean(_.REC_INC),
              ded_inc=mean(_.DED_INC)
              }) |> DataFrame
    return ya
end
#=
function savestr(p::cv.ModelParameters, custominsert="/", customstart="")
    datestr = (Dates.format(Dates.now(), dateformat"mmdd_HHMM"))
    ## setup folder name based on model parameters
    taustr = replace(string(p.τmild), "." => "")
    fstr = replace(string(p.fmild), "." => "")
    rstr = replace(string(p.β), "." => "")
    prov = replace(string(p.prov), "." => "")
    eldr = replace(string(p.eldq), "." => "")
    eldqag = replace(string(p.eldqag), "." => "")     
    fpreiso = replace(string(p.fpreiso), "." => "")
    tpreiso = replace(string(p.tpreiso), "." => "")
    fsev = replace(string(p.fsevere), "." => "")    
    frelasymp = replace(string(p.frelasymp), "." => "")
    strat = replace(string(p.ctstrat), "." => "")
    pct = replace(string(p.fctcapture), "." => "")
    cct = replace(string(p.fcontactst), "." => "")
    idt = replace(string(p.cidtime), "." => "") 
    tback = replace(string(p.cdaysback), "." => "")     
    fldrname = "/data/covid19abm/simresults/$(custominsert)/$(customstart)_$(prov)_strat$(strat)_pct$(pct)_cct$(cct)_idt$(idt)_tback$(tback)_fsev$(fsev)_tau$(taustr)_fmild$(fstr)_q$(eldr)_qag$(eldqag)_relasymp$(frelasymp)_tpreiso$(tpreiso)_preiso$(fpreiso)/"
    mkpath(fldrname)
end=#

function _calibrate(nsims, myp::cv.ModelParameters)
    myp.calibration != true && error("calibration parameter not turned on")
    vals = zeros(Int64, nsims)
    println("calibrating with beta: $(myp.β), total sims: $nsims, province: $(myp.prov)")
    println("calibration parameters:")
    dump(myp)
    cdr = pmap(1:nsims) do i 
        h,hh = cv.main(myp,i) ## gets the entire model. 
        val = sum(cv._get_column_incidence(h, covid19abm.LAT))            
        return val
    end
    return mean(cdr), std(cdr)
end

function calibrate(beta, nsims, herdi = 0, cali2 = false, fs = 0.0, prov=:usa, init_inf=1, size=10000)
    myp = cv.ModelParameters() # set up default parameters 
    myp.β = beta
    myp.prov = prov
    myp.popsize = size
    myp.modeltime = 30
    myp.calibration = true
    myp.calibration2 = cali2
    myp.fsevere = fs
    myp.initialinf = init_inf
    myp.herd = herdi
    m, sd = _calibrate(nsims, myp)
    println("mean R0: $(m) with std: $(sd)")
    myp.calibration = false       
    return m
end

function calibrate_robustness(beta, reps, prov=:usa)
    #[:ontario, :alberta, :bc, :manitoba, :newbruns, :newfdland, :nwterrito, :novasco, :nunavut, :pei, :quebec, :saskat, :yukon]
    # once a beta is found based on nsims simulations, 
    # see how robust it is. run calibration with same beta 100 times 
    # to see the variation in R0 produced. 
    nsims = [500, 1000]
    means = zeros(Float64, reps, length(nsims))
    for (i, ns) in enumerate(nsims)
        cd = map(1:reps) do x 
            println("iter: $x, sims: $ns")
            mval = calibrate(beta, ns, prov)         
            return mval
        end
        means[:, i] = cd
    end
    # for i in 2:nworkers()
    #     ## mf defined as: @everywhere mg() = covid19abm.p.β     
    #     rpr = remotecall_fetch(mf,  i+1).prov
    #     rpr != prov && error("province didn't get set in the remote workers")
    # end
    return means
end

function create_folder(ip::cv.ModelParameters)
    #strategy = ip.apply_vac_com == true ? "T" : "UT"
    strategy = ip.vaccinating == true ? "$(ip.days_before)" : "NV"
    #RF = string("heatmap/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vaccine_ef), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy","cov_$(replace(string(ip.cov_val)))") ## 
    main_folder = "/data/thomas-covid/new_vac_usa"
    if ip.set_g_cov
        RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vac_efficacy), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy","_cov_$(replace(string(ip.cov_val), "." => "_"))") ## 
    else
        RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vac_efficacy), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy") ##  
    end
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end

## now, running vaccine and herd immunity, focusing and not focusing in comorbidity, first  argument turns off vac
function run_param(b,herd_im_v = [0],fs=0.0,vaccinate = false,days_b = [0],vac_ef_v = [0.0],nsims=1000)
    for v_e = vac_ef_v, h_i = herd_im_v,days_b1 = days_b
        #bd = Dict(5=>0.074, 10=>0.076, 20=>0.089)
        #b = bd[h_i]
        @everywhere ip = cv.ModelParameters(β=$b,fsevere = $fs,vaccinating = $vaccinate, days_before = $days_b1,vac_efficacy = $v_e,herd = $(h_i))
        folder = create_folder(ip)

        #println("$v_e $(ip.vaccine_ef)")
        run(ip,nsims,folder)
    end
end


## now, running vaccine and herd immunity, focusing and not focusing in comorbidity, first  argument turns off vac
function run_param_fix(herd_im_v = [0],fs=0.0,vaccinate = false,days_b = [0],vac_ef_v = [0.0],sc = false,cov = 0.0,nsims=1000)
    for v_e = vac_ef_v, h_i = herd_im_v,days_b1 = days_b
        bd = Dict(3=>0.0425,5=>0.0425, 10=>0.0448, 20=>0.0515)
        b = bd[h_i]
        @everywhere ip = cv.ModelParameters(β=$b,fsevere = $fs,vaccinating = $vaccinate, days_before = $days_b1,vac_efficacy = $v_e,herd = $(h_i),set_g_cov = $sc,cov_val = $cov)
        folder = create_folder(ip)

        #println("$v_e $(ip.vaccine_ef)")
        run(ip,nsims,folder)
    end
end


## now, running vaccine and herd immunity, focusing and not focusing in comorbidity, first  argument turns off vac
function run_calibration(beta = 0.0345,herd_im_v = 0,cali = false,fs = 0.0,nsims=1000)
   
    @everywhere ip = cv.ModelParameters(β=$beta,herd = $herd_im_v,modeltime = 30,fsevere = $fs,calibration2 = $cali)
    folder = create_folder(ip)

    #println("$v_e $(ip.vaccine_ef)")
    run(ip,nsims,folder)

    R0 = readdlm(string(folder,"/R0.dat"),header=false)[:,1]
    m = mean(R0)
    sd = std(R0)
    println("mean R0: $(m) with std: $(sd)")
end