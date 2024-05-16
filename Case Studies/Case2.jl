using Pkg
Pkg.activate(".")
# Pkg.instantiate()

using JuMP
using Ipopt
using PowerModels
using Memento
using InfrastructureModels
using StochasticPowerModels
using PowerModelsACDC
using FlexPlan
using Plots
using Juniper
using XLSX

#Constants
const _PM = PowerModels
const _SPM = StochasticPowerModels
const _PMACDC = PowerModelsACDC
const _FP = FlexPlan

Memento.setlevel!(Memento.getlogger(StochasticPowerModels), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")

#Solver inputs
ipopt_solver1 = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>15000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-4,
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

ipopt_solver2 = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>15000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-6,
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver2, "atol" => 1e-6, "log_levels" => [])

#gPC degree input
deg  = 2



#RES inputs
pen_level_start = 0.70
pen_level_step = 0.05
pen_level_end = 0.70

#Report decision
Report_curt = true
global solve_SOTS = true
global solve_SOPF = false


#Necessary initializations
result_SOTS = Dict()
result_SOTS_bin = Dict()
result_SOPF = Dict()
result_ACDC_ots = Dict()
result_ACDC_opf = Dict()

obj_SOTS = Dict()
obj_SOTS_bin = Dict()
obj_SOPF = Dict()

stat_SOTS = Dict()
stat_SOTS_bin = Dict()
stat_SOPF = Dict()

time_SOTS = Dict()
time_SOTS_bin = Dict()
time_SOPF = Dict()

sample1 = Dict()
sample2 = Dict()

p_size_dict = Dict()



global feas_ctr1 = 0
global feas_ctr2 = 0
global feas_ctr3 = 0




#Case file and data reading
# case = "case5_ACDC_SPM_95cc.m"
case = "case67acdc_SPM_90cc_det.m"
case_det = "case67acdc_SPM_90cc_det.m"

file  = joinpath(BASE_DIR, "test/data/matpower", case)
file_det = joinpath(BASE_DIR, "test/data/matpower", case_det)


s = Dict("RES Curtailment" => true,
         "Load Curtailment" => false,
         )

#Penetration level iteration
for pen_level = pen_level_start:pen_level_step:pen_level_end
    
    data = _PM.parse_file(file)
    data_det_AC = _PM.parse_file(file_det) 
    data_ACDC = _PM.parse_file(file)
    data_OPF = _PM.parse_file(file)

    _PMACDC.process_additional_data!(data)
    _PMACDC.process_additional_data!(data_ACDC)
    _PMACDC.process_additional_data!(data_OPF)


    if solve_SOTS || solve_SOPF
        println("\n", "Penetration Level = $pen_level")

        global total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
        global p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus
        p_size_dict[pen_level] = p_size

        
        for i = 1:length(data_det_AC["load"])
            data_det_AC["load"]["$i"]["pd"] = data["load"]["$i"]["pd"] #- p_size #initial point belirlerken denemelik cikarttim pen_level'i
        end
        _PMACDC.process_additional_data!(data_det_AC)
         
        println("   Deterministic OTS - Solution progress: Solving...")
            result_ACDC_ots[pen_level] = _PMACDC.run_acdcots(data_det_AC, _PM.ACPPowerModel, minlp_solver)         
        println("   Deterministic OTS - Solution progress: Solved! (", string(result_ACDC_ots[pen_level]["primal_status"]), ")")
        
        println("   Deterministic OPF - Solution progress: Solving...")
            result_ACDC_opf[pen_level] = _PMACDC.run_acdcopf(data_det_AC, _PM.ACPPowerModel, minlp_solver)
        println("   Deterministic OPF - Solution progress: Solved! (", string(result_ACDC_opf[pen_level]["primal_status"]), ")")
         
        for (b, branch) in result_ACDC_ots[pen_level]["solution"]["branch"]

            if branch["br_status"] < 0.5
                data_ACDC["branch"][b]["br_status_initial"] = 0
                # data_OPF["branch"][b]["br_status"] = 0
            else
                data_ACDC["branch"][b]["br_status_initial"] = 1
            end
        end
            
        for (b, branch) in result_ACDC_ots[pen_level]["solution"]["branchdc"]
            
            if branch["br_status_dc"] < 0.5
                data_ACDC["branchdc"][b]["br_status_initial"] = 1
                # data_OPF["branchdc"][b]["status"] = 1
            else
                data_ACDC["branchdc"][b]["br_status_initial"] = 1
            end
            
        end

        if solve_SOTS

            sdata_OTS = _SPM.build_stochastic_data_ACDC_RES(data_ACDC, deg, p_size)

            PCE_data = Dict(
                  "T2" => sdata_OTS["T2"],
                  "T3" => sdata_OTS["T3"],
                  "T4" => sdata_OTS["T4"],
                  "mop" => sdata_OTS["mop"],
               )
            num_of_coeff = PCE_data["mop"].dim
            
            delete!(sdata_OTS, "T2")
            delete!(sdata_OTS, "T3")
            delete!(sdata_OTS, "T4")
            delete!(sdata_OTS, "mop")
            
            _FP.add_dimension!(sdata_OTS, :PCE_coeff, num_of_coeff; metadata = PCE_data)
            
            sdata_OTS["curtailment"] = s

            
            println("   SOTS - Solution progress: Solving...")
            global result_SOTS[pen_level] = _SPM.solve_sots_iv_acdc(sdata_OTS, _PM.IVRPowerModel, ipopt_solver1, deg=deg, p_size=p_size);
            
            #Store necessary values for reporting
            obj_SOTS[pen_level] = result_SOTS[pen_level]["objective"]
            stat_SOTS[pen_level] = string(result_SOTS[pen_level]["primal_status"])
            time_SOTS[pen_level] = result_SOTS[pen_level]["solve_time"]

            if sdata_OTS["curtailment"]["RES Curtailment"] == true
                for i = 1:num_of_coeff
                    result_SOTS[pen_level]["solution"]["nw"]["$i"]["RES"]["4"] = Dict()
                    result_SOTS[pen_level]["solution"]["nw"]["$i"]["RES"]["4"]["p_RES_curt"] = sum([result_SOTS[pen_level]["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:3])
                end
            end
            
            
            if string(result_SOTS[pen_level]["primal_status"]) != "FEASIBLE_POINT" && string(result_SOTS[pen_level]["primal_status"]) !="NEARLY_FEASIBLE_POINT"
                global feas_ctr1 += 1
            else
                global feas_ctr1 = 0
            end
            
            println("   SOTS - Solution progress: Solved! (", string(result_SOTS[pen_level]["primal_status"]), ")")

            if feas_ctr1 >= 2
                global solve_SOTS = false
                println("   SOTS reached infeasibility on penetration level of $pen_level.")
            end

            

        end

        if solve_SOPF

            for (b, branch) in result_SOTS[pen_level]["solution"]["nw"]["1"]["branch"]

                if branch["br_status"] < 0.5
                    data_OPF["branch"][b]["br_status"] = 0
                else
                    data_OPF["branch"][b]["br_status"] = 1
                end
            end

            sdata_OPF = _SPM.build_stochastic_data_ACDC_RES(data_OPF, deg, p_size)

            PCE_data = Dict(
                "T2" => sdata_OPF["T2"],
                "T3" => sdata_OPF["T3"],
                "T4" => sdata_OPF["T4"],
                "mop" => sdata_OPF["mop"],
            )
            num_of_coeff = PCE_data["mop"].dim

            delete!(sdata_OPF, "T2")
            delete!(sdata_OPF, "T3")
            delete!(sdata_OPF, "T4")
            delete!(sdata_OPF, "mop")

            _FP.add_dimension!(sdata_OPF, :PCE_coeff, num_of_coeff; metadata = PCE_data)

            sdata_OPF["curtailment"] = s

            println("   SOPF - Solution progress: Solving...")
            global result_SOPF[pen_level] = _SPM.solve_sopf_iv_acdc(sdata_OPF, _PM.IVRPowerModel, ipopt_solver1, deg=deg, p_size=p_size);
            
            #Store necessary values for reporting
            obj_SOPF[pen_level] = result_SOPF[pen_level]["objective"]
            stat_SOPF[pen_level] = string(result_SOPF[pen_level]["primal_status"])
            time_SOPF[pen_level] = result_SOPF[pen_level]["solve_time"]

            if sdata_OPF["curtailment"]["RES Curtailment"] == true
                for i = 1:num_of_coeff
                    result_SOPF[pen_level]["solution"]["nw"]["$i"]["RES"]["4"] = Dict()
                    result_SOPF[pen_level]["solution"]["nw"]["$i"]["RES"]["4"]["p_RES_curt"] = sum([result_SOPF[pen_level]["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:3])
                end
            end
            
            if string(result_SOPF[pen_level]["primal_status"]) != "FEASIBLE_POINT" && string(result_SOPF[pen_level]["primal_status"]) !="NEARLY_FEASIBLE_POINT"
                global feas_ctr3 += 1
            else
                global feas_ctr3 = 0
            end
            
            println("   SOPF - Solution progress: Solved! (", string(result_SOPF[pen_level]["primal_status"]), ")")

            if feas_ctr3 >= 2
                global solve_SOPF = false
                println("   SOPF reached infeasibility on penetration level of $pen_level.")
            end

        end

    end

end


# Reporting
if Report_curt

    if solve_SOTS

        # for i = 1:102
        #     sample_SOTS_branch[:,i] = _SPM.sample(result_SOTS[pen_level_start], "branch", i, "cr_fr"; sample_size=10000); 
        # end

        # for i = 1:12
        #     sample_SOTS_branchdc[:,i] = _SPM.sample(result_SOTS[pen_level_start], "branchdc", i, "pf"; sample_size=10000); 
        # end


        # sample1 = _SPM.sample(result_SOTS[pen_level_start], "branch", 14, "cr_fr"; sample_size=10000);    
        # sample2 = _SPM.sample(result_SOTS[pen_level_start], "branch", 15, "cr_fr"; sample_size=10000); 
        # sample3 = _SPM.sample(result_SOTS[pen_level_start], "branch", 17, "cr_fr"; sample_size=10000); 
        # sample4 = _SPM.sample(result_SOTS[pen_level_start], "branch", 32, "cr_fr"; sample_size=10000); 
        # sample5 = _SPM.sample(result_SOTS[pen_level_start], "branch", 41, "cr_fr"; sample_size=10000); 
        # sample6 = _SPM.sample(result_SOTS[pen_level_start], "branch", 52, "cr_fr"; sample_size=10000); 
        # sample7 = _SPM.sample(result_SOTS[pen_level_start], "branch", 69, "cr_fr"; sample_size=10000); 
        # sample8 = _SPM.sample(result_SOTS[pen_level_start], "branch", 71, "cr_fr"; sample_size=10000); 

        # sample9 = _SPM.sample(result_SOTS[pen_level_start], "branch", 14, "cmss"; sample_size=10000);    
        # sample10 = _SPM.sample(result_SOTS[pen_level_start], "branch", 15, "cmss"; sample_size=10000); 
        # sample11 = _SPM.sample(result_SOTS[pen_level_start], "branch", 17, "cmss"; sample_size=10000); 
        # sample12 = _SPM.sample(result_SOTS[pen_level_start], "branch", 32, "cmss"; sample_size=10000); 
        # sample13 = _SPM.sample(result_SOTS[pen_level_start], "branch", 41, "cmss"; sample_size=10000); 
        # sample14 = _SPM.sample(result_SOTS[pen_level_start], "branch", 52, "cmss"; sample_size=10000); 
        # sample15 = _SPM.sample(result_SOTS[pen_level_start], "branch", 69, "cmss"; sample_size=10000); 
        # sample16 = _SPM.sample(result_SOTS[pen_level_start], "branch", 71, "cmss"; sample_size=10000); 

        # sample17 = _SPM.sample(result_SOTS[pen_level_start], "RES", 3, "p_RES_curt"; sample_size=10000);
        # sample18 = _SPM.sample(result_SOTS[pen_level_start], "RES", 4, "p_RES_curt"; sample_size=10000);

        # sample19 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 11, "pf"; sample_size=10000);

        # sample20 = _SPM.sample(result_SOTS[pen_level_start], "gen", 10, "pg"; sample_size=10000);
        # sample21 = _SPM.sample(result_SOTS[pen_level_start], "gen", 11, "pg"; sample_size=10000);

        sample1 = _SPM.sample(result_SOTS[pen_level_start], "branch", 2, "cr_fr"; sample_size=10000);    
        sample2 = _SPM.sample(result_SOTS[pen_level_start], "branch", 3, "cr_fr"; sample_size=10000); 
        sample3 = _SPM.sample(result_SOTS[pen_level_start], "branch", 14, "cr_fr"; sample_size=10000); 
        sample4 = _SPM.sample(result_SOTS[pen_level_start], "branch", 15, "cr_fr"; sample_size=10000); 
        sample5 = _SPM.sample(result_SOTS[pen_level_start], "branch", 16, "cr_fr"; sample_size=10000); 
        sample6 = _SPM.sample(result_SOTS[pen_level_start], "branch", 93, "cr_fr"; sample_size=10000); 

        sample7 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 1, "pf"; sample_size=10000); 
        sample8 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 2, "pf"; sample_size=10000); 
        sample9 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 3, "pf"; sample_size=10000);    
        sample10 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 4, "pf"; sample_size=10000); 
        sample11 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 6, "pf"; sample_size=10000); 
        sample12 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 7, "pf"; sample_size=10000); 
        sample13 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 8, "pf"; sample_size=10000); 
        sample14 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 9, "pf"; sample_size=10000); 
        sample15 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 11, "pf"; sample_size=10000); 
        sample16 = _SPM.sample(result_SOTS[pen_level_start], "branchdc", 12, "pf"; sample_size=10000); 


    
        file_name1 = "Results\\Results - SOTS (1e-2) - (atol = 1e-4) - $case ($pen_level_start).xlsx"
        fid    = XLSX.openxlsx(file_name1, mode="w")
        header = ["1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9"; "10";
                    "11"; "12"; "13"; "14"; "15";
                    "16"]
        data   = [[sample1]; [sample2]; [sample3]; [sample4]; [sample5];
                    [sample6]; [sample7]; [sample8]; [sample9]; [sample10];
                    [sample11]; [sample12]; [sample13]; [sample14]; [sample15];
                    [sample16]]
        XLSX.writetable(file_name1, data, header)
    end

    if solve_SOPF

        # sample1 = _SPM.sample(result_SOPF[pen_level_start], "branch", 14, "cr_fr"; sample_size=10000);    
        # sample2 = _SPM.sample(result_SOPF[pen_level_start], "branch", 15, "cr_fr"; sample_size=10000); 
        # sample3 = _SPM.sample(result_SOPF[pen_level_start], "branch", 17, "cr_fr"; sample_size=10000); 
        # sample4 = _SPM.sample(result_SOPF[pen_level_start], "branch", 32, "cr_fr"; sample_size=10000); 
        # sample5 = _SPM.sample(result_SOPF[pen_level_start], "branch", 41, "cr_fr"; sample_size=10000); 
        # sample6 = _SPM.sample(result_SOPF[pen_level_start], "branch", 52, "cr_fr"; sample_size=10000); 
        # sample7 = _SPM.sample(result_SOPF[pen_level_start], "branch", 69, "cr_fr"; sample_size=10000); 
        # sample8 = _SPM.sample(result_SOPF[pen_level_start], "branch", 71, "cr_fr"; sample_size=10000); 

        # sample9 = _SPM.sample(result_SOPF[pen_level_start], "branch", 14, "cmss"; sample_size=10000);    
        # sample10 = _SPM.sample(result_SOPF[pen_level_start], "branch", 15, "cmss"; sample_size=10000); 
        # sample11 = _SPM.sample(result_SOPF[pen_level_start], "branch", 17, "cmss"; sample_size=10000); 
        # sample12 = _SPM.sample(result_SOPF[pen_level_start], "branch", 32, "cmss"; sample_size=10000); 
        # sample13 = _SPM.sample(result_SOPF[pen_level_start], "branch", 41, "cmss"; sample_size=10000); 
        # sample14 = _SPM.sample(result_SOPF[pen_level_start], "branch", 52, "cmss"; sample_size=10000); 
        # sample15 = _SPM.sample(result_SOPF[pen_level_start], "branch", 69, "cmss"; sample_size=10000); 
        # sample16 = _SPM.sample(result_SOPF[pen_level_start], "branch", 71, "cmss"; sample_size=10000); 

        # sample17 = _SPM.sample(result_SOPF[pen_level_start], "RES", 3, "p_RES_curt"; sample_size=10000);
        # sample18 = _SPM.sample(result_SOPF[pen_level_start], "RES", 4, "p_RES_curt"; sample_size=10000);

        # sample19 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 11, "pf"; sample_size=10000);

        # sample20 = _SPM.sample(result_SOPF[pen_level_start], "gen", 10, "pg"; sample_size=10000);
        # sample21 = _SPM.sample(result_SOPF[pen_level_start], "gen", 11, "pg"; sample_size=10000);

        sample1 = _SPM.sample(result_SOPF[pen_level_start], "branch", 2, "cr_fr"; sample_size=10000);    
        sample2 = _SPM.sample(result_SOPF[pen_level_start], "branch", 3, "cr_fr"; sample_size=10000); 
        sample3 = _SPM.sample(result_SOPF[pen_level_start], "branch", 14, "cr_fr"; sample_size=10000); 
        sample4 = _SPM.sample(result_SOPF[pen_level_start], "branch", 15, "cr_fr"; sample_size=10000); 
        sample5 = _SPM.sample(result_SOPF[pen_level_start], "branch", 16, "cr_fr"; sample_size=10000); 
        sample6 = _SPM.sample(result_SOPF[pen_level_start], "branch", 93, "cr_fr"; sample_size=10000); 

        sample7 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 1, "pf"; sample_size=10000); 
        sample8 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 2, "pf"; sample_size=10000); 
        sample9 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 3, "pf"; sample_size=10000);    
        sample10 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 4, "pf"; sample_size=10000); 
        sample11 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 6, "pf"; sample_size=10000); 
        sample12 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 7, "pf"; sample_size=10000); 
        sample13 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 8, "pf"; sample_size=10000); 
        sample14 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 9, "pf"; sample_size=10000); 
        sample15 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 11, "pf"; sample_size=10000); 
        sample16 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 12, "pf"; sample_size=10000); 

    
        file_name1 = "Results\\Results - SOPF (1e-8) - (atol = 1e-4) - $case ($pen_level_start).xlsx"
        fid    = XLSX.openxlsx(file_name1, mode="w")
        header = ["1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9"; "10";
                    "11"; "12"; "13"; "14"; "15";
                    "16"]
        data   = [[sample1]; [sample2]; [sample3]; [sample4]; [sample5];
                    [sample6]; [sample7]; [sample8]; [sample9]; [sample10];
                    [sample11]; [sample12]; [sample13]; [sample14]; [sample15];
                    [sample16]]
        XLSX.writetable(file_name1, data,header)
    end

end
    
    



