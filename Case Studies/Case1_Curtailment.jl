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

using Suppressor

Memento.setlevel!(Memento.getlogger(StochasticPowerModels), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")

#Solver inputs
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>10000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-8,
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "atol" => 1e-8, "log_levels" => [])
# atol 1e-4 iken MINLP local'e takildi

#gPC degree input
deg  = 2



#RES inputs
pen_level_start = 0.78
pen_level_step = 0.05
pen_level_end = 0.78

#Report decision
Report_curt = true
# Report_curt = false
global solve_SOTS = true
global solve_SOTS_bin = false
global solve_SOPF = true


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
case = "case5_ACDC_SPM_95cc.m"
file  = joinpath(BASE_DIR, "test/data/matpower", case)


s = Dict("RES Curtailment" => true,
         "Load Curtailment" => false,
         )

#Penetration level iteration
for pen_level = pen_level_start:pen_level_step:pen_level_end
    
    data = _PM.parse_file(file)
    data_det_AC = _PM.parse_file(file) 
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
                data_OPF["branch"][b]["br_status"] = 1
            else
                data_ACDC["branch"][b]["br_status_initial"] = 1
            end
        end
            
        for (b, branch) in result_ACDC_ots[pen_level]["solution"]["branchdc"]
            
            if branch["br_status_dc"] < 0.5
                data_ACDC["branchdc"][b]["br_status_initial"] = 1
                data_OPF["branchdc"][b]["status"] = 1
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

            if solve_SOTS
                println("   SOTS - Solution progress: Solving...")
                global result_SOTS[pen_level] = _SPM.solve_sots_iv_acdc(sdata_OTS, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
                
                #Store necessary values for reporting
                obj_SOTS[pen_level] = result_SOTS[pen_level]["objective"]
                stat_SOTS[pen_level] = string(result_SOTS[pen_level]["primal_status"])
                time_SOTS[pen_level] = result_SOTS[pen_level]["solve_time"]
                


                if sdata_OTS["curtailment"]["RES Curtailment"] == true
                    for i = 1:num_of_coeff
                        result_SOTS[pen_level]["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
                        result_SOTS[pen_level]["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_SOTS[pen_level]["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])
                    end
        
                    @suppress begin
                       global sample1[pen_level] = Dict()
                       for idx = 1:6
                          sample1[pen_level][idx] = _SPM.sample(result_SOTS[pen_level], "RES", idx, "p_RES_curt"; sample_size=10000);
                       end
                    end
                end
                
                if Report_curt
                    file_name1 = "Results\\Curtailment\\Curtailment - SOTS - $case - ($pen_level).xlsx"
                    fid    = XLSX.openxlsx(file_name1, mode="w")
                    header = ["RES 1";"RES 2";"RES 3";"RES 4";"RES 5";"RES Tot";]
                    data   = [[sample1[pen_level][1]];[sample1[pen_level][2]];[sample1[pen_level][3]];[sample1[pen_level][4]];[sample1[pen_level][5]];[sample1[pen_level][6]];]
                    XLSX.writetable(file_name1, data,header)
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

        end

        if solve_SOPF

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
            global result_SOPF[pen_level] = _SPM.solve_sopf_iv_acdc(sdata_OPF, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
            
            #Store necessary values for reporting
            obj_SOPF[pen_level] = result_SOPF[pen_level]["objective"]
            stat_SOPF[pen_level] = string(result_SOPF[pen_level]["primal_status"])
            time_SOPF[pen_level] = result_SOPF[pen_level]["solve_time"]
            
            if sdata_OPF["curtailment"]["RES Curtailment"] == true
                for i = 1:num_of_coeff
                    result_SOPF[pen_level]["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
                    result_SOPF[pen_level]["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_SOPF[pen_level]["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])
                end
    
                @suppress begin
                   global sample2[pen_level] = Dict()
                   for idx = 1:6
                    sample2[pen_level][idx] = _SPM.sample(result_SOPF[pen_level], "RES", idx, "p_RES_curt"; sample_size=10000);
                   end
                end
             end
             if Report_curt
                file_name1 = "Results\\Curtailment\\Curtailment - SOPF - $case - ($pen_level).xlsx"
                fid    = XLSX.openxlsx(file_name1, mode="w")
                header = ["RES 1";"RES 2";"RES 3";"RES 4";"RES 5";"RES Tot";]
                data   = [[sample2[pen_level][1]];[sample2[pen_level][2]];[sample2[pen_level][3]];[sample2[pen_level][4]];[sample2[pen_level][5]];[sample2[pen_level][6]];]
                XLSX.writetable(file_name1, data,header)
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


# println("\nPenetration Level = $pen_level")
# println("   Solution progress: Solving...")
# result_spm = _SPM.solve_sopf_iv_acdc(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
# println("   Solution progress: Solved! (", string(result_spm["primal_status"]), ")")



# # Show results on the terminal
# println("\n\n>>> SPMACDC Results >>>")
# println(result_spm["primal_status"])
# print("Objective: "
# print(result_spm["objective"])
# print("\nSolve Time: ")
# print(result_spm["solve_time"])



# total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
# p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus

# for i = 1:length(data_det_AC["load"])
#    data_det_AC["load"]["$i"]["pd"] = data["load"]["$i"]["pd"] - p_size
# end
# _PMACDC.process_additional_data!(data_det_AC)

# println("\n\n", "Solving Deterministic ACDC OTS Problem >>>>>", "\n")
#     result_ACDC_ots = _PMACDC.run_acdcots(data_det_AC, _PM.ACPPowerModel, minlp_solver)

# println("\n\n", "Solving Deterministic ACDC OPF Problem >>>>>", "\n")
#     result_ACDC_opf = _PMACDC.run_acdcopf(data_det_AC, _PM.ACPPowerModel, minlp_solver)

# println("\n\n>>> OTS Results >>>")
# println(result_ACDC_ots["primal_status"])
# print("Objective: ")
# print(result_ACDC_ots["objective"])
# print("\nSolve Time: ")
# print(result_ACDC_ots["solve_time"])

# println("\n\n>>> OPF Results >>>")
# println(result_ACDC_opf["primal_status"])
# print("Objective: ")
# print(result_ACDC_opf["objective"])
# print("\nSolve Time: ")
# print(result_ACDC_opf["solve_time"])


# for (b, branch) in result_ACDC_ots["solution"]["branch"]

#    if branch["br_status"] < 0.5
#        data_ACDC["branch"][b]["br_status_initial"] = 0
#        data_OPF["branch"][b]["br_status"] = 1
#    else
#        data_ACDC["branch"][b]["br_status_initial"] = 1
#    end
# end

# for (b, branch) in result_ACDC_ots["solution"]["branchdc"]

#    if branch["br_status_dc"] < 0.5
#        data_ACDC["branchdc"][b]["br_status_initial"] = 1
#        data_OPF["branchdc"][b]["status"] = 1
#    else
#        data_ACDC["branchdc"][b]["br_status_initial"] = 1
#    end

# end


# sdata_OPF = _SPM.build_stochastic_data_ACDC_RES(data_OPF, deg, p_size)

# PCE_data = Dict(
#       "T2" => sdata_OPF["T2"],
#       "T3" => sdata_OPF["T3"],
#       "T4" => sdata_OPF["T4"],
#       "mop" => sdata_OPF["mop"],
#    )
# num_of_coeff = PCE_data["mop"].dim

# delete!(sdata_OPF, "T2")
# delete!(sdata_OPF, "T3")
# delete!(sdata_OPF, "T4")
# delete!(sdata_OPF, "mop")

# _FP.add_dimension!(sdata_OPF, :PCE_coeff, num_of_coeff; metadata = PCE_data)

# sdata_OPF["curtailment"] = s

# println("\n\nPenetration Level = $pen_level")
# println("   Solution progress: Solving...")
# result_SOPF = _SPM.solve_sopf_iv_acdc(sdata_OPF, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
# println("   Solution progress: Solved! (", string(result_SOPF["primal_status"]), ")")



# sdata_OTS = _SPM.build_stochastic_data_ACDC_RES(data_ACDC, deg, p_size)

# PCE_data = Dict(
#       "T2" => sdata_OTS["T2"],
#       "T3" => sdata_OTS["T3"],
#       "T4" => sdata_OTS["T4"],
#       "mop" => sdata_OTS["mop"],
#    )
# num_of_coeff = PCE_data["mop"].dim

# delete!(sdata_OTS, "T2")
# delete!(sdata_OTS, "T3")
# delete!(sdata_OTS, "T4")
# delete!(sdata_OTS, "mop")

# _FP.add_dimension!(sdata_OTS, :PCE_coeff, num_of_coeff; metadata = PCE_data)

# sdata_OTS["curtailment"] = s

# println("\n\nPenetration Level = $pen_level")
# println("   Solution progress: Solving...")
# result_SOTS = _SPM.solve_sots_iv_acdc_deneme(sdata_OTS, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
# println("   Solution progress: Solved! (", string(result_SOTS["primal_status"]), ")")


# # println("\n   (MINLP) Solution progress: Solving...")
# # result_SOTS_bin = _SPM.solve_sots_iv_acdc_binary(sdata_OTS, _PM.IVRPowerModel, minlp_solver, deg=deg, p_size=p_size);
# # println("   (MINLP) Solution progress: Solved! (", string(result_spm["primal_status"]), ")")







# println("\n\n>>> SPMACDC SOPF Results >>>")
# println(result_SOPF["primal_status"])
# print("Objective: ")
# print(result_SOPF["objective"])
# print("\nSolve Time: ")
# print(result_SOPF["solve_time"])

# println("\n\n>>> SPMACDC SOTS Results >>>")
# println(result_SOTS["primal_status"])
# print("Objective: ")
# print(result_SOTS["objective"])
# print("\nSolve Time: ")
# print(result_SOTS["solve_time"])

# # println("\n\n>>> SPMACDC SOTS (MINLP) Results >>>")
# # println(result_SOTS_bin["primal_status"])
# # print("Objective: ")
# # print(result_SOTS_bin["objective"])
# # print("\nSolve Time: ")
# # print(result_SOTS_bin["solve_time"])

# bus = 6

# sample1 = _SPM.sample(result_SOTS, "branch", bus, "cmss"; sample_size=100000); 
# sample2 = _SPM.sample(result_SOPF, "branch", bus, "cmss"; sample_size=100000); 


# Plots.histogram(sample1)

# Plots.histogram!(sample2)