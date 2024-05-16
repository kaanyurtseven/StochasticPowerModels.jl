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
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>10000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-6,
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "atol" => 1e-6, "log_levels" => [])

#gPC degree input
deg  = 2

#RES input
# pen_level = 0.78
pen_level = 0.70

#Case file and data reading
# case = "case5_ACDC_SPM_95cc.m"
# case_det = "case5_ACDC_SPM_95cc.m"
case = "case67acdc_SPM_90cc_det.m"
case_det = "case67acdc_SPM_90cc_det.m"


file  = joinpath(BASE_DIR, "test/data/matpower", case)
file_det = joinpath(BASE_DIR, "test/data/matpower", case_det)


data = _PM.parse_file(file)
data_det_AC = _PM.parse_file(file_det) 
data_ACDC = _PM.parse_file(file)
data_OPF = _PM.parse_file(file)

_PMACDC.process_additional_data!(data)
_PMACDC.process_additional_data!(data_ACDC)
_PMACDC.process_additional_data!(data_OPF)

s = Dict("RES Curtailment" => true,
         "Load Curtailment" => false,
         )




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



total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus

for i = 1:length(data_det_AC["load"])
   data_det_AC["load"]["$i"]["pd"] = data["load"]["$i"]["pd"] #- p_size
end
_PMACDC.process_additional_data!(data_det_AC)

println("\n\n", "Solving Deterministic ACDC OTS Problem >>>>>", "\n")
    result_ACDC_ots = _PMACDC.run_acdcots(data_det_AC, _PM.ACPPowerModel, minlp_solver)

println("\n\n", "Solving Deterministic ACDC OPF Problem >>>>>", "\n")
    result_ACDC_opf = _PMACDC.run_acdcopf(data_det_AC, _PM.ACPPowerModel, minlp_solver)

println("\n\n>>> OTS Results >>>")
println(result_ACDC_ots["primal_status"])
print("Objective: ")
print(result_ACDC_ots["objective"])
print("\nSolve Time: ")
print(result_ACDC_ots["solve_time"])

println("\n\n>>> OPF Results >>>")
println(result_ACDC_opf["primal_status"])
print("Objective: ")
print(result_ACDC_opf["objective"])
print("\nSolve Time: ")
print(result_ACDC_opf["solve_time"])


for (b, branch) in result_ACDC_ots["solution"]["branch"]

   if branch["br_status"] < 0.5
       data_ACDC["branch"][b]["br_status_initial"] = 0
    #    data_OPF["branch"][b]["br_status"] = 0
   else
       data_ACDC["branch"][b]["br_status_initial"] = 1
   end
end

for (b, branch) in result_ACDC_ots["solution"]["branchdc"]

   if branch["br_status_dc"] < 0.5
       data_ACDC["branchdc"][b]["br_status_initial"] = 1
    #    data_OPF["branchdc"][b]["status"] = 1
   else
       data_ACDC["branchdc"][b]["br_status_initial"] = 1
   end

end


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

println("\n\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_SOTS = _SPM.solve_sots_iv_acdc(sdata_OTS, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_SOTS["primal_status"]), ")")

for (b, branch) in result_SOTS["solution"]["nw"]["1"]["branch"]

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

println("\n\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_SOPF = _SPM.solve_sopf_iv_acdc(sdata_OPF, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_SOPF["primal_status"]), ")")


# println("\n   (MINLP) Solution progress: Solving...")
# result_SOTS_bin = _SPM.solve_sots_iv_acdc_binary(sdata_OTS, _PM.IVRPowerModel, minlp_solver, deg=deg, p_size=p_size);
# println("   (MINLP) Solution progress: Solved! (", string(result_SOTS_bin["primal_status"]), ")")

if sdata_OTS["curtailment"]["RES Curtailment"] == true
    for i = 1:num_of_coeff
    result_SOPF["solution"]["nw"]["$i"]["RES"]["4"] = Dict()
    result_SOPF["solution"]["nw"]["$i"]["RES"]["4"]["p_RES_curt"] = sum([result_SOPF["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:3])


    result_SOTS["solution"]["nw"]["$i"]["RES"]["4"] = Dict()
    result_SOTS["solution"]["nw"]["$i"]["RES"]["4"]["p_RES_curt"] = sum([result_SOTS["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:3])

    # result_SOPF["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
    # result_SOPF["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_SOPF["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])


    # result_SOTS["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
    # result_SOTS["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_SOTS["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])

    end
end



println("\n\n>>> SPMACDC SOPF Results >>>")
println(result_SOPF["primal_status"])
print("Objective: ")
print(result_SOPF["objective"])
print("\nSolve Time: ")
print(result_SOPF["solve_time"])

println("\n\n>>> SPMACDC SOTS Results >>>")
println(result_SOTS["primal_status"])
print("Objective: ")
print(result_SOTS["objective"])
print("\nSolve Time: ")
print(result_SOTS["solve_time"])



# println("\n\n>>> SPMACDC SOTS (Binary) Results >>>")
# println(result_SOTS_bin["primal_status"])
# print("Objective: ")
# print(result_SOTS_bin["objective"])
# print("\nSolve Time: ")
# print(result_SOTS_bin["solve_time"])

# println("\n\n>>> SPMACDC SOTS (MINLP) Results >>>")
# println(result_SOTS_bin["primal_status"])
# print("Objective: ")
# print(result_SOTS_bin["objective"])
# print("\nSolve Time: ")
# print(result_SOTS_bin["solve_time"])

# sum([result_SOTS["solution"]["nw"]["1"]["branch"]["$i"]["br_status"] for i=1:102])

bus = 93

sample1 = _SPM.sample(result_SOTS, "branch", bus, "cr_fr"; sample_size=10000); 
sample2 = _SPM.sample(result_SOPF, "branch", bus, "cr_fr"; sample_size=10000); 

# bus = 11

# sample1 = _SPM.sample(result_SOTS, "branchdc", bus, "pf"; sample_size=10000); 
# sample2 = _SPM.sample(result_SOPF, "branchdc", bus, "pf"; sample_size=10000); 

# bus = 18

# sample1 = _SPM.sample(result_SOTS, "gen", bus, "pg"; sample_size=10000); 
# sample2 = _SPM.sample(result_SOPF, "gen", bus, "pg"; sample_size=10000); 

# bus = 67

# sample1 = _SPM.sample(result_SOTS, "bus", bus, "vm"; sample_size=10000); 
# sample2 = _SPM.sample(result_SOPF, "bus", bus, "vm"; sample_size=10000); 


Plots.histogram(sample1, bin=50)

Plots.histogram!(sample2, bin=50)