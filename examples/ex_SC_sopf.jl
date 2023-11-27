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
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>3000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-4,
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

#gPC degree input
deg  = 2

#RES input
pen_level = 0.5

#Case file and data reading
case = "case5_ACDC_SPM_95cc_RES.m"
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data = _PM.parse_file(file)
_PMACDC.process_additional_data!(data)

curt_settings = Dict("RES Curtailment" => false,
                     "Load Curtailment" => false,
                     )

action_settings = Dict("Preventive Redispatch" => true,
                       "Preventive HVDC" => true,
                     )

objective_settings = Dict("Minimize Dispatch" => false,
                          "Minimize Re-dispatch" => false,
                          "Minimize Both" => true,
                     )


total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus

sdata = _SPM.build_stochastic_data_ACDC_RES(data, deg, p_size)

contingencies = Dict(parse(Int, c) => cont for (c, cont) in data["contingencies"])
number_of_cont = length(contingencies)

PCE_data = Dict(
      "T2" => sdata["T2"],
      "T3" => sdata["T3"],
      "T4" => sdata["T4"],
      "mop" => sdata["mop"],
   )
num_of_coeff = PCE_data["mop"].dim

delete!(sdata, "T2")
delete!(sdata, "T3")
delete!(sdata, "T4")
delete!(sdata, "mop")


sdata_cont = Dict()
sdata_cont = deepcopy(sdata)
sdata_cont["nw"] = Dict{String, Any}()
sdata_cont["name"] = Dict()

for i = 1:number_of_cont
   for j = 1:num_of_coeff

      sdata_cont["nw"]["$((i-1)*num_of_coeff+j)"] = deepcopy(sdata["nw"]["$j"])

   end
end

_FP.add_dimension!(sdata_cont, :PCE_coeff, num_of_coeff; metadata = PCE_data)
_FP.add_dimension!(sdata_cont, :contingencies, contingencies)
sdata_cont["curtailment"] = curt_settings
sdata_cont["action"] = action_settings
sdata_cont["objective"] = objective_settings

sdata_cont["gen_buses"] = Dict()
dummy = zeros(length(data["gen"]),1)
dummy = [floor(Int,x) for x in dummy]
for i in 1:length(data["gen"])
    dummy[i] = data["gen"]["$i"]["gen_bus"]
end
sdata_cont["gen_buses"] = dummy



for i = 1:number_of_cont*num_of_coeff

   cont_id = _FP.coord(sdata_cont, i, :contingencies)

   if cont_id >=2 #skip the base case

      for b in 1:length(sdata_cont["nw"]["$i"]["branch"])

         if contingencies[cont_id]["branch_id1"] == b

            sdata_cont["nw"]["$i"]["branch"]["$b"]["br_status"] = 0

         end

      end

      for b in 1:length(sdata_cont["nw"]["$i"]["branchdc"])
         
         if contingencies[cont_id]["branchdc_id1"] == b

            sdata_cont["nw"]["$i"]["branchdc"]["$b"]["status"] = 0

         end

      end

      for c in 1:length(sdata_cont["nw"]["$i"]["convdc"])
         
         if contingencies[cont_id]["convdc_id1"] == c

            sdata_cont["nw"]["$i"]["convdc"]["$c"]["status"] = 0

         end

      end

   end

end




println("\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_spm = _SPM.solve_sopf_iv_acdc_SC(sdata_cont, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm["primal_status"]), ")")



# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm["primal_status"])
print("Objective: ")
print(result_spm["objective"])
print("\nSolve Time: ")
print(result_spm["solve_time"])


result_base = deepcopy(result_spm)
# result_base["solution"]["nw"] = Dict()
for i = 1:num_of_coeff*number_of_cont
   
   if i > num_of_coeff
      delete!(result_base["solution"]["nw"], "$i")
   end
end



bus = 2

sample1 = _SPM.sample(result_base, "gen", 1, "pg"; sample_size=100000); 
sample2 = _SPM.sample(result_base, "gen", 2, "pg"; sample_size=100000); 


Plots.histogram(sample1)
# Plots.histogram!(sample2)

Plots.histogram(sample2)
Plots.histogram!(sample1)



# for x in sample1
#    if x>0 && x<0.15
#       print("\n$x")
#    end
# end

minimum(sample1)
minimum(sample2)



# file_name = "Generation Data 0.90.xlsx"
# fid    = XLSX.openxlsx(file_name, mode="w")
# header = ["Gen1";"Gen2"]
# data   = [[sample1];[sample2]]
# XLSX.writetable(file_name, data,header)