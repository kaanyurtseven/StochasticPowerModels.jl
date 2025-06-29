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



using ExaModels, CUDA
using MadNLPGPU
using MadNLPMumps


Memento.setlevel!(Memento.getlogger(ExaModels), "error")
Memento.setlevel!(Memento.getlogger(CUDA), "error")
Memento.setlevel!(Memento.getlogger(MadNLPGPU), "error")

# using HSL_jll


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
                                                               "max_iter"=>1000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-4,
                                                               # "hsllib" => HSL_jll.libhsl_path,
                                                               # "linear_solver" => "ma27",
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

#PCE degree input
deg  = 2
#RES input
pen_level = 0.3

#Case file and data reading
case = "case5_ACDC_SPM_95cc.m"

# case = "pglib_opf_case588_sdet_acdc_SPM_90cc.m"

file  = joinpath(BASE_DIR, "test/data/matpower", case)
data = _PM.parse_file(file)
_PMACDC.process_additional_data!(data)


s = Dict("RES Curtailment" => true,
         "Load Curtailment" => false,
         )


total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus


sdata = _SPM.build_stochastic_data_ACDC_RES(data, deg, p_size)
_SPM.dimension_management_PCE(sdata)
sdata["curtailment"] = s



CPU_solver = JuMP.optimizer_with_attributes(() -> ExaModels.MadNLPOptimizer(print_level=MadNLP.NOTICE,
                                                        max_iter = 3000,
                                                        tol = 1e-4,
                                                        bound_relax_factor = 1e-4,
                                                        acceptable_tol = 1e-4
                                                        ))

GPU_solver = JuMP.optimizer_with_attributes(() -> ExaModels.MadNLPOptimizer(CUDABackend(),
                                                        print_level=MadNLP.NOTICE,
                                                        max_iter =5000,
                                                        tol = 1e-4,    
                                                        bound_relax_factor = 1e-4,
                                                        acceptable_tol = 1e-4, 
                                                        ))                # ensure you're using MadNLPGPU backend)


println("\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_spm_ipopt = _SPM.solve_sopf_iv_acdc(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm_ipopt["primal_status"]), ")")

println("\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_spm = _SPM.solve_sopf_iv_acdc(sdata, _PM.IVRPowerModel, CPU_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm["primal_status"]), ")")

println("\n\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_spm_GPU = _SPM.solve_sopf_iv_acdc(sdata, _PM.IVRPowerModel, GPU_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm_GPU["primal_status"]), ")")



# Show results on the terminal
println("\n\n>>> GPU Results >>>")
println(result_spm_ipopt["primal_status"])
print("Objective: ")
print(result_spm_ipopt["objective"])
print("\nSolve Time: ")
print(result_spm_ipopt["solve_time"])

# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm["primal_status"])
print("Objective: ")
print(result_spm["objective"])
print("\nSolve Time: ")
print(result_spm["solve_time"])


# Show results on the terminal
println("\n\n>>> GPU Results >>>")
println(result_spm_GPU["primal_status"])
print("Objective: ")
print(result_spm_GPU["objective"])
print("\nSolve Time: ")
print(result_spm_GPU["solve_time"])





# bus = 4

# sample1 = _SPM.sample(result_spm_GPU, "gen", bus, "pg"; sample_size=100000); 
# sample2 = _SPM.sample(result_spm, "gen", bus, "pg"; sample_size=100000); 


# Plots.histogram(sample1, bin=50)

# Plots.histogram!(sample2, bin=50)





                                                

# # ---- Setup model builder
# model_type = _PM.IVRPowerModel
# build_function = _SPM.build_sopf_iv_acdc
# extensions = [_PMACDC.add_ref_dcgrid!, _SPM.add_ref_RES!]

# # ---- Build JuMP model manually
# jump_model = JuMP.Model()

# _PM.instantiate_model(
#     sdata,
#     model_type,
#     build_function;
#     ref_extensions = extensions,
#     jump_model = jump_model
# );


# exa_model = ExaModels.ExaModel(jump_model; backend = CUDABackend())


# # Solve using native MadNLP call (this runs fully on GPU)
# result_gpu = madnlp(
#     exa_model;
#     # warm_start_init_point = true,
#     print_level = MadNLP.INFO,
#     max_iter = 5000,
#     tol = 1e-4,
#     acceptable_tol = 1e-4,
#     bound_relax_factor = 1e-4
# )


# # optimizer(exa_model, optimizer)
