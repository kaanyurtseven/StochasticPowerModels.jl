################################################################################
#  Copyright 2021, Frederik Geth                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

""
function run_sopf_iv_reduced(data, model_constructor, optimizer; deg::Int=1, kwargs...)
    sdata = build_stochastic_data(data, deg)

    return _PMs.run_model(sdata, model_constructor, optimizer, build_sopf_iv_reduced; multinetwork=true, kwargs...)
end

""
function build_sopf_iv_reduced(pm::AbstractPowerModel)
    for (n, network) in _PMs.nws(pm) 
        variable_bus_voltage(pm, nw=n)

        variable_branch_current_reduced(pm, nw=n)

        variable_gen_power(pm, nw=n, bounded=false) # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false) # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)
    end

    for i in _PMs.ids(pm, :bus, nw=1)
        constraint_bus_voltage_squared_cc_limit(pm, i, nw=1)
    end

    for g in _PMs.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PMs.ids(pm, :branch, nw=1)
        constraint_branch_series_current_squared_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PMs.nws(pm)
        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
            constraint_gp_bus_voltage_squared(pm, i, nw=n)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            _PMs.constraint_voltage_drop(pm, b, nw=n)

            constraint_gp_branch_series_current_squared(pm, b, nw=n)
        end

        for g in _PMs.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PMs.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    objective_min_expected_generation_cost(pm)
end