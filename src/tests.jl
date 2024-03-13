#=
    tests.jl

This file contains functions for ensuring the archetype building definitions,
as well as the necessary input data, are valid.
=#

"""
    run_parameter_tests(mod::Module = @__MODULE__; limit::Real = Inf)

Run tests for archetype building model definition `Parameters` for module `mod`.
The `limit` keyword can be used to limit the number of tests run.
"""
function run_parameter_tests(mod::Module=@__MODULE__; limit::Real=Inf)
    params_for_testing = [ # Object parameters
        mod.average_apparent_sky_temperature_difference_K =>
            (type=Real, min=0, max=20),
        mod.average_structural_solar_absorption_coefficient =>
            (type=Real, min=0, max=1),
        mod.building_frame_depth_m => (type=Real, min=0),
        mod.effective_thermal_capacity_of_interior_air_and_furniture_J_m2K =>
            (type=Real, min=0),
        mod.energy_efficiency_override_multiplier => (type=Real, min=0),
        mod.external_radiative_surface_heat_transfer_coefficient_W_m2K =>
            (type=Real, min=0, max=10),
        mod.external_shading_coefficient => (type=Real, min=0, max=1),
        mod.internal_heat_gain_convective_fraction => (type=Real, min=0, max=1),
        mod.number_of_storeys => (type=Real, min=0),
        mod.partition_wall_length_ratio_to_external_walls_m_m => (type=Real, min=0),
        mod.partition_wall_load_bearing_fraction => (type=Real, min=0, max=1),
        mod.room_height_m => (type=Real, min=0, max=10),
        mod.solar_heat_gain_convective_fraction => (type=Real, min=0, max=1),
        mod.volumetric_heat_capacity_of_interior_air_J_m3K =>
            (type=Real, min=0, max=10000),
        mod.weather_end => (type=Symbol,),
        mod.weather_start => (type=Symbol,),
        mod.window_area_distribution_towards_cardinal_directions => (type=Map,),
        mod.window_area_to_external_wall_ratio_m2_m2 => (type=Real, min=0, max=1),
        mod.window_area_thermal_bridge_surcharge_W_m2K => (type=Real,),
        mod.window_non_perpendicularity_correction_factor =>
            (type=Real, min=0, max=1),
        mod.domestic_hot_water_demand_base_W => (type=SpineDataType,),
        mod.domestic_hot_water_demand_gfa_scaling_W_m2 => (type=SpineDataType,),
        mod.internal_heat_loads_base_W => (type=SpineDataType,),
        mod.internal_heat_loads_gfa_scaling_W_m2 => (type=SpineDataType,),
        mod.domestic_hot_water_demand_weight => (type=Real, min=0, max=1),
        mod.effective_thermal_mass_base_J_K => (type=Real,),
        mod.effective_thermal_mass_gfa_scaling_J_m2K => (type=Real,),
        mod.is_interior_node => (type=Bool,),
        mod.maximum_permitted_temperature_K => (type=Real, min=0),
        mod.minimum_permitted_temperature_K => (type=Real, min=0),
        mod.self_discharge_rate_base_W_K => (type=Real, min=0),
        mod.self_discharge_rate_gfa_scaling_W_m2K => (type=Real, min=0),
        mod.period_end => (type=Real, min=0, max=2050),
        mod.period_start => (type=Real, min=0, max=2050),
        mod.coefficient_of_performance_base => (type=SpineDataType, min=0),
        mod.coefficient_of_performance_minimum_temperature_delta =>
            (type=Real, min=1, max=100),
        mod.coefficient_of_performance_mode => (type=Symbol,),
        mod.coefficient_of_performance_sink_temperature_K => (type=Map,),
        mod.coefficient_of_performance_source_temperature_K => (type=Map,),
        mod.scope_period_end_year => (type=Real, min=0, max=2100),
        mod.scope_period_start_year => (type=Real, min=0, max=2100),
        mod.building_stock_year => (type=Real, min=1900, max=2050),
        mod.shapefile_path => (type=Symbol,),
        mod.ambient_temperature_K => (type=SpineDataType,),
        mod.diffuse_solar_irradiation_W_m2 => (type=SpineDataType,),
        mod.direct_solar_irradiation_W_m2 => (type=Map,),
        mod.location_name => (type=Symbol,),
        mod.exterior_resistance_m2K_W => (type=Real, min=0),
        mod.interior_resistance_m2K_W => (type=Real, min=0),
        mod.is_internal => (type=Bool,),
        mod.is_load_bearing => (type=Bool,),
        mod.linear_thermal_bridge_W_mK => (type=Real, min=0),
        mod.structure_type_notes => (type=Symbol,),
        mod.grid_name => (type=Symbol,), # Relationship parameters from here on out.
        mod.node_name => (type=Symbol,),
        mod.heat_transfer_coefficient_base_W_K => (type=Real,),
        mod.heat_transfer_coefficient_gfa_scaling_W_m2K => (type=Real,),
        mod.structure_type_weight => (type=Real, min=0, max=1),
        mod.maximum_power_base_W => (type=Real, min=0),
        mod.maximum_power_gfa_scaling_W_m2 => (type=Real, min=0),
        mod.building_stock_weight => (type=Real, min=0), # Building scope weights allow upscaling.
        mod.building_type_weight => (type=Real, min=0), # Building scope weights allow upscaling.
        mod.heat_source_weight => (type=Real, min=0), # Building scope weights allow upscaling.
        mod.location_id_weight => (type=Real, min=0), # Building scope weights allow upscaling.
        mod.average_gross_floor_area_m2_per_building => (type=Real, min=0),
        mod.number_of_buildings => (type=Real, min=0),
        mod.design_U_value_W_m2K => (type=Real, min=0),
        mod.effective_thermal_mass_J_m2K => (type=Real, min=0),
        mod.external_U_value_to_ambient_air_W_m2K => (type=Real, min=0),
        mod.external_U_value_to_ground_W_m2K => (type=Real, min=0),
        mod.internal_U_value_to_structure_W_m2K => (type=Real, min=0),
        mod.linear_thermal_bridges_W_mK => (type=Real, min=0),
        mod.total_U_value_W_m2K => (type=Real, min=0),
        mod.HRU_efficiency => (type=Real, min=0, max=1),
        mod.infiltration_rate_1_h => (type=Real, min=0),
        mod.total_normal_solar_energy_transmittance => (type=Real, min=0, max=1),
        mod.ventilation_rate_1_h => (type=Real, min=0),
        mod.window_U_value_W_m2K => (type=Real, min=0),
    ]
    @testset "Testing `Parameters`." begin
        for (param, tup) in params_for_testing
            if !isnothing(tup)
                test_parameter(
                    param,
                    tup.type,
                    mod;
                    value_min=get(tup, :min, -Inf),
                    value_max=get(tup, :max, Inf),
                    limit=limit
                )
            end
        end
    end
end


"""
    run_object_class_tests(mod::Module = @__MODULE__; limit::Real = Inf)

Run tests for archetype building model definition `ObjectClasses` for module `mod`.
The `limit` keyword can be used to limit the number of tests run.
"""
function run_object_class_tests(mod::Module=@__MODULE__; limit::Real=Inf)
    obj_classes_for_testing = [
        mod.building_archetype =>
            (rel=mod.building_archetype__building_fabrics, min=1, max=1),
        mod.building_archetype =>
            (rel=mod.building_archetype__building_loads, min=1, max=1),
        mod.building_archetype =>
            (rel=mod.building_archetype__building_scope, min=1, max=1),
        mod.building_archetype =>
            (rel=mod.building_archetype__building_systems, min=1, max=1),
        mod.building_archetype =>
            (rel=mod.building_archetype__building_weather, min=0, max=1),
        mod.building_archetype =>
            (rel=mod.building_archetype__system_link_node, min=1),
        mod.building_fabrics => (rel=mod.building_fabrics__building_node, min=1),
        mod.building_process =>
            (rel=mod.building_process__direction__building_node, min=1),
        mod.building_scope => (rel=mod.building_scope__building_stock, min=1),
        mod.building_scope => (rel=mod.building_scope__building_type, min=1),
        mod.building_scope => (rel=mod.building_scope__heat_source, min=1),
        mod.building_scope => (rel=mod.building_scope__location_id, min=1),
        mod.building_systems => (rel=mod.building_systems__building_process, min=1),
    ]
    @testset "Testing `ObjectClasses`." begin
        for (oc, tup) in obj_classes_for_testing
            test_object_class(
                oc,
                tup.rel,
                mod;
                count_min=get(tup, :min, 0),
                count_max=get(tup, :max, Inf),
                limit=limit
            )
        end
    end
end


"""
    run_structure_type_tests(mod::Module = @__MODULE__)

Ensure that `structure_type` contains the correct objects.
"""
function run_structure_type_tests(mod::Module=@__MODULE__)
    valid_structure_type_names =
        Symbol.([
            "base_floor",
            "exterior_wall",
            "light_exterior_wall",
            "light_partition_wall",
            "partition_wall",
            "roof",
            "separating_floor",
        ])
    existing_structure_type_names = getfield.(mod.structure_type(), :name)
    if isempty(existing_structure_type_names)
        @warn "Testing `structure_type` skipped, as the object class is empty."
    else
        @testset "Testing that `structure_type` contains exactly the right entries." begin
            for valid_name in valid_structure_type_names
                @test _check(
                    valid_name in existing_structure_type_names,
                    "`$(valid_name)` isn't included in `structure_type`!",
                )
            end
            for existing_name in existing_structure_type_names
                @test _check(
                    existing_name in valid_structure_type_names,
                    "`$(existing_name)` is not a valid `structure_type`!",
                )
            end
        end
    end
end
