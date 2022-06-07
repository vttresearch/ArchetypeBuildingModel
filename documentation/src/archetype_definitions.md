# Archetype definition and data format

This section contains automatically generated descriptions of the
*Spine Datastore* structure used for the archetype definitions and building stock input data,
based on the `archetype_definitions.json` file.
If you are not familiar with *Spine Datastores*, see the
[Spine Toolbox](https://github.com/Spine-project/Spine-Toolbox) documentation.

## Contents
```@contents
Pages = ["archetype_definitions.md"]
Depth = 3
```


## Object classes

### `building_archetype`

>Represents a portion of the building stock, aggregated based on the related elements. Results in a lumped-capacitance thermal model of the defined archetype building.

### `building_fabrics`

>Represents the building fabrics like walls, floors, ceilings, interior air, etc., as well as how they are represented as a network of lumped-capacitance nodes.

### `building_loads`

>Contains information about the internal heat loads and demands of buildings.

### `building_node`

>Represents individual lumped-capacitances for building fabrics, as well as for building technical systems.

### `building_period`

>Represents a period of time between a start year and an end year for the statistical data.

### `building_process`

>Represents energy conversion and transfer processes between lumped-capacitances in the building technical systems.

### `building_scope`

>Defines aggregation for the building stock over time, geographical location, building type, and heat source.

### `building_stock`

>Represents a particular building stock dataset.

### `building_systems`

>Represents building technical systems via processes and lumped-capacitance nodes.

### `building_type`

>Classification for different types of buildings.

### `building_weather`

>Represents weather for the archetype buildings.

### `direction`

>Direction of process commodity flows, either into or from nodes.

### `heat_source`

>Classification for primary energy sources for heating/cooling systems.

### `location_id`

>Geographical location identifier.

### `structure_type`

>Classification for different building envelope structures.

## Object parameters

### `building_frame_depth_m`

>Assumed average depth of the building. Default value based on an industry rule-of-thumb regarding natural light conditions (Rakenteellinen energiatehokkuus - Opas, 2015).

Object class: `building_archetype`

Default value: `10.4`

Value list: `nothing`

### `effective_thermal_capacity_of_interior_air_and_furniture_J_m2K`

>Assumed effective thermal capacity [J/m2K gross-floor area]  of interior air and furniture. Default value based on SFS-EN ISO 52016-1:2017 Table B.17.

Object class: `building_archetype`

Default value: `10000.0`

Value list: `nothing`

### `external_shading_coefficient`

>Assumed average external shading coefficient, effectively a factor for the amount of direct solar radiation reaching the windows. Default value based on SFS-EN ISO 5216-1:2017 Table B.18 due to lack of a better source.

Object class: `building_archetype`

Default value: `0.5`

Value list: `nothing`

### `external_wall_load_bearing_fraction`

>Assumed average fraction of external walls that are load-bearing. Default value guesstimated.

Object class: `building_archetype`

Default value: `0.5`

Value list: `nothing`

### `indoor_air_cooling_set_point_override_K`

>Option for overriding nodal maximum indoor air temperature on an archetype level in K.

Object class: `building_archetype`

Default value: `nothing`

Value list: `nothing`

### `indoor_air_heating_set_point_override_K`

>Option for overriding nodal minimum indoor air temperature on an archetype level in K.

Object class: `building_archetype`

Default value: `nothing`

Value list: `nothing`

### `internal_heat_gain_convective_fraction`

>Assumed fraction of convective over radiative internal heat gains. Convective heat gains are applied directly to the interior air and furniture, while radiative heat gains are applied to surrounding structures. Default value based on SFS-EN ISO 52016-1:2017 Table B.11.

Object class: `building_archetype`

Default value: `0.4`

Value list: `nothing`

### `number_of_storeys`

>Assumed average number of storeys.

Object class: `building_archetype`

Default value: `2.0`

Value list: `nothing`

### `partition_wall_length_ratio_to_external_walls_m_m`

>Assumed average length ratio of internal partition walls compared to the external walls [m/m].

Object class: `building_archetype`

Default value: `0.5`

Value list: `nothing`

### `partition_wall_load_bearing_fraction`

>Assumed average fraction of internal partition walls that are load-bearing.

Object class: `building_archetype`

Default value: `0.5`

Value list: `nothing`

### `room_height_m`

>Assumed average room height in metres.

Object class: `building_archetype`

Default value: `2.5`

Value list: `nothing`

### `solar_heat_gain_convective_fraction`

>Assumed fraction of convective over radiative solar heat gains. Convective heat gains are applied directly to the interior air and furniture, while radiative heat gains are applied to surrounding structures. Default value based on IDA ESBO calibrations.

Object class: `building_archetype`

Default value: `0.6`

Value list: `nothing`

### `volumetric_heat_capacity_of_interior_air_J_m3K`

>Isobaric volumetric heat capacity of air at approximately room temperature ~20C. Default value based on Wikipedia.

Object class: `building_archetype`

Default value: `1210.0`

Value list: `nothing`

### `weather_year`

>The year for which weather data will be processed, unless a `building_weather` object is connected to the archetype manually.

Object class: `building_archetype`

Default value: `nothing`

Value list: `nothing`

### `window_area_distribution_towards_cardinal_directions`

>Assumed distribution of window area towards the cardinal directions. By default, window area is assumed to be distributed equally towards all cardinal directions.

Object class: `building_archetype`

Default value: `Dict{String, Any}("data" => Any[Any["north", 0.25], Any["east", 0.25], Any["south", 0.25], Any["west", 0.25]], "type" => "map", "index_name" => "cardinal_direction", "index_type" => "str")`

Value list: `nothing`

### `window_area_to_external_wall_ratio_m2_m2`

>Assumed average ratio between external wall surface area and the surface area of windows [m2/m2].

Object class: `building_archetype`

Default value: `0.2`

Value list: `nothing`

### `window_non_perpendicularity_correction_factor`

>A correction factor for estimating the total solar energy transmittance of windows from the normal solar energy transmittance. Default value based on SFS-EN ISO 52016-1:2017 Table B.43.

Object class: `building_archetype`

Default value: `0.9`

Value list: `nothing`

### `domestic_hot_water_demand_base_W`

>Base component of total domestic hot water demand [W].

Object class: `building_loads`

Default value: `0.0`

Value list: `nothing`

### `domestic_hot_water_demand_gfa_scaling_W_m2`

>Gross-floor-area-scaling component of total domestic hot water demand [W/m2 gross-floor area].

Object class: `building_loads`

Default value: `0.0`

Value list: `nothing`

### `internal_heat_loads_base_W`

>Base component of total internal heat loads [W].

Object class: `building_loads`

Default value: `0.0`

Value list: `nothing`

### `internal_heat_loads_gfa_scaling_W_m2`

>Gross-floor-area-scaling component of total internal heat loads [W/m2 gross-floor area].

Object class: `building_loads`

Default value: `0.0`

Value list: `nothing`

### `domestic_hot_water_demand_weight`

>Weight for how much of the domestic hot water demand is allocated for this node.

Object class: `building_node`

Default value: `0.0`

Value list: `nothing`

### `effective_thermal_mass_base_J_K`

>Define the base component of total effective thermal mass of the node by hand [J/K].

Object class: `building_node`

Default value: `0.0`

Value list: `nothing`

### `effective_thermal_mass_gfa_scaling_J_m2K`

>Define the gross-floor-area-scaling component of total effective thermal mass of the node by hand [J/m2K]

Object class: `building_node`

Default value: `0.0`

Value list: `nothing`

### `interior_air_and_furniture_weight`

>Weight for how much of the interior air and furniture is included in this node.

Object class: `building_node`

Default value: `0.0`

Value list: `nothing`

### `maximum_permitted_temperature_K`

>Maximum allowed temperature for the node [K].

Object class: `building_node`

Default value: `400.0`

Value list: `nothing`

### `minimum_permitted_temperature_K`

>Minimum allowed temperature for the node [K].

Object class: `building_node`

Default value: `200.0`

Value list: `nothing`

### `self_discharge_rate_base_W_K`

>Define the base component of total self-discharge rate of the node [W/K], where energy is lost outside the model scope.

Object class: `building_node`

Default value: `0.0`

Value list: `nothing`

### `self_discharge_rate_gfa_scaling_W_m2K`

>Define the gross-floor-area-scaling component of total self-discharge rate of the node [W/m2K], where energy is lost outside the model scope.

Object class: `building_node`

Default value: `0.0`

Value list: `nothing`

### `period_end`

>nothing

Object class: `building_period`

Default value: `nothing`

Value list: `nothing`

### `period_start`

>nothing

Object class: `building_period`

Default value: `nothing`

Value list: `nothing`

### `coefficient_of_performance_base`

>The base coefficient of performance for the `building_process`. For temperature-dependent COPs, this parameter needs to be the COP at known reference temperatures, divided by the Carnot COP at those same temperatures.

Object class: `building_process`

Default value: `1.0`

Value list: `nothing`

### `coefficient_of_performance_minimum_temperature_delta`

>The minimum assumed temperature raise/decrease in the heat pump process, limiting the maximum Carnot COP of the process. This parameter matters only if source and sink temperatures are defined.

Object class: `building_process`

Default value: `5.0`

Value list: `nothing`

### `coefficient_of_performance_mode`

>The mode of the heat pump process, affecting how the source and sink temperatures are interpreted. Set to `cooling` for heat pumpts used for cooling.

Object class: `building_process`

Default value: `heating`

Value list: `COP_modes`

### `coefficient_of_performance_sink_temperature_K`

>The sink temperature for the heat pump process, defined as a `Map`. Use `ambient` and `ground` for weather dependent `temperature_K`, while the rest of the properties can be used to tweak the assumed heating curve.

Object class: `building_process`

Default value: `Dict{String, Any}("data" => Any[Any["temperature_K", nothing], Any["heating_curve_control_temperature_min_K", nothing], Any["heating_curve_control_temperature_max_K", nothing], Any["heating_curve_output_temperature_min_K", nothing], Any["heating_curve_output_temperature_max_K", nothing]], "type" => "map", "index_name" => "property", "index_type" => "str")`

Value list: `nothing`

### `coefficient_of_performance_source_temperature_K`

>The source temperature for the heat pump process, defined as a `Map`. Use `ambient` and `ground` for weather dependent `temperature_K`, while the rest of the properties can be used to tweak the assumed heating curve.

Object class: `building_process`

Default value: `Dict{String, Any}("data" => Any[Any["temperature_K", nothing], Any["heating_curve_control_temperature_min_K", nothing], Any["heating_curve_control_temperature_max_K", nothing], Any["heating_curve_output_temperature_min_K", nothing], Any["heating_curve_output_temperature_max_K", nothing]], "type" => "map", "index_name" => "property", "index_type" => "str")`

Value list: `nothing`

### `scope_period_end_year`

>Last year [y] of the construction time period to be included in the statistical scope.

Object class: `building_scope`

Default value: `nothing`

Value list: `nothing`

### `scope_period_start_year`

>First year [y] of the construction time period to be included in the statistical scope.

Object class: `building_scope`

Default value: `nothing`

Value list: `nothing`

### `building_stock_year`

>The year this `building_stock` is supposed to represent, like a snapshot of the building stock during this year.

Object class: `building_stock`

Default value: `nothing`

Value list: `nothing`

### `raster_weight_path`

>An optional filepath to a geographical raster data file containing weighting information for the weather data, e.g. population density or the like.

Object class: `building_stock`

Default value: `nothing`

Value list: `nothing`

### `shapefile_path`

>The filepath to a shapefile containing the geographical information about the building stock. Required for weather data processing.

Object class: `building_stock`

Default value: `nothing`

Value list: `nothing`

### `ambient_temperature_K`

>Ambient air temperature [K].

Object class: `building_weather`

Default value: `nothing`

Value list: `nothing`

### `diffuse_solar_irradiation_W_m2`

>Diffuse solar irradiation [W/m2].

Object class: `building_weather`

Default value: `nothing`

Value list: `nothing`

### `direct_solar_irradiation_W_m2`

>Direct solar irradiation on walls facing towards the cardinal directions [W/m2].

Object class: `building_weather`

Default value: `Dict{String, Any}("data" => Any[Any["north", 0.0], Any["east", 0.0], Any["south", 0.0], Any["west", 0.0]], "type" => "map", "index_name" => "cardinal_direction", "index_type" => "str")`

Value list: `nothing`

### `location_name`

>Name of the location corresponding to the identifier, e.g. the name of the municipality, region, country, etc.

Object class: `location_id`

Default value: `nothing`

Value list: `nothing`

### `exterior_resistance_m2K_W`

>Exterior surface thermal resistance of a structure [m2K/W].

Object class: `structure_type`

Default value: `nothing`

Value list: `nothing`

### `interior_resistance_m2K_W`

>Interior surface thermal resistance of a structure [m2K/W].

Object class: `structure_type`

Default value: `nothing`

Value list: `nothing`

### `is_internal`

>A boolean flag for whether a structure type is internal, meaning inside the building envelope and not directly in contact with ambient conditions.

Object class: `structure_type`

Default value: `nothing`

Value list: `nothing`

### `is_load_bearing`

>A boolean flag for whether a structure type is load-bearing, meaning it is designed to bear the weight of structures on top of it.

Object class: `structure_type`

Default value: `true`

Value list: `nothing`

### `linear_thermal_bridge_W_mK`

>Linear thermal bridges for the structure type [W/mK] caused by seams between structures.

Object class: `structure_type`

Default value: `nothing`

Value list: `nothing`

### `structure_type_notes`

>Generic notes about the different structure types.

Object class: `structure_type`

Default value: `nothing`

Value list: `nothing`

## Objects

### `@system_link_node_1`

>Reserved special node intended to be used to link the archetype building systems to the desired connection point in the overarching energy system model, defined using the `building_archetype__system_link_node` relationship class.

Object class: `building_node`

### `@system_link_node_2`

>Reserved special node intended to be used to link the archetype building systems to the desired connection point in the overarching energy system model, defined using the `building_archetype__system_link_node` relationship class.

Object class: `building_node`

### `@system_link_node_3`

>Reserved special node intended to be used to link the archetype building systems to the desired connection point in the overarching energy system model, defined using the `building_archetype__system_link_node` relationship class.

Object class: `building_node`

### `from_node`

>Direction of energy flow from a node into a process.

Object class: `direction`

### `to_node`

>Direction of energy flow from a process and into a node.

Object class: `direction`

## Relationship classes

### `building_archetype__building_fabrics`

>Defines which building fabric representation is used for the chosen archetype building.

Object classes: building_archetype, and building_fabrics

### `building_archetype__building_loads`

>Defines the internal heat loads for the archetype building.

Object classes: building_archetype, and building_loads

### `building_archetype__building_scope`

>Defines the statistical scope of the archetype building.

Object classes: building_archetype, and building_scope

### `building_archetype__building_systems`

>Defines the heating/cooling systems used by the archetype building.

Object classes: building_archetype, and building_systems

### `building_archetype__building_weather`

>Defines the weather for the archetype building.

Object classes: building_archetype, and building_weather

### `building_archetype__system_link_node`

>Defines how the archetype building is connected to overarching models using the reserved `@system_link_node`s.

Object classes: building_archetype, and building_node

### `building_fabrics__building_node`

>Defines which nodes make up the representation of the building fabrics.

Object classes: building_fabrics, and building_node

### `building_node__building_node`

>Allows configuring heat transfer coefficients between building nodes by hand. Useful for e.g. heating/cooling system definitions.

Object classes: building_node, and building_node

### `building_node__structure_type`

>Defines which structures are included in which nodes, as well as their weights within those nodes.

Object classes: building_node, and structure_type

### `building_process__direction__building_node`

>Defines how processes interact with building nodes, as well as sets process-flow specific technical parameters.

Object classes: building_process, direction, and building_node

### `building_scope__building_stock`

>Defines which `building_stock` datasets are used in the statistical scope, as well as their weights.

Object classes: building_scope, and building_stock

### `building_scope__building_type`

>Defines the building types included in the statistical scope, as well as their weights.

Object classes: building_scope, and building_type

### `building_scope__heat_source`

>Defines the heat sources included in the statistical scope, as well as their weights.

Object classes: building_scope, and heat_source

### `building_scope__location_id`

>Defines the geographical locations included in a statistical scope, as well as their weights.

Object classes: building_scope, and location_id

### `building_stock_statistics`

>Contains the statistical data about the composition of the building stock.

Object classes: building_stock, building_type, building_period, location_id, and heat_source

### `building_systems__building_node`

>Defines which lumped-capacitance nodes are included in the building systems.

Object classes: building_systems, and building_node

### `building_systems__building_process`

>Defines which processes are included in the building systems.

Object classes: building_systems, and building_process

### `structure_statistics`

>Contains statistics about the properties of different structure types.

Object classes: building_type, building_period, location_id, and structure_type

### `ventilation_and_fenestration_statistics`

>Contains statistics about thje properties of fenestration and ventilation.

Object classes: building_type, building_period, and location_id

## Relationship parameters

### `grid_name`

>Name of the Backbone `grid` this `@system_link_node` is included in.

Relationship class: `building_archetype__system_link_node`

Default value: `nothing`

Value list: `nothing`

### `node_name`

>Name of the node this `@system_link_node` is representing in the overarching large-scale energy system model.

Relationship class: `building_archetype__system_link_node`

Default value: `nothing`

Value list: `nothing`

### `heat_transfer_coefficient_base_W_K`

>Set the base heat transfer coefficient between two nodes [W/K] by hand.

Relationship class: `building_node__building_node`

Default value: `0.0`

Value list: `nothing`

### `heat_transfer_coefficient_gfa_scaling_W_m2K`

>Set the gross-floor-area-scaling heat transfer coefficient between two nodes [W/m2K] by hand.

Relationship class: `building_node__building_node`

Default value: `0.0`

Value list: `nothing`

### `structure_type_weight`

>Weight or share of the structure type to be included in this node.

Relationship class: `building_node__structure_type`

Default value: `1.0`

Value list: `nothing`

### `maximum_power_base_W`

>Maximum base power of the process flow [W].

Relationship class: `building_process__direction__building_node`

Default value: `0.0`

Value list: `nothing`

### `maximum_power_gfa_scaling_W_m2`

>Gross-floor-area-scaling maximum power of a process flow [W/m2].

Relationship class: `building_process__direction__building_node`

Default value: `0.0`

Value list: `nothing`

### `building_stock_weight`

>Weight for sampling the `building_stock` within the `building_scope`.

Relationship class: `building_scope__building_stock`

Default value: `1.0`

Value list: `nothing`

### `building_type_weight`

>Weight for sampling the `building_type` within the `building_scope`.

Relationship class: `building_scope__building_type`

Default value: `1.0`

Value list: `nothing`

### `heat_source_weight`

>Weight for sampling the `heat_source` within the `building_scope`.

Relationship class: `building_scope__heat_source`

Default value: `1.0`

Value list: `nothing`

### `location_id_weight`

>Weight for sampling the `location_id` within the `building_scope`.

Relationship class: `building_scope__location_id`

Default value: `1.0`

Value list: `nothing`

### `average_gross_floor_area_m2_per_building`

>Statistical average gross-floor-area per building [m2].

Relationship class: `building_stock_statistics`

Default value: `nothing`

Value list: `nothing`

### `number_of_buildings`

>Statistical number of buildings.

Relationship class: `building_stock_statistics`

Default value: `nothing`

Value list: `nothing`

### `design_U_value_W_m2K`

>Mean original design U-value [W/m2K] of the structures corresponding to the statistics.

Relationship class: `structure_statistics`

Default value: `nothing`

Value list: `nothing`

### `effective_thermal_mass_J_m2K`

>Mean calculated effective thermal mass [J/m2K] of the structures corresponding to the statistics, per area of the structure.

Relationship class: `structure_statistics`

Default value: `nothing`

Value list: `nothing`

### `external_U_value_to_ambient_air_W_m2K`

>Mean calculated U-value [W/m2K] from the structure into the ambient air.

Relationship class: `structure_statistics`

Default value: `0.0`

Value list: `nothing`

### `external_U_value_to_ground_W_m2K`

>Mean calculated effective U-value [W/m2K] from the structure into the ground, according to KissockK2013.

Relationship class: `structure_statistics`

Default value: `0.0`

Value list: `nothing`

### `internal_U_value_to_structure_W_m2K`

>Mean calculated U-value [W/m2K] from the structure into the interior air.

Relationship class: `structure_statistics`

Default value: `nothing`

Value list: `nothing`

### `linear_thermal_bridges_W_mK`

>Mean linear thermal bridges [W/mK] of the seams between structures.

Relationship class: `structure_statistics`

Default value: `nothing`

Value list: `nothing`

### `total_U_value_W_m2K`

>Mean total effective U-value [W/m2K] of the structure, from the interior air into the ambient air/ground.

Relationship class: `structure_statistics`

Default value: `nothing`

Value list: `nothing`

### `HRU_efficiency`

>Mean heat-recovery efficiency of ventilation heat-recovery units corresponding to the statistics.

Relationship class: `ventilation_and_fenestration_statistics`

Default value: `nothing`

Value list: `nothing`

### `infiltration_rate_1_h`

>Mean infiltration air change rate [1/h] corresponding to the statistics.

Relationship class: `ventilation_and_fenestration_statistics`

Default value: `nothing`

Value list: `nothing`

### `total_normal_solar_energy_transmittance`

>Mean total normal solar energy transmittance of windows corresponding to the statistics, already including the effect of the frame-area fraction.

Relationship class: `ventilation_and_fenestration_statistics`

Default value: `nothing`

Value list: `nothing`

### `ventilation_rate_1_h`

>Mean ventilation air change rate [1/h] corresponding to the statistics.

Relationship class: `ventilation_and_fenestration_statistics`

Default value: `nothing`

Value list: `nothing`

### `window_U_value_W_m2K`

>Mean window U-value [W/m2K] corresponding to the statistics.

Relationship class: `ventilation_and_fenestration_statistics`

Default value: `nothing`

Value list: `nothing`

