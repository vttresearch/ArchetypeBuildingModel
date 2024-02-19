#=
    export.jl

Functions for exporting `ArBuMo.jl` outputs into different
formats, e.g. Spine Datastores and .csv files.
=#


"""
    write_to_url(url::String, input::ModelInput; alternative::String="")

Write model `input` data into the database at the `url`.

The `alternative` the parameter values are saved to depends on the given keyword
"""
function write_to_url(url::String, input::ModelInput; alternative::String="")
    # Convert object and relationship classes to their Dict representations.
    data =
        SpineInterface._to_dict.([
            getfield(input, field) for field in fieldnames(typeof(input))
        ])

    # Process `alternative` if given.
    if !isempty(alternative)
        @info "Processing `alternative` `$(alternative)` for import..."
        @time for d in data
            push!.(get(d, :object_parameter_values, [[]]), alternative)
            push!.(get(d, :relationship_parameter_values, [[]]), alternative)
        end
    end

    # Merge data into a single Dict and import.
    data = merge(append!, data...)
    @info "Importing `$(typeof(input))` to datastore at `$(url)`..."
    @time import_data(
        url,
        data,
        "Import `$(typeof(input))` from `ArBuMo` output.",
    )
end
