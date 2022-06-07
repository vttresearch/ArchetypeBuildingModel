#=
    export.jl

Functions for exporting `ArchetypeBuildingModel.jl` outputs into different
formats, e.g. Spine Datastores and .csv files.
=#


"""
    write_to_url(url::String, input::ModelInput)

Write model `input` data into the database at the `url`.
"""
function write_to_url(url::String, input::ModelInput)
    data = [getfield(input, field) for field in fieldnames(typeof(input))]
    @info "Importing `$(typeof(input))` to datastore at `$(url)`..."
    @time import_data(
        url,
        data,
        "Import `$(typeof(input))` from `ArchetypeBuildingModel` output.",
    )
end
