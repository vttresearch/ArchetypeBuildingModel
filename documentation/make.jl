using Documenter
using ArchetypeBuildingModel
using JSON


## Automatically parse together the input data description based on `archetype_definitions.json`.

filename = "archetype_definitions.md"
definitions = JSON.parsefile("archetype_definitions.json")
system_string = ["""
# Input data reference

This section contains automatically generated descriptions of the
*Spine Datastore* structure used for the archetype definitions and building stock input data,
based on the `archetype_definitions.json` file.
If you are not familiar with *Spine Datastores*, see the
[Spine Toolbox](https://github.com/Spine-project/Spine-Toolbox) documentation.

## Contents
```@contents
Pages = ["$(filename)"]
Depth = 3
```\n\n
"""]
# Add the object classes section.
push!(system_string, "## Object classes\n\n")
for (name, description) in definitions["object_classes"]
    push!(system_string, "### `$(name)`\n\n")
    push!(system_string, ">$(description)\n\n")
end
# Add the object parameters section.
push!(system_string, "## Object parameters\n\n")
for (object_class, name, default, value_list, description) in
    definitions["object_parameters"]
    push!(system_string, "### `$(name)`\n\n")
    push!(system_string, ">$(description)\n\n")
    push!(system_string, "Object class: `$(object_class)`\n\n")
    push!(system_string, "Default value: `$(default)`\n\n")
    push!(system_string, "Value list: `$(value_list)`\n\n")
end
# Add the object section.
push!(system_string, "## Objects\n\n")
for (object_class, name, description) in definitions["objects"]
    push!(system_string, "### `$(name)`\n\n")
    push!(system_string, ">$(description)\n\n")
    push!(system_string, "Object class: `$(object_class)`\n\n")
end
# Add the relationship_classes section
push!(system_string, "## Relationship classes\n\n")
for (name, object_classes, description, logo) in definitions["relationship_classes"]
    push!(system_string, "### `$(name)`\n\n")
    push!(system_string, ">$(description)\n\n")
    push!(system_string, "Object classes: $(join(object_classes, ", ", ", and "))\n\n")
end
# Add the relationship_parameters section
push!(system_string, "## Relationship parameters\n\n")
for (relationship_class, name, default, value_list, description) in
    definitions["relationship_parameters"]
    push!(system_string, "### `$(name)`\n\n")
    push!(system_string, ">$(description)\n\n")
    push!(system_string, "Relationship class: `$(relationship_class)`\n\n")
    push!(system_string, "Default value: `$(default)`\n\n")
    push!(system_string, "Value list: `$(value_list)`\n\n")
end

# Write the `system_string` into a markdown file.
open("documentation/src/$(filename)", "w") do file
    write(file, join(system_string))
end


## Make the documentation.

makedocs(
    sitename="ArchetypeBuildingModel",
    format=Documenter.HTML(size_threshold=300000),
    modules=[ArchetypeBuildingModel],
    pages=[
        "index.md",
        "defining_archetype_buildings.md",
        "model_workflow.md",
        "archetype_modelling.md",
        "archetypebuildingweather.md",
        "output_processing.md",
        "archetype_definitions.md",
        "library.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
