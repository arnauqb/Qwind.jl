using DrWatson
#@quickactivate "/cosma7/data/dp004/dc-quer1/Qwind.jl"
quickactivate("/cosma7/data/dp004/dc-quer1/Qwind.jl")
using Qwind
using ArgParse

s= ArgParseSettings()
@add_arg_table s begin
    "--config", "-c"
        help = "Config path"
        arg_type = String
        required = true
end
parsed_args = parse_args(ARGS, s)

config = parsed_args["config"]

create_models_folders(config)

