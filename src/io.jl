typedict(x) = Dict{Symbol, Any}(
    fn=>getfield(x, fn) for fn in fieldnames(typeof(x))
)

function saveinput(input, filename)
    filename = endswith(filename, ".json") ? filename : filename * ".json"
    open(filename, "w") do io
        JSON.print(io, typedict(input), 4)
    end
end

function loadinput(::Type{InputType}, filename) where InputType <: SimulationInput
    filename = endswith(filename, ".json") ? filename : filename * ".json"
    d = JSON.parsefile(filename)
    kwargs = Dict(Symbol(k) => v for (k, v) in d)
    return InputType(; kwargs...)
end
