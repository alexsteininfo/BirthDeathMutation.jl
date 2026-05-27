module BirthDeathMutation

using Distributions
using Statistics
using Random
using StatsBase
using DataFrames
using GLM
using Printf
using CSV
using JSON
using DelimitedFiles
using AbstractTrees
using MakieCore

export
SimpleTreeCell,
BinaryNode,
Subclone,
SimulationInput,
BranchingMoranInput,
BranchingInput,
MoranInput,
Simulation,
SampledData,
Population,
AbstractSelection,
SelectionPredefined,
SelectionDistribution,
NeutralSelection,

#functions for running simulations
runsimulation,
runsimulation_timeseries,
runsimulation_timeseries_returnfinalpop,
allcells,
initialize_population,
getsubclonesizes,

#statistics (tree-based)
pairwisedistance,
pairwisedistances,
pairwise_differences,
pairwise_fixed_differences,
pairwise_fixed_differences_matrix,
pairwise_fixed_differences_statistics,
pairwise_fixed_differences_clonal,
average_mutations,
mutations_per_cell,
clonal_mutations,
age,

#functions for tree cells
endtime,
celllifetime,
celllifetimes,
getalivecells,
leftchild!,
rightchild!,
time_to_MRCA,
coalescence_times,
getsingleroot,
popsize,
findMRCA,

#other
getclonetype,

#util
newinput,
saveinput,
get_simulation,
loadinput

include("input.jl")
include("selection.jl")
include("cells_modules.jl")
include("population.jl")
include("results.jl")
include("initialisation.jl")
include("cellupdates.jl")
include("simulations.jl")
include("run.jl")
include("process.jl")
include("util.jl")
include("statistics.jl")
include("io.jl")
include("simulation_trees.jl")


end
