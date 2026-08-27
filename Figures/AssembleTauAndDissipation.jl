using JLD2
using Printf

########################################################################
# Assemble Figures 4 and 5's raw data (tau_err.txt, dissipation.txt)   #
# from the per-(L, vDD_vT) outputs of the production runs:             #
#                                                                       #
#   Solvers/GetEigenvectorsWithDMRG/Data/mu1_Vdd<vDD_vT>_L<L>.jld2   (from ExcitedStateDMRG.jl)
#   Solvers/ComputeDissipation/Data/Dissipation_Vdd<vDD_vT>_L<L>.jld2 (from ComputeDissipation.jl)
#                                                                       #
# Run this from the Figures/ folder, after all 16 (L, vDD_vT) points  #
# have been run through SteadyStateDMRG.jl -> ExcitedStateDMRG.jl ->  #
# ComputeDissipation.jl. It writes:                                    #
#                                                                       #
#   "Figure 5 Data/tau_err.txt"                                        #
#   "Figure 5 Data/dissipation.txt"                                    #
#                                                                       #
# in the format Figure_4/prep_fig4.py and Figure_5/prep_fig5.py parse. #
#                                                                       #
########################################################################

const LS = [1, 3, 5, 7]
const VDDS = [1.1, 1.2, 1.3, 1.4]

const MU1_DIR = "../Solvers/GetEigenvectorsWithDMRG/Data/"
const DISS_DIR = "../Solvers/ComputeDissipation/Data/"
const OUT_DIR = "Figure 5 Data/"

function load_grid(dir, prefix, suffix)
    grid = fill(NaN, length(LS), length(VDDS))
    missing_points = Tuple{Int,Float64}[]
    for (i, L) in enumerate(LS), (j, V) in enumerate(VDDS)
        path = dir * prefix * "Vdd$(V)_L$(L)" * suffix
        if isfile(path)
            # mu1 and Dissipation are saved as whatever dmrg()/scalar() returned, which
            # can carry a numerically-tiny imaginary part even though both are physically
            # real (the master equation operator isn't Hermitian, so `ishermitian=false`
            # is passed throughout, and ITensors doesn't discard the imaginary part for
            # you). real(...) here is a no-op if the value was already real.
            grid[i, j] = real(load_object(path))
        else
            push!(missing_points, (L, V))
        end
    end
    return grid, missing_points
end

mu1_grid, mu1_missing = load_grid(MU1_DIR, "mu1_", ".jld2")
diss_grid, diss_missing = load_grid(DISS_DIR, "Dissipation_", ".jld2")

if !isempty(mu1_missing) || !isempty(diss_missing)
    println("Not all 16 production points are available yet -- missing:")
    for (L, V) in mu1_missing
        println("  mu1:          L=$L, vDD_vT=$V  (expected $(MU1_DIR)mu1_Vdd$(V)_L$(L).jld2)")
    end
    for (L, V) in diss_missing
        println("  dissipation:  L=$L, vDD_vT=$V  (expected $(DISS_DIR)Dissipation_Vdd$(V)_L$(L).jld2)")
    end
    error("Aborting -- run the missing (L, vDD_vT) points through the production " *
          "pipeline first, then re-run this script. No partial tau_err.txt/" *
          "dissipation.txt is written, since prep_fig4.py/prep_fig5.py both assert " *
          "the full 4x4 grid is present.")
end

tau_grid = 2.0 ./ mu1_grid   # <tau_err> = 2 / mu1

mkpath(OUT_DIR)

function write_table(path, header_note, grid)
    open(path, "w") do io
        println(io, "# ", header_note)
        println(io, "# L  ", join(VDDS, "  "))
        for (i, L) in enumerate(LS)
            @printf(io, "%d", L)
            for j in eachindex(VDDS)
                @printf(io, "  %.10g", grid[i, j])
            end
            println(io)
        end
    end
end

write_table(OUT_DIR * "tau_err.txt",
    "<tau_err>, in units of tau_0. Assembled by AssembleTauAndDissipation.jl from " *
    "Solvers/GetEigenvectorsWithDMRG/Data/mu1_Vdd*_L*.jld2 via <tau_err> = 2/mu1.",
    tau_grid)

write_table(OUT_DIR * "dissipation.txt",
    "Qdot, in units of k_B T / tau_0. Assembled by AssembleTauAndDissipation.jl from " *
    "Solvers/ComputeDissipation/Data/Dissipation_Vdd*_L*.jld2.",
    diss_grid)

println("Wrote $(OUT_DIR)tau_err.txt and $(OUT_DIR)dissipation.txt from all 16 production points.")
