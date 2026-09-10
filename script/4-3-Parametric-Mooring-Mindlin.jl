include("..\\src\\1-2-research-mooring-zhao.jl")

import .Mooring_Case: Mooring_case_params, run_mooring_case_freq_domain, collect_results, run_Mooring

using Plots
using Parameters
using DrWatson
using JLD2
using DataFrames
using CSV
using Printf
using Roots
using Statistics
using LaTeXStrings

# -------------------------------------------------------------
#   1  REFERENCES DATA  ---------------------------------------
# -------------------------------------------------------------

# -------------------------------------------------------------
#   2  RUN CASE  ----------------------------------------------
# -------------------------------------------------------------
h_mesh = 0.2                               # mesh size
ω_data = [0.5, 2.18546, 4.39691, 8.97598]  # From near 0 to 6 rad/s
c_damping = [8, 4, 4, 4]                   # Damping length coefficient - 8 for low frequency for accuracy
g = 9.81
h₀ = 1.1

function solve_dispersion(ω, g, h₀)
    f(k) = ω^2 - g * k * tanh(k * h₀)
    k_initial = ω^2 / g       # Initial guess for k - Deep water approximation plus small offset
    k_solution = find_zero(f, k_initial, Order2())
    return k_solution
end

k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
T_data = 2π ./ ω_data
L_data = 2π ./ k_data
@show L_data
mooring_stiff = [0, 1e3, 1e5]
# label_ks = ["0","1e3","1e5"]
label_ks = [L"\mathbf{0}", L"\mathbf{1\times10^3}", L"\mathbf{1\times10^5}"]
label_wave = ["a","b","c","d"]

path=datadir("4-3-Parametric-Mooring")

time_init = @elapsed begin
case =  Mooring_case_params(
  name="Warm-up",
  h=1, # mesh size
  order=2,k_opt=k_data[3], ω_opt=ω_data[3], ks=0, cλ=2)
produce_or_load(path,case,run_Mooring;digits=8)
end
# warmup 270.16 seconds ~ 4.5 minutes 

# Mooring Comparison
# comparison per frequencies - 3 mooring cases - use hp = 0.038 m
#   a, b, c - wave data
#   1, 2, 3, 4 - mooring data
nFreq = length(ω_data)
nMooring = length(mooring_stiff)
n_data = nFreq * nMooring
time = zeros(n_data)

for i = 1:nFreq
  for j = 1:nMooring
    n = (i-1)*nMooring + j
    labelling = label_wave[i]
    case = Mooring_case_params(name="5$labelling-$j", hₚ=0.038, ks=mooring_stiff[j], 
          k_opt=k_data[i], ω_opt=ω_data[i], h=h_mesh, vtk_output=false, cλ=c_damping[i])
    time[n] = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
end
end

# time = [3.7340205, 3.9920534, 3.7581324, 0.5114063, 0.5172937, 0.5199981, 0.3347736, 0.4065685, 0.820993, 0.2941373, 0.314209, 0.3146159]
# Low freq ~ 3 s
# Mean time = 1.293 s
# Mean time without low freq = 0.390 s

# -------------------------------------------------------------
#   3  OUTPUT DATA COLLECTION  --------------------------------
# -------------------------------------------------------------
result = collect_results(path)
df = DataFrame(result)
##
res1 = filter(row -> row.name == "5a-1" && row.h == h_mesh, df); xs1 = res1.x; η_xs1 = res1.η
res2 = filter(row -> row.name == "5a-2" && row.h == h_mesh, df); xs2 = res2.x; η_xs2 = res2.η
res3 = filter(row -> row.name == "5a-3" && row.h == h_mesh, df); xs3 = res3.x; η_xs3 = res3.η

res4 = filter(row -> row.name == "5b-1" && row.h == h_mesh, df); xs4 = res4.x; η_xs4 = res4.η
res5 = filter(row -> row.name == "5b-2" && row.h == h_mesh, df); xs5 = res5.x; η_xs5 = res5.η
res6 = filter(row -> row.name == "5b-3" && row.h == h_mesh, df); xs6 = res6.x; η_xs6 = res6.η

res7 = filter(row -> row.name == "5c-1" && row.h == h_mesh, df); xs7 = res7.x; η_xs7 = res7.η
res8 = filter(row -> row.name == "5c-2" && row.h == h_mesh, df); xs8 = res8.x; η_xs8 = res8.η
res9 = filter(row -> row.name == "5c-3" && row.h == h_mesh, df); xs9 = res9.x; η_xs9 = res9.η

res10 = filter(row -> row.name == "5d-1" && row.h == h_mesh, df); xs10 = res10.x; η_xs10 = res10.η
res11 = filter(row -> row.name == "5d-2" && row.h == h_mesh, df); xs11 = res11.x; η_xs11 = res11.η
res12 = filter(row -> row.name == "5d-3" && row.h == h_mesh, df); xs12 = res12.x; η_xs12 = res12.η

xa = [xs1, xs2, xs3]
ηa = [η_xs1, η_xs2, η_xs3]
xb = [xs4, xs5, xs6]
ηb = [η_xs4, η_xs5, η_xs6]
xc = [xs7, xs8, xs9]
ηc = [η_xs7, η_xs8, η_xs9]
xd = [xs10, xs11, xs12]
ηd = [η_xs10, η_xs11, η_xs12]

# -------------------------------------------------------------
#   4  PLOTTING  ----------------------------------------------
# -------------------------------------------------------------
marker_style = [:square, :o, :utriangle]
color_scale = palette(:rainbow,5)
n_marks = 21
idx_marker = round.(Int, range(1, length(xs1[1]), length=n_marks))
# xxx = xs1[1]
# @show xxx[idx_marker]

function plotting_mooring(x_data, η_data, ω_data, label_wave, label_ks)
  default(
    tickfontsize=13,
    legendfontsize=10,
    labelfontsize=16,
    titlefontsize=16)

  plot(legend=:best, size = (700,600), palette=:rainbow, 
        grid=true, gridlinewidth=1, gridalpha=0.2)
  for j = 1:3
    labelling_ks = label_ks[j]
    plot!(x_data[j], η_data[j],xlims=(0,10),line=(2, color_scale[j]),
          label=false)    
    plot!((x_data[j]...)[idx_marker], (η_data[j]...)[idx_marker],xlims=(0,10),
           line=false, marker=(4,marker_style[j], color_scale[j]), label=false)
  end

  plot!([], [], line=(2,color_scale[1]), marker=(marker_style[1],4,color_scale[1]), label=latexstring("ks= $(label_ks[1]) kg m⁻¹ s⁻²"))
  plot!([], [], line=(2,color_scale[2]), marker=(marker_style[2],4,color_scale[2]), label=latexstring("ks= $(label_ks[2]) kg m⁻¹ s⁻²"))
  plot!([], [], line=(2,color_scale[3]), marker=(marker_style[3],4,color_scale[3]), label=latexstring("ks= $(label_ks[3]) kg m⁻¹ s⁻²"))
  
  xlabel!("x [m]")
  ylabel!("|η|/κ₀ [-]")
  ylims!(0,1.5)
  xlims!(0,10)
  title!("ω = $(@sprintf("%.4g", ω_data)) rad/s")
  display(current())
  savefig(plotsdir("4-3-Parametric-Mooring","4-3 Parametric-mooring-$label_wave.png"))
end

plotting_mooring(xa, ηa, ω_data[1], label_wave[1], label_ks)
plotting_mooring(xb, ηb, ω_data[2], label_wave[2], label_ks)
plotting_mooring(xc, ηc, ω_data[3], label_wave[3], label_ks)
plotting_mooring(xd, ηd, ω_data[4], label_wave[4], label_ks)

@show time_init
@show time
@show mean(time)
@show mean(time[4:12])  # mean time without low frequency

# #
# Extract the VTK output
# We use the case with higher mooring stiffness
# case = Mooring_case_params(name="Extracting_VTK_4-3", hₚ=0.038, ks=mooring_stiff[3], 
#         k_opt=k_data[2], ω_opt=ω_data[2], h=h_mesh, vtk_output=true, cλ=c_damping[2])
# data, file = produce_or_load(path,case,run_Mooring;digits=8)

# case = Mooring_case_params(name="Extracting_VTK_4-3_no_mooring", hₚ=0.038, ks=mooring_stiff[1], 
#         k_opt=k_data[2], ω_opt=ω_data[2], h=h_mesh, vtk_output=true, cλ=c_damping[2])
# data, file = produce_or_load(path,case,run_Mooring;digits=8)