include("..\\src\\1-2-research-mooring-michele (4-1).jl")

import .Mooring_Case: Mooring_case_params, run_mooring_case_freq_domain, collect_results, run_Mooring

using Plots
using Parameters
using DrWatson
using JLD2
using DataFrames
using CSV
using Roots
using LaTeXStrings
using Printf

# -------------------------------------------------------------
#   1  REFERENCES DATA  ---------------------------------------
# -------------------------------------------------------------

# -------------------------------------------------------------
#   2  FUNCTIONS  ---------------------------------------------
# -------------------------------------------------------------
# Func 1 ------------------------------------------------------
function Producing_Output(range_vpto, path_case, wave_number, wave_freq, mesh_size, damping_length_coeff, beam_length, mooring_stiff)
wave_number = wave_number
wave_freq = wave_freq
range_vpto = range_vpto
nRange_vpto = length(range_vpto)
nMooring = length(mooring_stiff)
comp_time = zeros(nRange_vpto, nMooring)
for j = 1:nMooring
  ks_inp = mooring_stiff[j]
  for i = 1:nRange_vpto
    ω_tmp = wave_freq
    k_tmp = wave_number
    comp_time[i,j] = @elapsed begin  
    case = Mooring_case_params(name="Michele-vpto$i-ks$j", v_pto=range_vpto[i],k_opt=k_tmp, ω_opt=ω_tmp, vtk_output=false,
                                h=mesh_size, cλ=damping_length_coeff, Lb=beam_length, ε=1e32, J=0,
                                ks=ks_inp, ks2=ks_inp, ny=10)
    data, file = produce_or_load(path_case,case,run_Mooring;digits=8)
    end
    @show comp_time[i,j]
  end
end
return comp_time
end


# -------------------------------------------------------------
#   3  RUN CASE  ----------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Wave Properties
function solve_dispersion(ω, g, h₀)
    f(k) = ω^2 - g * k * tanh(k * h₀)
    k_initial = ω^2 / g       # Initial guess for k - Deep water approximation plus small offset
    k_solution = find_zero(f, k_initial, Order2())
    return k_solution
end
g = 9.81
h₀ = 10

v_pto = 10 .^ range(0, 7, length=50)  # logscale
n_vpto = length(v_pto)
# mooring_stiff = [0, 1e3, 1e5]
mooring_stiff = [0, 1e3, 1e4, 1e5]
label_ks = [L"\mathbf{0}", L"\mathbf{1\times10^3}", L"\mathbf{1\times10^4}", L"\mathbf{1\times10^5}"]
n_mooring = length(mooring_stiff)
path=datadir("5-3-Parametric-Cw-10m")

ω_data = collect(2: 0.05 : 2) 
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
L_data = 2*π ./ k_data
@show L_data

time_init = @elapsed begin
case =  Mooring_case_params(
  name="Warm-up - Michele",
  h=1, # mesh size
  Lb=10,
  order=2,k_opt=k_data[1], ω_opt=ω_data[1], vtk_output=false, 
  cλ=2, v_pto = v_pto[1], ks=0, ks2=0, ε=1e32, J=0)
produce_or_load(path,case,run_Mooring;digits=8)
end
# 350.63 sec ~ 5.75 min

# #
# Partial running (it can be combined for the same cλ, if required)
# batch 1
res_1 = Producing_Output(v_pto, path, k_data[1], ω_data[1], 0.2, 4, 10, mooring_stiff)
res_2 = Producing_Output([0], path, k_data[1], ω_data[1], 0.2, 4, 10, [0])

# Collect data 
result = collect_results(path)
df = DataFrame(result)

# Collecting Output
# ω_data = collect(0.2: 0.05 : 6)
# n_wavedata = length(ω_data)
result_Cw_10 = zeros(n_vpto,n_mooring) # row - col
result_η_10 = zeros(n_vpto,n_mooring)
for i = 1:n_vpto
  for j = 1:n_mooring
    # result_data_20 = filter(row -> row.v_pto == round(v_pto[i], digits=8) && row.ω_opt == round(ω_data[j],digits=8) && row.Lb == 20, df)
    # result_Cw_20[i,j] = result_data_20.Cw[1]

    result_data_10 = filter(row -> row.v_pto == round(v_pto[i], digits=8) && row.ks == round(mooring_stiff[j],digits=8) && row.Lb == 10, df)
    result_Cw_10[i,j] = result_data_10.Cw[1]
  end
end

Cw_1 = result_Cw_10[:,1]
Cw_2 = result_Cw_10[:,2]
Cw_3 = result_Cw_10[:,3]
Cw_4 = result_Cw_10[:,4]
idx_max = argmax(Cw_1)
pto_max = v_pto[idx_max]

# Plots
for i = 1:1
default(
  tickfontsize=13,
  legendfontsize=12,
  labelfontsize=16,
  titlefontsize=16)  
plot(xscale=:log10, xticks = (10 .^ (0:7), ["10⁰","10¹","10²","10³","10⁴","10⁵","10⁶","10⁷"]),
     yticks = (0:0.2:1.2, ["0.0","0.2","0.4","0.6","0.8","1.0","1.2"]),
     legend=:best, size = (750,600), palette=:rainbow, 
     grid=true, gridlinewidth=1, gridalpha=0.2)
plot!(v_pto, Cw_1, lw=3, label="ks = $(label_ks[1]) kg m⁻¹ s⁻²")
plot!(v_pto, Cw_2, lw=3, label="ks = $(label_ks[2]) kg m⁻¹ s⁻²")
plot!(v_pto, Cw_3, lw=3, label="ks = $(label_ks[3]) kg m⁻¹ s⁻²")
plot!(v_pto, Cw_4, lw=3, label="ks = $(label_ks[4]) kg m⁻¹ s⁻²")
vline!([pto_max], line=(2, :dash, :black), label=false)
xlims!(1e0, 2e7)
ylims!(0, 1)
ylabel!(L"C_w \: [-]")
xlabel!(L"v_{PTO} \: [\mathrm{kg \cdot m^{-1} \cdot s^{-1}}]")
title!("Capture-width ratio for ω = $(@sprintf("%.4g", ω_data[1])) rad/s")
display(current())
savefig(plotsdir("5-3-Parametric-Cw-10m","5-3a Mooring-effect-on-Cw.png"))
end



label_vpto = [L"\mathbf{0}", L"\mathbf{3.7\times10^3}"]
for i =1:1
default(
  tickfontsize=13,
  legendfontsize=12,
  labelfontsize=16,
  titlefontsize=16)    
plot(legend=:top, size = (750,600), palette=:rainbow,
     grid=true, gridlinewidth=1, gridalpha=0.2)

result_defl_0 = filter(row -> row.v_pto == round(0, digits=8) && row.ks == round(0, digits=8)  && row.Lb == 10, df)
η_def0  = result_defl_0.η[1]; x_def0  = result_defl_0.x[1]; Cw_def0 = result_defl_0.Cw[1]
plot!(x_def0, η_def0,lw=3, label="vₚₜₒ = $(label_vpto[1]), ks = $(label_ks[1])")

result_defl_1 = filter(row -> row.v_pto == round(pto_max, digits=8) && row.ks == round(0,digits=8) && row.Lb == 10, df)
η_defl = result_defl_1.η[1]; x_defl = result_defl_1.x[1]; Cw_def1 = result_defl_1.Cw[1]
plot!(x_defl, η_defl,lw=3, label="vₚₜₒ = $(label_vpto[2]), ks = $(label_ks[1])")

result_defl_2 = filter(row -> row.v_pto == round(pto_max, digits=8) && row.ks == round(1e4,digits=8) && row.Lb == 10, df)
η_def2 = result_defl_2.η[1]; x_def2 = result_defl_2.x[1]; Cw_def2 = result_defl_2.Cw[1]
plot!(x_def2, η_def2,lw=3, label="vₚₜₒ = $(label_vpto[2]), ks = $(label_ks[2])")

result_defl_3 = filter(row -> row.v_pto == round(pto_max, digits=8) && row.ks == round(1e5,digits=8) && row.Lb == 10, df)
η_def3 = result_defl_3.η[1]; x_def3 = result_defl_3.x[1]; Cw_def3 = result_defl_3.Cw[1]
plot!(x_def3, η_def3,lw=3, label="vₚₜₒ = $(label_vpto[2]), ks = $(label_ks[3])")

annotate!([
  (1.2, 1.25, text("Cw = $(@sprintf("%.4f", Cw_def0))", 14)),
  (1.2, 0.70, text("Cw = $(@sprintf("%.4f", Cw_def1))", 14)),
  (1.2, 0.38, text("Cw = $(@sprintf("%.4f", Cw_def2))", 14)),
  (1.2, 0.18, text("Cw = $(@sprintf("%.4f", Cw_def3))", 14))])
 
ylabel!("|η|/κ₀ [-]")
xlabel!("x [m]")
title!("Hydroelastic response for ω = $(@sprintf("%.4g", ω_data[1])) rad/s")
ylims!(0,1.5)
display(current())
savefig(plotsdir("5-3-Parametric-Cw-10m","5-3b Mooring-effect-on-hydroelastic response.png"))
end


