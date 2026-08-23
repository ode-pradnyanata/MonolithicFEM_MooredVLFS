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
using Statistics

# -------------------------------------------------------------
#   1  REFERENCES DATA  ---------------------------------------
# -------------------------------------------------------------

# -------------------------------------------------------------
#   2  FUNCTIONS  ---------------------------------------------
# -------------------------------------------------------------
# Func 1 ------------------------------------------------------
function Producing_Output(range_vpto, path_case, wave_number, wave_freq, mesh_size, damping_length_coeff, beam_length)
wave_number = wave_number
wave_freq = wave_freq
range_vpto = range_vpto
nRange_vpto = length(range_vpto)
nWave = length(wave_number)
comp_time = zeros(nRange_vpto, nWave)
for j = 1:nWave
  # ks = range_s[i]
  for i = 1:nRange_vpto
    ω_tmp = wave_freq[j]
    k_tmp = wave_number[j]
    comp_time[i,j] = @elapsed begin  
    case = Mooring_case_params(name="Michele-vpto$i-ω$j", v_pto=range_vpto[i], k_opt=k_tmp, ω_opt=ω_tmp,vtk_output=false, 
                                  h=mesh_size, cλ=damping_length_coeff, Lb=beam_length, ε=1e32, J=0, ny=10)
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

v_pto = 10 .^ range(2, 7, length=50)  # logscale
n_vpto = length(v_pto)

path=datadir("5-2-Comparison-Cw-Michele")

ω_data = collect(0.2: 0.05 : 6) 
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
L_data = 2*π ./ k_data
@show L_data
# #
time_init = @elapsed begin
case =  Mooring_case_params(
  name="Warm-up - Michele",
  h=1, # mesh size
  Lb=20,
  order=2,k_opt=k_data[35], ω_opt=ω_data[35], vtk_output=false, 
  cλ=2, v_pto = v_pto[8], ks=0, ks2=0, ε=1e32, J=0)
produce_or_load(path,case,run_Mooring;digits=8)
end
# 345.15 sec ~ 5.75 min
# #
# Partial running (it can be combined for the same cλ, if required)
# batch 1
ω_data = collect(0.20: 0.05 : 0.20)
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res_20m_1 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 6, 20)
res_10m_1 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 6, 10)  

ω_data = collect(0.25: 0.05 : 0.50) 
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data] 
res_20m_2 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 6, 20) 
res_10m_2 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 6, 10) 

ω_data = collect(0.55: 0.05 : 1.00) 
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data] 
res_20m_3 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 6, 20) 
res_10m_3 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 6, 10) 

# Regime 1 (0.20 - 1.00 rad/s)  
# 10m
# Total = 6.28 hr
# Range = (9.19 - 97.80 sec)
# Mean =  26.61 sec
total_10_r1 = (sum(res_10m_1) + sum(res_10m_2) + sum(res_10m_3))/3600
range_10_r1 = extrema(hcat(res_10m_1, res_10m_2, res_10m_3))
mean_10_r1  = mean(hcat(res_10m_1, res_10m_2, res_10m_3))

# 20m
# Total = 5.79 hr
# Range = (9.24 - 83.98 sec)
# Mean =  24.54 sec
total_20_r1 = (sum(res_20m_1) + sum(res_20m_2) + sum(res_20m_3))/3600
range_20_r1 = extrema(hcat(res_20m_1, res_20m_2, res_20m_3))
mean_20_r1  = mean(hcat(res_20m_1, res_20m_2, res_20m_3))


# ===============================================================================
# batch 2
ω_data = collect(1.05: 0.05 : 2.00) 
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res_20m_4 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 4, 20) 
res_10m_4 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 4, 10) 

ω_data = collect(2.05: 0.05 : 3.00)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data] 
res_20m_5 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 4, 20)  
res_10m_5 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 4, 10)  

# Regime 2 (1.05 - 3.00 rad/s)  
# 10m
# Total = 1.89 hr
# Range = (1.88 - 6.35 sec)
# Mean = 3.40 sec
total_10_r2 = (sum(res_10m_4) + sum(res_10m_5))/3600
range_10_r2 = extrema(hcat(res_10m_4, res_10m_5))
mean_10_r2  = mean(hcat(res_10m_4, res_10m_5))

# 20m
# Total = 1.99 hr
# Range = (2.06 - 6.54 sec)
# Mean = 3.58 sec
total_20_r2 = (sum(res_20m_4) + sum(res_20m_5))/3600
range_20_r2 = extrema(hcat(res_20m_4, res_20m_5))
mean_20_r2  = mean(hcat(res_20m_4, res_20m_5))

# ===============================================================================
# batch 3
ω_data = collect(3.05: 0.05 : 4.00)   
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data] 
res_20m_6 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 2, 20) 
res_10m_6 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 2, 10)

ω_data = collect(4.05: 0.05 : 5.00)
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data] 
res_20m_7 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 2, 20)
res_10m_7 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 2, 10) 

ω_data = collect(5.05: 0.05 : 6.00)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data] 
res_20m_8 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 2, 20)  
res_10m_8 = Producing_Output(v_pto, path, k_data, ω_data, 0.2, 2, 10) 

# Regime 3 (3.05 - 6.00 rad/s)  
# 10m
# Total = 1.43 hr
# Range = (1.35 - 2.28 sec)
# Mean = 1.72 sec
total_10_r3 = (sum(res_10m_6) + sum(res_10m_7) + sum(res_10m_8))/3600
range_10_r3 = extrema(hcat(res_10m_6, res_10m_7, res_10m_8))
mean_10_r3  = mean(hcat(res_10m_6, res_10m_7, res_10m_8))

# 20m
# Total = 2.04 hr
# Range = (2.04 - 2.60 sec)
# Mean = 2.44 sec
total_20_r3 = (sum(res_20m_6) + sum(res_20m_7) + sum(res_20m_8))/3600
range_20_r3 = extrema(hcat(res_20m_6, res_20m_7, res_20m_8))
mean_20_r3  = mean(hcat(res_20m_6, res_20m_7, res_20m_8))


# 8) 117 x 50 L=10m - 9.61hr
#    117 x 50 L=20m - 9.82hr
#    ω = 0.20 - 1.00, cλ = 6
#    ω = 1.05 - 3.00, cλ = 4
#    ω = 3.05 - 6.00, cλ = 2

total_10_hr = total_10_r1 + total_10_r2 + total_10_r3
total_20_hr = total_20_r1 + total_20_r2 + total_20_r3

# Collect data 
result = collect_results(path)
df = DataFrame(result)

# Collecting Output
ω_data = collect(0.2: 0.05 : 6)
n_wavedata = length(ω_data)
result_Cw_10 = zeros(n_vpto,n_wavedata) # row - col
result_Cw_20 = zeros(n_vpto,n_wavedata) # row - col
for i = 1:n_vpto
  for j = 1:n_wavedata
    result_data_20 = filter(row -> row.v_pto == round(v_pto[i], digits=8) && row.ω_opt == round(ω_data[j],digits=8) && row.Lb == 20, df)
    result_Cw_20[i,j] = result_data_20.Cw[1]

    result_data_10 = filter(row -> row.v_pto == round(v_pto[i], digits=8) && row.ω_opt == round(ω_data[j],digits=8) && row.Lb == 10, df)
    result_Cw_10[i,j] = result_data_10.Cw[1]
  end
end

# Plots
default(
  tickfontsize=13,
  legendfontsize=13,
  labelfontsize=16,
  titlefontsize=16)
plt_10= heatmap(
    ω_data,
    v_pto,
    result_Cw_10;
    yscale = :log10,
    yticks = (10 .^ (2:7), ["10²","10³","10⁴","10⁵","10⁶","10⁷"]),
    levels = 30,
    color = :turbo,
    clims = (0.0, 1.0),
    xlabel = L"\omega\;(\mathrm{rad\ s^{-1}})",
    ylabel = L"\nu_{PTO}\;(\mathrm{kg\ m^{-1}\ s^{-1}})",
    colorbar_title = L"C_w",
    colorbar_titlefontsize = 16,   
    title=L"\textbf{Capture-Width \: Ratio \: (C_w)} , Lb = 10m",
    size = (700,600),
    framestyle = :box,
)
savefig(plt_10, plotsdir("5-2-Comparison-Cw","Cw-heatmap-10m"))

# contourf(
plt_20= heatmap(
    ω_data,
    (v_pto),
    result_Cw_20;
    yscale = :log10,
    yticks = (10 .^ (2:7), ["10²","10³","10⁴","10⁵","10⁶","10⁷"]),
    levels = 30,
    color = :turbo,
    clims = (0.0, 1.0),
    xlabel = L"\omega\;(\mathrm{rad\ s^{-1}})",
    ylabel = L"\nu_{PTO}\;(\mathrm{kg\ m^{-1}\ s^{-1}})",
    colorbar_title = L"C_w",
    colorbar_titlefontsize = 16,  
    title=L"\textbf{Capture-Width \: Ratio \: (C_w)} , Lb = 20m",
    size = (700,600),
    framestyle = :box,
)
savefig(plt_20, plotsdir("5-2-Comparison-Cw","Cw-heatmap-20m"))

