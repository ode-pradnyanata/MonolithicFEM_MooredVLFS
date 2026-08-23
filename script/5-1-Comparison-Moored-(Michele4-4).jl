include("..\\src\\1-2-research-mooring-michele (4-4).jl")

import .Mooring_Case: Mooring_case_params, run_mooring_case_freq_domain, collect_results, run_Mooring

using Plots
using Parameters
using DrWatson
using JLD2
using DataFrames
using CSV
using Roots
using Statistics

# -------------------------------------------------------------
#   1  REFERENCES DATA  ---------------------------------------
# -------------------------------------------------------------
Michele_Left_1e3 = DataFrame(CSV.File(datadir("Ref_data/Michele","Michele_Left_ks_1e3_Flex.csv"); header=false))
Michele_Left_1e4 = DataFrame(CSV.File(datadir("Ref_data/Michele","Michele_Left_ks_1e4_Flex.csv"); header=false))
Michele_Left_1e5 = DataFrame(CSV.File(datadir("Ref_data/Michele","Michele_Left_ks_1e5_Flex.csv"); header=false))
Michele_Right_1e3 = DataFrame(CSV.File(datadir("Ref_data/Michele","Michele_Right_ks_1e3_Flex.csv"); header=false))
Michele_Right_1e4 = DataFrame(CSV.File(datadir("Ref_data/Michele","Michele_Right_ks_1e4_Flex.csv"); header=false))
Michele_Right_1e5 = DataFrame(CSV.File(datadir("Ref_data/Michele","Michele_Right_ks_1e5_Flex.csv"); header=false))
sort!(Michele_Left_1e5, :Column1)
sort!(Michele_Right_1e5, :Column1)

# -------------------------------------------------------------
#   2  FUNCTIONS  ---------------------------------------------
# -------------------------------------------------------------
# Func 1 ------------------------------------------------------
# function Producing_Output(range_s, path_case, hp, wave_number, wave_freq)
function Producing_Output(range_s, path_case, hp, wave_number, wave_freq, mesh_size, damping_length_coeff)
wave_number = wave_number
wave_freq = wave_freq
range_s = range_s
nRange_s = length(range_s)
nWave = length(wave_number)
comp_time = zeros(nWave, nRange_s)
for i = 1:nRange_s
  # ks = range_s[i]
  for j = 1:nWave
    ω_tmp = wave_freq[j]
    k_tmp = wave_number[j]
    comp_time[j,i] = @elapsed begin
    case = Mooring_case_params(name="Michele-ks$i-$hp-ω$j", ks=range_s[i], hₚ=hp, k_opt=k_tmp, ω_opt=ω_tmp,
            vtk_output=false, Lb=40, h=mesh_size, cλ=damping_length_coeff, ε = 1e32, J = 0.0, ny=10) 
    data, file = produce_or_load(path_case,case,run_Mooring;digits=8)
    end
    @show comp_time[j,i]
  end
end
return comp_time
end

# Func 2 ------------------------------------------------------
marker_style = [:square, :o, :utriangle]
color_scale = palette(:rainbow,5)
n_marks = 25

function plotting_freq(location,hp, wave_freq,range_s)
nWave = length(wave_freq)
nRange_s = length(range_s)
idx_marker = round.(Int, range(1, nWave, length=n_marks))

default(
  tickfontsize=13,
  legendfontsize=13,
  labelfontsize=16,
  titlefontsize=16)
plt1=plot(legend=:best, size = (700,600), palette=:rainbow, 
            grid=true, gridlinewidth=1, gridalpha=0.2)           
for i = 1:nRange_s
# for i = 2:3
  x_bin = zeros(nWave)
  η_bin = zeros(nWave)
  for j = 1:nWave
    name="Michele-ks$i-$hp-ω$j"
    ks_select = range_s[i]
    wave_select = wave_freq[j]
    result_data = filter(row -> row.ω_opt == wave_select && row.ks == ks_select, df)
    η = result_data.η
    if location == "start"
      η_bin[j]=η[1][1]
    elseif location == "end"
      η_bin[j]=η[1][end]
    end
  end
  plot!(plt1,wave_freq,η_bin, line=(2, color_scale[i]), #lw=3,palette=:rainbow,
      label=false)
  plot!(plt1,wave_freq[idx_marker],η_bin[idx_marker], marker=(4, color_scale[i], marker_style[i]),        #lw=3,palette=:rainbow,
      label=false, line=false)
  plot!([], [], line=(2,color_scale[i]), marker=(marker_style[i],4,color_scale[i]), label="ks=$(label_ks[i]) kg m⁻¹ s⁻²")
end

if location == "start"
  xx = "-L (Left Boundary)"
  #= 
    In Michele (2024), the wave model is propagating to the left
    Therefore, to match with this model and minimum adjustment, 
    we compare the left boundary (incoming wave) with Michele results 
    at the Right boundary, and vice versa
  =#
  plot!(plt1, Michele_Right_1e3.Column1, Michele_Right_1e3.Column2, line=(1,:black,:dash), 
        label="Michele et al. ks = 1e3 kg m⁻¹ s⁻²")
  plot!(plt1, Michele_Right_1e4.Column1, Michele_Right_1e4.Column2, line=(1,:black,:dash), 
        label="Michele et al. ks = 1e4 kg m⁻¹ s⁻²")
  plot!(plt1, Michele_Right_1e5.Column1, Michele_Right_1e5.Column2, line=(1,:black,:dash), 
        label="Michele et al. ks = 1e5 kg m⁻¹ s⁻²")
elseif location == "end"
  xx = "L (Right Boundary)"
  plot!(plt1, Michele_Left_1e3.Column1, Michele_Left_1e3.Column2, line=(1,:black,:dash), 
        label="Michele et al. ks = 1e3 kg m⁻¹ s⁻²")
  plot!(plt1, Michele_Left_1e4.Column1, Michele_Left_1e4.Column2, line=(1,:black,:dash), 
        label="Michele et al. ks = 1e4 kg m⁻¹ s⁻²")
  plot!(plt1, Michele_Left_1e5.Column1, Michele_Left_1e5.Column2, line=(1,:black,:dash), 
        label="Michele et al. ks = 1e5 kg m⁻¹ s⁻²")
end
# vline!([0.8, 1.0], color=:black, lw=0.5, linestyle=:dash, label=false)
plot!(plt1, xlims=(0,6), ylims=(0,2.5))
xticks!(0:1:6)
yticks!(0:0.25:2.5)
xlabel!("ω [rad/s]")
ylabel!("|η|/η₀ [-]")
title!("Deflection at x = $xx", titlefont = font(12))
display(current())
plot_title = "5-1 Michele-loc at $xx, h=02"
savefig(plt1, plotsdir("5-1-Comparison-Moored-plot",plot_title))
end

# Func 3 ------------------------------------------------------
function plotting_plate(hp, wave_freq_data,wave_freq_target,mooring_stiff)
plt1=plot(legendfontsize=10, legend=:best, grid=true, gridlinewidth=1, gridalpha=0.2) 
idx = argmin(abs.(wave_freq_data .- wave_freq_target))
closest_freq = wave_freq_data[idx]
max_η = zeros(length(mooring_stiff))
wave_select = wave_freq_data[idx]
n_mooring = length(mooring_stiff)

for i = 1:n_mooring
  ks_select = mooring_stiff[i]
  name="Michele-ks$i-$hp-ω$idx"
  result_data = filter(row -> row.ω_opt == wave_select && row.ks == ks_select, df)
  # result_data = filter(row -> row.name == name, df)
  x = result_data.x; η = result_data.η

  plot!(plt1,x,η,lw=2,palette=:rainbow,
          label="Mooring-ks=$(label_ks[i]) kg m⁻¹ s⁻²")  
  max_η[i] = maximum(η[1][:])
end
@show (max_η)
if maximum(max_η) < 1.5
  y_limit = 1.5
  xticks!(0:5:40);  yticks!(0:0.25:y_limit)
else
  range_lim = round(maximum(max_η) / 5, digits=2)
  y_limit = maximum(max_η) + range_lim + 0.1
  xticks!(0:5:40);  yticks!(0:range_lim :y_limit)
end
# vline!([round(40/3, digits = 1, base = Int(1/0.01)), round(40/3*2, digits = 1, base = Int(1/0.01))],color=:black, lw=1, linestyle=:dash, label="Joint locations")
plot!(plt1, xlims=(0,40), ylims=(0,y_limit))
xlabel!("x [m]")
ylabel!("|η|/η₀ [-]")
title!("Plate Deflection
hₚ = $hp m, ω = $(round(closest_freq,digits=4)) rad/sec", 
titlefont = font(12))
plot_title = "Michele-Plate Deflection-ω$idx"
savefig(plt1, plotsdir("5-1-Comparison-Moored-plot",plot_title))
end

# -------------------------------------------------------------
#   3  RUN CASE  ----------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Wave Properties
g = 9.81
h₀ = 10
ω_data = collect(0.01: 0.01 : 6)  # From near 0 to 6 rad/s
@show ω_data[1:10]

# Checking marked wave
# idx_check = round.(Int, range(1, length(ω_data), length=25))
# wave_xx = ω_data[idx_check]; @show wave_xx

function solve_dispersion(ω, g, h₀)
    f(k) = ω^2 - g * k * tanh(k * h₀)
    k_initial = ω^2 / g       # Initial guess for k - Deep water approximation plus small offset
    k_solution = find_zero(f, k_initial, Order2())
    return k_solution
end

k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
T_data = 2π ./ ω_data

mooring_stiff = [1e3, 1e4, 1e5]
label_ks = ["1e3","1e4","1e5"]

path=datadir("5-1-Michele-Moored-Freq")

time_init = @elapsed begin
case =  Mooring_case_params(
  name="Warm-up - Michele",
  h=1, # mesh size
  Lb=10,
  order=2,k_opt=k_data[200], ω_opt=ω_data[200], ks=0, vtk_output=false, cλ=4, ε=1e32, J=0.0)
produce_or_load(path,case,run_Mooring;digits=8)
end 
# Start with lower order warm-up case. We use h=1, Lb=10, and order=2
# Running time for this warm up is about 360 second ~ 6 minutes

# # -------------
# ====================================================================================
# Very low frequency case
ω_data = collect(0.1: 0.01: 0.2)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res1 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 8) 
total_time_res1_hr = sum(res1)/3600
avg_time_res1_sec = mean(res1)
# For very low frequency, the time is much longer (range from lower to higher freq, 490 - 140 sec)
# total   = 1.80 hr 
# average = 196.8 sec

ω_data = collect(0.21: 0.01: 0.50)
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res2 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 8)
total_time_res2_hr = sum(res2)/3600
avg_time_res2_sec = mean(res2)
# total   =  1.34 hr 
# average =  53.7 sec

ω_data = collect(0.51: 0.01: 0.80)
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res3 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 8)
total_time_res3_hr = sum(res3)/3600
avg_time_res3_sec = mean(res3)
# total   =  0.52 hr 
# average =  20.87 sec

# ====================================================================================
ω_data = collect(0.81: 0.01: 1.00)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res4 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 6)
total_time_res4_hr = sum(res4)/3600
avg_time_res4_sec = mean(res4)
# total   =  0.15 hr 
# average =  9.29 sec

# ====================================================================================
ω_data = collect(1.01: 0.01: 2.00)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res5 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 4)
total_time_res5_hr = sum(res5)/3600
avg_time_res5_sec = mean(res5)
# total   =  0.31 hr 
# average =  3.79 sec

ω_data = collect(2.01: 0.01: 3.00)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res6 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 4)
total_time_res6_hr = sum(res6)/3600
avg_time_res6_sec = mean(res6)
# total   =  0.20 hr 
# average =  2.45 sec

ω_data = collect(3.01: 0.01: 4.00)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res7 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 4)
total_time_res7_hr = sum(res7)/3600
avg_time_res7_sec = mean(res7)
# total   =  0.19 hr 
# average =  2.37 sec

ω_data = collect(4.01: 0.01: 4.99)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
res8 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 4)
total_time_res8_hr = sum(res8)/3600
avg_time_res8_sec = mean(res8)
# total   =  0.19 hr 
# average =  2.37 sec

ω_data = collect(5: 0.01: 6)  
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
L_data = 2π ./ k_data
res9 = Producing_Output(mooring_stiff, path, 0.1, k_data, ω_data, 0.2, 4) 
total_time_res9_hr = sum(res9)/3600
avg_time_res9_sec = mean(res9)
# total   =  0.20 hr 
# average =  2.37 sec



# ====================================================================================
# ====================================================================================
# Per Regime
# Regime 1 (0.00 - 0.80 rad/s) = res1 + res2 + res3
# Total = 3.66hr
# Range = (14.46 - 392.34 sec)
# Mean = 61.99 sec
total_r1 = total_time_res1_hr + total_time_res2_hr + total_time_res3_hr
range_r1 = extrema(vcat(res1, res2, res3))
mean_r1 = mean(vcat(res1, res2, res3))

# Regime 2 (0.81 - 1.00 rad/s) = res4
# Total = 0.15hr
# Range = (7.32 - 11.54 sec)  
# Mean = 9.29 sec
total_r2 = total_time_res4_hr
range_r2 = extrema(res4)
mean_r2 = mean(res4)

# Regime 3 (1.01 - 6.00 rad/s) = res5 + res6 + res7 + res8 + res9
# Total = 1.11hr
# Range = (1.84 - 6.14 sec)
# Mean = 2.67 sec
total_r3 = total_time_res5_hr + total_time_res6_hr + total_time_res7_hr + total_time_res8_hr + total_time_res9_hr
range_r3 = extrema(vcat(res5, res6, res7, res8, res9))
mean_r3 = mean(vcat(res5, res6, res7, res8, res9))


# #
result = collect_results(path) 
df = DataFrame(result)

ω_data = collect(0.1: 0.01: 6.00) 
plotting_freq("start",0.1,ω_data, mooring_stiff)
plotting_freq("end",0.1,ω_data, mooring_stiff)
plotting_plate(0.1, ω_data, 5.8, mooring_stiff)
