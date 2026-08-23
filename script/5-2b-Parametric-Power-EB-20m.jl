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

# 2
function Producing_Output_freq(range_vpto, path_case, wave_number, wave_freq, mesh_size, damping_length_coeff, beam_length, mooring_stiff)
wave_number = wave_number
wave_freq = wave_freq
range_vpto = range_vpto
nRange_vpto = length(range_vpto)
nMooring = length(mooring_stiff)
nWave = length(wave_freq)
comp_time = zeros(nRange_vpto, nWave)
for j = 1:nWave
  ω_tmp = wave_freq[j]
  k_tmp = wave_number[j]  
  ks_inp = mooring_stiff[1]
  for i = 1:nRange_vpto
    comp_time[i,j] = @elapsed begin  
    case = Mooring_case_params(name="Michele-vpto$i-ω$j", v_pto=range_vpto[i],k_opt=k_tmp, ω_opt=ω_tmp, vtk_output=false,
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
n_mooring = length(mooring_stiff)
path=datadir("5-2-Parametric-Power-20m")

ω = 2
ω_data = collect(ω: 0.05 : ω) 
Lb_length = 20
# ω_data = 6 
k_data = [solve_dispersion(ω, g, h₀) for ω in ω_data]
L_data = 2*π ./ k_data
@show k_data[1]

time_init = @elapsed begin
case =  Mooring_case_params(
  name="Warm-up - Michele",
  h=1, # mesh size
  Lb=Lb_length,
  order=2,k_opt=k_data[1], ω_opt=ω_data[1], vtk_output=false, 
  cλ=2, v_pto = v_pto[1], ks=0, ks2=0, ε=1e32, J=0)
produce_or_load(path,case,run_Mooring;digits=8)
end
# 350.63 sec ~ 5.75 min

# -----------------------------------------------------
# -----------------------------------------------------

# INPUT 
v_pto_check = [1e3, 2e6]
ks_check = [0]
ω_check = collect(1: 0.01 : 10) 
k_check = [solve_dispersion(ω, g, h₀) for ω in ω_check]

# Producing output for the specific case
for v_pto_eval in v_pto_check
  Producing_Output_freq(v_pto_eval, path, k_check, ω_check, 0.2, 4, Lb_length, ks_check)
end

result_check = collect_results(path)
df_check = DataFrame(result_check)

κ₀ = 1
H  = 10
ρ  = 1000
Cg = (ω_check ./ (2 .* k_check)) .* (1 .+ (2 .* k_check .* H)./(sinh.(2 .* k_check .* H))) 
Pw = 0.5 * ρ *g * κ₀^2 .* Cg 


# function Plotting_Power_Cw(v_pto_eval)
# for i =1:1
begin
idx_pto = [1, 85, 165, 250, 335, 415, 500]    # for Lb = 20, h = 0.2
nPTO    = length(idx_pto)
nWave   = length(ω_check)
nRange_pto = length(v_pto_check)
total_P_pto = zeros(Float64, nRange_pto, nWave)
Cw_calc = zeros(Float64, nRange_pto, nWave)
Cw_run = zeros(Float64, nRange_pto, nWave) 

n_marks = 20
idx_marker = round.(Int, range(1, nWave, length=n_marks))

default(
  tickfontsize=13,
  legendfontsize=13,
  labelfontsize=16,
  titlefontsize=16)


for i = 1:nRange_pto
  plt1=plot(legend=:best, size = (700,600), palette=:rainbow, 
            grid=true, gridlinewidth=1, gridalpha=0.2)

  x_bin = zeros(nWave)
  η_bin = zeros(nWave)
  η_pto = zeros(nPTO)
  vpto_select = v_pto_check[i]
  
  for j = 1:nWave
    name="Michele-vpto$i-ω$j"
    wave_select = ω_check[j]
    result_data = filter(row -> row.ω_opt == wave_select && row.v_pto == vpto_select, df_check)
    η_res  = result_data.η[1]
    Cw_res = result_data.Cw[1]
    for k = 1:nPTO
      idx = idx_pto[k]
      η_pto[k] = η_res[idx]
    end

    Cw_run[i,j] = Cw_res

    # Manual calculation to check the results
    P_pto = 0.5 .* vpto_select .* wave_select^2 .* η_pto.^2
    total_P_pto[i,j] = sum(P_pto)
    Cw_calc[i,j] = total_P_pto[i,j] ./ Pw[j]

  end
  plot!(plt1,ω_check,total_P_pto[i,:], line=(2, :rainbow, :solid), #lw=3,palette=:rainbow, 
      label="Pₚₜₒ, vₚₜₒ=$(vpto_select) kg m⁻¹s⁻¹")
  # plot!(plt1,ω_check,total_P_pto[i,:],lw=3,palette=:rainbow,
  #     label=latexstring("P_{pto}, v_{pto} = $(@sprintf("%.2e", vpto_select)) \\, \\mathrm{Nm^{-1}}"))


# end

  ytick_vals = range(0, maximum(Pw), length=6)
  ytick_labels = [@sprintf("%.2f", v/1e3) for v in ytick_vals]
  plot!(plt1,ω_check,Pw, line=(3, :dot), palette=:rainbow,   #lw=3,
      label="Wave energy flux", legend=(0.45, 0.90),
      yticks = (ytick_vals, ytick_labels),
      ylabel = L"Wave Power $(\times 10^3)$")
  annotate!(plt1,
    1.1,          # x position (left edge)
    ylims()[2],          # y position (top)
    text(L"\times 10^3", :left, 16))

  plt2= twinx(plt1)      
  # for i = 1:nRange_pto
  # vpto_select = v_pto_check[i]
  plot!(plt2,ω_check,Cw_calc[i,:],lw=2,color=:red,
      label=false)
  plot!(plt2,legend=(0.45, 0.750))  
  # end
  plot!(plt2,ω_check[idx_marker],Cw_calc[i, idx_marker],line=false, marker=(:square, :red, 4),
      label=false)
  plot!(plt2,[], [], line=(2,:red), marker=(:square, :red, 4), label="Cw, vₚₜₒ=$(v_pto_check[i]) kg m⁻¹s⁻¹")  


  plot!(plt1, xlims=(1,(ω_check[end])), ylims=(0,maximum(Pw))) 
  plot!(plt2, ylims=(0,1)) 
  xlabel!(plt1,"ω [rad/s]")
  ylabel!(plt1,"Wave Power [W/m]")
  ylabel!(plt2,"Cw")
  xlabel!(plt2,"")
  title!("Parametric Power")
  display(current())

  plot_title = "Parametric power calculation - vpto$(Int(v_pto_check[i]))"
  savefig(plotsdir("5-2-Parametric-Power-20m",plot_title))

end

end


