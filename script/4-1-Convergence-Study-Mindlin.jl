include("..\\src\\1-2-research-mooring-zhao-conv.jl")

import .Mooring_Case: Mooring_case_params, run_mooring_case_freq_domain, collect_results, run_Mooring

using Plots
using Parameters
using DrWatson
using JLD2
using DataFrames
using CSV
using Printf
using LaTeXStrings
using Gridap
using Gridap.CellData
using Interpolations
using Roots
using LaTeXStrings

# -------------------------------------------------------------
#   1  RUN CASE  ----------------------------------------------
# -------------------------------------------------------------
# beam_length = 10
k_data = [0.8044 9.0434]                            # [C. Zhao] - dimensionless wave number = k[1/m] * H [m]
T_data = [2.8750 0.7000]   
ω_data = 2π ./ T_data
hp_data = [0.038, 0.075, 0.100]                             # [C. Zhao] plate thickness

path=datadir("4-1-Convergence-Study-Mindlin")

# Warm-up2 case
time_init = @elapsed begin
case =  Mooring_case_params(
  name="Warm-up",
  h=2, # mesh size
  order=2,k_opt=k_data[2], ω_opt=ω_data[2], ks=0, cλ=2, Lb=10, vtk_output=true)
# ~, ~, ~ = run_mooring_case_freq_domain(case)
~, ~, ~, ~ = run_mooring_case_freq_domain(case)
end


# #
# function run_convergence_L2(hp, element_order, order_ref, ny_elem, beam_length, k_wave, ω_wave)
function run_convergence_L2(element_order)
  # Initialize bins
  L2_err_bin = zeros(ComplexF64, n_mesh_size) 
  Linf_err_bin = zeros(ComplexF64, n_mesh_size)
  time_bin = zeros(n_mesh_size)

  for i = 1:(n_mesh_size)
    # 1. Running the case for each mesh size
    time_init = @elapsed begin
      case = Mooring_case_params(name="Convergence-Study-coarse-$element_order-$i", hₚ = hp, H=water_depth, Lb=beam_length, cλ=4,
             order=element_order, h=mesh_size_var[i], ny=ny_data[i], k_opt=k_check, ω_opt=ω_check, vtk_output=true) 
      ΓbH, ηH, U_ΓηH, V_ΓηH = run_mooring_case_freq_domain(case) 
    end
  
  # Method 2: fine to coarse
    # interp_ηₕ_f = Interpolable(ηₕ_f)                        # wrap fine solution
    # η_ref_on_coarse = interpolate_everywhere(interp_ηₕ_f, U_Γη)
    # e = ηₕ - η_ref_on_coarse
    # degree = 2 * element_order
    # dΓb  = Measure(Γb, degree)
    # L2_error = sqrt((sum( ∫( e ⊙ conj(e) )dΓb )))
    # Linf_error = maximum(abs.(get_free_dof_values(ηₕ) .- get_free_dof_values(η_ref_on_coarse)))

  # Method 2: coarse to fine
    interp_ηH = Interpolable(ηH; tol=1e-12)                        # wrap coarser solution
    ηH_on_h = interpolate_everywhere(interp_ηH, U_Γη_f)
    # e = ηₕ_f - ηH_on_h
    e = (ηH_on_h - ηₕ_f)
    degree = 2 * r_ref
    dΓbf  = Measure(Γb_f, degree)
    L2_error = sqrt(real(sum( ∫( e ⊙ conj(e) )dΓbf )))
    #  L2_error = √((∑( ∫( (e ⋅ e) )dΓbf )))                 # L2 norm function  
    Linf_error = maximum(abs.(get_free_dof_values(ηₕ_f) .- get_free_dof_values(ηH_on_h)))

  # 5. Store the 
    L2_err_bin[i] = L2_error
    Linf_err_bin[i] = Linf_error  
    time_bin[i] = time_init
  end

  return L2_err_bin, Linf_err_bin, time_bin
end

# # STUDY
# Run convergence cases
function solve_dispersion(ω, H0)
    f(k) = ω^2 - 9.81 * k * tanh(k * H0)
    k_initial = ω^2 / 9.81       # Initial guess for k - Deep water approximation plus small offset
    k_solution = find_zero(f, k_initial, Order2())
    return k_solution
end

# -------------------------------------------------------------------------------
# PARAMETERS:

# 1 - Wave parameters
water_depth = 1.1
# water_depth = 10 
ω_check = 10
k_check =solve_dispersion(ω_check, water_depth) * water_depth
# ω_check = ω_data[1]; k_check = k_data[1]
λ_check = 2π/k_check * water_depth

# 2 - Beam parameters
hp = 0.038;  beam_length = 10
# hp = 1.5;  beam_length = 10
# hp = 0.038;beam_length = 0.1 * λ_check                    #2
# hp = 0.038;beam_length = round(0.1 * λ_check, digits=4)   #2
ratio_beam2wave = beam_length/λ_check
ratio_length2thk = beam_length/hp

# 3 - Mesh parameters
# a - element order
r2 = 2; r3 = 3; r4 = 4

# b - nx data
nx_element = [16, 24, 32, 48, 64, 128]
# nx_element = [4, 6, 8, 10, 12, 15]        # 2
# nx_element = [4, 8 ,16, 20, 24]           # 3
mesh_size_var = [beam_length/nx_elem for nx_elem in nx_element] 
# mesh_size_var = [round(beam_length/nx_elem, digits=10) for nx_elem in nx_element]
@show mesh_size_var
n_mesh_size = length(mesh_size_var)

# c - ny data 
ny_data = [2, 3, 4, 6, 8, 16]
# ny_data = [72, 108, 144, 180, 216, 270]     # 2
# ny_data = [12, 24, 48, 60, 72]     # 3

# d - Reference data
r_ref  = 4
nx_ref = 256
# nx_ref = 32        # 2
mesh_ref = beam_length/nx_ref
# mesh_ref = round(beam_length/nx_ref, digits=4)
ny_ref = 32     # Squared custom
# ny_ref = 96     # Squared custom - 2

# 2. Running the reference - finest mesh (f - fine reference)
case_ref = Mooring_case_params(name="Convergence-Study-fine", hₚ = hp, H=water_depth, Lb=beam_length, cλ=4,
           order=r_ref, h=mesh_ref, ny=ny_ref, k_opt=k_check, ω_opt=ω_check, vtk_output=true) 
Γb_f, ηₕ_f, U_Γη_f, V_Γ = run_mooring_case_freq_domain(case_ref)

# For ω = 10 rad/s - slender - ratio ~ 200
h3_rate = mesh_size_var.^(r2+1)*0.55      # [10m vs 20m] 0.20 or 0.30
h4_rate = mesh_size_var.^(r3+1)*0.15
h5_rate = mesh_size_var.^(r4+1)*0.05

# For ω = 10 rad/s - non slender - ratio ~ 6
# h3_rate = mesh_size_var.^(r2+1)*0.026     # [10m vs 20m] 0.20 or 0.30
# h4_rate = mesh_size_var.^(r3+1)*0.007
# h5_rate = mesh_size_var.^(r4+1)*0.0025

# For ω = 10 rad/s - non slender - ratio ~ 1.6
# h3_rate = mesh_size_var.^(r2+1)*2    # [10m vs 20m] 0.20 or 0.30
# h4_rate = mesh_size_var.^(r3+1)*4e2 
# h5_rate = mesh_size_var.^(r4+1)*2e4 

# Make sure that the coordinate of structure node boundaries are located correctly through the sanity check
L2_10m_r2, Linf_10m_r2, ~ = run_convergence_L2(r2) 
L2_10m_r3, Linf_10m_r3, ~ = run_convergence_L2(r3) 
L2_10m_r4, Linf_10m_r4, ~ = run_convergence_L2(r4)

for i=1:1
default(
  tickfontsize=13,
  legendfontsize=13,
  labelfontsize=16,
  titlefontsize=16)

color_scale = palette(:rainbow,5)

p1 = plot(nx_element, real.(L2_10m_r2), label="L2 error, r = $(r2)", lw=2, marker=(:square,4), color=color_scale[1], 
    yscale=:log10, xscale=:log10)
plot!(nx_element, real.(L2_10m_r3), label="L2 error, r = $(r3)", lw=2, marker=(:o,4),  color=color_scale[2], 
    yscale=:log10, xscale=:log10)
plot!(nx_element, real.(L2_10m_r4), label="L2 error, r = $(r4)", lw=2, marker=(:utriangle,4),  color=color_scale[3], 
    yscale=:log10, xscale=:log10)
  
# plot!(1:n_mesh_size, h5x_rate, label="h ~5 reference", linestyle=:dash, color=:black)
plot!(nx_element, h3_rate, label="h³ reference", linestyle=:dash, color=:black)
plot!(nx_element, h4_rate, label="h⁴ reference", linestyle=:dashdot, color=:black)
plot!(nx_element, h5_rate, label="h⁵ reference", linestyle=:dot, color=:black)


plot!(p1, size = (700,600), grid=true, gridlinewidth=1, gridalpha=0.2)
# plot!(p1, legend=:topright)
plot!(p1, legend=:bottomleft)
yticks!(10.0 .^ (-15:0))
ylabel!("Error")
xticks!(nx_element, string.((nx_element)))
# xlabel!("Characteristic mesh size, h")
xlabel!("Number of beam elements, " * L"n_x = L_b/h")
title!("ω = $(round(ω_check, digits=3)) rad/s - $(beam_length)m beam
hₚ = $(hp), ny = $(ny_data[1:n_mesh_size])")
# , ny = $(ny_elem)
# title!("  ")
display(p1)
fig_title1 = "L2_conv_10m"
savefig(p1, plotsdir("4-1-Convergence-Study",fig_title1))

end


