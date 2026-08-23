include("..\\src\\1-2-research-mooring-zhao.jl")

import .Mooring_Case: Mooring_case_params, run_mooring_case_freq_domain, collect_results, run_Mooring

using Plots
using Parameters
using DrWatson
using JLD2
using DataFrames
using CSV
using Printf
using Statistics

# -------------------------------------------------------------
#   1  REFERENCES DATA  ---------------------------------------
# -------------------------------------------------------------
Zhao_5a_38 = DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5a_tp=0.038.csv"); header=false))
Zhao_5a_75 = DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5a_tp=0.075.csv"); header=false))
Zhao_5a_100 = DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5a_tp=0.1.csv"); header=false))
Zhao_5b_38 = DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5b_tp=0.038.csv"); header=false))
Zhao_5b_75 = DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5b_tp=0.075.csv"); header=false))
Zhao_5b_100 = DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5b_tp=0.1.csv"); header=false))
Zhao_5c_38= DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5c_tp=0.038.csv"); header=false))
Zhao_5c_75= DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5c_tp=0.075.csv"); header=false))
Zhao_5c_100 = DataFrame(CSV.File(datadir("Ref_data/Zhao","Zhao_fig5c_tp=0.1.csv"); header=false))
sort!(Zhao_5b_75, :Column1)
sort!(Zhao_5b_100, :Column1)

Utsunomiya_5a_38 = DataFrame(CSV.File(datadir("Ref_data/Utsunomiya","Utsunomiya_BF35_5a.csv"); header=false))
Utsunomiya_5b_38 = DataFrame(CSV.File(datadir("Ref_data/Utsunomiya","Utsunomiya_BF06_5b.csv"); header=false))
Utsunomiya_5c_38 = DataFrame(CSV.File(datadir("Ref_data/Utsunomiya","Utsunomiya_BT07_5c.csv"); header=false))

Karperaki_5a_38 = DataFrame(CSV.File(datadir("Ref_data/Karperaki","Karperaki_w1_38.csv"); header=false))
Karperaki_5b_38 = DataFrame(CSV.File(datadir("Ref_data/Karperaki","Karperaki_w2_38.csv"); header=false))
Karperaki_5c_38 = DataFrame(CSV.File(datadir("Ref_data/Karperaki","Karperaki_w3_38.csv"); header=false))

# -------------------------------------------------------------
#   2  RUN CASE  ----------------------------------------------
# -------------------------------------------------------------
k_data = [0.8044 2.2216 9.0434]     # [Mohapatra / C. Zhao] - dimensionless wave number = k[1/m] * H [m]
T_data = [2.8750 1.4290 0.7000]   
ω_data = 2π ./ T_data
h_mesh = 0.2                            # mesh size
path=datadir("4-2-Zhao-free-floating-comparison")

time_init = @elapsed begin
case =  Mooring_case_params(
  name="Warm-up",
  h=1, # mesh size
  order=2,k_opt=k_data[3], ω_opt=ω_data[3], ks=0, cλ=2)
produce_or_load(path,case,run_Mooring;digits=8)
end
# warmup 267.79 seconds ~ 4.46 minutes 

# No Mooring Comparison 
# comparison per plate thickness - 3 frequencies
#   a, b, c - wave data
#   1, 2, 3 - plate thickness data
case = Mooring_case_params(name="5a-1", hₚ=0.038, ks=0, k_opt=k_data[1], ω_opt=ω_data[1], h=h_mesh, vtk_output=false, cλ=4)
t1 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
case = Mooring_case_params(name="5b-1", hₚ=0.038, ks=0, k_opt=k_data[2], ω_opt=ω_data[2], h=h_mesh, vtk_output=false, cλ=4)
t2 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
case = Mooring_case_params(name="5c-1", hₚ=0.038, ks=0, k_opt=k_data[3], ω_opt=ω_data[3], h=h_mesh, vtk_output=false, cλ=4)
t3 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)

case = Mooring_case_params(name="5a-2", hₚ=0.075, ks=0, k_opt=k_data[1], ω_opt=ω_data[1], h=h_mesh, vtk_output=false, cλ=4)
t4 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
case = Mooring_case_params(name="5b-2", hₚ=0.075, ks=0, k_opt=k_data[2], ω_opt=ω_data[2], h=h_mesh, vtk_output=false, cλ=4)
t5 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
case = Mooring_case_params(name="5c-2", hₚ=0.075, ks=0, k_opt=k_data[3], ω_opt=ω_data[3], h=h_mesh, vtk_output=true, cλ=4)
t6 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)

case = Mooring_case_params(name="5a-3", hₚ=0.100, ks=0, k_opt=k_data[1], ω_opt=ω_data[1], h=h_mesh, vtk_output=false, cλ=4)
t7 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
case = Mooring_case_params(name="5b-3", hₚ=0.100, ks=0, k_opt=k_data[2], ω_opt=ω_data[2], h=h_mesh, vtk_output=false, cλ=4)
t8 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
case = Mooring_case_params(name="5c-3", hₚ=0.100, ks=0, k_opt=k_data[3], ω_opt=ω_data[3], h=h_mesh, vtk_output=false, cλ=4)
t9 = @elapsed data, file = produce_or_load(path,case,run_Mooring;digits=8)
# time = [0.6750128, 0.8657694, 0.3319197, 0.5258941, 0.9416283, 0.4357963, 0.5237433, 0.3708801, 0.8363453]
# mean(time) = 0.6118877

# -------------------------------------------------------------
#   3  OUTPUT DATA COLLECTION  --------------------------------
# -------------------------------------------------------------
result = collect_results(path)
df = DataFrame(result)

res1 = filter(row -> row.name == "5a-1" && row.h == h_mesh, df); xs1 = res1.x; η_xs1 = res1.η
res2 = filter(row -> row.name == "5b-1" && row.h == h_mesh, df); xs2 = res2.x; η_xs2 = res2.η
res3 = filter(row -> row.name == "5c-1" && row.h == h_mesh, df); xs3 = res3.x; η_xs3 = res3.η

res4 = filter(row -> row.name == "5a-2" && row.h == h_mesh, df); xs4 = res4.x; η_xs4 = res4.η
res5 = filter(row -> row.name == "5b-2" && row.h == h_mesh, df); xs5 = res5.x; η_xs5 = res5.η
res6 = filter(row -> row.name == "5c-2" && row.h == h_mesh, df); xs6 = res6.x; η_xs6 = res6.η

res7 = filter(row -> row.name == "5a-3" && row.h == h_mesh, df); xs7 = res7.x; η_xs7 = res7.η
res8 = filter(row -> row.name == "5b-3" && row.h == h_mesh, df); xs8 = res8.x; η_xs8 = res8.η
res9 = filter(row -> row.name == "5c-3" && row.h == h_mesh, df); xs9 = res9.x; η_xs9 = res9.η

x5a = [xs1, xs2, xs3,Zhao_5a_38.Column1,Zhao_5b_38.Column1,Zhao_5c_38.Column1]
η5a = [η_xs1, η_xs2, η_xs3,Zhao_5a_38.Column2,Zhao_5b_38.Column2,Zhao_5c_38.Column2]
x5b = [xs4, xs5, xs6,Zhao_5a_75.Column1,Zhao_5b_75.Column1,Zhao_5c_75.Column1]
η5b = [η_xs4, η_xs5, η_xs6,Zhao_5a_75.Column2,Zhao_5b_75.Column2,Zhao_5c_75.Column2]
x5c = [xs7, xs8, xs9,Zhao_5a_100.Column1,Zhao_5b_100.Column1,Zhao_5c_100.Column1]
η5c = [η_xs7, η_xs8, η_xs9,Zhao_5a_100.Column2,Zhao_5b_100.Column2,Zhao_5c_100.Column2]
hp  = [0.038, 0.075, 0.100]

xKarperaki = [Karperaki_5a_38.Column1, Karperaki_5b_38.Column1, Karperaki_5c_38.Column1]
ηKarperaki = [Karperaki_5a_38.Column2, Karperaki_5b_38.Column2, Karperaki_5c_38.Column2]

xExp = [Utsunomiya_5a_38.Column1, Utsunomiya_5b_38.Column1, Utsunomiya_5c_38.Column1]
ηExp = [Utsunomiya_5a_38.Column2, Utsunomiya_5b_38.Column2, Utsunomiya_5c_38.Column2]
naming = ["ω1", "ω2", "ω3"]

# -------------------------------------------------------------
#   4  PLOTTING  ----------------------------------------------
# -------------------------------------------------------------
default(
    tickfontsize=13,
    legendfontsize=13,
    labelfontsize=16,
    titlefontsize=16)

# 1) Plot plate thickness 38mm
for i = 1:1
  plot(legend=:best, size = (700,600),grid=true, gridlinewidth=1, gridalpha=0.2)
  for i = 1:3
    plot!(x5a[i], η5a[i],xlims=(0,10),line=(2,:red),
          label=false) #"Monolithic C/DG")   
    plot!(x5a[i+3], η5a[i+3], line=(2,:dash,:green), marker=(:square,4,:green), 
          label=false) #"Zhao et al.")
    plot!(xKarperaki[i], ηKarperaki[i], line=(3,:dot,:blue), 
          label=false) #"Karperaki et al.") 
    # uncomment this for experimental data comparison
    # plot!(xExp[i], ηExp[i], marker=(:circle,5), line=:false,
            # label=false)
  end
  plot!([], [], line=(2,:red), label="Monolithic C/DG")
  plot!([], [], line=(2,:green), marker=(:square,4,:green), label="Zhao et al.")
  plot!([], [], line=(2,:dot,:blue), label="Karperaki et al.")
  # plot!([], [], marker=(:circle,5), line=:false, label="Experimental Utsunomiya hₚ = 0.038 m")
  annotate!([
    (3.0, 1.05, text("ω = $(@sprintf("%.4g", ω_data[1])) rad/sec", 14)),
    (3.0, 0.70, text("ω = $(@sprintf("%.4g", ω_data[2])) rad/sec", 14)),
    (3.0, 0.25, text("ω = $(@sprintf("%.4g", ω_data[3])) rad/sec", 14))])
  xlabel!("x [m]")
  ylabel!("|η|/κ₀ [-]")
  ylims!(0,1.5)
  xlims!(0,10)
  title!("hₚ = 0.038m")
  display(current())
  savefig(plotsdir("4-2-Comparison-Zhao-Mindlin","4-2 Zhao-comparison-38.png"))
end

# 2) Plot plate thickness 75mm
for i = 1:1
  plot(legend=:best, size = (700,600),grid=true, gridlinewidth=1, gridalpha=0.2)
  for i = 1:3
    plot!(x5b[i], η5b[i],xlims=(0,10),line=(2,:red),
          label=false) #"Monolithic C/DG"   
    plot!(x5b[i+3], η5b[i+3], line=(2,:dash,:green), marker=(:square,4,:green), 
          label=false) #"Zhao et al."
  end
  plot!([], [], line=(2,:red),         label="Monolithic C/DG")
  plot!([], [], line=(2,:green), marker=(:square,4,:green), label="Zhao et al.")
  annotate!([
    (3.0, 1.00, text("ω = $(@sprintf("%.4g", ω_data[1])) rad/sec", 14)),
    (3.0, 0.55, text("ω = $(@sprintf("%.4g", ω_data[2])) rad/sec", 14)),
    (3.0, 0.20, text("ω = $(@sprintf("%.4g", ω_data[3])) rad/sec", 14))])
  xlabel!("x [m]");  ylabel!("|η|/κ₀ [-]")
  ylims!(0,1.5);     xlims!(0,10)
  title!("hₚ = 0.075m")
  display(current())
  savefig(plotsdir("4-2-Comparison-Zhao-Mindlin","4-2 Zhao-comparison-75.png"))
end

# 3) Plot plate thickness 100mm
for i = 1:1
  plot(legend=:best, size = (700,600),grid=true, gridlinewidth=1, gridalpha=0.2)
  for i = 1:3
    plot!(x5c[i], η5c[i],xlims=(0,10),line=(2,:red),
          label=false) #"Monolithic C/DG"   
    plot!(x5c[i+3], η5c[i+3], line=(2,:dash, :green), marker=(:square,4,:green), 
          label=false) #"Zhao et al."
  end
  plot!([], [], line=(2,:red),         label="Monolithic C/DG")
  plot!([], [], line=(2,:green), marker=(:square,4,:green), label="Zhao et al.")
  annotate!([
    (3.0, 1.00, text("ω = $(@sprintf("%.4g", ω_data[1])) rad/sec", 14)),
    (3.0, 0.55, text("ω = $(@sprintf("%.4g", ω_data[2])) rad/sec", 14)),
    (7.0, 0.55, text("ω = $(@sprintf("%.4g", ω_data[3])) rad/sec", 14))])
  quiver!([7.0], [0.5], quiver=([0], [-0.4]), color=:black)
  quiver!([3.0], [0.5], quiver=([0], [-0.25]), color=:black)  
  xlabel!("x [m]");  ylabel!("|η|/κ₀ [-]")
  ylims!(0,1.5);     xlims!(0,10)
  title!("hₚ = 0.100m")
  display(current())
  savefig(plotsdir("4-2-Comparison-Zhao-Mindlin","4-2 Zhao-comparison-100.png"))
end

time = [t1, t2, t3, t4, t5, t6, t7, t8, t9]
@show time_init
@show time
@show mean(time)