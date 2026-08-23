module Mooring_Case 

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Plots
using Parameters
using DrWatson
using JLD2
using DataFrames
using CSV

export Mooring_case_params, run_mooring_case_freq_domain, collect_results, run_Mooring

# -------------------------------------------------------------
#   1  --------------------------------------------------------
# -------------------------------------------------------------

function collect_results(path::AbstractString)
  files = filter(f->endswith(f, ".jld2"), readdir(path))
  results = Vector{Dict{Symbol,Any}}()
  for f in sort(files)
    file = joinpath(path, f)
    obj = JLD2.load(file)      # loads dictionary saved by DrWatson.save
    dict_sym = Dict(Symbol(k) => v for (k,v) in obj)   # convert keys to Symbols
    dict_sym[:__file__] = f
    push!(results, dict_sym)
  end
  return results
end

# -------------------------------------------------------------
#   2  --------------------------------------------------------
# -------------------------------------------------------------

@with_kw struct Mooring_case_params
  name::String = "MooringCase"
  # nx::Int = 20
  ny::Int = 5
  order::Int  = 4          
  kr::Float64 = 0             # rotational springs
  ks::Float64 = 0             # translational springs
  ks2:: Float64 = 0           # internal springs
  vtk_output::Bool = true
  hₚ::Float64 = 0.1            # [C. Zhao] plate thickness
  H::Float64 = 10.0           # [Michele] water depth
  k_opt::Float64 = 1.0 
  ω_opt::Float64 = 1.0
  v_pto::Float64 = 0.0        # PTO damping 
  T₀::Float64 = 0.0           # Pre-tensioned [Agarwal]
  h::Float64 = 0.2            # grid size
  cλ::Float64 = 2             # wave damping domain length coefficient
  Lb:: Float64 = 40           # [Michele] plate length   
  ε::Float64 = π^2/12         # coefficient for shear rigidity, used to switch between EB and Timoshenko beam theories (EB, ε → ∞)
  J::Float64 = hₚ^3 / 12       # Angular moment of inertia for torsion, for EB beam ~ 0
end

# -------------------------------------------------------------
#   3  --------------------------------------------------------
# -------------------------------------------------------------

function run_Mooring(case::Mooring_case_params)
  case_name = savename(case)
  println("-------------")
  println("Case: ",case_name)
  x, η = run_mooring_case_freq_domain(case)
  case_name_suffix = savename(case,"jld2";digits=8)
  file = datadir("Data bin Michele Comparison Data", case_name_suffix)
  prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
  push!(data,"x"=>x, "η"=>η)  
  save(file,data)
  return data
end

# -------------------------------------------------------------
#   4  --------------------------------------------------------
# -------------------------------------------------------------

### The function is disabled for debuging, [FIXED] later
function run_mooring_case_freq_domain(params::Mooring_case_params) 
  # Unpack input parameters
  @unpack name, ny, order, vtk_output, kr, ks, ks2, hₚ, H, k_opt, ω_opt, v_pto, T₀, h, cλ, Lb, ε, J = params        # for running

  # Fixed parameters
  # Lb = 40                     # [Michele] plate length 
  E  = 103 * 1e6                # [C. Zhao] young modulus [Pa]
  νₚ = 0.3                       # [C. Zhao] poisson ratio
  G  =  E / (2*(1+νₚ))  
  EI = 6.9 * 10^4               # [Michele] flexural rigidity
  C = ε * G * hₚ                 # [C. Zhao] shear rigidity - c⋅GA, → ∞ for EB beam theory
  β = 0.2
  α = 0.249
   
  # Physics
  mass = 44.0                  # [Michele] plate mass (kg/m)
  g = 9.81
  ρ = 1000                     # [C. Zhao / Michele] water density
  ρ_ratio = 4.5455             # ρ/ρ₀ (water/plate density)
  ρ₀ = mass/hₚ                  # [Michele] plate density
  d₀ = mass / ρ
  Dᵨ = EI/ρ
  ksᵨ = ks / ρ                 # Mooring stiffness
  ks2ᵨ = ks2 / ρ 
  krᵨ = kr / ρ
  vᵨ = v_pto/ρ                 # PTO damping coefficient 
  Tᵨ = T₀/ρ

  # wave properties
  k = k_opt
  λ = 2π / k
  ω = ω_opt
  η₀ = 1.0                  # [Michele] wave amplitude
  ηᵢₙ(x) = η₀*exp(im*k*x[1])
  ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*x[2]) / sinh(k*H))*exp(im*k*x[1])
  vᵢₙ(x) = (η₀*ω)*(cosh(k*x[2]) / sinh(k*H))*exp(im*k*x[1])
  vzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])

  # Domain size
  if cλ*λ < Lb
    Ld = Lb
  else
    Ld = ceil(cλ*λ)
  end
  LΩ = 2Ld + 2Lb
  x₀ = 0.0
  xdᵢₙ = x₀ + Ld
  xb₀ = xdᵢₙ + Lb/2
  xb₁ = xb₀ + Lb
  xdₒᵤₜ = LΩ - Ld
  nSegments = 3
  nJoints = nSegments - 1
  segment_length = Lb / nSegments
  # digits and base are used to ensure that the joint location is approximated to the nearest nodes location
  xbⱼ₁ = round(xb₀ + segment_length, digits = 1, base = Int(1/h))
  xbⱼ₂ = round(xb₀ + 2*segment_length, digits = 1, base = Int(1/h))
  
  # Numerics constants
  nx_total = Int(ceil(LΩ/h))   # recompute nx_total based on h
  γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  αₕ = -im*ω/g * (1-βₕ)/βₕ
  # @show Ld
  # @show LΩ
  # @show nx_total
  # @show h

  # Damping
  μ₀ = 2.5
  μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1])/Ld))
  μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒᵤₜ)/Ld))
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
  ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
  ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzᵢₙ(x)

  # Fluid model
  domain = (x₀, LΩ, 0.0, H)
  partition = (nx_total,ny)
  function f_y(x)
    if x == H
        return H
    end
    i = x / (H/ny)
    return H-H/(2^i)
    # return H-H/(2.5^i)
  end
  map(x) = VectorValue(x[1], f_y(x[2]))
  # map(x) = VectorValue(x[1], x[2]) 
  𝒯_Ω = CartesianDiscreteModel(domain,partition,map=map)

  # Labelling
  labels_Ω = get_face_labeling(𝒯_Ω)
  add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
  add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
  add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
  add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)

  # Triangulations
  Ω = Interior(𝒯_Ω)
  Γ = Boundary(𝒯_Ω,tags="surface")
  Γin = Boundary(𝒯_Ω,tags="inlet")

  # Auxiliar functions
  function is_beam1(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xb₀ <= x[1] <= xbⱼ₁ ) * ( x[2] ≈ H)
  end
  function is_beam2(xs) # Check if an element is inside the beam2
    n = length(xs)
    x = (1/n)*sum(xs)
    (xbⱼ₁ <= x[1] <= xbⱼ₂ ) * ( x[2] ≈ H)
  end
  function is_beam3(xs) # Check if an element is inside the beam3
    n = length(xs)
    x = (1/n)*sum(xs)
    (xbⱼ₂ <= x[1] <= xb₁ ) * ( x[2] ≈ H)
  end
  function is_damping1(xs) # Check if an element is inside the damping zone 1
    n = length(xs)
    x = (1/n)*sum(xs)
    (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ H)
  end
  function is_damping2(xs) # Check if an element is inside the damping zone 2
    n = length(xs)
    x = (1/n)*sum(xs)
    (xdₒᵤₜ <= x[1] ) * ( x[2] ≈ H)
  end
  function is_joint1(xs) # Check if an element is a joint
    is_on_xbⱼ₁ = [x[1]≈xbⱼ₁ && x[2]≈H for x in xs] 
    element_on_xbⱼ₁ = minimum(is_on_xbⱼ₁) 
    element_on_xbⱼ₁
  end
  function is_joint2(xs) # Check if an element is a joint
    is_on_xbⱼ₂ = [x[1]≈xbⱼ₂ && x[2]≈H for x in xs] 
    element_on_xbⱼ₂ = minimum(is_on_xbⱼ₂) 
    element_on_xbⱼ₂
  end
  function is_beam(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xb₀ <= x[1] <= xb₁ ) * ( x[2] ≈ H)
  end  

  # Beam triangulations
  xΓ = get_cell_coordinates(Γ)
  Γb1_to_Γ_mask = lazy_map(is_beam1,xΓ)
  Γb2_to_Γ_mask = lazy_map(is_beam2,xΓ) 
  Γb3_to_Γ_mask = lazy_map(is_beam3,xΓ) 
  Γd1_to_Γ_mask = lazy_map(is_damping1,xΓ)
  Γd2_to_Γ_mask = lazy_map(is_damping2,xΓ)
  Γb1_to_Γ = findall(Γb1_to_Γ_mask)
  Γb2_to_Γ = findall(Γb2_to_Γ_mask)
  Γb3_to_Γ = findall(Γb3_to_Γ_mask)
  Γd1_to_Γ = findall(Γd1_to_Γ_mask)
  Γd2_to_Γ = findall(Γd2_to_Γ_mask)
  Γf_to_Γ = findall(!,Γb1_to_Γ_mask .| Γb2_to_Γ_mask .| Γb3_to_Γ_mask .| Γd1_to_Γ_mask .| Γd2_to_Γ_mask)
  Γη_to_Γ = findall(Γb1_to_Γ_mask .| Γb2_to_Γ_mask .| Γb3_to_Γ_mask )
  Γκ_to_Γ = findall(!,Γb1_to_Γ_mask .| Γb2_to_Γ_mask .| Γb3_to_Γ_mask )
  Γb1 = Triangulation(Γ,Γb1_to_Γ)
  Γb2 = Triangulation(Γ,Γb2_to_Γ)
  Γb3 = Triangulation(Γ,Γb3_to_Γ)
  Γd1 = Triangulation(Γ,Γd1_to_Γ)
  Γd2 = Triangulation(Γ,Γd2_to_Γ)
  Γfs = Triangulation(Γ,Γf_to_Γ)
  Γη = Triangulation(Γ,Γη_to_Γ)
  Γκ = Triangulation(Γ,Γκ_to_Γ)
  Λb1 = Skeleton(Γb1)
  Λb2 = Skeleton(Γb2)
  Λb3 = Skeleton(Γb3)

  Γ_mask_in_Ω_dim_0 = get_face_mask(labels_Ω,"surface",0)
  grid_dim_0_Γ = GridPortion(Grid(ReferenceFE{0},𝒯_Ω),Γ_mask_in_Ω_dim_0)
  xΓ_dim_0 = get_cell_coordinates(grid_dim_0_Γ)

# --------------------------------------------------------
# --------------------------------------------------------
# 1 Beam joints
  Λj1_to_Γ_mask = lazy_map(is_joint1,xΓ_dim_0)
  Λj2_to_Γ_mask = lazy_map(is_joint2,xΓ_dim_0)
  Λj1_to_Γ = findall(Λj1_to_Γ_mask)
  Λj2_to_Γ = findall(Λj2_to_Γ_mask)  
  Λj1 = Skeleton(Γ,Λj1_to_Γ_mask)
  Λj2 = Skeleton(Γ,Λj2_to_Γ_mask)

  println("Sanity check for the joints")
  @show num_dims(Λj1)
  @show num_cells(Λj1)
  println("Node Λj1 = ",Λj1_to_Γ)
  println("located at ", get_cell_coordinates(Λj1))
  @show num_dims(Λj2)
  @show num_cells(Λj2)
  println("Node Λj2 = ",Λj2_to_Γ)
  println("located at ", get_cell_coordinates(Λj2))

# --------------------------------------------------------
# --------------------------------------------------------
# 2 Beam boundary conditions, left-right ends
  Γb1_cells_coords = get_cell_coordinates(Γb1)
  Γb2_cells_coords = get_cell_coordinates(Γb2)
  Γb3_cells_coords = get_cell_coordinates(Γb3)

  function is_at_left_end(xs)
    return any(x -> x[1] ≈ xb₀ && x[2] ≈ H, xs)
  end
  function is_at_right_end(xs)
    is_on_xbᵣ = [x[1]≈xb₁ && x[2]≈H for x in xs] 
    element_on_xbᵣ = minimum(is_on_xbᵣ) 
    element_on_xbᵣ
  end
  
  Λm1_to_Γb_mask = lazy_map(is_at_left_end, Γb1_cells_coords)
  Λm2_to_Γ_mask = lazy_map(is_at_right_end, xΓ_dim_0)
  Γₗ_to_Γ = findall(Λm1_to_Γb_mask)
  Γᵣ_to_Γ = findall(Λm2_to_Γ_mask)
  Γₗ = Boundary(Γb1, Λm1_to_Γb_mask)
  Γᵣ = Boundary(Γ, Λm2_to_Γ_mask)

  println("Sanity check for the mooring nodes")
  @show num_dims(Γₗ)                      # Expected to be 0 for node boundary
  @show num_cells(Γₗ)                     # Expected to be 0 for node boundary  
  println("Node Γₗ = ",Γₗ_to_Γ)
  println("located at ", get_cell_coordinates(Γₗ))        # Expected to have a coordinate equal to (xb₀, H)
  @show num_dims(Γᵣ)
  @show num_cells(Γᵣ)
  println("Node Γᵣ = ",Γᵣ_to_Γ)
  println("located at ", get_cell_coordinates(Γᵣ))        # Expected to have a coordinate equal to (xb₁, H)
  
  if vtk_output == true
    filename = "data/VTKOutput/Mooring-Monolithic-VTK/"*name
    writevtk(Ω,filename*"_Ω")
    writevtk(Γ,filename*"_Γ")
    writevtk(Γb1,filename*"_Γb1")
    writevtk(Γb2,filename*"_Γb2")
    writevtk(Γb3,filename*"_Γb3")    
    writevtk(Γd1,filename*"_Γd1")
    writevtk(Γd2,filename*"_Γd2")
    writevtk(Γfs,filename*"_Γfs")
    writevtk(Λb1,filename*"_Λb1")
    writevtk(Λb2,filename*"_Λb2")
    writevtk(Λb3,filename*"_Λb3")        
    writevtk(Λj1,filename*"_Λj1")
    writevtk(Λj2,filename*"_Λj2")        
    writevtk(Γₗ,filename*"_ΓbLeft")
    writevtk(Γᵣ,filename*"_ΓbRight")
    end

  # Measures
  degree = 2*order        
  dΩ = Measure(Ω,degree)
  dΓb1 = Measure(Γb1,degree)
  dΓb2 = Measure(Γb2,degree)
  dΓb3 = Measure(Γb3,degree)
  dΓd1 = Measure(Γd1,degree)
  dΓd2 = Measure(Γd2,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΛb1 = Measure(Λb1,degree)
  dΛb2 = Measure(Λb2,degree)
  dΛb3 = Measure(Λb3,degree)
  dΛj1 = Measure(Λj1,degree)
  dΛj2 = Measure(Λj2,degree)
  dΓₗ = Measure(Γₗ,degree)
  dΓᵣ = Measure(Γᵣ,degree)

  # Normals
  nΓₗ = get_normal_vector(Γₗ)
  nΓᵣ = get_normal_vector(Γᵣ)
  nΛb1 = get_normal_vector(Λb1)
  nΛb2 = get_normal_vector(Λb2)
  nΛb3 = get_normal_vector(Λb3)
  nΛj1 = get_normal_vector(Λj1)
  nΛj2 = get_normal_vector(Λj2)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  reffeᵩ= ReferenceFE(lagrangian,VectorValue{1,Float64},order-1)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γηᵩ = TestFESpace(Γη, reffeᵩ, conformity=:H1, vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  U_Γη = TrialFESpace(V_Γη)
  U_Γηᵩ = TrialFESpace(V_Γηᵩ)

  X = MultiFieldFESpace([U_Ω, U_Γκ, U_Γη])
  Y = MultiFieldFESpace([V_Ω, V_Γκ, V_Γη])

  A2 = ρ₀*J*ω^2/ρ + d₀*ρ*Dᵨ*ω^2/C - ρ*g*Dᵨ/C
  A3 = d₀*ρ₀*J*ω^4/C - d₀*ω^2 + g - ρ₀*g*J*ω^2/C
  A4 = im*ρ*Dᵨ*ω/C
  A5 = im*ρ₀*J*ω^3/C - im*ω

  @show Dᵨ
  @show Tᵨ
  @show vᵨ
  @show A2
  @show A3
  @show A4
  @show A5

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  ∇ₓ(ϕ) = ∇(ϕ)⋅VectorValue(1.0,0.0)
  # CASE 1
  a((ϕ,κ,η),(w,u,v)) = (∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1 +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2 +
    ∫(  ( Dᵨ*Δ(v)*Δ(η) - A2*(∇(v)⋅∇(η)) + A3*v*η - A4*(∇(v)⋅∇(ϕ)) + A5*v*ϕ  + Tᵨ*(∇(v)⋅∇(η))) +  im*ω*w*η )dΓb1 +
    ∫(  ( Dᵨ*Δ(v)*Δ(η) - A2*(∇(v)⋅∇(η)) + A3*v*η - A4*(∇(v)⋅∇(ϕ)) + A5*v*ϕ  + Tᵨ*(∇(v)⋅∇(η))) +  im*ω*w*η )dΓb2 +
    ∫(  ( Dᵨ*Δ(v)*Δ(η) - A2*(∇(v)⋅∇(η)) + A3*v*η - A4*(∇(v)⋅∇(ϕ)) + A5*v*ϕ  + Tᵨ*(∇(v)⋅∇(η))) +  im*ω*w*η )dΓb3 +
    ∫(  (Dᵨ) * ( - jump(∇(v)⋅nΛb1) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb1) + γ*( jump(∇(v)⋅nΛb1) * jump(∇(η)⋅nΛb1) ) ) )dΛb1 +
    ∫(  (Dᵨ) * ( - jump(∇(v)⋅nΛb2) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb2) + γ*( jump(∇(v)⋅nΛb2) * jump(∇(η)⋅nΛb2) ) ) )dΛb2 +
    ∫(  (Dᵨ) * ( - jump(∇(v)⋅nΛb3) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb3) + γ*( jump(∇(v)⋅nΛb3) * jump(∇(η)⋅nΛb3) ) ) )dΛb3 +
    ∫(  (mean(v) * ks2ᵨ * mean(η)) + (mean(v) * (-im*ω*vᵨ) * mean(η)) + (-Dᵨ*jump(∇(v)⋅nΛj1)*mean(Δ(η))) + 
          (-Dᵨ*jump(∇(η)⋅nΛj1)*mean(Δ(v)))  - (-Dᵨ*γ*(jump(∇(v)⋅nΛj1) * jump(∇(η)⋅nΛj1)))  )dΛj1  + 
    ∫(  (mean(v) * ks2ᵨ * mean(η)) + (mean(v) * (-im*ω*vᵨ) * mean(η)) + (-Dᵨ*jump(∇(v)⋅nΛj2)*mean(Δ(η))) +
          (-Dᵨ*jump(∇(η)⋅nΛj2)*mean(Δ(v)))  - (-Dᵨ*γ*(jump(∇(v)⋅nΛj2) * jump(∇(η)⋅nΛj2)))  )dΛj2  +
    ∫(  (ksᵨ*(v*η) + krᵨ*(∇(v)⋅∇(η)) ) + (-im*ω*vᵨ* (v*η)) )dΓₗ  + 
    ∫(  (ksᵨ*(v*η) + krᵨ*(∇(v)⋅∇(η)) ) + (-im*ω*vᵨ* (v*η)) )dΓᵣ) 
  l((w,u,v)) =  ∫( w*vᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1

  # Solution
  op = AffineFEOperator(a,l,X,Y)
  (ϕₕ,κₕ,ηₕ) = solve(op)

 if vtk_output == true
    writevtk(Ω,filename * "_Ω_solution.vtu",cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ),"phi_abs"=>abs(ϕₕ)],nsubcells=10)
    writevtk(Γκ,filename * "_Γκ_solution.vtu",cellfields = ["kappa_re" => real(κₕ),"kappa_im" => imag(κₕ),"kappa_abs"=>abs(κₕ)],nsubcells=10)
    writevtk(Γη,filename * "_Γη_solution.vtu",cellfields = ["eta_re" => real(ηₕ),"eta_im" => imag(ηₕ),"eta_abs"=>abs(ηₕ)],nsubcells=10)
  end

# -------------------------------------------------------
  # Postprocess
  # [1] Beam Deflection Solutions
  xy_cp = get_cell_points(get_fe_dof_basis(V_Γη)).cell_phys_point
  x_cp = [[xy_ij[1] for xy_ij in xy_i] for xy_i in xy_cp]
  η_cdv = get_cell_dof_values(ηₕ)

  p = sortperm(x_cp[1])
  x_cp_sorted = [x_i[p] for x_i in x_cp]
  η_cdv_sorted = [η_i[p] for η_i in η_cdv]

  xs = [x_i-xb₀ for x_i in vcat(x_cp_sorted...)]
  η_rel_xs = [abs(η_i)/η₀ for η_i in vcat(η_cdv_sorted...)]

  return (xs, η_rel_xs)
end
end
