
using XLSX, YAML, ApproxOperator, LinearAlgebra

config = YAML.load_file("./yml/rectangular_rkgsi_nitsche.yml")

ndiv = 20
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)

nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

s = 3.1 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

sp = RegularGrid(nodes[:x], nodes[:y], nodes[:z], n = 2, γ = 5)
sp(elements["Ω"],elements["Γ"],elements["Γₚ₁"],elements["Γₚ₂"],elements["Γₚ₃"],elements["Γₚ₄"])
set_memory_𝝭!(elements["Ωˢ"])

set∇𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ωˢ"],elements["Ω"])
set∇³𝝭!(elements["Γ"])
set∇²𝝭!(elements["Γₚ₁"])
set∇²𝝭!(elements["Γₚ₂"])
set∇²𝝭!(elements["Γₚ₃"])
set∇²𝝭!(elements["Γₚ₄"])

D = 1.0
ν = 0.3

coefficient = (:D=>1.0,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫VgdΓ,:α=>0.0,coefficient...),
       Operator(:∫MₙₙθdΓ,:α=>0.0,coefficient...),
       Operator(:ΔMₙₛg,:α=>0.0,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

# ops[1](elements["Ωˢ"],k)
ops[2](elements["Γ"],k,f)
# ops[3](elements["Γ"],k,f)
# ops[4](elements["Γₚ₁"],k,f)
# ops[4](elements["Γₚ₂"],k,f)
# ops[4](elements["Γₚ₃"],k,f)
# ops[4](elements["Γₚ₄"],k,f)

F = eigen(k)
F.values[1]
# sign(F.values[1])*log10(abs(F.values[1]))