
using XLSX, YAML, ApproxOperator, LinearAlgebra

config = YAML.load_file("./yml/rectangular_rkgsi_nitsche.yml")

ndiv = 10
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 3.1 / ndiv * ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)

elements["Γₚ"] = elements["Γₚ₁"]∪elements["Γₚ₂"]∪elements["Γₚ₃"]∪elements["Γₚ₄"]
elements["Γ"] = elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]

set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇³𝝭!(elements["Γ₁"])
set∇³𝝭!(elements["Γ₂"])
set∇³𝝭!(elements["Γ₃"])
set∇³𝝭!(elements["Γ₄"])
set∇²₂𝝭!(elements["Γₚ₁"])
set∇²₂𝝭!(elements["Γₚ₂"])
set∇²₂𝝭!(elements["Γₚ₃"])
set∇²₂𝝭!(elements["Γₚ₄"])

D = 1.0
ν = 0.3

coefficient = (:D=>1.0,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫VgdΓ,:α=>1e0,coefficient...),
       Operator(:∫MₙₙθdΓ,:α=>1e8,coefficient...),
       Operator(:ΔMₙₛg,:α=>0e5,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)
# ops[2](elements["Γ₁"],k,f)
ops[2](elements["Γ"],k,f)
ops[3](elements["Γ"],k,f)
# ops[4](elements["Γₚ"],k,f)
# ops[4](elements["Γₚ₂"],k,f)
# ops[4](elements["Γₚ₃"],k,f)
# ops[4](elements["Γₚ₄"],k,f)

F = eigen(k)
F.values[1]
# sign(F.values[1])*log10(abs(F.values[1]))