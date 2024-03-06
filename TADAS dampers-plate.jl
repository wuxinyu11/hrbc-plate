
using Revise, YAML, ApproxOperator, CairoMakie

ndiv = 12

config = YAML.load_file("./yml/TADAS dampers-cubic.yml")
elements, nodes = importmsh("./msh/TADAS dampers.msh", config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 3.5*12/ ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γᵗ"],:𝝭)


elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]

set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set∇∇̃²𝝭!(elements["Γᵗ"],elements["Ω∩Γᵗ"])

set∇₂𝝭!(elements["Γᵍ"])
set∇₂𝝭!(elements["Γᵗ"])

set∇∇̄²𝝭!(elements["Γᵍ"],Γᵍ=elements["Γᵍ"],Γᶿ=elements["Γᵍ"],Γᴾ=elements["Γₚ"])


E = 2E11;
h = 3.0;
ν = 0.3;
D = E*h^3/(12*(1-ν^2));
P = 100.0;
w(x,y) = 0.0
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γᵗ"],:V=>(x,y,z)->-P)

coefficient = (:D=>D,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫ṼgdΓ,coefficient...),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫M̃ₙₙθdΓ,coefficient...),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔM̃ₙₛg,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]


k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)

ops[3](elements["Γᵍ"],k,f)
ops[5](elements["Γᵍ"],k,f)

ops[4](elements["Γᵗ"],f)

d = k\f


push!(nodes,:d=>d)

