
using Revise, YAML, ApproxOperator, CairoMakie, JLD

ndiv = 8

config = YAML.load_file("./yml/TADAS dampers-rkgsi-penalty-cubic.yml")
elements, nodes = importmsh("./msh/ADAS damper.msh", config)



# nₚ = length(nodes)
# nₑ = length(elements["Ω"])
nₚ = getnₚ(elements["Ω"])

s = 3.5*320/ ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)



set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])

set∇₂𝝭!(elements["Γᵍ"])
set∇₂𝝭!(elements["Γᵗ"])






E = 2E11;
h = 10;
ν = 0.3;
D = E*h^3/(12*(1-ν^2));
P = 1E5;
w(x,y) = 0.0
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γᵗ"],:V=>(x,y,z)->-P)

coefficient = (:D=>D,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫vgdΓ,coefficient...,:α=>1e23*E),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫∇𝑛vθdΓ,coefficient...,:α=>1e23*E),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)

ops[3](elements["Γᵍ"],k,f)
ops[5](elements["Γᵍ"],k,f)

ops[4](elements["Γᵗ"],f)

d = k\f


# push!(nodes,:d=>d)
# @save compress=true "png/ADAS_penalty.jld" d