using Revise, YAML, ApproxOperator, CairoMakie, JLD

ndiv = 9

config = YAML.load_file("./yml/ADAS dampers-rkgsi-nitsche-cubic.yml")
elements, nodes = importmsh("./msh/ADAS damper.msh", config)
nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 3.5*240/ ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)


set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])

set∇³𝝭!(elements["Γᵍ"])
set∇³𝝭!(elements["Γᵗ"])





E = 2E11;
h = 5;
ν = 0.3;
D = E*h^3/(12*(1-ν^2));
P = 1E5;
w(x,y) = 0.0
set𝒏!(elements["Γᵍ"])
set𝒏!(elements["Γᵗ"])
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γᵗ"],:V=>(x,y,z)->-P)

coefficient = (:D=>D,:ν=>0.3)


coefficient = (:D=>D,:ν=>ν)
ops = [Operator(:∫κ̃ᵢⱼM̃ᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫VgdΓ,coefficient...,:α=>1e5*E),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫MₙₙθdΓ,coefficient...,:α=>1e5*E),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔMₙₛg,coefficient...,:α=>1e1),
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
@save compress=true "jld/ADAS_nitsche.jld" d