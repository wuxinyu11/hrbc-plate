
using Revise, YAML, ApproxOperator, CairoMakie, JLD

ndiv = 12

config = YAML.load_file("./yml/TADAS dampers-gauss-penalty.yml")
elements, nodes = importmsh("./msh/TADAS dampers.msh", config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 3.5*120/ ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)


set∇²₂𝝭!(elements["Ω"])
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
       Operator(:∫vgdΓ,coefficient...,:α=>1e3),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫∇𝑛vθdΓ,coefficient...,:α=>1e3),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]


k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)

ops[3](elements["Γᵍ"],k,f)
ops[5](elements["Γᵍ"],k,f)

ops[4](elements["Γᵗ"],f)

d = k\f


push!(nodes,:d=>d)
@save compress=true "jld/TADAS_gauss_penalty_"*string(ndiv)*".jld" d