using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample

include("impor_damper.jl")
ndiv = 100
# elements, nodes = import_damper_fem("msh/honeycomb_damper_"*string(ndiv)*".msh");
elements, nodes = import_damper_fem("msh/slit damper_"*string(ndiv)*".msh");

nₚ = length(nodes)

E = 2E11;ν = 0.3;h = 10
D = E*h^3/(12*(1-ν^2));
P = 1E5;

eval(prescribeForFem)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢgᵢds}(:α=>1e5*E),
]
k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)


ops[1](elements["Ω"],k)
ops[2](elements["Γᵗ"],f)
ops[3](elements["Γᵍ"],k,f)

d = k\f

d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]
print(d₁[1])
print(d₂[1])