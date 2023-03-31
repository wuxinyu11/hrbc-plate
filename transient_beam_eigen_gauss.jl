using Revise, ApproxOperator, YAML, LinearAlgebra, CairoMakie, XLSX

ndiv = 80
𝒑 = "cubic"
config = YAML.load_file("./yml/beam_rkgsi_gauss_"*𝒑*".yml")
elements, nodes = importmsh("./msh/beam_"*string(ndiv)*".msh",config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

push!(getfield(elements["Γ"][1].𝓖[1],:data),:n₁=>(2,[-1.0,1.0]))
set∇²₁𝝭!(elements["Ω"])
set∇³₁𝝭!(elements["Γ"])

ρ = 2500.0
h = 1.0
A = 1.0
L = 10.0
EI = 1.0/6.0*1e6

ops = [
    Operator(:∫κMdx,:EI=>EI),
    Operator(:∫ρhvwdΩ,:ρ=>ρ,:h=>h),
    Operator(:∫wVdΓ),
    Operator(:Vg,:EI=>EI,:α=>1e3),
]
k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],m)
ops[4](elements["Γ"],kα,fα)

vals, vecs = LinearAlgebra.eigen(k+kα,m)  

ω(n) = (n*π)^2*(EI/ρ/A/L^4)^0.5
err1 = vals[3]^0.5/ω(1) - 1.0
err2 = vals[4]^0.5/ω(2) - 1.0
err3 = vals[5]^0.5/ω(3) - 1.0
err4 = vals[6]^0.5/ω(4) - 1.0

index = [10,20,40,80]
XLSX.openxlsx("./xlsx/beam_eigen_"*𝒑*".xlsx", mode="rw") do xf
    row = "B"
    ω₁ = xf[2]
    ω₂ = xf[3]
    ω₃ = xf[4]
    ω₄ = xf[5]

    ind = findfirst(n->n==ndiv,index)+1
    row = row*string(ind)
    ω₁[row] = log10(abs(err1))
    ω₂[row] = log10(abs(err2))
    ω₃[row] = log10(abs(err3))
    ω₄[row] = log10(abs(err4))
end