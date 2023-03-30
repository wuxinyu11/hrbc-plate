using Revise, ApproxOperator, YAML, LinearAlgebra, CairoMakie, XLSX

ndiv = 10
𝒑 = "cubic"
config = YAML.load_file("./yml/beam_rkgsi_penalty_"*𝒑*".yml")
elements, nodes = importmsh("./msh/beam_"*string(ndiv)*".msh",config)

data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=4,γ=6)
data = Dict([:x=>(2,[5.0]),:y=>(2,[0.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Cubic1D,:□,:QuinticSpline,:Poi1}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝝭!(elements["Ω̃"],:∂²𝝭∂x²)

nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

set∇₁𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γ"])

ρ = 2500.0
h = 1.0
A = 1.0
L = 10.0
EI = 1.0/6.0*1e6

ops = [
    Operator(:∫κMdx,:EI=>EI),
    Operator(:∫ρhvwdΩ,:ρ=>ρ,:h=>h),
    Operator(:∫wVdΓ),
    Operator(:∫vgdΓ,:α=>1e15),
]
k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],m)
ops[4](elements["Γ"],kα,fα)

vals, vecs = LinearAlgebra.eigen(k+kα,m)  

ω(n) = (n*π)^2*(EI/ρ/A/L^4)^0.5
err1 = vals[1]^0.5/ω(1) - 1.0
err2 = vals[2]^0.5/ω(2) - 1.0
err3 = vals[3]^0.5/ω(3) - 1.0
err4 = vals[4]^0.5/ω(4) - 1.0


index = [10,20,40,80]
XLSX.openxlsx("./xlsx/beam_eigen_"*𝒑*".xlsx", mode="rw") do xf
    row = "A"
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