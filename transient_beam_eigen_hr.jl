using Revise, ApproxOperator, YAML, LinearAlgebra, CairoMakie, XLSX

ndiv = 40
𝒑 = "cubic"
config = YAML.load_file("./yml/beam_rkgsi_hr_"*𝒑*".yml")
elements, nodes = importmsh("./msh/beam_"*string(ndiv)*".msh",config)

set_memory_𝗠!(elements["Ω̃"],:∇̃²,:𝝭)
set_memory_𝝭!(elements["Ω̃"],:𝝭,:∂²𝝭∂x²)
set_memory_𝗠!(elements["Γ"],:𝝭,:∂𝝭∂x,:∇̃²,:∂∇̃²∂ξ)

nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

set∇₁𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set𝝭!(elements["Ω"])
set𝝭!(elements["Ω̃"])
set∇₁𝝭!(elements["Γ"])
set∇∇̃²𝝭!(elements["Γ"],elements["Ω"][[1,nₑ]])
set∇∇̄²𝝭!(elements["Γ"])

e1 = 0.0
e2 = 0.0
e3 = 0.0
ξ = elements["Γ"][1].𝓖[1]
B = ξ[:∂∂²𝝭∂x²∂x]
for (i,xᵢ) in enumerate(elements["Γ"][1].𝓒)
    global e1 += B[i]*xᵢ.x
    global e2 += B[i]*xᵢ.x^2
    global e3 += B[i]*xᵢ.x^3
end
# e3 -= 6*ξ.x

ρ = 2500.0
h = 1.0
A = 1.0
L = 10.0
EI = 1.0/6.0*1e6

ops = [
    Operator(:∫κMdx,:EI=>EI),
    Operator(:∫ρhvwdΩ,:ρ=>ρ,:h=>h),
    Operator(:∫wVdΓ),
    Operator(:Ṽg,:EI=>EI),
]
k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
ops[1](elements["Ω̃"],k)
# ops[2](elements["Ω"],m)
ops[2](elements["Ω̃"],m)
ops[4](elements["Γ"],kα,fα)

vals, vecs = LinearAlgebra.eigen(k+kα,m)  
# vals, vecs = LinearAlgebra.eigen(k,m)  

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

cpn = zeros(nₚ)
kΔx = zeros(nₚ)
for n in 1:nₚ
    cpn[n] = vals[n]^0.5/ω(n)
    kΔx[n] = (n-1)/(nₚ-1)
end
# for n in 2:nₚ
#     cpn[n] = vals[n]^0.5/ω(n-1)
#     kΔx[n] = (n-1)/(nₚ-1)
# end
f = Figure()
ax = Axis(f[1,1])
lines!(kΔx,cpn)
f