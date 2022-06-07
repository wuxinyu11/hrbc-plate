
using Revise, YAML, ApproxOperator

config = YAML.load_file("./yml/rkgsi.yml")

ndiv = 10
elements, nodes = importmsh("./msh/patchtest_"*string(ndiv)*".msh", config)

nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

type = (SNode,:Cubic2D,:□,:QuinticSpline)
s = 3.5/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
n = 5

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 2,γ = 5)
sp(elements["Ω"],elements["Γ₁"],elements["Γ₂"],elements["Γ₃"],elements["Γ₄"],elements["Γₚ₁"],elements["Γₚ₂"],elements["Γₚ₃"],elements["Γₚ₄"])
set_memory_𝝭!(elements["Ωˢ"])

set∇𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ωˢ"],elements["Ω"])
set∇²𝝭!(elements["Γₚ₁"])
set∇²𝝭!(elements["Γₚ₂"])
set∇²𝝭!(elements["Γₚ₃"])
set∇²𝝭!(elements["Γₚ₄"])
set∇³𝝭!(elements["Γ₁"])
set∇³𝝭!(elements["Γ₂"])
set∇³𝝭!(elements["Γ₃"])
set∇³𝝭!(elements["Γ₄"])

# w(x,y) = (1+2x+3y)^n
# w₁(x,y) = 2n*(1+2x+3y)^abs(n-1)
# w₂(x,y) = 3n*(1+2x+3y)^abs(n-1)
# w₁₁(x,y) = 4n*(n-1)*(1+2x+3y)^abs(n-2)
# w₂₂(x,y) = 9n*(n-1)*(1+2x+3y)^abs(n-2)
# w₁₂(x,y) = 6n*(n-1)*(1+2x+3y)^abs(n-2)
# w₁₁₁(x,y) = 8n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
# w₁₁₂(x,y) = 12n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
# w₁₂₂(x,y) = 18n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
# w₂₂₂(x,y) = 27n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
# w₁₁₁₁(x,y) = 16n*(n-1)*(n-2)*(n-3)*(1+2x+3y)^abs(n-4)
# w₁₁₂₂(x,y) = 36n*(n-1)*(n-2)*(n-3)*(1+2x+3y)^abs(n-4)
# w₂₂₂₂(x,y) = 81n*(n-1)*(n-2)*(n-3)*(1+2x+3y)^abs(n-4)
w(x,y) = - sin(π*x)*sin(π*y)
w₁(x,y) = - π*cos(π*x)*sin(π*y)
w₂(x,y) = - π*sin(π*x)*cos(π*y)
w₁₁(x,y) = π^2*sin(π*x)*sin(π*y)
w₂₂(x,y) = π^2*sin(π*x)*sin(π*y)
w₁₂(x,y) = - π^2*cos(π*x)*cos(π*y)
w₁₁₁(x,y) = π^3*cos(π*x)*sin(π*y)
w₁₁₂(x,y) = π^3*sin(π*x)*cos(π*y)
w₁₂₂(x,y) = π^3*cos(π*x)*sin(π*y)
w₂₂₂(x,y) = π^3*sin(π*x)*cos(π*y)
w₁₁₁₁(x,y) = - π^4*sin(π*x)*sin(π*y)
w₁₁₂₂(x,y) = - π^4*sin(π*x)*sin(π*y)
w₂₂₂₂(x,y) = - π^4*sin(π*x)*sin(π*y)

D = 1.0
ν = 0.3
M₁₁(x,y) = - D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = - D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = - D*(1-ν)*w₁₂(x,y)
prescribe!(elements["Ω"],:q,(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
prescribe!(elements["Γₚ₁"],:ΔM,(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₂"],:ΔM,(x,y,z)->-2*M₁₂(x,y))
prescribe!(elements["Γₚ₃"],:ΔM,(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₄"],:ΔM,(x,y,z)->-2*M₁₂(x,y))
prescribe!(elements["Γₚ₁"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₂"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₃"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₄"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₁"],:Δn₁s₂n₂s₁,(x,y,z)->2.0)
prescribe!(elements["Γₚ₂"],:Δn₁s₂n₂s₁,(x,y,z)->-2.0)
prescribe!(elements["Γₚ₃"],:Δn₁s₂n₂s₁,(x,y,z)->2.0)
prescribe!(elements["Γₚ₄"],:Δn₁s₂n₂s₁,(x,y,z)->-2.0)
prescribe!(elements["Γ₁"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γ₂"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γ₃"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γ₄"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γ₁"],:V,(x,y,z)-> - D*(-(2-ν)*w₁₁₂(x,y)-w₂₂₂(x,y)))
prescribe!(elements["Γ₂"],:V,(x,y,z)-> - D*(w₁₁₁(x,y)+(2-ν)*w₁₂₂(x,y)))
prescribe!(elements["Γ₃"],:V,(x,y,z)-> - D*((2-ν)*w₁₁₂(x,y)+w₂₂₂(x,y)))
prescribe!(elements["Γ₄"],:V,(x,y,z)-> - D*(-w₁₁₁(x,y)-(2-ν)*w₁₂₂(x,y)))
prescribe!(elements["Γ₁"],:θ,(x,y,z)->-w₂(x,y))
prescribe!(elements["Γ₂"],:θ,(x,y,z)-> w₁(x,y))
prescribe!(elements["Γ₃"],:θ,(x,y,z)-> w₂(x,y))
prescribe!(elements["Γ₄"],:θ,(x,y,z)->-w₁(x,y))
prescribe!(elements["Γ₁"],:M,(x,y,z)->M₂₂(x,y))
prescribe!(elements["Γ₂"],:M,(x,y,z)->M₁₁(x,y))
prescribe!(elements["Γ₃"],:M,(x,y,z)->M₂₂(x,y))
prescribe!(elements["Γ₄"],:M,(x,y,z)->M₁₁(x,y))

coefficient = (:D=>1.0,:α=>1e4,:ν=>0.3)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫MₙₙθdΓ,coefficient...),
       Operator(:∫VgdΓ,coefficient...),
       Operator(:ΔMₙₛg,coefficient...),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ωˢ"],k)
ops[2](elements["Ω"],f)

# ops[3](elements["Γ₁"],k,f)
# ops[3](elements["Γ₂"],k,f)
# ops[3](elements["Γ₃"],k,f)
# ops[3](elements["Γ₄"],k,f)
# ops[6](elements["Γ₁"],f)
# ops[6](elements["Γ₂"],f)
# ops[6](elements["Γ₃"],f)
# ops[6](elements["Γ₄"],f)

ops[4](elements["Γ₁"],k,f)
ops[4](elements["Γ₂"],k,f)
ops[4](elements["Γ₃"],k,f)
ops[4](elements["Γ₄"],k,f)
# ops[7](elements["Γ₁"],f)
# ops[7](elements["Γ₂"],f)
# ops[7](elements["Γ₃"],f)
# ops[7](elements["Γ₄"],f)

ops[5](elements["Γₚ₁"],k,f)
ops[5](elements["Γₚ₂"],k,f)
ops[5](elements["Γₚ₃"],k,f)
ops[5](elements["Γₚ₄"],k,f)
# ops[8](elements["Γₚ₁"],f)
# ops[8](elements["Γₚ₂"],f)
# ops[8](elements["Γₚ₃"],f)
# ops[8](elements["Γₚ₄"],f)

# d = [w(nodes[:x][i],nodes[:y][i]) for i in 1:length(nodes[:x])]
# f .-= k*d

d = k\f

push!(nodes,:d=>d)
set𝓖!(elements["Ω"],:TriGI16)
set_memory_𝝭!(elements["Ω"])

set∇³𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u,(x,y,z)->w(x,y))
prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->w₁(x,y))
prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->w₂(x,y))
prescribe!(elements["Ω"],:∂²u∂x²,(x,y,z)->w₁₁(x,y))
prescribe!(elements["Ω"],:∂²u∂x∂y,(x,y,z)->w₁₂(x,y))
prescribe!(elements["Ω"],:∂²u∂y²,(x,y,z)->w₂₂(x,y))
prescribe!(elements["Ω"],:∂³u∂x³,(x,y,z)->w₁₁₁(x,y))
prescribe!(elements["Ω"],:∂³u∂x²∂y,(x,y,z)->w₁₁₂(x,y))
prescribe!(elements["Ω"],:∂³u∂x∂y²,(x,y,z)->w₁₂₂(x,y))
prescribe!(elements["Ω"],:∂³u∂y³,(x,y,z)->w₂₂₂(x,y))
h3,h2,h1,l2 = ops[9](elements["Ω"])

index = [10,20,40,80]

XLSX.openxlsx("./xlsx/rectangular.xlsx", mode="rw") do xf
    row = "E"
    𝐿₂ = xf[2]
    𝐻₁ = xf[3]
    𝐻₂ = xf[4]
    𝐻₃ = xf[5]
    ind = findfirst(n->n==ndiv,index)+1
    row = row*string(ind)
    𝐿₂[row] = log10(l2)
    𝐻₁[row] = log10(h1)
    𝐻₂[row] = log10(h2)
    𝐻₃[row] = log10(h3)
end
