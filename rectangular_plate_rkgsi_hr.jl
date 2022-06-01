
using XLSX, YAML, ApproxOperator

config = YAML.load_file("./yml/rectangular_hrrk.yml")

ndiv = 20
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)

nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

s = 3.5 / 20 * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

sp = RegularGrid(nodes[:x], nodes[:y], nodes[:z], n = 2, γ = 5)
sp(elements["Ω"])
set_memory_𝝭!(elements["Ωˢ"])
set_memory_𝝭!(elements["Γ̃"])

elements["Ω∩Γ̃"] = elements["Ω"]∩elements["Γ̃"]

set∇𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ωˢ"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γ̃"],elements["Ω∩Γ̃"])
set𝝭!(elements["Γ̃"])
set∇∇̄²𝝭!(elements["Γ̃"],Γᵍ=elements["Γ̃"])

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
prescribe!(elements["Γ̃"],:g,(x,y,z)->w(x,y))

coefficient = (:D=>1.0,:α=>1e3,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫ṼgdΓ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ωˢ"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γ̃"],k,f)

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
h3,h2,h1,l2 = ops[4](elements["Ω"])

index = [10,20,40,80]
row = "F"
XLSX.openxlsx("rectangular.xlsx", mode="rw") do xf
    𝐿₂ = xf[1]
    𝐻₁ = xf[2]
    𝐻₂ = xf[3]
    𝐻₃ = xf[4]
    ind = findfirst(n->n==ndiv,index)+1
    row = row*string(ind)
    𝐿₂[row] = log10(l2)
    𝐻₁[row] = log10(h1)
    𝐻₂[row] = log10(h2)
    𝐻₃[row] = log10(h3)
end
