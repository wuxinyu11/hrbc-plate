
using Revise, XLSX, YAML, ApproxOperator

config = YAML.load_file("./yml/rectangular_hrrk.yml")

ndiv = 80
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)

nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

s = 3.1 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

sp = RegularGrid(nodes[:x], nodes[:y], nodes[:z], n = 2, γ = 5)
sp(elements["Ω"])
set_memory_𝝭!(elements["Ωˢ"])
set_memory_𝝭!(elements["Γ̃"])
set_memory_𝝭!(elements["Γ̃ₚ₁"])
set_memory_𝝭!(elements["Γ̃ₚ₂"])
set_memory_𝝭!(elements["Γ̃ₚ₃"])
set_memory_𝝭!(elements["Γ̃ₚ₄"])

elements["Ω∩Γ̃"] = elements["Ω"]∩elements["Γ̃"]
elements["Ω∩Γ̃ₚ₁"] = elements["Ω"]∩elements["Γ̃ₚ₁"]
elements["Ω∩Γ̃ₚ₂"] = elements["Ω"]∩elements["Γ̃ₚ₂"]
elements["Ω∩Γ̃ₚ₃"] = elements["Ω"]∩elements["Γ̃ₚ₃"]
elements["Ω∩Γ̃ₚ₄"] = elements["Ω"]∩elements["Γ̃ₚ₄"]
elements["Γ̃ₚ"] = elements["Γ̃ₚ₁"]∪elements["Γ̃ₚ₂"]∪elements["Γ̃ₚ₃"]∪elements["Γ̃ₚ₄"]
elements["Γ̃∩Γ̃ₚ"] = elements["Γ̃"]∩elements["Γ̃ₚ"]

set∇𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ωˢ"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γ̃"],elements["Ω∩Γ̃"])
set𝝭!(elements["Γ̃"])
set∇̃²𝝭!(elements["Γ̃ₚ₁"],elements["Ω∩Γ̃ₚ₁"])
set∇̃²𝝭!(elements["Γ̃ₚ₂"],elements["Ω∩Γ̃ₚ₂"])
set∇̃²𝝭!(elements["Γ̃ₚ₃"],elements["Ω∩Γ̃ₚ₃"])
set∇̃²𝝭!(elements["Γ̃ₚ₄"],elements["Ω∩Γ̃ₚ₄"])
set𝝭!(elements["Γ̃ₚ₁"])
set𝝭!(elements["Γ̃ₚ₂"])
set𝝭!(elements["Γ̃ₚ₃"])
set𝝭!(elements["Γ̃ₚ₄"])
set∇∇̄²𝝭!(elements["Γ̃"],Γᵍ=elements["Γ̃"],Γᴾ=elements["Γ̃ₚ"])
set∇̄²𝝭!(elements["Γ̃ₚ"],Γᵍ=elements["Γ̃∩Γ̃ₚ"],Γᴾ=elements["Γ̃ₚ"])

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
prescribe!(elements["Ω"],:q,(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
# prescribe!(elements["Γ̃ₚ₁"],:Δn₁s₂n₂s₁,(x,y,z)->2.0)
# prescribe!(elements["Γ̃ₚ₂"],:Δn₁s₂n₂s₁,(x,y,z)->-2.0)
# prescribe!(elements["Γ̃ₚ₃"],:Δn₁s₂n₂s₁,(x,y,z)->2.0)
# prescribe!(elements["Γ̃ₚ₄"],:Δn₁s₂n₂s₁,(x,y,z)->-2.0)

coefficient = (:D=>1.0,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫ṼgdΓ,coefficient...),
       Operator(:ΔM̃ₙₛg,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ωˢ"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γ̃"],k,f)
ops[4](elements["Γ̃ₚ"],k,f)

d = k\f

push!(nodes,:d=>d)
set𝓖!(elements["Ω"],:TriGI16)
set_memory_𝝭!(elements["Ω"])

set∇³𝝭!(elements["Ω"])
# set∇̂³𝝭!(elements["Ω"])
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
h3,h2,h1,l2 = ops[5](elements["Ω"])

f = checkConsistency(elements["Ω"],ApproxOperator.get∇³𝝭,ApproxOperator.get∇³𝒑)

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
