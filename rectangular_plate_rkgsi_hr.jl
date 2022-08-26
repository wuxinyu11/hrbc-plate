
using YAML, ApproxOperator, XLSX, TimerOutputs
# @CPUtime begin
to = TimerOutput()
@timeit to "Total Time" begin
@timeit to "searching" begin

ndiv = 80
# 𝒑 = "cubic"
𝒑 = "quartic"
config = YAML.load_file("./yml/rectangular_rkgsi_hr_"*𝒑*".yml")
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)
end 
nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
       
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γ₁"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ₂"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ₃"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ₄"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γₚ₁"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₂"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₃"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₄"],:𝝭,:∇̃²)

elements["Ω∩Γ₁"] = elements["Ω"]∩elements["Γ₁"]
elements["Ω∩Γ₂"] = elements["Ω"]∩elements["Γ₂"]
elements["Ω∩Γ₃"] = elements["Ω"]∩elements["Γ₃"]
elements["Ω∩Γ₄"] = elements["Ω"]∩elements["Γ₄"]
elements["Ω∩Γₚ₁"] = elements["Ω"]∩elements["Γₚ₁"]
elements["Ω∩Γₚ₂"] = elements["Ω"]∩elements["Γₚ₂"]
elements["Ω∩Γₚ₃"] = elements["Ω"]∩elements["Γₚ₃"]
elements["Ω∩Γₚ₄"] = elements["Ω"]∩elements["Γₚ₄"]
elements["Γₚ"] = elements["Γₚ₁"]∪elements["Γₚ₂"]∪elements["Γₚ₃"]∪elements["Γₚ₄"]
elements["Γ"] = elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]
elements["Γ∩Γₚ"] = elements["Γ"]∩elements["Γₚ"]

 
@timeit to "shape functions " begin      
set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
@timeit to "shape functions Γᵍ " begin      
set∇∇̃²𝝭!(elements["Γ₁"],elements["Ω∩Γ₁"])
set∇∇̃²𝝭!(elements["Γ₂"],elements["Ω∩Γ₂"])
set∇∇̃²𝝭!(elements["Γ₃"],elements["Ω∩Γ₃"])
set∇∇̃²𝝭!(elements["Γ₄"],elements["Ω∩Γ₄"])
set∇̃²𝝭!(elements["Γₚ₁"],elements["Ω∩Γₚ₁"])
set∇̃²𝝭!(elements["Γₚ₂"],elements["Ω∩Γₚ₂"])
set∇̃²𝝭!(elements["Γₚ₃"],elements["Ω∩Γₚ₃"])
set∇̃²𝝭!(elements["Γₚ₄"],elements["Ω∩Γₚ₄"])
set∇₂𝝭!(elements["Γ₁"])
set∇₂𝝭!(elements["Γ₂"])
set∇₂𝝭!(elements["Γ₃"])
set∇₂𝝭!(elements["Γ₄"])
set𝝭!(elements["Γₚ₁"])
set𝝭!(elements["Γₚ₂"])
set𝝭!(elements["Γₚ₃"])
set𝝭!(elements["Γₚ₄"])

# set∇∇̄²𝝭!(elements["Γ₁"],Γᵍ=elements["Γ₁"],Γᶿ=elements["Γ₁"],Γᴾ=elements["Γₚ"])
# set∇∇̄²𝝭!(elements["Γ₂"],Γᵍ=elements["Γ₂"],Γᶿ=elements["Γ₂"],Γᴾ=elements["Γₚ"])
# set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"],Γᶿ=elements["Γ₃"],Γᴾ=elements["Γₚ"])
# set∇∇̄²𝝭!(elements["Γ₄"],Γᵍ=elements["Γ₄"],Γᶿ=elements["Γ₄"],Γᴾ=elements["Γₚ"])
# set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γ∩Γₚ"],Γᶿ=elements["Γ∩Γₚ"],Γᴾ=elements["Γₚ"])

set∇∇̄²𝝭!(elements["Γ₁"],Γᵍ=elements["Γ₁"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₂"],Γᵍ=elements["Γ₂"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₄"],Γᵍ=elements["Γ₄"],Γᴾ=elements["Γₚ"])
set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γ∩Γₚ"],Γᴾ=elements["Γₚ"])

end
end


# set∇∇̄²𝝭!(elements["Γ₁"],Γᵍ=elements["Γ₁"],Γᶿ=elements["Γ₁"])
# set∇∇̄²𝝭!(elements["Γ₂"],Γᵍ=elements["Γ₂"],Γᶿ=elements["Γ₂"])
# set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"],Γᶿ=elements["Γ₃"])
# set∇∇̄²𝝭!(elements["Γ₄"],Γᵍ=elements["Γ₄"],Γᶿ=elements["Γ₄"])

# set∇∇̄²𝝭!(elements["Γ₁"],Γᵍ=elements["Γ₁"])
# set∇∇̄²𝝭!(elements["Γ₂"],Γᵍ=elements["Γ₂"])
# set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"])
# set∇∇̄²𝝭!(elements["Γ₄"],Γᵍ=elements["Γ₄"])

# set∇̄²𝝭!(elements["Γ̃ₚ"],Γᴾ=elements["Γ̃ₚ"])

# set∇∇̄²𝝭!(elements["Γ₁"],Γᵍ=elements["Γ₁"],Γᴾ=elements["Γₚ"])
# set∇∇̄²𝝭!(elements["Γ₂"],Γᵍ=elements["Γ₂"],Γᴾ=elements["Γₚ"])
# set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"],Γᴾ=elements["Γₚ"])
# set∇∇̄²𝝭!(elements["Γ₄"],Γᵍ=elements["Γ₄"],Γᴾ=elements["Γₚ"])
# set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γ∩Γₚ"],Γᴾ=elements["Γₚ"])
# set∇̄²𝝭!(elements["Γ̃ₚ"],Γᵍ=elements["Γ̃₁"],Γᴾ=elements["Γ̃ₚ"])

# set∇∇̄²𝝭!(elements["Γ̃₁"],Γᶿ=elements["Γ̃₁"],Γᴾ=elements["Γ̃ₚ"])
# set∇̄²𝝭!(elements["Γ̃₁"],Γᶿ=elements["Γ̃₁"],Γᴾ=elements["Γ̃ₚ"])
# set∇̄²𝝭!(elements["Γ̃ₚ"],Γᶿ=elements["Γ̃₁"],Γᴾ=elements["Γ̃ₚ"])

# set∇̄²𝝭!(elements["Γ̃₁"],Γᶿ=elements["Γ̃₁"])

# set∇̄²𝝭!(elements["Γ̃₁"],Γᴾ=elements["Γ̃ₚ"])
# set∇̄²𝝭!(elements["Γ̃ₚ"],Γᶿ=elements["Γ̃₁"])

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
prescribe!(elements["Ω"],:q=>(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
prescribe!(elements["Γ₁"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γ₂"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γ₃"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γ₄"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γ₁"],:V=>(x,y,z)-> - D*(-(2-ν)*w₁₁₂(x,y)-w₂₂₂(x,y)))
prescribe!(elements["Γ₂"],:V=>(x,y,z)-> - D*(w₁₁₁(x,y)+(2-ν)*w₁₂₂(x,y)))
prescribe!(elements["Γ₃"],:V=>(x,y,z)-> - D*((2-ν)*w₁₁₂(x,y)+w₂₂₂(x,y)))
prescribe!(elements["Γ₄"],:V=>(x,y,z)-> - D*(-w₁₁₁(x,y)-(2-ν)*w₁₂₂(x,y)))
prescribe!(elements["Γ₁"],:θ=>(x,y,z)->-w₂(x,y))
prescribe!(elements["Γ₂"],:θ=>(x,y,z)-> w₁(x,y))
prescribe!(elements["Γ₃"],:θ=>(x,y,z)-> w₂(x,y))
prescribe!(elements["Γ₄"],:θ=>(x,y,z)->-w₁(x,y))
prescribe!(elements["Γ₁"],:M=>(x,y,z)->M₂₂(x,y))
prescribe!(elements["Γ₂"],:M=>(x,y,z)->M₁₁(x,y))
prescribe!(elements["Γ₃"],:M=>(x,y,z)->M₂₂(x,y))
prescribe!(elements["Γ₄"],:M=>(x,y,z)->M₁₁(x,y))
prescribe!(elements["Γₚ₁"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₂"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₃"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₄"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₁"],:ΔM=>(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₂"],:ΔM=>(x,y,z)->-2*M₁₂(x,y))
prescribe!(elements["Γₚ₃"],:ΔM=>(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₄"],:ΔM=>(x,y,z)->-2*M₁₂(x,y))


coefficient = (:D=>1.0,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫ṼgdΓ,coefficient...),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫M̃ₙₙθdΓ,coefficient...),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔM̃ₙₛg,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

@timeit to "assembly " begin       
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
@timeit to "assembly  Γᵍ" begin       

ops[3](elements["Γ₁"],k,f)
ops[3](elements["Γ₂"],k,f)
ops[3](elements["Γ₃"],k,f)
ops[3](elements["Γ₄"],k,f)
# # ops[6](elements["Γ₁"],f)
# # ops[6](elements["Γ₂"],f)
# # ops[6](elements["Γ₃"],f)
# # ops[6](elements["Γ₄"],f)

# ops[5](elements["Γ₁"],k,f)
# ops[5](elements["Γ₂"],k,f)
# ops[5](elements["Γ₃"],k,f)
# ops[5](elements["Γ₄"],k,f)
# # ops[7](elements["Γ₁"],f)
# # ops[7](elements["Γ₂"],f)
# # ops[7](elements["Γ₃"],f)
# # ops[7](elements["Γ₄"],f)

ops[7](elements["Γₚ"],k,f)
# ops[5](elements["Γ̃ₚ₁"],k,f)
# ops[5](elements["Γ̃ₚ₂"],k,f)
# ops[5](elements["Γ̃ₚ₃"],k,f)
# ops[5](elements["Γ̃ₚ₄"],k,f)
# ops[8](elements["Γₚ₁"],f)
# ops[8](elements["Γₚ₂"],f)
# ops[8](elements["Γₚ₃"],f)
# ops[8](elements["Γₚ₄"],f)
#
# # d = [w(nodes[:x][i],nodes[:y][i]) for i in 1:length(nodes[:x])]
# # f .-= k*d
end
end

d = k\f
end

push!(nodes,:d=>d)
set𝓖!(elements["Ω"],:TriGI16,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
set∇̂³𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u=>(x,y,z)->w(x,y))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->w₁(x,y))
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->w₂(x,y))
prescribe!(elements["Ω"],:∂²u∂x²=>(x,y,z)->w₁₁(x,y))
prescribe!(elements["Ω"],:∂²u∂x∂y=>(x,y,z)->w₁₂(x,y))
prescribe!(elements["Ω"],:∂²u∂y²=>(x,y,z)->w₂₂(x,y))
prescribe!(elements["Ω"],:∂³u∂x³=>(x,y,z)->w₁₁₁(x,y))
prescribe!(elements["Ω"],:∂³u∂x²∂y=>(x,y,z)->w₁₁₂(x,y))
prescribe!(elements["Ω"],:∂³u∂x∂y²=>(x,y,z)->w₁₂₂(x,y))
prescribe!(elements["Ω"],:∂³u∂y³=>(x,y,z)->w₂₂₂(x,y))
h3,h2,h1,l2 = ops[9](elements["Ω"])
show(to)

index = [10,20,40,80]
XLSX.openxlsx("./xlsx/rectangular_"*𝒑*".xlsx", mode="rw") do xf
    row = "G"
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
