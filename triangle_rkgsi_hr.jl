
using Revise, YAML, ApproxOperator
config = YAML.load_file("./yml/triangle_rkgsi_hr.yml")
ndiv = 10
elements, nodes = importmsh("./msh/triangle_"*string(ndiv)*".msh", config)


nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 4*10/ndiv*ones(nₚ)
#s = 4.5*10/ndiv*ones(nₚ)
#push!(nodes,:s₁=>3^(0.5)/2 .*s,:s₂=>s,:s₃=>s)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γ₁"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ₂"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ₃"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γₚ₁"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₂"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₃"],:𝝭,:∇̃²)


elements["Ω∩Γ₁"] = elements["Ω"]∩elements["Γ₁"]
elements["Ω∩Γ₂"] = elements["Ω"]∩elements["Γ₂"]
elements["Ω∩Γ₃"] = elements["Ω"]∩elements["Γ₃"]
elements["Ω∩Γₚ₁"] = elements["Ω"]∩elements["Γₚ₁"]
elements["Ω∩Γₚ₂"] = elements["Ω"]∩elements["Γₚ₂"]
elements["Ω∩Γₚ₃"] = elements["Ω"]∩elements["Γₚ₃"]
elements["Γₚ"] = elements["Γₚ₁"]∪elements["Γₚ₂"]∪elements["Γₚ₃"]
elements["Γ"] = elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]
elements["Γ∩Γₚ"] = elements["Γ"]∩elements["Γₚ"]

set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γ₁"],elements["Ω∩Γ₁"])
set∇∇̃²𝝭!(elements["Γ₂"],elements["Ω∩Γ₂"])
set∇∇̃²𝝭!(elements["Γ₃"],elements["Ω∩Γ₃"])

set∇̃²𝝭!(elements["Γₚ₁"],elements["Ω∩Γₚ₁"])
set∇̃²𝝭!(elements["Γₚ₂"],elements["Ω∩Γₚ₂"])
set∇̃²𝝭!(elements["Γₚ₃"],elements["Ω∩Γₚ₃"])

set∇₂𝝭!(elements["Γ₁"])
set∇₂𝝭!(elements["Γ₂"])
set∇₂𝝭!(elements["Γ₃"])

set𝝭!(elements["Γₚ₁"])
set𝝭!(elements["Γₚ₂"])
set𝝭!(elements["Γₚ₃"])


set∇∇̄²𝝭!(elements["Γ₁"],Γᵍ=elements["Γ₁"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₂"],Γᵍ=elements["Γ₂"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"],Γᴾ=elements["Γₚ"])
set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γ∩Γₚ"],Γᴾ=elements["Γₚ"])

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
w(x,y) = 1/640*(x^3-3y^2*x-10(x^2+y^2)+4000/27)*(400/9-x^2-y^2)
w₁(x,y) = 1/640*(3*x^2-3*y^2-20x)*(4/9*100-x^2-y^2)+1/640*(x^3-3y^2*x-10(x^2+y^2)+4/27*1000)*(-2*x)
w₂(x,y) = 1/640*(0-6y*x-20*y)*(4/9*100-x^2-y^2)+1/640*(x^3-3*y^2*x-10(x^2+y^2)+4/27*1000)*(-2*y)
w₁₁(x,y) = 1/640*(6*x-20)*(4/9*100-x^2-y^2)+1/640*(3*x^2-3*y^2-20*x)*(-2*x)*2-2/640*(x^3-3*y^2*x-10(x^2+y^2)+4/27*1000)
w₂₂(x,y) = 1/640*(-6*x-20)*(4/9*100-x^2-y^2)+1/640*(0-6*y*x-20*y)*(-2*y)*2-2/640*(x^3-3*y^2*x-10(x^2+y^2)+4/27*1000)
w₁₂(x,y) = 1/640*(-6*y)*(4/9*100-x^2-y^2)+1/640*(3*x^2-3*y^2-20*x)*(-2*y)+1/640*(0-6*y*x-20*y)*(-2*x)
w₁₁₁(x,y) = 1/640*6*(4/9*100-x^2-y^2)+1/640*(6*x-20)*(-2*x)*3+1/640*(3*x^2-3*y^2-20*x)*(-2)*3
w₁₁₂(x,y) = 0+1/640*(6*x-20)*(-2*y)+1/640*(0-6*y)*(-2*x)*2-2/640*(0-6*y*x-20*y)
w₁₂₂(x,y) = 1/640*(-6)*(4/9*100-x^2-y^2)+1/640*(-6*y)*(-2*y)*2+1/640*(3*x^2-3*y^2-20*x)*(-2)+1/640*(0-6*x-20)*(-2*x)
w₂₂₂(x,y) = 1/640*(-6*x-20)*(-2*y)*3+1/640*(0-6*y*x-20*y)*(-2)*3
w₁₁₁₁(x,y) = 1/640*6*(-2*x)+1/640*(6)*(-2*x)*3+1/640*(6*x-20)*(-2)*3+1/640*(6*x-20)*(-2)*3
w₁₁₂₂(x,y) = 1/640*(6*x-20)*(-2)+1/640*(0-6)*(-2*x)*2-2/640*(0-6*x-20)
w₂₂₂₂(x,y) = 1/640*(-6*x-20)*(-2)*3+1/640*(0-6*x-20)*(-2)*3

D = 1.0
ν = 0.3
M₁₁(x,y) = - D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = - D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = - D*(1-ν)*w₁₂(x,y)

prescribe!(elements["Ω"],:q=>(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
#prescribe!(elements["Γ₁"],:g=>(x,y,z)->w(x,y))
#prescribe!(elements["Γ₂"],:g=>(x,y,z)->w(x,y))
#prescribe!(elements["Γ₃"],:g=>(x,y,z)->w(x,y))

#prescribe!(elements["Γ₁"],:V=>(x,y,z)-> - D*(-(2-ν)*w₁₁₂(x,y)-w₂₂₂(x,y)))
#prescribe!(elements["Γ₂"],:V=>(x,y,z)-> - D*(w₁₁₁(x,y)+(2-ν)*w₁₂₂(x,y)))
#prescribe!(elements["Γ₃"],:V=>(x,y,z)-> - D*((2-ν)*w₁₁₂(x,y)+w₂₂₂(x,y)))

#prescribe!(elements["Γ₁"],:θ=>(x,y,z)->1/2*w₁(x,y)-3^(0.5)/2*w₂(x,y))
#prescribe!(elements["Γ₂"],:θ=>(x,y,z)->-1*w₁(x,y))
#prescribe!(elements["Γ₃"],:θ=>(x,y,z)->1/2*w₁(x,y)+3^(0.5)/2*w₂(x,y))

#prescribe!(elements["Γ₁"],:M=>(x,y,z)->1/4*M₁₁(x,y)+3/4*M₂₂(x,y)-3^(0.5)/2*M₁₂(x,y))
#prescribe!(elements["Γ₂"],:M=>(x,y,z)->M₁₁(x,y))
#prescribe!(elements["Γ₃"],:M=>(x,y,z)->1/4*M₁₁(x,y)+3/4*M₂₂(x,y)+3^(0.5)/2*M₁₂(x,y))

#prescribe!(elements["Γₚ₁"],:g=>(x,y,z)->w(x,y))
#prescribe!(elements["Γₚ₂"],:g=>(x,y,z)->w(x,y))
#prescribe!(elements["Γₚ₃"],:g=>(x,y,z)->w(x,y))

#prescribe!(elements["Γₚ₁"],:Δn₁s₁=>(x,y,z)->-3^(0.5)/2)
#prescribe!(elements["Γₚ₁"],:Δn₁s₂n₂s₁=>(x,y,z)->0)
#prescribe!(elements["Γₚ₁"],:Δn₂s₂=>(x,y,z)->3^(0.5)/2)
#prescribe!(elements["Γₚ₂"],:Δn₁s₁=>(x,y,z)->3^(0.5)/4)
#prescribe!(elements["Γₚ₂"],:Δn₁s₂n₂s₁=>(x,y,z)->-3/2)
#prescribe!(elements["Γₚ₂"],:Δn₂s₂=>(x,y,z)->-3^(0.5)/4)
#prescribe!(elements["Γₚ₃"],:Δn₁s₁=>(x,y,z)->3^(0.5)/4)
#prescribe!(elements["Γₚ₃"],:Δn₁s₂n₂s₁=>(x,y,z)->3/2)
#prescribe!(elements["Γₚ₃"],:Δn₂s₂=>(x,y,z)->-3^(0.5)/4)

#prescribe!(elements["Γₚ₁"],:ΔM=>(x,y,z)->0*M₁₂(x,y))
#prescribe!(elements["Γₚ₂"],:ΔM=>(x,y,z)->3/4*M₁₂(x,y))
#prescribe!(elements["Γₚ₃"],:ΔM=>(x,y,z)->-3/4*M₁₂(x,y))

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

ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)

ops[3](elements["Γ₁"],k,f)
ops[3](elements["Γ₂"],k,f)
ops[3](elements["Γ₃"],k,f)

# # ops[6](elements["Γ₁"],f)
# # ops[6](elements["Γ₂"],f)
# # ops[6](elements["Γ₃"],f)
# # ops[6](elements["Γ₄"],f)

#ops[5](elements["Γ₁"],k,f)
#ops[5](elements["Γ₂"],k,f)
#ops[5](elements["Γ₃"],k,f)

# # ops[7](elements["Γ₁"],f)
# # ops[7](elements["Γ₂"],f)
# # ops[7](elements["Γ₃"],f)
# # ops[7](elements["Γ₄"],f)

ops[7](elements["Γₚ"],k,f)

#ops[7](elements["Γₚ₁"],k,f)
#ops[7](elements["Γₚ₂"],k,f)
#ops[7](elements["Γₚ₃"],k,f)
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

d = k\f

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
H1=log10(h1)
H2=log10(h2)
H3=log10(h3)
L2=log10(l2)
h=log10(1/ndiv)