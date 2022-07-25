

using Revise, YAML, ApproxOperator

config = YAML.load_file("./yml/hollow_cylinder_rkgsi_hr.yml")

ndiv = 10
elements, nodes = importmsh("./msh/hollow_cylinder_"*string(ndiv)*".msh", config)
# elements, nodes = importmsh("./msh/patchtest_"*string(ndiv)*".msh", config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γ̃₁"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ̃₂"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ̃₃"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ̃₄"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γₚ₁"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₂"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₃"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₄"],:𝝭,:∇̃²)

elements["Ω∩Γ₁"] = elements["Ω"]∩elements["Γ̃₁"]
elements["Ω∩Γ₂"] = elements["Ω"]∩elements["Γ̃₂"]
elements["Ω∩Γ₃"] = elements["Ω"]∩elements["Γ̃₃"]
elements["Ω∩Γ₄"] = elements["Ω"]∩elements["Γ̃₄"]
elements["Ω∩Γₚ₁"] = elements["Ω"]∩elements["Γₚ₁"]
elements["Ω∩Γₚ₂"] = elements["Ω"]∩elements["Γₚ₂"]
elements["Ω∩Γₚ₃"] = elements["Ω"]∩elements["Γₚ₃"]
elements["Ω∩Γₚ₄"] = elements["Ω"]∩elements["Γₚ₄"]
elements["Γₚ"] = elements["Γₚ₁"]∪elements["Γₚ₂"]∪elements["Γₚ₃"]∪elements["Γₚ₄"]
elements["Γ"] = elements["Γ̃₁"]∪elements["Γ̃₂"]∪elements["Γ̃₃"]∪elements["Γ̃₄"]
elements["Γ∩Γₚ"] = elements["Γ"]∩elements["Γₚ"]


set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γ̃₁"],elements["Ω∩Γ₁"])
set∇∇̃²𝝭!(elements["Γ̃₂"],elements["Ω∩Γ₂"])
set∇∇̃²𝝭!(elements["Γ̃₃"],elements["Ω∩Γ₃"])
set∇∇̃²𝝭!(elements["Γ̃₄"],elements["Ω∩Γ₄"])
set∇₂𝝭!(elements["Γ₁"])
set∇₂𝝭!(elements["Γ₂"])
set∇₂𝝭!(elements["Γ₃"])
set∇₂𝝭!(elements["Γ₄"])
set∇̃²𝝭!(elements["Γₚ₁"],elements["Ω∩Γₚ₁"])
set∇̃²𝝭!(elements["Γₚ₂"],elements["Ω∩Γₚ₂"])
set∇̃²𝝭!(elements["Γₚ₃"],elements["Ω∩Γₚ₃"])
set∇̃²𝝭!(elements["Γₚ₄"],elements["Ω∩Γₚ₄"])
set∇₂𝝭!(elements["Γ̃₁"])
set∇₂𝝭!(elements["Γ̃₂"])
set∇₂𝝭!(elements["Γ̃₃"])
set∇₂𝝭!(elements["Γ̃₄"])

set𝝭!(elements["Γₚ₁"])
set𝝭!(elements["Γₚ₂"])
set𝝭!(elements["Γₚ₃"])
set𝝭!(elements["Γₚ₄"])
set∇∇̄²𝝭!(elements["Γ̃₁"],Γᶿ=elements["Γ̃₁"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ̃₂"],Γᵍ=elements["Γ̃₂"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ̃₃"],Γᶿ=elements["Γ̃₃"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ̃₄"],Γᵍ=elements["Γ̃₄"],Γᴾ=elements["Γₚ"])
set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γ̃₂"]∪elements["Γ̃₄"],Γᶿ=elements["Γ̃₁"]∪elements["Γ̃₃"],Γᴾ=elements["Γₚ"])

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
n = 2
w(x,y) = (1+2x+3y)^n
w₁(x,y) = 2n*(1+2x+3y)^abs(n-1)
w₂(x,y) = 3n*(1+2x+3y)^abs(n-1)
w₁₁(x,y) = 4n*(n-1)*(1+2x+3y)^abs(n-2)
w₂₂(x,y) = 9n*(n-1)*(1+2x+3y)^abs(n-2)
w₁₂(x,y) = 6n*(n-1)*(1+2x+3y)^abs(n-2)
w₁₁₁(x,y) = 8n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
w₁₁₂(x,y) = 12n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
w₁₂₂(x,y) = 18n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
w₂₂₂(x,y) = 27n*(n-1)*(n-2)*(1+2x+3y)^abs(n-3)
w₁₁₁₁(x,y) = 16n*(n-1)*(n-2)*(n-3)*(1+2x+3y)^abs(n-4)
w₁₁₂₂(x,y) = 36n*(n-1)*(n-2)*(n-3)*(1+2x+3y)^abs(n-4)
w₂₂₂₂(x,y) = 81n*(n-1)*(n-2)*(n-3)*(1+2x+3y)^abs(n-4)

# w(x,y) = 4/3/(1-ν)*log((x^2+y^2)/2)-1/3/*(1+ν)*(x^2+y^2-4)
# w₁(x,y) = 4/3/(1-ν)*(x^2+y^2)^(-1)*2*x-2*x/3/(1+ν)
# w₂(x,y) = 4/3/(1-ν)*(x^2+y^2)^(-1)*2*y-2*y/3/(1+ν)
# w₁₁(x,y) = 4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*x*2*x+2*4/3/(1-ν)*(x^2+y^2)^(-1)-2/3/(1+ν)
# w₂₂(x,y) = 4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*y*2*y+2*4/3/(1-ν)*(x^2+y^2)^(-1)-2/3/(1+ν)
# w₁₂(x,y) = 4/3/(1-ν)*(-1)(x^2+y^2)^(-2)*2*y*2*x
# w₁₁₁(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*2*x*2*x*2*x+8*x*4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*x
# w₁₁₂(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*2*y*4*x^2+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*y
# w₁₂₂(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*2*x*4*y^2+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*x
# w₂₂₂(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*8*y^3+4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*8*y+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*y
# w₁₁₁₁(x,y) = 4/3/(1-ν)*(-6)*(x^2+y^2)^(-4)*16*x^4+4/3/(1-ν)*2*(x^2+y^2)^(-3)*24*x^2-32/3/(1-ν)*(x^2+y^2)^(-2)+8/3*16*x^2*(x^2+y^2)^(-3)+16/3/(1-ν)*(x^2+y^2)^(-3)*4*x^2-16/3/(1-ν)*(x^2+y^2)^(-2)
# w₁₁₂₂(x,y) = 4/3/(1-ν)*(-6)*(x^2+y^2)^(-4)*4*y^2*4*x^2+8/3*(x^2+y^2)^(-3)*8*x^2+16/3/(1-ν)*(x^2+y^2)^(-3)*4*y^2-16/3/(1-ν)*(x^2+y^2)^(-2)
# w₂₂₂₂(x,y) = 4/3/(1-ν)*(-6)*(x^2+y^2)^(-4)*16*y^2+3/4/(1-ν)*8*(x^2+y^2)^(-3)*24*y^2*16*4+4/3/(1-ν)*16*2*(x^2+y^2)^(-2)
D = 1.0
ν = 0.3
M₁₁(x,y) = - D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = - D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = - D*(1-ν)*w₁₂(x,y)

prescribe!(elements["Ω"],:q=>(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
prescribe!(elements["Γ̃₂"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γ̃₄"],:g=>(x,y,z)->w(x,y))

# prescribe!(elements["Γ₁"],:V=>(x,y,z)->0.0)
# prescribe!(elements["Γ₁"],:V=>(x,y,z)-> - D*(-(2-ν)*w₁₁₂(x,y)-w₂₂₂(x,y)))
# prescribe!(elements["Γ₃"],:V=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:V=>(x,y,z)-> - D*((2-ν)*w₁₁₂(x,y)+w₂₂₂(x,y)))

# prescribe!(elements["Γ̃₁"],:θ=>(x,y,z)->0.0)
# prescribe!(elements["Γ̃₁"],:θ=>(x,y,z)->-w₂(x,y))
prescribe!(elements["Γ̃₁"],:θ=>(x,y,z)->-w₂(x,y))
# prescribe!(elements["Γ̃₃"],:θ=>(x,y,z)-> 0.0)
prescribe!(elements["Γ̃₃"],:θ=>(x,y,z)->-w₁(x,y))
# prescribe!(elements["Γ̃₃"],:θ=>(x,y,z)->w₂(x,y))

# prescribe!(elements["Γ₂"],:M=>(x,y,z)->1/2*M₁₁(x,y)+1/2*M₂₂(x,y)+M₁₂(x,y))
# prescribe!(elements["Γ₄"],:M=>(x,y,z)->1/2*M₁₁(x,y)+1/2*M₂₂(x,y)+M₁₂(x,y))
prescribe!(elements["Γ₂"],:M=>(x,y,z,n₁,n₂,s₁,s₂)->M₁₁(x,y)*n₁*n₁+2*M₁₂(x,y)*n₁*n₂+M₂₂(x,y)*n₂*n₂)
prescribe!(elements["Γ₄"],:M=>(x,y,z,n₁,n₂,s₁,s₂)->M₁₁(x,y)*n₁*n₁+2*M₁₂(x,y)*n₁*n₂+M₂₂(x,y)*n₂*n₂)
# prescribe!(elements["Γ₂"],:M=>(x,y,z)->M₁₁(x,y))
# prescribe!(elements["Γ₄"],:M=>(x,y,z)->M₁₁(x,y))

prescribe!(elements["Γₚ₁"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₂"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₃"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₄"],:g=>(x,y,z)->w(x,y))


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

ops[3](elements["Γ̃₂"],k,f)
ops[3](elements["Γ̃₄"],k,f)

 ops[4](elements["Γ₁"],f)
 ops[4](elements["Γ₃"],f)

ops[5](elements["Γ̃₁"],k,f)
ops[5](elements["Γ̃₃"],k,f)

ops[6](elements["Γ₂"],f)
ops[6](elements["Γ₄"],f)

ops[7](elements["Γₚ"],k,f)

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
# H1=log10(h1)
# H2=log10(h2)
# H3=log10(h3)
# L2=log10(l2)
# h=log10(1/ndiv)
