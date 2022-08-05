

using Revise, YAML, ApproxOperator

config = YAML.load_file("./yml/hollow_cylinder_rkgsi_hr.yml")

ndiv = 80
elements, nodes = importmsh("./msh/hollow_cylinder_"*string(ndiv)*".msh", config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 3.1*π/2/ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γᶿ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²)
set_memory_𝗠!(elements["Γᴾ"],:𝝭,:∇̃²)

elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]
elements["Ω∩Γᶿ"] = elements["Ω"]∩elements["Γᶿ"]
elements["Ω∩Γᴾ"] = elements["Ω"]∩elements["Γᴾ"]
elements["Γ"] = elements["Γᵍ"]∪elements["Γᶿ"]
elements["Γ∩Γᴾ"] = elements["Γ"]∩elements["Γᴾ"]


set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set𝝭!(elements["Γᵍ"])
set∇̃²𝝭!(elements["Γᶿ"],elements["Ω∩Γᶿ"])
set∇₂𝝭!(elements["Γᶿ"])
set∇₂𝝭!(elements["Γᴹ"])
set𝝭!(elements["Γⱽ"])
set∇̃²𝝭!(elements["Γᴾ"],elements["Ω∩Γᴾ"])
set𝝭!(elements["Γᴾ"])

set∇∇̄²𝝭!(elements["Γᵍ"],Γᵍ=elements["Γᵍ"],Γᴾ=elements["Γᴾ"])
set∇̄²𝝭!(elements["Γᶿ"],Γᶿ=elements["Γᶿ"],Γᴾ=elements["Γᴾ"])
set∇̄²𝝭!(elements["Γᴾ"],Γᵍ=elements["Γᵍ"],Γᶿ=elements["Γᶿ"],Γᴾ=elements["Γᴾ"])

# n = 4
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


w(x,y)=4/(3*(1-ν))*log((x^2+y^2)^(1/2)/2)-1/(3*(1+ν))*(x^2+y^2-4)
w₁(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-1)*x-2*x/(3*(1+ν))
w₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-1)*y-2*y/(3*(1+ν))
w₁₁(x,y)=-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*x^2+4/(3*(1-ν))*(x^2+y^2)^(-1)-1/(3*(1+ν))*2
w₁₂(x,y)=-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*y*x
w₂₂(x,y)=-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*y^2+4/(3*(1-ν))*(x^2+y^2)^(-1)-1/(3*(1+ν))*2
w₁₁₁(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*x^3-24*x/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₂(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*y*x^2-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*y
w₁₂₂(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*x*y^2-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*x
w₂₂₂(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*y^3-24*y/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₁₁(x,y)=-24/(3*(1-ν))*(x^2+y^2)^(-4)*8*x^4+192*x^2/(3*(1-ν))*(x^2+y^2)^(-3)-24/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₂₂(x,y)=-24/(3*(1-ν))*(x^2+y^2)^(-4)*8*x^2*y^2+32/(3*(1-ν))*(x^2+y^2)^(-3)*x^2+32/(3*(1-ν))*(x^2+y^2)^(-3)*y^2-8/(3*(1-ν))*(x^2+y^2)^(-2)
w₂₂₂₂(x,y)=-24/(3*(1-ν))*(x^2+y^2)^(-4)*8*y^4+192*y^2/(3*(1-ν))*(x^2+y^2)^(-3)-24/(3*(1-ν))*(x^2+y^2)^(-2)

D = 1.0
ν = 0.3
M₁₁(x,y) = - D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = - D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = - D*(1-ν)*w₁₂(x,y)
function Vₙ(x,y,n₁,n₂)
    s₁ = -n₂
    s₂ = n₁
    D₁₁₁ = -D*(n₁ + n₁*s₁*s₁ + ν*n₂*s₁*s₂)
    D₁₁₂ = -D*(n₂ + n₂*s₁*s₁ + 2*n₁*s₁*s₂ + (n₂*s₂*s₂ - n₂*s₁*s₁ - n₁*s₁*s₂)*ν)
    D₁₂₂ = -D*(n₁ + n₁*s₂*s₂ + 2*n₂*s₁*s₂ + (n₁*s₁*s₁ - n₁*s₂*s₂ - n₂*s₁*s₂)*ν)
    D₂₂₂ = -D*(n₂ + n₂*s₂*s₂ + ν*n₁*s₁*s₂)
    return D₁₁₁*w₁₁₁(x,y)+D₁₁₂*w₁₁₂(x,y)+D₁₂₂*w₁₂₂(x,y)+D₂₂₂*w₂₂₂(x,y)
end

prescribe!(elements["Ω"],:q=>(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γᶿ"],:θ=>(x,y,z,n₁,n₂)->w₁(x,y)*n₁+w₂(x,y)*n₂)
set𝒏!(elements["Γᴹ"])
prescribe!(elements["Γᴹ"],:M=>(x,y,z,n₁,n₂)->M₁₁(x,y)*n₁*n₁+2*M₁₂(x,y)*n₁*n₂+M₂₂(x,y)*n₂*n₂)
set𝒏!(elements["Γⱽ"])
prescribe!(elements["Γⱽ"],:V=>(x,y,z,n₁,n₂)->Vₙ(x,y,n₁,n₂))
prescribe!(elements["Γᴾ"],:g=>(x,y,z)->w(x,y))


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

ops[3](elements["Γᵍ"],k,f)
ops[4](elements["Γⱽ"],f)
ops[5](elements["Γᶿ"],k,f)
ops[6](elements["Γᴹ"],f)
ops[7](elements["Γᴾ"],k,f)

d = k\f

# d = [w(n.x,n.y) for n in nodes]
# f .-= k*d

push!(nodes,:d=>d)
set𝓖!(elements["Ω"],:TriGI16,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
set∇̂³𝝭!(elements["Ω"])
# set∇³𝝭!(elements["Ω"])
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
# h=log10(π/2/ndiv)
# H1=log10(h1)
# H2=log10(h2)
# H3=log10(h3)
L2=log10(l2)
