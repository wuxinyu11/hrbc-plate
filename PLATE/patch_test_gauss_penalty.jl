
using Revise, YAML, ApproxOperator

ndiv = 15
𝒑 = "quartic"
# 𝒑 = "cubic"
config = YAML.load_file("./yml/patch_test_gauss_penalty_"*𝒑*".yml")
elements,nodes = importmsh("./msh/patchtest.msh", config)
# elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)
nₚ = length(nodes)

# s = 3.5/ndiv*ones(nₚ)
s = 4.5/ndiv*ones(nₚ)

# nₚ = length(nodes)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇²₂𝝭!(elements["Ω"])
set∇₂𝝭!(elements["Γ₁"])
set∇₂𝝭!(elements["Γ₂"])
set∇₂𝝭!(elements["Γ₃"])
set∇₂𝝭!(elements["Γ₄"])
set𝝭!(elements["Γₚ₁"])
set𝝭!(elements["Γₚ₂"])
set𝝭!(elements["Γₚ₃"])
set𝝭!(elements["Γₚ₄"])

n = 4

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
prescribe!(elements["Γₚ₁"],:Δn₁s₂n₂s₁=>(x,y,z)->2.0)
prescribe!(elements["Γₚ₂"],:Δn₁s₂n₂s₁=>(x,y,z)->-2.0)
prescribe!(elements["Γₚ₃"],:Δn₁s₂n₂s₁=>(x,y,z)->2.0)
prescribe!(elements["Γₚ₄"],:Δn₁s₂n₂s₁=>(x,y,z)->-2.0)
prescribe!(elements["Γₚ₁"],:ΔM=>(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₂"],:ΔM=>(x,y,z)->-2*M₁₂(x,y))
prescribe!(elements["Γₚ₃"],:ΔM=>(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₄"],:ΔM=>(x,y,z)->-2*M₁₂(x,y))

coefficient = (:D=>D,:ν=>ν)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫vgdΓ,coefficient...,:α=>1e7),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫∇𝑛vθdΓ,coefficient...,:α=>1e7),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)

ops[3](elements["Γ₁"],k,f)
ops[3](elements["Γ₂"],k,f)
ops[3](elements["Γ₃"],k,f)
ops[3](elements["Γ₄"],k,f)

ops[5](elements["Γ₁"],k,f)
ops[5](elements["Γ₂"],k,f)
ops[5](elements["Γ₃"],k,f)
ops[5](elements["Γ₄"],k,f)
# ops[6](elements["Γ₁"],f)
# ops[6](elements["Γ₂"],f)
# ops[6](elements["Γ₃"],f)
# ops[6](elements["Γ₄"],f)

ops[3](elements["Γₚ₁"],k,f)
ops[3](elements["Γₚ₂"],k,f)
ops[3](elements["Γₚ₃"],k,f)
ops[3](elements["Γₚ₄"],k,f)
# ops[7](elements["Γₚ₁"],f)
# ops[7](elements["Γₚ₂"],f)
# ops[7](elements["Γₚ₃"],f)
# ops[7](elements["Γₚ₄"],f)

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
h3,h2,h1,l2 = ops[8](elements["Ω"])
