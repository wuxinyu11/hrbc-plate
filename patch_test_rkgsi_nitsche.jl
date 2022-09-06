
using Revise, YAML, ApproxOperator

ndiv = 10
# 𝒑 = "quartic"
𝒑 = "cubic"
config = YAML.load_file("./yml/patch_test_rkgsi_nitsche_"*𝒑*".yml")
# elements,nodes = importmsh("./msh/patchtest.msh", config)
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)
nₚ = length(nodes)
s = 4.5/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)

set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇³𝝭!(elements["Γ₁"])
set∇³𝝭!(elements["Γ₂"])
set∇³𝝭!(elements["Γ₃"])
set∇³𝝭!(elements["Γ₄"])
set∇²₂𝝭!(elements["Γₚ₁"])
set∇²₂𝝭!(elements["Γₚ₂"])
set∇²₂𝝭!(elements["Γₚ₃"])
set∇²₂𝝭!(elements["Γₚ₄"])

f1 = ApproxOperator.check∇̃²₂𝝭(elements["Ω̃"])
f2 = check∇³𝝭(elements["Γ₁"])

n = 3
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
w(x,y) = x^n
w₁(x,y) = n*x^abs(n-1)
w₂(x,y) = 0.
w₁₁(x,y) = n*(n-1)*x^abs(n-2)
w₂₂(x,y) = 0.
w₁₂(x,y) = 0.
w₁₁₁(x,y) = n*(n-1)*(n-2)*x^abs(n-3)
w₁₁₂(x,y) = 0.
w₁₂₂(x,y) = 0.
w₂₂₂(x,y) = 0.
w₁₁₁₁(x,y) = n*(n-1)*(n-2)*(n-3)*x^abs(n-4)
w₁₁₂₂(x,y) = 0.
w₂₂₂₂(x,y) = 0.
D = 1.0
ν = 0.0
M₁₁(x,y) = - D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = - D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = - D*(1-ν)*w₁₂(x,y)
set𝒏!(elements["Γ₁"])
set𝒏!(elements["Γ₂"])
set𝒏!(elements["Γ₃"])
set𝒏!(elements["Γ₄"])
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
# prescribe!(elements["Γₚ₁"],:Δn₁s₂n₂s₁=>(x,y,z)->2.0)
# prescribe!(elements["Γₚ₂"],:Δn₁s₂n₂s₁=>(x,y,z)->-2.0)
# prescribe!(elements["Γₚ₃"],:Δn₁s₂n₂s₁=>(x,y,z)->2.0)
# prescribe!(elements["Γₚ₄"],:Δn₁s₂n₂s₁=>(x,y,z)->-2.0)
prescribe!(elements["Γₚ₁"],:ΔM=>(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₂"],:ΔM=>(x,y,z)->-2*M₁₂(x,y))
prescribe!(elements["Γₚ₃"],:ΔM=>(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₄"],:ΔM=>(x,y,z)->-2*M₁₂(x,y))

coefficient = (:D=>D,:ν=>ν)
ops = [Operator(:∫κ̃ᵢⱼM̃ᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫VgdΓ,coefficient...,:α=>1e3*ndiv^2),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫MₙₙθdΓ,coefficient...,:α=>1e3*ndiv),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔMₙₛg,coefficient...,:α=>1e3*ndiv^2),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)

# ops[3](elements["Γ₁"],k,f)
# ops[3](elements["Γ₂"],k,f)
# ops[3](elements["Γ₃"],k,f)
# ops[3](elements["Γ₄"],k,f)
ops[4](elements["Γ₁"],f)
ops[4](elements["Γ₂"],f)
ops[4](elements["Γ₃"],f)
ops[4](elements["Γ₄"],f)

# ops[5](elements["Γ₁"],k,f)
# ops[5](elements["Γ₂"],k,f)
# ops[5](elements["Γ₃"],k,f)
# ops[5](elements["Γ₄"],k,f)
ops[6](elements["Γ₁"],f)
ops[6](elements["Γ₂"],f)
ops[6](elements["Γ₃"],f)
ops[6](elements["Γ₄"],f)

# ops[7](elements["Γₚ₁"],k,f)
# ops[7](elements["Γₚ₂"],k,f)
# ops[7](elements["Γₚ₃"],k,f)
# ops[7](elements["Γₚ₄"],k,f)
ops[8](elements["Γₚ₁"][1],f)
ops[8](elements["Γₚ₂"][1],f)
ops[8](elements["Γₚ₃"][1],f)
ops[8](elements["Γₚ₄"][1],f)

d = [w(n.x,n.y) for n in nodes]
f .-= k*d

# errk = d'*k*d
# errf = d'*f+216
# err = d'*k*d-d'*f
# d = k\f

# push!(nodes,:d=>d)
# set𝓖!(elements["Ω"],:TriGI16,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
# set∇̂³𝝭!(elements["Ω"])
# prescribe!(elements["Ω"],:u=>(x,y,z)->w(x,y))
# prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->w₁(x,y))
# prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->w₂(x,y))
# prescribe!(elements["Ω"],:∂²u∂x²=>(x,y,z)->w₁₁(x,y))
# prescribe!(elements["Ω"],:∂²u∂x∂y=>(x,y,z)->w₁₂(x,y))
# prescribe!(elements["Ω"],:∂²u∂y²=>(x,y,z)->w₂₂(x,y))
# prescribe!(elements["Ω"],:∂³u∂x³=>(x,y,z)->w₁₁₁(x,y))
# prescribe!(elements["Ω"],:∂³u∂x²∂y=>(x,y,z)->w₁₁₂(x,y))
# prescribe!(elements["Ω"],:∂³u∂x∂y²=>(x,y,z)->w₁₂₂(x,y))
# prescribe!(elements["Ω"],:∂³u∂y³=>(x,y,z)->w₂₂₂(x,y))
# h3,h2,h1,l2 = ops[9](elements["Ω"])
