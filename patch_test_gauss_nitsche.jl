
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

# type = (Node,:Quadratic2D,:□,:QuinticSpline)
type = (Node,:Cubic2D,:□,:QuinticSpline)
s = 3.5/20*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
n = 4

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 3,γ = 5)
elements["Ω"] = ReproducingKernel{type...,:Tri3}(elements["Ω"],sp)
elements["Γᵍ"] = ReproducingKernel{type...,:Seg2}(elements["Γᵍ"],sp)
elements["Γᶿ₁"] = ReproducingKernel{type...,:Seg2}(elements["Γᶿ₁"],sp)
elements["Γᶿ₂"] = ReproducingKernel{type...,:Seg2}(elements["Γᶿ₂"],sp)
elements["Γᶿ₃"] = ReproducingKernel{type...,:Seg2}(elements["Γᶿ₃"],sp)
elements["Γᶿ₄"] = ReproducingKernel{type...,:Seg2}(elements["Γᶿ₄"],sp)

set𝓖!(elements["Ω"],:TriGI7,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²)
set𝓖!(elements["Γᵍ"],:SegGI2,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³)
set𝓖!(elements["Γᶿ₁"],:SegGI2,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z)
set𝓖!(elements["Γᶿ₂"],:SegGI2,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z)
set𝓖!(elements["Γᶿ₃"],:SegGI2,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z)
set𝓖!(elements["Γᶿ₄"],:SegGI2,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z)

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
prescribe!(elements["Ω"],:q,(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
prescribe!(elements["Γᵍ"],:w,(x,y,z)->w(x,y))
# prescribe!(elements["Γᵍ"],:g,(x,y,z)->w(x,y))
prescribe!(elements["Γᶿ₁"],:θ,(x,y,z)->-w₂(x,y))
prescribe!(elements["Γᶿ₂"],:θ,(x,y,z)-> w₁(x,y))
prescribe!(elements["Γᶿ₃"],:θ,(x,y,z)-> w₂(x,y))
prescribe!(elements["Γᶿ₄"],:θ,(x,y,z)->-w₁(x,y))

coefficient = (:D=>1.0,:α=>1e3,:ν=>0.3)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫MθdΓ,coefficient...),
       Operator(:∫VwdΓ,coefficient...),
       Operator(:H₃),
       Operator(:∫∇𝑛vθdΓ,:α=>1e7),
       Operator(:∫vgdΓ,:α=>1e7)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γᶿ₁"],k,f)
ops[3](elements["Γᶿ₂"],k,f)
ops[3](elements["Γᶿ₃"],k,f)
ops[3](elements["Γᶿ₄"],k,f)
ops[4](elements["Γᵍ"],k,f)

# ops[6](elements["Γᶿ₁"],k,f)
# ops[6](elements["Γᶿ₂"],k,f)
# ops[6](elements["Γᶿ₃"],k,f)
# ops[6](elements["Γᶿ₄"],k,f)
# ops[7](elements["Γᵍ"],k,f)

d = k\f

push!(nodes,:d=>d)
set𝓖!(elements["Ω"],:TriGI16,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³)
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
