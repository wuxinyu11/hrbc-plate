
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

# type = (Node,:Quadratic2D,:□,:QuinticSpline)
type = (SNode,:Cubic2D,:□,:QuinticSpline)
s = 3.5/20*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
n = 2

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 3,γ = 5)
elements["Ω"] = ReproducingKernel{type...,:Tri3}(elements["Ω"],sp)
elements["Ωˢ"] = ReproducingKernel{type...,:Tri3}(elements["Ω"])
elements["Γₚ₁"] = ReproducingKernel{type...,:Poi1}(elements["Γₚ₁"],sp)
elements["Γₚ₂"] = ReproducingKernel{type...,:Poi1}(elements["Γₚ₂"],sp)
elements["Γₚ₃"] = ReproducingKernel{type...,:Poi1}(elements["Γₚ₃"],sp)
elements["Γₚ₄"] = ReproducingKernel{type...,:Poi1}(elements["Γₚ₄"],sp)
elements["Γ₁"] = ReproducingKernel{type...,:Seg2}(elements["Γ₁"],sp)
elements["Γ₂"] = ReproducingKernel{type...,:Seg2}(elements["Γ₂"],sp)
elements["Γ₃"] = ReproducingKernel{type...,:Seg2}(elements["Γ₃"],sp)
elements["Γ₄"] = ReproducingKernel{type...,:Seg2}(elements["Γ₄"],sp)

set𝓖!(elements["Ω"],:TriRK6,:∂1,:∂x,:∂y,:∂z)
set𝓖!(elements["Ωˢ"],:TriGI3,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z,:∂x∂z,:∂y∂z,:∂z²)
set𝓖!(elements["Γₚ₁"],:PoiGI1,:∂1)
set𝓖!(elements["Γₚ₂"],:PoiGI1,:∂1)
set𝓖!(elements["Γₚ₃"],:PoiGI1,:∂1)
set𝓖!(elements["Γₚ₄"],:PoiGI1,:∂1)
set𝓖!(elements["Γ₁"],:SegRK3,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z,:∂x∂z,:∂y∂z,:∂z²)
set𝓖!(elements["Γ₂"],:SegRK3,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z,:∂x∂z,:∂y∂z,:∂z²)
set𝓖!(elements["Γ₃"],:SegRK3,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z,:∂x∂z,:∂y∂z,:∂z²)
set𝓖!(elements["Γ₄"],:SegRK3,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂z,:∂x∂z,:∂y∂z,:∂z²)

set∇𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ωˢ"],elements["Ω"])
set𝝭!(elements["Γₚ₁"])
set𝝭!(elements["Γₚ₂"])
set𝝭!(elements["Γₚ₃"])
set𝝭!(elements["Γₚ₄"])
set∇²𝝭!(elements["Γ₁"])
set∇²𝝭!(elements["Γ₂"])
set∇²𝝭!(elements["Γ₃"])
set∇²𝝭!(elements["Γ₄"])

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
M₁₁(x,y) = D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = D*(1-ν)*w₁₂(x,y)
prescribe!(elements["Ω"],:q,(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
prescribe!(elements["Γₚ₁"],:ΔM,(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₂"],:ΔM,(x,y,z)->-2*M₁₂(x,y))
prescribe!(elements["Γₚ₃"],:ΔM,(x,y,z)->2*M₁₂(x,y))
prescribe!(elements["Γₚ₄"],:ΔM,(x,y,z)->-2*M₁₂(x,y))
prescribe!(elements["Γ₁"],:M,(x,y,z)->M₂₂(x,y))
prescribe!(elements["Γ₂"],:M,(x,y,z)->M₁₁(x,y))
prescribe!(elements["Γ₃"],:M,(x,y,z)->M₂₂(x,y))
prescribe!(elements["Γ₄"],:M,(x,y,z)->M₁₁(x,y))
prescribe!(elements["Γ₁"],:V,(x,y,z)->D*(-(2-ν)*w₁₁₂(x,y)-w₂₂₂(x,y)))
prescribe!(elements["Γ₂"],:V,(x,y,z)->D*(w₁₁₁(x,y)+(2-ν)*w₁₂₂(x,y)))
prescribe!(elements["Γ₃"],:V,(x,y,z)->D*((2-ν)*w₁₁₂(x,y)+w₂₂₂(x,y)))
prescribe!(elements["Γ₄"],:V,(x,y,z)->D*(-w₁₁₁(x,y)-(2-ν)*w₁₂₂(x,y)))
# prescribe!(elements["Γ₁"],:M₁,(x,y,z)->-M₁₂(x,y))
# prescribe!(elements["Γ₁"],:M₂,(x,y,z)->-M₂₂(x,y))
# prescribe!(elements["Γ₂"],:M₁,(x,y,z)->M₁₁(x,y))
# prescribe!(elements["Γ₂"],:M₂,(x,y,z)->M₁₂(x,y))
# prescribe!(elements["Γ₃"],:M₁,(x,y,z)->M₁₂(x,y))
# prescribe!(elements["Γ₃"],:M₂,(x,y,z)->M₂₂(x,y))
# prescribe!(elements["Γ₄"],:M₁,(x,y,z)->-M₁₁(x,y))
# prescribe!(elements["Γ₄"],:M₂,(x,y,z)->-M₁₂(x,y))

coefficient = (:D=>1.0,:α=>1e3,:ν=>ν)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫MθdΓ,coefficient...),
       Operator(:∫VwdΓ,coefficient...),
       Operator(:H₃),
       Operator(:∫∇𝑛vθdΓ,:α=>1e7),
       Operator(:∫vgdΓ,:α=>1e7),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫∇wMdΓ,coefficient...),
       Operator(:wΔMₙₛ,coefficient...)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ωˢ"],k)
ops[2](elements["Ω"],f)
# ops[3](elements["Γᶿ₁"],k,f)
# ops[3](elements["Γᶿ₂"],k,f)
# ops[3](elements["Γᶿ₃"],k,f)
# ops[3](elements["Γᶿ₄"],k,f)
# ops[8](elements["Γᶿ₁"],f)
# ops[8](elements["Γᶿ₂"],f)
# ops[8](elements["Γᶿ₃"],f)
# ops[8](elements["Γᶿ₄"],f)
# ops[4](elements["Γᵍ"],k,f)

# ops[6](elements["Γᶿ₁"],k,f)
# ops[6](elements["Γᶿ₂"],k,f)
# ops[6](elements["Γᶿ₃"],k,f)
# ops[6](elements["Γᶿ₄"],k,f)
# ops[7](elements["Γᵍ"],k,f)

ops[8](elements["Γ₁"],f)
ops[8](elements["Γ₂"],f)
ops[8](elements["Γ₃"],f)
ops[8](elements["Γ₄"],f)

ops[9](elements["Γ₁"],f)
ops[9](elements["Γ₂"],f)
ops[9](elements["Γ₃"],f)
ops[9](elements["Γ₄"],f)

# ops[10](elements["Γ₁"],f)
# ops[10](elements["Γ₂"],f)
# ops[10](elements["Γ₃"],f)
# ops[10](elements["Γ₄"],f)

ops[11](elements["Γₚ₁"],f)
ops[11](elements["Γₚ₂"],f)
ops[11](elements["Γₚ₃"],f)
ops[11](elements["Γₚ₄"],f)

d = [w(nodes[:x][i],nodes[:y][i]) for i in 1:length(nodes[:x])]

f .-= k*d
