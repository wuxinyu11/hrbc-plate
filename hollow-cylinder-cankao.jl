
using Revise, YAML, ApproxOperator

ndiv = 10
config = YAML.load_file("./yml/patch_test_rkgsi_nitsche.yml")
elements,nodes = importmsh("./msh/patchtest_"*string(ndiv)*".msh", config)
nₚ = length(nodes)
s = 3.5/ndiv*ones(nₚ)
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

w(x,y)=4/(3*(1-ν))*log((x^2+y^2)^(1/2)/2)-1/(3*(1+ν))*(x^2+y^2-4)
w₁(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-1)*2*x-1/(3*(1+ν))*2*x
w₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-1)*2*y-1/(3*(1+ν))*2*y
w₁₁(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-2)*4*x^2+8/(3*(1-ν))*(x^2+y^2)^(-1)-1/(3*(1+ν))*2
w₁₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-2)*4*y*x
w₂₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-2)*4*y^2+8/(3*(1-ν))*(x^2+y^2)^(-1)-1/(3*(1+ν))*2
w₁₁₁(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-3)*8*x^3+48*x/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-3)*8*y*x^2+8/(3*(1-ν))*(x^2+y^2)^(-2)*2*y
w₁₂₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-3)*8*x*y^2+8/(3*(1-ν))*(x^2+y^2)^(-2)*2*x
w₂₂₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-3)*8*y^3+48*y/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₁₁(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-4)*16*x^4+192*x/(3*(1-ν))*(x^2+y^2)^(-3)+48/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₂₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-4)*16*x^2*y^2+4/(3*(1-ν))*(x^2+y^2)^(-3)*8*x^3+8/(3*(1-ν))*(x^2+y^2)^(-3)*4*y^2+16/(3*(1-ν))*(x^2+y^2)^(-2)
w₂₂₂₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-4)*16*y^4+192*y/(3*(1-ν))*(x^2+y^2)^(-3)+48/(3*(1-ν))*(x^2+y^2)^(-2)



#w(x,y) = 4/3/(1-ν)*ln((x^2+y^2)/2)-1/3/*(1+ν)*(x^2+y^2-4)
#w₁(x,y) = 4/3/(1-ν)*(x^2+y^2)^(-1)*2*x-2*x/3/(1+ν)
#w₂(x,y) = 4/3/(1-ν)*(x^2+y^2)^(-1)*2*y-2*y/3/(1+ν)
#w₁₁(x,y) = 4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*x*2*x+2*4/3/(1-ν)*(x^2+y^2)^(-1)-2/3/(1+ν)
#w₂₂(x,y) = 4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*y*2*y+2*4/3/(1-ν)*(x^2+y^2)^(-1)-2/3/(1+ν)
#w₁₂(x,y) = 4/3/(1-ν)*(-1)(x^2+y^2)^(-2)*2*y*2*x
#w₁₁₁(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*2*x*2*x*2*x+8*x*4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*x
#w₁₁₂(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*2*y*4*x^2+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*y
#w₁₂₂(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*2*x*4*y^2+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*x
#w₂₂₂(x,y) = 4/3/(1-ν)*(-1)*(-2)*(x^2+y^2)^(-3)*8*y^3+4/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*8*y+8/3/(1-ν)*(-1)*(x^2+y^2)^(-2)*2*y
#w₁₁₁₁(x,y) = 4/3/(1-ν)*(-6)*(x^2+y^2)^(-4)*16*x^4+4/3/(1-ν)*2*(x^2+y^2)^(-3)*24*x^2-32/3/(1-ν)*(x^2+y^2)^(-2)+8/3*16*x^2*(x^2+y^2)^(-3)+16/3/(1-ν)*(x^2+y^2)^(-3)*4*x^2-16/3/(1-ν)*(x^2+y^2)^(-2)
#w₁₁₂₂(x,y) = 4/3/(1-ν)*(-6)*(x^2+y^2)^(-4)*4*y^2*4*x^2+8/3*(x^2+y^2)^(-3)*8*x^2+16/3/(1-ν)*(x^2+y^2)^(-3)*4*y^2-16/3/(1-ν)*(x^2+y^2)^(-2)
#w₂₂₂₂(x,y) = 4/3/(1-ν)*(-6)*(x^2+y^2)^(-4)*16*y^2+3/4/(1-ν)*8*(x^2+y^2)^(-3)*24*y^2*16*4+4/3/(1-ν)*16*2*(x^2+y^2)^(-2)
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

prescribe!(elements["Γ₁"],:V=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:V=>(x,y,z)->0.0)

prescribe!(elements["Γ₁"],:θ=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:θ=>(x,y,z)->0.0)

prescribe!(elements["Γ₂"],:M=>(x,y,z)->1)
prescribe!(elements["Γ₄"],:M=>(x,y,z)->2)

prescribe!(elements["Γₚ₁"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₂"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₃"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₄"],:g=>(x,y,z)->w(x,y))

prescribe!(elements["Γₚ₁"],:Δn₁s₁=>(x,y,z)->1/2)
prescribe!(elements["Γₚ₁"],:Δn₁s₂n₂s₁=>(x,y,z)->-1)
prescribe!(elements["Γₚ₁"],:Δn₂s₂=>(x,y,z)->-1/2)
prescribe!(elements["Γₚ₂"],:Δn₁s₁=>(x,y,z)->-1/2)
prescribe!(elements["Γₚ₂"],:Δn₁s₂n₂s₁=>(x,y,z)->-1)
prescribe!(elements["Γₚ₂"],:Δn₂s₂=>(x,y,z)->1/2)
prescribe!(elements["Γₚ₃"],:Δn₁s₁=>(x,y,z)->1/2)
prescribe!(elements["Γₚ₃"],:Δn₁s₂n₂s₁=>(x,y,z)->1)
prescribe!(elements["Γₚ₃"],:Δn₂s₂=>(x,y,z)->-1/2)
rescribe!(elements["Γₚ₄"],:Δn₁s₁=>(x,y,z)->-1/2)
prescribe!(elements["Γₚ₄"],:Δn₁s₂n₂s₁=>(x,y,z)->1)
prescribe!(elements["Γₚ₄"],:Δn₂s₂=>(x,y,z)->1/2)

prescribe!(elements["Γₚ₁"],:ΔM=>(x,y,z)->-1*M₁₂(x,y))
prescribe!(elements["Γₚ₂"],:ΔM=>(x,y,z)->-1*M₁₂(x,y))
prescribe!(elements["Γₚ₃"],:ΔM=>(x,y,z)->1*M₁₂(x,y))
prescribe!(elements["Γₚ₄"],:ΔM=>(x,y,z)->1*M₁₂(x,y))

coefficient = (:D=>D,:ν=>ν)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
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

#ops[3](elements["Γ₁"],k,f)
ops[3](elements["Γ₂"],k,f)
#ops[3](elements["Γ₃"],k,f)
ops[3](elements["Γ₄"],k,f)

 ops[4](elements["Γ₁"],f)
# ops[4](elements["Γ₂"],f)
 ops[4](elements["Γ₃"],f)
# ops[4](elements["Γ₄"],f)

ops[5](elements["Γ₁"],k,f)
#ops[5](elements["Γ₂"],k,f)
ops[5](elements["Γ₃"],k,f)
#ops[5](elements["Γ₄"],k,f)

# ops[6](elements["Γ₁"],f)
 ops[6](elements["Γ₂"],f)
# ops[6](elements["Γ₃"],f)
 ops[6](elements["Γ₄"],f)

ops[7](elements["Γₚ₁"],k,f)
ops[7](elements["Γₚ₂"],k,f)
ops[7](elements["Γₚ₃"],k,f)
ops[7](elements["Γₚ₄"],k,f)
# ops[8](elements["Γₚ₁"],f)
# ops[8](elements["Γₚ₂"],f)
# ops[8](elements["Γₚ₃"],f)
# ops[8](elements["Γₚ₄"],f)

# d = [w(n.x,n.y) for n in nodes]
# f .-= k*d

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