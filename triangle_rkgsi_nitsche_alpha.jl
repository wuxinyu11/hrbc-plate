

using YAML, ApproxOperator, XLSX

ndiv = 10
# 𝒑 = "cubic"
𝒑 = "quartic"
config = YAML.load_file("./yml/triangle_rkgsi_nitsche_alpha_"*𝒑*".yml")
elements,nodes = importmsh("./msh/triangle_"*string(ndiv)*".msh", config)

nₚ = getnₚ(elements["Ω"])

s = 5*10/ndiv*ones(nₚ)
#s = 4.5*10/ndiv*ones(nₚ)
#push!(nodes,:s₁=>3^(0.5)/2 .*s,:s₂=>s,:s₃=>s)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)


set_memory_𝗠!(elements["Ω̃"],:∇̃²)

set∇₂𝝭!(elements["Ω"])
set𝝭!(elements["Ω̄"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇³𝝭!(elements["Γ₁"])
set∇³𝝭!(elements["Γ₂"])
set∇³𝝭!(elements["Γ₃"])

set∇²₂𝝭!(elements["Γₚ₁"])
set∇²₂𝝭!(elements["Γₚ₂"])
set∇²₂𝝭!(elements["Γₚ₃"])

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

set𝒏!(elements["Γ₁"])
set𝒏!(elements["Γ₂"])
set𝒏!(elements["Γ₃"])
prescribe!(elements["Γₚ₁"],:Δn₁s₁=>(x,y,z)->-3^(0.5)/2)
prescribe!(elements["Γₚ₁"],:Δn₂s₂=>(x,y,z)->3^(0.5)/2)
prescribe!(elements["Γₚ₂"],:Δn₁s₁=>(x,y,z)->3^(0.5)/4)
prescribe!(elements["Γₚ₂"],:Δn₁s₂n₂s₁=>(x,y,z)->-3/2)
prescribe!(elements["Γₚ₂"],:Δn₂s₂=>(x,y,z)->-3^(0.5)/4)
prescribe!(elements["Γₚ₃"],:Δn₁s₁=>(x,y,z)->3^(0.5)/4)
prescribe!(elements["Γₚ₃"],:Δn₁s₂n₂s₁=>(x,y,z)->3/2)
prescribe!(elements["Γₚ₃"],:Δn₂s₂=>(x,y,z)->-3^(0.5)/4)
prescribe!(elements["Ω̄"],:u=>(x,y,z)->w(x,y))


# cubic
# α = 1e3*ndiv^2
# quartic
# ndiv = 10, 20, 40, 80, α = 1e3*ndiv^2
# ndiv = 80, α = 1e2*ndiv^3
coefficient = (:D=>D,:ν=>ν)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫VgdΓ,coefficient...,:α=>1e3*ndiv^2),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫MₙₙθdΓ,coefficient...,:α=>1e3*ndiv),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔMₙₛg,coefficient...,:α=>1e3),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:L₂)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
d = zeros(nₚ)
push!(nodes,:d=>d)

αs = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16]
for (i,α) in enumerate(αs)
    println(i)

    fill!(k,0.0)
    fill!(f,0.0)
    ops[1](elements["Ω̃"],k)
    ops[2](elements["Ω"],f)

    opv = Operator(:∫VgdΓ,coefficient...,:α=>α)
    opv(elements["Γ₁"],k,f)
    opv(elements["Γ₂"],k,f)
    opv(elements["Γ₃"],k,f)

    ops[7](elements["Γₚ₁"],k,f)
    ops[7](elements["Γₚ₂"],k,f)
    ops[7](elements["Γₚ₃"],k,f)

    d .= k\f

    l2 = ops[9](elements["Ω̄"])
    println(log10(l2))

    XLSX.openxlsx("./xlsx/alpha.xlsx", mode="rw") do xf
        𝐿₂_row = "C"*string(i+1)
        xf_ = 𝒑 == "cubic" ? xf[3] : xf[4]
        xf_[𝐿₂_row] = log10(l2)
    end
end