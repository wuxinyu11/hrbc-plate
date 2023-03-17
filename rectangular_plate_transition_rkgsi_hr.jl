
using YAML, ApproxOperator, XLSX,LinearAlgebra

ndiv = 10
𝒑 = "cubic"
# 𝒑 = "quartic"
config = YAML.load_file("./yml/rectangular_rkgsi_hr_"*𝒑*".yml")
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)

# naturall bc
sp = ApproxOperator.RegularGrid(nodes,n=2,γ=5)
data = Dict([:x=>(2,[5.0]),:y=>(2,[5.0],:z=>(2,[0.0],:𝑤=>(2,[1.0])))])
ξ = SNode((1,1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Quadratic2D,:□,:QuinticSpline,:Tri3}(𝓒,[ξ])]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)
 
nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
       
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)

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

 
    
set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set𝝭!(elements["Ω"])
 
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



set∇∇̄²𝝭!(elements["Γ₁"],Γᵍ=elements["Γ₁"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₂"],Γᵍ=elements["Γ₂"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"],Γᴾ=elements["Γₚ"])
set∇∇̄²𝝭!(elements["Γ₄"],Γᵍ=elements["Γ₄"],Γᴾ=elements["Γₚ"])
set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γ∩Γₚ"],Γᴾ=elements["Γₚ"])



function w(x,y)
    w_ = 0.0
    max_iter = 5
    for m in 1:max_iter
        for n in 1:max_iter
            w_ += W(x,y,m,n)*η(t,m,n)
        end
    end
    return w_
end
W(x,y,m,n) = 2/a/(ρh)^0.5*sin(m*π*x/a)*sin(n*π*y/a)
η(t,m,n) = 2*F₀/(ω(m,n)^2-Θ^2)/a/(ρh)^0.5*sin(m*π/a)*sin(n*π/a)*(sin(Θ*t)-Θ/ω(m,n)*sin(ω(m,n)*t))

prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₄"],:g=>(x,y,z)->0.0)


coefficient = (:D=>1.0,:ν=>0.3,:ρ=>8000.0,:h=>0.05)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫ρhvwdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫ṼgdΓ,coefficient...),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫M̃ₙₙθdΓ,coefficient...),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔM̃ₙₛg,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
                                   # ----施加边界：集中力
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
f = zeros(nₚ)
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],m)

ops[4](elements["Γ₁"],k,f)
ops[4](elements["Γ₂"],k,f)
ops[4](elements["Γ₃"],k,f)
ops[4](elements["Γ₄"],k,f)
ops[7](elements["Γ₁"],f)
ops[7](elements["Γ₂"],f)
ops[7](elements["Γ₃"],f)
ops[7](elements["Γ₄"],f)
# ops[9](elements["Γₚ"],k,f)


A=eigvals(m,k)

Θ = π
β = 0.0
γ = 0.5
Δt = 0.1
total_time = 1.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
deflection = zeros(length(times))
                #  ---时间从d=0到第n步
v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)
                            # --第几个t
    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->sin(Θ*t))   
                        #    ----把这个点上的集中力设到这个单元上
    f = zeros(nₚ)
    ops[5](elements["Γᵗ"],f)

    a = (m + β*Δt^2*k)\(f-k*d)
                     #       ----这个m是前面的m(kₚ,kₚ)的结果
    # predictor phase
    d .+= Δt*v + Δt^2/2.0*(1.0-2.0*β)*aₙ
    v .+= Δt*(1.0-γ)*aₙ

    # Corrector phase
    d .+= β*Δt^2*a
    v .+= γ*Δt*a

    # cal deflection
    ξ = elements["Γᵗ"].𝓖[1]
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(elements["Γᵗ"].𝓒)
        I = xᵢ.𝐼
        deflection[n] += N[i]*d[I]
    end
end


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
# h3,h2,h1,l2 = ops[10](elements["Ω"])


# index = [10,20,40,80]
# XLSX.openxlsx("./xlsx/rectangular_"*𝒑*".xlsx", mode="rw") do xf
#     row = "G"
#     𝐿₂ = xf[2]
#     𝐻₁ = xf[3]
#     𝐻₂ = xf[4]
#     𝐻₃ = xf[5]
#     ind = findfirst(n->n==ndiv,index)+1
#     row = row*string(ind)
#     𝐿₂[row] = log10(l2)
#     𝐻₁[row] = log10(h1)
#     𝐻₂[row] = log10(h2)
#     𝐻₃[row] = log10(h3)
# end

