

using YAML, ApproxOperator, XLSX, TimerOutputs

to = TimerOutput()
@timeit to "Total Time" begin
@timeit to "searching" begin

ndiv = 10
𝒑 = "cubic"
# 𝒑 = "quartic"
config = YAML.load_file("./yml/triangle_rkgsi_nitsche_"*𝒑*".yml")
elements,nodes = importmsh("./msh/triangle_"*string(ndiv)*".msh", config)
end

nₚ = length(nodes)
s = 5*10/ndiv*ones(nₚ)
# s = 4.5*10/ndiv*ones(nₚ)
#push!(nodes,:s₁=>3^(0.5)/2 .*s,:s₂=>s,:s₃=>s)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)


set_memory_𝗠!(elements["Ω̃"],:∇̃²)

@timeit to "shape functions " begin      
set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
@timeit to "shape functions Γᵍ " begin      
set∇³𝝭!(elements["Γ₁"])
set∇³𝝭!(elements["Γ₂"])
set∇³𝝭!(elements["Γ₃"])

set∇²₂𝝭!(elements["Γₚ₁"])
set∇²₂𝝭!(elements["Γₚ₂"])
set∇²₂𝝭!(elements["Γₚ₃"])
end
end

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

set𝒏!(elements["Γ₁"])
set𝒏!(elements["Γ₂"])
set𝒏!(elements["Γ₃"])
prescribe!(elements["Ω"],:q=>(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))

prescribe!(elements["Γₚ₁"],:Δn₁s₁=>(x,y,z)->-3^(0.5)/2)
prescribe!(elements["Γₚ₁"],:Δn₂s₂=>(x,y,z)->3^(0.5)/2)
prescribe!(elements["Γₚ₂"],:Δn₁s₁=>(x,y,z)->3^(0.5)/4)
prescribe!(elements["Γₚ₂"],:Δn₁s₂n₂s₁=>(x,y,z)->-3/2)
prescribe!(elements["Γₚ₂"],:Δn₂s₂=>(x,y,z)->-3^(0.5)/4)
prescribe!(elements["Γₚ₃"],:Δn₁s₁=>(x,y,z)->3^(0.5)/4)
prescribe!(elements["Γₚ₃"],:Δn₁s₂n₂s₁=>(x,y,z)->3/2)
prescribe!(elements["Γₚ₃"],:Δn₂s₂=>(x,y,z)->-3^(0.5)/4)


# cubic
# α = 1e3*ndiv^2
# quartic
# ndiv = 10, 20, 40, 80, α = 1e3*ndiv^2
# ndiv = 80, α = 1e2*ndiv^3

    #    quartic-0909
       # ndiv = 10, α = 1e4
       # ndiv = 20, α = 1e6
       # ndiv = 40, α = 1e6
       # ndiv = 80, α = 5e8

       #    cubic-0909
       # ndiv = 10, α = 1e4
       # ndiv = 20, α = 1e5
       # ndiv = 40, α = 1e6
       # ndiv = 80, α = 1e7

coefficient = (:D=>D,:ν=>ν)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫VgdΓ,coefficient...,:α=>1e4),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫MₙₙθdΓ,coefficient...,:α=>1e3*ndiv),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔMₙₛg,coefficient...,:α=>1e1),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

@timeit to "assembly" begin       
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)

@timeit to "assembly Γᵍ" begin       
ops[3](elements["Γ₁"],k,f)
ops[3](elements["Γ₂"],k,f)
ops[3](elements["Γ₃"],k,f)


#ops[5](elements["Γ₁"],k,f)
#ops[5](elements["Γ₂"],k,f)
#ops[5](elements["Γ₃"],k,f)

# ops[6](elements["Γ₁"],f)
# ops[6](elements["Γ₂"],f)
# ops[6](elements["Γ₃"],f)
# ops[6](elements["Γ₄"],f)

ops[7](elements["Γₚ₁"],k,f)
ops[7](elements["Γₚ₂"],k,f)
ops[7](elements["Γₚ₃"],k,f)

# ops[8](elements["Γₚ₁"],f)
# ops[8](elements["Γₚ₂"],f)
# ops[8](elements["Γₚ₃"],f)
# ops[8](elements["Γₚ₄"],f)

end
end

d = k\f
end

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
show(to)

index = [10,20,40,80]
XLSX.openxlsx("./xlsx/triangle_"*𝒑*".xlsx", mode="rw") do xf
    row = "F"
    𝐿₂ = xf[2]
    𝐻₁ = xf[3]
    𝐻₂ = xf[4]
    𝐻₃ = xf[5]
    ind = findfirst(n->n==ndiv,index)+1
    row = row*string(ind)
    𝐿₂[row] = log10(l2)
    𝐻₁[row] = log10(h1)
    𝐻₂[row] = log10(h2)
    𝐻₃[row] = log10(h3)
end