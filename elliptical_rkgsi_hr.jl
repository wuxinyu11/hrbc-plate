
using YAML, ApproxOperator, XLSX, TimerOutputs

to = TimerOutput()
@timeit to "Total Time" begin
@timeit to "searching" begin

𝒑 = "cubic"
# 𝒑 = "quartic"
config = YAML.load_file("./yml/elliptical_rkgsi_hr_"*𝒑*".yml")

ndiv = 8
elements, nodes = importmsh("./msh/elliptical_"*string(ndiv)*".msh", config)
end

nₚ = length(nodes)
nₑ = length(elements["Ω"])

# s = 3.1*π/2/ndiv * ones(nₚ)
  s = 3.5/ndiv * ones(nₚ)
  push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

# s = zeros(nₚ)
# push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
# for node in nodes
#     x = node.x
#     y = node.y
#     #  r = (x^2+y^2)^0.5
#     # quartic, ndiv = 32, s = 4.05
#     # sᵢ = 3.05*r*π/4/ndiv
#     sᵢ = 3.5/ndiv

#     node.s₁ = sᵢ
#     node.s₂ = sᵢ
#     node.s₃ = sᵢ
# end

set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
# set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
# set_memory_𝗠!(elements["Γᴹ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
# set_memory_𝗠!(elements["Γⱽ"],:𝝭)
set_memory_𝗠!(elements["Γᶿ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²)
# set_memory_𝗠!(elements["Γᴹ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²)
set_memory_𝗠!(elements["Γᴾ"],:𝝭,:∇̃²)
set_memory_𝗠!(elements["Γₚ₁"],:𝝭,:∇̃²)

elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]
elements["Ω∩Γᶿ"] = elements["Ω"]∩elements["Γᶿ"]
# elements["Ω∩Γᴹ"] = elements["Ω"]∩elements["Γᴹ"]
elements["Ω∩Γᴾ"] = elements["Ω"]∩elements["Γᴾ"]
elements["Ω∩Γₚ₁"] = elements["Ω"]∩elements["Γₚ₁"]
elements["Γ"] = elements["Γᵍ"]∪elements["Γᶿ"]
elements["Γₚ"] = elements["Γₚ₁"]∪elements["Γᴾ"]
elements["Γ∩Γₚ"] = elements["Γ"]∩elements["Γₚ"]



@timeit to "shape functions " begin      
set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
# set∇₂𝝭!(elements["Γᴹ"])
set𝝭!(elements["Γⱽ"])
@timeit to "shape functions Γᵍ " begin      
set∇∇̃²𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set𝝭!(elements["Γᵍ"])
set∇̃²𝝭!(elements["Γᶿ"],elements["Ω∩Γᶿ"])
set∇₂𝝭!(elements["Γᶿ"])
# set∇̃²𝝭!(elements["Γᴹ"],elements["Ω∩Γᴹ"])
# set∇₂𝝭!(elements["Γᴹ"])
set∇̃²𝝭!(elements["Γᴾ"],elements["Ω∩Γᴾ"])
set𝝭!(elements["Γᴾ"])
set∇̃²𝝭!(elements["Γₚ₁"],elements["Ω∩Γₚ₁"])
set𝝭!(elements["Γₚ₁"])


set∇∇̄²𝝭!(elements["Γᵍ"],Γᵍ=elements["Γᵍ"],Γᶿ=elements["Γᵍ"],Γᴾ=elements["Γᴾ"])
set∇̄²𝝭!(elements["Γᶿ"],Γᶿ=elements["Γᶿ"],Γᴾ=elements["Γₚ"],)
# set∇̄²𝝭!(elements["Γᴹ"],Γᵍ=elements["Γᵍ"],Γᶿ=elements["Γᴹ"],Γᴾ=elements["Γₚ"],)
set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γᵍ"],Γᶿ=elements["Γᶿ"]∪elements["Γᵍ"],Γᴾ=elements["Γₚ"])
  # set∇̄²𝝭!(elements["Γₚ₁"],Γᵍ=elements["Γₚ₁"],Γᴾ=elements["Γₚ₁"])
end
end

n = 1
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
 a = 2^(0.5) 
 b = 1.0
 p₀= 1.0
 D = 1.0
 ν = 0.3
 C=p₀*(a^4)*(b^4)/(8*D*(3*a^4+3*b^4+2*a^2*b^2))

# w(x,y) = (x^2/a^2+y^2/b^2-1)^2*C
# w₁(x,y) = 4*x/a^2*(x^2/a^2+y^2/b^2-1)*C
# w₂(x,y) = 4*y/b^2*(x^2/a^2+y^2/b^2-1)*C
# w₁₁(x,y) = 8*x^2/a^4*C+4*C/a^2*(x^2/a^2+y^2/b^2-1)
# w₂₂(x,y) = 8*y^2/b^4*C+4*C/b^2*(x^2/a^2+y^2/b^2-1)
# w₁₂(x,y) = 8*x*y/a^2/b^2*C
# w₁₁₁(x,y) = 24*x*C/a^4
# w₁₁₂(x,y) = 8*y*C/a^2/b^2
# w₁₂₂(x,y) = 8*x*C/a^2/b^2
# w₂₂₂(x,y) = 24*y*C/b^4
# w₁₁₁₁(x,y) = 24*C/a^4
# w₁₁₂₂(x,y) = 8*C/a^2/b^2
# w₂₂₂₂(x,y) = 24*C/b^4


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
# prescribe!(elements["Γᴹ"],:θ=>(x,y,z,n₁,n₂)->w₁(x,y)*n₁+w₂(x,y)*n₂)
prescribe!(elements["Γᶿ"],:θ=>(x,y,z,n₁,n₂)->w₁(x,y)*n₁+w₂(x,y)*n₂)
# set𝒏!(elements["Γᴹ"])
# prescribe!(elements["Γᴹ"],:M=>(x,y,z,n₁,n₂)->M₁₁(x,y)*n₁*n₁+2*M₁₂(x,y)*n₁*n₂+M₂₂(x,y)*n₂*n₂)
set𝒏!(elements["Γⱽ"])
 prescribe!(elements["Γⱽ"],:V=>(x,y,z,n₁,n₂)->Vₙ(x,y,n₁,n₂))
prescribe!(elements["Γᴾ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γₚ₁"],:g=>(x,y,z)->w(x,y))
# prescribe!(elements["Γₚ₁"],:ΔM=>(x,y,z)->0*M₁₂(x,y))


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

@timeit to "assembly" begin       
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
ops[4](elements["Γⱽ"],f)
# ops[6](elements["Γᴹ"],f)
@timeit to "assembly Γᵍ" begin       

ops[3](elements["Γᵍ"],k,f)
# ops[5](elements["Γᴹ"],k,f)
 ops[5](elements["Γᵍ"],k,f)
ops[5](elements["Γᶿ"],k,f)
ops[7](elements["Γₚ"],k,f)
# ops[7](elements["Γₚ₁"],k,f)

end
end

# F = eigen(k)
# F.values[1]
d = k\f
end

# # d = [w(n.x,n.y) for n in nodes]
# # f .-= k*d

push!(nodes,:d=>d)
set𝓖!(elements["Ω"],:TriGI16,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
set∇̂³𝝭!(elements["Ω"])
# set_memory_𝗠!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
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
show(to)
h3,h2,h1,l2 = ops[9](elements["Ω"])



  # index = [10,20,40,80]
#    index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/circular_"*𝒑*".xlsx", mode="rw") do xf
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