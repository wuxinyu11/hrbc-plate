

using YAML, ApproxOperator, XLSX, TimerOutputs

to = TimerOutput()
@timeit to "Total Time" begin
@timeit to "searching" begin

# 𝒑 = "cubic"
𝒑 = "quartic"
config = YAML.load_file("./yml/hollow_cylinder_rkgsi_hr_"*𝒑*".yml")

ndiv = 64
elements, nodes = importmsh("./msh/hollow_cylinder_"*string(ndiv)*".msh", config)
# elements, nodes = importmsh("./msh/hollow_cylinder_regular_"*string(ndiv)*".msh", config)
end

nₚ = length(nodes)
nₑ = length(elements["Ω"])

# s = 3.1*π/2/ndiv * ones(nₚ)
s = zeros(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
for node in nodes
    x = node.x
    y = node.y
    r = (x^2+y^2)^0.5
    # quartic, ndiv = 8, s = 4.1
    # quartic, ndiv = 16, s = 4.1
    # quartic, ndiv = 32, s = 4.05
    # quartic, ndiv = 64, s = 4.1
    sᵢ = 4.1*r*π/4/ndiv
    node.s₁ = sᵢ
    node.s₂ = sᵢ
    node.s₃ = sᵢ
end

set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γᶿ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²)
set_memory_𝗠!(elements["Γᴾ"],:𝝭,:∇̃²)

elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]
elements["Ω∩Γᶿ"] = elements["Ω"]∩elements["Γᶿ"]
elements["Ω∩Γᴾ"] = elements["Ω"]∩elements["Γᴾ"]
elements["Γ"] = elements["Γᵍ"]∪elements["Γᶿ"]
elements["Γ∩Γᴾ"] = elements["Γ"]∩elements["Γᴾ"]


@timeit to "shape functions " begin      
set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇₂𝝭!(elements["Γᴹ"])
set𝝭!(elements["Γⱽ"])
@timeit to "shape functions Γᵍ " begin      
set∇∇̃²𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set𝝭!(elements["Γᵍ"])
set∇̃²𝝭!(elements["Γᶿ"],elements["Ω∩Γᶿ"])
set∇₂𝝭!(elements["Γᶿ"])
set∇̃²𝝭!(elements["Γᴾ"],elements["Ω∩Γᴾ"])
set𝝭!(elements["Γᴾ"])

set∇∇̄²𝝭!(elements["Γᵍ"],Γᵍ=elements["Γᵍ"],Γᴾ=elements["Γᴾ"])
set∇̄²𝝭!(elements["Γᶿ"],Γᶿ=elements["Γᶿ"],Γᴾ=elements["Γᴾ"])
set∇̄²𝝭!(elements["Γᴾ"],Γᵍ=elements["Γᵍ"],Γᶿ=elements["Γᶿ"],Γᴾ=elements["Γᴾ"])

end
end

# n = 3
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

c1 = 4/D/(1-ν)/3
c2 = -2/D/(1+ν)/3
r²(x,y) = x^2+y^2
w(x,y)=c1*log(r²(x,y)^0.5/2)+c2/2*(x^2+y^2-4)
w₁(x,y)=c1*x/r²(x,y)+c2*x
w₂(x,y)=c1*y/r²(x,y)+c2*y
w₁₁(x,y)=c1*(-2*x^2/r²(x,y)^2+1/r²(x,y))+c2
w₁₂(x,y)=c1*(-2*x*y/r²(x,y)^2)
w₂₂(x,y)=c1*(-2*y^2/r²(x,y)^2+1/r²(x,y))+c2
w₁₁₁(x,y)=c1*(8*x^3/r²(x,y)^3-6*x/r²(x,y)^2)
w₁₁₂(x,y)=c1*(8*x^2*y/r²(x,y)^3-2*y/r²(x,y)^2)
w₁₂₂(x,y)=c1*(8*x*y^2/r²(x,y)^3-2*x/r²(x,y)^2)
w₂₂₂(x,y)=c1*(8*y^3/r²(x,y)^3-6*y/r²(x,y)^2)
w₁₁₁₁(x,y)=c1*(-48*x^4/r²(x,y)^4+48*x^2/r²(x,y)^3-6/r²(x,y)^2)
w₁₁₂₂(x,y)=c1*(-48*x^2*y^2/r²(x,y)^4+8*x^2/r²(x,y)^3+8*y^2/r²(x,y)^3-2/r²(x,y)^2)
w₂₂₂₂(x,y)=c1*(-48*y^4/r²(x,y)^4+48*y^2/r²(x,y)^3-6/r²(x,y)^2)

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

@timeit to "assembly" begin       
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
ops[4](elements["Γⱽ"],f)
ops[6](elements["Γᴹ"],f)
@timeit to "assembly Γᵍ" begin       

ops[3](elements["Γᵍ"],k,f)
ops[5](elements["Γᶿ"],k,f)
ops[7](elements["Γᴾ"],k,f)
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

# index = [10,20,40,80]
index = [8,16,32,64]
XLSX.openxlsx("./xlsx/hollow_cylinder_"*𝒑*".xlsx", mode="rw") do xf
    row = "G"
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
XLSX.openxlsx("./xlsx/hollow_cylinder_contour.xlsx", mode="rw") do xf
    sheet = xf[1]
    row = "E"
    sheet[row*string(1)] = "rkgsi-hr"
    for (i,node) in enumerate(nodes)
        sheet[row*string(i+1)] = node.d
    end
end

# log10(l2)

# set𝓖!(elements["Ω"],:TriRK6,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
# using GLMakie, CairoMakie
# f = Figure()
# ax = Axis(f[1, 1])
# ax.aspect = 1
# hidespines!(ax)
# hidedecorations!(ax)
# vertices = zeros(nₚ,2)
# colors = zeros(nₚ)
# faces = zeros(Int,nₑ,3)

# set∇²₂𝝭!(elements["Ω"])
# for (c,elm) in enumerate(elements["Ω"])
#     𝓒 = elm.𝓒
#     𝓖 = elm.𝓖
#     faces[c,1] = 𝓒[1].𝐼
#     faces[c,2] = 𝓒[2].𝐼
#     faces[c,3] = 𝓒[3].𝐼
#     for (j,ξ) in enumerate(𝓖[1:3])
#         r = (ξ.x^2+ξ.y^2)^0.5
#         n₁ = ξ.x/r
#         n₂ = ξ.y/r
#         B₁₁ = ξ[:∂²𝝭∂x²]
#         B₁₂ = ξ[:∂²𝝭∂x∂y]
#         B₂₂ = ξ[:∂²𝝭∂y²]
#         M₁₁ = 0.0
#         M₁₂ = 0.0
#         M₂₂ = 0.0
#         for (i,xᵢ) in enumerate(𝓒)
#             M₁₁ -= D*(B₁₁[i]+ν*B₂₂[i])*xᵢ.d
#             M₂₂ -= D*(ν*B₁₁[i]+B₂₂[i])*xᵢ.d
#             M₁₂ -= D*(1-ν)*B₁₂[i]*xᵢ.d
#         end
#         vertices[𝓒[j].𝐼,1] = ξ.x
#         vertices[𝓒[j].𝐼,2] = ξ.y
#         colors[𝓒[j].𝐼] = M₁₁*n₁*n₁+2*M₁₂*n₁*n₂+M₂₂*n₂*n₂ 
#     end
# end
# mesh!(
#     vertices,
#     faces,
#     color = colors,
#     colormap = :balance,
#     shading = false
#     )
# contour!(vertices[:,1],vertices[:,2],colors,levels=1:0.1:2,color=:black)
# contour!()

# f
# save("./figure/rkgsi_hr.png",f)
