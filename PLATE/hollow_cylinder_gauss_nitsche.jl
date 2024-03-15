
using YAML, ApproxOperator, XLSX, TimerOutputs

to = TimerOutput()
@timeit to "Total Time" begin
@timeit to "searching" begin

# 𝒑 = "cubic"
𝒑 = "quartic"
ndiv = 64
config = YAML.load_file("./yml/hollow_cylinder_gauss_nitsche_"*𝒑*".yml")
elements,nodes = importmsh("./msh/hollow_cylinder_"*string(ndiv)*".msh", config)
nₚ = length(nodes)
end

# s = 5/ndiv*ones(nₚ)
s = zeros(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
for node in nodes
    x = node.x
    y = node.y
    r = (x^2+y^2)^0.5
    sᵢ = 4.1*r*π/4/ndiv
    node.s₁ = sᵢ
    node.s₂ = sᵢ
    node.s₃ = sᵢ
end

@timeit to "shape functions " begin      
set∇²₂𝝭!(elements["Ω"])
set∇₂𝝭!(elements["Γᴹ"])
set𝝭!(elements["Γⱽ"])
# set∇³𝝭!(elements["Γᴹ"])
# set∇³𝝭!(elements["Γⱽ"])
@timeit to "shape functions Γᵍ " begin      
set∇³𝝭!(elements["Γᵍ"])
set∇²₂𝝭!(elements["Γᶿ"])
# set∇³𝝭!(elements["Γᶿ"])
set∇²₂𝝭!(elements["Γᴾ"])

end
end

w(x,y)=4/(3*(1-ν))*log((x^2+y^2)^(1/2)/2)-1/(3*(1+ν))*(x^2+y^2-4)
w₁(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-1)*x-2*x/(3*(1+ν))
w₂(x,y)=4/(3*(1-ν))*(x^2+y^2)^(-1)*y-2*y/(3*(1+ν))
w₁₁(x,y)=-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*x^2+4/(3*(1-ν))*(x^2+y^2)^(-1)-1/(3*(1+ν))*2
w₁₂(x,y)=-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*y*x
w₂₂(x,y)=-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*y^2+4/(3*(1-ν))*(x^2+y^2)^(-1)-1/(3*(1+ν))*2
w₁₁₁(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*x^3-24*x/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₂(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*y*x^2-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*y
w₁₂₂(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*x*y^2-4/(3*(1-ν))*(x^2+y^2)^(-2)*2*x
w₂₂₂(x,y)=8/(3*(1-ν))*(x^2+y^2)^(-3)*4*y^3-24*y/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₁₁(x,y)=-24/(3*(1-ν))*(x^2+y^2)^(-4)*8*x^4+192*x^2/(3*(1-ν))*(x^2+y^2)^(-3)-24/(3*(1-ν))*(x^2+y^2)^(-2)
w₁₁₂₂(x,y)=-24/(3*(1-ν))*(x^2+y^2)^(-4)*8*x^2*y^2+32/(3*(1-ν))*(x^2+y^2)^(-3)*x^2+32/(3*(1-ν))*(x^2+y^2)^(-3)*y^2-8/(3*(1-ν))*(x^2+y^2)^(-2)
w₂₂₂₂(x,y)=-24/(3*(1-ν))*(x^2+y^2)^(-4)*8*y^4+192*y^2/(3*(1-ν))*(x^2+y^2)^(-3)-24/(3*(1-ν))*(x^2+y^2)^(-2)

D =  1.0
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
set𝒏!(elements["Γᵍ"])
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->w(x,y))
set𝒏!(elements["Γᶿ"])
prescribe!(elements["Γᶿ"],:θ=>(x,y,z,n₁,n₂)->w₁(x,y)*n₁+w₂(x,y)*n₂)
set𝒏!(elements["Γᴹ"])
prescribe!(elements["Γᴹ"],:M=>(x,y,z,n₁,n₂)->M₁₁(x,y)*n₁*n₁+2*M₁₂(x,y)*n₁*n₂+M₂₂(x,y)*n₂*n₂)
set𝒏!(elements["Γⱽ"])
prescribe!(elements["Γⱽ"],:V=>(x,y,z,n₁,n₂)->Vₙ(x,y,n₁,n₂))
prescribe!(elements["Γᴾ"],:g=>(x,y,z)->w(x,y))

 coefficient = (:D=>D,:ν=>ν)
 ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
        Operator(:∫wqdΩ,coefficient...),
        # cubic
       # ndiv = 10, α = 1e4*ndiv^3
       # ndiv = 20, α = 1e3*ndiv^3
       # ndiv = 40, α = 1e3*ndiv^3
       # ndiv = 80, α = 1e3*ndiv^4
    #    quartic
       # ndiv = 8, α = 1e4*ndiv^3
       # ndiv = 16, α = 1e6*ndiv^3
       # ndiv = 32, α = 1e4*ndiv^3
       # ndiv = 64, α = 1e4*ndiv^3
        Operator(:∫VgdΓ,coefficient...,:α=>1e8),
        Operator(:∫wVdΓ,coefficient...),
       # ndiv = 10, α = 1e3*ndiv
       # ndiv = 80, α = 1e2*ndiv
        Operator(:∫MₙₙθdΓ,coefficient...,:α=>1e6),
        Operator(:∫θₙMₙₙdΓ,coefficient...),
       # ndiv = 10, α = 1e1*ndiv^2
       # ndiv = 80, α = 1e0*ndiv^2
        Operator(:ΔMₙₛg,coefficient...,:α=>1e5),
        Operator(:wΔMₙₛ,coefficient...),
        Operator(:H₃),
        Operator(:L₂),]
 
 k = zeros(nₚ,nₚ)
 f = zeros(nₚ)

@timeit to "assembly" begin       
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[4](elements["Γⱽ"],f)
ops[6](elements["Γᴹ"],f)
@timeit to "assembly Γᵍ" begin       
ops[3](elements["Γᵍ"],k,f)
ops[5](elements["Γᶿ"],k,f)
ops[7](elements["Γᴾ"],k,f)
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

# index = [10,20,40,80]
index = [8,16,32,64]
XLSX.openxlsx("./xlsx/hollow_cylinder_"*𝒑*".xlsx", mode="rw") do xf
    # row = "B"
    row = "D"
    𝐿₂ = xf[2]
    𝐻₁ = xf[3]
    𝐻₂ = xf[4]
    𝐻₃ = xf[5]
    ind = findfirst(n->n==ndiv,index)+1
    row = row*string(ind)
    𝐻₁[row] = log10(h1)
    𝐻₂[row] = log10(h2)
    𝐻₃[row] = log10(h3)
    𝐿₂[row] = log10(l2)
end

XLSX.openxlsx("./xlsx/hollow_cylinder_contour.xlsx", mode="rw") do xf
    sheet = xf[1]
    row = "B"
    sheet[row*string(1)] = "gauss-nitsche"
    for (i,node) in enumerate(nodes)
        sheet[row*string(i+1)] = node.d
    end
end