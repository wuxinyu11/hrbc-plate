
using YAML, ApproxOperator, XLSX

# 𝒑 = "cubic"
𝒑 = "quartic"
ndiv = 32
config = YAML.load_file("./yml/hollow_cylinder_gauss_penalty_alpha_"*𝒑*".yml")
elements,nodes = importmsh("./msh/hollow_cylinder_"*string(ndiv)*".msh", config)
nₚ = length(nodes)

s = zeros(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
for node in nodes
    x = node.x
    y = node.y
    r = (x^2+y^2)^0.5
    sᵢ = 4.05*r*π/4/ndiv
    node.s₁ = sᵢ
    node.s₂ = sᵢ
    node.s₃ = sᵢ
end

# set∇₂𝝭!(elements["Ω"])
set𝝭!(elements["Ω̄"])
set∇²₂𝝭!(elements["Ω"])
set∇₂𝝭!(elements["Γᴹ"])
set𝝭!(elements["Γⱽ"])
set𝝭!(elements["Γᵍ"])
set∇₂𝝭!(elements["Γᶿ"])
set𝝭!(elements["Γᴾ"])

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
set𝒏!(elements["Γᵍ"])
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->w(x,y))
set𝒏!(elements["Γᶿ"])
prescribe!(elements["Γᶿ"],:θ=>(x,y,z,n₁,n₂)->w₁(x,y)*n₁+w₂(x,y)*n₂)
set𝒏!(elements["Γᴹ"])
prescribe!(elements["Γᴹ"],:M=>(x,y,z,n₁,n₂)->M₁₁(x,y)*n₁*n₁+2*M₁₂(x,y)*n₁*n₂+M₂₂(x,y)*n₂*n₂)
set𝒏!(elements["Γⱽ"])
prescribe!(elements["Γⱽ"],:V=>(x,y,z,n₁,n₂)->Vₙ(x,y,n₁,n₂))
prescribe!(elements["Γᴾ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Ω̄"],:u=>(x,y,z)->w(x,y))

coefficient = (:D=>D,:ν=>ν)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫vgdΓ,coefficient...,:α=>1e8),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫∇𝑛vθdΓ,coefficient...,:α=>1e5),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:L₂)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
d = zeros(nₚ)
push!(nodes,:d=>d)
       
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[4](elements["Γⱽ"],f)
ops[6](elements["Γᴹ"],f)
ops[3](elements["Γᴾ"],k,f)

αs = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16]
for (i,α) in enumerate(αs)
    for (j,β) in enumerate(αs)
# for (i,α) in enumerate([1e12,1e13])
#     for (j,β) in enumerate([1e5,1e5,1e7])
        println(i,j)

        fill!(kα,0.0)
        fill!(fα,0.0)


        opv = Operator(:∫vgdΓ,coefficient...,:α=>α)
        opm = Operator(:∫∇𝑛vθdΓ,coefficient...,:α=>β)

        opv(elements["Γᵍ"],kα,fα)
        opm(elements["Γᶿ"],kα,fα)

        d .= (k+kα)\(f+fα)

        l2 = ops[8](elements["Ω̄"])
        println(log10(l2))
    
        XLSX.openxlsx("./xlsx/alpha.xlsx", mode="rw") do xf
            𝐿₂_row = Char(64+j)*string(i)
            xf[5][𝐿₂_row] = log10(l2)
        end
    end
end