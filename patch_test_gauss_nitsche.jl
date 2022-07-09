
using Revise, YAML, ApproxOperator

config = YAML.load_file("./yml/patch_test_gauss_nitsche.yml")
elements = importmsh("./msh/patchtest_10.msh",config)
nₚ = getnₚ(elements["Ω"])

s = 3.5/20*ones(nₚ)
push!(elements["Ω"],:s₁=>s,:s₂=>s,:s₃=>s)
set𝓖!(elements["Ω"],:TriGI3)
n = 4

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
prescribe!(elements["Ω"],:q=>(x,y,z)->w₁₁₁₁(x,y)+2*w₁₁₂₂(x,y)+w₂₂₂₂(x,y))
# prescribe!(elements["Γᵍ"],:w,(x,y,z)->w(x,y))
# # prescribe!(elements["Γᵍ"],:g,(x,y,z)->w(x,y))
# prescribe!(elements["Γᶿ₁"],:θ,(x,y,z)->-w₂(x,y))
# prescribe!(elements["Γᶿ₂"],:θ,(x,y,z)-> w₁(x,y))
# prescribe!(elements["Γᶿ₃"],:θ,(x,y,z)-> w₂(x,y))
# prescribe!(elements["Γᶿ₄"],:θ,(x,y,z)->-w₁(x,y))

coefficient = (:D=>1.0,:ν=>0.3)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫VwdΓ,coefficient...,:α=>1e3),
       Operator(:∫MθdΓ,coefficient...,:α=>1e3),
       Operator(:ΔMₙₛg,coefficient...,:α=>1e3),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)

d = k\f

push!(elements["Ω"],:d=>d)
# set𝓖!(elements["Ω"],:TriGI16,:∂1,:∂x,:∂y,:∂x²,:∂x∂y,:∂y²,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³)
# prescribe!(elements["Ω"],:u,(x,y,z)->w(x,y))
# prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->w₁(x,y))
# prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->w₂(x,y))
# prescribe!(elements["Ω"],:∂²u∂x²,(x,y,z)->w₁₁(x,y))
# prescribe!(elements["Ω"],:∂²u∂x∂y,(x,y,z)->w₁₂(x,y))
# prescribe!(elements["Ω"],:∂²u∂y²,(x,y,z)->w₂₂(x,y))
# prescribe!(elements["Ω"],:∂³u∂x³,(x,y,z)->w₁₁₁(x,y))
# prescribe!(elements["Ω"],:∂³u∂x²∂y,(x,y,z)->w₁₁₂(x,y))
# prescribe!(elements["Ω"],:∂³u∂x∂y²,(x,y,z)->w₁₂₂(x,y))
# prescribe!(elements["Ω"],:∂³u∂y³,(x,y,z)->w₂₂₂(x,y))
# h3,h2,h1,l2 = ops[5](elements["Ω"])
