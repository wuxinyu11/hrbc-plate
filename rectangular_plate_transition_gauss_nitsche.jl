using YAML, ApproxOperator,LinearAlgebra

ndiv = 10
𝒑 = "cubic"
# 𝒑 = "quartic"
config = YAML.load_file("./yml/rectangular_gauss_nitsche_"*𝒑*".yml")
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh",config)

# naturall bc
# sp = ApproxOperator.RegularGrid(nodes,n=2,γ=5)

data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=2,γ=5)
data = Dict([:x=>(2,[5.0]),:y=>(2,[5.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Quadratic2D,:□,:QuinticSpline,:Tri3}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)


 
nₚ = length(nodes)

s = 3.5*10/ ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

set_memory_𝝭!(elements["Γᵗ"],:𝝭)

set∇²₂𝝭!(elements["Ω"])
set∇³𝝭!(elements["Γ₁"])
set∇³𝝭!(elements["Γ₂"])
set∇³𝝭!(elements["Γ₃"])
set∇³𝝭!(elements["Γ₄"])
set∇²₂𝝭!(elements["Γₚ₁"])
set∇²₂𝝭!(elements["Γₚ₂"])
set∇²₂𝝭!(elements["Γₚ₃"])
set∇²₂𝝭!(elements["Γₚ₄"])
set𝝭!(elements["Γᵗ"])

a=10
ρ=8000
F₀=100
Θ=π

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

set𝒏!(elements["Γ₁"])
set𝒏!(elements["Γ₂"])
set𝒏!(elements["Γ₃"])
set𝒏!(elements["Γ₄"])

coefficient = (:D=>1.0,:ν=>0.3,:ρ=>8000.0,:h=>0.05)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫ρhvwdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       # ndiv = 10, α = 1e3*ndiv^3
       Operator(:∫VgdΓ,coefficient...,:α=>1e8),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫MₙₙθdΓ,coefficient...,:α=>1e3*ndiv),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔMₙₛg,coefficient...,:α=>1e1),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
f = zeros(nₚ)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],m)
ops[3](elements["Ω"],f)
ops[4](elements["Γ₁"],k,f)
ops[4](elements["Γ₂"],k,f)
ops[4](elements["Γ₃"],k,f)
ops[4](elements["Γ₄"],k,f)

ops[8](elements["Γₚ₁"],k,f)
ops[8](elements["Γₚ₂"],k,f)
ops[8](elements["Γₚ₃"],k,f)
ops[8](elements["Γₚ₄"],k,f)




# A=eigvals(m,k)

Θ = π
β = 0.0
γ = 0.5
Δt = 0.1
total_time = 1
times = 0.0:Δt:total_time
d = zeros(nₚ)
deflection = zeros(length(times))

v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)
                           
    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->sin(Θ*t))   
                       
    f = zeros(nₚ)
    ops[5](elements["Γᵗ"],f)

    a = (m + β*Δt^2*k)\(f-k*d)
                    
    # predictor phase
    d .+= Δt*v + Δt^2/2.0*(1.0-2.0*β)*aₙ
    v .+= Δt*(1.0-γ)*aₙ

    # Corrector phase
    d .+= β*Δt^2*a
    v .+= γ*Δt*a

    # cal deflection
    ξ = elements["Γᵗ"][1].𝓖[1]
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(elements["Γᵗ"][1].𝓒)
        I = xᵢ.𝐼
        deflection[n] += N[i]*d[I]
    end
end
deflection