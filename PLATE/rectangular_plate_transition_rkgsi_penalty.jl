using YAML, ApproxOperator,LinearAlgebra,CairoMakie

ndiv = 10
𝒑 = "cubic"
# 𝒑 = "quartic"
config = YAML.load_file("./yml/rectangular_rkgsi_penalty_"*𝒑*".yml")
elements,nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)
nₚ = getnₚ(elements["Ω"])

# naturall bc
data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=2,γ=5)
data = Dict([:x=>(2,[5.0]),:y=>(2,[5.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Cubic2D,:□,:QuinticSpline,:Tri3}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)


s = 3.5*10/ndiv*ones(nₚ)

push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²,:𝝭)
set_memory_𝝭!(elements["Ω̃"],:∂²𝝭∂x²)

set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])


set∇₂𝝭!(elements["Γ₁"])
set∇₂𝝭!(elements["Γ₂"])
set∇₂𝝭!(elements["Γ₃"])
set∇₂𝝭!(elements["Γ₄"])
set𝝭!(elements["Γₚ₁"])
set𝝭!(elements["Γₚ₂"])
set𝝭!(elements["Γₚ₃"])
set𝝭!(elements["Γₚ₄"])
set𝝭!(elements["Γᵗ"])

a =10
ρ = 8000
h = 0.05
F₀ = 100
Θ = π
E = 2*10e11
ν = 0.3
D = 1.0
# ω(m,n) = π^2*(D/ρ/h)^0.5*((m/a)^2+(n/b)^2)
function w(x,y,t)
    w_ = 0.0
    max_iter = 5
    for m in 1:max_iter
        for n in 1:max_iter
            ω(m,n) = π^2*(D/ρ/h)^0.5*((m/a)^2+(n/a)^2)
            W(x,y,m,n) = 2/a/(ρ*h)^0.5*sin(m*π*x/a)*sin(n*π*y/a)
            # η(t,m,n) = 2*F₀/(ω(m,n)^2-Θ^2)/a/(ρ*h)^0.5*sin(m*π/2)*sin(n*π/2)*(sin(Θ*t)-Θ/ω(m,n)*sin(ω(m,n)*t))
            η(t,m,n) = 2*F₀/(ω(m,n)^2-Θ^2)/a/(ρ*h)^0.5*sin(m*π/2)*sin(n*π/2)*(ω(m,n)*sin(Θ*t)-Θ*sin(ω(m,n)*t))
            w_ += W(x,y,m,n)*η(t,m,n)
        end
    end
    return w_
end

prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γₚ₄"],:g=>(x,y,z)->0.0)


coefficient = (:D=>1.0,:ν=>0.3,:ρ=>8000.0,:h=>0.05)

# cubic
# ndiv = 10, α = 1e7
# ndiv = 20, α = 1e8
# ndiv = 40, α = 1e9
# ndiv = 10, α = 1e10
# quartic
# ndiv = 10, α = 1e7
# ndiv = 20, α = 1e9
# ndiv = 40, α = 1e11
# ndiv = 10, α = 1e14

# cubic-wxy
# ndiv = 10, α = 1e6
# ndiv = 20, α = 1e8
# ndiv = 40, α = 1e9
# ndiv = 80, α = 1e11

# quartic-wxy
# ndiv = 10, α = 1e8
# ndiv = 20, α = 1e10
# ndiv = 40, α = 1e12
# ndiv = 80, α = 1e11

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫ρhvwdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫vgdΓ,coefficient...,:α=>1e7),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫∇𝑛vθdΓ,coefficient...,:α=>1e7),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)


ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],m)

# ops[3](elements["Ω"],fα)


# ops[4](elements["Γ₁"],kα,fα)
# ops[4](elements["Γ₂"],kα,fα)
# ops[4](elements["Γ₃"],kα,fα)
# ops[4](elements["Γ₄"],kα,fα)

# ops[5](elements["Γ₁"],k,f)
# ops[5](elements["Γ₂"],k,f)
# ops[5](elements["Γ₃"],k,f)
# ops[5](elements["Γ₄"],k,f)
# ops[6](elements["Γ₁"],f)
# ops[6](elements["Γ₂"],f)
# ops[6](elements["Γ₃"],f)
# ops[6](elements["Γ₄"],f)

# ops[4](elements["Γₚ₁"],kα,fα)
# ops[4](elements["Γₚ₂"],kα,fα)
# ops[4](elements["Γₚ₃"],kα,fα)
# ops[4](elements["Γₚ₄"],kα,fα)
# ops[7](elements["Γₚ₁"],f)
# ops[7](elements["Γₚ₂"],f)
# ops[7](elements["Γₚ₃"],f)
# ops[7](elements["Γₚ₄"],f)


Θ = π
# β = 0.25
# γ = 0.5
β = 0.25
γ = 0.5
Δt = 0.01
total_time = 5.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
deflection = zeros(length(times))
dexact = zeros(length(times))

v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)
                           
    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->100.0*sin(Θ*t))   
                       
    fₙ = zeros(nₚ)
    ops[5](elements["Γᵗ"],fₙ)

    # predictor phase
    d .+= Δt*v + Δt^2/2.0*(1.0-2.0*β)*aₙ
    v .+= Δt*(1.0-γ)*aₙ
    a = (m + β*Δt^2*(k+kα))\(fₙ+fα-(k+kα)*d)
    # Corrector phase
    d .+= β*Δt^2*a
    v .+= γ*Δt*a
    aₙ .= a

    # cal deflection
    ξ = elements["Γᵗ"][1].𝓖[1]
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(elements["Γᵗ"][1].𝓒)
        I = xᵢ.𝐼
        deflection[n] += N[i]*d[I]
    end
      # cal exact solution
      dexact[n] = w(5.0,5.0,t)
end
f = Figure()
ax = Axis(f[1,1])

scatterlines!(times,deflection)
lines!(times,dexact)

f






