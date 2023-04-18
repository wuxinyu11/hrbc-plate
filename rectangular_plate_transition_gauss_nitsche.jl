using YAML, ApproxOperator,LinearAlgebra,CairoMakie

ndiv = 10
𝒑 = "cubic"
# 𝒑 = "quartic"
config = YAML.load_file("./yml/rectangular_gauss_nitsche_"*𝒑*"_GI13.yml")
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh",config)
nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

# naturall bc
# sp = ApproxOperator.RegularGrid(nodes,n=2,γ=5)

data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=3,γ=5)
data = Dict([:x=>(2,[5.0]),:y=>(2,[5.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Cubic2D,:□,:QuinticSpline,:Tri3}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
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



       

𝑎 = 10
ρ = 8000.0
h = 0.05
F₀ = 100.0
Θ = π
E = 2e11
ν = 0.3
D = E*h^3/12/(1-ν^2)
ω(m,n) = π^2*(D/ρ/h)^0.5*((m/𝑎)^2+(n/𝑎)^2)
W(x,y,m,n) = 2/𝑎/(ρ*h)^0.5*sin(m*π*x/𝑎)*sin(n*π*y/𝑎)
# η(t,m,n) = 2*F₀/(ω(m,n)^2-Θ^2)/𝑎/(ρ*h)^0.5*sin(m*π*5/𝑎)*sin(n*π*5/𝑎)*(ω(m,n)*sin(Θ*t)-Θ*sin(ω(m,n)*t))
η(t,m,n) = 2*F₀/(ω(m,n)^2-Θ^2)/𝑎/(ρ*h)^0.5*sin(m*π*5/𝑎)*sin(n*π*5/𝑎)*(sin(Θ*t)-Θ/ω(m,n)*sin(ω(m,n)*t))
function w(x,y,t)
    w_ = 0.0
    max_iter = 100
    for m in 1:max_iter
        for n in 1:max_iter
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
set𝒏!(elements["Γ₁"])
set𝒏!(elements["Γ₂"])
set𝒏!(elements["Γ₃"])
set𝒏!(elements["Γ₄"])

coefficient = (:D=>D,:ν=>ν,:ρ=>ρ,:h=>h)
ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫ρhvwdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫VgdΓ,coefficient...,:α=>1e8),   
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫MₙₙθdΓ,coefficient...,:α=>1e3*ndiv),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔMₙₛg,coefficient...,:α=>1e1),      
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],m)
ops[3](elements["Ω"],fα)

ops[4](elements["Γ₁"],kα,fα)
ops[4](elements["Γ₂"],kα,fα)
ops[4](elements["Γ₃"],kα,fα)
ops[4](elements["Γ₄"],kα,fα)
ops[8](elements["Γₚ₁"],kα,fα)
ops[8](elements["Γₚ₂"],kα,fα)
ops[8](elements["Γₚ₃"],kα,fα)
ops[8](elements["Γₚ₄"],kα,fα)

β = 0.25
γ = 0.5
# β = 0.0
# γ = 0.5
Δt = 0.01
total_time = 5.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
deflection = zeros(length(times))
dexact = zeros(length(times))
error = zeros(length(times))

v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)
                           
    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->F₀*sin(Θ*t))   
                       
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
error = deflection - dexact

f = Figure()
ax = Axis(f[1,1])
xlims!(ax, 1,5)
ax.xlabel = "time"
ax.ylabel = "deflection"

# scatter!(times[1:10:500],deflection[1:10:500],markersize = 15,color = "#C00E0E",
#    label = "gauss_nitsche")
# lines!(times[1:10:500],dexact[1:10:500],linewidth = 4,color = :black,
#     label = "exact")
lines!(times[1:10:500],error[1:10:500],linewidth = 4,color = "#C00E0E",
label = "rkgsi_hr")


    axislegend(position = :lb)

f
