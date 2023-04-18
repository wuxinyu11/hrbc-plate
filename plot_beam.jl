# rkgsi-hr
using Revise, ApproxOperator, YAML, CairoMakie
ndiv = 10
𝒑 = "cubic"
config = YAML.load_file("./yml/beam_rkgsi_hr_"*𝒑*".yml")
elements, nodes = importmsh("./msh/beam_"*string(ndiv)*".msh",config)

data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=4,γ=6)
data = Dict([:x=>(2,[5.0]),:y=>(2,[0.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Cubic1D,:□,:QuinticSpline,:Poi1}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)

set_memory_𝗠!(elements["Ω̃"],:∇̃²,:𝝭)
set_memory_𝝭!(elements["Ω̃"],:𝝭,:∂²𝝭∂x²)
set_memory_𝗠!(elements["Γ"],:𝝭,:∂𝝭∂x,:∇̃²,:∂∇̃²∂ξ)


nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

set∇₁𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set𝝭!(elements["Ω"])
set𝝭!(elements["Ω̃"])
set∇₁𝝭!(elements["Γ"])
set∇∇̃²𝝭!(elements["Γ"],elements["Ω"][[1,nₑ]])
set∇∇̄²𝝭!(elements["Γ"])
set𝝭!(elements["Γᵗ"])

F₀ = 10.0
ρ = 2500.0
h = 1.0
A = 1.0
L = 10.0
ω = π
EI = 1.0/6.0*1e6
function w(x,t)
    w_ = 0.0
    max_iter = 5
    for i in 1:2:max_iter
        ωᵢ = (i*π)^2/L^2*((EI)/(ρ*A))^0.5    
        w_ += sin((i*π)/2)*sin(i*π*x/L)/(ωᵢ^2-ω^2)*(sin(ω*t)-(ω/ωᵢ)*sin(ωᵢ*t))
    end
    w_ *= 2.0*F₀/(ρ*A*L)
    return w_    
end

ops = [
    Operator(:∫κMdx,:EI=>EI),
    Operator(:∫ρhvwdΩ,:ρ=>ρ,:h=>h),
    Operator(:∫wVdΓ),
    Operator(:Ṽg,:EI=>EI),
]
k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)

ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],m)
# ops[2](elements["Ω̃"],m)
ops[4](elements["Γ"],kα,fα)


Θ = π
β = 0.25
γ = 0.5
Δt = 0.01
total_time = 10.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
x = zeros(length(times))
deflection_hr = zeros(length(times))
dexact_hr = zeros(length(times))
error_hr= zeros(length(times))
v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)

    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->F₀*sin(Θ*t))   
                       
    fₙ = zeros(nₚ)
    ops[3](elements["Γᵗ"],fₙ)

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
        deflection_hr[n] += N[i]*d[I]
    end 

    # cal exact solution
    dexact_hr[n] = w(5.0,t)

end

error_hr = deflection_hr - dexact_hr

# gauss-GI2
using Revise, ApproxOperator, YAML, CairoMakie
ndiv = 10
𝒑 = "cubic"
config = YAML.load_file("./yml/beam_gauss_"*𝒑*"_GI2.yml")
elements, nodes = importmsh("./msh/beam_"*string(ndiv)*".msh",config)
data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=4,γ=6)
data = Dict([:x=>(2,[5.0]),:y=>(2,[0.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Cubic1D,:□,:QuinticSpline,:Poi1}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)
nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

# push!(getfield(elements["Γ"][1].𝓖[1],:data),:n₁=>(2,[-1.0,1.0]))
set∇²₁𝝭!(elements["Ω"])
set∇³₁𝝭!(elements["Γ"])
set𝝭!(elements["Γᵗ"])


F₀ = 10.0
ρ = 2500.0
h = 1.0
A = 1.0
L = 10.0
ω = π
EI = 1.0/6.0*1e6
function w(x,t)
    w_ = 0.0
    max_iter = 5
    for i in 1:2:max_iter
        ωᵢ = (i*π)^2/L^2*((EI)/(ρ*A))^0.5    
        w_ += sin((i*π)/2)*sin(i*π*x/L)/(ωᵢ^2-ω^2)*(sin(ω*t)-(ω/ωᵢ)*sin(ωᵢ*t))
    end
    w_ *= 2.0*F₀/(ρ*A*L)
    return w_    
end
ops = [
    Operator(:∫κMdx,:EI=>EI),
    Operator(:∫ρhvwdΩ,:ρ=>ρ,:h=>h),
    Operator(:∫wVdΓ),
    Operator(:Vg,:EI=>EI,:α=>1e8),
]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],m)
ops[4](elements["Γ"],kα,fα)
Θ = π
β = 0.25
γ = 0.5
Δt = 0.01
total_time = 10.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
x = zeros(length(times))
deflection_gi2 = zeros(length(times))
dexact_gi2 = zeros(length(times))
error_gi2 = zeros(length(times))
v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)

    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->F₀*sin(Θ*t))   
                       
    fₙ = zeros(nₚ)
    ops[3](elements["Γᵗ"],fₙ)

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
        deflection_gi2[n] += N[i]*d[I]
    end 

    # cal exact solution
    dexact_gi2[n] = w(5.0,t)

end
error_gi2 = deflection_gi2 - dexact_gi2
# gauss-GI5
using Revise, ApproxOperator, YAML, CairoMakie
ndiv = 10
𝒑 = "cubic"
config = YAML.load_file("./yml/beam_gauss_"*𝒑*"_GI5.yml")
elements, nodes = importmsh("./msh/beam_"*string(ndiv)*".msh",config)
data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=4,γ=6)
data = Dict([:x=>(2,[5.0]),:y=>(2,[0.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Cubic1D,:□,:QuinticSpline,:Poi1}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)
nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

# push!(getfield(elements["Γ"][1].𝓖[1],:data),:n₁=>(2,[-1.0,1.0]))
set∇²₁𝝭!(elements["Ω"])
set∇³₁𝝭!(elements["Γ"])
set𝝭!(elements["Γᵗ"])


F₀ = 10.0
ρ = 2500.0
h = 1.0
A = 1.0
L = 10.0
ω = π
EI = 1.0/6.0*1e6
function w(x,t)
    w_ = 0.0
    max_iter = 5
    for i in 1:2:max_iter
        ωᵢ = (i*π)^2/L^2*((EI)/(ρ*A))^0.5    
        w_ += sin((i*π)/2)*sin(i*π*x/L)/(ωᵢ^2-ω^2)*(sin(ω*t)-(ω/ωᵢ)*sin(ωᵢ*t))
    end
    w_ *= 2.0*F₀/(ρ*A*L)
    return w_    
end
ops = [
    Operator(:∫κMdx,:EI=>EI),
    Operator(:∫ρhvwdΩ,:ρ=>ρ,:h=>h),
    Operator(:∫wVdΓ),
    Operator(:Vg,:EI=>EI,:α=>1e8),
]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],m)
ops[4](elements["Γ"],kα,fα)
Θ = π
β = 0.25
γ = 0.5
Δt = 0.01
total_time = 10.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
x = zeros(length(times))
deflection_gi5 = zeros(length(times))
dexact_gi5 = zeros(length(times))
error_gi5 = zeros(length(times))
v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)

    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->F₀*sin(Θ*t))   
                       
    fₙ = zeros(nₚ)
    ops[3](elements["Γᵗ"],fₙ)

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
        deflection_gi5[n] += N[i]*d[I]
    end 

    # cal exact solution
    dexact_gi5[n] = w(5.0,t)

end
error_gi5 = deflection_gi5 - dexact_gi5





f = Figure()
# ytickformat = :string
ax = Axis(f[1,1])
xlims!(ax, 1,10)
# ax=Axis(f[1, 1],aspect = 1)
ax.xlabel = "time"
# ax.ylabel = "deflection"

ax.ylabel = "deflection error"
# lines!(times[1:20:1000],dexact_hr[1:20:1000],linewidth = 4,color = :black,
#     label = "exact")
# scatter!(times[1:20:1000],deflection_hr[1:20:1000],markersize = 15,color = "#C00E0E",
#    label = "rkgsi_hr")
# scatter!(times[1:20:1000],deflection_gi2[1:20:1000],marker = :utriangle,markersize = 15,color = "#7D2D89",
#     label = "gauss-GI2")
# scatter!(times[1:20:1000],deflection_gi5[1:20:1000],marker = :utriangle,markersize = 15,color = "#114A97",
#     label = "gauss-GI5")  


lines!(times[1:20:1000],error_hr[1:20:1000],linewidth = 4,color = "#C00E0E",
label = "rkgsi_hr")
lines!(times[1:20:1000],error_gi2[1:20:1000], linestyle = :dot,linewidth = 6,color = "#7D2D89",
label = "gauss-GI2")
lines!(times[1:20:1000],error_gi5[1:20:1000],linestyle = :dot,linewidth = 6,color = "#114A97",
label = "gauss-GI5")

axislegend(position = :lb)

# Legend(f[1, 2], ax)
f


