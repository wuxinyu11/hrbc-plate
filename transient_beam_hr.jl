using Revise, ApproxOperator, YAML, CairoMakie

ndiv = 10
# 𝒑 = "cubic"
𝒑 = "quartic"

config = YAML.load_file("./yml/beam_rkgsi_hr_"*𝒑*".yml")
elements, nodes = importmsh("./msh/beam_"*string(ndiv)*".msh",config)

data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2];n=4,γ=6)
data = Dict([:x=>(2,[5.0]),:y=>(2,[0.0]),:z=>(2,[0.0]),:𝑤=>(2,[1.0])])
ξ = ApproxOperator.SNode((1,1,0),data)
𝓒 = [nodes[i] for i in sp(ξ)]
𝗠 = Dict{Symbol,ApproxOperator.SymMat}()
elements["Γᵗ"] = [ApproxOperator.ReproducingKernel{:Quartic1D,:□,:QuinticSpline,:Poi1}(𝓒,[ξ],𝗠)]
set_memory_𝗠!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)

set_memory_𝗠!(elements["Ω̃"],:∇̃²,:𝝭)
set_memory_𝝭!(elements["Ω̃"],:𝝭,:∂²𝝭∂x²)
set_memory_𝗠!(elements["Γ"],:𝝭,:∂𝝭∂x,:∇̃²,:∂∇̃²∂ξ)


nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 4.5*10 / ndiv * ones(nₚ)
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
# β = 0.25
# γ = 0.5
β = 0.25
γ = 0.5
Δt = 0.01
total_time = 10.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
x = zeros(length(times))
deflection = zeros(length(times))
dexact = zeros(length(times))
error = zeros(length(times))
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
        deflection[n] += N[i]*d[I]
    end 

    # cal exact solution
    dexact[n] = w(5.0,t)

end

f = Figure()
ax = Axis(f[1,1])
error = deflection - dexact
scatterlines!(times,deflection)
lines!(times,dexact)
scatterlines!(times,error)
f