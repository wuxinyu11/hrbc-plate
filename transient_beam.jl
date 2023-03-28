using Revise, ApproxOperator, YAML,XLSX

ndiv = 10
𝒑 = "cubic"
config = YAML.load_file("./yml/beam_rkgsi_penalty_"*𝒑*".yml")
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
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝝭!(elements["Ω̃"],:∂²𝝭∂x²)

nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 3.5*10 / ndiv * ones(nₚ)
# s = 4.5 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

set∇₁𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γ"])

F₀ = 10
ρ = 2500
t = 1.0
A = 1.0
L = 10
ω = π
E = 2*10e6
I = 1/12
EI = 1.0/6.0*1e6
# function w(x,t)
#     w_ = 0.0
#     max_iter = 5
#     for i in 1:max_iter
#         ωᵢ = (i*i*π*π)/(L*L)*((E*I)/(ρ*A))^abs(1/2)    
#         # w_ += W(x,t)
#         w_ += 2*F₀/(ρ*A*L)*(sin((i*π)/2)*sin(i*π*x/L)/(ωᵢ²-ω²))*(sin(ω*t)-(ω/ωᵢ)*sin(ωᵢ*t))
#         # ωᵢ = (i*i*π*π)/L*L*((E*I)/(ρ*A))^abs(1/2)    
#     end
#     return w_    
# end
# W(x,t)= 2*F₀/(ρ*A*L)*(sin((i*π)/2)*sin(i*π*x/L)/(ωᵢ²-ω²))*(sin(ω*t)-(ω/ωᵢ)*sin(ωᵢ*t))

ops = [
    Operator(:∫κMdx,:EI=>1.0/6.0*1e6),
    Operator(:∫ρhvwdΩ,:ρ=>2500.0,:h=>1.0),
    Operator(:∫wVdΓ),
    Operator(:∫vgdΓ,:α=>1e8),
]
k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],m)
ops[4](elements["Γ"],kα,fα)

Θ = π
# β = 0.25
# γ = 0.5
β = 0.25
γ = 0.5
Δt = 0.01
total_time = 5.0
times = 0.0:Δt:total_time
d = zeros(nₚ)
x = zeros(length(times))
deflection = zeros(length(times))
v = zeros(nₚ)
aₙ = zeros(nₚ)
for (n,t) in enumerate(times)
    ωₜ = (t*t*π*π)/(L*L)*((E*I)/(ρ*A))^abs(1/2)    
    w(x,t)= 2*F₀/(ρ*A*L)*(sin((t*π)/2)*sin(t*π*x/L)/(ωₜ*ωₜ-ω*ω))*(sin(ω)-(ω/ωₜ)*sin(ωₜ))                      
    prescribe!(elements["Γᵗ"],:V=>(x,y,z)->10.0*sin(Θ*t))   
                       
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
    ξ = elements["Γᵗ"][1].𝓖[1]
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(elements["Γᵗ"][1].𝓒)
        I = xᵢ.𝐼
        x[n] += N[i]*d[I]
    end
end

x
  






# index = [10,20,40,80]
# XLSX.openxlsx("./xlsx/transient_"*𝒑*".xlsx", mode="rw") do xf
#     row = "A"
# #     row = "C"
#     D = xf[1]
#     # T = xf[3]

#     ind = findfirst(n->n==ndiv,index)+1
#     row = row*string(ind)
#     D[row] = deflection
#     # T[row] = times

# end