using YAML, ApproxOperator, XLSX
ndiv = 5
config = YAML.load_file("./yml/cantilever_rkgsi_hr.yml")
elements, nodes = importmsh("./msh/TADAS dampers.msh", config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])
 s = 3.5/ ndiv * ones(nₚ)
 push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

 set_memory_𝗠!(elements["Ω̃"],:∇̃)
 set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃)
 elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]

 set∇₂𝝭!(elements["Ω"])
set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
set∇̃𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set∇₂𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])
E = 2E11;ν = 0.3;P = 1000;

prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0)
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

coefficient = (:E=>E,:ν=>ν)
ops = [Operator{:∫∫εᵢⱼσᵢⱼdxdy}(coefficient...),
       Operator{:∫vᵢtᵢds}(coefficient...),
       Operator{:∫σᵢⱼnⱼgᵢds}(coefficient...),
       Operator{:∫σ̄ᵢⱼnⱼgᵢds}(coefficient...),
       Operator{:Hₑ_PlaneStress}(coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

 ops[1](elements["Ω̃"],k)
 ops[2](elements["Γᵗ"],f)
 ops[3](elements["Γᵍ"],k,f)
 ops[4](elements["Γᵍ"],k,f)

 d = k\f