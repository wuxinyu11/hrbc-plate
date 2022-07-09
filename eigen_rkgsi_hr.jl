using YAML, ApproxOperator

config = YAML.load_file("./yml/rectangular_hrrk.yml")

ndiv = 20
elements, nodes = importmsh("./msh/rectangular_"*string(ndiv)*".msh", config)

nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

s = 3.1 / ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

sp = RegularGrid(nodes[:x], nodes[:y], nodes[:z], n = 2, γ = 5)
sp(elements["Ω"])
set_memory_𝝭!(elements["Ωˢ"])
set_memory_𝝭!(elements["Γ̃"])
set_memory_𝝭!(elements["Γ̃ₚ₁"])
set_memory_𝝭!(elements["Γ̃ₚ₂"])
set_memory_𝝭!(elements["Γ̃ₚ₃"])
set_memory_𝝭!(elements["Γ̃ₚ₄"])

elements["Ω∩Γ̃"] = elements["Ω"]∩elements["Γ̃"]
elements["Ω∩Γ̃ₚ₁"] = elements["Ω"]∩elements["Γ̃ₚ₁"]
elements["Ω∩Γ̃ₚ₂"] = elements["Ω"]∩elements["Γ̃ₚ₂"]
elements["Ω∩Γ̃ₚ₃"] = elements["Ω"]∩elements["Γ̃ₚ₃"]
elements["Ω∩Γ̃ₚ₄"] = elements["Ω"]∩elements["Γ̃ₚ₄"]
elements["Γ̃ₚ"] = elements["Γ̃ₚ₁"]∪elements["Γ̃ₚ₂"]∪elements["Γ̃ₚ₃"]∪elements["Γ̃ₚ₄"]
elements["Γ̃∩Γ̃ₚ"] = elements["Γ̃"]∩elements["Γ̃ₚ"]

set∇𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ωˢ"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γ̃"],elements["Ω∩Γ̃"])
set∇𝝭!(elements["Γ̃"])
set∇̃²𝝭!(elements["Γ̃ₚ₁"],elements["Ω∩Γ̃ₚ₁"])
set∇̃²𝝭!(elements["Γ̃ₚ₂"],elements["Ω∩Γ̃ₚ₂"])
set∇̃²𝝭!(elements["Γ̃ₚ₃"],elements["Ω∩Γ̃ₚ₃"])
set∇̃²𝝭!(elements["Γ̃ₚ₄"],elements["Ω∩Γ̃ₚ₄"])
set𝝭!(elements["Γ̃ₚ₁"])
set𝝭!(elements["Γ̃ₚ₂"])
set𝝭!(elements["Γ̃ₚ₃"])
set𝝭!(elements["Γ̃ₚ₄"])
set∇∇̄²𝝭!(elements["Γ̃"],Γᵍ=elements["Γ̃"])
# set∇∇̄²𝝭!(elements["Γ̃"],Γᶿ=elements["Γ̃"])
# set∇∇̄²𝝭!(elements["Γ̃"],Γᵍ=elements["Γ̃"],Γᴾ=elements["Γ̃ₚ"])
# set∇̄²𝝭!(elements["Γ̃ₚ"],Γᵍ=elements["Γ̃∩Γ̃ₚ"],Γᴾ=elements["Γ̃ₚ"])

D = 1.0
ν = 0.3

coefficient = (:D=>1.0,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫ṼgdΓ,coefficient...),
       Operator(:∫M̃ₙₙθdΓ,coefficient...),
       Operator(:ΔM̃ₙₛg,coefficient...)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

# ops[1](elements["Ωˢ"],k)
ops[2](elements["Γ̃"],k,f)
# ops[3](elements["Γ̃"],k,f)
# ops[4](elements["Γ̃ₚ"],k,f)

F = eigen(k)
F.values[1]
# sign(F.values[1])*log10(abs(F.values[1]))