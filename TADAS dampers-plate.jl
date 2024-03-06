
using Revise, YAML, ApproxOperator

ndiv = 16

config = YAML.load_file("./yml/TADAS dampers-cubic.yml")
elements, nodes = importmsh("./msh/TADAS dampers.msh", config)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 3.5*12/ ndiv * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
set_memory_𝗠!(elements["Ω̃"],:∇̃²)
set_memory_𝗠!(elements["Γ9"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ2"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ3"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ4"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ5"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ10"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ7"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γ8"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
set_memory_𝗠!(elements["Γᵗ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)


# set_memory_𝗠!(elements["Γₚ₁"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₂"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₃"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₄"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₅"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₆"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₇"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₈"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₉"],:𝝭,:∇̃²)
# set_memory_𝗠!(elements["Γₚ₁₀"],:𝝭,:∇̃²)

elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]
elements["Ω∩Γᵗ"] = elements["Ω"]∩elements["Γᵗ"]
elements["Ω∩Γ2"] = elements["Ω"]∩elements["Γ2"]
elements["Ω∩Γ3"] = elements["Ω"]∩elements["Γ3"]
elements["Ω∩Γ4"] = elements["Ω"]∩elements["Γ4"]
elements["Ω∩Γ5"] = elements["Ω"]∩elements["Γ5"]
elements["Ω∩Γ7"] = elements["Ω"]∩elements["Γ7"]
elements["Ω∩Γ8"] = elements["Ω"]∩elements["Γ8"]
elements["Ω∩Γ9"] = elements["Ω"]∩elements["Γ9"]
elements["Ω∩Γ10"] = elements["Ω"]∩elements["Γ10"]

# elements["Ω∩Γₚ₁"] = elements["Ω"]∩elements["Γₚ₁"]
# elements["Ω∩Γₚ₂"] = elements["Ω"]∩elements["Γₚ₂"]
# elements["Ω∩Γₚ₃"] = elements["Ω"]∩elements["Γₚ₃"]
# elements["Ω∩Γₚ₄"] = elements["Ω"]∩elements["Γₚ₄"]
# elements["Ω∩Γₚ₅"] = elements["Ω"]∩elements["Γₚ₅"]
# elements["Ω∩Γₚ₆"] = elements["Ω"]∩elements["Γₚ₆"]
# elements["Ω∩Γₚ₇"] = elements["Ω"]∩elements["Γₚ₇"]
# elements["Ω∩Γₚ₈"] = elements["Ω"]∩elements["Γₚ₈"]
# elements["Ω∩Γₚ₉"] = elements["Ω"]∩elements["Γₚ₉"]
# elements["Ω∩Γₚ₁₀"] = elements["Ω"]∩elements["Γₚ₁₀"]



# elements["Γₚ"] = elements["Γₚ₁"]∪elements["Γₚ₂"]∪elements["Γₚ₃"]∪elements["Γₚ₄"]∪elements["Γₚ₅"]∪elements["Γₚ₆"]∪elements["Γₚ₇"]∪elements["Γₚ₈"]∪elements["Γₚ₉"]∪elements["Γₚ₁₀"]
elements["Γ"] = elements["Γᵍ"]∪elements["Γᵗ"]∪elements["Γ2"]∪elements["Γ3"]∪elements["Γ4"]∪elements["Γ5"]∪elements["Γ7"]∪elements["Γ8"]∪elements["Γ9"]∪elements["Γ10"]
# elements["Γ∩Γₚ"] = elements["Γ"]∩elements["Γₚ"]

set∇₂𝝭!(elements["Ω"])
set∇̃²𝝭!(elements["Ω̃"],elements["Ω"])
set∇∇̃²𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set∇∇̃²𝝭!(elements["Γᵗ"],elements["Ω∩Γᵗ"])
set∇∇̃²𝝭!(elements["Γ2"],elements["Ω∩Γ2"])
set∇∇̃²𝝭!(elements["Γ3"],elements["Ω∩Γ3"])
set∇∇̃²𝝭!(elements["Γ4"],elements["Ω∩Γ4"])
set∇∇̃²𝝭!(elements["Γ5"],elements["Ω∩Γ5"])
set∇∇̃²𝝭!(elements["Γ7"],elements["Ω∩Γ7"])
set∇∇̃²𝝭!(elements["Γ8"],elements["Ω∩Γ8"])
set∇∇̃²𝝭!(elements["Γ9"],elements["Ω∩Γ9"])
set∇∇̃²𝝭!(elements["Γ10"],elements["Ω∩Γ10"])




# set∇̃²𝝭!(elements["Γₚ₁"],elements["Ω∩Γₚ₁"])
# set∇̃²𝝭!(elements["Γₚ₂"],elements["Ω∩Γₚ₂"])
# set∇̃²𝝭!(elements["Γₚ₃"],elements["Ω∩Γₚ₃"])
# set∇̃²𝝭!(elements["Γₚ₄"],elements["Ω∩Γₚ₄"])
# set∇̃²𝝭!(elements["Γₚ₅"],elements["Ω∩Γₚ₅"])
# set∇̃²𝝭!(elements["Γₚ₆"],elements["Ω∩Γₚ₆"])
# set∇̃²𝝭!(elements["Γₚ₇"],elements["Ω∩Γₚ₇"])
# set∇̃²𝝭!(elements["Γₚ₈"],elements["Ω∩Γₚ₈"])
# set∇̃²𝝭!(elements["Γₚ₉"],elements["Ω∩Γₚ₉"])
# set∇̃²𝝭!(elements["Γₚ₁₀"],elements["Ω∩Γₚ₁₀"])

set∇₂𝝭!(elements["Γᵍ"])
set∇₂𝝭!(elements["Γᵗ"])
set∇₂𝝭!(elements["Γ2"])
set∇₂𝝭!(elements["Γ3"])
set∇₂𝝭!(elements["Γ4"])
set∇₂𝝭!(elements["Γ5"])
set∇₂𝝭!(elements["Γ7"])
set∇₂𝝭!(elements["Γ8"])
set∇₂𝝭!(elements["Γ9"])
set∇₂𝝭!(elements["Γ10"])




# set𝝭!(elements["Γₚ₁"])
# set𝝭!(elements["Γₚ₂"])
# set𝝭!(elements["Γₚ₃"])
# set𝝭!(elements["Γₚ₄"])
# set∇∇̄²𝝭!(elements["Γ₃"],Γᵍ=elements["Γ₃"],Γᶿ=elements["Γ₃"],Γᴾ=elements["Γₚ"])
# set∇∇̄²𝝭!(elements["Γ₄"],Γᵍ=elements["Γ₄"],Γᶿ=elements["Γ₄"],Γᴾ=elements["Γₚ"])
# set∇̄²𝝭!(elements["Γₚ"],Γᵍ=elements["Γ∩Γₚ"],Γᶿ=elements["Γ∩Γₚ"],Γᴾ=elements["Γₚ"])


set∇∇̄²𝝭!(elements["Γᵍ"],Γᵍ=elements["Γᵍ"],Γᶿ=elements["Γᵍ"])
set∇∇̄²𝝭!(elements["Γᵗ"],Γᵍ=elements["Γᵗ"],Γᶿ=elements["Γᵗ"])


E = 2E11;
ν = 0.3;
P = 100.0;
w(x,y) = 0.0
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->w(x,y))

prescribe!(elements["Γᵗ"],:V=>(x,y,z)->-P)




coefficient = (:D=>1.0,:ν=>0.3)

ops = [Operator(:∫κᵢⱼMᵢⱼdΩ,coefficient...),
       Operator(:∫wqdΩ,coefficient...),
       Operator(:∫ṼgdΓ,coefficient...),
       Operator(:∫wVdΓ,coefficient...),
       Operator(:∫M̃ₙₙθdΓ,coefficient...),
       Operator(:∫θₙMₙₙdΓ,coefficient...),
       Operator(:ΔM̃ₙₛg,coefficient...),
       Operator(:wΔMₙₛ,coefficient...),
       Operator(:H₃)]


k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)
# ops[2](elements["Ω"],f)

ops[3](elements["Γᵍ"],k,f)

# ops[4](elements["Γ₁"],f)
# ops[4](elements["Γ₂"],f)
# ops[4](elements["Γ₃"],f)
# ops[4](elements["Γ₄"],f)

ops[4](elements["Γᵗ"],f)

# ops[6](elements["Γ₁"],f)
# ops[6](elements["Γ₂"],f)
# ops[6](elements["Γ₃"],f)
# ops[6](elements["Γ₄"],f)

# d = k\f



