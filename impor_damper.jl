
using Tensors, BenchmarkExample
import Gmsh: gmsh

function import_damper_fem(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationOrder = 4
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"], integrationOrder)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],  integrationOrder)
    # gmsh.finalize()
    return elements, nodes
end

prescribeForFem = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵍ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵍ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵍ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(ξ¹,ξ²,ξ³)->1.0)
    
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->-P)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)
end
