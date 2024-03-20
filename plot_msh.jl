
using ApproxOperator, GLMakie

import Gmsh: gmsh

gmsh.initialize()
gmsh.open("./msh/ADAS damper.msh")
# gmsh.open("./msh/TADAS dampers.msh")
# gmsh.open("./msh/slit damper.msh")
# gmsh.open("./msh/honeycomb_damper.msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"])
elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
elements["∂Ω"] = elements["Γᵍ"]∪elements["Γᵗ"]
# gmsh.finalize()

f = Figure()

# axis
ax = Axis3(f[1, 1], perspectiveness = 0.8, aspect = :data, azimuth = -0.5*pi, elevation = 0.5*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
hidespines!(ax)
hidedecorations!(ax)

x =  nodes.x
y = nodes.y
z = 0
ps = Point3f.(x,y,z)
scatter!(ps, 
    marker=:circle,
    markersize = 5,
    color = :black
)

# elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]

    lines!(x,y,linestyle = :dash, linewidth = 0.5, color = :black)
end

# # boundaries
for elm in elements["∂Ω"]
    ξ¹ = [x.x for x in elm.𝓒]
    ξ² = [x.y for x in elm.𝓒]
    x =  [x.x for x in elm.𝓒]
    y =  [x.y for x in elm.𝓒]
    lines!(x,y,linewidth = 1.5, color = :black)
end
save("./png/ADAS damper_msh.png",f)
# save("./png/TADAS dampers_msh.png",f)
# save("./png/slit damper_msh.png",f)
# save("./png/honeycomb_damper_msh.png",f)
f