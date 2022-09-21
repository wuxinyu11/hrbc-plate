using ApproxOperator,Plots

## 1D style
# linestyle = (:solid,2)
# markstyle = (:circle,6,Plots.stroke(0, :black))
# xlims_ = (0,1)
# ylims_ = (-1,1)

## 2D style
# patch_test:
#   linewidth: 2->1.5, 12->1
#   markersize: 2->6, 6->5

xlims_ = (-0.05,1.05)
ylims_ = (-0.05,1.05)

# plate_with_hole:
#   linewidth: 3,6->2, 12->1, 24->0.5
#   markersize: 3->6, 6->5, 12->3.5, 24->1

# xlims_ = (-0.1,5.1)
# ylims_ = (-0.1,5.1)

# cantilever:
#   linewidth: 4->2, 8->1.5, 16->1, 32->0.5
#   markersize: 4->5, 8->4, 16->3, 32->1

# xlims_ = (-0.5,48.5)
# ylims_ = (-6.5,6.5)

linestyle = (:solid,1.25)
markstyle = (:circle,5,Plots.stroke(0, :black))

function plotmesh(as::Vector{T}) where T<:ApproxOperator.AbstractElement
    p = plot(;framestyle=:none,legend = false,ylims = ylims_,background_color=:transparent,aspect_ratio=:equal)
    for a in as
        plotmesh(a,p)
    end
    return p
end

function plotmesh(a::T,p::Plots.Plot) where T<:ApproxOperator.AbstractElement{:Seg2}
    plot!([a.𝓒[i].x for i in 1:2],zeros(2),color=:black,line=linestyle,marker=markstyle)
end

function plotmesh(a::T,p::Plots.Plot) where T<:ApproxOperator.AbstractElement{:Tri3}
    plot!([a.𝓒[i].x for i in (1,2,3,1)],[a.𝓒[i].y for i in (1,2,3,1)],color=:black,line=linestyle,marker=markstyle)
end

# elements, nodes = importmsh("./msh/bar.msh")
# p = plotmesh(elements["Ω"])
# p;savefig(p,"./figure/bar_mesh.svg")

elements, nodes = importmsh("./msh/patchtest.msh,config")
p = plotmesh(elements["Ω"])
savefig(p,"./figure/patchtest.svg")

# elements, nodes = importmsh("./msh/cantilever_32.msh")
# p = plotmesh(elements["Ω"])
# savefig(p,"./figure/cantilever_32.svg")

# elements, nodes = importmsh("./msh/plate_with_hole_24.msh")
# p = plotmesh(elements["Ω"])
# savefig(p,"./figure/plate_with_hole_24.svg")

plot(p)
