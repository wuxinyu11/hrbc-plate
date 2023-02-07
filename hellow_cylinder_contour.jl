using ApproxOperator, GLMakie, CairoMakie, XLSX, YAML

ndiv = 16
config = YAML.load_file("./yml/hollow_cylinder_rkgsi_hr_quartic.yml")
elements, nodes = importmsh("./msh/hollow_cylinder_"*string(ndiv)*".msh", config)
nₚ = length(nodes)
s = zeros(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
for node in nodes
    x = node.x
    y = node.y
    r = (x^2+y^2)^0.5
    sᵢ = 4.1*r*π/4/ndiv
    node.s₁ = sᵢ
    node.s₂ = sᵢ
    node.s₃ = sᵢ
end

data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2],n=3,γ=5)

sheet = XLSX.readtable("./xlsx/hollow_cylinder_contour.xlsx", "Sheet1")

inte = 300
x = [2.0/inte*i for i in 0:inte]
y = [2.0/inte*i for i in 0:inte]


D = 1.0
ν = 0.3
c1 = 4/D/(1-ν)/3
c2 = -2/D/(1+ν)/3
r²(x,y) = x^2+y^2
w₁₁(x,y)=c1*(-2*x^2/r²(x,y)^2+1/r²(x,y))+c2
w₁₂(x,y)=c1*(-2*x*y/r²(x,y)^2)
w₂₂(x,y)=c1*(-2*y^2/r²(x,y)^2+1/r²(x,y))+c2
M₁₁(x,y) = - D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = - D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = - D*(1-ν)*w₁₂(x,y)
function Mₙ(x,y)
    r = (x^2+y^2)^0.5
    n₁ = x/r
    n₂ = y/r
    return M₁₁(x,y)*n₁*n₁+2*M₁₂(x,y)*n₁*n₂+M₂₂(x,y)*n₂*n₂
end

Mᵣᵣ = Array{Union{Missing,Float64}}(missing,inte+1,inte+1)
for j in 1:inte+1
    for i in 1:inte+1
        xᵢ = x[i]
        yᵢ = y[j]
        r = (xᵢ^2+yᵢ^2)^0.5
        if r>=1.0 && r<=2.0
            Mᵣᵣ[i,j] = Mₙ(xᵢ,yᵢ)
        end
    end
end

f = Figure(fontsize=24,fonts = (; regular = "Times New Roman"))
ax = Axis(f[1, 1])
ax.aspect = 1
hidespines!(ax)
hidedecorations!(ax)
surface!(x,y,Mᵣᵣ,
    # colormap = :balance,
    # colormap = :bluesreds,
    colormap = :haline,
    shading = false,
    colorrange=(1.,2.)
    )
contour!(x,y,Mᵣᵣ,color=:whitesmoke,levels=0.9:0.2:2.1)
arc!(Point2f(0), 1, 0, π/2,color=:black)
arc!(Point2f(0), 2, 0, π/2,color=:black)
lines!([1.,2.],[0.,0.],color=:black)
lines!([0.,0.],[1.,2.],color=:black)
Colorbar(f[1,2], colormap=:haline,limits = (1, 2),ticks = 1:0.2:2, size=50 )
f
save("./figure/exact.png",f)

# Mᵣᵣ = Array{Union{Missing,Float64}}(missing,inte+1,inte+1)
# n = 5
# d = sheet[1][n]
# for j in 1:inte+1
#     for i in 1:inte+1
#         xᵢ = x[i]
#         yᵢ = y[j]
#         r = (xᵢ^2+yᵢ^2)^0.5
#         if r>=1.0 && r<=2.0
#             Mᵣᵣ[i,j] = 0.0
#             n₁ = xᵢ/r
#             n₂ = yᵢ/r
#             indices = sp(xᵢ,yᵢ,0.0)
#             ξ = ApproxOperator.SNode((1,1,0),Dict([:x=>(2,[xᵢ]),:y=>(2,[yᵢ]),:z=>(2,[0.])]))
#             ap = ApproxOperator.ReproducingKernel{:Quartic2D,:□,:QuinticSpline,:Tri3}([nodes[ind] for ind in indices],[ξ],Dict{Symbol,ApproxOperator.SymMat}())
#             set_memory_𝗠!(ap,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²)
#             set_memory_𝝭!([ap],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²)
#             set∇²₂𝝭!(ap,ξ)
#             B₁₁ = ξ[:∂²𝝭∂x²]
#             B₁₂ = ξ[:∂²𝝭∂x∂y]
#             B₂₂ = ξ[:∂²𝝭∂y²]
#             for (ind,I) in enumerate(indices)
#                 M₁₁ = -D*(B₁₁[ind]+ν*B₂₂[ind])*d[I]
#                 M₂₂ = -D*(ν*B₁₁[ind]+B₂₂[ind])*d[I]
#                 M₁₂ = -D*(1-ν)*B₁₂[ind]*d[I]
#                 Mᵣᵣ[i,j] += (M₁₁*n₁*n₁ + 2*M₁₂*n₁*n₂ + M₂₂*n₂*n₂)
#             end
#         end
#     end
# end
# f = Figure()
# ax = Axis(f[1, 1])
# ax.aspect = 1
# hidespines!(ax)
# hidedecorations!(ax)
# surface!(x,y,Mᵣᵣ,colormap = :haline, shading = false, colorrange=(1.,2.))
# contour!(x,y,Mᵣᵣ,color=:whitesmoke,levels=0.9:0.2:2.1)
# arc!(Point2f(0), 1, 0, π/2,color=:black)
# arc!(Point2f(0), 2, 0, π/2,color=:black)
# lines!([1.,2.],[0.,0.],color=:black)
# lines!([0.,0.],[1.,2.],color=:black)
# save("./figure/"*string(n)*".png",f)