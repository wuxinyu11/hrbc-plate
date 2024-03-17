using ApproxOperator, JLD, GLMakie 

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

ndiv = 8
α = 2e3
gmsh.initialize()
gmsh.open("msh/slit damper.msh")
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
nₚ = length(nodes)
s = 3.5*240/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
# gmsh.finalize()

type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
𝗠 = zeros(55)
∂𝗠∂y = zeros(55)
∂𝗠∂x = zeros(55)
∂²𝗠∂x² = zeros(55)
∂²𝗠∂y² = zeros(55)
∂²𝗠∂x∂y = zeros(55)
h = 5
E = 2e11
ν = 0.3
Dᵢᵢᵢᵢ = E*h^3/(12*(1-ν^2))
Dᵢⱼᵢⱼ = E*h^3/(24*(1+ν))
# ds = Dict(load("jld/slit_hr.jld"))
# ds = Dict(load("jld/slit_nitsche.jld"))
ds = Dict(load("jld/slit_penalty.jld"))

push!(nodes,:d=>ds["d"])

ind = 41
xs1 = zeros(ind,ind)
ys1 = zeros(ind,ind)
zs1 = zeros(ind,ind)
cs1 = zeros(ind,ind)
ys1s = zeros(ind,ind)
zs1s = zeros(ind,ind)
ys1x = zeros(ind,ind)
zs1x = zeros(ind,ind)
xs2 = zeros(ind,ind)
ys2 = zeros(ind,ind)
zs2 = zeros(ind,ind)
cs2 = zeros(ind,ind)
ys2x = zeros(ind,ind)
zs2x = zeros(ind,ind)
ys2s = zeros(ind,ind)
zs2s = zeros(ind,ind)
xs3 = zeros(ind,ind)
ys3 = zeros(ind,ind)
zs3 = zeros(ind,ind)
cs3 = zeros(ind,ind)
xs4 = zeros(ind,ind)
ys4 = zeros(ind,ind)
zs4 = zeros(ind,ind)
xs5 = zeros(ind,ind)
ys5 = zeros(ind,ind)
zs5 = zeros(ind,ind)
xs6 = zeros(ind,ind)
ys6 = zeros(ind,ind)
zs6 = zeros(ind,ind)
xs7 = zeros(ind,ind)
ys7 = zeros(ind,ind)
zs7 = zeros(ind,ind)
xs8 = zeros(ind,ind)
ys8 = zeros(ind,ind)
zs8 = zeros(ind,ind)
xs9 = zeros(ind,ind)
ys9 = zeros(ind,ind)
zs9 = zeros(ind,ind)
xs10 = zeros(ind,ind)
ys10 = zeros(ind,ind)
zs10 = zeros(ind,ind)
xs11 = zeros(ind,ind)
ys11 = zeros(ind,ind)
zs11 = zeros(ind,ind)
ys3s = zeros(ind,ind)
zs3s = zeros(ind,ind)
ys3x = zeros(ind,ind)
zs3x = zeros(ind,ind)
ys4s = zeros(ind,ind)
zs4s = zeros(ind,ind)
ys4x = zeros(ind,ind)
zs4x = zeros(ind,ind)
ys5s = zeros(ind,ind)
zs5s = zeros(ind,ind)
ys5x = zeros(ind,ind)
zs5x = zeros(ind,ind)
ys6s = zeros(ind,ind)
zs6s = zeros(ind,ind)
ys6x = zeros(ind,ind)
zs6x = zeros(ind,ind)
ys7s = zeros(ind,ind)
zs7s = zeros(ind,ind)
ys7x = zeros(ind,ind)
zs7x = zeros(ind,ind)
ys8s = zeros(ind,ind)
zs8s = zeros(ind,ind)
ys8x = zeros(ind,ind)
zs8x = zeros(ind,ind)
ys9s = zeros(ind,ind)
zs9s = zeros(ind,ind)
ys9x = zeros(ind,ind)
zs9x = zeros(ind,ind)
ys10s = zeros(ind,ind)
zs10s = zeros(ind,ind)
ys10x = zeros(ind,ind)
zs10x = zeros(ind,ind)
ys11s = zeros(ind,ind)
zs11s = zeros(ind,ind)
ys11x = zeros(ind,ind)
zs11x = zeros(ind,ind)
z = zeros(ind,ind)
y = zeros(ind,ind)
xl₁ = zeros(ind)
yl₁ = zeros(ind)
zl = zeros(ind)
xl₂ = zeros(ind)
yl₂ = zeros(ind)
xl₃ = zeros(ind)
yl₃ = zeros(ind)
xl₄ = zeros(ind)
yl₄ = zeros(ind)
xl₅ = zeros(ind)
yl₅ = zeros(ind)
xl₆ = zeros(ind)
yl₆ = zeros(ind)
xl₇ = zeros(ind)
yl₇ = zeros(ind)
xl₈ = zeros(ind)
yl₈ = zeros(ind)
xl₉ = zeros(ind)
yl₉ = zeros(ind)
xl₁₀ = zeros(ind)
yl₁₀ = zeros(ind)
xl₁₁ = zeros(ind)
yl₁₁ = zeros(ind)
xl₁₂ = zeros(ind)
yl₁₂ = zeros(ind)
xl₁₃ = zeros(ind)
yl₁₃ = zeros(ind)
xl₁₄ = zeros(ind)
yl₁₄ = zeros(ind)
xl₁₅ = zeros(ind)
yl₁₅ = zeros(ind)
xl₁₆ = zeros(ind)
yl₁₆ = zeros(ind)
xl₁₇ = zeros(ind)
yl₁₇ = zeros(ind)
xl₁₈ = zeros(ind)
yl₁₈ = zeros(ind)
xl₁₉ = zeros(ind)
yl₁₉ = zeros(ind)
xl₂₀ = zeros(ind)
yl₂₀ = zeros(ind)
xl₂₁ = zeros(ind)
yl₂₁ = zeros(ind)
xl₂₂ = zeros(ind)
yl₂₂ = zeros(ind)
xl₂₃ = zeros(ind)
yl₂₃ = zeros(ind)
xl₂₄ = zeros(ind)
yl₂₄ = zeros(ind)


yl1ₛ = zeros(ind)
yl1ₓ = zeros(ind)
zl1ₛ = zeros(ind)
zl1ₓ = zeros(ind)
yl2ₛ = zeros(ind)
yl2ₓ = zeros(ind)
zl2ₛ = zeros(ind)
zl2ₓ = zeros(ind)
yl3ₛ = zeros(ind)
yl3ₓ = zeros(ind)
zl3ₛ = zeros(ind)
zl3ₓ = zeros(ind)
yl4ₛ = zeros(ind)
yl4ₓ = zeros(ind)
zl4ₛ = zeros(ind)
zl4ₓ = zeros(ind)
yl5ₛ = zeros(ind)
yl5ₓ = zeros(ind)
zl5ₛ = zeros(ind)
zl5ₓ = zeros(ind)
yl6ₛ = zeros(ind)
yl6ₓ = zeros(ind)
zl6ₛ = zeros(ind)
zl6ₓ = zeros(ind)
yl7ₛ = zeros(ind)
yl7ₓ = zeros(ind)
zl7ₛ = zeros(ind)
zl7ₓ = zeros(ind)
yl8ₛ = zeros(ind)
yl8ₓ = zeros(ind)
zl8ₛ = zeros(ind)
zl8ₓ = zeros(ind)
yl9ₛ = zeros(ind)
yl9ₓ = zeros(ind)
zl9ₛ = zeros(ind)
zl9ₓ = zeros(ind)
yl10ₛ = zeros(ind)
yl10ₓ = zeros(ind)
zl10ₛ = zeros(ind)
zl10ₓ = zeros(ind)
yl11ₛ = zeros(ind)
yl11ₓ = zeros(ind)
zl11ₛ = zeros(ind)
zl11ₓ = zeros(ind)
yl12ₛ = zeros(ind)
yl12ₓ = zeros(ind)
zl12ₛ = zeros(ind)
zl12ₓ = zeros(ind)
yl13ₛ = zeros(ind)
yl13ₓ = zeros(ind)
zl13ₛ = zeros(ind)
zl13ₓ = zeros(ind)
yl14ₛ = zeros(ind)
yl14ₓ = zeros(ind)
zl14ₛ = zeros(ind)
zl14ₓ = zeros(ind)
yl15ₛ = zeros(ind)
yl15ₓ = zeros(ind)
zl15ₛ = zeros(ind)
zl15ₓ = zeros(ind)
yl16ₛ = zeros(ind)
yl16ₓ = zeros(ind)
zl16ₛ = zeros(ind)
zl16ₓ = zeros(ind)
yl17ₛ = zeros(ind)
yl17ₓ = zeros(ind)
zl17ₛ = zeros(ind)
zl17ₓ = zeros(ind)
yl18ₛ = zeros(ind)
yl18ₓ = zeros(ind)
zl18ₛ = zeros(ind)
zl18ₓ = zeros(ind)
yl19ₛ = zeros(ind)
yl19ₓ = zeros(ind)
zl19ₛ = zeros(ind)
zl19ₓ = zeros(ind)
yl20ₛ = zeros(ind)
yl20ₓ = zeros(ind)
zl20ₛ = zeros(ind)
zl20ₓ = zeros(ind)
yl21ₛ = zeros(ind)
yl21ₓ = zeros(ind)
zl21ₛ = zeros(ind)
zl21ₓ = zeros(ind)
yl22ₛ = zeros(ind)
yl22ₓ = zeros(ind)
zl22ₛ = zeros(ind)
zl22ₓ = zeros(ind)
yl23ₛ = zeros(ind)
yl23ₓ = zeros(ind)
zl23ₛ = zeros(ind)
zl23ₓ = zeros(ind)
yl24ₛ = zeros(ind)
yl24ₓ = zeros(ind)
zl24ₛ = zeros(ind)
zl24ₓ = zeros(ind)
xld1 = zeros(ind)
yld1 = zeros(ind)
zld1 = zeros(ind)
xld2 = zeros(ind)
yld2 = zeros(ind)
zld2 = zeros(ind)
xld3 = zeros(ind)
yld3 = zeros(ind)
zld3 = zeros(ind)
xld4 = zeros(ind)
yld4 = zeros(ind)
zld4 = zeros(ind)
xld5 = zeros(ind)
yld5 = zeros(ind)
zld5 = zeros(ind)
xld6 = zeros(ind)
yld6 = zeros(ind)
zld6 = zeros(ind)
xld7 = zeros(ind)
yld7 = zeros(ind)
zld7 = zeros(ind)
xld8 = zeros(ind)
yld8 = zeros(ind)
zld8 = zeros(ind)
xld9 = zeros(ind)
yld9 = zeros(ind)
zld9 = zeros(ind)
xld10 = zeros(ind)
yld10 = zeros(ind)
zld10 = zeros(ind)
M₁₂ = zeros(ind,ind)
M₁₁ = zeros(ind,ind)
M₂₂ = zeros(ind,ind)
M2₁₂ = zeros(ind,ind)
M2₁₁ = zeros(ind,ind)
M2₂₂ = zeros(ind,ind)
M3₁₂ = zeros(ind,ind)
M3₁₁ = zeros(ind,ind)
M3₂₂ = zeros(ind,ind)
M4₁₂ = zeros(ind,ind)
M4₁₁ = zeros(ind,ind)
M4₂₂ = zeros(ind,ind)
M5₁₂ = zeros(ind,ind)
M5₁₁ = zeros(ind,ind)
M5₂₂ = zeros(ind,ind)
M6₁₂ = zeros(ind,ind)
M6₁₁ = zeros(ind,ind)
M6₂₂ = zeros(ind,ind)
M7₁₂ = zeros(ind,ind)
M7₁₁ = zeros(ind,ind)
M7₂₂ = zeros(ind,ind)
M8₁₂ = zeros(ind,ind)
M8₁₁ = zeros(ind,ind)
M8₂₂ = zeros(ind,ind)
M9₁₂ = zeros(ind,ind)
M9₁₁ = zeros(ind,ind)
M9₂₂ = zeros(ind,ind)
M10₁₂ = zeros(ind,ind)
M10₁₁ = zeros(ind,ind)
M10₂₂ = zeros(ind,ind)
M11₁₂ = zeros(ind,ind)
M11₁₁ = zeros(ind,ind)
M11₂₂ = zeros(ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0,240, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, 20, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        # set𝝭!(ap)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
           xs1[I,J] = ξ¹
           ys1[I,J] = ξ²
           zs1[I,J] = α*u₃
           ys1s[I,J] = ξ²+α*u₃/(180-ξ²)*h
           zs1s[I,J] = α*u₃+5
           ys1x[I,J] = ξ²-α*u₃/(180-ξ²)*h
           zs1x[I,J] = α*u₃-5
           xl₂[J] = 0
           yl₂[J] = ξ²
           xl₃[J] = 240
           yl₃[J] = ξ²
           yl2ₛ[J] = ys1s[1,J]
           zl2ₛ[J] = zs1s[1,J]
           yl2ₓ[J] = ys1x[1,J]
           zl2ₓ[J] = zs1x[1,J]
           yl3ₛ[J] = ys1s[ind,J]
           zl3ₛ[J] = zs1s[ind,J]
           yl3ₓ[J] = ys1x[ind,J]
           zl3ₓ[J] = zs1x[ind,J]
    end
    xl₁[I] = ξ¹
    yl₁[I] = 0
    zl[I] = 0
    yl1ₛ[I] = ys1s[I,1]
    zl1ₛ[I] = zs1s[I,1]
    yl1ₓ[I] = ys1x[I,1]
    zl1ₓ[I] = zs1x[I,1]
end

for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*20/(ind-1)
        Δx=(j-(ind+1)/2)*(40-(40*Δy-Δy^2)^0.5)/(ind-1)*2
        x₂ = 40 + Δx
        y₂ = 20 + Δy
        indices = sp(x₂,y₂,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[x₂]),:y=>(2,[y₂]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M2₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M2₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M2₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs2[i,j] = x₂
        ys2[i,j] = y₂
        zs2[i,j] = α*u₃
        ys2s[i,j] = y₂+α*u₃/(180-y₂)*h
        zs2s[i,j] = α*u₃+5
        ys2x[i,j] = y₂-α*u₃/(180-y₂)*h
        zs2x[i,j] = α*u₃-5
        xl₄[i]=x₂
        yl₄[i]=y₂
        xl₅[i]=40-(40-(40*Δy-Δy^2)^0.5)
        yl₅[i]=y₂
    end
    yl4ₛ[i] = ys2s[i,1]
    zl4ₛ[i] = zs2s[i,1]
    yl4ₓ[i] = ys2x[i,1]
    zl4ₓ[i] = zs2x[i,1]   
    yl5ₛ[i] = ys2s[i,ind]
    zl5ₛ[i] = zs2s[i,ind]
    yl5ₓ[i] = ys2x[i,ind]
    zl5ₓ[i] = zs2x[i,ind]
end
for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*20/(ind-1)
        Δx=(j-(ind+1)/2)*(40-(40*Δy-Δy^2)^0.5)/(ind-1)*2
        x₃ = 120 + Δx
        y₃ = 20 + Δy
        indices = sp(x₃,y₃,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[x₃]),:y=>(2,[y₃]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M3₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M3₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M3₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs3[i,j] = x₃
        ys3[i,j] = y₃
        zs3[i,j] = α*u₃
        ys3s[i,j] = y₃+α*u₃/(180-y₃)*h
        zs3s[i,j] = α*u₃+5
        ys3x[i,j] = y₃-α*u₃/(180-y₃)*h
        zs3x[i,j] = α*u₃-5
        xl₆[i]=120+(40-(40*Δy-Δy^2)^0.5)
        yl₆[i]=y₃
        xl₇[i]=120-(40-(40*Δy-Δy^2)^0.5)
        yl₇[i]=y₃
    end
    yl6ₛ[i] = ys3s[i,1]
    zl6ₛ[i] = zs3s[i,1]
    yl6ₓ[i] = ys3x[i,1]
    zl6ₓ[i] = zs3x[i,1]   
    yl7ₛ[i] = ys3s[i,ind]
    zl7ₛ[i] = zs3s[i,ind]
    yl7ₓ[i] = ys3x[i,ind]
    zl7ₓ[i] = zs3x[i,ind]
end

for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*20/(ind-1)
        Δx=(j-(ind+1)/2)*(40-(40*Δy-Δy^2)^0.5)/(ind-1)*2
        x₄ = 200 + Δx
        y₄ = 20 + Δy
        indices = sp(x₄,y₄,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[x₄]),:y=>(2,[y₄]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M4₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M4₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M4₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs4[i,j] = x₄
        ys4[i,j] = y₄
        zs4[i,j] = α*u₃
        ys4s[i,j] = y₄+α*u₃/(180-y₄)*h
        zs4s[i,j] = α*u₃+5
        ys4x[i,j] = y₄-α*u₃/(180-y₄)*h
        zs4x[i,j] = α*u₃-5
        xl₈[i]=200+(40-(40*Δy-Δy^2)^0.5)
        yl₈[i]=y₄
        xl₉[i]=200-(40-(40*Δy-Δy^2)^0.5)
        yl₉[i]=y₄
    end
    yl8ₛ[i] = ys4s[i,1]
    zl8ₛ[i] = zs4s[i,1]
    yl8ₓ[i] = ys4x[i,1]
    zl8ₓ[i] = zs4x[i,1]   
    yl9ₛ[i] = ys4s[i,ind]
    zl9ₛ[i] = zs4s[i,ind]
    yl9ₓ[i] = ys4x[i,ind]
    zl9ₓ[i] = zs4x[i,ind]
end

for (I,ξ¹) in enumerate(LinRange(20,60, ind))
    for (J,ξ²) in enumerate(LinRange(40, 140, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        # set𝝭!(ap)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M5₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M5₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M5₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
        xs5[I,J] = ξ¹
        ys5[I,J] = ξ²
        zs5[I,J] = α*u₃
        ys5s[I,J] = ξ²+α*u₃/(180-ξ²)*h
        zs5s[I,J] = α*u₃+5
        ys5x[I,J] = ξ²-α*u₃/(180-ξ²)*h
        zs5x[I,J] = α*u₃-5
        xl₁₀[J]=20
        yl₁₀[J]=ξ²
        xl₁₁[J]=60
        yl₁₁[J]=ξ²
        yl10ₛ[J] = ys5s[1,J]
        zl10ₛ[J] = zs5s[1,J]
        yl10ₓ[J] = ys5x[1,J]
        zl10ₓ[J] = zs5x[1,J]
        yl11ₛ[J] = ys5s[ind,J]
        zl11ₛ[J] = zs5s[ind,J]
        yl11ₓ[J] = ys5x[ind,J]
        zl11ₓ[J] = zs5x[ind,J]
    end
end

for (I,ξ¹) in enumerate(LinRange(100,140, ind))
    for (J,ξ²) in enumerate(LinRange(40, 140, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        # set𝝭!(ap)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M6₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M6₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M6₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
        xs6[I,J] = ξ¹
        ys6[I,J] = ξ²
        zs6[I,J] = α*u₃
        ys6s[I,J] = ξ²+α*u₃/(180-ξ²)*h
        zs6s[I,J] = α*u₃+5
        ys6x[I,J] = ξ²-α*u₃/(180-ξ²)*h
        zs6x[I,J] = α*u₃-5
        xl₁₂[J]=100
        yl₁₂[J]=ξ²
        xl₁₃[J]=140
        yl₁₃[J]=ξ²
        yl12ₛ[J] = ys6s[1,J]
        zl12ₛ[J] = zs6s[1,J]
        yl12ₓ[J] = ys6x[1,J]
        zl12ₓ[J] = zs6x[1,J]
        yl13ₛ[J] = ys6s[ind,J]
        zl13ₛ[J] = zs6s[ind,J]
        yl13ₓ[J] = ys6x[ind,J]
        zl13ₓ[J] = zs6x[ind,J]
    end
end

for (I,ξ¹) in enumerate(LinRange(180,220, ind))
    for (J,ξ²) in enumerate(LinRange(40, 140, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        # set𝝭!(ap)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M7₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M7₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M7₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
        xs7[I,J] = ξ¹
        ys7[I,J] = ξ²
        zs7[I,J] = α*u₃
        ys7s[I,J] = ξ²+α*u₃/(180-ξ²)*h
        zs7s[I,J] = α*u₃+5
        ys7x[I,J] = ξ²-α*u₃/(180-ξ²)*h
        zs7x[I,J] = α*u₃-5
        xl₁₄[J]=180
        yl₁₄[J]=ξ²
        xl₁₅[J]=220
        yl₁₅[J]=ξ²
        yl14ₛ[J] = ys7s[1,J]
        zl14ₛ[J] = zs7s[1,J]
        yl14ₓ[J] = ys7x[1,J]
        zl14ₓ[J] = zs7x[1,J]
        yl15ₛ[J] = ys7s[ind,J]
        zl15ₛ[J] = zs7s[ind,J]
        yl15ₓ[J] = ys7x[ind,J]
        zl15ₓ[J] = zs7x[ind,J]
    end
end

for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*20/(ind-1)
        Δx=(j-(ind+1)/2)*(40-(400-Δy^2)^0.5)/(ind-1)*2
        x₅ = 40 + Δx
        y₅ = 140 + Δy
        indices = sp(x₅,y₅,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[x₅]),:y=>(2,[y₅]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M8₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M8₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M8₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs8[i,j] = x₅
        ys8[i,j] = y₅
        zs8[i,j] = α*u₃
        ys8s[i,j] = y₅+α*u₃/(180-y₅)*h
        zs8s[i,j] = α*u₃+5
        ys8x[i,j] = y₅-α*u₃/(180-y₅)*h
        zs8x[i,j] = α*u₃-5
        xl₁₆[i]=40+(40-(400-Δy^2)^0.5)
        yl₁₆[i]=y₅
        xl₁₇[i]=40-(40-(400-Δy^2)^0.5)
        yl₁₇[i]=y₅
    end
    yl16ₛ[i] = ys8s[i,1]
    zl16ₛ[i] = zs8s[i,1]
    yl16ₓ[i] = ys8x[i,1]
    zl16ₓ[i] = zs8x[i,1]   
    yl17ₛ[i] = ys8s[i,ind]
    zl17ₛ[i] = zs8s[i,ind]
    yl17ₓ[i] = ys8x[i,ind]
    zl17ₓ[i] = zs8x[i,ind]
end

for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*20/(ind-1)
        Δx=(j-(ind+1)/2)*(40-(400-Δy^2)^0.5)/(ind-1)*2
        x₆ = 120 + Δx
        y₆ = 140 + Δy
        indices = sp(x₆,y₆,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[x₆]),:y=>(2,[y₆]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M9₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M9₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M9₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs9[i,j] = x₆
        ys9[i,j] = y₆
        zs9[i,j] = α*u₃
        ys9s[i,j] = y₆+α*u₃/(180-y₆)*h
        zs9s[i,j] = α*u₃+5
        ys9x[i,j] = y₆-α*u₃/(180-y₆)*h
        zs9x[i,j] = α*u₃-5
        xl₁₈[i]=120+(40-(400-Δy^2)^0.5)
        yl₁₈[i]=y₆
        xl₁₉[i]=120-(40-(400-Δy^2)^0.5)
        yl₁₉[i]=y₆
    end
    yl18ₛ[i] = ys9s[i,1]
    zl18ₛ[i] = zs9s[i,1]
    yl18ₓ[i] = ys9x[i,1]
    zl18ₓ[i] = zs9x[i,1]   
    yl19ₛ[i] = ys9s[i,ind]
    zl19ₛ[i] = zs9s[i,ind]
    yl19ₓ[i] = ys9x[i,ind]
    zl19ₓ[i] = zs9x[i,ind]
end

for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*20/(ind-1)
        Δx=(j-(ind+1)/2)*(40-(400-Δy^2)^0.5)/(ind-1)*2
        x₇ = 200 + Δx
        y₇ = 140 + Δy
        indices = sp(x₇,y₇,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[x₇]),:y=>(2,[y₇]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M10₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M10₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M10₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs10[i,j] = x₇
        ys10[i,j] = y₇
        zs10[i,j] = α*u₃
        ys10s[i,j] = y₇+α*u₃/(180-y₇)*h
        zs10s[i,j] = α*u₃+5
        ys10x[i,j] = y₇-α*u₃/(180-y₇)*h
        zs10x[i,j] = α*u₃-5
        xl₂₀[i]=200+(40-(400-Δy^2)^0.5)
        yl₂₀[i]=y₇
        xl₂₁[i]=200-(40-(400-Δy^2)^0.5)
        yl₂₁[i]=y₇
    end
    yl20ₛ[i] = ys10s[i,1]
    zl20ₛ[i] = zs10s[i,1]
    yl20ₓ[i] = ys10x[i,1]
    zl20ₓ[i] = zs10x[i,1]   
    yl21ₛ[i] = ys10s[i,ind]
    zl21ₛ[i] = zs10s[i,ind]
    yl21ₓ[i] = ys10x[i,ind]
    zl21ₓ[i] = zs10x[i,ind]
end

for (I,ξ¹) in enumerate(LinRange(0.0,240, ind))
    for (J,ξ²) in enumerate(LinRange(160, 180, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,∂𝗠∂x),:∂𝗠∂y=>(0,∂𝗠∂y),:∂²𝗠∂x²=>(0,∂²𝗠∂x²),:∂²𝗠∂y²=>(0,∂²𝗠∂y²),:∂²𝗠∂x∂y=>(0,∂²𝗠∂x∂y)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇²𝝭!(ap)
        u₃ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
            κ₁₁ += -B₁₁[i]*xᵢ.d
            κ₁₂ += -B₁₂[i]*xᵢ.d
            κ₂₂ += -B₂₂[i]*xᵢ.d
        end
        M11₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M11₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M11₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
        xs11[I,J] = ξ¹
        ys11[I,J] = ξ²
        zs11[I,J] = α*u₃
        ys11s[I,J] = ξ²+α*u₃/(320-ξ²)*h
        zs11s[I,J] = α*u₃+5
        ys11x[I,J] = ξ²-α*u₃/(320-ξ²)*h
        zs11x[I,J] = α*u₃-5
        xl₂₂[J]=0
        yl₂₂[J]=ξ²
        xl₂₃[J]=240
        yl₂₃[J]=ξ²
        yl22ₛ[J] = ys11s[1,J]
        zl22ₛ[J] = zs11s[1,J]
        yl22ₓ[J] = ys11x[1,J]
        zl22ₓ[J] = zs11x[1,J]
        yl23ₛ[J] = ys11s[ind,J]
        zl23ₛ[J] = zs11s[ind,J]
        yl23ₓ[J] = ys11x[ind,J]
        zl23ₓ[J] = zs11x[ind,J]
    end
    xl₂₄[I]=ξ¹
    yl₂₄[I]=180
    yl24ₛ[I] = ys11s[I,ind]
    zl24ₛ[I] = zs11s[I,ind]
    yl24ₓ[I] = ys11x[I,ind]
    zl24ₓ[I] = zs11x[I,ind]
end

fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data, azimuth = -0.25*pi, elevation = 0.10*pi)

hidespines!(ax)
hidedecorations!(ax)
# M₁₂ colorrange = (-100000,100000) M₁₁ colorrange = (-800000,200000) M₂₂ colorrange = (-100000,3800000)
s = surface!(ax,zs1,xs1,ys1, color=M₂₂, colormap=:haline,colorrange = (-100000,3800000))
# s = surface!(ax,zs1,xs1,ys1, color=M₂₂, colormap=:haline, colorrange = (-110000,4730000))
# s = surface!(ax,xs1s,ys1s,zs1s, color=cs1, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs1x,ys1x,zs1x, color=cs1, colormap=:redsblues, colorrange = (-0.11,0))
s = surface!(ax,zs2,xs2,ys2, color=M2₂₂, colormap=:haline,colorrange = (-100000,3800000))
# # s = surface!(ax,xs2s,ys2s,zs2s, color=cs2, colormap=:redsblues, colorrange = (-0.11,0))
# # s = surface!(ax,xs2x,ys2x,zs2x, color=cs2, colormap=:redsblues, colorrange = (-0.11,0))
s = surface!(ax,zs3,xs3,ys3, color=M3₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs4,xs4,ys4, color=M4₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs5,xs5,ys5, color=M5₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs6,xs6,ys6, color=M6₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs7,xs7,ys7, color=M7₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs8,xs8,ys8, color=M8₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs9,xs9,ys9, color=M9₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs10,xs10,ys10, color=M10₂₂, colormap=:haline,colorrange = (-100000,3800000))
s = surface!(ax,zs11,xs11,ys11, color=M11₂₂, colormap=:haline,colorrange = (-100000,3800000))
# # s = surface!(ax,xs3s,ys3s,zs3s, color=cs3, colormap=:redsblues, colorrange = (-0.11,0))
# # s = surface!(ax,xs3x,ys3x,zs3x, color=cs3, colormap=:redsblues, colorrange = (-0.11,0))
lines!(ax,zl,xl₁,yl₁,color=:black,linestyle = :dash)
lines!(ax,zl,xl₂,yl₂,color=:black,linestyle = :dash)
lines!(ax,zl,xl₃,yl₃,color=:black,linestyle = :dash)
lines!(ax,zl,xl₄,yl₄,color=:black,linestyle = :dash)
lines!(ax,zl,xl₅,yl₅,color=:black,linestyle = :dash)
lines!(ax,zl,xl₆,yl₆,color=:black,linestyle = :dash)
lines!(ax,zl,xl₇,yl₇,color=:black,linestyle = :dash)
lines!(ax,zl,xl₈,yl₈,color=:black,linestyle = :dash)
lines!(ax,zl,xl₉,yl₉,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₀,yl₁₀,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₁,yl₁₁,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₂,yl₁₂,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₃,yl₁₃,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₄,yl₁₄,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₅,yl₁₅,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₆,yl₁₆,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₇,yl₁₇,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₈,yl₁₈,color=:black,linestyle = :dash)
lines!(ax,zl,xl₁₉,yl₁₉,color=:black,linestyle = :dash)
lines!(ax,zl,xl₂₀,yl₂₀,color=:black,linestyle = :dash)
lines!(ax,zl,xl₂₁,yl₂₁,color=:black,linestyle = :dash)
lines!(ax,zl,xl₂₂,yl₂₂,color=:black,linestyle = :dash)
lines!(ax,zl,xl₂₃,yl₂₃,color=:black,linestyle = :dash)
lines!(ax,zl,xl₂₄,yl₂₄,color=:black,linestyle = :dash)

lines!(ax,zl1ₛ,xl₁,yl1ₛ,color=:gray)
lines!(ax,zl1ₓ,xl₁,yl1ₓ,color=:gray)
lines!(ax,zl2ₛ,xl₂,yl2ₛ,color=:gray)
lines!(ax,zl2ₓ,xl₂,yl2ₓ,color=:gray)
lines!(ax,zl3ₛ,xl₃,yl3ₛ,color=:gray)
lines!(ax,zl3ₓ,xl₃,yl3ₓ,color=:gray)
lines!(ax,zl4ₛ,xl₄,yl4ₛ,color=:gray)
lines!(ax,zl4ₓ,xl₄,yl4ₓ,color=:gray)
lines!(ax,zl5ₛ,xl₅,yl5ₛ,color=:gray)
lines!(ax,zl5ₓ,xl₅,yl5ₓ,color=:gray)
lines!(ax,zl6ₛ,xl₆,yl6ₛ,color=:gray)
lines!(ax,zl6ₓ,xl₆,yl6ₓ,color=:gray)
lines!(ax,zl7ₛ,xl₇,yl7ₛ,color=:gray)
lines!(ax,zl7ₓ,xl₇,yl7ₓ,color=:gray)
lines!(ax,zl8ₛ,xl₈,yl8ₛ,color=:gray)
lines!(ax,zl8ₓ,xl₈,yl8ₓ,color=:gray)
lines!(ax,zl9ₛ,xl₉,yl9ₛ,color=:gray)
lines!(ax,zl9ₓ,xl₉,yl9ₓ,color=:gray)
lines!(ax,zl10ₛ,xl₁₀,yl10ₛ,color=:gray)
lines!(ax,zl10ₓ,xl₁₀,yl10ₓ,color=:gray)
lines!(ax,zl11ₛ,xl₁₁,yl11ₛ,color=:gray)
lines!(ax,zl11ₓ,xl₁₁,yl11ₓ,color=:gray)
lines!(ax,zl12ₛ,xl₁₂,yl12ₛ,color=:gray)
lines!(ax,zl12ₓ,xl₁₂,yl12ₓ,color=:gray)
lines!(ax,zl13ₛ,xl₁₃,yl13ₛ,color=:gray)
lines!(ax,zl13ₓ,xl₁₃,yl13ₓ,color=:gray)
lines!(ax,zl14ₛ,xl₁₄,yl14ₛ,color=:gray)
lines!(ax,zl14ₓ,xl₁₄,yl14ₓ,color=:gray)
lines!(ax,zl15ₛ,xl₁₅,yl15ₛ,color=:gray)
lines!(ax,zl15ₓ,xl₁₅,yl15ₓ,color=:gray)
lines!(ax,zl16ₛ,xl₁₆,yl16ₛ,color=:gray)
lines!(ax,zl16ₓ,xl₁₆,yl16ₓ,color=:gray)
lines!(ax,zl17ₛ,xl₁₇,yl17ₛ,color=:gray)
lines!(ax,zl17ₓ,xl₁₇,yl17ₓ,color=:gray)
lines!(ax,zl18ₛ,xl₁₈,yl18ₛ,color=:gray)
lines!(ax,zl18ₓ,xl₁₈,yl18ₓ,color=:gray)
lines!(ax,zl19ₛ,xl₁₉,yl19ₛ,color=:gray)
lines!(ax,zl19ₓ,xl₁₉,yl19ₓ,color=:gray)
lines!(ax,zl20ₛ,xl₂₀,yl20ₛ,color=:gray)
lines!(ax,zl20ₓ,xl₂₀,yl20ₓ,color=:gray)
lines!(ax,zl21ₛ,xl₂₁,yl21ₛ,color=:gray)
lines!(ax,zl21ₓ,xl₂₁,yl21ₓ,color=:gray)
lines!(ax,zl22ₛ,xl₂₂,yl22ₛ,color=:gray)
lines!(ax,zl22ₓ,xl₂₂,yl22ₓ,color=:gray)
lines!(ax,zl23ₛ,xl₂₃,yl23ₛ,color=:gray)
lines!(ax,zl23ₓ,xl₂₃,yl23ₓ,color=:gray)
lines!(ax,zl24ₛ,xl₂₄,yl24ₛ,color=:gray)
lines!(ax,zl24ₓ,xl₂₄,yl24ₓ,color=:gray)

# lines!(ax,zld1,xld1,yld1,color=:gray)
# lines!(ax,zld2,xld2,yld2,color=:gray)
# lines!(ax,zld3,xld3,yld3,color=:gray)
# lines!(ax,zld4,xld4,yld4,color=:gray)
# lines!(ax,zld5,xld5,yld5,color=:gray)
# lines!(ax,zld6,xld6,yld6,color=:gray)
# lines!(ax,zld7,xld7,yld7,color=:gray)
# lines!(ax,zld8,xld8,yld8,color=:gray)
# lines!(ax,zld9,xld9,yld9,color=:gray)
# lines!(ax,zld10,xld10,yld10,color=:gray)


# Colorbar(fig[1, 2],s)
save("./png/slit_penalty_M22.png",fig)
# save("./png/slit_hr_M22.png",fig)
# save("./png/slit_nitsche_M22.png",fig)
fig
