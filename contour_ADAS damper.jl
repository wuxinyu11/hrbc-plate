using ApproxOperator, JLD, GLMakie 

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

ndiv = 8
α = 5e2
gmsh.initialize()
gmsh.open("msh/ADAS damper.msh")
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
nₚ = length(nodes)
s = 3.5*320/ndiv*ones(nₚ)
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
ds = Dict(load("jld/ADAS_hr.jld"))
# ds = Dict(load("jld/ADAS_penalty.jld"))
# ds = Dict(load("jld/ADAS_nitsche.jld"))

push!(nodes,:d=>ds["d"])

ind = 21
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
for (I,ξ¹) in enumerate(LinRange(0.0,320, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, 40, ind))
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
        M₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
           xs1[I,J] = ξ¹
           ys1[I,J] = ξ²
           zs1[I,J] = α*u₃
           ys1s[I,J] = ξ²+α*u₃/(320-ξ²)*h
           zs1s[I,J] = α*u₃+5
           ys1x[I,J] = ξ²-α*u₃/(320-ξ²)*h
           zs1x[I,J] = α*u₃-5
           xl₂[J] = 0
           yl₂[J] = ξ²
           xl₃[J] = 320
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
    for (i,x) in enumerate(LinRange(0.0,40, ind))
        xl₄[i]=x
        yl₄[i]=40
        xl₅[i]=x+280
        yl₅[i]=40
    end
    yl4ₛ[I] = ys1s[I,ind]
    zl4ₛ[I] = zs1s[I,ind]
    yl4ₓ[I] = ys1x[I,ind]
    zl4ₓ[I] = zs1x[I,ind]
    yl5ₛ[I] = ys1s[I,ind]
    zl5ₛ[I] = zs1s[I,ind]
    yl5ₓ[I] = ys1x[I,ind]
    zl5ₓ[I] = zs1x[I,ind]
end

for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*30/(ind-1)
        Δx=(j-(ind+1)/2)*(120-(60*Δy-Δy^2)^0.5)/(ind-1)*2
        x₂ = 160 + Δx
        y₂ = 40 + Δy
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
        ys2s[i,j] = y₂+α*u₃/(320-y₂)*h
        zs2s[i,j] = α*u₃+5
        ys2x[i,j] = y₂-α*u₃/(320-y₂)*h
        zs2x[i,j] = α*u₃-5
        xl₆[i]=x₂
        yl₆[i]=y₂
        xl₇[i]=160-(120-(60*Δy-Δy^2)^0.5)
        yl₇[i]=y₂
    end
    yl6ₛ[i] = ys2s[i,1]
    zl6ₛ[i] = zs2s[i,1]
    yl6ₓ[i] = ys2x[i,1]
    zl6ₓ[i] = zs2x[i,1]   
    yl7ₛ[i] = ys2s[i,ind]
    zl7ₛ[i] = zs2s[i,ind]
    yl7ₓ[i] = ys2x[i,ind]
    zl7ₓ[i] = zs2x[i,ind]
end

for i in 1:ind
    for j in 1:ind
        Δy=90/(ind-1)
        Δx=(-7*((ind-i)Δy)/9+90)/(ind-1)*2
        x₃ = 160 + (j-(ind+1)/2)*Δx
        y₃ = 160 - (i-1)*Δy
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
        ys3s[i,j] = y₃+α*u₃/(320-y₃)*h
        zs3s[i,j] = α*u₃+5
        ys3x[i,j] = y₃-α*u₃/(320-y₃)*h
        zs3x[i,j] = α*u₃-5
        xl₈[i]=160+(-7*((ind-i)Δy)/9+90)
        yl₈[i]=y₃
        xl₉[i]=160-(-7*((ind-i)Δy)/9+90)
        yl₉[i]=y₃
    end
    yl8ₛ[i] = ys3s[i,1]
    zl8ₛ[i] = zs3s[i,1]
    yl8ₓ[i] = ys3x[i,1]
    zl8ₓ[i] = zs3x[i,1]   
    yl9ₛ[i] = ys3s[i,ind]
    zl9ₛ[i] = zs3s[i,ind]
    yl9ₓ[i] = ys3x[i,ind]
    zl9ₓ[i] = zs3x[i,ind]
end

for i in 1:ind
    for j in 1:ind
        Δy=90/(ind-1)
        Δx=(7*((ind-i)Δy)/9+20)/(ind-1)*2
        x₄ = 160 + (j-(ind+1)/2)*Δx
        y₄ = 250 - (i-1)*Δy
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
        M4₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M4₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M4₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs4[i,j] = x₄
        ys4[i,j] = y₄
        zs4[i,j] = α*u₃
        ys4s[i,j] = y₄+α*u₃/(320-y₄)*h
        zs4s[i,j] = α*u₃+5
        ys4x[i,j] = y₄-α*u₃/(320-y₄)*h
        zs4x[i,j] = α*u₃-5
        xl₁₀[i]=160+(7*((ind-i)Δy)/9+20)
        yl₁₀[i]=y₄
        xl₁₁[i]=160-(7*((ind-i)Δy)/9+20)
        yl₁₁[i]=y₄
    end
    yl10ₛ[i] = ys4s[i,1]
    zl10ₛ[i] = zs4s[i,1]
    yl10ₓ[i] = ys4x[i,1]
    zl10ₓ[i] = zs4x[i,1]   
    yl11ₛ[i] = ys4s[i,ind]
    zl11ₛ[i] = zs4s[i,ind]
    yl11ₓ[i] = ys4x[i,ind]
    zl11ₓ[i] = zs4x[i,ind]
end

for i in 1:ind
    for j in 1:ind
        Δy=(ind-i)*30/(ind-1)
        Δx=(j-(ind+1)/2)*(120-(900-Δy^2)^0.5)/(ind-1)*2
        x₅ = 160 + Δx
        y₅ = 250 + Δy
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
#         # set𝝭!(ap)
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
        M5₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M5₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M5₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs5[i,j] = x₅
        ys5[i,j] = y₅
        zs5[i,j] = α*u₃
        ys5s[i,j] = y₅+α*u₃/(320-y₅)*h
        zs5s[i,j] = α*u₃+5
        ys5x[i,j] = y₅-α*u₃/(320-y₅)*h
        zs5x[i,j] = α*u₃-5
        xl₁₂[i]=x₅
        yl₁₂[i]=y₅
        xl₁₃[i]=160-(120-(900-Δy^2)^0.5)
        yl₁₃[i]=y₅
    end
    yl12ₛ[i] = ys5s[i,1]
    zl12ₛ[i] = zs5s[i,1]
    yl12ₓ[i] = ys5x[i,1]
    zl12ₓ[i] = zs5x[i,1]   
    yl13ₛ[i] = ys5s[i,ind]
    zl13ₛ[i] = zs5s[i,ind]
    yl13ₓ[i] = ys5x[i,ind]
    zl13ₓ[i] = zs5x[i,ind]
end

for (I,ξ¹) in enumerate(LinRange(0.0,320, ind))
    for (J,ξ²) in enumerate(LinRange(280, 320, ind))
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
        for (i,x) in enumerate(LinRange(0.0,40, ind))
            xl₁₄[i]=x
            yl₁₄[i]=280
            xl₁₅[i]=x+280
            yl₁₅[i]=280
        end
        ys6s[I,J] = ξ²+α*u₃/(320-ξ²)*h
        zs6s[I,J] = α*u₃+5
        ys6x[I,J] = ξ²-α*u₃/(320-ξ²)*h
        zs6x[I,J] = α*u₃-5
        ys6s[I,ind] = 320
        zs6s[I,ind] = 5
        ys6x[I,ind] = 320
        zs6x[I,ind] = -5
        yl14ₛ[J] = ys6s[I,1]
        zl14ₛ[J] = zs6s[I,1]
        yl14ₓ[J] = ys6x[I,1]
        zl14ₓ[J] = zs6x[I,1]
        yl15ₛ[J] = ys6s[I,1]
        zl15ₛ[J] = zs6s[I,1]
        yl15ₓ[J] = ys6x[I,1]
        zl15ₓ[J] = zs6x[I,1]
        xl₁₆[J] = 0
        yl₁₆[J] = ξ²
        xl₁₇[J] = 320
        yl₁₇[J] = ξ²
        yl16ₛ[J] = ys6s[1,J]
        zl16ₛ[J] = zs6s[1,J]
        yl16ₓ[J] = ys6x[1,J]
        zl16ₓ[J] = zs6x[1,J]
        yl17ₛ[J] = ys6s[1,J]
        zl17ₛ[J] = zs6s[1,J]
        yl17ₓ[J] = ys6x[1,J]
        zl17ₓ[J] = zs6x[1,J]
    end
    xl₁₈[I] = ξ¹
    yl₁₈[I] = 320
    yl18ₛ[I] = ys6s[I,ind]
    zl18ₛ[I] = zs6s[I,ind]
    yl18ₓ[I] = ys6x[I,ind]
    zl18ₓ[I] = zs6x[I,ind]
end

fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data, azimuth = -0.25*pi, elevation = 0.10*pi)

hidespines!(ax)
hidedecorations!(ax)
# M₁₂ colorrange = (-700000,700000) M₁₁ colorrange = (-1500000,1000000) M₂₂ colorrange = (-800000,8000000)
s = surface!(ax,zs1,xs1,ys1, color=M₁₂, colormap=:haline,colorrange = (-700000,700000))
# s = surface!(ax,zs1,xs1,ys1, color=M₂₂, colormap=:haline, colorrange = (-300000,300000))
# s = surface!(ax,xs1s,ys1s,zs1s, color=cs1, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs1x,ys1x,zs1x, color=cs1, colormap=:redsblues, colorrange = (-0.11,0))
s = surface!(ax,zs2,xs2,ys2, color=M2₁₂, colormap=:haline,colorrange = (-700000,700000))
# # s = surface!(ax,xs2s,ys2s,zs2s, color=cs2, colormap=:redsblues, colorrange = (-0.11,0))
# # s = surface!(ax,xs2x,ys2x,zs2x, color=cs2, colormap=:redsblues, colorrange = (-0.11,0))
s = surface!(ax,zs3,xs3,ys3, color=M3₁₂, colormap=:haline,colorrange = (-700000,700000))
s = surface!(ax,zs4,xs4,ys4, color=M4₁₂, colormap=:haline,colorrange = (-700000,700000))
s = surface!(ax,zs5,xs5,ys5, color=M5₁₂, colormap=:haline,colorrange = (-700000,700000))
s = surface!(ax,zs6,xs6,ys6, color=M6₁₂, colormap=:haline,colorrange = (-700000,700000))
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


# Colorbar(fig[1, 2],s)
# save("./png/ADAS_nitsche_M12.png",fig)
# save("./png/ADAS_hr_M12.png",fig)
# save("./png/ADAS_penalty_M12.png",fig)
fig
