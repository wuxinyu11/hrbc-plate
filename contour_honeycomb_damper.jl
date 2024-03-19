using ApproxOperator, JLD, GLMakie 

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

ndiv = 19
α = 4e3
gmsh.initialize()
gmsh.open("msh/honeycomb_damper.msh")
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 2,γ = 5)
nₚ = length(nodes)
s = 2.5*560/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
# gmsh.finalize()

type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
𝗠 = zeros(21)
∂𝗠∂y = zeros(21)
∂𝗠∂x = zeros(21)
∂²𝗠∂x² = zeros(21)
∂²𝗠∂y² = zeros(21)
∂²𝗠∂x∂y = zeros(21)
h = 5
E = 2e11
ν = 0.3
Dᵢᵢᵢᵢ = E*h^3/(12*(1-ν^2))
Dᵢⱼᵢⱼ = E*h^3/(24*(1+ν))
# ds = Dict(load("jld/honeycomb_damper_hr.jld"))
ds = Dict(load("jld/honeycomb_damper_nitsche.jld"))
# ds = Dict(load("jld/honeycomb_damper_penalty.jld"))

push!(nodes,:d₁=>ds["d₁"],:d₂=>ds["d₂"])

ind = 31
xs1 = zeros(ind,ind)
ys1 = zeros(ind,ind)
xs2 = zeros(ind,ind)
ys2 = zeros(ind,ind)
xs3 = zeros(ind,ind)
ys3 = zeros(ind,ind)
zs3 = zeros(ind,ind)
xs4 = zeros(ind,ind)
ys4 = zeros(ind,ind)
xs5 = zeros(ind,ind)
ys5 = zeros(ind,ind)
xs6 = zeros(ind,ind)
ys6 = zeros(ind,ind)
xs7 = zeros(ind,ind)
ys7 = zeros(ind,ind)
xs8 = zeros(ind,ind)
ys8 = zeros(ind,ind)

xl₁ = zeros(ind)
yl₁ = zeros(ind)
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
xl₂₅ = zeros(ind)
yl₂₅ = zeros(ind)
xl₂₆ = zeros(ind)
yl₂₆ = zeros(ind)
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
for (I,ξ¹) in enumerate(LinRange(0.0,560, ind))
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
        xs1[I,J] = ξ¹+α*u₁
        ys1[I,J] = ξ²+α*u₂
        xl₂[J] = 0
        yl₂[J] = ξ²
        xl₃[J] = 560
        yl₃[J] = ξ²
    end
    xl₁[I] = ξ¹
    yl₁[I] = 0
    for (i,x) in enumerate(LinRange(0.0,60, ind))
        xl₄[i]=x
        yl₄[i]=40
        xl₅[i]=x+500
        yl₅[i]=40
    end
    for (i,x) in enumerate(LinRange(180,220, ind))
        xl₆[i]=x
        yl₆[i]=40
        xl₇[i]=x+160
        yl₇[i]=40
    end
end

for i in 1:ind
    for j in 1:ind
        Δy= 140/(ind-1)
        Δx=(-2*((ind-i)Δy)/7+60)/(ind-1)*2
        x₂ = 120  + (j-(ind+1)/2)*Δx
        y₂ = 180 - (i-1)*Δy
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M2₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M2₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M2₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs2[i,j] = x₂+α*u₁
        ys2[i,j] = y₂+α*u₂
        xl₈[i]=120+(-2*((ind-i)Δy)/7+60)
        yl₈[i]=y₂
        xl₉[i]=120-(-2*((ind-i)Δy)/7+60)
        yl₉[i]=y₂
    end
end

for i in 1:ind
    for j in 1:ind
        Δy= 140/(ind-1)
        Δx=(-2*((ind-i)Δy)/7+60)/(ind-1)*2
        x₃ = 280  + (j-(ind+1)/2)*Δx
        y₃ = 180 - (i-1)*Δy
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M3₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M3₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M3₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs3[i,j] = x₃+α*u₁
        ys3[i,j] = y₃+α*u₂
        xl₁₀[i]=280+(-2*((ind-i)Δy)/7+60)
        yl₁₀[i]=y₃
        xl₁₁[i]=280-(-2*((ind-i)Δy)/7+60)
        yl₁₁[i]=y₃
    end
end
for i in 1:ind
    for j in 1:ind
        Δy= 140/(ind-1)
        Δx=(-2*((ind-i)Δy)/7+60)/(ind-1)*2
        x₄ = 440 + (j-(ind+1)/2)*Δx
        y₄ = 180 - (i-1)*Δy
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M4₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M4₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M4₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs4[i,j] = x₄+α*u₁
        ys4[i,j] = y₄+α*u₂
        xl₁₂[i]=440+(-2*((ind-i)Δy)/7+60)
        yl₁₂[i]=y₄
        xl₁₃[i]=440-(-2*((ind-i)Δy)/7+60)
        yl₁₃[i]=y₄
    end
end

for i in 1:ind
    for j in 1:ind
        Δy= 140/(ind-1)
        Δx=(2*((ind-i)Δy)/7+20)/(ind-1)*2
        x₅ = 120  + (j-(ind+1)/2)*Δx
        y₅ = 320 - (i-1)*Δy
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M5₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M5₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M5₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs5[i,j] = x₅+α*u₁
        ys5[i,j] = y₅+α*u₂
        xl₁₄[i]=120+(2*((ind-i)Δy)/7+20)
        yl₁₄[i]=y₅
        xl₁₅[i]=120-(2*((ind-i)Δy)/7+20)
        yl₁₅[i]=y₅
    end
end
for i in 1:ind
    for j in 1:ind
        Δy= 140/(ind-1)
        Δx=(2*((ind-i)Δy)/7+20)/(ind-1)*2
        x₆ = 280  + (j-(ind+1)/2)*Δx
        y₆ = 320 - (i-1)*Δy
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M6₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M6₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M6₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs6[i,j] = x₆+α*u₁
        ys6[i,j] = y₆+α*u₂
        xl₁₆[i]=280+(2*((ind-i)Δy)/7+20)
        yl₁₆[i]=y₆
        xl₁₇[i]=280-(2*((ind-i)Δy)/7+20)
        yl₁₇[i]=y₆
    end
end

for i in 1:ind
    for j in 1:ind
        Δy= 140/(ind-1)
        Δx=(2*((ind-i)Δy)/7+20)/(ind-1)*2
        x₇ = 440 + (j-(ind+1)/2)*Δx
        y₇ = 320 - (i-1)*Δy
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M7₁₁[i,j] = Dᵢᵢᵢᵢ*κ₁₁
        M7₁₂[i,j] = Dᵢⱼᵢⱼ*κ₁₂
        M7₂₂[i,j] = Dᵢᵢᵢᵢ*κ₂₂
        xs7[i,j] = x₇+α*u₁
        ys7[i,j] = y₇+α*u₂
        xl₁₈[i]=440+(2*((ind-i)Δy)/7+20)
        yl₁₈[i]=y₇
        xl₁₉[i]=440-(2*((ind-i)Δy)/7+20)
        yl₁₉[i]=y₇
    end
end

for (I,ξ¹) in enumerate(LinRange(0.0,560, ind))
    for (J,ξ²) in enumerate(LinRange(320, 360, ind))
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
        u₁ = 0.0
        u₂ = 0.0
        κ₁₁ = 0.0
        κ₁₂ = 0.0
        κ₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            κ₁₁ += -B₁₁[i]*xᵢ.d₁
            κ₁₂ += -(B₁₂[i]*xᵢ.d₁+B₁₂[i]*xᵢ.d₂)
            κ₂₂ += -B₂₂[i]*xᵢ.d₂
        end
        M8₁₁[I,J] = Dᵢᵢᵢᵢ*κ₁₁
        M8₁₂[I,J] = Dᵢⱼᵢⱼ*κ₁₂
        M8₂₂[I,J] = Dᵢᵢᵢᵢ*κ₂₂
        xs8[I,J] = ξ¹+α*u₁
        ys8[I,J] = ξ²+α*u₂
        xl₂₀[J] = 0
        yl₂₀[J] = ξ²
        xl₂₁[J] = 560
        yl₂₁[J] = ξ²
    end
    xl₂₂[I] = ξ¹
    yl₂₂[I] = 360
    for (i,x) in enumerate(LinRange(0.0,60, ind))
        xl₂₃[i]=x
        yl₂₃[i]=320
        xl₂₄[i]=x+500
        yl₂₄[i]=320
    end
    for (i,x) in enumerate(LinRange(180,220, ind))
        xl₂₅[i]=x
        yl₂₅[i]=320
        xl₂₆[i]=x+160
        yl₂₆[i]=320
    end
end

fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data, azimuth =-0.5*pi, elevation = 0.5*pi)

hidespines!(ax)
hidedecorations!(ax)
# M₁₂ colorrange = (-200000,200000) M₁₁ colorrange = (-200000,200000) M₂₂ colorrange = (-200000,200000)
s = surface!(ax,xs1,ys1,  color=M₂₂, colormap=:haline,colorrange = (-200000,200000))
s = surface!(ax,xs2,ys2, color=M2₂₂, colormap=:haline,colorrange = (-200000,200000))
s = surface!(ax,xs3,ys3, color=M3₂₂, colormap=:haline,colorrange = (-200000,200000))
s = surface!(ax,xs4,ys4, color=M4₂₂, colormap=:haline,colorrange = (-200000,200000))
s = surface!(ax,xs5,ys5, color=M5₂₂, colormap=:haline,colorrange = (-200000,200000))
s = surface!(ax,xs6,ys6, color=M6₂₂, colormap=:haline,colorrange = (-200000,200000))
s = surface!(ax,xs7,ys7, color=M7₂₂, colormap=:haline,colorrange = (-200000,200000))
s = surface!(ax,xs8,ys8, color=M8₂₂, colormap=:haline,colorrange = (-200000,200000))

lines!(ax,xl₁,yl₁,color=:black,linestyle = :dash)
lines!(ax,xl₂,yl₂,color=:black,linestyle = :dash)
lines!(ax,xl₃,yl₃,color=:black,linestyle = :dash)
lines!(ax,xl₄,yl₄,color=:black,linestyle = :dash)
lines!(ax,xl₅,yl₅,color=:black,linestyle = :dash)
lines!(ax,xl₆,yl₆,color=:black,linestyle = :dash)
lines!(ax,xl₇,yl₇,color=:black,linestyle = :dash)
lines!(ax,xl₈,yl₈,color=:black,linestyle = :dash)
lines!(ax,xl₉,yl₉,color=:black,linestyle = :dash)
lines!(ax,xl₁₀,yl₁₀,color=:black,linestyle = :dash)
lines!(ax,xl₁₁,yl₁₁,color=:black,linestyle = :dash)
lines!(ax,xl₁₂,yl₁₂,color=:black,linestyle = :dash)
lines!(ax,xl₁₃,yl₁₃,color=:black,linestyle = :dash)
lines!(ax,xl₁₄,yl₁₄,color=:black,linestyle = :dash)
lines!(ax,xl₁₅,yl₁₅,color=:black,linestyle = :dash)
lines!(ax,xl₁₆,yl₁₆,color=:black,linestyle = :dash)
lines!(ax,xl₁₇,yl₁₇,color=:black,linestyle = :dash)
lines!(ax,xl₁₈,yl₁₈,color=:black,linestyle = :dash)
lines!(ax,xl₁₉,yl₁₉,color=:black,linestyle = :dash)
lines!(ax,xl₂₀,yl₂₀,color=:black,linestyle = :dash)
lines!(ax,xl₂₁,yl₂₁,color=:black,linestyle = :dash)
lines!(ax,xl₂₂,yl₂₂,color=:black,linestyle = :dash)
lines!(ax,xl₂₃,yl₂₃,color=:black,linestyle = :dash)
lines!(ax,xl₂₄,yl₂₄,color=:black,linestyle = :dash)
lines!(ax,xl₂₅,yl₂₅,color=:black,linestyle = :dash)
lines!(ax,xl₂₆,yl₂₆,color=:black,linestyle = :dash)

# Colorbar(fig[1, 2],s)
# save("./png/honeycomb_damper_penalty_M22.png",fig)
# save("./png/honeycomb_damper_hr_M22.png",fig)
# save("./png/honeycomb_damper_nitsche_M22.png",fig)
fig
