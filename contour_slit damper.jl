using ApproxOperator, JLD, GLMakie 

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

ndiv = 19
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
Cᵢᵢᵢᵢ = E/(1-ν^2)
Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
Cᵢⱼᵢⱼ = E/2/(1+ν)
ds = Dict(load("jld/slit damper_hr.jld"))
# ds = Dict(load("jld/slit damper_nitsche.jld"))
# ds = Dict(load("jld/slit damper_penalty.jld"))

push!(nodes,:d₁=>ds["d₁"],:d₂=>ds["d₂"])

ind = 21
xs1 = zeros(3*ind,ind)
ys1 = zeros(3*ind,ind)
xs2 = zeros(ind,ind)
ys2 = zeros(ind,ind)
xs3 = zeros(ind,ind)
ys3 = zeros(ind,ind)
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
xs9 = zeros(ind,ind)
ys9 = zeros(ind,ind)
xs10 = zeros(ind,ind)
ys10 = zeros(ind,ind)
xs11 = zeros(3*ind,ind)
ys11 = zeros(3*ind,ind)

xl₁ = zeros(3*ind)
yl₁ = zeros(3*ind)
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
xl₂₄ = zeros(3*ind)
yl₂₄ = zeros(3*ind)
σ₁₂ = zeros(3*ind,ind)
σ₁₁ = zeros(3*ind,ind)
σ₂₂ = zeros(3*ind,ind)
σ2₁₂ = zeros(ind,ind)
σ2₁₁ = zeros(ind,ind)
σ2₂₂ = zeros(ind,ind)
σ3₁₂ = zeros(ind,ind)
σ3₁₁ = zeros(ind,ind)
σ3₂₂ = zeros(ind,ind)
σ4₁₂ = zeros(ind,ind)
σ4₁₁ = zeros(ind,ind)
σ4₂₂ = zeros(ind,ind)
σ5₁₂ = zeros(ind,ind)
σ5₁₁ = zeros(ind,ind)
σ5₂₂ = zeros(ind,ind)
σ6₁₂ = zeros(ind,ind)
σ6₁₁ = zeros(ind,ind)
σ6₂₂ = zeros(ind,ind)
σ7₁₂ = zeros(ind,ind)
σ7₁₁ = zeros(ind,ind)
σ7₂₂ = zeros(ind,ind)
σ8₁₂ = zeros(ind,ind)
σ8₁₁ = zeros(ind,ind)
σ8₂₂ = zeros(ind,ind)
σ9₁₂ = zeros(ind,ind)
σ9₁₁ = zeros(ind,ind)
σ9₂₂ = zeros(ind,ind)
σ10₁₂ = zeros(ind,ind)
σ10₁₁ = zeros(ind,ind)
σ10₂₂ = zeros(ind,ind)
σ11₁₂ = zeros(3*ind,ind)
σ11₁₁ = zeros(3*ind,ind)
σ11₂₂ = zeros(3*ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0,240, 3*ind))
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ₁₁[I,J] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ₁₂[I,J] = Cᵢⱼᵢⱼ*ε₁₂
        σ₂₂[I,J] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs1[I,J] = ξ¹+α*u₁
        ys1[I,J] = ξ²+α*u₂
        xl₂[J] = 0
        yl₂[J] = ξ²
        xl₃[J] = 240
        yl₃[J] = ξ²
    end
    xl₁[I] = ξ¹
    yl₁[I] = 0
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ2₁₁[i,j] =Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ2₁₂[i,j] =Cᵢⱼᵢⱼ*ε₁₂
        σ2₂₂[i,j] =Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs2[i,j] = x₂+α*u₁
        ys2[i,j] = y₂+α*u₂
        xl₄[i]=40+(40-(40*Δy-Δy^2)^0.5)
        yl₄[i]=y₂
        xl₅[i]=40-(40-(40*Δy-Δy^2)^0.5)
        yl₅[i]=y₂
    end
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ3₁₁[i,j] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ3₁₂[i,j] = Cᵢⱼᵢⱼ*ε₁₂
        σ3₂₂[i,j] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs3[i,j] = x₃+α*u₁
        ys3[i,j] = y₃+α*u₂
        xl₆[i]=120+(40-(40*Δy-Δy^2)^0.5)
        yl₆[i]=y₃
        xl₇[i]=120-(40-(40*Δy-Δy^2)^0.5)
        yl₇[i]=y₃
    end
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ4₁₁[i,j] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ4₁₂[i,j] = Cᵢⱼᵢⱼ*ε₁₂
        σ4₂₂[i,j] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs4[i,j] = x₄+α*u₁
        ys4[i,j] = y₄+α*u₂
        xl₈[i]=200+(40-(40*Δy-Δy^2)^0.5)
        yl₈[i]=y₄
        xl₉[i]=200-(40-(40*Δy-Δy^2)^0.5)
        yl₉[i]=y₄
    end
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ5₁₁[I,J] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ5₁₂[I,J] = Cᵢⱼᵢⱼ*ε₁₂
        σ5₂₂[I,J] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs5[I,J] = ξ¹+α*u₁
        ys5[I,J] = ξ²+α*u₂
        xl₁₀[J]=20
        yl₁₀[J]=ξ²
        xl₁₁[J]=60
        yl₁₁[J]=ξ²
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ6₁₁[I,J] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ6₁₂[I,J] = Cᵢⱼᵢⱼ*ε₁₂
        σ6₂₂[I,J] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs6[I,J] = ξ¹+α*u₁
        ys6[I,J] = ξ²+α*u₂
        xl₁₂[J]=100
        yl₁₂[J]=ξ²
        xl₁₃[J]=140
        yl₁₃[J]=ξ²
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ7₁₁[I,J] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ7₁₂[I,J] = Cᵢⱼᵢⱼ*ε₁₂
        σ7₂₂[I,J] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs7[I,J] = ξ¹+α*u₁
        ys7[I,J] = ξ²+α*u₂
        xl₁₄[J]=180
        yl₁₄[J]=ξ²
        xl₁₅[J]=220
        yl₁₅[J]=ξ²
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ8₁₁[i,j] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ8₁₂[i,j] = Cᵢⱼᵢⱼ*ε₁₂
        σ8₂₂[i,j] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs8[i,j] = x₅+α*u₁
        ys8[i,j] = y₅+α*u₂
        xl₁₆[i]=40+(40-(400-Δy^2)^0.5)
        yl₁₆[i]=y₅
        xl₁₇[i]=40-(40-(400-Δy^2)^0.5)
        yl₁₇[i]=y₅
    end
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ9₁₁[i,j] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ9₁₂[i,j] = Cᵢⱼᵢⱼ*ε₁₂
        σ9₂₂[i,j] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs9[i,j] = x₆+α*u₁
        ys9[i,j] = y₆+α*u₂
        xl₁₈[i]=120+(40-(400-Δy^2)^0.5)
        yl₁₈[i]=y₆
        xl₁₉[i]=120-(40-(400-Δy^2)^0.5)
        yl₁₉[i]=y₆
    end
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ10₁₁[i,j] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ10₁₂[i,j] = Cᵢⱼᵢⱼ*ε₁₂
        σ10₂₂[i,j] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs10[i,j] = x₇+α*u₁
        ys10[i,j] = y₇+α*u₂
        xl₂₀[i]=200+(40-(400-Δy^2)^0.5)
        yl₂₀[i]=y₇
        xl₂₁[i]=200-(40-(400-Δy^2)^0.5)
        yl₂₁[i]=y₇
    end
end

for (I,ξ¹) in enumerate(LinRange(0.0,240, 3*ind))
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
        set∇𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        ε₁₁ = 0.0
        ε₁₂ = 0.0
        ε₂₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₂ 
            ε₁₂ += (B₂[i]*xᵢ.d₂+B₁[i]*xᵢ.d₁)/2
            ε₂₂ += B₂[i]*xᵢ.d₁
        end
        σ11₁₁[I,J] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ11₁₂[I,J] = Cᵢⱼᵢⱼ*ε₁₂
        σ11₂₂[I,J] = Cᵢᵢᵢᵢ*ε₂₂+Cᵢᵢⱼⱼ*ε₁₁
        xs11[I,J] = ξ¹+α*u₁
        ys11[I,J] = ξ²+α*u₂
        xl₂₂[J]=0
        yl₂₂[J]=ξ²
        xl₂₃[J]=240
        yl₂₃[J]=ξ²
    end
    xl₂₄[I]=ξ¹
    yl₂₄[I]=180
end

fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data, azimuth = -0.5*pi, elevation = 0.5*pi)

hidespines!(ax)
hidedecorations!(ax)
# M₁₂ colorrange = (-150000,150000) M₁₁ colorrange = (-3000000,1000000) M₂₂ colorrange = (-100000,4000000)
s = surface!(ax,xs1,ys1,   color=  σ₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs2,ys2,   color= σ2₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs3,ys3,   color= σ3₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs4,ys4,   color= σ4₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs5,ys5,   color= σ5₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs6,ys6,   color= σ6₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs7,ys7,   color= σ7₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs8,ys8,   color= σ8₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs9,ys9,   color= σ9₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs10,ys10, color=σ10₂₂, colormap=:haline,colorrange = (-100000,4000000))
s = surface!(ax,xs11,ys11, color=σ11₂₂, colormap=:haline,colorrange = (-100000,4000000))
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

Colorbar(fig[1, 2],s,ticklabelsize = 25,height = 300,ticks =0:2000000:4000000)
# save("./png/slit damper_penalty_M22.png",fig)
save("./png/slit damper_hr_M22.png",fig)
# save("./png/slit damper_nitsche_M22.png",fig)
fig
