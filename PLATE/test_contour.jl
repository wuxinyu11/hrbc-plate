
using GLMakie, CairoMakie

f = Figure()
ax = Axis(f[1, 1])

x = zeros(10,10)
y = zeros(10,10)
for i in 1:10
    x[i,:] .= 1/9*(i-1)
    y[:,i] .= 1/9*(i-1)
end
contour!(x,y,rand(10,10))

f