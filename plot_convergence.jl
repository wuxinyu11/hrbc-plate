using CairoMakie, XLSX

xf = XLSX.readxlsx("./xlsx/rectangular.xlsx")
h = xf["h!C2:C5"]
L₂ = xf["L2!A2:E5"]
L₂_slope = xf["L2!F3:J5"]

L₂_figure = Figure()
Axis(L₂_figure[1, 1])
scatterlines!(Base.vect(h...),Base.vect(L₂[:,5]...), color = :red)

save("./png/test.png",L₂_figure)
