function Contour(rangeX, rangeY, step, f, _title, xPlot, yPlot, solPlotX, solPlotY, finalSolX, finalSolY)
    plotf(x,y) = f([x, y])
    _x = -10.0:0.03:10.0
    _y = -10.0:0.03:10.0
    X = repeat(reshape(_x, 1, :), length(_y), 1)
    Y = repeat(_y, 1, length(_x))
    Z = map(plotf, X, Y)
    p1 = Plots.contour(_x,_y, plotf, fill = true)
    plot(p1, legend = false, title = _title)
    plot!(xPlot, yPlot, color = "white")
    scatter!(xPlot, yPlot, markersize = 2, color = "red")
    scatter!(solPlotX, solPlotY, color = "green")
    scatter!(finalSolX, finalSolY, color = "white")
end
