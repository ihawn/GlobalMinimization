include("graphing.jl")
include("test_objectives.jl")
include("cloud_search.jl")


x0 = [-10.0,10.0]
ϵ = 1e-8
η = 1e-3
α = 0.5
β = 0.8
κ = 1
ℓ = 1
γ = 0.5
ρ = 5
searchWidth = 10


xPlot = []
yPlot = []
solPlotX = []
solPlotY = []
searchX = []
searchY = []
finalSolX = []
finalSolY = []

var = x0
maxIterations = 150

#f(x) = x[1]^2 + x[2]^2 + 7*sin(x[1] + x[2]) + 10*sin(5x[1])
#f(x) = (x[2] - 0.129*x[1]^2 + 1.6*x[1] - 6)^2 + 6.07*cos(x[1]) + 10
#f(x) = Rastrigin(x, 2)
f(x) = Ackley(x)
#f(x) = Bukin(x)
#f(x) = Holder_Table(x)
#f(x) = Schaffer_N2(x)
#f(x) = Styblinski_Tang(x,2)

@time minimum = Grad_Cloud_Search(f, x0, α, β, η, ϵ, κ, ℓ, γ, ρ, searchWidth, maxIterations)

Contour([-5, 5], [-5, 5], 0.05, f, "Global Minimization With Gradient Cloud Search", xPlot, yPlot, solPlotX, solPlotY, finalSolX, finalSolY)
