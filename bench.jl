using Plots
include("global_minimization.jl")


function Bench(maxDimension, f, startRange, testsPerDimension)
    cloud_search_times = zeros(maxDimension - 1)
    cloud_search_accuracy_sum = 0.0

    for i = 2:maxDimension
        println("Dimension ", i)
        x0 = (rand(i) .- 0.5)*startRange

        #cloud search
        println("Testing Cloud Search")
        sum = 0
        for k = 1:testsPerDimension
            sum += @elapsed sol = Cloud_Search(f, x0, α, β, η, ϵ, κ, ℓ, γ, ρ, searchWidth, maxIterations)[1]
            if norm(sol - zeros(i)) < 1e-5
                cloud_search_accuracy_sum += 1.0
            end
            println(sol)
        end
        cloud_search_times[i-1] = sum/testsPerDimension

    end

    println("Cloud Search Accuracy: ", 100*cloud_search_accuracy_sum/(testsPerDimension*(maxDimension - 1)), "%")
    plot(cloud_search_times)
end
