include("convex_algorithms.jl")

function Generate_ℓ_Vector(_x, _ℓ, _k, _ρ)
    return (rand(length(_x)) .- 0.5) * _ℓ
end


function Search(local_sol, _f, _ℓ, _γ, _η, _ρ, _α, _β, _κ)

    ℓ_start = _ℓ

    for k in 1:_ρ

        _ℓ = ℓ_start
        vec = Generate_ℓ_Vector(local_sol[1], _ℓ, k, _ρ)

        while Grad_Descent(_f, local_sol[1] + vec, _α, _β, _η, _κ)[2] >= local_sol[2] + _η && Grad_Descent(_f, local_sol[1] - vec, _α, _β, _η, _κ)[2] >= local_sol[2] + _η && _ℓ > _η
            _ℓ *= _γ
            vec = Generate_ℓ_Vector(local_sol[1], _ℓ, k, _ρ)


            push!(searchX, (local_sol[1] + vec)[1])
            push!(searchY, (local_sol[1] + vec)[2])
            push!(searchX, (local_sol[1] - vec)[1])
            push!(searchY, (local_sol[1] - vec)[2])
        end


        if Grad_Descent(_f, local_sol[1] + vec, _α, _β, _η, _κ)[2] < local_sol[2] - _η
            return local_sol[1] + vec, _ℓ
        elseif Grad_Descent(_f, local_sol[1] - vec, _α, _β, _η, _κ)[2] < local_sol[2] - _η
            return local_sol[1] - vec, _ℓ
        end
    end

    return local_sol[1], _ℓ
end


function Grad_Cloud_Search(_f, _x, _α, _β, _η, _ϵ, _κ, _ℓ, _γ, _ρ, width, maxIt)

    sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)
    s = Search(sol, _f, _ℓ, _γ, _η, _ρ, _α, _β, _κ)
    x_prev = s[1]
    it = 0

    while abs(norm(_x) - norm(x_prev)) > _η && it < maxIt

        s = Search(sol, _f, _ℓ, _γ, _η, _ρ, _α, _β, _κ)
        x_prev = s[1]

        sol = Grad_Descent(_f, x_prev, _α, _β, _η, _κ)
        s = Search(sol, _f, _ℓ, _γ, _η, _ρ, _α, _β, _κ)
        _x = s[1]

        it +=1
    end

    sol = Unconstrained_Newton(_f, _x, _α, _β, _κ, _ϵ, maxIt)

    println(sol)

    return sol
end
