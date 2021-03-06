include("convex_algorithms.jl")


function m(_p, _fx, _∇f, _∇2f)
    return _fx + transpose(_∇f)*_p + transpose(_p)*_∇2f*_p/2.0
end


function Rho(_fx, _f, _x, _p, _m, _∇f, _∇2f)
    return (_fx - _f(_x + _p))/(_m(zeros(length(_x)), _fx, _∇f, _∇2f) - _m(_p, _fx, _∇f, _∇2f))
end


function Subproblem_Cauchy_Point(_Δk, _∇f, _∇2f)

    τ = 1.0
    nrm_Δf = norm(_∇f)
    val = transpose(_∇f) * _∇2f * _∇f

    if val > 0
        τ = min(nrm_Δf^3 / (_Δk * val), 1.0)
    end

    return -τ * _Δk/nrm_Δf * _∇f
end


function Conj_Grad_Steihaug(_Δk, _∇f, _∇2f, _ϵ, _itt)
    z = 0
    r = _∇f
    d = -r

    if norm(r) < _ϵ
        return 0
    end

    for j = 0:_itt
        if transpose(d)*_∇2f*d <= 0

        end
    end

end



function Trust_Region(_f, _x, _Δk, _Δm, _η1, _η2, _η3, _t1, _t2, _ϵ, _δ, itt)

    println("\n\n")

    for k = 1:itt
        fx = _f(_x)
        ∇f = Compute_Gradient(_f, _x)
        ∇2f = Compute_Hessian(_f, _x)


        p = Subproblem_Cauchy_Point(_Δk, ∇f, ∇2f) #Solve trust region subproblem

        ρ = Rho(fx, _f, _x, p, m, ∇f, ∇2f)

        if ρ < _η2
            _Δk *= _t1
        elseif ρ > _η3 && abs(norm(p) - _Δk) <= _δ
            _Δk = min(_t2*_Δk, _Δm)
        else
            #do nothing i.e. _Δk remains the same
        end

        if ρ > _η1
            _x += p
        else
            #_x stays the same i.e. model is poor and we need to solve another subproblem
        end

        println("Iteration: ", k)
        println("x = ", _x)
        println("f(x) = ", _f(_x))

        push!(xPlot, _x[1])
        push!(yPlot, _x[2])

        if norm(∇f) <= _ϵ
            break
        end
    end

    return _x
end
