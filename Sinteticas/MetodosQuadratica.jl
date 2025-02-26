# Implementação: 04 de novembro de 2024

# Principais pacotes
using LinearAlgebra, SparseArrays

# Switch to Intel MKL BLAS
using MKL

# Float128 type (interface to libquadmath system package)
# using Quadmath

# VERSÃO 07/12/23
# algumas modifcações de performance inspiradas em https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl/blob/c34fc10afa854ba605aec4722414f9e4f21ac0a9/src/cg.jl

# Conjugate Gradient method (CG) (see Algorithm 6.17 of Saad's book)
function cg(
    A,
    b;
    x0 = nothing,
    eps::Real = 1e-8,
    kmax::Int = 10000,
    history::Bool = false,
    T::DataType=Float64
    )

    # in this function is supposed that A is SPD and A, b have compatible dimensions
    # A can be sparse or dense
    A = Symmetric(A)

    n = length(b)

    # pre-allocate vectors
    x  = zeros(T,n)
    r  = similar(x)
    d  = similar(x)
    Ad = similar(x)

    # initialize objects
    if !isnothing(x0)
        copyto!(x,x0)
    end

    # r = b - A*x
    copyto!(r,b)
    mul!(r,A,x,-1.0,1.0)

    copyto!(d,r)
    rsupn = norm(r, Inf)

    # r^t * r
    rtr = dot(r,r)

    # next r^t * r
    new_rtr = rtr

    k  = 0
    st = 1

    if history
        histx = [x]
        histg = [rsupn]
    else
        histx = nothing
        histg = nothing
    end

    while true

        # test whether gradient norm is small enough
        if rsupn <= eps
            st = 0
            break
        end

        # test whether maximum number of iterations is reached
        if k >= kmax
            st = 1
            break
        end

        # Ad = A*d
        mul!(Ad,A,d)

        t = rtr / dot(d,Ad)

        x .= x .+ t .* d
        r .= r .- t .* Ad

        new_rtr = dot(r,r)

        beta = new_rtr/rtr

        d .= r .+ beta .* d

        rtr = new_rtr

        rsupn = norm(r, Inf)

        # Save the iterate for debugging
        if history
            push!(histx, x)
            push!(histg, rsupn)
        end

        k += 1
    end

    return x, rsupn, k, st, histx, histg
end

# Delayed Weighted Gradient Method (DWGM) (from Oviedo Leon)
function dwgm(
    A,
    b;
    x0 = nothing,
    eps::Real = 1e-8,
    kmax::Int = 10000,
    history::Bool = false,
    T::DataType=Float64
    )


    # in this function is supposed that A is SPD and A, b have compatible dimensions
    # A can be sparse or dense
    A = Symmetric(A)
    
    n = length(b)

    # pre-allocate vectors
    x   = zeros(T, n)
    x_1 = similar(x)

    g   = similar(x)
    g_1 = similar(x)

    y = similar(x)
    w = similar(x)

    aux = similar(x)

    # initialize objects
    if !isnothing(x0)
        copyto!(x,x0)
    end

    copyto!(x_1,x)

    copyto!(g,b)
    mul!(g, A, x, 1.0, -1.0)

    copyto!(g_1,g)

    k  = 0
    st = 1

    gsupn = norm(g, Inf)

    if history
        histx = [x]
        histg = [gsupn]
    else
        histx = nothing
        histg = nothing
    end

    while true

        # test whether gradient norm is small enough
        if gsupn <= eps
            st = 0
            break
        end

        # test whether maximum number of iterations is reached
        if k >= kmax
            st = 1
            break
        end

        mul!(w, A, g)
        t = dot(g, w)/dot(w, w)

        y   .= (g_1 .- g) .+ t .* w
        beta = dot(g_1, y)/dot(y, y)

        copyto!(aux,x)
        x .= (1.0 - beta).*x_1 .+ beta .* (x .- t .* g)
        copyto!(x_1,aux)

        copyto!(aux,g)
        g .= g_1 .- beta .* y
        copyto!(g_1,aux)

        gsupn = norm(g, Inf)

        # Save the iterate for debugging
        if history
            push!(histx, x)
            push!(histg, gsupn)
        end

        k += 1
    end

    return x, gsupn, k, st, histx, histg
end

# # Conjugate Residual Method (CR) (see Algorithm 6.19 of Saad's book)
function cr(
    A,
    b;
    x0 = nothing,
    eps::Real = 1e-8,
    kmax::Int = 10000,
    history::Bool = false,
    T::DataType=Float64
    )

    # in this function is supposed that A is SPD and A, b have compatible dimensions
    # A can be sparse or dense
    A = Symmetric(A)
    
    n = length(b)

    # pre-allocate vectors
    x  = zeros(T,n)
    r  = similar(x)
    Ar = similar(x)
    d  = similar(x)
    Ad = similar(x)

    # initialize objects
    if !isnothing(x0)
        copyto!(x,x0)
    end

    # r = b - A*x
    copyto!(r,b)
    mul!(r, A, x, -1.0, 1.0)

    copyto!(d,r)
    rsupn = norm(r, Inf)

    mul!(Ad,A,d)
    mul!(Ar,A,r)

    rAr = dot(r,Ar)

    new_rAr = rAr

    k  = 0
    st = 1

    if history
        histx = [x]
        histg = [rsupn]
    else
        histx = nothing
        histg = nothing
    end

    while true

        # test whether gradient norm is small enough
        if rsupn <= eps
            st = 0
            break
        end

        # test whether maximum number of iterations is reached
        if k >= kmax
            st = 1
            break
        end

        t = rAr/dot(Ad,Ad)

        x .= x .+ t .* d
        r .= r .- t .* Ad

        mul!(Ar,A,r)

        new_rAr = dot(r,Ar)

        beta = new_rAr/rAr

        d .= r .+ beta .* d

        Ad .= Ar .+ beta .* Ad

        rAr = new_rAr

        rsupn = norm(r, Inf)

        # Save the iterate for debugging
        if history
            push!(histx, x)
            push!(histg, rsupn)
        end

        k += 1
    end

    return x, rsupn, k, st, histx, histg
end
