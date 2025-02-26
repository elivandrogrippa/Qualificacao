# Principais pacotes
using LinearAlgebra, SparseArrays

# Switch to Intel MKL BLAS
using MKL

##################################### código pra plotagem de gráfico de erro #####################

function cr_dwgm(
    A,
    b;
    x0 = nothing,
    prec::Real = 1e-8,
    kmax::Int = 10000,
    T::DataType=Float64
    )
    
    # in this function is supposed that A is SPD and A, b have compatible dimensions
    # A can be sparse or dense
    A = Symmetric(A)

    n = length(b)

    # ------------------------------------------------
    # ITERAÇÃO DO MÉTODO CR
    # pre-allocate vectors
    cr_x  = zeros(T,n)
    cr_r  = similar(cr_x)
    cr_Ar = similar(cr_x)
    cr_d  = similar(cr_x)
    cr_Ad = similar(cr_x)

    # initialize objects
    if !isnothing(x0)
        copyto!(cr_x,x0)
    end

    # r = b - A*x
    copyto!(cr_r,b)
    mul!(cr_r, A, cr_x, -1.0, 1.0)

    copyto!(cr_d, cr_r)
    cr_rsupn = norm(cr_r, Inf)

    mul!(cr_Ad, A, cr_d)
    mul!(cr_Ar, A, cr_r)

    cr_rAr = dot(cr_r,cr_Ar)

    cr_new_rAr = cr_rAr

    cr_k  = 0
    cr_st = 1

    # ------------------------------------------------

    # ------------------------------------------------
    # ITERAÇÃO DO DWGM 
    # pre-allocate vectors
    dwgm_x   = zeros(T, n)
    dwgm_x_1 = similar(dwgm_x)

    dwgm_g   = similar(dwgm_x)
    dwgm_g_1 = similar(dwgm_x)

    dwgm_y = similar(dwgm_x)
    dwgm_w = similar(dwgm_x)

    dwgm_aux = similar(dwgm_x)

    # initialize objects
    if !isnothing(x0)
        copyto!(dwgm_x,x0)
    end

    copyto!(dwgm_x_1, dwgm_x)

    copyto!(dwgm_g,b)
    mul!(dwgm_g, A, dwgm_x, 1.0, -1.0)

    copyto!(dwgm_g_1,dwgm_g)

    dwgm_k  = 0
    dwgm_st = 1

    dwgm_gsupn = norm(dwgm_g, Inf)
    # ------------------------------------------------

    # ARMAZENAMENTO DE INFORMAÇÕES ÚTEIS PARA ANÁLISE

    # diferenças com respeito aos iterandos |x_cr - x_dwgm|
    difx = T[0.]
    difg = T[0.]

    # diferenças com respeito ao ótimo |x_k - x^*|
    dwgm_dif_opt = T[norm(dwgm_x .- 1)]
    cr_dif_opt   = T[norm(cr_x .- 1)]

    dwgm_normx = [0.] # ||x||
    cr_normx   = [0.]

    while true
        # test whether gradient norm is small enough
        if cr_rsupn <= prec && dwgm_gsupn <= prec # &&  cr_dif_opt[cr_k+1] <= prec && dwgm_dif_opt[dwgm_k+1] <= prec
            cr_st   = 0
            dwgm_st = 0
            break
        end

        # # test whether the distance of the optimal is small enough
        # if cr_dif_opt[cr_k+1] <= prec && dwgm_dif_opt[dwgm_k+1] <= prec
        #     cr_st   = 0
        #     dwgm_st = 0
        #     break
        # end

        # test whether maximum number of iterations is reached
        if cr_k >= kmax || dwgm_k >= kmax
            cr_st = 1
            dwgm_st = 1
            break
        end        

        if cr_rsupn > prec
            # ------------------------------------------------------
            # CR ITERATION 
            cr_t = cr_rAr/dot(cr_Ad,cr_Ad)

            cr_x .= cr_x .+ cr_t .* cr_d
            cr_r .= cr_r .- cr_t .* cr_Ad

            mul!(cr_Ar, A, cr_r)

            cr_new_rAr = dot(cr_r,cr_Ar)

            cr_beta = cr_new_rAr/cr_rAr

            cr_d .= cr_r .+ cr_beta .* cr_d

            cr_Ad .= cr_Ar .+ cr_beta .* cr_Ad

            cr_rAr = cr_new_rAr

            cr_rsupn = norm(cr_r, Inf)

            cr_k += 1
            # ------------------------------------------------------

            push!(cr_normx, norm(cr_x, 2)/norm(ones(n), 2))
        end

        if dwgm_gsupn > prec
            # ------------------------------------------------------
            # DWGM ITERATION
            mul!(dwgm_w, A, dwgm_g)

            dwgm_alpha = dot(dwgm_g, dwgm_w)/dot(dwgm_w, dwgm_w)

            dwgm_y .= (dwgm_g_1 .- dwgm_g) .+ dwgm_alpha .* dwgm_w

            dwgm_beta = dot(dwgm_g_1, dwgm_y)/dot(dwgm_y, dwgm_y)

            copyto!(dwgm_aux,dwgm_x)
            dwgm_x .= (1.0 - dwgm_beta).*dwgm_x_1 .+ dwgm_beta .* (dwgm_x .- dwgm_alpha .* dwgm_g)
            copyto!(dwgm_x_1,dwgm_aux)

            copyto!(dwgm_aux,dwgm_g)
            dwgm_g .= dwgm_g_1 .- dwgm_beta .* dwgm_y
            copyto!(dwgm_g_1,dwgm_aux)

            dwgm_gsupn = norm(dwgm_g, Inf)

            dwgm_k += 1
            # ------------------------------------------------------

            push!(dwgm_normx, norm(dwgm_x, 2)/norm(ones(n), 2))
        end

        # Armazenando informações importantes
        push!(difx, norm(cr_x - dwgm_x))
        push!(difg, abs(cr_rsupn - dwgm_gsupn))

        push!(dwgm_dif_opt, norm(dwgm_x .- 1))
        push!(cr_dif_opt, norm(cr_x .- 1))
    end

    return dwgm_k, cr_k, difx, difg, dwgm_dif_opt, cr_dif_opt, dwgm_normx, cr_normx
end
