
# Nome: Elivandro oliveira grippa
# 23 de abril de 2024

# principais pacotes utilizados
using SparseArrays, Printf, DataFrames, LinearAlgebra, JLD2, BenchmarkTools

# Carregando os métodos para quadráticas
include("MetodosQuadratica.jl")

# Switch to Intel MKL BLAS
using MKL

# Configurações adicionais
BLAS.set_num_threads(1)

# Apenas para confirmar as configurações realizadas
# println(LinearAlgebra.BLAS.get_num_threads())


# Instâncias propostas no paper inicial (oviedo leon)
function GerandoInstancia(n, ncond)
    # eye = I(n)

    eye = diagm(ones(n))
    Q   = similar(eye)

    copyto!(Q, eye)

    # Determinando vetores e normalizando.
    v1 = randn(n, 1); v1 = v1 ./ norm(v1, 2)
    v2 = randn(n, 1); v2 = v2 ./ norm(v2, 2)
    v3 = randn(n, 1); v3 = v3 ./ norm(v3, 2)

    # Computando a matriz Q
    Q = (eye - 2 * v1 * v1') * (eye - 2 * v2 * v2') * (eye - 2 * v3 * v3')

    d = zeros(n)

    # cond(A) = exp(ncond)
    d = [exp(((i - 1)/(n - 1)) * ncond) for i in 1:n]
    D = diagm(d)

    A = Q * D * Q'

    x = 2 * rand(n, 1) - ones(n, 1)
    b = A * x

    return A, b
end

# TESTES PARA A DISSERTAÇÃO DE MESTRADO...

# Função teste -  principal
function teste_quadraticas()

    # Carrega o dataframe para salvar os resultados
    resultados = DataFrame(metodo=[], n=[], ncond=[], it=[], normg=[], st=[], tempo=[])

    # máximo de iterações para os métodos
    maxit = 50_000 # 5 * 1e4
    eps = 1e-8             # Precisão para a norma do gradiente

    # metodos a serem testados
    metodos =[
        "dwgm"
        "cg"
        "cr"
    ]

    for n in [1000; 5000; 10000; 20000]
        for ncond in [5; 10; 15]

            A, b = GerandoInstancia(n, ncond)
            x0   = zeros(n)
            
            # println("Resolvendo o problema para n = $n e ncond = $ncond")
            
            for metodo in metodos
                # Chave para a função do método
                mf = getfield(Main, Symbol(metodo))

                # println("*Com o método: $metodo*")

                # Chamando o problema e o método
                tempo = @elapsed ~, normg, it, st, ~, ~ = mf(A, b, kmax= maxit, x0 = x0, eps = eps, history = false)
                
                # Se o tempo for muito pequena faça um Benchmark...
                if tempo < 30.0
                    tempo = @belapsed mf(A, b, kmax= maxit, x0 = x0, eps = eps) setup=(A = $A; b = $b; maxit = $maxit; x0 = $x0; eps= $eps; mf = $mf)
                end

                push!(resultados, (metodo, n, ncond, it, normg, st, tempo))

                # salva resultados para acompanharmos...
                txt = open("test3_resultados_desktop_small.txt", "w"); write(txt, @sprintf("%s", resultados)); close(txt)
            end
        end
    end

    @save "test3_dataframe_desktop_small.jdl2" resultados
    
    txt = open("test3_resultados_desktop_small.txt", "w")
    write(txt, @sprintf("%s", resultados))
    close(txt)

    return resultados
end

teste_quadraticas()
