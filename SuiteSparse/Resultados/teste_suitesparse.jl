# Esse código corresponde aos testes para funções quadráticas. 
# Tenho como objetivo testar o desempenho de alguns métodos (DWGM, CR e CG)
# em matflrizes da coletânea da Universidade da órida (SuiteSparseMatrixCollection).
# As matrizes devem ser simétricas e definidas positivas. 

# Os resultados serão armazenados em "resultados.txt" e "resultadosjdl.jld2"

# Carregando pacotes
using LinearAlgebra, MatrixDepot, DataFrames, SparseArrays, JLD2, Printf, CSV, BenchmarkTools

# Switch to Intel MKL BLAS
using MKL

# Configurações adicionais
BLAS.set_num_threads(1)

# Apenas para confirmar as configurações realizadas
# println(LinearAlgebra.BLAS.get_num_threads())

# Carregando os métodos
include("MetodosQuadratica_02-11-2024.jl")


function main()
    # metodos a serem testados
    metodos =["cg"; "cr"; "dwgm"]

    # Precisão para a norma do gradiente (estou utilizando a norma do máximo)
    ϵ = 1e-7

    # Retomada de testes passados
    if isfile("resultados_desktop_ssmc.jdl2")
        @load "resultados_desktop_ssmc.jdl2" resultados
        
        # não sei se isso funciona..(provavelmente...)
        while true
            if size(resultados, 1)%length(metodos) != 0
                delete!(resultados, nrow(resultados))
            else
                break
            end
        end
    else
        # Construção do dataframe resultados
         resultados = DataFrame(Nome=[], n=[], nnz=[], metodo=[], it=[], normg=[], normg_c=[], xnorm=[], st=[], tempo=[])
    end

    # As matrizes abaixo foram classificadas no site da coletânea
    # com as seguintes especificações: posdef, simétrica, real, n >= 20.000 e
    # nnz <= 400.000.000 totalizando 91 matrizes (foram modificados, agora valem 57 = 86 - 24(non solve))
    data = CSV.read("info_matrizes_ssmc.csv", DataFrame)

    sort!(data, :n) # organiza o dataframe em ordem crescente com respeito a coluna n (dimensão)

    # Aplicar filtro nas instâncias (removendo instâncias já resolvidas em teste anterior.)
    # ?? Precisa de modificação caso pare na metade um método (continuação a partir da próxima instância)
    if isfile("resultados_desktop_ssmc.jdl2")
        data = filter(row -> !(row[:grupo_nome] in resultados[:, :Nome]), eachrow(data))
    end
    
    # Pega as informações do dataframe data
    # list_dir = @. data.grupo * "/" * data.nome
    list_dir = data.grupo_nome
    list_n   = data.n
    list_nnz = data.nnz

    # Número de problemas
    nprob = length(list_dir)
    
    try

    for i in 1:nprob
        # Organizando os dados do problema
        name = list_dir[i]
        A    = matrixdepot(name)

        n    = list_n[i]
        b    = A*ones(Float64, n)

        x0   = zeros(Float64, n)
        nz   = list_nnz[i]

        # Número de iterações máximo 
        maxit = 5*n

        for metodo in metodos
            # Chave para a função do método
            mf = getfield(Main, Symbol(metodo))

            # roda uma vez para inicializar e carregar tudo...
            if i == 1
                mf(A, b, kmax=1);
            end

            # Chamando o problema e o método
            tempo = @elapsed ~, normg, normg_c, xnorm, it, st = mf(A, b, kmax= maxit, x0 = x0, eps = ϵ)
    
            # Se o tempo for pequeno faça um Benchmarking
            if tempo < 30.0
                tempo = @belapsed mf(A, b, kmax= maxit, x0 = x0, eps = ϵ) setup=(A = $A; b = $b; maxit = $maxit; x0 = $x0; ϵ = $ϵ; mf = $mf)
            end

            # Salvando informações no DataFrame resultados
            push!(resultados, (name, n, nz, metodo, it, normg, normg_c, xnorm, st, tempo))
            
            # salva resultados em TXT
            txt = open("resultados_desktop_ssmc.txt", "w")
            write(txt, @sprintf("%s", resultados[:, 1:end]))
            close(txt)

            # salva resultados em CSV
            CSV.write("resultados_desktop_ssmc.csv", resultados)

            # você pode carregá-lo executando "@load result_ptoint.jdl2".
            @save "resultados_desktop_ssmc.jdl2" resultados
        end
    end

    catch err
        # Informa o erro e logo após salva as informações
        text = open("erros_execucao.txt", "a"); write(text, "Ocorreu um erro ", err,"\n"); close(text)
    finally
        # salva resultados em TXT
        txt = open("resultados_desktop_ssmc.txt", "w")
        write(txt, @sprintf("%s", resultados[:, 1:end]))
        close(txt)

        # salva resultados em CSV
        CSV.write("resultados_desktop_ssmc.csv", resultados)

        # você pode carregá-lo executando "@load result_ptoint.jdl2".
        @save "resultados_desktop_ssmc.jdl2" resultados
    end

    return 0
end

main()
