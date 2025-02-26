# Autor: Elivandro Grippa
# Data: 26 de Abril de 2024

# estudando erro associado as sequências em um problema da biblioteca SSMC

# Carregando pacotes
using LinearAlgebra, MatrixDepot, DataFrames, SparseArrays, JLD2, Printf, CSV, BenchmarkTools

using MKL
# Configurações adicionais
# BLAS.set_num_threads(1)

# instâncias menores
instancias = "info_matrizes_ssmc.csv"

# Carregando os métodos
include("MetodosQuadratica.jl")

function main()
    # máximo de iterações para os métodos
    ϵ      = 1e-7  # Precisão para a norma do gradiente

    # Carrega o dataframe para salvar os resultados
    resultados = DataFrame(nome=[], n=[], nnz=[], dwgm_k=[], cr_k=[], difx=[], difg=[], 
                           dwgm_dif_opt=[], cr_dif_opt=[], dwgm_normx=[], cr_normx = [])


    # Informações da Matriz selecionada
    data = CSV.read(instancias, DataFrame)
    
    sort!(data, :nnz) # organiza o dataframe em ordem crescente com respeito a coluna nnz

    # Pega as informações do dataframe data
    list_dir = @. data.grupo * "/" * data.nome
    list_n   = data.n
    list_nnz = data.nnz
    
    i = 1
    # Organizando os dados do problema
    name = list_dir[i]
    A    = matrixdepot(name)

    n    = list_n[i]
    b    = A*ones(n)

    # x0   = zeros(n)
    nz   = list_nnz[i]

    maxit = 3_000

    dwgm_k, cr_k, difx, difg, dwgm_dif_opt, cr_dif_opt, dwgm_normx, cr_normx = cr_dwgm(A, b, kmax= maxit, prec = ϵ)

    push!(resultados, (name, n, nz, dwgm_k, cr_k, difx, difg, dwgm_dif_opt, cr_dif_opt, dwgm_normx, cr_normx))

    # você pode carregá-lo executando "@load result_ptoint.jdl2".
    @save "resultadosjdl_desktop_$(instancias[15:end-4]).jdl2" resultados

    return resultados
end


# Rotina para cálculo das distâncias
main()


# Construção das Figuras
include("pp.jl")
pp_normx(save = true)
pp_difx(save = true)
