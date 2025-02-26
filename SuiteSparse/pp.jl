# Construção de perfis de desempenho (Performance Profile)

# As versões recentes do pacote BenchmarkProfiles apresenta problemas
# quando utilizado em conjunto com pgfplotsx().
# Utilizar: BenchmarkProfiles v0.3.4 (a escala fica melhor)

using LinearAlgebra, DataFrames, JLD2, Plots, PGFPlotsX, BenchmarkProfiles, LaTeXStrings
using CSV

pgfplotsx()
# gr()

# Separa problemas irrelevantes
function Separa_irrelevantes(resultados; saida = true)
    n = size(resultados, 1)

    # Problemas desconsiderados - contador
    num_prob = 0
    list_linha = Int64[]

    for i in 1:3:n
        cg_st = resultados[i, :st]
        cr_st = resultados[i+1, :st] 
        dwgm_st =  resultados[i+2, :st]

        if dwgm_st * cg_st * cr_st >= 1
            push!(list_linha, i, i+1, i+2)

            # println("problema: ", resultados[i, :Nome])
            num_prob = num_prob + 1
        end
    end

    delete!(resultados, list_linha)

    if saida
        println("\nTotal de problemas: ", Int64(n/3))
        println("Problemas desconsiderados: ",num_prob)
    end
end

# Corrige o Status do problema (caso necessário)
function corrigindoStCriterio(resultados)
	n = size(resultados, 1)

	for i in 1:3:n
        prec = 1e-7
		linha_cg = resultados[i, :]
		linha_cr = resultados[i+1, :]
        linha_dwgm = resultados[i+2, :]
		
		if linha_cg[:normg] < prec && linha_cg[:xnorm] < prec
 			resultados[i, :st] = 0
		else
			resultados[i, :st] = 1
		end
		
		if linha_cr[:normg] < prec && linha_cr[:xnorm] < prec
			resultados[i+1, :st] = 0
		else
			resultados[i+1, :st] = 1
		end

        if linha_dwgm[:normg] < prec && linha_dwgm[:xnorm] < prec
			resultados[i+2, :st] = 0
		else
			resultados[i+2, :st] = 1
        end
	end
end

# Separa as informações do DataFrame resultados 
function result_data(resultados)
	n = size(resultados, 1)
	
	tempo = Inf * ones(Int(n/3), 3)
	iter  = Inf * ones(Int(n/3), 3)

	k = 1
	for i in 1:3:n
		li_cg   = resultados[i,:]
		li_cr   = resultados[i+1,:]
		li_dwgm = resultados[i+2,:]

		
		if li_cg[:st] == 0
			tempo[k, 1] = li_cg[:tempo]
			iter[k, 1]  = li_cg[:it]
		end

		if li_cr[:st] == 0
			tempo[k, 2] = li_cr[:tempo]
			iter[k, 2]  = li_cr[:it]
		end

		if li_dwgm[:st] == 0
			tempo[k, 3] = li_dwgm[:tempo]
			iter[k, 3]  = li_dwgm[:it]
		end
		
		k = k + 1
	end

	return tempo, iter
end

# Construção dos perfis com carregamento e organização das informações - todas as instancias
function pp_CriterioHibrido(;save=false)

    # Carregando informações do Dataframe com os resultados
    @load "Resultados/resultados_desktop_ssmc.jdl2"

    # Remove os problemas que nenhum dos métodos resolveu...
    Separa_irrelevantes(resultados)

    # Construção dos dados para os perfis (iterações e tempo de CPU)
    data_tempo, data_iter = result_data(resultados);

    fig = performance_profile(
        PlotsBackend(),
        data_tempo,
        ["CG","CR" ,"DWGM"],
        lw = 1.7,
        # title = "Perfil de desempenho - tempo de CPU",
        xlabel = latexstring("\$\$Dentro deste fator do melhor (escala logarítmica)"),
        ylabel = latexstring("Proporção de problemas \$\$"),
        legend = :bottomright,
        # titlefont=font(17),
        guidefont=font(15),
        tickfont=font(11),
        legendfont=font(12),
        ytick=0:0.2:1.0,
        # ytick=(0:0.2:1, ["\$ $i \\% \$" for i in 0:20:100])
        linestyles=[:solid, :dash, :dot]
    )
	
    println("\nObs.: Relações apresentadas da escala logarítmica!!\n")
    println("Com relação ao perfil de tempo: ")

    (ratios, max_ratio) = performance_ratios(data_tempo, logscale=true, drawtol=0.0)
    (np, ns) = size(ratios)
    
    function rho(s, tau)
        cont =  sum(Float64.(ratios[:, s] .<= tau)) / np
        return round(cont, digits=2)    
    end
    
 
    println("\nCG   : τ(0) = $(rho(1, 0)) e τ(max) = ", rho(1, max_ratio))
    println("CR   : τ(0) = $(rho(2, 0)) e τ(max) = ", rho(2, max_ratio))
    println("DWGM : τ(0) = $(rho(3, 0)) e τ(max) = ", rho(3, max_ratio))

	if save
        savefig(fig, "Figuras/pp_tempo-CGxCRxDWGM.pdf")
    else 
        return fig
	end
    

    println("\nCom relação ao perfil de iterações: ")
    
    (ratios, max_ratio) = performance_ratios(data_iter, logscale=true, drawtol=0.0)
    
    println("\nCG   : τ(0) = $(rho(1, 0)) e τ(max) = ", rho(1, max_ratio))
    println("CR   : τ(0) = $(rho(2, 0)) e τ(max) = ", rho(2, max_ratio))
    println("DWGM : τ(0) = $(rho(3, 0)) e τ(max) = ", rho(3, max_ratio))
    
    fig = performance_profile(
        PlotsBackend(),
        data_iter,
        ["CG","CR" ,"DWGM"],
        lw = 1.7,
        # title = latexstring("Perfil de desempenho - número de iterações\$\$"),
        xlabel = "Dentro deste fator do melhor (escala logarítmica)",
        ylabel = latexstring("Proporção de problemas \$\$"),
        legend = :bottomright,
        # titlefont=font(17),
        guidefont=font(15),
        tickfont=font(11),
        legendfont=font(12),
        ytick=0:0.2:1.0,
        linestyles=[:solid, :dash, :dot]
    )
    
	if save
  		savefig(fig, "Figuras/pp_iter-CGxCRxDWGM.pdf")
    else
        return fig
	end

    return 0
end

# Computa média geométrica para os dois métodos mais competitivos (seguindo os perfis)
function Media_geometrica_CGxCR()

    # Carregando informações do Dataframe com os resultados
    @load "Resultados/resultados_desktop_ssmc.jdl2"

    # Remove os problemas que nenhum dos métodos resolveu...
    Separa_irrelevantes(resultados, saida = false)

    p = size(resultados, 1)

	cost  = zeros(Int(p/3))

	k = 1
	for i in 1:3:p
		li_cg = resultados[i,:]
		li_cr = resultados[i+1,:]

		tempo_cg = li_cg[:tempo]
		tempo_cr = li_cr[:tempo]

        cost[k]  = tempo_cg/tempo_cr
		
		k = k + 1
	end
	
	return println("\nMG das razões (t_cg/t_cr) = ",prod(cost)^(1/length(cost)))
end

# função auxiliar, retorna o tipo do problema dado o nome da instância
function info_matrizes(nome)
    data = CSV.read("Resultados/info_matrizes_ssmc.csv", DataFrame)

    result = filter(row -> row.grupo_nome == nome, data)

    return result[1, :tipo]
end

# Filtro em problemas estruturais utilizados
function filtro_estruturais_resultados(resultados)
    resultados[:, :tipo] = info_matrizes.(resultados[:, :Nome])

    # filtra com relação a problemas estruturais
    filter!(row -> row.tipo == "Structural Problem" || row.tipo == "Subsequent Structural Problem", resultados)
end

# Construção dos perfis com relação aos problemas de classe estrutural
function pp_Structural(;save=false)

    # Carregando informações do Dataframe com os resultados
    @load "Resultados/resultados_desktop_ssmc.jdl2"

    # Remove os problemas que nenhum dos métodos resolveu...
    Separa_irrelevantes(resultados)
    filtro_estruturais_resultados(resultados)

    # Construção dos dados para os perfis (iterações e tempo de CPU)
    data_tempo, data_iter = result_data(resultados);

    fig = performance_profile(
        PlotsBackend(),
        data_tempo,
        ["CG","CR" ,"DWGM"],
        lw = 1.7,
        # title = "Perfil de desempenho - tempo de CPU",
        xlabel = latexstring("\$\$Dentro deste fator do melhor (escala logarítmica)"),
        ylabel = latexstring("Proporção de problemas \$\$"),
        legend = :bottomright,
        # titlefont=font(17),
        guidefont=font(15),
        tickfont=font(11),
        legendfont=font(12),
        ytick=0:0.2:1.0,
        # ytick=(0:0.2:1, ["\$ $i \\% \$" for i in 0:20:100])
        linestyles=[:solid, :dash, :dot]
    )
	
    println("\nObs.: Relações apresentadas da escala logarítmica!!\n")
    println("Com relação ao perfil de tempo: ")

    (ratios, max_ratio) = performance_ratios(data_tempo, logscale=true, drawtol=0.0)
    (np, ns) = size(ratios)
    
    function rho(s, tau)
        cont =  sum(Float64.(ratios[:, s] .<= tau)) / np
        return round(cont, digits=2)    
    end
    
    println("\nCG   : τ(0) = $(rho(1, 0)) e τ(max) = ", rho(1, max_ratio))
    println("CR   : τ(0) = $(rho(2, 0)) e τ(max) = ", rho(2, max_ratio))
    println("DWGM : τ(0) = $(rho(3, 0)) e τ(max) = ", rho(3, max_ratio))

	if save
        savefig(fig, "Figuras/pp_tempo-CGxCRxDWGM-Structural.pdf")
    else 
        return fig
	end

    println("\nCom relação ao perfil de iterações: ")
    
    (ratios, max_ratio) = performance_ratios(data_iter, logscale=true, drawtol=0.0)
    
    println("\nCG   : τ(0) = $(rho(1, 0)) e τ(max) = ", rho(1, max_ratio))
    println("CR   : τ(0) = $(rho(2, 0)) e τ(max) = ", rho(2, max_ratio))
    println("DWGM : τ(0) = $(rho(3, 0)) e τ(max) = ", rho(3, max_ratio))
    
    fig = performance_profile(
        PlotsBackend(),
        data_iter,
        ["CG","CR" ,"DWGM"],
        lw = 1.7,
        # title = latexstring("Perfil de desempenho - número de iterações\$\$"),
        xlabel = "Dentro deste fator do melhor (escala logarítmica)",
        ylabel = latexstring("Proporção de problemas \$\$"),
        legend = :bottomright,
        # titlefont=font(17),
        guidefont=font(15),
        tickfont=font(11),
        legendfont=font(12),
        ytick=0:0.2:1.0,
        linestyles=[:solid, :dash, :dot]
    )
    
	if save
  		savefig(fig, "Figuras/pp_iter-CGxCRxDWGM-Structural.pdf")
    else
        return fig
	end
end

# Computa média geométrica para os dois métodos mais competitivos (segundo os perfis)
function Media_geometrica_Structural_CGxCR()

    # Carregando informações do Dataframe com os resultados
    @load "Resultados/resultados_desktop_ssmc.jdl2"

    # Remove os problemas que nenhum dos métodos resolveu...
    Separa_irrelevantes(resultados, saida = false)
    filtro_estruturais_resultados(resultados)

    p = size(resultados, 1)

	cost  = zeros(Int(p/3))

	k = 1
	for i in 1:3:p
		li_cg = resultados[i,:]
		li_cr = resultados[i+1,:]

		tempo_cg = li_cg[:tempo]
		tempo_cr = li_cr[:tempo]

        cost[k]  = tempo_cr/tempo_cg
		
		k = k + 1
	end
	
	return println("\nMG das razões (t_cr/t_cg) (Structural) = ",prod(cost)^(1/length(cost)))
end

# Media_geometrica_Structural_CGxCR()

# pp_CriterioHibrido_Structural(save=true)

# pp_CriterioHibrido(save=true)

# Media_geometrica_CGxCR()

# Filtro em problemas não estruturais utilizados
function filtro_NonEstruturais_resultados(resultados)
    resultados[:, :tipo] = info_matrizes.(resultados[:, :Nome])

    # filtra com relação a problemas estruturais
    filter!(row -> row.tipo != "Structural Problem" && row.tipo != "Subsequent Structural Problem", resultados)
end

# Construção dos pefis com carregamento e organização das informações
# Nessa função analisamos os problemas não estruturais
# Busca por padrões
function pp_NonStructural(;save=false)

    # Carregando informações do Dataframe com os resultados
    @load "Resultados/resultados_desktop_ssmc.jdl2"

    # Remove os problemas que nenhum dos métodos resolveu...
    Separa_irrelevantes(resultados)
    filtro_NonEstruturais_resultados(resultados)

    # Construção dos dados para os perfis (iterações e tempo de CPU)
    data_tempo, data_iter = result_data(resultados);

    fig = performance_profile(
        PlotsBackend(),
        data_tempo,
        ["CG","CR" ,"DWGM"],
        lw = 1.7,
        # title = "Perfil de desempenho - tempo de CPU",
        xlabel = latexstring("\$\$Dentro deste fator do melhor (escala logarítmica)"),
        ylabel = latexstring("Proporção de problemas \$\$"),
        legend = :bottomright,
        # titlefont=font(17),
        guidefont=font(15),
        tickfont=font(11),
        legendfont=font(12),
        ytick=0:0.2:1.0,
        # ytick=(0:0.2:1, ["\$ $i \\% \$" for i in 0:20:100])
        linestyles=[:solid, :dash, :dot]
    )
	
    println("\nObs.: Relações apresentadas da escala logarítmica!!\n")
    println("Com relação ao perfil de tempo: ")

    (ratios, max_ratio) = performance_ratios(data_tempo, logscale=true, drawtol=0.0)
    (np, ns) = size(ratios)
    
    function rho(s, tau)
        cont =  sum(Float64.(ratios[:, s] .<= tau)) / np
        return round(cont, digits=2)    
    end
    
 
    println("\nCG   : τ(0) = $(rho(1, 0)) e τ(max) = ", rho(1, max_ratio))
    println("CR   : τ(0) = $(rho(2, 0)) e τ(max) = ", rho(2, max_ratio))
    println("DWGM : τ(0) = $(rho(3, 0)) e τ(max) = ", rho(3, max_ratio))

	if save
        savefig(fig, "Figuras/pp_tempo-CGxCRxDWGM-NonStructural.pdf")
    else 
        return fig
	end
    

    println("\nCom relação ao perfil de iterações: ")
    
    (ratios, max_ratio) = performance_ratios(data_iter, logscale=true, drawtol=0.0)
    
    println("\nCG   : τ(0) = $(rho(1, 0)) e τ(max) = ", rho(1, max_ratio))
    println("CR   : τ(0) = $(rho(2, 0)) e τ(max) = ", rho(2, max_ratio))
    println("DWGM : τ(0) = $(rho(3, 0)) e τ(max) = ", rho(3, max_ratio))
    
    fig = performance_profile(
        PlotsBackend(),
        data_iter,
        ["CG","CR" ,"DWGM"],
        lw = 1.7,
        # title = latexstring("Perfil de desempenho - número de iterações\$\$"),
        xlabel = "Dentro deste fator do melhor (escala logarítmica)",
        ylabel = latexstring("Proporção de problemas \$\$"),
        legend = :bottomright,
        # titlefont=font(17),
        guidefont=font(15),
        tickfont=font(11),
        legendfont=font(12),
        ytick=0:0.2:1.0,
        linestyles=[:solid, :dash, :dot]
    )
    
	if save
  		savefig(fig, "Figuras/pp_iter-CGxCRxDWGM-NonStructural.pdf")
    else
        return fig
	end
end

# Computa média geométrica para os dois métodos mais competivos (segundo os pefis)
function Media_geometrica_NonStructural_CGxCR()

    # Carregando informações do Dataframe com os resultados
    @load "Resultados/resultados_desktop_ssmc.jdl2"

    # Remove os problemas que nenhum dos métodos resolveu...
    Separa_irrelevantes(resultados, saida = false)
    filtro_NonEstruturais_resultados(resultados)

    p = size(resultados, 1)

	cost  = zeros(Int(p/3))

	k = 1
	for i in 1:3:p
		li_cg = resultados[i,:]
		li_cr = resultados[i+1,:]

		tempo_cg = li_cg[:tempo]
		tempo_cr = li_cr[:tempo]

        cost[k]  = tempo_cg/tempo_cr
		
		k = k + 1
	end
	
	return println("\nMG das razões (t_cg/t_cr) (NonStructural) = ",prod(cost)^(1/length(cost)))
end

# pp_NonStructural(save=true)
# Media_geometrica_NonStructural_CGxCR()


# instâncias não consideradas na tabela da dissertação - auxilio na escrita
function InstanciasNonConsideredInfo()
    data = CSV.read("Resultados/info_matrizes_ssmc.csv", DataFrame)

    sort!(data, :n) # organiza o dataframe em ordem crescente com respeito a coluna n (dimensão)

    # Aplicar filtro nas instâncias (removendo instâncias já resolvidas em teste anterior.)
    # ?? Precisa de modificação caso pare na metade um método (continuação a partir da próxima instância)
    if isfile("Resultados/resultados_desktop_ssmc.jdl2")
        
        @load "Resultados/resultados_desktop_ssmc.jdl2" resultados

        data = filter(row -> !(row[:grupo_nome] in resultados[:, :Nome]), eachrow(data))
    end
end

# Analisaremos o quanto os métodos erram o gradiente na última iteração
# Desconsiderei as instâncias que não foram resolvidas por nenhum dos métodos,
# mas podem ser consideradas a depender da análise. 
function erro_gradiente()
    # Carregando informações do Dataframe com os resultados
    @load "Resultados/resultados_desktop_ssmc.jdl2"

    # Remove os problemas que nenhum dos métodos resolveu...
    Separa_irrelevantes(resultados)

    resultados[:, :erro_grad] .= abs.(resultados[:, :normg] - resultados[:, :normg_c])

    resul_cg = filter(row -> row.metodo == "cg", resultados)
    resul_cr = filter(row -> row.metodo == "cr", resultados)
    resul_dwgm = filter(row -> row.metodo == "dwgm", resultados)

    filter!(row -> !isnan(row.erro_grad), resul_dwgm)
    
    data_cg, data_cr, data_dwgm = describe(resul_cg),  describe(resul_cr),  describe(resul_dwgm)

    data = DataFrame(metodo=[], media=[], min=[], mediana=[], max=[])

    push!(data, ("cg", data_cg[11, 2], data_cg[11, 3], data_cg[11, 4], data_cg[11, 5]))
    push!(data, ("cr", data_cr[11, 2], data_cr[11, 3], data_cr[11, 4], data_cr[11, 5]))
    push!(data, ("dwgm", data_dwgm[11, 2], data_dwgm[11, 3], data_dwgm[11, 4], data_dwgm[11, 5]))

    return data 

end

# data = erro_gradiente()



