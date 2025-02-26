# Autor: Elivandro
# Data: 17 de março de 2024

# Construção de perfis de desempenho (Perfomance Profile)

using LinearAlgebra, DataFrames, JLD2, Plots, BenchmarkProfiles, LaTeXStrings, PGFPlotsX

pgfplotsx()
# gr()

# Carregando informações do Dataframe com os resultados
@load "resultadosjdl_desktop_ssmc.jdl2"

# Gráfico: |g_cr - g_dwgm|
function pp_difg(; save = false)
    pl   = resultados[1, :difg]
    name = split(resultados[1, :nome], "/")[2]

    fig = plot(xlabel=LaTeXString("\$ k \$"),
               ylabel=LaTeXString("\$ \\||g_k^{\\text{CR}}| - |g_k^{\\text{DWGM}}| \\| \$"),
               legend=:topright,
               lw=1.3,
               guidefont=font(15),
               tickfont=font(11),
               legendfont=font(12)
            )

    K1 = length(pl)
    x = 3:K1
    
    plot!(x, pl[3:K1], label= nothing)

    	
	if save
        savefig(fig, "Figuras/plot_difg-DWGMxCR - $name.pdf")
    else
        return fig
    end
end

# Gráfico ||x_cr - x_dwgm||
function pp_difx(; save = false)
    pl   = resultados[1, :difx]
    name = resultados[1, :nome]

    k_cr   = resultados[1, :cr_k]
    k_dwgm = resultados[1, :dwgm_k]

    name = split(resultados[1, :nome], "/")[2]

    fig = plot(xlabel=LaTeXString("\$ k \$"),
               ylabel=LaTeXString("\$ \\|x_k^{CR} - x_k^{DWGM} \\|_2 \$"),
            #    title = LaTeXString("Problema $name"),
               legend=:topright,
               lw=1.7,
               yscale=:log10,
               guidefont=font(15),
               tickfont=font(11),
               legendfont=font(12)
            #    xscale=:log10
            #    size = (700,500)
            )


    K1 = max(k_cr, k_dwgm)
    k = 3:K1
    plot!(k, pl[3:K1], lw=1.5, label=nothing)


	if save
        savefig(fig, "Figuras/plot_difx-DWGMxCR - $(name).pdf")
    else
        return fig
    end
end

# gráfico: ||x_k - x^*||
function pp_normx(; save = false)
    pl_cr   = resultados[1, :cr_dif_opt]
    pl_dwgm = resultados[1, :dwgm_dif_opt]

    k_cr   = resultados[1, :cr_k]
    k_dwgm = resultados[1, :dwgm_k]

    name = split(resultados[1, :nome], "/")[2]

    lab = [
        "DWGM"
        "CR"    
    ]

    fig = plot(xlabel=LaTeXString("\$ k \$"),
               ylabel=LaTeXString("\$ \\| x_k - \\bar x\\|_2 \$"),
            #    title = LaTeXString("Problema $name"),
            legend=:topright,
            yscale=:log10,
            lw  =1.7,
            guidefont=font(15),
            tickfont=font(11),
            legendfont=font(12)
            #    xscale=:log10
            #    size = (700,500)
        )


        # para o dwgm
        K1 = k_dwgm
        x = 1:K1
        plot!(x, pl_dwgm[1:K1], yscale=:log10, label = lab[1], lw=1.5)
    
        # para o cr 
        K2 = k_cr
        x = 1:K2
        plot!(x, pl_cr[1:K2], yscale=:log10, label = lab[2], lw=1.5)

    	
	if save
        savefig(fig, "Figuras/plot_normx-DWGMxCR - $name.pdf")
    end

    return fig
end

# gráfico: ||x_k||/||x^*||
function pp_normxtwo(; save = false)
    pl_cr   = resultados[1, :cr_normx]
    pl_dwgm = resultados[1, :dwgm_normx]

    k_cr   = resultados[1, :cr_k]
    k_dwgm = resultados[1, :dwgm_k]

    name = split(resultados[1, :nome], "/")[2]

    lab = [
        "DWGM"
        "CR"    
    ]

    fig = plot(xlabel=LaTeXString("\$ k \$"),
                ylabel=LaTeXString("\$ \\| x_k \\| / \\| \\bar x \\| \$"),
                legend=:bottomright,#:topright,
                yscale=:log10,
                lw  =1.3,
                guidefont=font(15),
                tickfont=font(11),
                legendfont=font(12)
        )


        # para o dwgm
        K1 = min(5000, k_dwgm)
        x = 2:K1

        plot!(x, pl_dwgm[2:K1],
            yscale=:log10,
            xscale=:log10,
            label = lab[1],
        )
    
        # para o cr 
        K2 = min(5000, k_cr)
        x = 2:K2

        plot!(x, pl_cr[2:K2],
            yscale=:log10,
            xscale=:log10,
            label = lab[2],
        )

    	
	if save
        savefig(fig, "Figuras/plot_normx2-DWGMxCR - $name.pdf")
    end

    return fig
end

# pp_normx(save = true)
# pp_difx(save = true)

# pp_normx()
# pp_normxtwo()

