# Comparação de Métodos Iterativos para a minimização de quadráticas estritamente convexas

Este repositório contém um compilado dos principais testes comparando três métodos iterativos para problemas quadráticos:

- **Método dos Gradientes Conjugados (CG)**
- **Método dos Resíduos Conjugados (CR)**  
- **Método do Gradiente Ponderado com Atraso**  

Os resultados incluem análises detalhadas e perfis de desempenho para diferentes classes de matrizes.

## 📌 Resultados Apresentados  
Os experimentos foram realizados em duas classes de matrizes:  
1. **Matrizes Sintéticas**  
2. **Matrizes da coletânea SuiteSparse**  

Os resultados foram apresentados na defesa de mestrado em **17 de fevereiro de 2025**.  

## 📊 Perfis de Desempenho  
Os perfis de desempenho foram gerados para avaliar a eficiência dos métodos em diferentes cenários.  

## 📁 Organização do Repositório  
- `dados/` - Contém as instâncias das matrizes utilizadas.  
- `scripts/` - Códigos para execução dos experimentos.  
- `resultados/` - Saídas dos testes e gráficos dos perfis de desempenho.  
- `README.md` - Este documento.  

## 🚀 Como Reproduzir os Testes  
1. Clone este repositório:  
   ```bash
   git clone https://github.com/seu-usuario/seu-repositorio.git
   cd seu-repositorio
