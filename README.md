# Comparação entre Métodos Iterativos para a minimização de quadráticas estritamente convexas

Este repositório contém um compilado dos principais testes comparando três métodos iterativos para problemas quadráticos:

- **Método de Gradientes Conjugados (CG)** (com parâmetro de conjugação [Fletcher-Reeves](https://doi.org/10.1093/comjnl/7.2.149))
- **Método de Resíduos Conjugados (CR)**  (proposto por [Stiefel](https://dx.doi.org/10.1007/bf02564277))
- **Método do Gradiente Ponderado com Atraso** (proposto por [Luenberguer](https://doi.org/10.1137/0707032)-[Oviedo](https://dx.doi.org/10.1007/s10589-019-00125-6))

Os resultados foram resumidos em tabelas e perfis de desempenho para as diferentes classes de matrizes.

## 📌 Resultados Apresentados  
Os experimentos foram realizados em duas classes de matrizes:  
1. **Matrizes Sintéticas** (utilizadas por [Oviedo](https://dx.doi.org/10.1007/s10589-019-00125-6))
2. **Matrizes da coletânea [SuiteSparse](https://sparse.tamu.edu/)**  

Os resultados foram apresentados na defesa de mestrado em **17 de fevereiro de 2025**.  

## 📊 Perfis de Desempenho  
Os perfis de desempenho foram gerados para avaliar a eficiência dos métodos em diferentes cenários.  

## 📁 Organização do Repositório  
- `Sintéticas/` - Contém os códigos e resultados correspondentes as matrizes sintéticas.  
- `SuiteSparse/`- Contém os códigos e resultados correspondentes as matrizes da coletânea _SuiteSparse_.   
- `Figuras/` - Contém os perfis de desempenho para os problemas gerais, estruturais e não estruturais.
  
## 🚀 Como Reproduzir os Testes  
1. Clone este repositório:  
   ```bash
   git clone https://github.com/elivandrogrippa/Qualificacao.git
   cd Qualificacao]
2. Execute os códigos iniciados em "Testes" no prompt Julia. 
