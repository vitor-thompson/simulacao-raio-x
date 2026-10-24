# Simulação Computacional de Raios-X em C++ (Monte Carlo)

Ferramenta de simulação em física computacional desenvolvida em **C++** para modelagem espectral de raios-X e transporte estocástico de fótons utilizando o **Método de Monte Carlo**.

---

## 🚀 O que o código faz?

1. **Leitura de Bases de Dados Atômicas**: Lê arquivos de configuração de elementos químicos (como Prata, Carbono, etc.) contendo propriedades fundamentais como número atômico ($Z$), massa atômica ($A$), densidade, energias de ligação dos subníveis eletrônicos e coeficientes de interação.
2. **Modelagem de Tubos de Raios-X (`tubo_rx`)**: 
   - Calcula a componente de **radiação contínua (*Bremsstrahlung*)**.
   - Calcula a componente de **radiação discreta (picos característicos de fluorescência XRF do anodo)**.
   - Gera o espectro completo de emissão do tubo e constrói uma Função de Distribuição Acumulada (CDF) para amostragem estocástica.
3. **Simulação de Interações da Radiação na Matéria (Monte Carlo)**:
   - Sorteia a energia dos fótons emitidos pelo tubo respeitando o espectro real.
   - Calcula as seções de choque de **espalhamento elástico (Rayleigh)**, **espalhamento inelástico (Compton)** e **absorção fotoelétrica**.
   - Determina probabilisticamente o destino de cada fóton ao atravessar uma espessura de material alvo (Transmissão, Rayleigh, Compton ou Fotoelétrica).
4. **Exportação de Dados**: Gera arquivos de saída com o espectro simulado (`espectro.txt`).

---

## 💻 Ambiente de Execução

O código é otimizado e configurado para compilação e execução nativa no ambiente **Microsoft Visual Studio (MSVC)** no Windows.

---

## ⚙️ Configuração do Visual Studio (Passo a Passo para o `dbase`)

Como o programa utiliza argumentos de linha de comando para localizar a pasta com os dados atômicos e definir o anodo e o alvo, siga os passos abaixo para configurar o projeto no Visual Studio:

### 1. Organização dos Arquivos
1. Abra a pasta raiz do seu projeto no Explorador de Arquivos do Windows (onde está o arquivo `.cpp` principal).
2. Cole a pasta **`dbase`** (contendo os arquivos `.txt` dos elementos, ex: `Ag.txt`, `C.txt`) diretamente dentro dessa pasta do projeto.

### 2. Configuração dos Argumentos no Visual Studio
Para que o programa saiba onde encontrar a pasta `dbase` e quais elementos simular:

1. No menu superior do Visual Studio, vá em **Projeto** > **Propriedades** (ou pressione `Alt + F7`).
2. No menu à esquerda, expanda **Propriedades de Configuração** e clique em **Depuração** (*Debugging*).
3. No campo **Argumentos do Comando** (*Command Arguments*), insira os parâmetros:
   ```text
   dbase Ag C
(Onde dbase é a pasta de dados, Ag é o anodo e C é o material alvo).
4. Clique em Aplicar e depois em OK.

Pronto! Agora basta compilar (Ctrl + Shift + B) e executar com depuração (F5).
