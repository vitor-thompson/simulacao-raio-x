#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

//---------------------------------------------------class fluorescencia-----------------------------

// Classe que armazena os dados de transição de fluorescência de raios X (XRF)
class fluorescencia {
private:
  char   linha[6];      // Nome da linha de emissão XRF (ex: K-L3)
  double prob;          // Probabilidade individual da transição
  double proba;         // Probabilidade acumulada da transição
  double energia;       // Energia do fóton emitido pela transição (em keV)

public:
  // Métodos Getters e Setters para manipulação dos atributos privados da fluorescência
  void          Linha (const char* nome)   { strcpy(linha, nome); } // Copia o nome do parâmetro para a variável membro linha
  char*         Linha (void)               { return linha; }        // Retorna o ponteiro para a string com o nome da linha
  
  void          Prob  (double e)           { prob = e; }            // Atribui o valor da probabilidade individual
  double        Prob  (void)               { return prob; }         // Retorna o valor da probabilidade individual
  
  void          ProbA (double e)           { proba = e; }           // Atribui o valor da probabilidade acumulada
  double        ProbA (void)               { return proba; }        // Retorna o valor da probabilidade acumulada
  
  void          E     (double e)           { energia = e; }         // Atribui a energia da linha de fluorescência
  double        E     (void)               { return energia; }      // Retorna a energia da linha de fluorescência
};                        

//----------------------------------------------------class subnivel-----------------------------------

// Classe que representa um subnível eletrônico atômico (ex: K, L1, L2, L3...)
class subnivel {
private:
  char      simb[6];    // Nome/símbolo do subnível
  double    Eb;         // Energia de ligação (Binding Energy - Eb) em keV
  double    rs;         // Raio de espalhamento / fator de forma atômico
  double    w;          // Rendimento de fluorescência ou largura de linha
  double    coef[6];    // Coeficientes do polinômio para a seção de choque fotoelétrica
  
  int Nf;               // Número total de linhas de fluorescência associadas a este subnível
  fluorescencia *xrf;   // Ponteiro para o array dinâmico de objetos da classe fluorescencia
  
public: 
  // Getters e Setters das propriedades do subnível
  void          Simb   (const char *nome)  { strcpy(simb, nome); } // Copia o nome do subnível
  const char*   Simb   (void)              { return simb; }        // Retorna o nome do subnível
  
  void          GetEb  (double e)          { Eb = e; }             // Atribui a energia de ligação Eb
  double        GetEb  (void)              { return Eb; }          // Retorna a energia de ligação Eb
  
  void          Rs     (double r)          { rs = r; }             // Atribui o raio de espalhamento/fator de forma
  double        Rs     (void)              { return rs; }          // Retorna o raio de espalhamento/fator de forma
  
  void          W      (double e)          { w = e; }              // Atribui a largura/rendimento de fluorescência
  double        W      (void)              { return w; }           // Retorna a largura/rendimento de fluorescência
  
  // Define a quantidade de linhas de fluorescência e aloca dinamicamente a memória necessária
  void          GetNf (int n)              { Nf = n; xrf = new fluorescencia[Nf]; }
  int           GetNf (void)               { return Nf; }          // Retorna a quantidade de linhas alocadas
  
  // Preenche os dados de uma linha específica dentro do array de fluorescências
  void   XRF(int pos, const char* l, double p, double pa, double e) {
    xrf[pos].Linha(l);  // Define a nomenclatura da linha na posição 'pos'
    xrf[pos].Prob(p);   // Define a probabilidade individual na posição 'pos'
    xrf[pos].ProbA(pa); // Define a probabilidade acumulada na posição 'pos'
    xrf[pos].E(e);      // Define a energia em keV na posição 'pos'
  }
  
  // Métodos de acesso direto aos elementos do vetor dinâmico de fluorescência
  char   *Linha (int pos) { return xrf[pos].Linha(); } // Retorna o nome da linha do índice informado
  double Prob  (int pos)  { return xrf[pos].Prob(); }  // Retorna a probabilidade do índice informado
  double ProbA (int pos)  { return xrf[pos].ProbA(); } // Retorna a probabilidade acumulada do índice informado
  double E     (int pos)  { return xrf[pos].E(); }     // Retorna a energia do índice informado
  
  // Copia o vetor de 6 coeficientes passado por parâmetro para o atributo membro da classe
  void   Coef (double *c) {
    for (int i = 0; i < 6; i++) { // Loop executado 6 vezes para copiar cada termo do ajuste polinomial
      coef[i] = c[i]; 
    }
  }
                                            
  double Coef (int pos)  { return coef[pos]; } // Retorna um coeficiente específico pela sua posição no vetor
  
  // Calcula a seção de choque fotoelétrica parcial deste subnível para uma dada energia incident (hv)
  double SecaoChoqueAbsFotoeletrica(double hv) {
    double s = 0, p = 1, lnE = log(hv); // Inicializa a soma acumulada 's', a potência 'p' (hv^0 = 1) e calcula ln(hv)
    
    // Avalia o polinômio: s = c0 + c1*ln(E) + c2*(ln(E))^2 + ... + c5*(ln(E))^5
    for (int i = 0; i < 6; i++) {
      s += coef[i] * p; // Acumula o termo atual do polinômio na variável 's'
      p *= lnE;         // Eleva lnE à próxima potência para a seguinte iteração
    }
    
    return exp(s); // Retorna a seção de choque na escala linear aplicando a função exponencial
  }
};

//-------------------------------------------------Class atomo-----------------------------------

// Classe principal que representa o elemento atômico e suas propriedades físico-químicas
class atomo {
private: 
  char   simb[3];       // Símbolo químico do elemento (ex: "Ag", "B")
  int    Na;            // Número Atômico Z (quantidade de prótons)
  double Ma;            // Massa Atômica A (em u.m.a ou g/mol)
  double rho;           // Densidade do elemento puro (em g/cm³)
  double coer[6];       // Coeficientes do espalhamento coerente (Rayleigh)
  double incoer[6];     // Coeficientes do espalhamento incoerente (Compton)  
  int Ns;               // Número total de subníveis eletrônicos cadastrados
  subnivel *ene;        // Ponteiro para o array dinâmico contendo os subníveis
    
public:  
  // Encapsulamento / Getters e Setters do átomo
  void          Simb(const char *S) { strcpy(simb, S); } // Define o símbolo químico
  const char*   Simb(void)          { return simb; }     // Retorna o símbolo químico
  
  void          Z(int n)            { Na = n; }          // Define o número atômico Z
  int           Z(void)             { return Na; }       // Retorna o número atômico Z
  
  void          A(double m)         { Ma = m; }          // Define a massa atômica A
  double        A(void)             { return Ma; }       // Retorna a massa atômica A
  
  void          Rho(double r)       { rho = r; }         // Define a densidade
  double        Rho(void)           { return rho; }      // Retorna a densidade
  
  // Aloca dinamicamente o vetor de subníveis eletrônicos
  void          GetNs (int n)       { Ns = n; ene = new subnivel[Ns]; }
  int           GetNs (void)        { return Ns; }       // Retorna a quantidade de subníveis
  
  // Protótipo da função que lê e analisa o arquivo de dados
  void LerArquivo(const char *Nome_arq);
  
  // Exibe todos os parâmetros carregados do arquivo formatados no terminal
  void ImprimirDados(void) {
    // Imprime os dados gerais do elemento químico
    printf("Simb: %s\n", simb);
    printf("Z: %d\tA: %lf\trho: %lf\n", Na, Ma, rho);
    
    // Imprime o vetor de coeficientes do espalhamento coerente
    printf("Coer:\n");
    for(int i = 0; i < 6; i++) printf("|%lf ", coer[i]);
    
    // Imprime o vetor de coeficientes do espalhamento incoerente
    printf("\nIncoer:\n");
    for(int i = 0; i < 6; i++) printf("|%lf ", incoer[i]);
    printf("\n\n");

    // Itera sobre todos os subníveis alocados no átomo
    for (int i = 0; i < Ns; i++) {
      printf("Subnivel: %s\n", ene[i].Simb()); // Imprime a identificação do subnível
      printf("Eb: %lf | rs: %lf | w: %lf\n", ene[i].GetEb(), ene[i].Rs(), ene[i].W()); // Imprime as propriedades físicas do subnível
      
      // Imprime os coeficientes fotoelétricos do subnível atual
      printf("coef: ");
      for (int j = 0; j < 6; j++) printf("|%lf ", ene[i].Coef(j));
      printf("\n\n");

      // Se existirem linhas de fluorescência cadastradas para este subnível, imprime a tabela
      if (ene[i].GetNf() > 0) {
        printf("\nlinha\t  prob\t     proba      energia\n");
        printf("-----------------------------------------\n");
        for (int j = 0; j < ene[i].GetNf(); j++) {
          printf("%s   |%lf   |%lf   |%lf\n", 
                 ene[i].Linha(j), 
                 ene[i].Prob(j), 
                 ene[i].ProbA(j), 
                 ene[i].E(j));
          printf("-----------------------------------------\n");
        }
        printf("\n"); 
      }
    }
  }
  
  // --- Métodos de Cálculo de Seção de Choque ---
  
  // Calcula a seção de choque de espalhamento elástico (Rayleigh) para a energia hv
  double SecaoChoqueElastico(double hv) {
    double s = 0, p = 1, lnE = log(hv); // Inicializa acumulador, fator de potência e calcula o logaritmo natural de hv
    
    // Avalia o polinômio logarítmico para o espalhamento coerente
    for (int i = 0; i < 6; i++) {
      s += coer[i] * p; // Soma o produto do coeficiente pela potência atual de ln(hv)
      p *= lnE;         // Atualiza a potência para a próxima iteração
    }
    return exp(s); // Eleva 'e' à potência resultante do ajuste polinomial
  }
  
  // Calcula a seção de choque de espalhamento inelástico (Compton) para a energia hv
  double SecaoChoqueIneElastico(double hv) {
    double s = 0, p = 1, lnE = log(hv); // Inicializa variáveis para o ajuste da curva Compton
    
    // Avalia o polinômio logarítmico para o espalhamento incoerente
    for (int i = 0; i < 6; i++) {
      s += incoer[i] * p; // Soma o produto do coeficiente pela potência de ln(hv)            
      p *= lnE;          // Incrementa a ordem da potência
    }
    return exp(s); // Retorna o valor calculado aplicando a exponencial
  }
  
  // Calcula a seção de choque fotoelétrica total somando a contribuição dos subníveis fotoeletricamente ativos
  double SecaoChoqueAbsFotoeletrica(double hv) {
    double s = 0; // Inicializa a seção de choque total acumulada com zero
    
    for (int i = 0; i < Ns; i++) { 
      // Condição física: a absorção no subnível só ocorre se a energia do fóton incident for maior ou igual à energia de ligação
      if (ene[i].GetEb() <= hv) {
        s += ene[i].SecaoChoqueAbsFotoeletrica(hv); // Acumula a seção de choque do subnível
      }
    }
    return s; // Retorna o total acumulado para a absorção fotoelétrica
  }
};

// --- IMPLEMENTAÇÃO DO MÉTODO DE LEITURA DO ARQUIVO ---
// Lê a estrutura do arquivo de dados espectrais e carrega as tabelas na memória
void atomo::LerArquivo(const char *Nome_arq) {
  FILE *arq = fopen(Nome_arq, "rt"); // Abre o arquivo especificado no modo de leitura de texto ("rt")
  if (!arq) { // Verifica se houve falha na abertura do arquivo (ex: caminho incorreto)
    printf("Erro ao abrir arquivo %s\n", Nome_arq);
    return; // Interrompe a execução do método
  }

  // Declaração de variáveis auxiliares temporárias para leitura do arquivo
  char lixo[64], linha[10];
  int i, j, pos, Nf;
  double e, r, w;
  double prob, proba, energia;
  double coef[6];

  // --- Lendo Dados Principais do Átomo ---
  fscanf(arq, "%s\n", simb);                    // Lê o símbolo químico do elemento
  fgets(lixo, 64, arq);                         // Lê e descarta a linha de cabeçalho do arquivo
  fscanf(arq, "%d %lf %lf ", &Na, &Ma, &rho);   // Lê o número atômico Z, a massa A e a densidade rho
  fgets(lixo, 64, arq);                         // Lê e descarta linha de comentários
  
  // Lê os 6 coeficientes do ajuste polinomial para espalhamento coerente
  for(i = 0; i < 6; i++) fscanf(arq, "%lf ", &coer[i]);
  fgets(lixo, 64, arq);                         // Descarta a linha divisória
  
  // Lê os 6 coeficientes do ajuste polinomial para espalhamento incoerente
  for(i = 0; i < 6; i++) fscanf(arq, "%lf ", &incoer[i]);
  fgets(lixo, 64, arq);                         // Descarta a linha divisória

  // --- Primeira Passagem no Arquivo: Contagem do total de Subníveis ---
  pos = ftell(arq); // Armazena a posição atual do ponteiro de leitura no arquivo
  Ns = 0;           // Inicializa o contador de subníveis
  
  // Percorre o arquivo contando quantas linhas iniciam com '*' (marcador do início do bloco de subnível)
  while (fgets(lixo, 64, arq) != NULL) {
    if (lixo[0] == '*') Ns++; // Incrementa o total de subníveis encontrados
  }

  GetNs(Ns);                 // Executa a alocação dinâmica do vetor de subníveis 'ene' com tamanho Ns
  fseek(arq, pos, SEEK_SET); // Reposiciona o ponteiro de leitura do arquivo de volta ao início dos subníveis

  // --- Segunda Passagem no Arquivo: Leitura detalhada dos Subníveis ---
  for (i = 0; i < Ns; i++) {
    fscanf(arq, "%s ", lixo); // Lê a string identificadora do subnível
    ene[i].Simb(lixo);        // Armazena a identificação no objeto do subnível
    
    fgets(lixo, 64, arq);     // Descarta a linha de cabeçalho
    
    // Lê as propriedades do subnível: energia de ligação (Eb), raio (rs) e largura de linha (w)
    fscanf(arq, "%lf %lf %lf ", &e, &r, &w); 
    ene[i].GetEb(e), ene[i].Rs(r), ene[i].W(w); // Atribui os dados lidos
    
    fgets(lixo, 64, arq);     // Descarta linha de comentário intermediária
    
    // Lê os 6 coeficientes fotoelétricos associados ao subnível atual
    for(j = 0; j < 6; j++) {
      fscanf(arq, "%lf ", &coef[j]);
    }
    ene[i].Coef(coef); // Salva o array de coeficientes no objeto do subnível
    
    // --- Lendo os Dados de Fluorescência do Subnível Atual ---
    fgets(lixo, 64, arq);     // Descarta linha de separação
    fgets(lixo, 64, arq);     // Descarta o cabeçalho das colunas de fluorescência
    
    pos = ftell(arq);         // Salva a posição do ponteiro antes da contagem das linhas de fluorescência
    
    Nf = 0;                   // Reseta o contador de linhas de fluorescência do subnível
    
    // Conta o número de transições de fluorescência até alcançar o marcador '*' de fim do bloco
    while(fgets(lixo, 64, arq) != NULL) {
      if (lixo[0] == '*') break; // Encontrou o final das linhas de fluorescência deste subnível
      Nf++;                      // Incrementa o número de transições encontradas
    }

    ene[i].GetNf(Nf);          // Aloca dinamicamente a memória do vetor de fluorescência do subnível
       
    fseek(arq, pos, SEEK_SET); // Restaura o ponteiro para o início da lista de linhas de fluorescência

    // Lê os valores de cada transição de fluorescência do subnível
    for (j = 0; j < ene[i].GetNf(); j++) {
      fscanf(arq, "%s %lf %lf %lf ", 
             linha, 
             &prob, 
             &proba, 
             &energia);
      ene[i].XRF(j, linha, prob, proba, energia); // Insere os dados lidos no objeto de fluorescência
    }
    fgets(lixo, 64, arq);     // Lê e descarta a linha que contém o caractere '*'
  }

  fclose(arq);                 // Fecha o ponteiro do arquivo aberto
}


// --- FUNÇÃO PRINCIPAL ---
int main(void) {
  
  // Declaração de variáveis para varredura de energia e coeficientes de atenuação
  double hv, e = 20.0; // Define o limite máximo da energia incidente 'e' em keV (20.0 keV)
  double sigma_e, sigma_i, tau; // Variáveis para seções de choque: elástica (Rayleigh), inelástica (Compton) e fotoelétrica (tau)
  double mu; // Coeficiente de atenuação total (soma das seções de choque)
  
  atomo atm; // Instancia o objeto 'atm' referente à classe 'atomo'
  
  // Executa o carregamento dos dados do elemento Prata (Ag) a partir da base de dados física
  atm.LerArquivo("/home/administrador/Documentos/vitor_modelagem/dbase/Ag.txt");
  
  // Imprime os parâmetros estruturais do elemento lido
  atm.ImprimirDados();
  
  // Imprime o cabeçalho da tabela de resultados físicos calculados
  printf("    hv |  sigma_e|  sigma_i|    tau   |     mu\n");
  printf("-------|---------|---------|----------|--------------\n");
  
  // Laço de repetição: Varia a energia incidente (hv) de 1.0 keV até 20.0 keV em passos de 1.0 keV
  for (hv = 1.0; hv <= e; hv += 1.0) {
  
    sigma_e = atm.SecaoChoqueElastico(hv);       // Calcula a seção de choque elástica (Rayleigh) em cm²/g
    sigma_i = atm.SecaoChoqueIneElastico(hv);     // Calcula a seção de choque inelástica (Compton) em cm²/g
    tau     = atm.SecaoChoqueAbsFotoeletrica(hv);  // Calcula a seção de choque fotoelétrica total em cm²/g
    mu      = sigma_e + sigma_i + tau;           // Calcula o coeficiente de atenuação total pela soma dos três processos
    
    // Imprime no console os valores calculados formatados em colunas alinhadas
    printf("%7.3lf| %8.5lf| %8.5lf| %9.3lf| %9.3lf\n", hv, sigma_e, sigma_i, tau, mu);
    printf("-------|---------|---------|----------|--------------\n");

  }

  return 0; // Retorna 0 para o sistema operacional informando que o programa finalizou sem erros
}
