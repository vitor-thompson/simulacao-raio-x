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
  double prob;          // Probabilidade de transição
  double proba;         // Probabilidade acumulada
  double energia;       // Energia do fóton emitido (hv)

public: // Tornado público para permitir leitura e impressão direta
  // Métodos Getters e Setters para manipulação dos atributos privados da fluorescência
  void          Linha (const char* nome)   {strcpy (linha, nome);} // Define o nome da linha
  char*         Linha (void)               {return linha;}         // Retorna o nome da linha
  
  void          Prob  (double e)           { prob = e ; }          // Define a probabilidade
  double        Prob  (void)               { return prob; }        // Retorna a probabilidade
  
  void          ProbA (double e)           { proba = e ; }         // Define a probabilidade acumulada
  double        ProbA (void)               { return proba; }       // Retorna a probabilidade acumulada
  
  void          E     (double e)           { energia = e ; }       // Define a energia do fóton
  double        E     (void)               { return energia; }     // Retorna a energia do fóton
};                        

//----------------------------------------------------class subnivel-----------------------------------

// Classe que representa um subnível eletrônico do átomo (ex: K, L1, L2...)
class subnivel {
private:
  char      simb[6];    // Nome do subnivel
  double    Eb;         // Energia de ligacao (Binding Energy - Eb)
  double    rs;         // Raio de espalhamento / fator de forma
  double    w;          // Largura de linha / peso estatistico
  double    coef[6];    // Coeficientes para secao de choque fotoeletrica
  
  int Nf;               // Número de linhas de fluorescência associadas a este subnível
  fluorescencia *xrf;   // Ponteiro para o array dinâmico de linhas de fluorescência
  
public: // Tornado público para ser acessado pela classe atomo 
  // Getters e Setters do subnível
  void          Simb   (const char *nome)  { strcpy(simb, nome); } // Define o símbolo do subnível
  const char*   Simb   (void)              { return simb; }        // Retorna o símbolo do subnível
  
  void          GetEb  (double e)          { Eb = e ; }            // Define a energia de ligação
  double        GetEb  (void)              { return Eb; }          // Retorna a energia de ligação
  
  void          Rs     (double r)          { rs = r ; }            // Define o raio de espalhamento
  double        Rs     (void)              { return rs; }          // Retorna o raio de espalhamento
  
  void          W      (double e)          { w = e ; }             // Define a largura de linha
  double        W      (void)              { return w; }           // Retorna a largura de linha
  
  // Aloca dinamicamente o array de fluorescências para o tamanho Nf fornecido
  void          GetNf (int n)              { Nf = n ; xrf = new fluorescencia[Nf];}
  
  int           GetNf (void)               { return Nf; }          // Retorna a quantidade de linhas de fluorescência
  
  // Preenche os dados de um elemento específico do vetor de fluorescência
  void   XRF(int pos, const char* l, double p, double pa, double e) {
    xrf[pos].Linha(l);
    xrf[pos].Prob(p);
    xrf[pos].ProbA(pa);
    xrf[pos].E(e);
  }
  
  // Métodos para obter dados de uma determinada posição do vetor de fluorescência
  char   *Linha (int pos) {return xrf[pos].Linha(); } // Retorna o nome da linha na posição pos
  double Prob  (int pos)  {return xrf[pos].Prob(); }  // Retorna a probabilidade na posição pos
  double ProbA (int pos)  {return xrf[pos].ProbA(); } // Retorna a probabilidade acumulada na posição pos
  double E     (int pos)  {return xrf[pos].E(); }     // Retorna a energia na posição pos
  
  // Define o vetor de 6 coeficientes fotoelétricos do subnível
  void   Coef (double *c) {
    int i;
    for (i = 0; i < 6; i++)
      coef[i] = c[i]; 
  }
                                            
  double Coef (int pos)  { return coef[pos]; } // Retorna um coeficiente específico pela posição
  
  // Calcula a seção de choque fotoelétrica do subnível para a energia hv dada por um polinômio em ln(hv)
  double SecaoChoqueAbsFotoeletrica(double hv) {
    double s = 0, p = 1, lnE = log(hv); 
    
    for (int i = 0; i < 6; i++) {
      s += coef[i] * p;
      p *= lnE;          
    }
    
    return exp(s); // Retorna a seção de choque calculada em escala exponencial
  }
  
};

//-------------------------------------------------Class atomo-----------------------------------

// Classe principal que representa o átomo e suas propriedades físicas e de interação
class atomo {
private: 
  // Atributos privados
  char   simb[3];       // Simbolo quimico
  int    Na;            // Numero Atomico (Z)
  double Ma;            // Massa Atomica (A)
  double rho;           // Densidade
  double coer[6];       // Coeficientes Coerente
  double incoer[6];     // Coeficientes Incoerente  
  int Ns;               // Numero de subniveis
  subnivel *ene;        // Ponteiro para o array dinamico de subniveis
    
public:  
  // Encapsulamento / Getters e Setters
  void          Simb(const char *S) { strcpy(simb, S); } // Define o símbolo químico
  const char*   Simb(void)          { return simb; }     // Retorna o símbolo químico
  
  void          Z(int n)            { Na = n; }          // Define o número atômico Z
  int           Z(void)             { return Na; }       // Retorna o número atômico Z
  
  void          A(double m)         { Ma = m; }          // Define a massa atômica A
  double        A(void)             { return Ma; }       // Retorna a massa atômica A
  
  void          Rho(double r)       { rho = r; }         // Define a densidade rho
  double        Rho(void)           { return rho; }      // Retorna a densidade rho
  
  // Aloca dinamicamente o array de subníveis para o tamanho Ns fornecido
  void          GetNs (int n)       { Ns = n ; ene = new subnivel[Ns];}
  int           GetNs (void)        { return Ns; }       // Retorna o número de subníveis
  
  // Declaração do método de leitura de arquivo estruturado
  void LerArquivo(const char *Nome_arq);
  
  // Exibe todos os dados carregados na memória de forma formatada no console
  void ImprimirDados(void) {
  
    int i, j;
   
    // Imprime os dados gerais do átomo
    printf("Simb: %s\n", simb);
    printf("Z: %d\tA: %lf\trho: %lf\n", Na, Ma, rho);
    
    // Imprime os coeficientes de espalhamento
    printf("Coer:\n");
    for(int i = 0; i < 6; i++) printf("%lf ", coer[i]);
    printf("\nIncoer:\n");
    for(int i = 0; i < 6; i++) printf("%lf ", incoer[i]);
    printf("\n\n");

    // Loop para imprimir cada um dos subníveis cadastrados
    for (int i = 0; i < Ns; i++) {
      printf("Subnivel: %s\n", ene[i].Simb());
      printf("Eb: %lf | rs: %lf | w: %lf\n", ene[i].GetEb(), ene[i].Rs(), ene[i].W());
      
      printf("coef: ");
      for (int j = 0; j < 6; j++) printf("%lf ", ene[i].Coef(j));
      printf("\n");

      // Se o subnível possuir linhas de fluorescência, imprime a tabela de transição
      if (ene[i].GetNf() > 0) {
        printf("\nlinha\tprob\tproba\tenergia\n");
        for (int j = 0; j < ene[i].GetNf(); j++) {
          printf("%s %lf %lf %lf\n", 
                 ene[i].Linha(j), 
                 ene[i].Prob(j), 
                 ene[i].ProbA(j), 
                 ene[i].E(j));
        }
        printf("\n"); 
      }
    }
  }
  
  // Metodos de Calculo de Secao de Choque
  
  // Calcula a seção de choque elástica (Coerente - Rayleigh)
  double SecaoChoqueElastico(double hv) {
    double s = 0, p = 1, lnE = log(hv);
    for (int i = 0; i < 6; i++) {
      s += coer[i] * p;
      p *= lnE;          
    }
    return exp(s);
  }
  
  // Calcula a seção de choque inelástica (Incoerente - Compton)
  double SecaoChoqueIneElastico(double hv) {
    double s = 0, p = 1, lnE = log(hv);
    for (int i = 0; i < 6; i++) {
      s += incoer[i] * p;             
      p *= lnE;          
    }
    return exp(s);
  }
  
  // Calcula a seção de choque fotoelétrica total somando as contribuições de subníveis válidos
  double SecaoChoqueAbsFotoeletrica(double hv) {
    double s = 0;
    for (int i = 0; i < Ns; i++) { 
      // Apenas subníveis com energia de ligação menor ou igual à energia do fóton incidentes participam
      if (ene[i].GetEb() <= hv) {
        s += ene[i].SecaoChoqueAbsFotoeletrica(hv);
      }
    }
    return s;
  }
};

// --- IMPLEMENTACAO DO METODO LER ARQUIVO ---
// Método responsável por realizar a leitura do arquivo texto e popular os objetos
void atomo::LerArquivo(const char *Nome_arq) {
  FILE *arq = fopen(Nome_arq, "rt"); // Abre o arquivo em modo texto para leitura
  if (!arq) {
    printf("Erro ao abrir arquivo %s\n", Nome_arq);
    return;
  }

  // Variáveis auxiliares para armazenamento temporário dos dados do arquivo
  char lixo[64], linha[10];
  int i, j, pos, Nf;
  double e, r, w;
  double prob, proba, energia;
  double coef[6];

  // --- Lendo Dados Principais do Atomo ---
  fscanf(arq, "%s\n", simb);                    // Lê o símbolo do átomo
  fgets(lixo, 64, arq);                         // Descarta a linha de cabeçalho
  fscanf(arq, "%d %lf %lf ", &Na, &Ma, &rho);   // Lê Z, A e a densidade rho
  fgets(lixo, 64, arq);                         // Descarta linha de comentário
  
  // Lê os 6 coeficientes do espalhamento coerente
  for(i = 0; i < 6; i++) fscanf(arq, "%lf ", &coer[i]);
  fgets(lixo, 64, arq);                         // Descarta linha de comentário
  
  // Lê os 6 coeficientes do espalhamento incoerente
  for(i = 0; i < 6; i++) fscanf(arq, "%lf ", &incoer[i]);
  fgets(lixo, 64, arq);                         // Descarta linha de comentário

  // --- Primeira Passagem: Contagem de Subníveis ---
  pos = ftell(arq); // Salva a posição atual do ponteiro do arquivo
  Ns = 0;
  
  // Conta a quantidade de marcadores '*' no arquivo para determinar o total de subníveis
  while (fgets(lixo, 64, arq) != NULL) {
    if (lixo[0] == '*') Ns++;
  }

  GetNs(Ns);               // Aloca o vetor dinâmico de subníveis com tamanho Ns
  fseek(arq, pos, SEEK_SET); // Restaura o ponteiro do arquivo para a posição salva

  // --- Segunda Passagem: Leitura Individual de cada Subnível ---
  for (i = 0; i < Ns; i++) {
    fscanf(arq, "%s ", lixo); // Lê o nome do subnível
    ene[i].Simb(lixo);
    
    fgets(lixo, 64, arq);     // Descarta o cabeçalho
    
    // Lê as constantes físicas do subnível (Eb, rs, w)
    fscanf(arq, "%lf %lf %lf ", &e, &r, &w); 
    ene[i].GetEb(e), ene[i].Rs(r), ene[i].W(w);
    
    fgets(lixo, 64, arq);     // Descarta a linha de comentário
    
    // Lê os 6 coeficientes fotoelétricos do subnível
    for(j = 0; j < 6; j++) {
      fscanf(arq, "%lf ", &coef[j]);
    }
    ene[i].Coef(coef);
    
    // --- Lendo Dados da Fluorescência do Subnível Atual ---
    fgets(lixo, 64, arq);     // Descarta linha divisória
    fgets(lixo, 64, arq);     // Descarta cabeçalho
    
    pos = ftell(arq);         // Salva a posição do arquivo antes de contar as linhas de fluorescência
    
    Nf = 0;
    
    // Conta quantas linhas de fluorescência existem até encontrar o '*' de encerramento
    while(fgets(lixo, 64, arq) != NULL) {
      if (lixo[0] == '*') break;
      Nf++;
    }

    ene[i].GetNf(Nf);          // Aloca o vetor de fluorescências para este subnível
       
    fseek(arq, pos, SEEK_SET); // Retorna o ponteiro para o início do bloco de fluorescência

    // Lê os valores de cada linha de fluorescência do subnível
    for (j = 0; j < ene[i].GetNf(); j++) {
      fscanf(arq, "%s %lf %lf %lf ", 
             linha, 
             &prob, 
             &proba, 
             &energia);
      ene[i].XRF(j, linha, prob, proba, energia); // Armazena no objeto de fluorescência
    }
    fgets(lixo, 64, arq);     // Descarta a linha do '*' de fechamento
  }

  fclose(arq);                 // Fecha o arquivo texto aberto
}


// --- FUNCAO PRINCIPAL ---
int main(void) {
  atomo atm; // Instancia o objeto principal da classe atomo
  
  // Chama a leitura do arquivo apontando para o caminho correto
  atm.LerArquivo("/home/administrador/Documentos/vitor_modelagem/dbase/B.txt");
  
  // Exibe os dados utilizando o metodo da classe
  atm.ImprimirDados();

  return 0; // Finaliza a execução do programa com sucesso
}
