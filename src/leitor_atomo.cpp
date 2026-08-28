#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
  char   linha[6];      // Nome da transicao de Raio-X (ex: K-L3, KL3)
  double prob;          // Probabilidade de transicao (taxa/rate)
  double proba;         // Probabilidade acumulada (irate)
  double energia;       // Energia do foton emitido hv (em keV)
} fluorescencia;        // Estrutura de dados para linhas de fluorescencia

typedef struct {
  char      simb[6];     // Nome do subnivel (ex: K, L1, L2, L3)
  double    Eb;          // Energia de ligacao (Binding Energy, em keV)
  double    rs;          // Raio de espalhamento / fator de forma
  double    w;           // Largura de linha / rendimento de fluorescencia
  double    coef[6];     // Coeficientes de absorcao fotoeletrica (Array fixo de 6 posicoes)
  
  int Nf;                // Quantidade de linhas de fluorescencia deste subnivel
  fluorescencia *xrf;    // Ponteiro para o array dinamico de fluorescencias
} subnivel;              // Estrutura de dados do subnivel atomico

typedef struct {
  char   simb[3];        // Simbolo quimico (ex: "B", "C", "Fe")
  int    Z;              // Numero atomico (Z)
  double A;              // Massa atomica (A)
  double rho;            // Densidade (g/cm³)
  double coer[6];        // Coeficientes de espalhamento coerente (Rayleigh)
  double incoer[6];      // Coeficientes de espalhamento incoerente (Compton)
  
  int Ns;                // Numero total de subniveis do atomo
  subnivel *ene;         // Ponteiro para o array dinamico de subniveis
} atomo;                 // Estrutura de dados principal do atomo

int main(void) {
  FILE *arq;             // Ponteiro que guarda o endereco da estrutura FILE
  atomo atm;             // Variavel principal que armazena o atomo
  char lixo[64];         // Buffer temporario para descartar linhas/cabecalhos
  int i, j, pos;         // Variaveis de controle de loops e posicao no arquivo
  
  // Abertura do arquivo da base de dados (use caminho relativo na entrega final)
  arq = fopen("data/B.txt", "rt"); 
  if (arq == NULL) {
    printf("Erro ao abrir o arquivo!\n");
    return 1;
  }

  // --- LEITURA DOS DADOS DO ATOMO ---
  fscanf(arq, "%s\n", atm.simb);                              // Le o simbolo quimico
  fgets(lixo, 64, arq);                                       // Descarta cabecalho de Z, A, rho
  fscanf(arq, "%d %lf %lf ", &atm.Z, &atm.A, &atm.rho);       // Le Z, A e rho
  fgets(lixo, 64, arq);                                       // Descarta cabecalho de coerentes
  
  for(i = 0; i < 6; i++) {
    fscanf(arq, "%lf ", &atm.coer[i]);                        // Le os 6 coeficientes coerentes
  }
  
  fgets(lixo, 64, arq);                                       // Descarta cabecalho de incoerentes
  
  for(i = 0; i < 6; i++) {
    fscanf(arq, "%lf ", &atm.incoer[i]);                      // Le os 6 coeficientes incoerentes
  }
  
  fgets(lixo, 64, arq);                                       // Descarta linha separadora
  
  // --- CONTAGEM DE SUBNIVEIS (NS) ---
  pos = ftell(arq);                                           // Salva a posicao inicial dos subniveis
  atm.Ns = 0;                                                 // Inicializa o contador de subniveis
  
  while (fgets(lixo, 64, arq) != NULL) {                      // Le linha por linha ate o fim do arquivo
    if (lixo[0] == '*') {
      atm.Ns++;                                               // Se a linha comeca com '*', conta +1 subnivel
    }
  }

  // Alocacao dinamica da memoria para os subniveis (requer compilador C++)
  atm.ene = new subnivel[atm.Ns];
  
  // Retorna o ponteiro do arquivo para onde comecam os subniveis
  fseek(arq, pos, SEEK_SET);
  
  // --- LEITURA DOS SUBNIVEIS ---
  for (i = 0; i < atm.Ns; i++) {
    fscanf(arq, "%s ", atm.ene[i].simb);                      // Le o nome do subnivel (ex: K, L1)
    fgets(lixo, 64, arq);                                     // Descarta cabecalho de Eb, rs, w
    fscanf(arq, "%lf %lf %lf ", &atm.ene[i].Eb, &atm.ene[i].rs, &atm.ene[i].w); // Le Eb, rs, w
    fgets(lixo, 64, arq);                                     // Descarta cabecalho dos coeficientes
    
    for(j = 0; j < 6; j++) {
      fscanf(arq, "%lf ", &atm.ene[i].coef[j]);               // Le os 6 coeficientes de absorcao
    }
    
    // --- LEITURA DA FLUORESCENCIA DO SUBNIVEL ---
    fgets(lixo, 64, arq);                                     // Descarta linha '#'
    fgets(lixo, 64, arq);                                     // Descarta cabecalho '# xray rate irate hv'
    
    pos = ftell(arq);                                         // Salva posicao onde iniciam as linhas XRF
    atm.ene[i].Nf = 0;                                        // Inicializa contador de fluorescencias
    
    while(fgets(lixo, 64, arq) != NULL) {
      if (lixo[0] == '*') {                                   // O caractere '*' indica o fim deste subnivel
        break;
      }
      atm.ene[i].Nf++;                                        // Incrementa contador de linhas de fluorescencia
    }
    
    // Aloca a tabela dinamica de fluorescencia para o subnivel i
    atm.ene[i].xrf = new fluorescencia[atm.ene[i].Nf];
    
    // Retorna para o inicio dos dados de fluorescencia
    fseek(arq, pos, SEEK_SET);
    
    for (j = 0; j < atm.ene[i].Nf; j++) {
      fscanf(arq, "%s %lf %lf %lf ", 
             atm.ene[i].xrf[j].linha, 
             &atm.ene[i].xrf[j].prob, 
             &atm.ene[i].xrf[j].proba, 
             &atm.ene[i].xrf[j].energia);                     // Le dados de cada transicao XRF
    }
    
    fgets(lixo, 64, arq);                                     // Consome a linha contendo '*' no final do bloco
  }

  fclose(arq);                                                // Fecha o arquivo de dados

  // --- IMPRESSAO DOS RESULTADOS ---
  printf("Simb: %s\n", atm.simb);
  printf("Z: %d\tA: %lf\trho: %lf\n", atm.Z, atm.A, atm.rho);
  
  printf("Coer:\n");
  for(i = 0; i < 6; i++) printf("%lf ", atm.coer[i]);
  printf("\nIncoer:\n");
  for(i = 0; i < 6; i++) printf("%lf ", atm.incoer[i]);
  printf("\n\n");

  for (i = 0; i < atm.Ns; i++) {
    printf("Subnivel: %s\n", atm.ene[i].simb);
    printf("Eb: %lf | rs: %lf | w: %lf\n", atm.ene[i].Eb, atm.ene[i].rs, atm.ene[i].w);
    
    printf("Coef: ");
    for (j = 0; j < 6; j++) printf("%lf ", atm.ene[i].coef[j]);
    printf("\n");
    
    if (atm.ene[i].Nf > 0) {
      printf("\n  linha\t\tprob\t\tproba\t\tenergia\n");
      for (j = 0; j < atm.ene[i].Nf; j++) {
        printf("  %s\t\t%lf\t%lf\t%lf\n", 
               atm.ene[i].xrf[j].linha, 
               atm.ene[i].xrf[j].prob, 
               atm.ene[i].xrf[j].proba, 
               atm.ene[i].xrf[j].energia);
      }
      printf("\n");
    }
  }
