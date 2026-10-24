#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// ==========================================
// CLASSE: fluorescencia
// Armazena os dados de transição de raios X (XRF) de um subnível
// ==========================================
class fluorescencia {
private:
  char linha[6];          // Nome da linha de emissão (ex: K-L3)
  double prob, proba, energia; // Probabilidade, probabilidade acumulada e energia (keV)
  
public:
  // Métodos de acesso (Getters e Setters)
  void  Linha(const char* nome) { strcpy(linha,nome); }
  char* Linha(void)             { return linha; }
  
  void   Prob(double e) { prob = e; }
  double Prob(void)     { return prob; }

  void   Proba(double e) { proba = e; }
  double Proba(void)     { return proba; }

  void   E(double e) { energia = e; }
  double E(void)     { return energia; }
};

// ==========================================
// CLASSE: subnivel
// Representa um subnível eletrônico do átomo (ex: K, L1, L2...)
// ==========================================
class subnivel {
private:
  char simb[6];           // Nome/símbolo do subnível
  double eb, rs, w;       // Eb: energia de ligação, rs: razão de salto, w: rendimento/largura
  double coef[6];         // Coeficientes do polinômio para a seção de choque fotoelétrica
  int nf;                 // Número de linhas de fluorescência associadas
  fluorescencia *xrf;     // Ponteiro para o vetor dinâmico de fluorescências

public:
  // Métodos Getters e Setters
  void  Simb(const char *nome) { strcpy(simb,nome); }
  char* Simb(void)             { return simb; }

  void   Eb(double e) { eb = e; }
  double Eb(void)     { return eb; }

  void   RS(double r) { rs = r; }
  double RS(void)     { return rs; }

  void   W(double e) { w = e; }
  double W(void)     { return w; }
  
  // Aloca dinamicamente o array de fluorescências com tamanho 'n'
  void   NF(int n)   { nf = n; xrf = new fluorescencia[nf]; }
  int    NF(void)    { return nf; }
  
  // Preenche os dados de uma fluorescência específica na posição 'pos'
  void   XRF(int pos, const char* l, double p, double pa, double e ) {
    xrf[pos].Linha(l); xrf[pos].Prob(p); xrf[pos].Proba(pa); xrf[pos].E(e);  
  }
  
  char*  Linha(int pos)  { return xrf[pos].Linha(); }
  double Prob(int pos)   { return xrf[pos].Prob(); }
  double Proba(int pos)  { return xrf[pos].Proba(); }
  double E(int pos)      { return xrf[pos].E(); }
  
  // Copia os 6 coeficientes fotoelétricos
  void   Coef(double *c) {
    int i; 
    for ( i=0; i<6; i++ ) coef[i] = c[i];
  }
  double Coef(int pos) { return coef[pos]; }

  // Calcula a seção de choque fotoelétrica parcial deste subnível usando expansão logarítmica
  double SecaoChoqueAbsFotoeletrica(double hv) {
    double s = 0, p = 1, lnE = log(hv);
    for (int i=0; i<6; i++) {
      s += coef[i]*p;
      p *= lnE;
    }
    return exp(s);
  }
};

// ==========================================
// CLASSE: atomo
// Representa o elemento químico alvo e suas propriedades de interação com radiação
// ==========================================
class atomo {
private:
  char simb[3];           // Símbolo químico do átomo (ex: Ag, C)
  int  na;                // Número atômico (Z)
  double ma, rho;         // Massa atômica (A) e densidade (g/cm³)
  double coer[6], incoer[6]; // Coeficientes de espalhamento elástico (Rayleigh) e inelástico (Compton)
  subnivel *ene;          // Ponteiro para o vetor dinâmico de subníveis
  int ns;                 // Número total de subníveis

public:
  void   Simb(const char *s) { strcpy(simb,s); }
  char  *Simb(void)          { return simb; }
  
  void          Z(unsigned char n) { na = n; }
  unsigned char Z(void)            { return na; }
  
  void   A(double m) { ma = m; }
  double A(void)     { return ma; }
  
  void   Rho(double d) { rho = d; }
  double Rho(void)     { return rho; }
  
  // Aloca o vetor de subníveis com tamanho 'n'
  void   NS(int n)   { ns = n; ene = new subnivel[ns]; }
  int    NS(void)    { return ns; }
  
  // Métodos de atalho para acessar propriedades internas dos subníveis
  double Eb(int i)       { return ene[i].Eb(); }
  double W(int i)        { return ene[i].W(); }
  double NF(int i)       { return ene[i].NF(); }
  double E(int i, int j)    { return ene[i].E(j); }
  double Prob(int i, int j) { return ene[i].Prob(j); }

  void   LerArquivo(const char *nome);
  
  // Calcula seção de choque de espalhamento elástico (Rayleigh)
  double SecaoChoqueElastico(double hv) {
    double s = 0, p = 1, lnE = log(hv);
    for (int i=0; i<6; i++) {
      s += coer[i]*p;
      p *= lnE;
    }
    return exp(s);
  }
  
  // Calcula seção de choque de espalhamento inelástico (Compton)
  double SecaoChoqueInelastico(double hv) {
    double s = 0, p = 1, lnE = log(hv);
    for (int i=0; i<6; i++) {
      s += incoer[i]*p;
      p *= lnE;
    }
    return exp(s);
  }
  
  // Calcula a seção de choque fotoelétrica total somando os subníveis válidos
  double SecaoChoqueAbsFotoeletrica(double hv)  {
    double f, s = 0;
    for ( int i=0; i<ns; i++ ) {
      if (ene[i].Eb() <= hv) { // Só absorve se a energia do fóton for maior que a energia de ligação
        f = ene[i].SecaoChoqueAbsFotoeletrica(hv);
        s += f;
      }
    }
    return s;
  }
  
  // Imprime no console todas as propriedades estruturais do átomo
  void Imprimir(void) {
    int i, j;
    printf("Simb: %s\n",simb);
    printf("Z: %d\tA: %lf\trho: %lf\n", na, ma, rho);
    for ( i=0; i<6; i++ ) printf("%lf ", coer[i]);
    printf("\n");
    for ( i=0; i<6; i++ ) printf("%lf ", incoer[i]);
    printf("\n");

    for ( i=0; i<ns; i++ ) {
      printf("%s\n",ene[i].Simb());
      printf("%lf\t%lf\t%lf\n", ene[i].Eb(), ene[i].RS(), ene[i].W());
      for ( j=0; j<6; j++ ) printf("%lf ",ene[i].Coef(j));
      printf("\n");
      if (ene[i].NF()>0) {
        printf("Linha\tProb\tProbA\tEnergia\n");
        for ( j=0; j<ene[i].NF(); j++ )
          printf("%s %lf %lf %lf\n", ene[i].Linha(j), ene[i].Prob(j),
                                     ene[i].Proba(j), ene[i].E(j));
      }
    }
  }
};

// Lê a base de dados em arquivo texto para preencher o objeto atomo
void atomo::LerArquivo(const char *nome) {
  FILE *fp;
  char  lixo[64], linha[10];
  int   i, j, pos, nf;
  double e, r, w;
  double prob, proba, energia;
  double coef[6];
  
  fp = fopen(nome,"rt");

  // Lê o símbolo químico do elemento
  fscanf(fp,"%s ", simb);
  fgets(lixo,64,fp); // Descarta cabeçalho
  // Lê número atômico, massa e densidade
  fscanf(fp,"%d %lf %lf ", &na, &ma, &rho);
  fgets(lixo,64,fp); // Descarta cabeçalho
  for ( i=0; i<6; i++ ) fscanf(fp,"%lf ", &coer[i]);
  fgets(lixo,64,fp); // Descarta cabeçalho
  for ( i=0; i<6; i++ ) fscanf(fp,"%lf ", &incoer[i]);
  fgets(lixo,64,fp); // Descarta cabeçalho
  
  pos = ftell(fp); // Salva a posição atual do ponteiro do arquivo

  // Primeira passada: conta quantos subníveis existem (marcador '*')
  ns = 0;
  while (1) {
    fgets(lixo,64,fp);
    if (feof(fp)!=0) break; 
    if (lixo[0]=='*') ns++;
  }

  NS(ns); // Aloca dinamicamente o espaço para os subníveis
  fseek(fp,pos,SEEK_SET); // Retorna o cursor para o início dos subníveis

  // Segunda passada: lê detalhadamente cada subnível e suas fluorescências
  for ( i=0; i<ns; i++ ) {
    fscanf(fp,"%s ",lixo);
    ene[i].Simb(lixo);
    
    fgets(lixo,64,fp);
    
    // Lê energia de ligação, razão de salto e rendimento de fluorescência
    fscanf(fp,"%lf %lf %lf ", &e, &r, &w);
    ene[i].Eb(e); ene[i].RS(r); ene[i].W(w);
    
    fgets(lixo,64,fp);
    for ( j=0; j<6; j++ ) fscanf(fp,"%lf ", &coef[j]);
    ene[i].Coef(coef);
    
    fgets(lixo,64,fp); // Descarta divisor #
    fgets(lixo,64,fp); // Descarta cabeçalho de raios X

    pos = ftell(fp); // Salva posição das linhas de fluorescência
    
    nf = 0;
    // Conta quantas linhas de fluorescência o subnível possui
    while (1) {
      fgets(lixo,64,fp);
      if (lixo[0]=='*') break;
      nf++;
    }
    
    ene[i].NF(nf); // Aloca o vetor de fluorescência do subnível
    fseek(fp,pos,SEEK_SET); // Posiciona cursor antes da tabela de fluorescência
    
    for ( j=0; j<ene[i].NF(); j++ ) {
      fscanf(fp,"%s %lf %lf %lf ", linha, &prob, &proba, &energia);
      ene[i].XRF(j, linha, prob, proba, energia);
    }
    fgets(lixo,64,fp); // Descarta o asterisco '*' de fechamento do bloco
  }
  
  fclose(fp);
}

// ==========================================
// CLASSE: tubo_rx
// Simula a geração do espectro contínuo e discreto de raios X de um tubo com anodo metálico
// ==========================================
class tubo_rx {
private:
   double E0; // Energia máxima do tubo (kV)
   double dE; // Passo de energia por canal
   int    Z;  // Número atômico do anodo
   double A;  // Massa atômica do anodo
   
   double m, J, eta, logZ, rho_zm, x, f_corr;
   double omega, t, i, phi, epsilon;

   double *spe, *speacum; // Espectro de emissão e espectro acumulado para Monte Carlo
   int    nc;             // Número total de canais do espectro
   atomo  atm;            // Objeto átomo correspondente ao material do anodo
   
   // Calcula a componente de radiação contínua (Bremsstrahlung)
   void Continuo(int n) {
      int i;
      double E;
      double cte = 1.35e9; 
      double U0, logU0, rho_z, tau_E_j, fa;

      nc = n;
      dE = E0/nc;
      spe = new double[nc];
      
      for ( i=0; i<nc; i++ ) {
         E = (i+1)*dE;
         U0 = E0/E;
         logU0 = log(U0);

         rho_z = rho_zm*logU0*(0.49269 - 1.09870*eta + 0.78557*eta*eta)/
                              (0.70256 - 1.09865*eta + 1.00460*eta*eta + logU0);

         tau_E_j = atm.SecaoChoqueAbsFotoeletrica(E);
         fa = U0>1 ? (1 - exp(-tau_E_j*2*rho_z*sin(phi)/sin(epsilon)))/
                              (tau_E_j*2*rho_z*sin(phi)/sin(epsilon)) : 1.0;

         spe[i] = t*cte*Z*pow(E0/E-1,x)*fa;
      }
      spe[0] = 0;
   }
   
   // Calcula a componente de raios X característicos (espectro discreto / linhas de fluorescência do anodo)
   void Discreto(void) {
      int k, l;
      double w_jk, p_jkl, zk, bk;
      double E_jk, E_jkl, tau_jkl;
      double U0_jk, logU0, rho_z, fa_jkl, S_jk_1, C_jkl, R;
      int    ch;
      
      for ( k=0; k <= (atm.NS()<=8 ? atm.NS() : 8); k++ ) {
         E_jk = atm.Eb(k); // Energia de ligação do subnível k
         if (E_jk < E0) {
            w_jk = atm.W(k);  // Rendimento de fluorescência

            U0_jk = E0/E_jk;
            logU0 = log(U0_jk);
            rho_z = rho_zm*logU0*(0.49269 - 1.09870*eta + 0.78557*eta*eta)/
                                 (0.70256 - 1.09865*eta + 1.00460*eta*eta + logU0);
            R = 1 - 0.0081517*Z + 3.613e-5*Z*Z + 0.009583*Z*exp(-U0_jk) + 0.001141*E0;

            // Constantes empíricas por subnível para o cálculo de intensidade
            switch (k){
               case 0: C_jkl = 1e13*(4.697 + 0.134*U0_jk - 0.00268*U0_jk*U0_jk)/
                                    (1 - exp(-7*(U0_jk - 1))); 
                       zk = 2; bk = 0.35;
                       break;
               case 1: C_jkl = f_corr*0.71e13; zk = 8; bk = 0.25; break;
               case 2: C_jkl = f_corr*2.7e13; zk = 8; bk = 0.25; break;
               case 3: C_jkl = 4.94e13; zk = 8; bk = 0.25; break;
               case 4: 
               case 5: C_jkl = 0; break;
               case 6: C_jkl = 2.32e13; zk = 1.2; bk = 0.7; break;
               case 7: C_jkl = 15.8e13; zk = 1.2; bk = 0.7; break;
               case 8: C_jkl = 20.5e13; zk = 1.2; bk = 0.7; break;
            }
            
            // Itera pelas linhas de fluorescência do subnível
            for ( l=0; l<atm.NF(k); l++ ) {
               E_jkl   = atm.E(k,l);     // Energia da linha XRF
               p_jkl   = atm.Prob(k,l);  // Probabilidade da transição

               tau_jkl = atm.SecaoChoqueAbsFotoeletrica(E_jkl);
               fa_jkl  = (1 - exp(-tau_jkl*2*rho_z*sin(phi)/sin(epsilon)))/
                                 (tau_jkl*2*rho_z*sin(phi)/sin(epsilon));

               S_jk_1 = (zk*bk/Z)*
                        (1.0 + U0_jk*(logU0 - 1.0))*
                        (1 + 16.05*sqrt(J/E_jk)*(sqrt(U0_jk)*logU0 + 2*(1-sqrt(U0_jk)))/
                                                     (1.0 + U0_jk*(logU0 - 1.0)));

               ch = int(E_jkl/dE); // Canal correspondente no espectro
               double e = t*C_jkl*R*S_jk_1*w_jk*p_jkl*fa_jkl;
               spe[ch] += e; // Acumula intensidade no canal do espectro discreto
            }
         }
      }
   }

public:
   // Configura o tubo de raios X, carrega o anodo, gera os espectros contínuo e discreto e cria a CDF acumulada
   void Anodo(double e, int n, const char *str) { 
      int idx;
      double sum, isum;
      E0 = e;
      atm.LerArquivo(str); 

      Z = atm.Z();
      x = 1.109 - 0.00435*Z + 0.00175*E0;
      A = atm.A();
      m = 0.1382 - 0.9211/sqrt(Z);
      J = 0.0135*Z;
      logZ = log(Z);
      eta = pow(E0,m)*
            (0.1904 - 0.2236*logZ + 0.1292*logZ*logZ - 0.0149*logZ*logZ*logZ);
      rho_zm = (A/Z)*
               (0.787e-5*sqrt(J)*pow(E0,1.5) + 0.735e-6*E0*E0);
      f_corr = (Z<=80) ? -0.4814 + 0.03781*Z - 2.413e-4*Z*Z : 1;

      omega = 0.00367; // Ângulo sólido (sr)
      t = 1e-6;        // Tempo de exposição (s)
      i = 1;           // Corrente (mA)

      phi = M_PI_2;          // Ângulo de decolagem (rad)
      epsilon = 35*M_PI/180; // Ângulo de incidência (rad)

      Continuo(n);
      Discreto();

      // Normaliza o espectro para transformar em probabilidades de emissão
      for (sum = 0, idx=0; idx<nc; idx++) sum += spe[idx];
      isum = 1/sum;
      for (idx=0; idx<nc; idx++) spe[idx] *= isum;
      
      // Cria o vetor de distribuição acumulada (CDF) para uso no método de Monte Carlo
      speacum = (double*)new double[nc];
      speacum[0] = spe[0];
      for (idx=1; idx<nc; idx++) 
         speacum[idx] = speacum[idx-1]+spe[idx];
    }
   
   void Imprimir(void) {
      for ( int i=0; i<nc; i++ )
         printf("%5.3lf\t%5.3e\n", (i+1)*dE, spe[i]);
   }
   
   void Imprimir(const char *nome) {
      FILE *fp;
      fp = fopen(nome,"wt");
      for ( int i=0; i<nc; i++ )
         fprintf(fp,"%5.3lf\t%5.3e\n", (i+1)*dE, spe[i]);
      fclose(fp);
   }
   
   // Sorteia uma energia de fóton via Monte Carlo baseada no espectro acumulado do tubo
   double Emissao(void) {
      int idx;
      double eps = rand()/(double)RAND_MAX; // Número randômico entre 0 e 1
      for (idx = 0; idx < nc; idx++ ) 
         if (speacum[idx]>eps) break;
      return (idx/(double)nc)*E0; // Retorna a energia correspondente sorteada
   }
};

// ==========================================
// FUNÇÃO PRINCIPAL (main)
// Executa a simulação completa: configura o tubo, carrega o alvo e simula interações estocásticas
// ==========================================
int main(int argc, char *argv[]) {

   char dir[128], arq[128];
   int     i;
   double  hv = 5.0; // Energia inicial padrão (keV)
   double  sigma_e, sigma_i, tau, mu;
   double  x = 0.01; // Espessura do material alvo (cm)
   double  Ptrans, Pee, Pei, Pabs;
   double  eps;
   
   srand(time(NULL)); // Inicializa a semente geradora de números aleatórios

   atomo   atm;
   tubo_rx f;

   // Validação dos argumentos passados via linha de comando
   if (argc != 4) {
      printf("1) Passar o caminho completo para o diretório dbase.\n");
      printf("2) Passar o elemento químico do anodo.\n");
      printf("3) Passar o elemento químico alvo.\n");
      printf("Exemplo:\n");
      printf("\t<programa> \"/home/user/work/dbase\" Ag C\n");
      exit(0);
   }

   // Monta o caminho do diretório da base de dados
   strcpy(dir,argv[1]);
   strcat(dir,"/");

   // Configura o arquivo de dados do anodo do tubo de raios X
   strcpy(arq,dir);
   strcat(arq,argv[2]);
   strcat(arq,".txt");
   
   printf("Configurando tubo de raio-x com anodo de %s.\n",argv[2]);
   
   f.Anodo(25, 1000, arq); // Inicializa tubo com 25 kV e 1000 canais
   f.Imprimir("espectro.txt"); // Salva o espectro gerado em arquivo
   
   // Configura o arquivo de dados do material alvo
   strcpy(arq,dir);
   strcat(arq,argv[3]);
   strcat(arq,".txt");
   
   atm.LerArquivo(arq);
   atm.Imprimir();

   printf("Carregando informação sobre material: %s.\n",argv[3]);

   // Loop de simulação de Monte Carlo para 10 fótons gerados pelo tubo
   for ( i=0; i<10; i++ ) {
   
      // Sorteia a energia do fóton emitida pelo tubo de raios X
      hv = f.Emissao();
      printf("hv = %lf\n",hv);
      
      // Cálculo das seções de choque para o fóton com energia 'hv'
      sigma_e = atm.SecaoChoqueElastico(hv);     // Rayleigh
      sigma_i = atm.SecaoChoqueInelastico(hv);   // Compton
      tau     = atm.SecaoChoqueAbsFotoeletrica(hv); // Fotoelétrica

      mu = sigma_e + sigma_i + tau; // Coeficiente de atenuação total

      printf("\tSigma_e: %lf\tSigma_i: %lf\tTau: %lf\n",sigma_e, sigma_i, tau);

      // Cálculo das probabilidades de interação física na espessura 'x'
      Ptrans = exp(-mu*atm.Rho()*x); // Probabilidade de Transmissão (Beer-Lambert)
      Pee    = sigma_e*(1-Ptrans)/mu; // Probabilidade de Espalhamento Elástico
      Pei    = sigma_i*(1-Ptrans)/mu; // Probabilidade de Espalhamento Inelástico
      Pabs   = tau*(-Ptrans + 1)/mu;  // Probabilidade de Absorção Fotoelétrica (ajustado ao original)

      printf("\tPtrans = %e\n",Ptrans);
      printf("\tPee    = %e\n",Pee);
      printf("\tPei    = %e\n",Pei);
      printf("\tPabs   = %e\n",Pabs);
      
      // Sorteio estatístico (Monte Carlo) para determinar o destino do fóton na matéria
      eps = rand()/(double)RAND_MAX;
      if (eps < Ptrans) {
         printf("\t\tTransmissão\n");
      } else if (eps < Ptrans + Pee) {
         printf("\t\tEspalhamento Elástico\n");
      } else if (eps < Ptrans + Pee + Pei) {
         printf("\t\tEspalhamento Inelástico\n");
      } else {
         printf("\t\tAbsorção Fotoelétrica\n");
      }      
   }
   return 0;
}
