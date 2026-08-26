#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {

  char      simb[6];                                                // Nome do subnivel;
  double    Eb;                                                     // Energia de ligacao (Binding Energy - Eb, em keV ou eV);                                                                         
  double    rs;                                                     // Raio de espalhamento / fator de forma ou parametro de escala (rs);
  double    w;                                                      // Largura de linha / peso estatistico de transicao (w);
  double    coef[6];                                                // Vetor dinâmico para coeficientes adicionais de atenuacao/absorcao;

} subnivel;                                                         // Estrutura de Dados;                                                                    
                                                                                                                            
typedef struct {                                                                                                              
                                                                                                                              
  char   simb[3];                                                   // Simbolo quimico;  
  int    Z;                                                         // Z   (Numero Atomico);   
  double A;                                                         // A   (Massa Atomica);
  double rho;                                                       // rho (Densidade);
  double coer[6];                                                   // Coeficientes Coerente[6];
  double incoer[6];                                                 // Coeficientes Incoerente[6];
  
  subnivel *ene;                                                    // Ponteiro;
  
  int Ns;                                                           // Numero de subniveis;                                                                 

} atomo;                                                            // Estrutura de Dados;
                                                                                                                                                                                   
                                                                                        
int main(void){                                                                       
                                                                              
  FILE *arq;                                                        // Ponteiro;
  atomo atm;                                                        // variavel do tipo atomo;
  char lixo[64];                                                    // descarte de string
  int i, pos;                                                       // variaveis;                                                              
  
  arq = fopen("/home/administrador/Documentos/vitor_modelagem/dbase/U.txt","rt"); // abre arquivo;
                                                                      
//----------------------------------------------------------------------------------------------------// 
  fscanf(arq,"%s\n", atm.simb);                                     // leitura do simbolo quimico;
                                                                              
  fgets(lixo, 64, arq);                                             // descarte da linha 2;
                                                                            
  fscanf(arq, "%d %lf %lf ", &atm.Z, &atm.A, &atm.rho);             // leitura de Z, A, rho;
                                                                            
  fgets(lixo, 64, arq);                                             // descarte a linha 4
                                                                        
  for(i = 0; i < 6; i++){                                             
    fscanf(arq,"%lf ", &atm.coer[i]);                                         
  }                                                                 // leitura dos 6 coeficientes coerentes; 
                                                                            
  fgets(lixo, 64, arq);                                             // descarte da linha 6;  
                                                                    
  for(i = 0; i < 6; i++){                                           
    fscanf(arq,"%lf ", &atm.incoer[i]);                               
  }                                                                 // leitura dos 6 coeficientes incoerentes;
  
  fgets(lixo, 64, arq);                                             // descarte da linha 8;  
  
  pos = ftell(arq);                                                 // Marca a possição do marcador;
  printf("pos= %d\n", pos);
  
  atm.Ns = 0;                                                       // Numero iniciais de subniveis
  while (1) {                           
    
    fgets(lixo, 64, arq);                                           // Descarte de estrings;
    
    if(feof(arq)!=0){
    
      break;
    
    }                                                               // Verifica se chegou ao final do arquivo;
    
    if (lixo[0]=='*') {
    
      atm.Ns++;
    
    }                                                               // Se ler '*' soma +1
  
  }                                                                 // conta quantos '*' tem em um arquivo enquanto não chegar no final do arquivo            
  
  printf("atm.Ns = %d\n", atm.Ns);
    
//----------------------------------------------------------------------------------------------------//  
                                                                    
  fclose(arq);                                                      //fecha arquivo
  
  
  printf("Simb: %s\n", atm.simb);
  printf("Z: %d\tA: %lf\trho: %lf\n", atm.Z, atm.A, atm.rho);
  printf("Coer:\n");
  for(i = 0; i < 6; i++){
    printf("%lf ", atm.coer[i]);} 
  printf("\n");
  printf("Incoer:\n");
  for(i = 0; i < 6; i++){
  printf("%lf ", atm.incoer[i]);}
  printf("\n");  
  return 0;
}






























