#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
  char   simb[3];         // simbolo quimico string 3 caracteres char;  
  int    Z;               // Z                                   int;   
  double A;               // A                                   double;
  double rho;             // P                                   double;
  double coer[6];         // Coeficientes Coerente[6]            double;
  double incoer[6];       // Coeficientes Incoerente[6]          double;
  
} atomo;



int main(void){
  
  FILE *arq;                                                    // ponteiro;
  atomo atm;                                                    // variavel do atomo;
  char lixo[64];  
  int i;                                               // variavel;
  
  arq = fopen("~/Documentos/vitor_modelagem/dbase/C.txt","rt"); // abre arquivo;
//----------------------------------------------------------------------------------------------------// 
  fscanf(arq,"%s\n", atm.simb);                                 // lendo o simbolo;
  
  fgets(lixo, 64, arq);                                         // descarte da linha 2;
  
  fscanf(arq, "%d %lf %lf ", &atm.Z, &atm.A, &atm.rho);         // lendo Z, A, rho;
  
  fgets(lixo, 64, arq); 
  
  for(i = 0; i < 6; i++){
    fscanf(arq,"%lf ", &atm.coer[i]);
    }     // lendo Coeficientes Coerente[6]  
   
  fgets(lixo, 64, arq);                                         // descarte da linha 6;  
  
  for(i = 0; i < 6; i++){
  fscanf(arq,"%lf ", &atm.incoer[i]);
  }   // lendo Coeficientes inCoerente[6] 
//----------------------------------------------------------------------------------------------------//  
  fclose(arq);                                                  //fecha arquivo
  
  
  printf("Simb: %s\n", atm.simb);
  printf("Z: %d\tA: %lf\trho: %lf\n", atm.Z, atm.A, atm.rho);
  for(i = 0; i < 6; i++){
    printf("Coer: %lf ", atm.coer[i]);} 
  printf("\n");
  for(i = 0; i < 6; i++){
  printf("Incoer: %lf ", atm.incoer[i]);}
  printf("\n");  
  
  return 0;
}
