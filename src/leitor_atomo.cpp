#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

//---------------------------------------------------class fluorescencia-----------------------------

class fluorescencia {
public: // Tornado público para permitir leitura e impressão direta
    char   linha[6];       // Xray (char[6])
    double prob;          // Probabilidade
    double proba;         // Probabilidade acumulada
    double energia;       // Energia hv
};

//----------------------------------------------------class subnivel-----------------------------------

class subnivel {
public: // Tornado público para ser acessado pela classe atomo
    char      simb[6];    // Nome do subnivel
    double    Eb;         // Energia de ligacao (Binding Energy - Eb)
    double    rs;         // Raio de espalhamento / fator de forma
    double    w;          // Largura de linha / peso estatistico
    double    coef[6];    // Coeficientes para secao de choque fotoeletrica

    int Nf;
    fluorescencia* xrf;

    // Metodos
    double SecaoChoqueAbsFotoeletrica(double hv) {
        double s = 0, p = 1, lnE = log(hv); // Corrigido de lh(hv) para log(hv)

        for (int i = 0; i < 6; i++) {
            s += coef[i] * p;
            p *= lnE;
        }

        return exp(s);
    }
};

//-------------------------------------------------Class atomo-----------------------------------

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
    subnivel* ene;        // Ponteiro para o array dinamico de subniveis

public:
    // Encapsulamento / Getters e Setters
    void          Simb(const char* S) { strcpy(simb, S); } // Corrigido parâmetro S
    const char* Simb(void) { return simb; }
    void          Z(int n) { Na = n; }           // Removido parêntese extra
    int           Z(void) { return Na; }
    void          A(double m) { Ma = m; }           // Removido parêntese extra
    double        A(void) { return Ma; }
    void          Rho(double r) { rho = r; }          // Removido parêntese extra
    double        Rho(void) { return rho; }
    int           GetNs(void) { return Ns; }
    subnivel* GetSubniveis(void) { return ene; }

    void LerArquivo(const char* Nome_arq);
    void ImprimirDados(void);

    // Metodos de Calculo de Secao de Choque
    double SecaoChoqueElastico(double hv) {
        double s = 0, p = 1, lnE = log(hv);
        for (int i = 0; i < 6; i++) {
            s += coer[i] * p;
            p *= lnE;
        }
        return exp(s);
    }

    double SecaoChoqueIneElastico(double hv) {
        double s = 0, p = 1, lnE = log(hv);
        for (int i = 0; i < 6; i++) {
            s += incoer[i] * p;
            p *= lnE;
        }
        return exp(s);
    }

    double SecaoChoqueAbsFotoeletrica(double hv) {
        double s = 0;
        for (int i = 0; i < Ns; i++) { // Corrigido a inicialização de i
            if (ene[i].Eb <= hv) {
                s += ene[i].SecaoChoqueAbsFotoeletrica(hv);
            }
        }
        return s;
    }
};

// --- IMPLEMENTACAO DO METODO LER ARQUIVO ---
void atomo::LerArquivo(const char* Nome_arq) {
    FILE* arq = fopen(Nome_arq, "rt");
    if (!arq) {
        printf("Erro ao abrir arquivo %s\n", Nome_arq);
        return;
    }

    char lixo[64];
    int i, j, pos;

    // Dados do Atomo
    fscanf(arq, "%s\n", simb);
    fgets(lixo, 64, arq);
    fscanf(arq, "%d %lf %lf ", &Na, &Ma, &rho); // Ajustado para as variaveis da classe (Na, Ma, rho)
    fgets(lixo, 64, arq);

    for (i = 0; i < 6; i++) fscanf(arq, "%lf ", &coer[i]);
    fgets(lixo, 64, arq);
    for (i = 0; i < 6; i++) fscanf(arq, "%lf ", &incoer[i]);
    fgets(lixo, 64, arq);

    // Dados do Subnivel
    pos = ftell(arq);
    Ns = 0;

    while (fgets(lixo, 64, arq) != NULL) {
        if (lixo[0] == '*') Ns++;
    }

    ene = new subnivel[Ns];
    fseek(arq, pos, SEEK_SET);

    for (i = 0; i < Ns; i++) {
        fscanf(arq, "%s ", ene[i].simb);
        fgets(lixo, 64, arq);
        fscanf(arq, "%lf %lf %lf ", &ene[i].Eb, &ene[i].rs, &ene[i].w);
        fgets(lixo, 64, arq);

        for (j = 0; j < 6; j++) {
            fscanf(arq, "%lf ", &ene[i].coef[j]);
        }

        // Fluorescencia
        fgets(lixo, 64, arq);
        fgets(lixo, 64, arq);

        pos = ftell(arq);
        ene[i].Nf = 0;

        while (fgets(lixo, 64, arq) != NULL) {
            if (lixo[0] == '*') break;
            ene[i].Nf++;
        }

        ene[i].xrf = new fluorescencia[ene[i].Nf];
        fseek(arq, pos, SEEK_SET);

        for (j = 0; j < ene[i].Nf; j++) {
            fscanf(arq, "%s %lf %lf %lf ",
                ene[i].xrf[j].linha,
                &ene[i].xrf[j].prob,
                &ene[i].xrf[j].proba,
                &ene[i].xrf[j].energia);
        }
        fgets(lixo, 64, arq);
    }

    fclose(arq);
}

// --- METODO DE IMPRESSAO DENTRO DA CLASSE ATOMO ---
void atomo::ImprimirDados(void) {
    printf("Simb: %s\n", simb);
    printf("Z: %d\tA: %lf\trho: %lf\n", Na, Ma, rho);

    printf("Coer:\n");
    for (int i = 0; i < 6; i++) printf("%lf ", coer[i]);
    printf("\nIncoer:\n");
    for (int i = 0; i < 6; i++) printf("%lf ", incoer[i]);
    printf("\n\n");

    for (int i = 0; i < Ns; i++) {
        printf("Subnivel: %s\n", ene[i].simb);
        printf("Eb: %lf | rs: %lf | w: %lf\n", ene[i].Eb, ene[i].rs, ene[i].w);

        printf("coef: ");
        for (int j = 0; j < 6; j++) printf("%lf ", ene[i].coef[j]);
        printf("\n");

        if (ene[i].Nf > 0) {
            printf("\nlinha\tprob\tproba\tenergia\n");
            for (int j = 0; j < ene[i].Nf; j++) {
                printf("%s %lf %lf %lf\n",
                    ene[i].xrf[j].linha,
                    ene[i].xrf[j].prob,
                    ene[i].xrf[j].proba,
                    ene[i].xrf[j].energia);
            }
            printf("\n");
        }
    }
}

// --- FUNCAO PRINCIPAL ---
int main(void) {
    atomo atm;

    // Chama a leitura do arquivo apontando para o caminho correto
    atm.LerArquivo("data/B.txt");

    // Exibe os dados utilizando o metodo da classe
    atm.ImprimirDados();

    return 0;
}
