#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

//---------------------------------------------------class fluorescencia-----------------------------

class fluorescencia {
private:
    char   linha[16];   // Aumentado para evitar estouro de buffer
    double prob;        // Probabilidade individual
    double proba;       // Probabilidade acumulada
    double energia;     // Energia (keV)

public:
    fluorescencia() : prob(0.0), proba(0.0), energia(0.0) {
        linha[0] = '\0';
    }

    void        Linha(const char* nome) { strncpy(linha, nome, sizeof(linha) - 1); linha[sizeof(linha) - 1] = '\0'; }
    char*       Linha(void)             { return linha; }
    
    void        Prob(double e)          { prob = e; }
    double      Prob(void) const        { return prob; }
    
    void        ProbA(double e)         { proba = e; }
    double      ProbA(void) const       { return proba; }
    
    void        E(double e)             { energia = e; }
    double      E(void) const           { return energia; }
};                     

//----------------------------------------------------class subnivel-----------------------------------

class subnivel {
private:
    char          simb[16];   
    double        Eb;         // Energia de ligação
    double        rs;         // Raio de espalhamento
    double        w;          // Largura / Rendimento
    double        coef[6];    // Coeficientes do polinômio fotoelétrico
    
    int           Nf;         // Quantidade de linhas
    fluorescencia *xrf;       // Array dinâmico de fluorescência
    
public: 
    subnivel() : Eb(0.0), rs(0.0), w(0.0), Nf(0), xrf(nullptr) {
        simb[0] = '\0';
        for (int i = 0; i < 6; i++) coef[i] = 0.0;
    }

    // Destrutor para liberação da memória
    ~subnivel() {
        if (xrf != nullptr) {
            delete[] xrf;
            xrf = nullptr;
        }
    }

    void        Simb(const char *nome)  { strncpy(simb, nome, sizeof(simb) - 1); simb[sizeof(simb) - 1] = '\0'; }
    const char* Simb(void) const        { return simb; }
    
    void        GetEb(double e)         { Eb = e; }
    double      GetEb(void) const       { return Eb; }
    
    void        Rs(double r)            { rs = r; }
    double      Rs(void) const          { return rs; }
    
    void        W(double e)             { w = e; }
    double      W(void) const           { return w; }
    
    void GetNf(int n) { 
        if (xrf != nullptr) delete[] xrf;
        Nf = n; 
        xrf = (Nf > 0) ? new fluorescencia[Nf] : nullptr; 
    }
    int  GetNf(void) const              { return Nf; }
    
    void XRF(int pos, const char* l, double p, double pa, double e) {
        if (pos >= 0 && pos < Nf && xrf != nullptr) {
            xrf[pos].Linha(l);
            xrf[pos].Prob(p);
            xrf[pos].ProbA(pa);
            xrf[pos].E(e);
        }
    }
    
    char*  Linha(int pos) const { return xrf[pos].Linha(); }
    double Prob (int pos) const { return xrf[pos].Prob(); }
    double ProbA(int pos) const { return xrf[pos].ProbA(); }
    double E    (int pos) const { return xrf[pos].E(); }
    
    void Coef(double *c) {
        for (int i = 0; i < 6; i++) coef[i] = c[i];
    }
                                            
    double Coef(int pos) const { return coef[pos]; }
    
    double SecaoChoqueAbsFotoeletrica(double hv) const {
        double s = 0.0, p = 1.0, lnE = log(hv);
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
    char      simb[8];     
    int       Na;          // Z
    double    Ma;          // A
    double    rho;         // Densidade
    double    coer[6];     // Coeficientes Coerentes (Rayleigh)
    double    incoer[6];   // Coeficientes Incoerentes (Compton)  
    int       Ns;          // Número de subníveis
    subnivel  *ene;        // Array dinâmico de subníveis
    
public:  
    atomo() : Na(0), Ma(0.0), rho(0.0), Ns(0), ene(nullptr) {
        simb[0] = '\0';
        for(int i = 0; i < 6; i++) {
            coer[i] = 0.0;
            incoer[i] = 0.0;
        }
    }

    // Destrutor da classe atomo
    ~atomo() {
        if (ene != nullptr) {
            delete[] ene;
            ene = nullptr;
        }
    }

    void        Simb(const char *S) { strncpy(simb, S, sizeof(simb) - 1); simb[sizeof(simb) - 1] = '\0'; }
    const char* Simb(void) const    { return simb; }
    
    void        Z(int n)            { Na = n; }
    int         Z(void) const       { return Na; }
    
    void        A(double m)         { Ma = m; }
    double      A(void) const       { return Ma; }
    
    void        Rho(double r)       { rho = r; }
    double      Rho(void) const     { return rho; }
    
    void GetNs(int n) { 
        if (ene != nullptr) delete[] ene;
        Ns = n; 
        ene = (Ns > 0) ? new subnivel[Ns] : nullptr; 
    }
    int  GetNs(void) const          { return Ns; }
    
    void LerArquivo(const char *Nome_arq);
    
    double SecaoChoqueElastico(double hv) const {
        double s = 0.0, p = 1.0, lnE = log(hv);
        for (int i = 0; i < 6; i++) {
            s += coer[i] * p;
            p *= lnE;
        }
        return exp(s);
    }
    
    double SecaoChoqueIneElastico(double hv) const {
        double s = 0.0, p = 1.0, lnE = log(hv);
        for (int i = 0; i < 6; i++) {
            s += incoer[i] * p;            
            p *= lnE;
        }
        return exp(s);
    }
    
    double SecaoChoqueAbsFotoeletrica(double hv) const {
        double s = 0.0;
        for (int i = 0; i < Ns; i++) { 
            if (ene[i].GetEb() <= hv) {
                s += ene[i].SecaoChoqueAbsFotoeletrica(hv);
            }
        }
        return s;
    }
};

// --- IMPLEMENTAÇÃO DA LEITURA DO ARQUIVO ---

void atomo::LerArquivo(const char *Nome_arq) {
    FILE *arq = fopen(Nome_arq, "rt");
    if (!arq) {
        std::printf("Erro ao abrir arquivo: %s\n", Nome_arq);
        return;
    }

    char lixo[128], linha[16];
    int i, j, Nf;
    long pos;
    double e, r, w;
    double prob, proba, energia;
    double coef[6];

    if (fscanf(arq, "%7s\n", simb) != 1) { fclose(arq); return; }
    if (!fgets(lixo, sizeof(lixo), arq)) { fclose(arq); return; }
    if (fscanf(arq, "%d %lf %lf ", &Na, &Ma, &rho) != 3) { fclose(arq); return; }
    if (!fgets(lixo, sizeof(lixo), arq)) { fclose(arq); return; }
    
    for (i = 0; i < 6; i++) {
        if (fscanf(arq, "%lf ", &coer[i]) != 1) break;
    }
    if (!fgets(lixo, sizeof(lixo), arq)) { fclose(arq); return; }
    
    for (i = 0; i < 6; i++) {
        if (fscanf(arq, "%lf ", &incoer[i]) != 1) break;
    }
    if (!fgets(lixo, sizeof(lixo), arq)) { fclose(arq); return; }

    // Primeira Passagem: Contagem dos Subníveis
    pos = ftell(arq);
    Ns = 0;
    while (fgets(lixo, sizeof(lixo), arq) != NULL) {
        if (lixo[0] == '*') Ns++;
    }

    GetNs(Ns);
    fseek(arq, pos, SEEK_SET);

    // Segunda Passagem: Carregamento
    for (i = 0; i < Ns; i++) {
        if (fscanf(arq, "%127s ", lixo) != 1) break;
        ene[i].Simb(lixo);
        
        if (!fgets(lixo, sizeof(lixo), arq)) break;
        if (fscanf(arq, "%lf %lf %lf ", &e, &r, &w) != 3) break;
        
        ene[i].GetEb(e);
        ene[i].Rs(r);
        ene[i].W(w);
        
        if (!fgets(lixo, sizeof(lixo), arq)) break;
        
        for (j = 0; j < 6; j++) {
            if (fscanf(arq, "%lf ", &coef[j]) != 1) break;
        }
        ene[i].Coef(coef);
        
        if (!fgets(lixo, sizeof(lixo), arq)) break;
        if (!fgets(lixo, sizeof(lixo), arq)) break;
        
        pos = ftell(arq);
        Nf = 0;
        
        while (fgets(lixo, sizeof(lixo), arq) != NULL) {
            if (lixo[0] == '*') break;
            Nf++;
        }

        ene[i].GetNf(Nf);
        fseek(arq, pos, SEEK_SET);

        for (j = 0; j < ene[i].GetNf(); j++) {
            if (fscanf(arq, "%15s %lf %lf %lf ", linha, &prob, &proba, &energia) == 4) {
                ene[i].XRF(j, linha, prob, proba, energia);
            }
        }
        if (!fgets(lixo, sizeof(lixo), arq)) break;
    }

    fclose(arq);
}

// --- FUNÇÃO PRINCIPAL ---

int main(void) {
    double hv = 5.0;              // Energia do fóton incidente (keV)
    double x = 0.01;              // Espessura do meio (cm)
    
    atomo atm;
    atm.LerArquivo("/home/administrador/Documentos/vitor_modelagem/dbase/C.txt");
    
    // Cálculos de seções de choque (cm²/g)
    double sigma_e = atm.SecaoChoqueElastico(hv);
    double sigma_i = atm.SecaoChoqueIneElastico(hv);
    double tau     = atm.SecaoChoqueAbsFotoeletrica(hv);
    
    // Coeficiente de atenuação mássico total (\mu/\rho em cm²/g)
    double mu = sigma_e + sigma_i + tau;
    
    // Cálculo das probabilidades de interação no caminho x
    double Ptrans = exp(-mu * atm.Rho() * x);
    double Pee    = (mu > 0) ? (sigma_e * (1.0 - Ptrans) / mu) : 0.0;
    double Pei    = (mu > 0) ? (sigma_i * (1.0 - Ptrans) / mu) : 0.0;
    double Pabs   = (mu > 0) ? (tau     * (1.0 - Ptrans) / mu) : 0.0;

    std::printf("sigma_e: %lf\tsigma_i: %lf\ttau: %lf\n\n", sigma_e, sigma_i, tau);
    std::printf("Ptrans = %e\n", Ptrans);
    std::printf("Pee    = %e\n", Pee);
    std::printf("Pei    = %e\n", Pei);
    std::printf("Pabs   = %e\n", Pabs);

    return 0;
}
