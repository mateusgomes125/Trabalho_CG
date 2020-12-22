
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip> 
#include <string> 

#include <fstream>
#include <Eigen/Dense>
#include "matplotlibcpp.h"
#include "PGM.h"
#include "PPM.h"
#include "Coord2D.h"
#include "Regiao.h"

using namespace std;
using namespace Eigen;

namespace plt = matplotlibcpp;

float graus2radianos(float g) {
	float r = g * (float)(M_PI / 180.0);

	return r;
}

void getMax(vector<float> c, int &mx){
	
	float f = *max_element(c.begin(), c.end());
	
	int v = (int) abs(round(f));
	if(v> mx)
		mx = v;
	
}

void transform2D( Matrix3f M, vector<Vector3f> coordsIn, vector<Vector3f> &coordsOut){
	coordsOut.clear();
	coordsOut = vector<Vector3f>(coordsIn.size());
	for(size_t i=0; i<coordsIn.size(); i++)
		coordsOut[i] = M * coordsIn[i];
	
}

void transformacao2D(Matrix3f M, PGM *entrada, PGM *saida) {    
    int idxSaida, idxEntrada;		
		for (int x = 0; x < entrada->getL(); x++) {
			for (int y = 0; y < entrada->getA(); y++) {
				Vector3f pEntrada(x, y, 1);
				Vector3f p = M * pEntrada;
				Vector2f pSaida = p.hnormalized();

				int xSaida = (int) round( pSaida(0) );
				int ySaida = (int) round( pSaida(1) );

				if (saida->verificarCoordenada(xSaida, ySaida)) {
					idxEntrada = y * entrada->getL() + x;
					idxSaida = ySaida * saida->getL() + xSaida;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];
				}
			}
		}
	} 

void transformacao2DPPM(Matrix3f M, PPM *entrada, PPM *saida) {    
    int idxSaida, idxEntrada;		
		for (int x = 0; x < entrada->getL(); x++) {
			for (int y = 0; y < entrada->getA(); y++) {
				Vector3f pEntrada(x, y, 1);
				Vector3f p = M * pEntrada;
				Vector2f pSaida = p.hnormalized();

				int xSaida = (int) round( pSaida(0) );
				int ySaida = (int) round( pSaida(1) );
				if (saida->verificarCoordenada(xSaida, ySaida)) {
					idxEntrada = y * 3 * entrada->getL() + x * 3;
					idxSaida = ySaida * 3 * saida->getL() + xSaida * 3;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];

					idxEntrada++;
					idxSaida++;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];

					idxEntrada++;
					idxSaida++;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];
				}
			}
		}
	}   

void transformacao2Dinv(Matrix3f M, PGM *entrada, PGM *saida) {   
    int idxSaida, idxEntrada;		
		for (int x = 0; x < saida->getL(); x++) {
			for (int y = 0; y < saida->getA(); y++) {
				Vector3f pSaida(x, y, 1);
				Vector3f p = M * pSaida;
				Vector2f pEntrada = p.hnormalized();
				int xEntrada = (int) round( pEntrada(0) );
				int yEntrada = (int) round( pEntrada(1) );
				if (entrada->verificarCoordenada(xEntrada, yEntrada)) {
					idxSaida = y * saida->getL() + x;
					idxEntrada = yEntrada * entrada->getL() + xEntrada;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];
				}
			}
		}
	} 


	void transformacao2DinvPPM(Matrix3f M, PPM *entrada, PPM *saida) {   
    int idxSaida, idxEntrada;		
		for (int x = 0; x < saida->getL(); x++) {
			for (int y = 0; y < saida->getA(); y++) {
				Vector3f pSaida(x, y, 1);
				Vector3f p = M * pSaida;
				Vector2f pEntrada = p.hnormalized();
				int xEntrada = (int) round( pEntrada(0) );
				int yEntrada = (int) round( pEntrada(1) );
				if (entrada->verificarCoordenada(xEntrada, yEntrada)) {
					idxSaida = y * 3 * saida->getL() + x * 3;
					idxEntrada = yEntrada* 3 * entrada->getL() + xEntrada * 3;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];

					idxEntrada++;
					idxSaida++;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];

					idxEntrada++;
					idxSaida++;
					saida->pixels[idxSaida] = entrada->pixels[idxEntrada];

				}
			}
		}
	}   

void hnormalized(vector<Vector3f>coords, vector<float> &x, vector<float> &y){
	for(size_t i=0; i<coords.size(); i++){
		
		Vector2f p = coords[i].hnormalized();
		x.push_back(p[0]);
		y.push_back(p[1]);
	}	
	
}
//----------------------------------FUNÇÕES DE TRANSFORMAÇÕES------------------------------------------------
Matrix3f computarMRotacao(float angulo){

	float rad = graus2radianos(angulo);
	Matrix3f T = Matrix3f::Identity();
	float c = cos(rad);
	float s = sin(rad);
		T(0, 0) = c;
		T(0, 1) = -s;
		T(1, 0) = s;
		T(1, 1) = c;
	return T;
}
Matrix3f computarMRotacaoInv(float angulo){
	Matrix3f T = computarMRotacao(angulo);
	Matrix3f RInv = T.transpose();
	return RInv;

}
Matrix3f computarMEscala(float x, float y){
	Matrix3f T = Matrix3f::Identity();
	T(0, 0) = x;
	T(1, 1) = y;
	return T;
}
Matrix3f computarMEscalaInv(float x, float y){
	Matrix3f T = Matrix3f::Identity();
	T(0, 0) = 1/x;
	T(1, 1) = 1/y;
	return T;
}
Matrix3f computarMTranslacao(int x, int y){
	Matrix3f T = Matrix3f::Identity();
	T(0, 2) = x;
	T(1, 2) = y;
	return T;
}
Matrix3f computarMTranslacaoInv(int x, int y){
	Matrix3f T = Matrix3f::Identity();
	T(0, 2) = -x;
	T(1, 2) = -y;
	return T;
}
Matrix3f computarMReflexaoX(){
	Matrix3f T = Matrix3f::Identity();
	T(1, 1) = -1;
	return T;		
}
Matrix3f computarMReflexaoY(){
	Matrix3f T = Matrix3f::Identity();
	T(0, 0) = -1;
	return T;		
}
Matrix3f computarMCizalhamentoH(float ciFator){	
	Matrix3f T = Matrix3f::Identity();
	T(0, 1) = ciFator;
	return T;
}	
Matrix3f computarMCizalhamentoHInv(float ciFator){	
	Matrix3f T = Matrix3f::Identity();
	T(0, 1) = -ciFator;
	return T;
}		
Matrix3f computarMCizalhamentoV(float ciFator){	
	Matrix3f T = Matrix3f::Identity();
	T(1, 0) = ciFator;
	return T;
}
Matrix3f computarMCizalhamentoVInv(float ciFator){	
	Matrix3f T = Matrix3f::Identity();
	T(1, 0) = -ciFator;
	return T;
}		

int main(void){

	std::streambuf *cinbuf = std::cin.rdbuf(); //salvar o buffer atual!
    ifstream file("entrada.txt");
    if (file.is_open()){
		cin.rdbuf(file.rdbuf()); //redirecionar std::cin para entrada.txt!	
	}	

	//VARIAVEIS
	string img, map, codigo, ci, re, rot;
	float angulo, tx, ty, rx, ry, ciFator;

	
	Matrix3f T;

	
	vector<Matrix3f> transformacoes;
	Matrix3f M = Matrix3f::Identity();

cout << "identidade: " << M;

	while(true){
		cin >> codigo;
		cout << "\n" << codigo;
			if (!std::cin)
						break;
			if (codigo.size() > 0 && codigo[0] == '#') {
						string aux;
						getline(cin,aux);	
			}	
			else if(codigo == "IMG"){//TIPO DE MAPEAMENTO
				cin >> img;				
			}
			else if(codigo == "MAP"){//TIPO DE MAPEAMENTO
				cin >> map;				
			}
			//-----------------------------TRANSFORMAÇÕES-----------------------------------------
			else if(codigo == "R"){				
				cin >> angulo;
				if(map == "INV"){
					T = computarMRotacaoInv(angulo);
					transformacoes.push_back(T);
				}					
				else{
						T = computarMRotacao(angulo);
					transformacoes.push_back(T);
				}	
				cout << "\nentrei na rotaçao!\n";				
			}
			else if(codigo == "T"){
				cin >> tx;
				cin >> ty;
				if(map == "INV")
					transformacoes.push_back(computarMTranslacaoInv(tx, ty));
				else{
					transformacoes.push_back(computarMTranslacao(tx, ty));				
				}					
			}
			else if(codigo == "S"){
				cin >> rx;
				cin >> ry;
				if(map == "INV"){
					T = computarMEscalaInv(rx, ry);
					transformacoes.push_back(T);
				}
					
				else{
					T = computarMEscala(rx, ry);
					transformacoes.push_back(T);
				}
			}
			else if(codigo == "CI"){
					cin >> ci;
					if(ci == "H"){
						cout << "entrou h!";
						cin >> ciFator;
						if(map == "INV")
							transformacoes.push_back(computarMCizalhamentoHInv(ciFator));
						else{
							transformacoes.push_back(computarMCizalhamentoH(ciFator));
						}
					}
					else if(ci == "V"){
						cout << "entrou v!";
						cin >> ciFator;
						if(map == "INV")
							transformacoes.push_back(computarMCizalhamentoVInv(ciFator));
						else{
							transformacoes.push_back(computarMCizalhamentoV(ciFator));
						}
							
					}
			}
			else if(codigo == "RE"){
					cin >> re;
					if(re == "H"){
						T = computarMReflexaoX();
						transformacoes.push_back(T);
						M = M * T;
					}
					else if(re == "V"){
						T = computarMReflexaoY();
						transformacoes.push_back(T);
						M = M * T;
					}
			}
	}			
	file.close(); //fechar o arquivo!
	std::cin.rdbuf(cinbuf);   //restabelecer a entrada padrão novamente!

//-----------------------------------VERIFICAÇÃO DO TIPO DA IMAGEM----------------------------------------------------------
PGM* pgmE = NULL, *pgmS = NULL;
PPM* ppmE = NULL, *ppmS = NULL;
Coord2D c;

// ----------------------------------VERIFICAÇÃO DE IMAGEM PGM -----------------------------------------------------------------
if(PGM::verificarTipoP2(img)){
	pgmE = new PGM();
	if(!pgmE->ler(img)){//-----------LENDO IMAGEM-------------------------------------------------------------------------------
		cout << "Imagem não existe!"<< endl;
		return EXIT_SUCCESS;
	}
	pgmS = new PGM(pgmE->getL(), pgmE->getA());
	c = pgmE->computarCentro();

	if(map == "INV"){//--------------VERICANDO MAPEAMENTO INVERSO------------------------------------------------------------------------
		
		cout << "\n MAPEAMENTO INVERSO EM IMAGEM PGM!";		
		Matrix3f TX = computarMTranslacaoInv(c.x, c.y);
		transformacoes.push_back(TX);

		for(int i = transformacoes.size()-1; i >=0 ; i--)M = transformacoes[i] * M;	//COLETANDO TRANSFORMAÇÕES
		M = computarMTranslacao(c.x, c.y) * M;
		cout << pgmE->getL();
		transformacao2Dinv(M, pgmE, pgmS);	//APLICANDO A TRANSFORMAÇÃO		
		cout << "\nMATRIZ DE TRANSFORMACAO: " << M;
	}
	 else{//------------------------VERIFICANDO MAPEAMENTO DIRETO---------------------------------------------------------------------

		 cout << "\n MAPEAMENTO DIRETO EM IMAGEM PGM!";	
		 M = M * computarMTranslacao(c.x, c.y);
		 for(int i = 0; i < transformacoes.size() ; i++)M = M * transformacoes[i];	//COLETANDO TRANSFORMAÇÕES
		 
		 T = computarMTranslacaoInv(c.x, c.y);
		 M = M * T;
	    cout << "Matriz " << M;
		transformacao2D(M, pgmE, pgmS); //APLICANDO A TRANSFORMAÇÃO
		cout << "\nMATRIZ DE TRANSFORMACAO: " << M;
	}	
pgmS->gravar("lenatransformada.pgm"); //GRAVANDO IMAGEM
}else if(PPM::verificarTipoP3(img)){//--------VERIFICANDO TIPO DE IMAGEM ---------------------------------------------------------------
	ppmE = new PPM();
	if(!ppmE->ler(img)){//--------------------LENDO IAMGEM------------------------------------------------------------------------
		cout << "Imagem não existe!"<< endl;
		return EXIT_SUCCESS;
	}
	c = ppmE->computarCentro();
	ppmS = new PPM(ppmE->getL(), ppmE->getA());

	if(map == "INV"){//----------------------VERIFICANDO MAPEAMENTO INVERSO----------------------------------------------------------------------		
		
		cout << "\n MAPEAMENTO INVERSO EM IMAGEM PPM!";	
		Matrix3f TX = computarMTranslacaoInv(c.x, c.y);
		transformacoes.push_back(TX);
		
		for(int i = transformacoes.size()-1; i >=0 ; i--)M = transformacoes[i] * M;	//COLETANDO TRANSFORMAÇÕES
		
		M = computarMTranslacao(c.x, c.y) * M;
		transformacao2DinvPPM(M, ppmE, ppmS);	//APLICANDO A TRANSFORMAÇÃO	
		cout << "\nMATRIZ DE TRANSFORMACAO: " << M;
	}
	 else{//--------------------------------VERIFICANDO MAPEAMENTO DIRETO--------------------------------------------------------------------------------
		 cout << "\n MAPEAMENTO DIRETO EM IMAGEM PPM!";	
		 M = M * computarMTranslacao(c.x, c.y);
		
		 for(int i = 0; i < transformacoes.size() ; i++)M = M * transformacoes[i];	//COLETANDO TRANSFORMAÇÕES
		 
		 T = computarMTranslacaoInv(c.x, c.y);
		 M = M * T;
		transformacao2DPPM(M, ppmE, ppmS);//APLICANDO A TRANSFORMAÇÃO
		cout << "\nMATRIZ DE TRANSFORMACAO: " << M;
	}	
	ppmS->gravar("spidertransformado.ppm"); //GRAVANDO IMAGEM
}

    return EXIT_SUCCESS;
}