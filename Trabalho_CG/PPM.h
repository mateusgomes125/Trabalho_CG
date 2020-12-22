//#pragma once

#ifndef PPM_H_INCLUDED
#define PPM_H_INCLUDED

#include<cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

#include <Eigen/Dense>
#include "Coord2D.h"
#include "Regiao.h"
using namespace std;

class PPM
{
public:
	PPM(int l, int a) {
		this->larg = l;
		this->alt = a;
		this->pixels = new unsigned char[this->larg * this->alt * 3];
		
		this->tipo = "P3";
		this->vmax = 1;

		int tam = this->larg * this->alt * 3;
		for (int i = 0; i < tam; i++) 
			this->pixels[i] = 0;
	}
	PPM() {
		this->pixels = NULL;
		this->larg = 0;
		this->alt = 0;
		this->tipo = "";
	}

	~PPM() {
		if (this->pixels)
			delete this->pixels;
	}
	
	static bool verificarTipoP3(string caminho){
		ifstream arq;
		string aux;
		arq.open(caminho);
			if(!arq.is_open()){
				return false;
			}
			if(!PPM::lerDado(arq, &aux)){
				return false;
			}
			if(aux == "P3")
				return true;
		return false;
	}

	/*void preencherRegiao(Coord2D pmin, Coord2D pmax, unsigned char cor) {
		if (!this->pixels){
			//cout << ...
			return;
		}
		for (int x = max(pmin.x, 0); x< min(pmax.x, this->larg-1); x++) {
			for (int y = max(pmin.y, 0); y < min(pmax.y, this->alt - 1); y++) {
				this->pixels[y * this->alt + x] = cor;
			}
		}

	}*/

	int calcMax(){
		if(!this->pixels)
			return 1;

			int tam = this->larg * this->alt * 3;
			int mx = this->pixels[0];
			for(int i = 0; i < tam; i++){
				if(this->pixels[i] > mx)
					mx = this->pixels[i];				
			}
		if(mx == 0)
			mx = 1;

			return mx;
			
		
	}

	static void map(PPM* entrada, PPM* saida, Regiao rEntrada, Regiao rSaida) {
		
		//percorrendo a entrada!
		for (int x = max(rEntrada.cMin.x, 0); x < min(rEntrada.cMax.x, entrada->larg - 1); x++) {
			for (int y = max(rEntrada.cMin.y, 0); y < min(rEntrada.cMax.y, entrada->alt - 1); y++) {
				
				float px = (x - rEntrada.cMin.x) / (float)(rEntrada.cMax.x - rEntrada.cMin.x);
				int xSaida = (int) round(px * (rSaida.cMax.x - rSaida.cMin.x) + rSaida.cMin.x);

				float py = (y - rEntrada.cMin.y) / (float)(rEntrada.cMax.y - rEntrada.cMin.y);
				int ySaida = (int)round(py * (rSaida.cMax.y - rSaida.cMin.y) + rSaida.cMin.y);

				saida->pixels[ySaida * saida->larg + xSaida] = entrada->pixels[y* entrada->larg + x];
			}
		}
	}

	static void mapInv(PPM* entrada, PPM* saida, Regiao rEntrada, Regiao rSaida) {

		//percorrendo a entrada!
		for (int x = max(rSaida.cMin.x, 0); x < min(rSaida.cMax.x, saida->larg - 1); x++) {
			for (int y = max(rSaida.cMin.y, 0); y < min(rSaida.cMax.y, saida->alt - 1); y++) {

				float px = (x - rSaida.cMin.x) / (float)(rSaida.cMax.x - rSaida.cMin.x);
				int xEntrada = (int)round(px * (rEntrada.cMax.x - rEntrada.cMin.x) + rEntrada.cMin.x);

				float py = (y - rSaida.cMin.y) / (float)(rSaida.cMax.y - rSaida.cMin.y);
				int yEntrada = (int)round(py * (rEntrada.cMax.y - rEntrada.cMin.y) + rEntrada.cMin.y);

				saida->pixels[y * saida->larg + x] = entrada->pixels[yEntrada * entrada->larg + xEntrada];
			}
		}
	}

	Coord2D computarCentro(){
		Coord2D centro;
		centro.x = this->alt/2;
		centro.y = this->larg/2;
		return centro;
	}

	int calcVmax() {
		if (!this->pixels)
			return 1;

		int tam = this->larg * this->alt;
		int mx = this->pixels[0];
		for (int i = 0; i < tam; i++) {
			if (this->pixels[i] > mx)
				mx = this->pixels[i];
		}

		if (mx == 0)
			mx = 1;

		return mx;
	}
	
	bool gravar(string caminho) {

		if (!this->pixels)
			return false;

		ofstream arq;

		arq.open(caminho);
		if (!arq.is_open()) {
			cout << "PPM: endereco do arquivo invalido\n";
			return false;
		}

		arq << this->tipo << endl;

		arq << this->larg << " " << this->alt << endl;

		arq << calcVmax() << endl;

		int tam = this->larg * this->alt * 3;
		for (int i = 0; i < tam; i++)
			arq << (int) this->pixels[i] << endl;

		arq.close();
		return true;
	}


	bool ler(string caminho) {
		ifstream arq;
		string aux;

		if (this->pixels)
			delete this->pixels;

		arq.open(caminho);
		if (!arq.is_open()) {
			cout << "PPM: endereco do arquivo invalido\n";
			return false;
		}

		if (!lerDado(arq, &aux)) {
			cout << "PPM: erro ao ler o tipo da imagem\n";
			return false;
		}
		this->tipo = aux;


		if (!lerDado(arq, &aux)) {
			cout << "PPM: erro ao ler a largura da imagem\n";
			return false;
		}
		this->larg = atoi(aux.c_str());

		if (!lerDado(arq, &aux)) {
			cout << "PPM: erro ao ler a largura da imagem\n";
			return false;
		}
		this->alt = atoi(aux.c_str());


		if (!lerDado(arq, &aux)) {
			cout << "PPM: erro ao ler a largura da imagem\n";
			return false;
		}
		this->vmax = atoi(aux.c_str());

		this->pixels = new unsigned char[this->larg * this->alt * 3];

		int i = 0;
		while (!arq.eof()) {
			if (!lerDado(arq, &aux)) {
				break;
			}

			this->pixels[i] = (unsigned char) atoi(aux.c_str());
			i++;

		}

		//cout << "i: " << i << endl;
		if (i != this->larg * this->alt * 3) {
			cout << "PPM: erro ao ler os pixels\n";
			return false;
		}

		//cout << this->tipo << endl;
		//cout << this->larg << endl;
		//cout << this->alt << endl;
		//cout << this->vmax << endl;


		return true;
	}

	int getL() {
		return larg;
	}
	int getA() {
		return alt;
	}
	bool verificarCoordenada(int x, int y) {
		if (x < 0 || y < 0 || x >= this->larg || y >= this->alt)
			return false;

		return true;
	}

	unsigned char* getPixels(){
		return pixels;
	}
unsigned char* pixels;	
private:
static	bool lerDado(ifstream &arq, string *dado) {
			string aux;
			while (!arq.eof()) {
				arq >> aux;
				if (aux.size() > 0 && aux[0] == '#') {
					std::getline(arq, aux);
				}else if(aux.size() > 0){
					*dado = aux;
					return true;
				}

				aux = "";
			}
			return false;
		}
	
	string tipo;
	int larg;
	int alt;
	int vmax;
	


};

#endif