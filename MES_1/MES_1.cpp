#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstddef>      
#include <vector>
#include <sstream>
#include <iomanip>

using namespace std;

const double eps = 1e-12; // stała przybliżenia zera
bool gauss(int n, double** AB, double* X);

struct  GlobalData {
	int SimulationTime;
	int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
	int NodesNumber;
	int ElementsNumber;

	GlobalData() {
		SimulationTime = 0;
		SimulationStepTime = 0;
		Conductivity = 0;
		Alfa = 0;
		Tot = 0;
		InitialTemp = 0;
		Density = 0;
		SpecificHeat = 0;
		NodesNumber = 0;
		ElementsNumber = 0;
	}

	GlobalData(vector <int> data) {
		SimulationTime = data[0];
		SimulationStepTime = data[1];
		Conductivity = data[2];
		Alfa = data[3];
		Tot = data[4];
		InitialTemp = data[5];
		Density = data[6];
		SpecificHeat = data[7];
		NodesNumber = data[8];
		ElementsNumber = data[9];
	}
};

struct  Node {
	double x;
	double y;
	int BC;

	Node(vector<double>node) {
		x = node[1];
		y = node[2];
		BC = 0;
	}
};

struct  Element {
	vector <int> ID;

	Element(vector<int> id) {
		for (int i = 1; i < id.size(); i++) {
			ID.push_back(id[i]);
		};
	}
};

struct Side {
	vector <int> pc;
	vector <double> ksi;
	vector <double> eta;
	double** N;

	double pomoc[2] = { -1,1 };
	// double  ksiTable[4] = {-1,1,-1,1};
	// double etaTable[4] = {-1,-1,1,1};

	Side() {}
	Side(int pc_numbers, double srt, vector <int> ksiTable, vector <int> etaTable, int numerBoku) {
		N = new double* [pc_numbers];
		for (int i = 0; i < pc_numbers; i++) {
			N[i] = new double[4];
			pc.push_back(i);
			if (numerBoku % 2 == 0) {
				ksi.push_back(ksiTable[i] * srt);
				eta.push_back(pomoc[numerBoku / 2]);
			}
			else {
				ksi.push_back(pomoc[1 / numerBoku]);
				eta.push_back(ksiTable[i] * srt);
			}
			if (numerBoku == 0) {
				N[i][0] = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
				N[i][1] = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
				N[i][2] = 0;
				N[i][3] = 0;
			}
			else if (numerBoku == 1) {
				N[i][0] = 0;
				N[i][1] = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
				N[i][2] = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
				N[i][3] = 0;
			}
			else if (numerBoku == 2) {
				N[i][0] = 0;
				N[i][1] = 0;
				N[i][2] = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
				N[i][3] = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
			}
			else {
				N[i][0] = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
				N[i][1] = 0;
				N[i][2] = 0;
				N[i][3] = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
			}
		}
		/*
		cout << "Bok nr: " << numerBoku << endl;
		for (int j = 0; j < pc_numbers; j++) {
			cout << "Nr PC: " << pc[j] << " ksi: " << ksi[j] << " eta: " << eta[j] << " N1: " << N[j][0] << " N2: " << N[j][1] << " N3: " << N[j][2] << " N4: " << N[j][3] << endl;
		}
		*/
	}
};

struct  Grid {
	vector <Node> ND;
	vector <Element> EL;
	int nEL;
	int nNd;
};


int main()
{
	GlobalData data;
	Grid grid;
	vector <int> globalData;
	vector <double> node;
	vector <int> element;
	Side boki[4];

	bool timeForNodes = false;
	bool timeForElements = false;
	bool timeForBC = false;

	fstream file;
	string line;
	size_t found;
	string word;

	file.open("Test1_4_4.txt", ios::in);
	//file.open("Test2_4_4_MixGrid.txt", ios::in);
	//file.open("Test3_31_31_kwadrat.txt", ios::in);

	
	if (file.good() == true)
	{
		while (!file.eof())
		{
			getline(file, line);

			stringstream sLine(line);

			if (line == "*Node")
				timeForNodes = true;
			else if (line == "*Element, type=DC2D4")
				timeForElements = true;
			else if (line == "*BC")
				timeForBC = true;
			else if (timeForBC == true)
			{
				while (getline(sLine, word, ',')) {
					for (int i = 0; i < grid.ND.size(); i++)
					{
						if (i == stoi(word) - 1)
							grid.ND[i].BC = 1;
					}
				}
			}
			else if (timeForElements == true && timeForBC == false)
			{
				while (getline(sLine, word, ',')) {
					element.push_back(stoi(word) - 1);
				}
				grid.EL.push_back(Element(element));
				element.clear();
			}
			else if (timeForNodes == true)
			{
				while (getline(sLine, word, ',')) {
					node.push_back(stod(word));
				}
				grid.ND.push_back(Node(node));
				node.clear();
			}
			else {
				found = line.find_last_of(" ");
				globalData.push_back(stoi(line.substr(found + 1)));
			}

		}
		data = GlobalData(globalData);
		grid.nEL = grid.EL.size();
		grid.nNd = grid.ND.size();

		file.close();

		cout << "NodesNumber: " << grid.nNd << endl;

		for (int i = 0; i < grid.ND.size(); i++)
		{
			cout << i << "\tx: " << grid.ND[i].x << " y: " << grid.ND[i].y << " BC: " << grid.ND[i].BC << endl;
		}

		cout << "ElementsNumber: " << grid.nEL << endl;

		for (int i = 0; i < grid.EL.size(); i++)
		{
			cout << i << "\t " << grid.EL[i].ID[0] << ", " << grid.EL[i].ID[1] << ", " << grid.EL[i].ID[2] << ", " << grid.EL[i].ID[3] << endl;
		}
	}
	else {
		cout << "Nie udało się otworzyć pliku!" << endl;
	}

	double** macierz_H_Globalna = new double* [data.NodesNumber];
	for (int i = 0; i < data.NodesNumber; ++i) {
		macierz_H_Globalna[i] = new double[data.NodesNumber];
		for (int j = 0; j < data.NodesNumber; ++j)
			macierz_H_Globalna[i][j] = 0.0;
	}

	double** macierz_C_Globalna = new double* [data.NodesNumber];
	for (int i = 0; i < data.NodesNumber; ++i) {
		macierz_C_Globalna[i] = new double[data.NodesNumber];
		for (int j = 0; j < data.NodesNumber; ++j)
			macierz_C_Globalna[i][j] = 0.0;
	}

	double* wektor_P_Globalny = new double[data.NodesNumber];
	for (int j = 0; j < data.NodesNumber; j++)
		wektor_P_Globalny[j] = 0.0;

	double* wektor_P_Globalny_kopia = new double[data.NodesNumber];
	for (int j = 0; j < data.NodesNumber; j++)
		wektor_P_Globalny_kopia[j] = 0.0;

	double** macierz_Gaussa = new double* [data.NodesNumber];
	for (int i = 0; i < data.NodesNumber; i++) {
		macierz_Gaussa[i] = new double[data.NodesNumber + 1];
		for (int j = 0; j < data.NodesNumber + 1; j++)
			macierz_Gaussa[i][j] = 0.0;
	}

	double* temp = new double[data.NodesNumber];
	for (int j = 0; j < data.NodesNumber; j++)
		temp[j] = data.InitialTemp;

	double* MinTemp = new double[20];
	double* MaxTemp = new double[20];
	for (int j = 0; j < 10; j++) {
		MinTemp[j] = data.InitialTemp;
		MaxTemp[j] = 0;
	}

	////////////////////////////////////////
	const int liczbaPunktow = 4;
	////////////////////////////////////////

	vector<double> waga;
	vector <int> etaTable;
	vector <int>  ksiTable;
	double srt;

	if (liczbaPunktow == 4) {
		waga = { 1,1 };
		srt = 1 / sqrt(3);
		ksiTable = { -1,1,1,-1 };
		etaTable = { -1,-1,1,1 };
	}
	else if (liczbaPunktow == 9) {
		waga = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };
		ksiTable = { -1,0,1,-1,0,1, -1,0,1 };
		etaTable = { -1,-1,-1,0,0,0,1,1, 1 };
		srt = sqrt(3.0 / 5.0);
	}

	for (int i = 0; i < 4; i++) {
		boki[i] = Side(sqrt(liczbaPunktow), srt, ksiTable, etaTable, i);
	}

	double tab_dN_dksi[liczbaPunktow][4] = {};
	double tab_dN_deta[liczbaPunktow][4] = {};
	double tab_N[liczbaPunktow][4] = {};

	//wartość pochodnych funkcji kształtu po ksi i eta dla każdego punktu całkowania
	for (int i = 0; i < liczbaPunktow; i++) {
		tab_dN_dksi[i][0] = -0.25 * (1 - srt * etaTable[i]);
		tab_dN_dksi[i][1] = 0.25 * (1 - srt * etaTable[i]);
		tab_dN_dksi[i][2] = 0.25 * (1 + srt * etaTable[i]);
		tab_dN_dksi[i][3] = -0.25 * (1 + srt * etaTable[i]);

		tab_dN_deta[i][0] = -0.25 * (1 - srt * ksiTable[i]);
		tab_dN_deta[i][1] = -0.25 * (1 + srt * ksiTable[i]);
		tab_dN_deta[i][2] = 0.25 * (1 + srt * ksiTable[i]);
		tab_dN_deta[i][3] = 0.25 * (1 - srt * ksiTable[i]);

		//wartość funkcji kształtu dla każdego punktu całkowania
		tab_N[i][0] = 0.25 * (1 - srt * ksiTable[i]) * (1 - srt * etaTable[i]);
		tab_N[i][1] = 0.25 * (1 + srt * ksiTable[i]) * (1 - srt * etaTable[i]);
		tab_N[i][2] = 0.25 * (1 + srt * ksiTable[i]) * (1 + srt * etaTable[i]);
		tab_N[i][3] = 0.25 * (1 - srt * ksiTable[i]) * (1 + srt * etaTable[i]);
	}

	/*
	cout << "dNi / ksi" << endl;
	for (int i = 0; i < liczbaPunktow; i++) {
		for (int j = 0; j < 4; j++)
			cout << tab_dN_dksi[i][j] << " ";
		cout << endl;
	}
	cout << endl;

	cout << "dNi / eta" << endl;
	for (int i = 0; i < liczbaPunktow; i++) {
		for (int j = 0; j < 4; j++)
			cout << tab_dN_deta[i][j] << " ";
		cout << endl;
	}
	cout << endl;
   */

	int numerElementu = 0;

	for (Element element : grid.EL) {

		numerElementu++;

		//cout << "\nWezly w elemencie nr " << numerElementu << ":\n" << element.ID[0] + 1 << ", " << grid.ND[element.ID[0]].x << ", " << grid.ND[element.ID[0]].y << "\n" << element.ID[1] + 1 << ", " << grid.ND[element.ID[1]].x << ", " << grid.ND[element.ID[1]].y << "\n" << element.ID[2] + 1 << ", " << grid.ND[element.ID[2]].x << ", " << grid.ND[element.ID[2]].y << "\n" << element.ID[3] + 1 << ", " << grid.ND[element.ID[3]].x << ", " << grid.ND[element.ID[3]].y << endl;

		double ID[2][4] = { {grid.ND[element.ID[0]].x,grid.ND[element.ID[1]].x,grid.ND[element.ID[2]].x,grid.ND[element.ID[3]].x},{grid.ND[element.ID[0]].y,grid.ND[element.ID[1]].y,grid.ND[element.ID[2]].y,grid.ND[element.ID[3]].y} }; //współrzędne węzłów

		//MACIERZ H I C

		double tab_jakobianow[liczbaPunktow][4] = {}; //macierz jakobianu odwrotnego
		double det[liczbaPunktow] = {}; // wyznacznik macierzy jacobiego

		//składowe jakobianu dla każdego punktu całkowania
		for (int i = 0; i < liczbaPunktow; i++) {
			tab_jakobianow[i][0] = tab_dN_deta[i][0] * ID[1][0] + tab_dN_deta[i][1] * ID[1][1] + tab_dN_deta[i][2] * ID[1][2] + tab_dN_deta[i][3] * ID[1][3];
			tab_jakobianow[i][1] = -1 * (tab_dN_dksi[i][0] * ID[1][0] + tab_dN_dksi[i][1] * ID[1][1] + tab_dN_dksi[i][2] * ID[1][2] + tab_dN_dksi[i][3] * ID[1][3]);
			tab_jakobianow[i][2] = -1 * (tab_dN_deta[i][0] * ID[0][0] + tab_dN_deta[i][1] * ID[0][1] + tab_dN_deta[i][2] * ID[0][2] + tab_dN_deta[i][3] * ID[0][3]);
			tab_jakobianow[i][3] = tab_dN_dksi[i][0] * ID[0][0] + tab_dN_dksi[i][1] * ID[0][1] + tab_dN_dksi[i][2] * ID[0][2] + tab_dN_dksi[i][3] * ID[0][3];

			det[i] = tab_jakobianow[i][3] * tab_jakobianow[i][0] - tab_jakobianow[i][1] * tab_jakobianow[i][2];
		}

		double tab_dN_dx[liczbaPunktow][4] = {}; //dN[1,..9]/dx 
		double tab_dN_dy[liczbaPunktow][4] = {}; //dN[1,..9]/dy

		for (int i = 0; i < liczbaPunktow; i++) {
			tab_dN_dx[i][0] = 1 / det[i] * tab_jakobianow[i][0] * tab_dN_dksi[i][0] + 1 / det[i] * tab_jakobianow[i][1] * tab_dN_deta[i][0];
			tab_dN_dx[i][1] = 1 / det[i] * tab_jakobianow[i][0] * tab_dN_dksi[i][1] + 1 / det[i] * tab_jakobianow[i][1] * tab_dN_deta[i][1];
			tab_dN_dx[i][2] = 1 / det[i] * tab_jakobianow[i][0] * tab_dN_dksi[i][2] + 1 / det[i] * tab_jakobianow[i][1] * tab_dN_deta[i][2];
			tab_dN_dx[i][3] = 1 / det[i] * tab_jakobianow[i][0] * tab_dN_dksi[i][3] + 1 / det[i] * tab_jakobianow[i][1] * tab_dN_deta[i][3];

			tab_dN_dy[i][0] = 1 / det[i] * tab_jakobianow[i][2] * tab_dN_dksi[i][0] + 1 / det[i] * tab_jakobianow[i][3] * tab_dN_deta[i][0];
			tab_dN_dy[i][1] = 1 / det[i] * tab_jakobianow[i][2] * tab_dN_dksi[i][1] + 1 / det[i] * tab_jakobianow[i][3] * tab_dN_deta[i][1];
			tab_dN_dy[i][2] = 1 / det[i] * tab_jakobianow[i][2] * tab_dN_dksi[i][2] + 1 / det[i] * tab_jakobianow[i][3] * tab_dN_deta[i][2];
			tab_dN_dy[i][3] = 1 / det[i] * tab_jakobianow[i][2] * tab_dN_dksi[i][3] + 1 / det[i] * tab_jakobianow[i][3] * tab_dN_deta[i][3];
		}

		double tab_1[4][4] = {}; //pomoc przy mnożeniu wektorów
		double tab_2[4][4] = {}; //pomoc przy mnożeniu wektorów 
		double tab_punktow_H[liczbaPunktow][4][4] = {}; //macierze H dla punktów całkowania
		double tab_punktow_C[liczbaPunktow][4][4] = {}; //macierze C dla punktów całkowania

		for (int j = 0; j < liczbaPunktow; j++) {
			for (int i = 0; i < 4; i++) {
				tab_1[i][0] = tab_dN_dx[j][i] * tab_dN_dx[j][0];
				tab_1[i][1] = tab_dN_dx[j][i] * tab_dN_dx[j][1];
				tab_1[i][2] = tab_dN_dx[j][i] * tab_dN_dx[j][2];
				tab_1[i][3] = tab_dN_dx[j][i] * tab_dN_dx[j][3];

				tab_2[i][0] = tab_dN_dy[j][i] * tab_dN_dy[j][0];
				tab_2[i][1] = tab_dN_dy[j][i] * tab_dN_dy[j][1];
				tab_2[i][2] = tab_dN_dy[j][i] * tab_dN_dy[j][2];
				tab_2[i][3] = tab_dN_dy[j][i] * tab_dN_dy[j][3];

				//obliczanie macierzy H w poszczegółnych punktów całkowania
				tab_punktow_H[j][i][0] = (tab_1[i][0] + tab_2[i][0]) * data.Conductivity * det[j];
				tab_punktow_H[j][i][1] = (tab_1[i][1] + tab_2[i][1]) * data.Conductivity * det[j];
				tab_punktow_H[j][i][2] = (tab_1[i][2] + tab_2[i][2]) * data.Conductivity * det[j];
				tab_punktow_H[j][i][3] = (tab_1[i][3] + tab_2[i][3]) * data.Conductivity * det[j];

				//obliczenie macierzy C w poszczególnych punktach całkowania
				tab_punktow_C[j][i][0] = tab_N[j][i] * tab_N[j][0] * data.Density * data.SpecificHeat * det[j];
				tab_punktow_C[j][i][1] = tab_N[j][i] * tab_N[j][1] * data.Density * data.SpecificHeat * det[j];
				tab_punktow_C[j][i][2] = tab_N[j][i] * tab_N[j][2] * data.Density * data.SpecificHeat * det[j];
				tab_punktow_C[j][i][3] = tab_N[j][i] * tab_N[j][3] * data.Density * data.SpecificHeat * det[j];
			}
		}

		double macierzH[4][4] = {};
		double macierzC[4][4] = {};

		//zsumowanie macierzy H i C z poszczególnych punktów całkowania wymnożonych przez wagi
		for (int k = 0; k < 4; k++) {
			int liczba_punktów = 0;
			for (int i = 0; i < waga.size(); i++) {
				for (int j = 0; j < waga.size(); j++) {
					macierzH[k][0] += tab_punktow_H[liczba_punktów][k][0] * waga[j] * waga[i];
					macierzH[k][1] += tab_punktow_H[liczba_punktów][k][1] * waga[j] * waga[i];
					macierzH[k][2] += tab_punktow_H[liczba_punktów][k][2] * waga[j] * waga[i];
					macierzH[k][3] += tab_punktow_H[liczba_punktów][k][3] * waga[j] * waga[i];

					macierzC[k][0] += tab_punktow_C[liczba_punktów][k][0] * waga[j] * waga[i];
					macierzC[k][1] += tab_punktow_C[liczba_punktów][k][1] * waga[j] * waga[i];
					macierzC[k][2] += tab_punktow_C[liczba_punktów][k][2] * waga[j] * waga[i];
					macierzC[k][3] += tab_punktow_C[liczba_punktów][k][3] * waga[j] * waga[i];

					liczba_punktów++;
				}
			}
		}

		// MACIERZ HBC I WEKTOR P

		double det_BC[4]; //wyznacznik jacobianu
		double tab_pomocniczna[4][4][4] = {};
		double tab_powierzchni[4][4][4] = {};
		double wektorP[4] = {}; //wektor obciążeń

		//sprawdzenie flag na węzłach ścian
		for (int i = 0; i < 4; i++) {
			if (i != 3 && grid.ND[element.ID[i]].BC == 1 && grid.ND[element.ID[i + 1]].BC == 1) {
				det_BC[i] = sqrt(pow((ID[0][i + 1] - ID[0][i]), 2) + pow((ID[1][i + 1] - ID[1][i]), 2)) / 2;
			}
			else if (i == 3 && grid.ND[element.ID[0]].BC == 1 && grid.ND[element.ID[i]].BC == 1) {
				det_BC[i] = sqrt(pow((ID[0][0] - ID[0][i]), 2) + pow((ID[1][0] - ID[1][i]), 2)) / 2;
			}
			else {
				det_BC[i] = 0;
			}
			for (int k = 0; k < sqrt(liczbaPunktow); k++) {
				for (int j = 0; j < 4; j++) {
					tab_pomocniczna[k][j][0] = boki[i].N[k][j] * boki[i].N[k][0];
					tab_pomocniczna[k][j][1] = boki[i].N[k][j] * boki[i].N[k][1];
					tab_pomocniczna[k][j][2] = boki[i].N[k][j] * boki[i].N[k][2];
					tab_pomocniczna[k][j][3] = boki[i].N[k][j] * boki[i].N[k][3];

					tab_powierzchni[i][j][0] += tab_pomocniczna[k][j][0] * waga[k] * data.Alfa * det_BC[i];
					tab_powierzchni[i][j][1] += tab_pomocniczna[k][j][1] * waga[k] * data.Alfa * det_BC[i];
					tab_powierzchni[i][j][2] += tab_pomocniczna[k][j][2] * waga[k] * data.Alfa * det_BC[i];
					tab_powierzchni[i][j][3] += tab_pomocniczna[k][j][3] * waga[k] * data.Alfa * det_BC[i];

					//zsumowanie wartości wektora obciążeń
					wektorP[j] += boki[i].N[k][j] * data.Tot * data.Alfa * det_BC[i]  * waga[k];
				}
			}
		}

		double macierzHBC[4][4] = {}; //macierz HBC

		//zsumowanie macierzy HBC każdej ściany
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				macierzHBC[i][0] += tab_powierzchni[j][i][0];
				macierzHBC[i][1] += tab_powierzchni[j][i][1];
				macierzHBC[i][2] += tab_powierzchni[j][i][2];
				macierzHBC[i][3] += tab_powierzchni[j][i][3];
			}
		}
		
		//dodanie do macierzy H macierzy HBC
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				macierzH[i][j] += macierzHBC[i][j];
			}
		}
		/*
		cout << "\nMacierz H:" << endl;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++)
				cout << macierzH[j][k] << " ";
			cout << endl;
		}

		cout << "\nMacierz HBC:" << endl;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++)
				cout << macierzHBC[j][k] << " ";
			cout << endl;
		}
		*/
		// AGREGACJA MACIERZY H, C I WEKTORA P

		for (int j = 0; j < 4; j++) {
			for (int i = 0; i < 4; i++) {
				macierz_H_Globalna[element.ID[j]][element.ID[i]] += macierzH[j][i];
				macierz_C_Globalna[element.ID[j]][element.ID[i]] += macierzC[j][i];
			}
			wektor_P_Globalny[element.ID[j]] += wektorP[j];
		}
	}

	//wypisywanie macierzy H i C po agregacjach
	/*
	cout << "\nMacierz [H]:" << endl;
	for (int j = 0; j < data.NodesNumber; j++) {
		for (int i = 0; i < data.NodesNumber; i++)
			cout << macierz_H_Globalna[j][i] << " ";
		cout << endl;
	}

	cout << "\nMacierz [C]:" << endl;
	for (int j = 0; j < data.NodesNumber; j++) {
		for (int i = 0; i < data.NodesNumber; i++)
			cout << macierz_C_Globalna[j][i] << " ";
		cout << endl;
	}
	*/
	//cout << "\nMacierz [H] = [H]+[C]/dT:" << endl;
	for (int j = 0; j < data.NodesNumber; j++) {
		for (int i = 0; i < data.NodesNumber; i++) {
			macierz_H_Globalna[j][i] += (macierz_C_Globalna[j][i] / data.SimulationStepTime);
			//cout << macierz_H_Globalna[j][i] << " ";
		}
		//cout << endl;
	}

	//cout << "\nWektor {P} = {P}+{[C]/dT}*{T0}" << endl;

	// ROZWIĄZANIE UKŁADU RÓWNAŃ

	for (int k = 0; (k+1)*data.SimulationStepTime<=data.SimulationTime; k++) {

		for (int j = 0; j < data.NodesNumber; j++) {
			wektor_P_Globalny_kopia[j] = wektor_P_Globalny[j];
		}

		for (int j = 0; j < data.NodesNumber; j++) {
			for (int i = 0; i < data.NodesNumber; i++) {
				wektor_P_Globalny_kopia[j] += (macierz_C_Globalna[j][i] / data.SimulationStepTime) * temp[i];
			}
			//if (k < 1)
				//cout << wektor_P_Globalny_kopia[j] << " ";
		}
		//if (k < 1)
			//cout << endl;

		//cout << endl;

		for (int j = 0; j < data.NodesNumber; j++) {
			for (int i = 0; i <= data.NodesNumber; i++) {
				if (i < data.NodesNumber)
					macierz_Gaussa[j][i] = macierz_H_Globalna[j][i];
				else
					macierz_Gaussa[j][i] = wektor_P_Globalny_kopia[j];
			}
		}

		if (gauss(data.NodesNumber, macierz_Gaussa, temp))
		{
			//for (int i = 0; i < data.NodesNumber; i++)
				//cout << "t" << i + 1 << " = " << temp[i] << endl;

			MinTemp[k] = temp[0];

			for (int j = 0; j < data.NodesNumber; j++) {
				if (MinTemp[k] > temp[j])
					MinTemp[k] = temp[j];
				if (MaxTemp[k] < temp[j])
					MaxTemp[k] = temp[j];
			}

			cout << "Czas: " << (k + 1) * data.SimulationStepTime << " MinTemp: " << MinTemp[k] << " MaxTemp: " << MaxTemp[k] << endl;
		}
		else
			cout << "Dzielenie przez zero!\n";
	}
	cout << endl;

	// CZYSZCZENIE PAMIĘCI

	for (int j = 0; j < data.NodesNumber; j++) {
		delete[] macierz_H_Globalna[j];
		delete[] macierz_C_Globalna[j];
		delete[] macierz_Gaussa[j];
	}
	delete[] macierz_H_Globalna;
	delete[] macierz_C_Globalna;
	delete[] macierz_Gaussa;
	delete[] wektor_P_Globalny;
	delete[] wektor_P_Globalny_kopia;
	delete[] temp;
	delete[] MinTemp;
	delete[] MaxTemp;

	system("PAUSE");
	return(0);
}

bool gauss(int n, double** AB, double* X)
{

	int i, j, k;
	double m, s;

	// eliminacja współczynników
	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych
	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}

	return true;
}