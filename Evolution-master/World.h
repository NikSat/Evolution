/*
 * World.h
 *
 *  Created on: Apr 26, 2011
 *      Author: paola
 */

#ifndef WORLD_H_
#define WORLD_H_
#include "Host.h"//Host.h includes Virus.h and Gene.h which includes Bitstrings.h and Mathfunctions.h

class World {
public:
	World(); //Default constructor
	virtual ~World(){}; //destructor that so far does nothing at all!

	void LoadParameterFile(const string& fileName); //works
	bool Initialize();

	bool Birth(int index, Host& host,Contain& cont); 	// Birth() function takes a host at index (index) and picks a random host at a random index and create a "child" that is copied to (host).
										// It returns true or false according if the child should be added to the population or not
	bool Death (int index);
	void Infect(int index);

	void ShuffleHosts();
	void IntroduceVirus(Virus& Vir);
	void RemoveDeadHosts();
	void Updatemap(int id,string pepseq);
	void Initializepeppools();			//fills self, nonself pools;
	//void LoadPools();

	void CreateVirus();	//Creates one Virus
	//void LoadVirusPool(); //Restoring the virus

	void CheckAssociations();
	void Simulate();

	void SaveGenes(); //works
	void SavePopulationSize(int babies,int dead_people);//works
	void SaveParameters();//works
	bool SaveAgeDyingHosts(double lastOutfileTime, double lastStopOutfileTime);

	void SaveBackupFile();//works
	//void SaveVirusPool();
	//void SavePools();
	void SaveMap();

	void LoadBackupFile(const string& backupName);//works
	void LoadMap();
	void Loadmatrix(string& mfile);

protected:
	double birthRate;
	double deathRate;
	double infectionRate;
	double escapeRate;

	double timeStep;
	double timeEnd;
	double simulationTime;
	double timeIntroducingInfection;
	double outfileRate;
	double backupRate;
	double timeRecording;
	int maxHostPop;

	int contactRate;
	double timeMHCDownregulation;
	double timeDecoy;

	double deltaVirus;
	double timeInfection;

	double transmissionRateAcute;
	double transmissionRateChronic;

	vector<int> shuffledHosts;
	vector<int> virtualHosts;

	int acute_infected;
	int chronic_infected;
	int immune;
	int virnum;

	bool isFileOpen;
	fstream populationSize;
	fstream genesFile;
	fstream parameterFile;
	fstream backupFile;
	fstream notBornChildren;
	fstream dyingHosts;

	Contain container;
	vector<Host> hosts; //vector of hosts
	MHCGenePool MHCPool; //vector of integers representing MHC pools
	GenePool KIRPool; //vector of integers representing KIR pools
	Virus Wildtype; //The only virus in the simulation

	vector<uint32_t> self;		//common for all humans
	vector<uint32_t> nonself;	//pool from which the virus fills its pool
	double mutationRate;
	bool education;
	int expressionExtraKIRs;
	bool HLA_C;
	double KIRspecificity;
	int KIRLoci;
	int MHCLoci;
	unsigned int initHostPop;
	bool printoption;
	int selfpoolsize;
	int nonselfpoolsize;
	double specificityPEP;
	int peplength;


};

#endif /* WORLD_H_ */
