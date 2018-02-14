/*
 * Host.h
 *
 *  Created on: Apr 21, 2011
 *      Author: paola
 */

#ifndef HOST_H_
#define HOST_H_
#include "Genes.h"
#define TWO 2

const double YEAR = 31207680.0;
const double MONTH = 2600640.0;
const double WEEK = 604800.0;

class Virus
{
public:
	Virus(){};
	//	Virus(vector<uint32_t>& _nonself,Contain& cont,int _size);
	virtual ~Virus(){};
	void SetViralParameters(uint32_t _ID,vector<uint32_t>& ns,Contain& cont,double _viralLoad, double _lifeTimeVirus);
	Virus& Copy(Virus& rhsVirus);//works
	double GetLifeTimeVirus()const{return lifeTimeVirus;}
	double GetViralLoad()const{return viralLoad;}
	uint32_t GiveVirusID(){return virusID;};
	void SetVirusID(uint32_t num);
	void Fillpool(vector<uint32_t>& ns,Contain& cont);
	vector<uint32_t> GivePep(Contain& cont);
	void SavePep(vector<uint32_t> pool,Contain& cont);
	void SaveBackupVirus(fstream& backupFile);//works
	string RestoreVirus(stringstream& svline);//works

	void SaveParametersVirus(fstream& outfile);

protected:
	double viralLoad; //as the increase of the intrinsic death rate
	double lifeTimeVirus; // time that a virus can live in one's organism
	uint32_t virusID; 			//This virus ID


};


class Host {
public:
	Host(){}; //for creating memory
	Host(int loci_kir, int loci_mhc , double specificity, double _mutationRate, bool _tuning,MHCGenePool& mhcPool, bool hla,
			GenePool& kirPool, int numberOfExtraKirs, Contain& cont); //for initialization of the population

	//This is where all my troubles come from
	Host(int loci_kir,int loci_mhc, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist, GenePool& kirFromTheGenePool,
			vector<Gene>& kirGenesMother, vector<Gene>& kirGenesFather, double specificity, double _mutationRate,
			bool _tuning,int numberOfExtraKirs, Contain& cont);


	//Host(int loci_kir, Gene* mhcGenesParent, int mhcGenePool, Gene* kirGenesMother, Gene* kirGenesFather, int specificity, double _mutationRate, bool _tuning); //for birth function
	virtual ~Host(){};
	void SetDead(){dead = true;}
	bool IsDead()const{return dead;}
//	bool IsHostToBeTuned()const{return tuning;}

	void InitializeHostParameters(double mutationRate, bool _tuning, int kirloci,int mhcloci);
	void SetHostParameters(bool t, double mut, int inftyp, double inftim, double viraldeathm ,double clrtim);
	void EducateKIRs(Contain& cont);
	void EducateKIRs(double specificity,Contain& cont);//works
	void ExpressKIRs(int numberOfExpressedKirs);
	double GetAge()const{return age;}
	void SetAge(double number){age = number;}
	int CountFunctionalKIRs();//works
	int CountExpressedKIRs();
	Host& Copy(Host& rightHandSideHost);//works

	double GetIntrinsicDeathRate();
	double GetAgeDependentBirthRate();

	bool IsSusceptible();
	bool IsImmune();
	bool IsChronicInfected();
	bool IsAcuteInfected();
	bool IsIncubating();

	void InfectWith(Virus& nastyVirus, double simulationTime, double spec,Contain& cont);//works
	void ResetInfection(double simulationTime);
	void ClearInfection(double simulationTime,Contain& cont);
	double ProbabilityOfClearingTheInfection(int number_of_protective_receptors);

	void SetInfectionType(double simulationTime, double lifeTimeVirus);
	double GetInfectionTime(){return infectionTime;};
	double GetClearanceTime(){return clearanceTime;};
	int GetInfectionType(){return infectionType;};
	double GetMutationRate(){return mutationRateHost;};
	double GetViralDeathRate(){return viralDeathRate;}

	void UpdateParameters(double timeStep, double simulationTime,double lifeTimeVirus);//works
	void SaveGenes(fstream& outfile,Contain& cont);//works
	void SaveGenesID (fstream& outfile);
	void SaveParameters(fstream& outfile);//works
	void SaveAgeDyingHost(fstream& outfile);
	vector<Gene> GetmhcGenes(){return mhcGenes;};
	vector<Gene> GetkirGenes(){return kirGenes;};

	void SaveBackupHost(fstream& file);//works
	void RestoreHost(const string& sline);//works

	//Gene mhcGenes[LOCI_MHC][TWO];
	//Gene kirGenes[LOCI_KIR][TWO];

	enum state{susceptible, incubating, acute, chronic, immune};

	Virus pathogen;
protected:

	vector<Gene> mhcGenes;
	vector<Gene> kirGenes;
	int LOCI_KIR;
	int LOCI_MHC;
	double age;
	bool dead;
	bool tuning;
//	double tuning_rate;
	double mutationRateHost;

	double intrinsicDeathRate;
	double ageDependentBirthrate;
	double viralDeathRate;
    state infectionType;

    double infectionTime;
    double clearanceTime;
    int immunityTime;
    double ageInfection;
    double ageClearance;
    double PepSpecificity;
    //int MHCsearchinPool;
    //int functionalKIRs;

};

#endif /* HOST_H_ */
