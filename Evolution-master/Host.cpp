/*
 * Host.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: paola
 */

#include "Host.h"
/* FUNCTIONS OF CLASS VIRUSSPECIES */


/* FUNTIONS OF CLASS VIRUS
 * 	Constructs the virus*/

void Virus::SetVirusID(uint32_t num)
{
	virusID=num;
}

void Virus::SavePep(vector<uint32_t> pool,Contain& cont)
{
	cont.UpdateViralPool(virusID,pool);
}

/* Fills the Virus individual non self pool from the greater nonself pool, ****ONE VIRUS*** */

void Virus::Fillpool(vector<uint32_t>& ns,Contain& cont)
{
	vector<uint32_t> nspeps=ns;
	cont.UpdateViralPool(virusID,nspeps);
}

void Virus::SetViralParameters(uint32_t _ID,vector<uint32_t>& ns,Contain& cont,double _viralLoad, double _lifeTimeVirus)
{
	virusID=_ID;
	this->Fillpool(ns,cont);
	viralLoad = _viralLoad;
	lifeTimeVirus = _lifeTimeVirus;

}

/*This function performs a deep copy*/
Virus & Virus::Copy(Virus& rhsVirus)
{
	// checking if it is a self-assignment
	if(this == &rhsVirus)
		return *this;
	//copying member variables
	this->virusID=rhsVirus.virusID;
	this->lifeTimeVirus = rhsVirus.lifeTimeVirus;
	this->viralLoad = rhsVirus.viralLoad;
	return *this;

}

void Virus::SaveBackupVirus(fstream& file)
{
	file << virusID << lifeTimeVirus << "\t" << viralLoad <<"\t";
}

string Virus::RestoreVirus(stringstream& svline)
{
	svline >> virusID;
	svline >> lifeTimeVirus;
	svline >> viralLoad;
	string vstring = svline.str();
	return vstring;
}

void Virus :: SaveParametersVirus(fstream& outfile)
{
	outfile << virusID << "\t" << viralLoad <<"\t" ;
}


vector<uint32_t> Virus::GivePep(Contain& cont)
{
	vector<uint32_t> pool=cont.GetViralPool(virusID);
	return pool;
}


/*FUNCTIONS OF CLASS HOST
 *Constructs a host for the initialization of the population: it fills the MHC genes with a randomly picked allele from the population
 * and creates KIRs that match their own MHC according to the specificity */
Host::Host(int loci_kir, int loci_mhc, double specificity, double _mutationRate, bool _tuning, MHCGenePool& mhcPool, bool hla,GenePool& kirPool,int numberOfExtraKirs, Contain& cont )
{

	InitializeHostParameters(_mutationRate,_tuning, loci_kir, loci_mhc);
	//fill the mhc Genes
	for(int i = 0; i <LOCI_MHC*TWO; i ++)
		{
		Gene firstGene;
		int mhc1 = mhcPool.RandomlyPickGene(hla);
		firstGene.SetGeneID(mhc1);
		mhcGenes.push_back(firstGene);
		}


	//create kirs that recognize their own mhcs / or mhc from the population (not necessarily their own)
	for(int i=0; i<LOCI_KIR; i++)
		{
//		int kir_id = kirPool.RandomlyPickGene(false); this is if I want to start with more KIRs
		uint32_t kir_id = kirPool.GetGenes().at(10*i+6); //i will start with 10 KIR!
		Gene kir;
		kir.SetGeneID(kir_id);
		kir.SetGeneFunctionality(false);
		kir.SetGeneExpression(false);
		kirGenes.push_back(kir);
		}
	//copy the 5 first kirs into the other haplotype
	for(int i=0; i<LOCI_KIR; i++)
		{
		Gene kir = kirGenes.at(i);
		kirGenes.push_back(kir);
		}

	if(tuning == true)
	EducateKIRs(cont);
	ExpressKIRs(numberOfExtraKirs);

	age = RandomNumber(1,70); //population initialized with a random age between 1 and 70

}

/*Constructs a baby host out of two parents*/
Host::Host(int loci_kir,int loci_mhc, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist, GenePool& kirPool, vector<Gene>& kirGenesMother, vector<Gene>& kirGenesFather,  double pepspecificity,double _mutationRate, bool _tuning, int numberOfExtraKirs, Contain& cont)
{	//to create a NEW host: the haplotypes of KIR of BOTH parents are needed. Besides one MHC haplotype of one parent plus one of the pool

	//int hhhhaaaap = 0;
	InitializeHostParameters(_mutationRate, _tuning, loci_kir, loci_mhc);
	//uint32_t mhcFromTheGenePool = mhcPool.RandomlyPickGene(dist);
//	int matchingGenes = 0;

	int KIR_init_mum = 0;
	int KIR_end_mum = 0;
	int MHC_init_mum = 0;
	int MHC_end_mum = 0;

	int KIR_init_dad = 0;
	int KIR_end_dad = 0;
	//pick haplotype 1/0 of each parent!
	int hap_mum=(RandomNumberDouble()<0.5);
	int hap_dad= (RandomNumberDouble()<0.5);
	//int k= 0;
	if(hap_mum)
		{
		KIR_init_mum = LOCI_KIR;
		KIR_end_mum = LOCI_KIR*TWO;
		}
	else
		{
		KIR_init_mum = 0;
 		KIR_end_mum = LOCI_KIR;
		}

	if(hap_dad)
		{
		KIR_init_dad = LOCI_KIR;
		KIR_end_dad = LOCI_KIR*TWO;
		}
	else
		{
		KIR_init_dad = 0;
		KIR_end_dad = LOCI_KIR;
		}

	int hap_mhc = (RandomNumberDouble()<0.5);
	if(hap_mhc)
		{
		MHC_init_mum = LOCI_MHC;
		MHC_end_mum = LOCI_MHC*TWO;
		}

	else
		{
		MHC_init_mum = 0;
		MHC_end_mum = LOCI_MHC;
		}


	//copy the KIR haplotype into the new host (mutation occurs!)
	int hap_child = (RandomNumberDouble()<0.5);
	if(hap_child)
		{
		for(int i=KIR_init_dad; i<KIR_end_dad; i++)
			{
			Gene kir_hap2;
			kir_hap2.SetGeneID(kirGenesFather.at(i).GetGeneID());
			kir_hap2.SetGeneFunctionality(false);
			kir_hap2.SetGeneExpression(false);
//			kir_hap2.Mutate(mutationRateHost); //this is point mutation... i will try to use another trick now.
											   //which is why I will just steal another KIR from the pool as a mutation operation

			if(RandomNumberDouble() < mutationRateHost)
				{
				kir_hap2.Mutate(kirPool.RandomlyPickGene(false));
//				kir_hap2.CheckCross(mhcPool,pool,specificity,pepspecificity,cont);
				}
			kirGenes.push_back(kir_hap2);

			}
		for(int j=KIR_init_mum; j<KIR_end_mum; j++)
			{
			Gene kir_hap1;
			kir_hap1.SetGeneID(kirGenesMother.at(j).GetGeneID());
			kir_hap1.SetGeneFunctionality(false);
			kir_hap1.SetGeneExpression(false);
			if(RandomNumberDouble() < mutationRateHost)
				{
				kir_hap1.Mutate(kirPool.RandomlyPickGene(false));
//				kir_hap1.CheckCross(mhcPool,pool,specificity,pepspecificity,cont);
				}

			kirGenes.push_back(kir_hap1);

			}
		//copy the MHC haplotype into the new host
		for(int m = MHC_init_mum; m < MHC_end_mum; m++)
			{
			Gene mhc1;
			mhc1.Copy(mhcGenesParent.at(m));
			mhcGenes.push_back(mhc1);
			//cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
			}

		//generate new mhc genes (from the pool) for the new born
		for(int mm = 0; mm<LOCI_MHC; mm++)
			{
			Gene mhc2;
			int mhcFromTheGenePool = mhcPool.RandomlyPickGene(dist);
			mhc2.SetGeneID(mhcFromTheGenePool);
			mhcGenes.push_back(mhc2);
			}
		}
	else
		{
		for(int i=KIR_init_mum; i<KIR_end_mum; i++)
			{
			Gene kir_hap1;
			kir_hap1.SetGeneID(kirGenesMother.at(i).GetGeneID());
			kir_hap1.SetGeneFunctionality(false);
			kir_hap1.SetGeneExpression(false);
			if(RandomNumberDouble() < mutationRateHost)
				{
				kir_hap1.Mutate(kirPool.RandomlyPickGene(false));
//				kir_hap1.CheckCross(mhcPool,pool,specificity,pepspecificity,cont);
				}
			kirGenes.push_back(kir_hap1);

			}
		for(int j=KIR_init_dad; j<KIR_end_dad; j++)
			{
			Gene kir_hap2;
			kir_hap2.SetGeneID(kirGenesFather.at(j).GetGeneID());
			kir_hap2.SetGeneFunctionality(false);
			kir_hap2.SetGeneExpression(false);
			if(RandomNumberDouble() < mutationRateHost)
				{
				kir_hap2.Mutate(kirPool.RandomlyPickGene(false));
//				kir_hap2.CheckCross(mhcPool,pool,specificity,pepspecificity,cont);
				}
			kirGenes.push_back(kir_hap2);

			}


		//generate new mhc genes (from the pool) for the new born
		for(int mm = 0; mm<LOCI_MHC; mm++)
			{
			Gene mhc2;
			int mhcFromTheGenePool = mhcPool.RandomlyPickGene(dist);
			mhc2.SetGeneID(mhcFromTheGenePool);
			mhcGenes.push_back(mhc2);
			}
			//copy the MHC haplotype into the new host
		for(int m = MHC_init_mum; m < MHC_end_mum; m++)
			{
			Gene mhc1;
			mhc1.Copy(mhcGenesParent.at(m));
			mhcGenes.push_back(mhc1);
			//cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
			}

	}



	if(tuning == true)
		{
		EducateKIRs(cont);
		}
	ExpressKIRs(numberOfExtraKirs);
	age = 1.0; //newborns are given the age of 15
}

void Host::InitializeHostParameters(double mutationRate, bool _tuning, int loci_kir,int loci_mhc)
{
	dead = false;
	mutationRateHost = mutationRate;
	tuning = _tuning;
	LOCI_KIR = loci_kir;
	LOCI_MHC = loci_mhc;
	infectionType = susceptible;
	infectionTime = 0.0;
	immunityTime = 0.0;
	viralDeathRate = 0.0;
	clearanceTime = 0.0;
	ageInfection = 0.0;
	ageClearance = 0.0;
	age = 0.0;
	//MHCsearchinPool = 0;

}

/*This functions tunes the KIR repertoire according to the self MHC repertoire*/

void Host :: EducateKIRs(Contain& cont)
{
	for(unsigned int i=0; i<kirGenes.size(); i++)  	//For how many MHCs each KIR is licensed is already calculated and stored in the association map, education is just retrieval of that data
		{
		int count=0;
		std::pair <std::multimap<uint32_t,uint32_t>::iterator, std::multimap<uint32_t,uint32_t>::iterator> ret;  //For each KIR, iterate through the association map and get the MHCs this KIR is licensed for
		ret = cont.GiveAssoc().equal_range(kirGenes.at(i).GetGeneID());
		for (std::multimap<uint32_t,uint32_t>::iterator it=ret.first; it!=ret.second; ++it)
		   	{
		   	for (unsigned int i=0; i<mhcGenes.size(); i++)    //For each MHC, check whether if it is one of the collected MHCs the KIR is licensed for
		   		{
		   		if (it->second==mhcGenes.at(i).GetGeneID()) count++;
	    		}
	    	}
	    	if (count>0)				//Set KIRs as functional and expressed (expressed does not make a difference, it is a remnant of one of Paola's projects)
	    		{
	    		kirGenes.at(i).SetGeneFunctionality(true);
	    		kirGenes.at(i).SetGeneExpression(true);
	    		}

		}
}



void Host :: ExpressKIRs(int numberOfExtraKirs)
{
	int counter = 0;
	if(tuning == false)
		{
		for(unsigned int i=0; i<kirGenes.size(); i++)
			kirGenes.at(i).SetGeneExpression(true);
		// if there is no education, all KIRs should be expressed!
		}
	else //otherwise express more KIRs besides those which are already functional
		{
		for(unsigned int i=0; i<kirGenes.size(); i++)
			{
			if(kirGenes.at(i).IsFunctional()) //KIRs that are functional are already expressed (see EducateKIRs!)
				continue;
			else
				{
				if(counter<numberOfExtraKirs)
					{
					kirGenes.at(i).SetGeneExpression(true);
					counter++;
					}
				else
					return;
				}
			}
		}
}

/* This function counts the UNIQUE KIRs within one host that are functional, i.e. that are able to recognize their own MHC*/
int Host::CountFunctionalKIRs()
{
	int kirsThatAreFunctional = 0;

	for(unsigned int i= 0; i<kirGenes.size(); i++)
	{
		if(kirGenes.at(i).IsFunctional())
		{
			int flagNotUnique = 1;
			for(unsigned int j= 0; j<i; j++)
			{
				if(kirGenes.at(i).GetGeneID()==kirGenes.at(j).GetGeneID())
					flagNotUnique++;
			}
			if (flagNotUnique == 1)
				kirsThatAreFunctional ++;
		}
	}

	return kirsThatAreFunctional;
}

/* This function counts ALL KIRs within one host that are expressed*/
int Host::CountExpressedKIRs()
{
	int kirsThatAreExpressed = 0;

	for(unsigned int i= 0; i<kirGenes.size(); i++)
	{
		if(kirGenes.at(i).IsExpressed())
			{
			kirsThatAreExpressed ++;
			}
	}

	return kirsThatAreExpressed;
}
/*This function performs a deep copy*/
Host& Host:: Copy(Host& rightHandSideHost)
{
	// checking if it is a self-assignment
	if(this == &rightHandSideHost)
		return *this;

	//copying member variables
	this->dead = rightHandSideHost.dead;
	this->tuning = rightHandSideHost.tuning;
	this->age = rightHandSideHost.age;
	this->LOCI_KIR = rightHandSideHost.LOCI_KIR;
	this->infectionType = rightHandSideHost.infectionType;
	this->infectionTime = rightHandSideHost.infectionTime;
	this->immunityTime = rightHandSideHost.immunityTime;
	this->mutationRateHost = rightHandSideHost.mutationRateHost;
	this->viralDeathRate = rightHandSideHost.viralDeathRate;
	this->clearanceTime = rightHandSideHost.clearanceTime;
	this->ageInfection = rightHandSideHost.ageInfection;
	this->ageClearance = rightHandSideHost.ageClearance;

	//copy genes
	for(unsigned int i=0; i<rightHandSideHost.mhcGenes.size(); i++)
		{
		this->mhcGenes.push_back(rightHandSideHost.mhcGenes.at(i));
		//cout << "mhc genes: \t"<<mhcGenes.at(i).GetGeneID() << "\t" << rightHandSideHost.mhcGenes.at(i).GetGeneID() <<endl;
		}
	for(unsigned int i=0; i<rightHandSideHost.kirGenes.size(); i++)
		{
		this->kirGenes.push_back(rightHandSideHost.kirGenes.at(i));
		}

	return *this; // returns self-reference so cascaded assignment works
}

double Host :: GetIntrinsicDeathRate()
{
	intrinsicDeathRate = exp(0.1*age-10.5)+ exp(-0.4*age-8)+ viralDeathRate;
	return intrinsicDeathRate;
}

double Host :: GetAgeDependentBirthRate()
{
	ageDependentBirthrate = -1/(1+ exp(age-20)) + 1/(1+exp(age-45));
	return ageDependentBirthrate;
}

/*This functions return the infection type of the host*/
bool Host :: IsSusceptible()
{
	if (infectionType == susceptible)
		return true;
	else
		return false;
	//return (infection == susceptible);
}
bool Host :: IsImmune()
{
	if (infectionType == immune)
		return true;
	else
		return false;
}
bool Host:: IsIncubating()
{
	if(infectionType == incubating)
		return true;
	else
		return false;
}
bool Host :: IsAcuteInfected()
{
	if (infectionType == acute)
		return true;
	else
		return false;
}
bool Host :: IsChronicInfected()
{
	if(infectionType == chronic)
		return true;
	else
		return false;
}

/*This function sets the pathogen parameters as infectious*/
void Host :: InfectWith(Virus& nastyVirus, double simulationTime, double spec,Contain& cont)
{
	if(infectionType == susceptible||infectionType == incubating)
		{
		pathogen.Copy(nastyVirus);
		infectionTime = simulationTime;
		ageInfection = age;
		immunityTime = 10;
		}
}

/*This function changes the infection type of the host according to the time of infection*/
void Host::SetInfectionType(double simulationTime, double lifeTimeVirus)
{
	//infection has started but has not been cleared yet:
	if(infectionTime > 0.0 && clearanceTime == 0.0)
		{
		//check first if the host is in the infection period (i.e. between 1-4 weeks)
		if((simulationTime-infectionTime)<=1.0*WEEK)
			infectionType = incubating;

		if((simulationTime-infectionTime)>1.0*WEEK && (simulationTime-infectionTime)<1.0*WEEK+4.0*lifeTimeVirus)
			{
			viralDeathRate = pathogen.GetViralLoad();
			infectionType = acute;
			}
		if((simulationTime -infectionTime) > 1.0*WEEK+ 4.0*lifeTimeVirus)
			{
			viralDeathRate = 0.6*pathogen.GetViralLoad();
			infectionType = chronic;
			}
		}

	if(infectionTime == 0.0 && clearanceTime > 0.0)
		{
		if(simulationTime - clearanceTime<ImmunityTime(immunityTime,0.5)*YEAR)
			{
			infectionType = immune;
			}
		else
			{
			infectionType = susceptible;
			clearanceTime = 0.0;
			}
		}

	if(infectionTime == 0.0 && clearanceTime == 0.0)
		infectionType = susceptible;
}

void Host :: UpdateParameters(double timeStep, double simulationTime, double lifeTimeVirus)
{
	age += (timeStep/YEAR);
	SetInfectionType(simulationTime, lifeTimeVirus);

}


/*This function gives the probability of clearing the infection for linear contribution of each KIR to clearance */

double Host :: ProbabilityOfClearingTheInfection(int number_of_protective_receptors)
{
	double p = 0.25;
	if(number_of_protective_receptors) //if the number of protective receptors is > than 0
		p = ((number_of_protective_receptors *0.60)/CountFunctionalKIRs())+0.25;
	return p;
}


/*
//Linear Contribution Rule: The clearance probability is variable, each licensed protective KIR contributes linearly to it to a maximum of 0.85 (if all licensed KIRs are protective, then the clearance probability is 0.85)

///************NOT USED IN THIS SIMULATION*****************

void Host::ClearInfection(double simulationTime,Contain& cont)
{
	// Rule 1: Linear Contribution Rule. Each KIR contributes linearly to the infection clearance, starting from 0.25 when there are no protective KIRs
	//to a maximum of 0.85
	int kirsnotrecognizingallpmhcs=0;
	int functionalkirs=0;



	for (unsigned int k=0; k<kirGenes.size(); k++)
	{
		if (kirGenes.at(k).IsFunctional())
		{

			functionalkirs++;
			Gene KIR=kirGenes.at(k);
			uint32_t KIRID=KIR.GetGeneID();


//			Sequence KIRSeq(KIRID,cont);

//			cout << "KIRID\t"<<KIRID <<"\n";

//			cout << "KIRSequence\t"<<KIRSeq.GetPep() <<"\n";

			int sum=0;

			for (unsigned int l=0;l<mhcGenes.size();l++)
				{
				Gene MHC=mhcGenes.at(l);
				uint32_t MHCID=MHC.GetGeneID();



//				Sequence MHCSeq(MHCID,cont);

//  			cout <<"MHCID\t"<< MHCID <<"\n";

//				cout << "MHCSequence"<<MHCSeq.GetPep() <<"\n";




				int subsum=0;
				for (unsigned int m=0;m<pathogen.GivePep(cont).size();m++)
					{
					uint32_t VirusID=pathogen.GivePep(cont).at(m);

//					Sequence VirusSeq(VirusID,cont);

//					cout << "VirusID\t"<<VirusID <<"\n";

//					cout << "VirusSequence\t"<<VirusSeq.GetPep() <<"\n";

//					if (cont.GetbindValue(MHCID,KIRID,VirusID)==1)
//					{
//					cout << "Binds?\tYes\n";
//					}
					int one=cont.GetbindValue(MHCID,KIRID,VirusID);
//					cout << one << "\t";
					subsum+=cont.GetbindValue(MHCID,KIRID,VirusID);
					}
				sum+=subsum;
//				cout << subsum << "\n";
				}
			if (sum==0)
				{
					kirsnotrecognizingallpmhcs++;
				}
//				cout << sum << "\t"<< kirsnotrecognizingallpmhcs <<"\n";
		}

	}



	// We now have the number of functional kirs and the number of KIRs that do not recognize all the pMHCs

	//If we have functional kirs:

	if(CountFunctionalKIRs())
	{
		if(RandomNumberDouble()<ProbabilityOfClearingTheInfection(kirsnotrecognizingallpmhcs))
		{
			ResetInfection(simulationTime);
			return;
		}
	}
//	cout <<"\n\n\n";
}
*/



//All or nothing Rule: ALL licensed KIRs should not recognize any viral peptides to have 0.85 clearance probability, if not the clearance is 0.25


void Host::ClearInfection(double simulationTime,Contain& cont)
{
	// Rule 2: All or nothing rule. If all the functional KIRs fail to recognize all the pMHC complexes then there is a 85% probability
	// of clearing the infection
	int kirsnotrecognizingallpmhcs=0;
	int functionalkirs=0;

//Again all the interactions are pre-calculated, and it is just a matter of retrieving the information from the 3D vector

	for (unsigned int k=0; k<kirGenes.size(); k++)		//Select all the KIRs
		{
		if (kirGenes.at(k).IsFunctional())			//Check if they are functional (licensed)
			{
			functionalkirs++;
			Gene KIR=kirGenes.at(k);
			uint32_t KIRID=KIR.GetGeneID();

//			Sequence KIRSeq(KIRID,cont);

//			cout << "KIRID\t"<<KIRID <<"\n";

//			cout << "KIRSequence\t"<<KIRSeq.GetPep() <<"\n";

			int sum=0;

			for (unsigned int l=0;l<mhcGenes.size();l++)		//For each MHC check how many viral pMHC complexes the KIR recognizes
				{
				Gene MHC=mhcGenes.at(l);
				uint32_t MHCID=MHC.GetGeneID();



//				Sequence MHCSeq(MHCID,cont);

//  			cout <<"MHCID\t"<< MHCID <<"\n";

//				cout << "MHCSequence"<<MHCSeq.GetPep() <<"\n";




				int subsum=0;													//Check for all viral peptides whether pMHC complexes are formed and recognized by the KIR
				for (unsigned int m=0;m<pathogen.GivePep(cont).size();m++)
					{
					uint32_t VirusID=pathogen.GivePep(cont).at(m);

//					Sequence VirusSeq(VirusID,cont);

//					cout << "VirusID\t"<<VirusID <<"\n";

//					cout << "VirusSequence\t"<<VirusSeq.GetPep() <<"\n";

//					if (cont.GetbindValue(MHCID,KIRID,VirusID)==1)
//					{
//					cout << "Binds?\tYes\n";
//					}
					int one=cont.GetbindValue(MHCID,KIRID,VirusID);
//					cout << one << "\t";
					subsum+=cont.GetbindValue(MHCID,KIRID,VirusID);
					}
				sum+=subsum;
//				cout << subsum << "\n";
				}
			if (sum==0)
				{
				kirsnotrecognizingallpmhcs++;
				}
//				cout << sum << "\t"<< kirsnotrecognizingallpmhcs <<"\n";
		}

	}



	double P=RandomNumberDouble();

	// We now have the number of functional kirs and the number of KIRs that do not recognize all the pMHCs

	//If we have functional kirs:

	if (functionalkirs>0)
		{

	//If *all* of them do not recognize the pMHC we have maximum protection 0.85

		if (kirsnotrecognizingallpmhcs==functionalkirs)
			{
			if (P<=0.85)
				{
				ResetInfection(simulationTime);
//				cout <<"Yes\n";
				}

	//in all other cases we have minimum protection 0.25
			}
		else
			{
			if (P<=0.25)
				{
				ResetInfection(simulationTime);
//				cout <<"No\n";
				}
			}
		}
	else
		{
		if (P<=0.25) ResetInfection(simulationTime);
		//	cout <<"No\n";
		}
//exit(0);
//	cout <<"\n\n\n";
}

/* Resets the infection  */

void Host::ResetInfection(double simulationTime)
{
	Virus deadVirus;
	deadVirus.SetVirusID(0);
	viralDeathRate = 0.0;
	pathogen.Copy(deadVirus);
	infectionTime = 0.0;
	clearanceTime = simulationTime;
	ageClearance = age;

}

//functions for SAVING PARAMETERS, BACKUP, etc
/*this function saves each locus -> to keep track of KIR and MHC diversity*/
void Host::SaveGenes(fstream& outfile,Contain& cont) // This function prints the Sequences
{
	for(unsigned int i=0; i<mhcGenes.size(); i++)
		{
		outfile << cont.getseq(mhcGenes.at(i).GetGeneID()) << "\t"<<"\t";
		}
	for(unsigned int i=0; i<kirGenes.size(); i++)
		{
		outfile << cont.getseq(kirGenes.at(i).GetGeneID()) << "\t"<<"\t";
		}

	outfile << CountFunctionalKIRs() << "\t"<<"\t"<< CountExpressedKIRs() <<"\n";
}

void Host::SaveGenesID(fstream& outfile) // This funtion prints the ID's
{
	for(unsigned int i=0; i<mhcGenes.size(); i++)
		{
		outfile << mhcGenes.at(i).GetGeneID()<< "\t";
		}
	for(unsigned int i=0; i<kirGenes.size(); i++)
		{
		outfile << kirGenes.at(i).GetGeneID() << "\t" ; //Prints ID
		}

	outfile << CountFunctionalKIRs() << "\t"<<"\t"<< CountExpressedKIRs() <<"\n";

}

void Host ::SaveParameters(fstream& outfile)
{
	outfile << "\t"<<age << "\t" << infectionTime/YEAR <<"\t" << infectionType << "\t"<< viralDeathRate<<"\t"<<ageInfection << "\t"<< ageClearance<<"\t";
	pathogen.SaveParametersVirus(outfile);
	outfile <<  "\n";
}


void Host :: SaveAgeDyingHost(fstream& outfile)
{
	outfile << age << "\t";

	for(unsigned int i=0; i<kirGenes.size(); i++)
		{
		kirGenes.at(i).SaveGenes(outfile);
		}
	pathogen.SaveParametersVirus(outfile);
	outfile << infectionType << "\t"<< ageInfection << "\t"<< ageClearance<<"\t";
}


void Host::SaveBackupHost(fstream& backupFile)
{
	backupFile <<LOCI_MHC<< LOCI_KIR << "\t"<< age << "\t"<< infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<<immunityTime << "\t"<< tuning << "\t"
	<< dead << "\t"<< mutationRateHost << "\t"<< viralDeathRate << "\t"<<ageInfection <<"\t"<<ageClearance <<"\t";

	pathogen.SaveBackupVirus(backupFile);

	for(unsigned int i=0; i<mhcGenes.size(); i++)
		{
		mhcGenes.at(i).SaveBackupGenes(backupFile);
		}
	for(unsigned int i=0; i<kirGenes.size(); i++)
		{
		kirGenes.at(i).SaveBackupGenes(backupFile);
		}
	backupFile << "\n";
}

void Host::RestoreHost(const string& sline)
{
	stringstream ssline(sline);
	ssline >> LOCI_MHC;
	ssline >> LOCI_KIR;
	ssline >>age;
	//Ouss ssline >> infectionType;
	int inf;
	ssline >> inf;
	//cout <<"inf:  "<<inf <<endl;
	switch(inf)
		{
		case 0:{infectionType = susceptible;}break;
		case 1:(infectionType = incubating);break;
		case 2:{infectionType = acute;}break;
		case 3:{infectionType = chronic;}break;
		case 4:{infectionType = immune;}break;
		}
	//cout << "infection Type: "<<infectionType <<endl;
	ssline >> infectionTime;
	ssline >> clearanceTime;
	ssline >> immunityTime;
	ssline >> tuning;
	ssline >> dead;
	ssline >> mutationRateHost;
	ssline >> viralDeathRate;
	ssline >> ageInfection;
	ssline >> ageClearance;
//	cout << LOCI_KIR<<"\t"<<age << "\t"<< infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<< tuning << "\t"<< dead << "\t"<< mutationRateHost << "\t"<< viralDeathRate << "\t"<<endl;
	string geneString = pathogen.RestoreVirus(ssline);
	ssline.str() = geneString;

	for(int i=0; i<LOCI_MHC*TWO; i++)
		{
		Gene mhc;
		geneString = mhc.RestoreGenes(ssline);
		mhcGenes.push_back(mhc);
		ssline.str() = geneString;
		}

	for(int i=0; i<LOCI_KIR*TWO; i++)
		{
		Gene kir;
		geneString = kir.RestoreGenes(ssline);
		kirGenes.push_back(kir);
//		cout <<"id:\t"<<kir.GetGeneID()<<" func:\t"<<kir.IsFunctional()<<" kirSize: "<< kirGenes.size()<<endl;
//		cout <<ssline.str()<<endl;
		ssline.str() = geneString;
		}
//	cout << "halloooooooooooooooooooooooooooo"<<endl;

}
