/*
 * World.cpp
 *
 *  Created on: Apr 26, 2011
 *      Author: paola
 */

#include "World.h"
#include "backup.h"

World::World(){}

void World :: LoadParameterFile(const string& fileName)
{
	double t_outfile;
	double t_backup;
	kaBackup paramFile(fileName,false);
	paramFile.Load("Maximal population size", &maxHostPop);
	paramFile.Load("Initial population size", &initHostPop);
	paramFile.Load("Time Step", &timeStep);
	paramFile.Load("End of Simulation", &timeEnd);
	paramFile.Load("Time introducing infection", &timeIntroducingInfection);
	paramFile.Load("Time outfile", &t_outfile);
	//paramFile.Load("Time recording", &timeRecording);
	paramFile.Load("Time backup", &t_backup);
	paramFile.Load("Mutation rate host", &mutationRate);
	paramFile.Load("Contacts per week", &contactRate);
	paramFile.Load("KIR Loci", &KIRLoci);
	paramFile.Load("MHC Loci", &MHCLoci);
	paramFile.Load("HLA-C alleles distribution", &HLA_C);
	paramFile.Load("MHC-KIR specificity", &KIRspecificity);
	paramFile.Load("Tuning", &education);
	paramFile.Load("Extra Number KIRs to be expressed", &expressionExtraKIRs);
	paramFile.Load("Viral load", &deltaVirus);
	paramFile.Load("Life time of the virus", &timeInfection);
	paramFile.Load("Transmission rate acute infection", &transmissionRateAcute);
	paramFile.Load("Transmission rate chronic infection", &transmissionRateChronic);
	paramFile.Load("Gene output", &printoption);
	paramFile.Load("Self Pool Size", &selfpoolsize);
	paramFile.Load("Non Self Pool Size", &nonselfpoolsize);
	paramFile.Load("Presentation Threshold", &specificityPEP);
	paramFile.Load("Presented Peptide Length", &peplength);


	const string populationFile("PopulationSize.txt");

	populationSize.open(populationFile.c_str(), ios::out);
	populationSize << "#time\tpopSize\tbabies\tdeath\tacute\tchronich\timmune\n";
	timeStep= timeStep*WEEK;
	simulationTime = 0.0;
	timeEnd = timeEnd*YEAR;
	//timeRecording = timeRecording * YEAR;
	timeIntroducingInfection = timeIntroducingInfection*YEAR;

	outfileRate = 1.0/(t_outfile*YEAR);
	backupRate = 1.0/(t_backup*YEAR);

	birthRate = 0.5*timeStep/YEAR; // every birth event will happen once every four years
	deathRate = timeStep/YEAR; // every death event will happen once every year
	infectionRate = 52.0*timeStep/YEAR; //every infection event will happen every week (pop concerts removed!)
	escapeRate = 0.0001*timeStep/YEAR;


	acute_infected = 0;
	chronic_infected = 0;
	immune = 0;

	cout << "Beginning simulation with: \n";
	cout << "End Simulation after " << timeEnd/YEAR << " years \n";
	cout << "Maximal population size " << maxHostPop << " individuals \n";
	cout << "Initial population size " << initHostPop << " individuals \n";
	cout << "Introducing virus after " << timeIntroducingInfection/YEAR << " years \n";
	cout << "Mutation rate host " << mutationRate << "\n";
	cout << "Contacts per week: " << contactRate << "\n";
	cout << "HLA-C alleles distribution " << HLA_C<<"\n";
	cout << "MHC-KIR specificity " << KIRspecificity<<"\n";
	cout << "Tuning "<< education << "\n";
	cout << "Expression extra KIRs "<< expressionExtraKIRs << "\n";
	cout << "Viral load: " << deltaVirus << "\n";
	cout << "Life time virus: " << timeInfection << "\n";
	cout << "Tranmission rate acute infection "<< transmissionRateAcute << "\n";

}

// initialize host population
bool World::Initialize()
{



	MHCPool.FillGenePoolFromContainer(container); //Creates 14 completely random molecules

	//MHCPool.FillGenePoolFromContainer(container,self,specificityPEP,peplength); //Make sure it presents at least 4 - 6 % of self peptides

	this->Initializepeppools();		//Initializes the self / non self pools
	this->CreateVirus(); 			//Creates a virus

	//KIRGenePool poolKIRGenes(MHCPool,self,KIRspecificity,specificityPEP,container);	//Completely random molecules
	KIRGenePool poolKIRGenes(MHCPool,self,HLA_C, KIRspecificity,specificityPEP,container); // Creates KIRs that recognize at least one self pMHC molecules


	//copy the kir genes in the KIRPool(to store the initial KIR ID and be able to keep track which ones are getting lost)
	for(unsigned int i= 0; i<poolKIRGenes.GetPoolSize(); i++)
		{
		KIRPool.GetGenes().push_back(poolKIRGenes.GetGenes().at(i));
		}


	for(unsigned int i= 0; i<KIRPool.GetPoolSize(); i++)
		{
		Gene Dummy;
		Dummy.SetGeneID(KIRPool.GetGenes().at(i));
		Dummy.CheckCross(MHCPool,self, KIRspecificity,specificityPEP,container);  //Checks if one particular KIR can also become licensed for other MHCs
		}

	container.Backupassoc(); //Creates a text file which contains the KIRs and all the MHCs each KIR is licensed for.

// Now create a 3D array of all the interaction between peptides, MHCs and KIRS

	this->CheckAssociations();
//	container.Printbind(); //Prints the 3D array (for testing)


	//initialize the population with the genes of the pools
	for (unsigned int i = 0; i< initHostPop; i++)
		{
		Host dummyhost(KIRLoci, MHCLoci, KIRspecificity, mutationRate, education,MHCPool,HLA_C, KIRPool,expressionExtraKIRs,container);
		hosts.push_back(dummyhost);
		}
	return true;

}



void World::Initializepeppools()
{
	vector<uint32_t> allsequences;										//A pool that contains all the peptides self and non self is created
	while (allsequences.size()<(selfpoolsize+nonselfpoolsize+1))
		{
		SequenceVar dummy(peplength,container);					// A peptide candidate is created
		uint32_t a=dummy.PepHash();
		if (std::find(allsequences.begin(),allsequences.end(),a) !=allsequences.end()) //Make sure it doesn't already exist in the pool, if so next candidate
				continue;
		else
			{
			int counter=0;											//Check if and by how many MHCs this peptide can be presented
			for (unsigned int o=0;o<MHCPool.GetPoolSize();o++)		//All MHC genes are created and checked;
				{
				Sequence anotherDummy(MHCPool.GetGenes().at(o),container);
				int pos=floor((anotherDummy.GetPepSize()-dummy.GetPepSize())/2);
				if (anotherDummy.Bind(dummy,container,pos)>=specificityPEP)
					{
					counter++;
					}
				else continue;
				}
															//Keep it if it can be presented by at least one else create a new candidate
			if (counter>0)
				{
				uint32_t a=dummy.Pep2Int(container);
				allsequences.push_back(a);
				}
			else continue;
			}
		}
	std::random_shuffle (allsequences.begin(),allsequences.end());		//It is shuffled, just to be sure
	for (int i=0;i<selfpoolsize;i++)										//And divided into the 2 pools self
		{
		self.push_back(allsequences.at(i));
		}
	for (int i=selfpoolsize;i<selfpoolsize+nonselfpoolsize+1;i++)				//and non self
		{
		nonself.push_back(allsequences.at(i));
		}

}


void World::CreateVirus()	//Creates the virus
{
	int i=1;
	Wildtype.SetViralParameters(i,nonself,container,deltaVirus,timeInfection);
}

void World::CheckAssociations()
{
	vector<uint32_t> MHC=MHCPool.GetGenes();
	vector<uint32_t> KIR=KIRPool.GetGenes();
	vector<uint32_t> peps=nonself;

	container.Createbind(MHC,KIR,peps);  //Creates a 3D vector accordingly with respect to the number of MHCs, KIRs and non self peptides

	//Now fill the vector

	//For each MHC:

	for (unsigned int k=0;k<MHC.size();k++)
		{
		//Fist create the MHC
		Sequence Dummy(MHCPool.GetGenes().at(k),container);

		//cout << Dummy.GetPep() <<"\n";

		//Check all the viral peptides

		for (unsigned int l=0;l<peps.size();l++)
			{
			Sequence thirddummy=Dummy;
			//Now the peptide
			int ispresented=0;
			Sequence applicant;
			applicant.Int2Pep(peps.at(l),container);

			//cout << applicant.GetPep() <<"\n";

			//Check if the peptide is presented
			int pos=floor((Dummy.GetPepSize()-applicant.GetPepSize())/2);

			//cout << pos << "\n";

			if (Dummy.Bind(applicant,container,pos)>=specificityPEP)
				{
				ispresented=1;
				thirddummy.Present(applicant,pos); //Also create the pMHC complex
				}

			//cout << ispresented <<"\n";

			//cout << thirddummy.GetPep() <<"\n";

			//Check all the KIRs

			for (unsigned int m=0;m<KIR.size();m++)
				{
				//Check, if it is not presented at all no KIR can bind with the complex
				if (ispresented==0)
					{
					container.Updatebind(k,l,m,0);
					}
				else	//If it is presented we will have to check if it binds with the pMHC complex
					{
					//Create the KIR
					Sequence anotherdummy(KIRPool.GetGenes().at(m),container);

					if (anotherdummy.Bind(thirddummy,container)>KIRspecificity)
						{
						container.Updatebind(k,l,m,1);

//						cout << anotherdummy.GetPep() << "\n";

//						cout << anotherdummy.Bind(thirddummy,container) <<"\n";

//						cout << "Yes\n";
						}
					else
						{
						container.Updatebind(k,l,m,0);

//						cout << anotherdummy.GetPep() << "\n";

//						cout << anotherdummy.Bind(thirddummy,container) <<"\n";

//						cout << "No\n";
						}
					}
				}
			}

	}
}





/*EVENT functions*/
/*Birth function: creates a child with the haplotype of one parent and another randomly chosen host*/
bool World::Birth(int index, Host& baby_host,Contain& cont)
{
	//check if host's age allows him to become a parent
	double ageDependentBirth = hosts.at(index).GetAgeDependentBirthRate();
	if(RandomNumberDouble()< birthRate*ageDependentBirth*(1-(hosts.size()/(maxHostPop*0.99753))))
		{
		int randomindex = RandomNumber(0,shuffledHosts.size()-1);
		//check if the potential parent is himself
		while(randomindex == index)
			{
			randomindex = RandomNumber(0,shuffledHosts.size()-1);
			}

	//choose randomly which parent is going to "donate" his/her mhc molecule
	int parent;
	if((RandomNumberDouble()<0.5))
		parent = index;
	else
		parent = randomindex;

	//int mhcFromPool = MHCPool.RandomlyPickGene(HLA_C);
	//int kirFromPool = KIRPool.RandomlyPickGeneHost newHost(KIRLoci, hosts.at(parent).mhcGenes ,mhcFromPool,KIRPool, hosts.at(index).kirGenes,hosts.at(randomindex).kirGenes,KIRspecificity,mutationRate,education);(false);
	vector<Gene> mchparent=hosts.at(parent).GetmhcGenes();
	vector<Gene> kirparent=hosts.at(parent).GetkirGenes();
	vector<Gene> kirrandom=hosts.at(randomindex).GetkirGenes();

//	Host testHost(KIRLoci, hosts.at(parent).mhcGenes ,mhcFromPool,KIRPool, hosts.at(index).kirGenes,hosts.at(randomindex).kirGenes,KIRspecificity,mutationRate,education,expressionExtraKIRs);
	Host testHost(KIRLoci,MHCLoci, mchparent ,MHCPool, HLA_C,KIRPool, kirparent ,kirrandom ,KIRspecificity,mutationRate,education,expressionExtraKIRs,cont);

/*		if(education==false)
		{
			baby_host.Copy(testHost);
			return true;
		}
		else
		{
			//be sure to create an offspring!: i.e. to have at least ONE functional KIR
			int score=testHost.CountFunctionalKIRs();

			if(score>0)
			{
				baby_host.Copy(testHost); //if there is already at least one functional KIR leave the loop
				return true;
			}
			else //but if it's not:'
			{
				while(score== 0) //pick an MHC allele until you get at least one functional KIR and the baby can be born
				{
					mhcFromPool = MHCPool.RandomlyPickGene(HLA_C);
					Host newBaby(KIRLoci, hosts.at(parent).mhcGenes ,mhcFromPool,KIRPool, hosts.at(index).kirGenes,hosts.at(randomindex).kirGenes,KIRspecificity,mutationRate,education,expressionExtraKIRs);
					if(newBaby.CountFunctionalKIRs()>0)
					{
						score= newBaby.CountFunctionalKIRs();
						baby_host.Copy(newBaby);
					}
				}
				return true;
			}
		}*/

	baby_host.Copy(testHost);
	return true;
	}
	return false;
}

/*Death function: according to an age-dependent rate, a host will be removed from the population*/
bool World::Death(int index)
{
	double intrinsicDeath;
	intrinsicDeath = deathRate*hosts.at(index).GetIntrinsicDeathRate();
	if(RandomNumberDouble()<intrinsicDeath)
		{
		return true;
		}
	return false;
}
/*Infect function: upon contact between two hosts, virus can spread*/

void World::Infect(int index)
{
	if(RandomNumberDouble()<infectionRate)
		{
		for(int i=0; i<contactRate; i++)
			{
			//pick random partner that is NOT yourself!
			int randomindex = RandomNumber(0,shuffledHosts.size()-1);
			while (randomindex == index)
				randomindex=RandomNumber(0,shuffledHosts.size()-1);
			//check whether partner is infectious
			if(hosts.at(randomindex).IsSusceptible() || hosts.at(randomindex).IsImmune())
				continue;
			//transmit the virus according to the type of infection
			else
				{
				if(hosts.at(randomindex).IsAcuteInfected())
					{
					if(RandomNumberDouble()<transmissionRateAcute)
						{
						hosts.at(index).InfectWith(hosts.at(randomindex).pathogen, simulationTime,specificityPEP,container);
						return;
						}
					}
				if(hosts.at(randomindex).IsChronicInfected())
					{
					if(RandomNumberDouble()<transmissionRateChronic)
						{
						hosts.at(index).InfectWith(hosts.at(randomindex).pathogen, simulationTime,specificityPEP,container);
						return;
						}
					}
				}
			}
		}
}


/*SIMULATION functions*/

void World::ShuffleHosts()
{
	//virtualHosts.clear();
	shuffledHosts.clear();

	for(unsigned int i=0; i<hosts.size(); i++)
		{
		shuffledHosts.push_back(i);
		}

   for(unsigned int i = 0 ; i < hosts.size() ; i ++ )
   	   {
      int j = RandomNumber(i,hosts.size()-1);
      // swap a[i] and a[j]
      int t = shuffledHosts.at(j);
      shuffledHosts.at(j) = shuffledHosts.at(i);
      shuffledHosts.at(i) = t;
   	   }

   /*
	while(virtualHosts.size()!=0)
	{
		int randomIndex = RandomNumber(0,virtualHosts.size()-1);
		shuffledHosts.push_back(virtualHosts.at(randomIndex));
		vector<int>:: iterator it = virtualHosts.begin();
		virtualHosts.erase(it + randomIndex);
	}*/
}

void World::Simulate()
{
	double lastOutfileTime = 0.0;
	//double lastStopOutfileTime = timeRecording;
	double lastBackupTime = 0.0;
	double lastAcuteInfectionTime = 0.0;
	isFileOpen = false;

	while(simulationTime <=timeEnd)
		{
		// printing out the backup files
		if(floor((simulationTime-lastBackupTime)*backupRate)>0 || simulationTime==0)
			{
			cout <<"\tSaving Backup\n";
			SaveBackupFile();
			SaveMap();
			//container.Backupassoc();
			lastBackupTime = simulationTime;
			}


		int number_babies = 0;
		int number_dead_people = 0;
		//cout<< "Time: " << simulationTime/YEAR <<endl;
		ShuffleHosts();
		vector<int>::iterator shuffledHostsit;

		//introduce infection


		if(floor(timeIntroducingInfection -simulationTime) > 0 && floor(timeIntroducingInfection -simulationTime) < 3.0*WEEK)
				IntroduceVirus(Wildtype);

		for(shuffledHostsit = shuffledHosts.begin(); shuffledHostsit!=shuffledHosts.end(); shuffledHostsit++)
			{
			//let event happen for every random host
			int index = *shuffledHostsit;
			Host babyHost;
			if(Birth(index, babyHost,container))
				{
				hosts.push_back(babyHost);
				number_babies ++;
				}
			if(hosts.at(index).IsSusceptible())
				Infect(index);
			if(Death(index))
				{
				hosts.at(index).SetDead();
				number_dead_people++;
				}
			//clear the infection
			if(hosts.at(index).IsAcuteInfected())
				{
				if((simulationTime -hosts.at(index).GetInfectionTime()) == (1.0 + 4.0*timeInfection)*WEEK)
					{
					hosts.at(index).ClearInfection(simulationTime,container);
					}
				}

			hosts.at(index).UpdateParameters(timeStep,simulationTime,timeInfection*WEEK);
			}
		//record those who are dying!
		/*	if(SaveAgeDyingHosts(lastOutfileTime, lastStopOutfileTime))
			lastStopOutfileTime = simulationTime;*/

		RemoveDeadHosts();
		SavePopulationSize(number_babies,number_dead_people);

		// printing out the gene files
		bool timeToPrintOut = floor((simulationTime-lastOutfileTime)*outfileRate)>0;
		if(simulationTime == 0.0 || timeToPrintOut)
			{
			cout <<"\tPrinting parameters\n";
			cout<< "Time: " << simulationTime/YEAR <<endl;
			SaveGenes();
			SaveParameters();
			populationSize.flush();
			if(populationSize.bad())
				cout << "something bad happened \n";
			lastOutfileTime = simulationTime;
			}

		simulationTime+=timeStep;
		}
	populationSize.close();

	cout << "The simulation was carried out with: \n";
	cout << "Maximal population size " << maxHostPop << " individuals \n";
	cout << "Initial population size " << initHostPop << " individuals \n";
	cout << "Introducing virus after " << timeIntroducingInfection/YEAR << " years \n";
	cout << "Mutation rate host " << mutationRate << "\n";
	cout << "Contacts per week: " << contactRate << "\n";
	cout << "HLA-C alleles distribution " << HLA_C<<"\n";
	cout << "MHC-KIR specificity " << KIRspecificity<<"\n";
	cout << "Tuning "<< education << "\n";
	cout << "Viral load: " << deltaVirus << "\n";
	cout << "Life time virus: " << timeInfection << "\n";
	cout << "Transmission rate acute infection "<< transmissionRateAcute << "\n";
}

void World ::IntroduceVirus(Virus& Vir)
{

	/*  This function introduces the viruses  */


	cout <<"\t Introducing the infection"<<endl;
	for(int i= 0; i<0.05* hosts.size(); i++)
		{
		int randomindex = RandomNumber(0,hosts.size()-1);
		while(!hosts.at(randomindex).IsSusceptible())
			{
			randomindex = RandomNumber(0,hosts.size()-1);
			}
		hosts.at(randomindex).InfectWith(Vir, simulationTime,specificityPEP,container);
		}

}


void World::RemoveDeadHosts()
{
	//deleting the hosts
	acute_infected = 0;
	chronic_infected = 0;
	immune = 0;
	int pos = 0;
	int max = hosts.size();
	while(pos < max)
	{
		if(hosts.at(pos).IsDead())
			{
			//remove the host
			hosts.erase(hosts.begin() + pos);
			max = hosts.size();
			}
		else
			{
			if(hosts.at(pos).IsAcuteInfected())
				acute_infected++;
			if(hosts.at(pos).IsChronicInfected())
				chronic_infected++;
			if(hosts.at(pos).IsImmune())
				immune ++;
			}
		pos ++;
	}
}

/*Functions which have to do the OUTPUT files*/

/*this function saves the genes of each host -> to keep track of MHC and KIR diversity*/
void World:: SaveGenes()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Genes.txt";
	
	genesFile.open(s.c_str(), ios::out);

	genesFile << "#Host ID \t Mhc 1 \t Mhc2 \t Mhc 3 \t Mhc4 \t Kir1 \t Kir2 \t Kir3 \t Kir4 \t Kir5 \t Kir6 \t Kir7 \t Kir8 \t Kir9 \t Kir10\t#functional Kirs\n";
	int index = 0;
	vector<Host>::iterator hostIt;
	if(printoption)
		{
		for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
			{
			genesFile << index << "\t";
			hostIt->SaveGenes(genesFile,container);
			index ++;
			}
		}
	else
		{
		for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
			{
			genesFile << index << "\t";
			hostIt->SaveGenesID(genesFile);
			index ++;
			}
		}

	genesFile.close();
}

/*This function keeps track (and saves) of the population size*/
void World::SavePopulationSize(int babies,int dead_people)
{
	populationSize <<simulationTime/YEAR <<"\t"<< hosts.size() <<"\t"<<babies <<"\t"<<dead_people<<"\t"<<acute_infected << "\t"<< chronic_infected << "\t"<< immune << "\t"<<"\n";
}

/*This function keeps track (ans saves) several parameters of the host and virus*/
void World::SaveParameters()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".txt";
	parameterFile.open(s.c_str(), ios::out);

	parameterFile << "#Host ID\t" << "Age \t" << "Infection Time \t" << "Infection type \t" << "Viral Death Rate\t"<<"Age Infection\t"<<"Age Clearance\t"<<"Virus Type\t"<<"Base Viral Load\n";
	vector<Host>::iterator hostIt;
	int index = 0;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
		{
		parameterFile << index <<"\t";
		hostIt->SaveParameters(parameterFile);
		index ++;
		}

	parameterFile.close();
}

bool World::SaveAgeDyingHosts(double lastOutfileTime, double lastStopOutfileTime)
{

	bool timeToPrintOut = floor((simulationTime-lastOutfileTime)*outfileRate);
	bool timeToStopOut = floor((simulationTime-lastStopOutfileTime)*outfileRate);

	if(timeToPrintOut>0)
		{
		isFileOpen = true;
		stringstream ss;
		ss << simulationTime/YEAR;
		string s(ss.str());
		s+=".Age.txt";
		dyingHosts.open(s.c_str(), ios::out);
		dyingHosts << "#Age \t KIR genes\t virus type\t viral load\t decoy ID\t onlyAcute?\t infection type\t infection time\t clearance time\n";
		}

	if(isFileOpen)
		{
		vector<Host>::iterator hostIt;
		//int index = 0;
		for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
			{
			if (hostIt->IsDead())
				{
				//write in a file their age, gene content, and infection type
				hostIt->SaveAgeDyingHost(dyingHosts);
				dyingHosts << simulationTime/YEAR <<"\n";
				}
			//index++;
			}
		if(timeToStopOut>0)
			{
			isFileOpen = false;
			dyingHosts.close();
			return true;
			}
		else
			return false;
		}
	return false;
}


/*Functions which do the BACKUP files*/

void World::SaveBackupFile()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Backup.data";
	backupFile.open(s.c_str(), ios::out);

	if(!backupFile.is_open())
		{
		throw OussException("Backup: Could not Open File");
		}
	else
		{
		cout <<"Saving backup file\n"<<endl;
		}

	backupFile << simulationTime<< "\t "<< KIRPool.GetPoolSize() << "\t";
	//save mhc gene pool
	for(unsigned int i = 0; i<MHCPool.GetPoolSize(); i++)
		{
		backupFile << MHCPool.GetGenes().at(i) << "\t";
		}
	//save kir gene pool
	for(unsigned int i = 0; i<KIRPool.GetPoolSize(); i++)
		{
		backupFile << KIRPool.GetGenes().at(i) << "\t";
		}
	backupFile << "\n";
	//Now for the self / non-self pool (line 2)
	for(unsigned int i=0;i< self.size();i++) 	//self
		{
		backupFile << self.at(i) << "\t";
		}
	for(unsigned int i=0;i< nonself.size();i++)	//non self
		{
		backupFile << nonself.at(i) << "\t";
		}
	backupFile << "\n";		//Another line
	//Here go the viruses and their individual pools (line 3)
	backupFile << Wildtype.GiveVirusID() << "\t";
		for (int k=0;k<nonselfpoolsize;k++)
			{
			backupFile << Wildtype.GivePep(container).at(k) << "\t";
			}
	backupFile << "\n";		//Another line and we start with the iteration
	vector<Host>::iterator hostIt;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
		{
		hostIt->SaveBackupHost(backupFile);
		}
	backupFile.close();
}

void World ::SaveMap()
{
	system("rm -f Sequence.map");   //Erases the previous map file
	container.backupmap();			// Saves the file contents to the new Sequence.map file
	cout << "A backup of the map is made";
};


void World :: LoadMap()
{
	container.restoremap();
}

void World ::LoadBackupFile(const string& fileName)
{
	int l = 0;
	MHCPool.FillGenePoolFromContainer(container);
	int poolSize = MHCPool.GetGenes().size();
	//int poolSize=14;  //JUST TO CHECK IT WILL BE CHANGED
	int kirPoolSize = 0;
	MHCPool.GetGenes().clear();
	string sline;
	string slineIt;
	vector<Host>::iterator hostIt;
	backupFile.open(fileName.c_str(), ios::in);
	if(!backupFile.is_open())
		{
		throw OussException("Backup: Could not Open File");
		}
	while(!backupFile.eof())
		{
		if(l==0)
			{
			getline(backupFile,sline);
			stringstream ssline(sline);
			ssline >>simulationTime;
			ssline >>kirPoolSize;
			//restore gene pools
			for(int i = 0; i<poolSize; i++)
				{
				uint32_t mhc_gene;
				ssline>>mhc_gene;
				MHCPool.GetGenes().push_back(mhc_gene);
				}
			for(int i = 0; i<kirPoolSize; i++)
				{
				uint32_t kir_gene;
				ssline>>kir_gene;
				KIRPool.GetGenes().push_back(kir_gene);
				}
			}
		else if (l==1)
			{
			getline(backupFile,sline);
			stringstream ssline(sline);
			for(int i=0;i< selfpoolsize;i++) 	//self
				{
				uint32_t selfp;
				ssline>>selfp;
				self.push_back(selfp);
				}
			for(int i=0;i< nonselfpoolsize;i++)	//non self
				{
				uint32_t nselfp;
				ssline>>nselfp;
				nonself.push_back(nselfp);
				}
			}
		else if (l==2)
			{
			getline(backupFile,sline);
			stringstream ssline(sline);
			uint32_t virID;
			ssline>>virID;
			Wildtype.SetVirusID(virID);
			vector<uint32_t> nspeps;
			for (int k=0;k<nonselfpoolsize;k++)
				{
				uint32_t pnum;
				ssline>>pnum;
				nspeps.push_back(pnum);
				}
			Wildtype.SavePep(nspeps,container);
			}
		else if (l>2)
			{
			string sline;
			getline(backupFile, sline);
			if(sline.size()!=0)
				{
				//cout <<"sline1:  "<<sline <<endl;
				Host tempHost;
				tempHost.RestoreHost(sline);

				//cout <<"restored:  "<<tempHost.GetAge() <<"\t" << tempHost.GetInfectionType() << "\t" << tempHost.GetInfectionTime() << "\t" <<tempHost.GetClearanceTime() <<endl;
				hosts.push_back(tempHost);

				}
			}
		l++;
	}
	backupFile.close();
	this->CheckAssociations();
}

void World::Loadmatrix(string& mfile)
{
	container.loadfile(mfile);
};
