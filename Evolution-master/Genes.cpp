/*
 * Genes.cpp
 *
 *  Created on: Apr 20, 2011
 *      Author: paola
 */

#include "Genes.h"


/* FUNCTIONS OF CLASS GenePool*/

/* MHC Related functions */

/* FUNCTION 1: Constructs a pool with HLA C Alleles (14) */
void MHCGenePool::FillGenePoolFromContainer(Contain& cont)
{
	cerr << "Filling MHC pool with random generated strings" << endl;
	cerr << "Working";
	for(int i=0; i<14; i++)
		{
		SequenceMHC dummy(cont);
		uint32_t geneID = dummy.Pep2Int(cont);
		if(!GenePool::GeneAlreadyInPool(geneID)){
			GenePool::genes.push_back(geneID);
		cerr << ".";
		};
	}
	cerr << endl;
}


/* KIR Related functions*/

/* This function fills the KIR Gene pool, and prints the Sequences and ID's*/

KIRGenePool :: KIRGenePool(GenePool& mhcPool,vector<uint32_t> pool, bool dist, double specificity,double pepspecificity,Contain& cont)
{
		/*  NEW FUNCTION  SERIES*/

	cerr << "Filling initial KIR pool" <<endl;
	cerr <<"Working";



	/*FUNCTION 3: Each KIR should recognize at least one self pMHC from the pool for each MHC
	 * 			  There should be an initial pool of 140 KIRs, 10 licensed for each MHC   */

	for (unsigned int o=0;o<mhcPool.GetPoolSize();o++)		//Start with the pool of MHCs
		{
			Sequence Dummy(mhcPool.GetGenes().at(o),cont); 	//The MHC gene is created;
			vector<uint32_t> presented;						//All the peptides it can present go here
			for (unsigned int i=0;i<pool.size();i++)
				{
				Sequence applicant;
				applicant.Int2Pep(pool.at(i),cont);
				int pos=floor((Dummy.GetPepSize()-applicant.GetPepSize())/2);
				if (Dummy.Bind(applicant,cont,pos)>=pepspecificity)
					{
					presented.push_back(pool.at(i));
					}
				else continue;
				}								//We have all the presented peptides for this MHC molecule
			int count=0;
			while (count<10)				// Ten KIRs that each one recognizes most of the presented peptides must be created by this process
				{
				SequenceMHC anotherdummy(cont);	//A random KIR applicant is created
				if (std::find(genes.begin(),genes.end(),anotherdummy.PepHash()) !=genes.end())	// First check whether the applicant is already in the pool
				continue;
				else
					{
					int anothercounter=0;  //  The number of pMHCs the applicant recognizes
					for (unsigned int i=0;i<presented.size();i++) 	//Each one of the pMHC complexes is created and the applicant is cross checked
						{
						Sequence tempdummy(presented.at(i),cont);
						Sequence thirddummy=Dummy;
						int pos=floor((Dummy.GetPepSize()-tempdummy.GetPepSize())/2);
						thirddummy.Present(tempdummy,pos);
						double tempsco=anotherdummy.Bind(thirddummy,cont);
						if (tempsco>=specificity) anothercounter++;
						else continue;
						}
/*******************
			//Different Education Processes: Possible implementations of different education rules here
//					if (anothercounter>=floor(0.75*presented.size()))
//					if (anothercounter>=floor(0.50*presented.size()))
//					if (anothercounter>=floor(0.25*presented.size()))
 	 	 	 	 	 	 	 	 	 	 ***********************************/
					if (anothercounter>=1)
						{
						genes.push_back(anotherdummy.Pep2Int(cont));		//If it recognizes at least one of the presented peptides keep it
						count++;											//Update the counter
						cont.Updateassoc(anotherdummy.PepHash(),mhcPool.GetGenes().at(o)); //The Association map is also updated
						}
					else continue;
					}
				}
			cerr << ".";
		}						//This should give 140 KIRs



	cerr <<endl;

	// Now the equivalent printing/storing part of the function


	genePoolFile.open("GenePool.data", ios::out);
	genePoolFile << "# MHC \t KIR \n";
	for(int i=0; i<mhcPool.GetPoolSize(); i++)
		{
		for(unsigned int j=0; j<10; j++)
			{
			Sequence mhc_molecule(mhcPool.GetGenes().at(i),cont);
			Sequence kir_molecule(genes.at((i*10)+j),cont);
			uint32_t kir_id = genes.at((i*10)+j);
			genePoolFile << mhc_molecule.Pep2Int(cont) << "\t" << kir_id <<endl;
			cout << mhc_molecule.Pep2Int(cont) << "\t" ;
			cout << mhc_molecule.GetPep() <<"\t";
			cout << "\t"<< kir_id << "\t";
			cout <<kir_molecule.GetPep() <<endl;
			}
		}



	genePoolFile.close();

	cout <<endl;

}



/* This function assures that every single gene in the pool is unique*/
bool GenePool:: GeneAlreadyInPool(uint32_t geneID)
{
	for(unsigned int i=0; i<genes.size(); i++)
		{
		uint32_t number = genes.at(i);
		if(number == geneID)
			return true;
		}
	return false;
}


/* This function sets each MHC allele with a predefined frequency. (adopted from dbMHC Project)*/
float GenePool :: GetAlleleFreq(int alleleIndex)
{
	static float allele_freq=0;
	allele_freq = 0.2*exp(-0.22*alleleIndex); // probability distribution for hla-c alleles in european population
	return allele_freq;
}


/*This function returns an allele "randomly". i.e. either according to the HLA -C distribution or with equal probability*/
uint32_t GenePool :: RandomlyPickGene(bool distribution)
{
	int j;
	if(distribution==true)
		{
		static double probabilities[14];
		for (int i=0; i<14; i++)
			{
			probabilities[i] = GetAlleleFreq(i);
			}
		j = RandomNumber(probabilities, 13);
		}
	else
		{
		j = RandomNumber(0,genes.size()-1);
		}
	return genes.at(j);
}

void GenePool::AddGene(uint32_t gid)
{
	genes.push_back(gid);
}

/* FUNCTIONS OF CLASS Gene*/
/*constructs default genes*/
Gene::Gene ()
{
	functional = true;
	isExpressed = true;
	geneID = 0;
}

/*This function checks whether two genes are equal*/
bool Gene :: operator == (Gene& rhs)
{
	if(geneID == rhs.geneID)
		return true;
	else
		return false;
}

void Gene :: SetGeneFunctionality(bool functionality)
{
	functional = functionality;
}

void Gene:: SetGeneExpression(bool expression)
{
	isExpressed = expression;
}

/*This function performs the different mutational operators for a Gene*/

/* Point mutation - not used*/

void Gene :: Mutate(Contain& cont,GenePool& kirPool)
{
	Sequence molecule(geneID,cont);
	molecule.PointMutation(cont);
	geneID=molecule.Pep2Int(cont); //A new gene is created and appended to the map
	kirPool.AddGene(geneID);// Now this gene is also appended to the gene pool
}

void Gene :: Mutate(uint32_t new_gene_id)
{
	geneID = new_gene_id;
	//Sequence molecule(geneID);
	//molecule.PointMutation();
}



/*This function counts unique genes within a given pool*/
bool Gene :: IsGeneUnique(vector<Gene>& genePool, int counter)
{
	cout <<"this is the counter : "<<counter << "\n";
	for (int i = 0; i<counter; i++)
		{
		Gene currentGene = genePool.at(i);
		if(geneID == currentGene.geneID)
			return false;
		}
	return true;
}

int Gene::BindMolecule(Gene& anotherMolecule, Contain& cont)
{

	Sequence anotherSequence(anotherMolecule.geneID,cont);
	Sequence currentSequence(geneID,cont);

	int bindingStrength = currentSequence.Bind(anotherSequence,cont);
	return bindingStrength;
}

/* Specific binding at a certain position */

int Gene::SpecBindMol(Sequence& anotherSequence,int pos,Contain& cont)
{
Sequence currentSequence(geneID,cont);
int a=currentSequence.Bind(anotherSequence,cont,pos);
return a;
}

void Gene::PresentPeptide(Sequence& anotherSeq,int pos,Contain& cont)
{
	Sequence currentSequence(geneID,cont);			//The Sequence is made
	currentSequence.Present(anotherSeq,pos);		// The fussion is made
	geneID=currentSequence.Pep2Int(cont);			//Hashed, appended and new ID given
}

void Gene::SaveBackupGenes(fstream& backupFile)
{
	backupFile << geneID << "\t" << functional << "\t"<<isExpressed <<"\t";

}

void Gene::SaveGenes(fstream& outfile)
{
	outfile << geneID << "\t";

}

Gene& Gene::Copy(Gene& rhsGene)
{
	this->functional = rhsGene.functional;
	this->geneID = rhsGene.geneID;
	this->isExpressed = rhsGene.isExpressed;
	return *this;
}

string Gene::RestoreGenes(stringstream& sgline)
{
	sgline >> geneID;
	sgline >> functional;
	sgline >> isExpressed;
	//cout <<geneID <<"\t" <<functional << "\t";
	string gstring = sgline.str();
	return gstring;
}

/* Checks if KIR genes are also licensed for other MHC molecules */

void Gene::CheckCross(GenePool& mhcPool,vector<uint32_t> pool, double specificity,double pepspecificity,Contain& cont)
{
	Sequence anotherdummy(geneID,cont);
	//anotherdummy.PrintSeq();
	//cerr << endl;


	for (unsigned int o=0;o<mhcPool.GetPoolSize();o++)
		{
		Sequence Dummy(mhcPool.GetGenes().at(o),cont); 	//The MHC gene is created;
		vector<uint32_t> presented;						//All the presented peptides go here
		for (unsigned int i=0;i<pool.size();i++)
			{
			Sequence applicant;
			applicant.Int2Pep(pool.at(i),cont);
			int pos=floor((Dummy.GetPepSize()-applicant.GetPepSize())/2);
			if (Dummy.Bind(applicant,cont,pos)>=pepspecificity)
				{
				presented.push_back(pool.at(i));
				}
			else continue;
			}								//We have all the presented peptides for this HLA C molecule

		int anothercounter=0;
		for (unsigned int i=0;i<presented.size();i++) 	//Each one of the pMHC complexes is created and the KIR is cross checked
			{
			Sequence tempdummy(presented.at(i),cont);
			Sequence thirddummy=Dummy;
			int pos=floor((Dummy.GetPepSize()-tempdummy.GetPepSize())/2);
			thirddummy.Present(tempdummy,pos);
			int score=anotherdummy.Bind(thirddummy,cont);
			if (score>=specificity)
				{
				anothercounter++;
				}
			else continue;

			}

	//	cerr << anothercounter;
/*******************
//Different Education Processes: Possible implementations of different education rules here
//			if (anothercounter>=floor(0.75*presented.size()))
//			if (anothercounter>=floor(0.50*presented.size()))
//			if (anothercounter>=floor(0.25*presented.size()))
 	 	 	 	 	 	 	 	 ***********************************/

		if (anothercounter>=1)
			{
			cont.Updateassoc(anotherdummy.PepHash(),mhcPool.GetGenes().at(o)); //If it recognizes most of the presented peptides (>80%) update the Association map
			//	cerr << " " << "yes" << endl;
			}
			//else cerr << " " << "no" << endl;
		}
}

