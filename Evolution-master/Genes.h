/*
 * Genes.h
 *
 *  Created on: Apr 20, 2011
 *      Author: paola
 */

#ifndef GENES_H_
#define GENES_H_

#include "Sequence.h"



class GenePool
{
public:
	GenePool(){};
	virtual ~GenePool(){};
	float GetAlleleFreq(int alleleIndex); //works
	bool GeneAlreadyInPool(uint32_t geneID); //works
	unsigned int GetPoolSize(){return genes.size();}
	uint32_t RandomlyPickGene(bool distribution); //works
	vector<uint32_t>& GetGenes(){return genes;}
	void AddGene(uint32_t gid);

protected:
	vector<uint32_t> genes;//a vector of integers
	fstream genePoolFile;
};



class Gene {


public:
	Gene();
	virtual ~Gene(){};
	bool operator ==(Gene& rhs);
	void SetGeneFunctionality(bool functionality);
	void SetGeneExpression(bool expression);
	bool IsFunctional(){return functional;};
	bool IsExpressed(){return isExpressed;}
	Gene& Copy(Gene& rhsGene);//works
	void SetGeneID(uint32_t _ID){geneID = _ID;};
	uint32_t GetGeneID(){return geneID;}
	void Mutate(Contain& cont,GenePool& kirPool);//works
	void Mutate(uint32_t new_id);
	bool IsGeneUnique(vector<Gene>& genePool, int counter);//works
	int BindMolecule(Gene& anotherMolecule, Contain& cont);
	int SpecBindMol(Sequence& anotherSeq,int pos,Contain& cont);
	void PresentPeptide (Sequence& anotherSeq,int pos,Contain& cont);
	void CheckCross(GenePool& mhcPool,vector<uint32_t> pool, double specificity,double pepspecificity,Contain& cont);
//	void CheckCross(GenePool& mhcPool,vector<uint32_t> pool, double specificity,double pepspecificity,Contain& cont,fstream& temp);

	void SaveBackupGenes(fstream& backupFile);//works
	void SaveGenes (fstream& outfile);
	string RestoreGenes(stringstream& sgline);


protected:
	uint32_t geneID;
	bool functional;
	bool isExpressed;
};


class MHCGenePool: public GenePool
{
public:
	MHCGenePool(){};
	void FillGenePoolFromContainer(Contain& cont);
//	void FillGenePoolFromContainer(Contain& cont,vector<uint32_t> pool,double cesp,int le);
	virtual ~MHCGenePool(){};
};

class KIRGenePool: public GenePool
{
public:
    KIRGenePool(GenePool& mhcPool,vector<uint32_t> pool, bool dist, double specificity,double pepspecificity,Contain& cont);
//	KIRGenePool(GenePool& mhcPool,vector<uint32_t> pool, double specificity,double pepspecificity,Contain& cont);
	virtual ~KIRGenePool(){};

};

#endif /* GENES_H_ */
