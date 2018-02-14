/*
 * Sequence.h
 *
 *  Created on: Apr 25, 2013
 *      Author: nikolaos
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_


#include "MathFunctions.h"
#include "Contain.h"
#include <cmath>
#include <utility>
#include <string.h>
#define MOLECULE_LENGTH 16
using namespace std;

class Sequence {
	public:
	Sequence (){};
	Sequence (uint32_t molecule_number,Contain& cont){Int2Pep(molecule_number,cont);}
	virtual ~Sequence (){};
	string& GetPep(){return peps;};
	unsigned int GetPepSize(){return peps.size();}
	bool operator==(Sequence& rhs);
	uint32_t Pep2Int(Contain& cont); //to test
	uint32_t PepHash(); //to test
	void Int2Pep(uint32_t molecule_number,Contain& cont);//to test
	/*int AdjacentMatch(Sequence& kirBits); //works
	int HammingMatch(Sequence& kirBits); // We will have to see about this one// works*/
	int Bind (Sequence& anotherPeps,Contain& cont);
	int Bind (Sequence& anotherPeps,Contain& cont,int pos);
	void PointMutation(Contain& cont); //works
	void Present(Sequence& anotherPeps,int pos);
	void PrintSeq();
	void Express(Sequence& anotherBits);

protected:
	string peps;
};

/* Random 16 character length sequence */

class SequenceMHC : public Sequence
{
public:
	SequenceMHC(Contain& cont);
	virtual ~SequenceMHC(){};
};


/* Sequence of variable length */

class SequenceVar: public Sequence
{
public:
	SequenceVar(int length,Contain& cont);
	virtual ~SequenceVar(){};
};

#endif /* SEQUENCE_H_ */
