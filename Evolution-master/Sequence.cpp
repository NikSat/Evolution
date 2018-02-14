/*
 * Sequence.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: nikolaos
 */


#include "Sequence.h"

/* Constructor: creates MHC molecules as a string of random characters, length 16*/
SequenceMHC::SequenceMHC(Contain& cont)
{
	string library=cont.getaminolibrary();
//	cout << library<< endl;
	for (int i=0; i< MOLECULE_LENGTH; i++){
		char e=library[RandomNumber(0,19)];
		peps.insert(peps.begin()+i,1,e);
		}
}


/* This creates peptide sequences of variable length */

SequenceVar::SequenceVar(int length,Contain& cont)
{
	string library=cont.getaminolibrary();
	for (int i=0; i< length; i++)
		{
		char e=library[RandomNumber(0,19)];
		peps.insert(peps.begin()+i,1,e);
		}

}

/* This functions checks whether two strings are equal*/

bool Sequence:: operator == (Sequence& rhs)
{
	if(rhs.GetPepSize() != peps.size())
		return false;
	else if (strcmp ((char*)rhs.GetPep().c_str(),(char*)peps.c_str()) !=0)
		return false;
	else
		return true;

}

/* This function gives each string an integer "ID" (Jenkins hash) and appends it to the map*/

uint32_t Sequence :: Pep2Int(Contain& cont)
{
	    uint32_t moleculeID, i;
	    for(moleculeID= i = 0; i < peps.size(); ++i)
	    {
	    	moleculeID += peps[i];
	    	moleculeID += (moleculeID << 10);
	    	moleculeID ^= (moleculeID >> 6);
	    }
	    moleculeID += (moleculeID << 3);
	    moleculeID ^= (moleculeID >> 11);
	    moleculeID += (moleculeID << 15);

	    cont.Updatemap(moleculeID,peps);	// Stores them in the map
	    return moleculeID;
}

/* This function does just the Hash and does not append to the map */

uint32_t Sequence :: PepHash()
{
	    uint32_t moleculeID, i;
	    for(moleculeID= i = 0; i < peps.size(); ++i)
	    {
	    	moleculeID += peps[i];
	    	moleculeID += (moleculeID << 10);
	    	moleculeID ^= (moleculeID >> 6);
	    }
	    moleculeID += (moleculeID << 3);
	    moleculeID ^= (moleculeID >> 11);
	    moleculeID += (moleculeID << 15);

 return moleculeID;
}

/*This function fetches the sequence that has this "ID" from the map */
void Sequence :: Int2Pep(uint32_t moleculeID,Contain& cont)
{
	peps=cont.getseq(moleculeID);
}

/*This is the new function based on the pseudo-Blosum Matrix  */

// Update: The library is now loaded from the first line of the matrix (to avoid errors from both Blosum and JM matrix).

int Sequence :: Bind (Sequence& anotherPeps,Contain& cont)
{
	string library=cont.getaminolibrary();
	int maxStrength=0;
	for (unsigned int i=0;i<peps.size();i++){
	int m=cont.getscore(library.find(peps[i]),library.find(anotherPeps.GetPep()[i]));
	// int m=myWorld.affinity[library.find(peps[i])][library.find(anotherPeps.GetPep()[i])][i]; A 3D affinity matrix that I will not use after all
	maxStrength=maxStrength+m;}
	return maxStrength;
}

int Sequence :: Bind (Sequence& anotherPeps,Contain& cont,int pos)
{
	string library=cont.getaminolibrary();
	int maxStrength=0;
	for (unsigned int i=0;i<anotherPeps.GetPep().size();i++)
		{
		int m=cont.getscore(library.find(peps[pos+i]),library.find(anotherPeps.GetPep()[i]));
		maxStrength=maxStrength+m;
		}
	return maxStrength;
}


void Sequence :: PointMutation(Contain& cont) //This function makes a random point mutation on a molecule (either MHC or KIR)
{
	string library=cont.getaminolibrary();
	int position = RandomNumber(0,peps.size()-1);
	peps.at(position)=library[RandomNumber(0,library.size()-1)];
	cout <<"point mutation just happened\n";
	}

void Sequence :: Express(Sequence& anotherPeps) //This function expresses the 3 char peptide at a random place on the MHC between positions 3 and 8
{
	int k=RandomNumber(3,8);
	for (int n=0;n<3;n++){
	peps[k+n]=anotherPeps.GetPep()[n];
	}
}

void Sequence :: PrintSeq()

	{
		cerr << peps << " ";
	}
void Sequence :: Present (Sequence& anotherPeps,int pos)	//This function takes a peptide and fusses it with the already existent MHC
{
	for (unsigned int n=0;n<anotherPeps.GetPepSize();n++)
	{
	peps[pos+n]=anotherPeps.GetPep()[n];
	}
}


