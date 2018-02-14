/*
 * Contain.cpp
 *
 *  Created on: May 7, 2013
 *      Author: nikolaos
 */

#include "Contain.h"


// This function loads the matrix. In this simulation we use the MJ matrix but it can accept other matrixes as long as the first line and column contain the amino acid library
//It assumes the first line and column contain the one letter abbreviation of the amino acid name.
// The amino acid library i.e. the names of the amino acids and their order in the matrix is taken from the first line

void Contain :: loadfile(string& mfile)
{
	std::ifstream mat (mfile.c_str()); //Import all the data from the file and on the stream lib
	if(mat)
		{					// Check the iostate of the stream
		string line;
		int n=-1;			// Counter of lines
		while(getline( mat,line ))
			{
			stringstream ss(line);
			if (n<0)
				{
				char charlib[20];
				for (int i=0;i<20;i++)
					{
					ss >> charlib[i];
					}
				string amin(charlib);
				aminolibrary=amin;
				}
//			cerr << aminolibrary << "\n";
			if(n>=0)
				{				//Ignores the first line that contains the names
				int m=-1;	//Counter of columns
				while(!ss.eof())
					{
					char buff[10];
					ss.getline(buff,10,' ');
					if(m>=0)
						{						//Ignores the first column (also names)
						affinity[n][m]=atoi(buff); //Gives the values to the matrix
						m++;
						}
					else m++;
					}
				n++;
				}
			else n++;
			}

		}
	else
		{
		cout << "Error: no data or unreadable matrix file"; 	//If the state is bad inform and
		exit (0);
		}; 											//exit the program

};

//Updates the map

void Contain ::Updatemap(uint32_t mID,string mStr)
{
	libr.insert(std::pair <uint32_t,string> (mID,mStr) );
};


//Gives the sequence from the map

string Contain :: getseq(uint32_t ID)
{
	string seq=libr[ID];
	return seq;
}

//Gives binding score extracted from the matrix

int Contain ::getscore(int a,int b)
{
	return affinity[a][b];
};

//For testing purposes (to test different matrixes eg MJ or PAM or BLOSUM) the maximum, minimum element as well as the maximum binding score for 16 char sequences from the matrix can be calculated with these functions

int Contain::getmaxelement()
{
	maxelement=*max_element(affinity[0],affinity[0]+400);
	return maxelement;
};

int Contain::getminelement(){

	minelement=*min_element(affinity[0],affinity[0]+400);
	return minelement;
}

int Contain::getmaxscore(){
	maxelement=*max_element(affinity[0],affinity[0]+400);
	maxmatrix=16*maxelement;
	return maxmatrix;
};

//Backup and restore functions

void Contain :: backupmap()
{
	std::ofstream outf ("Sequence.map", std::ofstream::out);	// Stores both hash values and stings in the file Sequence.map
	for (std::map<uint32_t,string>::iterator it=libr.begin(); it!=libr.end(); ++it){
		outf << it->first << " " << it->second <<endl ;}
	outf.close();
}

void Contain :: restoremap()
{
	uint32_t hashnum;
	string seqstring;

	ifstream load ("Sequence.map");		// Loads the backup map
	if(load)
		{
		string line;
		while(getline (load,line))
			{
			stringstream ss(line);
			ss >> hashnum >> seqstring;			//Gets the hash number and the string
			libr.insert (std::pair<uint32_t,string>(hashnum,seqstring));	// Stores them in the map
			}
		load.close();
		cout << aminolibrary;
		}
		else
			{
			cout << "Error:Non existent or unreadable map file" << endl;
			cout << "The file containing the hashed values named Sequence.map is essential" << endl;
			exit(0);
			}	//Inform and force quit
}

// Copy functions, again for some tests

void Contain::copymatrix(Contain& anothercont)
{
	for (int i=0;i<20;i++)
		{
		for (int j=0;j<20;j++)
			{
			affinity[i][j]=anothercont.getscore(i,j);
			}
		}
}

void Contain::copylib(Contain& anothercont)
{
	aminolibrary=anothercont.getaminolibrary();
}

void Contain::Updateexpected(double num)
{
	expected=num;
}

//Updates the association map, which contains for which MHC each KIR is licensed for

void Contain ::Updateassoc(uint32_t mID,uint32_t mMID)
{

	assoc.insert(std::pair <uint32_t,uint32_t> (mID,mMID) );

};

//Association file: a print of the association map

void Contain :: Backupassoc()
{
	fstream genePoolFile;
	genePoolFile.open("Association.data", ios::out);
	genePoolFile << "# KIR \t MHC \n";
	for (std::multimap<uint32_t,uint32_t>::iterator it=assoc.begin(); it!=assoc.end(); ++it)
		{
		uint32_t alpha=it->first;
		uint32_t beta=it->second;
		genePoolFile << libr[alpha] << " "<< libr[beta] << endl;
		}
	genePoolFile.close();

}

//Creates the bind 3D vector which contains all the peptide - MHC - KIR interactions

void Contain::Createbind(vector<uint32_t> _MHC,vector<uint32_t> _KIR,vector<uint32_t> _peps)
{
	// Store the sizes of the MHC, KIR pool
	int MHCsize=_MHC.size();
	int KIRsize=_KIR.size();
	// Keep the vectors for lookup
	MHC=_MHC;
	KIR=_KIR;
	int Pepsize=_peps.size();
	Peptide=_peps;

	// Now make the 3D vector

	bind.resize(MHCsize);
	for (int i = 0; i < MHCsize; ++i)
		{
	    bind[i].resize(Pepsize);
	    for (int j = 0; j < Pepsize; ++j)
	    	{
	    	bind[i][j].resize(KIRsize);
	    	}
		}
}

//Puts values on the vector

void Contain::Updatebind(int a,int b,int c,int value)
{
	bind[a][b][c]=value;
}

//Uses the MHC, KIR and peptide vectors to get the index, this is the position in the 3D vector of each MHC, KIR or peptide

int Contain::Getposition(vector<uint32_t>& pool,uint32_t x)
{
	int i = 0;
	for( vector<uint32_t>::iterator it = pool.begin(); it != pool.end() ; ++ it )
		{
		if( *it == x ){ return i; }
		i ++;
		}
	return -1;
	/*
	vector<uint32_t>::iterator it;
	it=find(pool.begin(), pool.end(), x);
	int pos=distance(pool.begin(),it);
	return pos;*/
}

//Retrieves the interaction value: 1 if pMHC-KIR bind, 0 in all other cases (pMHC - KIR not binding or peptide not presented at all)

int Contain::GetbindValue(uint32_t a,uint32_t b,uint32_t c)
{

	//First get the indexes (positions) in the 3D vector
	int alpha=Getposition(MHC,a);

	int beta=Getposition(KIR,b);

	int ci=Getposition(Peptide,c);

//	cout <<"MHC position\t "<<alpha<<"\n";
//	cout <<"KIR position\t "<< beta<<"\n";
//	cout << "Peptide position\t"<<ci<<"\n";


	int value=bind[alpha][ci][beta];

//	cout <<"Binds?\t" <<value <<"\n";

	//Return the value
	return value;
}

//Not used - just one virus

void Contain::UpdateViralPool(uint32_t _ID,vector<uint32_t> _pool)
{
	viralpools.insert(std::pair <uint32_t,vector<uint32_t> > (_ID,_pool) );
}

vector<uint32_t> Contain::GetViralPool(uint32_t _ID)
{
	vector<uint32_t> pool=viralpools[_ID];
	return pool;
}


// This function is just for testing: prints 140 files (one for each KIR) with the KIRs interactions with the peptides and MHC molecules

void Contain::Printbind()
{
	for (int i=0;i<KIR.size();i++)
		{
		string indKIR=libr[KIR.at(i)];
		fstream bindfile;
		string s=indKIR;
		s+=".txt";
		bindfile.open(s.c_str(), ios::out);
		//Print all the Peptides at the first line
		for (int o=0;o<Peptide.size();o++)
			{
			string title=libr[Peptide.at(o)];
			bindfile << "\t"<< title;
			}
		bindfile << "\n";

		//Now for each MHC print the recognition

		for (int k=0;k<MHC.size();k++)
			{
			bindfile << libr[MHC.at(k)] ;
			for (int l=0;l<Peptide.size();l++)
				{
				bindfile <<"\t" << bind[k][l][i];
				}
			bindfile << "\n";
			}
		bindfile.close();
		}
}





