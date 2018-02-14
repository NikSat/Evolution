/*
 * Contain.h
 *
 *  Created on: May 7, 2013
 *      Author: nikolaos
 */

#ifndef CONTAIN_H_
#define CONTAIN_H_

#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <iostream>
using namespace std;

class Contain {
public:


// Functions of the class Contain
	void loadfile(string& mfile);		//loads the affinity matrix
	void Updatemap(uint32_t mID,string mStr); //adds elements to the map
	string getseq(uint32_t);   //gets the corresponding sting
	int getscore(int a,int b);	 //gets the binding score between two aminoacids
	int getmaxscore();  		// gives the maximum score two sequences can get (this can also be set)
	int getmaxelement(); 		//gives the maximum possible binding score between two aminoacids
	int getminelement();		// gives the least possible binding score between two aminoacids
	void backupmap();
	void restoremap();
	string getaminolibrary(){return aminolibrary;};
	void copymatrix(Contain& anothercont);
	void copylib(Contain& anothercont);
	void Updateassoc(uint32_t mID,uint32_t mStr); //Updates the association map
	void Updateexpected(double num); 			//Updates the expected presentation threshold
	double Getexpected(){return expected;};							//Gives the expected threshold
	void Backupassoc();
	multimap<uint32_t,uint32_t> GiveAssoc(){return assoc;};
	void Createbind(vector<uint32_t> _MHC,vector<uint32_t> _KIR,vector<uint32_t> _peps);
	void Updatebind(int a,int b,int c,int value);
	int GetbindValue(uint32_t a,uint32_t b,uint32_t c);
	int FastGetValue(int a,int b,int c){return bind[a][c][b];};
	void UpdateViralPool(uint32_t _ID,vector<uint32_t> _pool);
	vector<uint32_t> GetViralPool(uint32_t _ID);
	int Getposition(vector<uint32_t>& pool,uint32_t x);
	void Printbind();


protected:			 //the two containers are now protected

	int affinity[20][20];             //the matrix that contains the interaction scores
	map<uint32_t,string> libr; 		// The map containing all the strings
	map<uint32_t,vector<uint32_t> > viralpools; 		// The map containing all the viruses NOT used
	int maxmatrix; 					// The maximum binding strength possible
	int maxelement;
	int minelement;
	string aminolibrary;
	multimap<uint32_t,uint32_t> assoc; //This map associates KIRs with MHCs. The initial pool as well as the newly created KIRs are scanned
								// in order to find for how many MHCs they can get licensed for. KIRs are the keys, MHCs are the values.
	double expected;			//The expected percentage of presentation
	vector<vector<vector<int> > > bind; //The 3D array that stores all the interactions between peptides, MHCs and KIRs
	vector<uint32_t> MHC;
	vector<uint32_t> KIR;
	vector<uint32_t> Peptide;

};



#endif /* CONTAIN_H_ */
