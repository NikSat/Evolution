//============================================================================
// Name        : CoevolutionKIR.cpp
// Author      : P.Carrillo-Bustamante and N. Satravelas
// Version     : 3.0 (String usage, Peptide Presentation)
// Copyright   : Your copyright notice
// Description : Simulation of the effects of MHC peptide presentation
//			   : on KIR evolution
//============================================================================


#include "World.h"
#include "Contain.h"
#include <string.h>
using namespace std;

int main(int argc, char*argv[])
{


	//ACTUAL SIMULATION STARTING!
	cerr << "Testing if stderr is redirected to stdoutput" << endl;
	if (argc<3)
		{
		cerr << "Usage: "<<argv[0] << " <Parameter file> <Matrix file> <-b Loading Backup file> \n CoevolutionKIR simulates the evolution of the complex KIR system. For correct usage you have to indicate a parameter file, the affinity matrix file and a backup file name. If you are not loading any backup, just give any random name"<< endl;
		exit(-1);
		}

	cout <<"deleting old files ..."<<endl;


	string parameterFile(argv[1]);  //Parameters.data
	string matrixfile(argv[2]);    //JMMATRIXCM.csv or another type of matrix
	string backupFile;
	bool loadingBackup=false;
	for(int i=3; i< argc;i++)
		{
		if(strcmp(argv[i],"-b")==0 && argc > i)
			{
			loadingBackup = true;
			backupFile.assign(argv[++i]);
			}
		}

	World theWorld;
	theWorld.Loadmatrix(matrixfile);
	theWorld.LoadParameterFile(parameterFile);
	cout <<"welcome, hello" << endl;
	if(!loadingBackup)
		{
		// initialize host population
		theWorld.Initialize();
		}
	else
		{
		cout << "\n Loading data from backup file: "  << backupFile << endl;
		theWorld.LoadMap();
		theWorld.LoadBackupFile(backupFile);
		}
	theWorld.Simulate();
	cout << "bye bye \n";
//*/
	return 0;
}

/*
 *	12/6/2013 : ** RECODING **
 *
 */
