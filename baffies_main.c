/*
 *		Bechelor program, 100% from snapshot:
 *
 *	B inary
 *	A
 *	F raction
 *	F or
 *	I
 *	E
 *	S tar
 *	.
 *	C lusters
 *
 *	HARD-CODED:
 *		-> snapshot file path in fopen()
 *		-> COL and ROW, number of columns and rows in snapshot file, also the name of file containing snapshot data
 *		-> magnitudo limitation mag_lim, stars fainter than that ar considered as unseeable
 *		-> mass ratio limitation q_lim, border between single stars and binaries
 *		-> define content of sequent columns 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "baffies_fn.h"

#define STEPS 100

#define mag_lim 9
#define q_lim 0.5

#define rad 1
#define ikb 6
#define ik1 7
#define ik2 8
#define sm1 9
#define sm2 10
#define mv 21
#define mi 23
#define mv1 24
#define mi1 26
#define mv2 27
#define mi2 29

int COL;
int ROW; 


int main(){

	int c,r;

	double **snap;

	char model[30];
	char filename[90];

	FILE * fsnap;

new_model:

		printf("Choose model:\n");
		printf("\taccessible models:\n");
//		printf("\t---> model1\n");
//		printf("\t---> model2\n");

		printf("\t---> diogo-ibp-modified-1\n");
		printf("\t---> diogo-ibp-modified-2\n");
		printf("\t---> diogo-ibp-modified-3\n");
		printf("\t---> diogo-ibp-modified-4\n");
		printf("\t---> diogo-ibp-modified-5\n");
		printf("\t---> diogo-ibp-modified-6\n");

		printf("\t---> diogo-ibp-original-1\n");
		printf("\t---> diogo-ibp-original-2\n");
		printf("\t---> diogo-ibp-original-3\n");
		printf("\t---> diogo-ibp-original-4\n");
		printf("\t---> diogo-ibp-original-5\n");
		printf("\t---> diogo-ibp-original-6\n");

		printf("\n\tOverall Analysis (OLD)---> old\n");
		printf("\tOverall Analysis (NEW)---> new\n");
		printf("\tFcomparison (NEW)---> fcomp\n");
		printf("\tMagnitude Dependance (NEW)---> magdep\n");
		
		printf("\nstop program\t\t-> stop\n");

		scanf("%s",model);

//			if( 0 == strcmp(model,"model1")){COL = 36;ROW = 964007;} 
//			else if( 0 == strcmp(model,"model2")){COL = 36;ROW = 514064;} 
			if( 0 == strcmp(model,"diogo-ibp-modified-1")){COL = 30;ROW = 559221;} 
			else if( 0 == strcmp(model,"diogo-ibp-modified-2")){COL = 30;ROW = 242457;} 
			else if( 0 == strcmp(model,"diogo-ibp-modified-3")){COL = 30;ROW = 12691;} 
			else if( 0 == strcmp(model,"diogo-ibp-modified-4")){COL = 30;ROW = 1565455;} 
			else if( 0 == strcmp(model,"diogo-ibp-modified-5")){COL = 30;ROW = 683393;} 
			else if( 0 == strcmp(model,"diogo-ibp-modified-6")){COL = 30;ROW = 193538;} 

			else if( 0 == strcmp(model,"diogo-ibp-original-1")){COL = 30;ROW = 552773;} 
			else if( 0 == strcmp(model,"diogo-ibp-original-2")){COL = 30;ROW = 242971;} 
			else if( 0 == strcmp(model,"diogo-ibp-original-3")){COL = 30;ROW = 11956;} 
			else if( 0 == strcmp(model,"diogo-ibp-original-4")){COL = 30;ROW = 1547011;} 
			else if( 0 == strcmp(model,"diogo-ibp-original-5")){COL = 30;ROW = 710530;} 
			else if( 0 == strcmp(model,"diogo-ibp-original-6")){COL = 30;ROW = 212316;} 

			else if( 0 == strcmp(model,"old")){system("clear");UltimatePlot();goto new_model;} 
			else if( 0 == strcmp(model,"new")){system("clear");UltimatePlot_New();goto new_model;} 
			else if( 0 == strcmp(model,"fcomp")){system("clear");Ultimate_Fcomparison();goto new_model;}
			else if( 0 == strcmp(model,"magdep")){system("clear");Plot_MagnitudeDependance();goto new_model;} 
			else if( 0 == strcmp(model,"stop") ){goto end_program;}

			sprintf(filename,"data/%s/new_snap.dat",model);

		if(NULL==(snap=malloc(ROW*sizeof(double*)))){perror("malloc\n");exit(1);}
		for(r=0;r!=ROW;r++){
			if(NULL==(snap[r]=malloc(COL*sizeof(double)))){perror("malloc\n");exit(1);}
		}


			if(NULL==(fsnap=fopen(filename,"r"))){perror("fopen\n");exit(1);}

			for(r=0;r!=ROW;r++){
				for(c=0;c!=COL;c++){
					fscanf(fsnap,"%lf",&snap[r][c]);
				}
			}

//______________________PROGRAM:

		system("clear");
		char condition[4];
once_more:
		printf("Tell me what to do...\n");
		printf("binaries statistics\t-> stt\n");
		printf("binary fraction\t\t-> bf\n");
		printf("bin. frac. distribution\t-> bfd\n");
		printf("mass-ratio distribution\t-> qd\n");
		printf("radial distribution\t-> rdd\n");
		printf("distribution of binary systems\t-> snd\n");
		printf("plotting\t\t-> plot\t[work in progress]\n");
		printf("radii distro rest\t\t-> radi");
			
		printf("\n\nstop program\t\t-> stop\n");
		printf("change model\t\t-> cz\n");
		printf("already using model\t-> mdl\n");	
		scanf("%s",condition);

		if( 0 == strcmp(condition,"stt") ){system("clear");BinariesStatistics(snap);goto once_more;}
		else if( 0 == strcmp(condition,"stop") ){goto end_program;}
		else if( 0 == strcmp(condition,"cz") ){
			fclose(fsnap);
			for(c=0;c!=COL;c++){
				free(snap[c]);
			}
			free(snap);
			goto new_model;
		}
		else if( 0 == strcmp(condition,"mdl") ){system("clear");ShowModel(filename);goto once_more;}
		else if( 0 == strcmp(condition,"bf") ){system("clear");BinaryFraction(snap);goto once_more;}
		else if( 0 == strcmp(condition,"bfd") ){system("clear");BinaryFractionDistribution(snap,model);goto once_more;}
		else if( 0 == strcmp(condition,"snd") ){system("clear");SysNumberDistribution(snap,model);goto once_more;}
		else if( 0 == strcmp(condition,"qd") ){system("clear");MassRatioDistribution(snap,model);goto once_more;}
		else if( 0 == strcmp(condition,"rdd") ){system("clear");RadialDistribution(snap,model);goto once_more;}
		else if( 0 == strcmp(condition,"plot") ){system("clear");Plotting(snap,model,filename);goto once_more;}
		else if( 0 == strcmp(condition,"radi")){system("clear");RadialBinFracFunction_New(snap,model);}//New
		else{system("clear");printf("I don't get it. :c\n");goto once_more;}
		
end_program:
		fclose(fsnap);

	for(c=0;c!=COL;c++){
		free(snap[c]);
	}
	free(snap);

return EXIT_SUCCESS;
}
