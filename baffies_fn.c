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
#define XLIM 4	//half-mass radius limit during plotting

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

double delta_bin_frac[12];
int g = 0; 

double LinearInterpolation(double *x,double *y,int n,double x_p){

	double y_p;
	int i;
		for(i=0;i!=n-1;i++){
			if(x[i]<=x_p && x_p<=x[i+1]){
				y_p=(x_p-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i])+y[i];
			}
		}

	return y_p;
}

double CubicSpline(double *x_tab_in, double *y_tab_in, int n, double p){
        int i,j;
        double b,c,d,e,sum,temp;
        double *h,*s,*F,**m;
if((h=calloc(n,sizeof(double)))==NULL){
        perror("calloc failure c^:\n");
        exit(EXIT_FAILURE);
}
if((s=calloc(n,sizeof(double)))==NULL){
        perror("calloc failure c^:\n");
        exit(EXIT_FAILURE);
}
if((F=calloc(n,sizeof(double)))==NULL){
        perror("calloc failure c^:\n");
        exit(EXIT_FAILURE);
}
if((m=calloc(n,sizeof(double*)))==NULL){
        perror("calloc failure c^:\n");
        exit(EXIT_FAILURE);
}
        for(i=0;i!=n;i++){
                if((m[i]=calloc(n,sizeof(double)))==NULL){
                        perror("calloc failure c^:\n");
                        exit(EXIT_FAILURE);
                }
        }

        for(i=n-1;i>0;i--){
        F[i]=(y_tab_in[i]-y_tab_in[i-1])/(x_tab_in[i]-x_tab_in[i-1]);
        h[i-1]=x_tab_in[i]-x_tab_in[i-1];
        }
        for(i=1;i<n-1;i++){
        m[i][i]=2*(h[i-1]+h[i]);
                if(i!=1){
                m[i][i-1]=h[i-1];
                m[i-1][i]=h[i-1];
                }
        m[i][n-1]=6*(F[i+1]-F[i]);
        }
       for(i=1;i<n-2;i++){
        temp=(m[i+1][i]/m[i][i]);
                for(j=1;j<=n-1;j++){m[i+1][j]-=temp*m[i][j];}
        }
        for(i=n-2;i>0;i--){
        sum=0;
                for(j=i;j<=n-2;j++){
                sum+=m[i][j]*s[j];
                s[i]=(m[i][n-1]-sum)/m[i][i];
                }
        }
sum=-20;
        for(i=0;i<n-1;i++){
                if(x_tab_in[i]<=p && p<=x_tab_in[i+1]){
                                e=(s[i+1]-s[i])/(6*h[i]);
                                b=s[i]/2;
                                c=(y_tab_in[i+1]-y_tab_in[i])/h[i]-(2*h[i]*s[i]+s[i+1]*h[i])/6;
                                d=y_tab_in[i];
                                sum=e*pow((p-x_tab_in[i]),3)+b*pow((p-x_tab_in[i]),2)+c*(p-x_tab_in[i])+d;
				printf("sum=%lf\n",sum);
                }
        }
free(h);free(s);free(F);
for(i=0;i!=n;i++){
        free(m[i]);
}
free(m);

return sum;
}

void BinariesStatistics(double **snap){
//______________________Do statistics of objects in cluster
//			MS, WD, COMPACT, OTHER
//			Binaries:
//				MS	WD	COMPACT	OTHER
//			MS	x	x	x	x

			int i,r;
			int binaries[4]={0};
			double temp;

			for(r=0;r!=ROW;r++){
				if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]<2 ){ binaries[0]++; }  

				if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]>=11 && snap[r][ik2]<=13 ){ binaries[1]++; } 
				if( snap[r][ikb]!=0 && snap[r][ik2]<2 && snap[r][ik1]>=11 && snap[r][ik1]<=13 ){ binaries[1]++; } 

 				if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]>=13 && snap[r][ik2]<=14  ){ binaries[2]++; }  
 				if( snap[r][ikb]!=0 && snap[r][ik2]<2 && snap[r][ik1]>=13 && snap[r][ik1]<=14  ){ binaries[2]++; }  

				if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]>=2 && snap[r][ik2]<=9 ){ binaries[3]++; }  
				if( snap[r][ikb]!=0 && snap[r][ik2]<2 && snap[r][ik1]>=2 && snap[r][ik1]<=9 ){ binaries[3]++; }  
			}
			printf("Binaries, no limit:\n");
			printf("\tMS\tWD\tCOMPACT\tOTHER\t\tWD/(MS+WD)\n");
			temp=(double)binaries[1]/(binaries[1]+binaries[0]);
			printf("MS\t%i\t%i\t%i\t%i\t\t%lf\n",binaries[0],binaries[1],binaries[2],binaries[3],temp);
			
//			Number of binaries with magnitude limit:
			
			int lim;

			for(lim=13;lim!=3;lim--){		

				for(i=0;i!=4;i++){binaries[i]=0;}

				for(r=0;r!=ROW;r++){
					if( snap[r][mi]<=lim ){
						if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]<2 ){ binaries[0]++; }  
	
						if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]>=11 && snap[r][ik2]<=13 ){ binaries[1]++; } 
						if( snap[r][ikb]!=0 && snap[r][ik2]<2 && snap[r][ik1]>=11 && snap[r][ik1]<=13 ){ binaries[1]++; } 
		
		 				if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]>=13 && snap[r][ik2]<=14  ){ binaries[2]++; }  
		 				if( snap[r][ikb]!=0 && snap[r][ik2]<2 && snap[r][ik1]>=13 && snap[r][ik1]<=14  ){ binaries[2]++; }  
		
						if( snap[r][ikb]!=0 && snap[r][ik1]<2 && snap[r][ik2]>=2 && snap[r][ik2]<=9 ){ binaries[3]++; }  
						if( snap[r][ikb]!=0 && snap[r][ik2]<2 && snap[r][ik1]>=2 && snap[r][ik1]<=9 ){ binaries[3]++; }  
					}			
				}
			printf("%i mag limit:\n",lim);
			printf("\tMS\tWD\tCOMPACT\tOTHER\t\tWD/(MS+WD)\n");
			temp=(double)binaries[1]/(binaries[1]+binaries[0]);
			printf("MS\t%i\t%i\t%i\t%i\t\t%lf\n",binaries[0],binaries[1],binaries[2],binaries[3],temp);
			}
}

void ObservationalBinaryFraction(double **snap){

	int r;
	int N_single=0;
	int N_binary=0;
	double bin_frac;
printf("\n\tObservational binary fraction:\n");
//______________SINGLES:
		for(r=0;r!=ROW;r++){
			if( mag_lim >= snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_single++;}
			}
		}
		printf("\tSingle stars:\t%i\n",N_single);
		
//______________BINARY:
		for(r=0;r!=ROW;r++){
			if( mag_lim >= snap[r][mi] ){
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_binary++;}
			}
		}
		printf("\tBinaries:\t%i\n",N_binary);
		
		bin_frac=(double)N_binary/(N_binary+N_single);
		printf("\tBinary fraction:\t%lf\n",bin_frac);
		printf("\n");
}

double RealBinaryFraction(double **snap){

	int r;
	
	int N_ms = 0;
	int N_msms_sq = 0;
	int N_msms_lq = 0;
	int N_msc = 0;
	int N_all;

	double bin_frac;
//	double limit;

//	printf("Type magnitude limit (default 9 mag):");
//	scanf("%lf",&limit);

	printf("\n\tReal binary fraction:\n");

//______________SINGLES:
		for(r=0;r!=ROW;r++){	
//				if(  limit >= snap[r][mi] && 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single++;}
			if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms++;}
		}
//		printf("\tSingle stars:\t%i\n",N_single);
		
//______________BINARY:
		for(r=0;r!=ROW;r++){
//				if( limit >= snap[r][mi] && 1 == snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] ){N_binary++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_msms_lq++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc++;}
		}
//		printf("\tBinaries:\t%i\n",N_binary);

		N_all = N_ms + N_msms_sq + N_msms_lq + N_msc;
		
		bin_frac=(double)(N_msms_sq + N_msms_lq + N_msc)/(double)N_all;
		printf("\tBinary fraction:\t%lf\n",bin_frac);
		printf("\n");

		return bin_frac;
}

double TrueBinaryFraction(double **snap){

	int r;
	
	int sin = 0;
	int bin = 0;

	double bin_frac;

	printf("\n\tTrue binary fraction:\n");

		for(r=0;r!=ROW;r++){	
			if( 0 == snap[r][ikb] ){sin++;}	//Single
			if( 0 != snap[r][ikb] ){bin++;}	//Binary
		}		
		
		bin_frac=(double)(bin)/(double)(sin+bin);
		printf("\tBinary fraction:\t%lf\n",bin_frac);
		printf("\n");

		return bin_frac;
}

void BinaryFraction(double **snap){

	char string[4];
binfrac_more:
	printf("What would you like to calculate?\n");
	printf("Observational binary fraction\t-> obs\n");
	printf("Real binary fraction\t\t-> real\n");
	printf("Both\t\t\t\t-> both\n");

	printf("\nReturn\t-> R\n");

	scanf("%s",string);

		if( 0 == strcmp(string,"obs") ){ObservationalBinaryFraction(snap);goto binfrac_more;}
		else if( 0 == strcmp(string,"real") ){RealBinaryFraction(snap);goto binfrac_more;}
		else if( 0 == strcmp(string,"both") ){ObservationalBinaryFraction(snap);RealBinaryFraction(snap);goto binfrac_more;}
		else if( 0 == strcmp(string,"R") ){system("clear");}
		else{system("clear");printf("I don't get it. :c\n");goto binfrac_more;}

}
//______________________________________________________________New magnitude distribution for 3 bin_frac parameters
double BinaryFractionDistribution_New(double **snap, char *model){

	printf("\n\t5 intervals in magnitude range [4:9] mag\n");
	printf("BinFracObs  = N_ms-ms_q>0.5/N_all\n");
	printf("BinFracMS   = N_ms-ms_all/N_all\n");
	printf("BinFracReal = N_ms-ms+N_ms-compact/N_all\n");
	printf("\tInterval[mag]\tBinFracObs\tBinFracMS\tBinFracReal\n");

	int i,r;
	int N_ms;		//single MS stars
	int N_msms_sq;		//systems MS-MS with small q (q < q_limit)
	int N_msms_lq;		//systems MS-MS with large q (q > q_limit)
	int N_msc;		//systems of MS stars with compact objects (WD, NS, BH)
	int N_all;

	int N_msms_lq_47 = 0, N_all_47 = 0; //for [4:7] magnitude

	double bin_frac_o,bin_frac_r,bin_frac_ms;
	double bin_frac_o_47 = 0;	//Observational binary fraction in magnitude [4:7]

	FILE *out;

	char filename[60];
	sprintf(filename,"bfm_new_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("bfm_new fopen failure");exit(1);}

	for(i=4;i!=10;i++){

	N_ms = 0;
	N_msms_sq = 0;
	N_msms_lq = 0;
	N_msc = 0;
	N_all = 0;

		for(r=0;r!=ROW;r++){
		if( 9 != i ){
				if( (double)i <= snap[r][mi] &&  (double)i+1 > snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_msms_lq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc++;}
				}
		}else{
				if( (double)i <= snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_msms_lq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc++;}
				}
		}

		}

		N_all = N_ms + N_msms_sq + N_msms_lq + N_msc;
//		printf("N_ms=%i\nN_msms_sq=%i\nN_msms_lq=%i\nN_msc=%i\nN_all=%i\n",N_ms,N_msms_sq,N_msms_lq,N_msc,N_all);
		bin_frac_o  = (double)N_msms_lq/(double)N_all;
		bin_frac_ms = (double)(N_msms_lq + N_msms_sq)/(double)N_all;
		bin_frac_r  = (double)(N_msms_lq + N_msms_sq + N_msc)/(double)N_all;
		
		if( 9 != i ){
			printf("\t[%i,%i]\t\t%lf\t%lf\t%lf\n",i,i+1,bin_frac_o,bin_frac_ms,bin_frac_r);
			fprintf(out,"[%i,%i]\t\t%lf\t%lf\t%lf\n",i,i+1,bin_frac_o,bin_frac_ms,bin_frac_r);
		}else{
			printf("\t[%i,max]\t\t%lf\t%lf\t%lf\n",i,bin_frac_o,bin_frac_ms,bin_frac_r);
			fprintf(out,"[%i,max]\t\t%lf\t%lf\t%lf\n",i,bin_frac_o,bin_frac_ms,bin_frac_r);
		
		}
		
		if( 4 <= i && 7 > i ){
			N_msms_lq_47 += N_msms_lq;
			N_all_47 += N_all;
		}

	}
	
	bin_frac_o_47 = (double)N_msms_lq_47/(double)N_all_47;

	fclose(out);
	
	return bin_frac_o_47; //returns observational binary fraction for [4:7] magnitudes
}
//Same as BinaryFractionDistribution_New but produce one column more to check magnitude dependance
void ObservationalMagnitudeDependance(double **snap, char *model){

	printf("\n\t5 intervals in magnitude range [4:9] mag\n");
	printf("BinFracObs  = N_ms-ms_q>0.5/N_all\n");
	printf("BinFracMS   = N_ms-ms_all/N_all\n");
	printf("BinFracReal = N_ms-ms+N_ms-compact/N_all\n");
	printf("\tInterval[mag]\tBinFracObs\tfobs/ftot\n");

	int i,r;
	int N_ms;		//single MS stars
	int N_msms_sq;		//systems MS-MS with small q (q < q_limit)
	int N_msms_lq;		//systems MS-MS with large q (q > q_limit)
	int N_msc;		//systems of MS stars with compact objects (WD, NS, BH)
	int N_all;

	double bin_frac_o,magdep;
	
	FILE *out;

	char filename[60];
	sprintf(filename,"magdep_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("magdep fopen failure");exit(1);}

	for(i=4;i!=10;i++){

	N_ms = 0;
	N_msms_sq = 0;
	N_msms_lq = 0;
	N_msc = 0;
	N_all = 0;

	double bin_frac_true = TrueBinaryFraction(snap);

		for(r=0;r!=ROW;r++){
		if( 9 != i ){
				if( (double)i <= snap[r][mi] &&  (double)i+1 > snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_msms_lq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc++;}
				}
		}else{
				if( (double)i <= snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_msms_lq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc++;}
				}
		}

		}

		N_all = N_ms + N_msms_sq + N_msms_lq + N_msc;

		bin_frac_o  = (double)N_msms_lq/(double)N_all;
		magdep	= bin_frac_o/bin_frac_true;
		
		if( 9 != i ){
			printf("\t[%i,%i]\t%lf\n",i,i+1,magdep);
			fprintf(out,"[%i,%i]\t%lf\n",i,i+1,magdep);
		}else{
			printf("\t[%i,max]\t%lf\n",i,magdep);
			fprintf(out,"[%i,max]\t%lf\n",i,magdep);
		
		}
	}

	fclose(out);
}

void BinaryFractionDistribution(double **snap, char *model){
	
	printf("\n\t9 intervals in magnitude range [4:13] mag\n");
	printf("\tInterval[mag]\tBin Frac[OBS]\tBin Frac[REAL]\n");

	int i,r;
	int N_single;
	int N_binary;

	double bin_frac_o,bin_frac_r;

	FILE *out;

	char filename[60];
	sprintf(filename,"bfm_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("bfm fopen failure");exit(1);}

	for(i=4;i!=13;i++){
	N_single=0;
	N_binary=0;

//______________SINGLES:
		for(r=0;r!=ROW;r++){
			if( (double)i <= snap[r][mi] &&  (double)i+1 > snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_single++;}
			}
		}		
//______________BINARY:
		for(r=0;r!=ROW;r++){
			if( (double)i <= snap[r][mi] &&  (double)i+1 > snap[r][mi] ){
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_binary++;}
			}
		}
		
		bin_frac_o=(double)N_binary/(N_binary+N_single);

	N_single=0;
	N_binary=0;

//______________SINGLES:
		for(r=0;r!=ROW;r++){	
			if( (double)i <= snap[r][mi] &&  (double)i+1 > snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single++;}
			}
		}
		
//______________BINARY:
		for(r=0;r!=ROW;r++){
			if( (double)i <= snap[r][mi] &&  (double)i+1 > snap[r][mi] ){
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] ){N_binary++;}
			}
		}

		bin_frac_r=(double)N_binary/(N_binary+N_single);

		printf("\t[%i,%i]\t\t%lf\t%lf\n",i,i+1,bin_frac_o,bin_frac_r);
		fprintf(out,"[%i,%i]\t\t%lf\t%lf\n",i,i+1,bin_frac_o,bin_frac_r);

	}
	fclose(out);
}

double MassRatioDistribution(double **snap, char *model){

	int i,r;

	int N_binary = 0;
	int N_total = 0;
	int INT = 10; //number of intervals

	double bin_frac_r;
	double bin_frac_tot=0;
	double q_a,q_b;

	FILE *out;

	char filename[60];
	sprintf(filename,"bfq_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("bfq fopen failure");exit(1);}


	printf("\n\tI provide values of real binary fraction in certain mass-ratio intervals\n");

	printf("\tInterval\tv\tN_binaries\n");

//______________SINGLES:
		for(r=0;r!=ROW;r++){	
			if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_total++;}
		}	
//______________BINARY:
		for(r=0;r!=ROW;r++){
			if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] ){N_total++;}
		}
		printf("\tN_total=%i\n",N_total);
	for(i=0;i!=INT;i++){

		q_a=(double)i/INT;
		q_b=((double)i+1)/INT;	

		N_binary=0;

//______________BINARY:
		for(r=0;r!=ROW;r++){
			if( q_b >= snap[r][sm2]/snap[r][sm1] && q_a < snap[r][sm2]/snap[r][sm1] ){
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] ){N_binary++;}
			}
		}
		bin_frac_r=(double)N_binary/N_total;
		bin_frac_tot+=bin_frac_r;
		printf("\t[%.1lf,%.1lf]\t%lf\t%i\n",q_a,q_b,bin_frac_r,N_binary);
		fprintf(out,"[%.1lf,%.1lf]\t%lf\t%i\n",q_a,q_b,bin_frac_r,N_binary);
	}
	printf("\n\tbin_frac_total=%lf\n",bin_frac_tot);
fclose(out);

return bin_frac_tot;
}

double MaxRadius(double **snap){
	
	int r;

	double r_MAX=0;
	double r_TMP=0;

		r_MAX=snap[0][rad];

		for(r=1;r!=ROW;r++){

			r_TMP=snap[r][rad];

			if(r_TMP > r_MAX){
				r_MAX = r_TMP;
			}	
		}
//	printf("r_max=%lf\n",r_MAX);
	return r_MAX;
}

double TotalMass(double **snap, double from, double to,double r_MAX){

//	{from, to} - in percents of max radius

	int r;
	
	double m_TOT=0;
	double r_TMP;

		for(r=0;r!=ROW;r++){

			r_TMP=100*snap[r][rad]/r_MAX; //in percents of max radius

			if(to == 100){
				if(r_TMP>=from && r_TMP<=to){
					m_TOT+=snap[r][sm1];
					m_TOT+=snap[r][sm2];	
				}
			}else{
				if(r_TMP>=from && r_TMP<to){
					m_TOT+=snap[r][sm1];
					m_TOT+=snap[r][sm2];	
				}		
			}
		}
	return m_TOT; 
}
double **returnRadialMassFunction(double **snap,int INTERVALS){
	
	int i;

	double m_TOT;
	double r_MAX;

	double **mass_fun;

	//	mass_fun[r][total_mass]

	if(NULL == (mass_fun = malloc((INTERVALS+1)*sizeof(double*)))){perror("malloc\n");exit(1);}
	for(i=0;i!=INTERVALS+1;i++){
		if(NULL == (mass_fun[i] = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	}

		r_MAX = MaxRadius(snap);
		m_TOT = TotalMass(snap,0,100,r_MAX);
		printf("r_MAX=%lf\n",r_MAX);

		double r_TMP;

		for(i=1;i!=INTERVALS+1;i++){
			r_TMP = (double)(i)/(double)INTERVALS*r_MAX;	//same units as in snapshot

			if(1 == i){
				mass_fun[i-1][0] = 0;
				mass_fun[i-1][1] = 0;
				
				mass_fun[i][0] = r_TMP/r_MAX;
				mass_fun[i][1] = TotalMass(snap,(i-1)*100/(double)INTERVALS,i*100/(double)INTERVALS,r_MAX)/m_TOT;							
			}else{
				mass_fun[i][0] = r_TMP/r_MAX;
				mass_fun[i][1] = TotalMass(snap,(i-1)*100/(double)INTERVALS,i*100/(double)INTERVALS,r_MAX)/m_TOT+mass_fun[i-1][1];				
			}
		}
	return mass_fun;
}

void RadialMassFunction(double **snap){
	
	int i;
	int INTERVALS;

	double m_TOT;
	double r_MAX;

	double **mass_fun;

	//	mass_fun[r][total_mass]

	printf("How many intervals should I consider?\n");
	scanf("%i",&INTERVALS);

	if(NULL == (mass_fun = malloc((INTERVALS+1)*sizeof(double*)))){perror("malloc\n");exit(1);}
	for(i=0;i!=INTERVALS+1;i++){
		if(NULL == (mass_fun[i] = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	}

		r_MAX = MaxRadius(snap);
		m_TOT = TotalMass(snap,0,100,r_MAX);
		printf("r_MAX=%lf\n",r_MAX);

		double r_TMP;

		for(i=1;i!=INTERVALS+1;i++){
			r_TMP = (double)(i)/(double)INTERVALS*r_MAX;	//same units as in snapshot

			if(1 == i){
				mass_fun[i-1][0] = 0;
				mass_fun[i-1][1] = 0;
				
				mass_fun[i][0] = r_TMP/r_MAX;
				mass_fun[i][1] = TotalMass(snap,(i-1)*100/(double)INTERVALS,i*100/(double)INTERVALS,r_MAX)/m_TOT;							
			}else{
				mass_fun[i][0] = r_TMP/r_MAX;
				mass_fun[i][1] = TotalMass(snap,(i-1)*100/(double)INTERVALS,i*100/(double)INTERVALS,r_MAX)/m_TOT+mass_fun[i-1][1];				
			}
		}

	printf("\tr/r_MAX\t\tm/m_TOT\n");
	for(i=0;i!=INTERVALS+1;i++){
		printf("\t%lf\t%lf\n",mass_fun[i][0],mass_fun[i][1]);
	}

	for(i=0;i!=INTERVALS+1;i++){
		free(mass_fun[i]);
	}
	free(mass_fun);
}

void RadialBinFracFunction_New(double **snap, char *model){

	int i,r;
	int N_ms;
	int N_msms_sq;
	int N_msms_lq;
	int N_msc;
	int N_all;
//Parameters below only for bin_frac_r. Case without magnitude limits.
	int N_ms_r;
	int N_msms_sq_r;
	int N_msms_lq_r;
	int N_msc_r;
	int N_all_r;

	double bin_frac_o, bin_frac_r, bin_frac_ms;

	double **mass_fun;
	if(NULL == (mass_fun = malloc((STEPS)*sizeof(double*)))){perror("malloc\n");exit(1);}
	for(i=0;i!=STEPS;i++){
		if(NULL == (mass_fun[i] = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	}

	double *radi,*mass;
	if(NULL == (radi=malloc(STEPS*sizeof(double)))){perror("RadialRealBinFracFunction malloc fail");exit(1);}
	if(NULL == (mass=malloc(STEPS*sizeof(double)))){perror("RadialRealBinFracFunction malloc fail");exit(1);}

	double r_MAX;
	double r_TMP;
	double L[6];
/*	L[0] = 0
 *	L[1] = L10%
 *	L[2] = L30%
 *	L[3] = L50%
 *	L[4] = L70%
 *	L[5] = L100%
 */
	r_MAX=MaxRadius(snap);

	FILE *fmassfun;
	fmassfun=fopen("fmassfun.dat","w");

	mass_fun = returnRadialMassFunction(snap,STEPS);
	for(i=0;i!=STEPS;i++){
		radi[i]=mass_fun[i][0];
		mass[i]=mass_fun[i][1];
		fprintf(fmassfun,"%lf\t%lf\n",radi[i],mass[i]);
	}
	fclose(fmassfun);

	FILE *out;

	char filename[60];
	sprintf(filename,"bfr_new_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("bfr fopen failure");exit(1);}

	L[0]=0;						printf("\trL0=%lf\n",L[0]);
	L[1]=LinearInterpolation(mass,radi,STEPS,0.1);	printf("\trL10=%lf\n",L[1]);
	L[2]=LinearInterpolation(mass,radi,STEPS,0.3);	printf("\trL30=%lf\n",L[2]);
	L[3]=LinearInterpolation(mass,radi,STEPS,0.5);	printf("\trL50=%lf\n",L[3]);
	L[4]=LinearInterpolation(mass,radi,STEPS,0.7);	printf("\trL70=%lf\n",L[4]);
	L[5]=1;						printf("\trL100=%lf\n",L[5]);

	printf("\n\tBinary fraction:\n");

	printf("Magnitude range [4-7] mag (only for bin_frac_ms and bin_frac_o)\n");
	printf("\tLagrange Radius\tbin_frac_o\tbin_frac_ms\tbin_frac_r\n");

	for(i=0;i!=5;i++){

	N_ms = N_msms_sq = N_msms_lq = N_msc = N_all = N_ms_r = N_msms_sq_r = N_msms_lq_r = N_msc_r = N_all_r = 0;
	
	if( 0 == i )	{fprintf(out,"0-10\t");printf("0-10\t");}
	else if( 1 == i){fprintf(out,"\n10-30\t");printf("\n10-30\t");}
	else if( 2 == i){fprintf(out,"\n30-50\t");printf("\n30-50\t");}
	else if( 3 == i){fprintf(out,"\n50-70\t");printf("\n50-70\t");}
	else if( 4 == i){fprintf(out,"\n70-100\t");printf("\n70-100\t");}

		for(r=0;r!=ROW;r++){	

			//As you can see, in that case there are no magnitude limits. It's only for bin_frac_r parameter. 

			r_TMP=snap[r][rad]/r_MAX; //as a part of max radius

			if(r_TMP>L[i] && r_TMP<=L[i+1]){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms_r++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq_r++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_msms_lq_r++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc_r++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc_r++;}
				}
		
			if( 4 <= snap[r][mi] && 7 >= snap[r][mi] ){ //Arbitrary limits between 4 mag and 7 mag for bin_frac_o and bin_frac_ms

			r_TMP=snap[r][rad]/r_MAX; //as a part of max radius

			if(r_TMP>L[i] && r_TMP<=L[i+1]){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_msms_lq++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc++;}
				}
			}
		}
		N_all = N_ms + N_msms_sq + N_msms_lq + N_msc;
		N_all_r = N_ms_r + N_msms_sq_r + N_msms_lq_r + N_msc_r;
//		printf("N_ms=%i\nN_msms_sq=%i\nN_msms_lq=%i\nN_msc=%i\nN_all=%i\n",N_ms,N_msms_sq,N_msms_lq,N_msc,N_all);
		bin_frac_o  = (double)N_msms_lq/(double)N_all;
		bin_frac_ms = (double)(N_msms_lq + N_msms_sq)/(double)N_all;

		bin_frac_r  = (double)(N_msms_lq_r + N_msms_sq_r + N_msc_r)/(double)N_all_r;

		printf("%lf\t%lf\t%lf",bin_frac_o,bin_frac_ms,bin_frac_r);
		fprintf(out,"%lf\t%lf\t%lf",bin_frac_o,bin_frac_ms,bin_frac_r);
	}

	fclose(out);
}

double RadialRealBinFracFunction(double **snap, char *model){

	int i,r;
	int N_single=0;
	int N_binary=0;
	int N_binary_i=0;
	int N_single_i=0;

	double **mass_fun;
	if(NULL == (mass_fun = malloc((STEPS)*sizeof(double*)))){perror("malloc\n");exit(1);}
	for(i=0;i!=STEPS;i++){
		if(NULL == (mass_fun[i] = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	}

	double *radi,*mass;
	if(NULL == (radi=malloc(STEPS*sizeof(double)))){perror("RadialRealBinFracFunction malloc fail");exit(1);}
	if(NULL == (mass=malloc(STEPS*sizeof(double)))){perror("RadialRealBinFracFunction malloc fail");exit(1);}

	double bin_frac[5],bin_frac_i[5];
	double limit;
	double r_MAX;
	double r_TMP;
	double L[6];
/*	L[0] = 0
 *	L[1] = L10%
 *	L[2] = L30%
 *	L[3] = L50%
 *	L[4] = L70%
 *	L[5] = L100%
 */
	r_MAX=MaxRadius(snap);

	FILE *fmassfun;
	fmassfun=fopen("fmassfun.dat","w");

	mass_fun = returnRadialMassFunction(snap,STEPS);
	for(i=0;i!=STEPS;i++){
		radi[i]=mass_fun[i][0];
		mass[i]=mass_fun[i][1];
		fprintf(fmassfun,"%lf\t%lf\n",radi[i],mass[i]);
	}
	fclose(fmassfun);

	FILE *out;

	char filename[60];
	sprintf(filename,"bfr_r_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("real bfr fopen failure");exit(1);}

	L[0]=0;						printf("\trL0=%lf\n",L[0]);
	L[1]=LinearInterpolation(mass,radi,STEPS,0.1);	printf("\trL10=%lf\n",L[1]);
	L[2]=LinearInterpolation(mass,radi,STEPS,0.3);	printf("\trL30=%lf\n",L[2]);
	L[3]=LinearInterpolation(mass,radi,STEPS,0.5);	printf("\trL50=%lf\n",L[3]);
	L[4]=LinearInterpolation(mass,radi,STEPS,0.7);	printf("\trL70=%lf\n",L[4]);
	L[5]=1;						printf("\trL100=%lf\n",L[5]);

	printf("\n\tReal binary fraction:\n");

	printf("Type magnitude limit (default 9 mag):");
	scanf("%lf",&limit);
//______________SINGLES:
		for(r=0;r!=ROW;r++){	
			if(  limit >= snap[r][mi] && 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single++;}
		}
printf("Singles:%i\n",N_single);
//______________BINARY:
		for(r=0;r!=ROW;r++){
			if( limit >= snap[r][mi] && 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] ){N_binary++;}
		}
printf("Binaries:%i\n",N_binary);

	for(i=0;i!=5;i++){
	
	if( 0 ==i )	{fprintf(out,"0-10\t");}
	else if( 1 == i){fprintf(out,"\n10-30\t");}
	else if( 2 == i){fprintf(out,"\n30-50\t");}
	else if( 3 == i){fprintf(out,"\n50-70\t");}
	else if( 4 == i){fprintf(out,"\n70-100\t");}

	N_binary_i=0;
	N_single_i=0;
	
//______________SINGLES:
		for(r=0;r!=ROW;r++){	

			r_TMP=snap[r][rad]/r_MAX; //as a part of max radius

			if(r_TMP>L[i] && r_TMP<=L[i+1]){

				if(  limit >= snap[r][mi] && 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single_i++;}
			}
		}

//______________BINARY:
		for(r=0;r!=ROW;r++){

			r_TMP=snap[r][rad]/r_MAX; //as a part of max radius

			if(r_TMP>L[i] && r_TMP<=L[i+1]){
				if( limit >= snap[r][mi] && 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] ){N_binary_i++;}
			}

		}
printf("--->Binaries_%i=%i\n",i,N_binary_i);

		bin_frac[i]=(double)N_binary_i/(N_binary+N_single);
		bin_frac_i[i]=(double)N_binary_i/(N_binary_i+N_single_i);

		fprintf(out,"%lf\t%i",bin_frac_i[i],N_binary_i);

	}
	printf("\tLagrange Radius\tBi/(B+S)\tBi/(Bi+Si)\n");

	printf("\t0-10\t\t%lf\t%lf\n",bin_frac[0],bin_frac_i[0]);
	printf("\t10-30\t\t%lf\t%lf\n",bin_frac[1],bin_frac_i[1]);
	printf("\t30-50\t\t%lf\t%lf\n",bin_frac[2],bin_frac_i[2]);
	printf("\t50-70\t\t%lf\t%lf\n",bin_frac[3],bin_frac_i[3]);
	printf("\t70-100\t\t%lf\t%lf\n",bin_frac[4],bin_frac_i[4]);
	printf("\n");

	fclose(out);

	double sum=0;
	for(i=0;i!=5;i++){
		sum+=bin_frac[i];
	}
	printf("\tSummary Bin. Frac.=%lf\n\n",sum);

return limit;
}

double RadialObsBinFracFunction(double **snap, char *model){
	
	int i,r;
	int N_single=0;
	int N_binary=0;
	int N_binary_i=0;
	int N_single_i=0;

	double **mass_fun;
	if(NULL == (mass_fun = malloc((STEPS)*sizeof(double*)))){perror("malloc\n");exit(1);}
	for(i=0;i!=STEPS;i++){
		if(NULL == (mass_fun[i] = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	}

	double *radi,*mass;
	if(NULL == (radi=malloc(STEPS*sizeof(double)))){perror("RadialObsBinFracFunction malloc fail");exit(1);}
	if(NULL == (mass=malloc(STEPS*sizeof(double)))){perror("RadialObsBinFracFunction malloc fail");exit(1);}

	double bin_frac[5],bin_frac_i[5];
	double limit;
	double r_MAX;
	double r_TMP;
	double L[6];
/*	L[0] = 0
 *	L[1] = L10%
 *	L[2] = L30%
 *	L[3] = L50%
 *	L[4] = L70%
 *	L[5] = L100%
 */
	r_MAX=MaxRadius(snap);

	FILE *fmassfun;
	fmassfun=fopen("fmassfun.dat","w");

	mass_fun = returnRadialMassFunction(snap,STEPS);
	for(i=0;i!=STEPS;i++){
		radi[i]=mass_fun[i][0];
		mass[i]=mass_fun[i][1];
		fprintf(fmassfun,"%lf\t%lf\n",radi[i],mass[i]);
	}
	fclose(fmassfun);

	FILE *out;

	char filename[60];
	sprintf(filename,"bfr_o_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("obs bfr fopen failure");exit(1);}

	L[0]=0;						printf("\trL0=%lf\n",L[0]);
	L[1]=LinearInterpolation(mass,radi,STEPS,0.1);	printf("\trL10=%lf\n",L[1]);
	L[2]=LinearInterpolation(mass,radi,STEPS,0.3);	printf("\trL30=%lf\n",L[2]);
	L[3]=LinearInterpolation(mass,radi,STEPS,0.5);	printf("\trL50=%lf\n",L[3]);
	L[4]=LinearInterpolation(mass,radi,STEPS,0.7);	printf("\trL70=%lf\n",L[4]);
	L[5]=1;						printf("\trL100=%lf\n",L[5]);

	printf("\n\tObservational binary fraction:\n");

	printf("Type magnitude limit (default 9 mag):");
	scanf("%lf",&limit);
//______________SINGLES:
		for(r=0;r!=ROW;r++){
			if(limit >= snap[r][mi]){	
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_single++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_single++;}
			}
		}
printf("Singles:%i\n",N_single);
//______________BINARY:
		for(r=0;r!=ROW;r++){
			if( limit >= snap[r][mi] && 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_binary++;}
		}
printf("Binaries:%i\n",N_binary);

	for(i=0;i!=5;i++){

	if( 0 ==i )	{fprintf(out,"0-10\t");}
	else if( 1 == i){fprintf(out,"\n10-30\t");}
	else if( 2 == i){fprintf(out,"\n30-50\t");}
	else if( 3 == i){fprintf(out,"\n50-70\t");}
	else if( 4 == i){fprintf(out,"\n70-100\t");}

	N_binary_i=0;
	N_single_i=0;
//______________SINGLES:
		for(r=0;r!=ROW;r++){	

			r_TMP=snap[r][rad]/r_MAX; //as a part of max radius

			if(r_TMP>L[i] && r_TMP<=L[i+1]){
				if(limit >= snap[r][mi]){	
					if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_single_i++;}
					if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_single_i++;}
					if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_single_i++;}
					if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_single_i++;}
				}
			}
		}

//______________BINARY:
		for(r=0;r!=ROW;r++){

			r_TMP=snap[r][rad]/r_MAX; //as a part of max radius

			if(r_TMP>L[i] && r_TMP<=L[i+1]){
		if( limit >= snap[r][mi] && 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim < snap[r][sm2]/snap[r][sm1] ){N_binary_i++;}
			}

		}
printf("--->Binaries_%i=%i\n",i,N_binary_i);

		bin_frac[i]=(double)N_binary_i/(N_binary+N_single);
		bin_frac_i[i]=(double)N_binary_i/(N_binary_i+N_single_i);

		fprintf(out,"%lf\t%i",bin_frac_i[i],N_binary_i);

	}
	fclose(out);

	printf("\tLagrange Radius\tBi/(B+S)\tBi/(Bi+Si)\n");

	printf("\t0-10\t\t%lf\t%lf\n",bin_frac[0],bin_frac_i[0]);
	printf("\t10-30\t\t%lf\t%lf\n",bin_frac[1],bin_frac_i[1]);
	printf("\t30-50\t\t%lf\t%lf\n",bin_frac[2],bin_frac_i[2]);
	printf("\t50-70\t\t%lf\t%lf\n",bin_frac[3],bin_frac_i[3]);
	printf("\t70-100\t\t%lf\t%lf\n",bin_frac[4],bin_frac_i[4]);
	printf("\n");
	double sum=0;
	for(i=0;i!=5;i++){
		sum+=bin_frac[i];
	}
	printf("\tSummary Bin. Frac.=%lf\n\n",sum);

return limit;
}

void RadialBinFracFunction(double **snap, char *model){
	
	char string[4];
bf_radii_more:
	printf("What would you like to calculate?\n");
	printf("Observational binary fraction\t-> obs\n");
	printf("Real binary fraction\t\t-> real\n");
	printf("Both\t\t\t\t-> both\n");

	printf("\nReturn\t-> R\n");

	scanf("%s",string);

		if( 0 == strcmp(string,"obs") ){RadialObsBinFracFunction(snap,model);goto bf_radii_more;}
		else if( 0 == strcmp(string,"real") ){RadialRealBinFracFunction(snap,model);goto bf_radii_more;}
		else if( 0 == strcmp(string,"both") ){RadialObsBinFracFunction(snap,model);RadialRealBinFracFunction(snap,model);goto bf_radii_more;}
		else if( 0 == strcmp(string,"R") ){system("clear");}
		else{system("clear");printf("I don't get it. :c\n");goto bf_radii_more;}
}

void RadialDistribution(double **snap, char *model){

	char string[4];
radii_more:
	printf("Radial Mass Distribution\t-> rmd\n");
	printf("Binary Fraction\t\t\t-> bf\n");

	printf("\nReturn\t-> R\n");

	scanf("%s",string);

	if( 0 == strcmp(string,"rmd") ){system("clear");RadialMassFunction(snap);goto radii_more;}
	else if( 0 == strcmp(string,"bf") ){system("clear");RadialBinFracFunction(snap,model);goto radii_more;}
	else if( 0 == strcmp(string,"R") ){system("clear");}
	else{system("clear");printf("I don't get it. :c\n");goto radii_more;}

}

int ** SysNrRadialDistribution_New(double **snap, char *model){

	int steps = 100;

	int i,r,m;
	double r_av,r_tmp,r_max;

	int N_ms[5] = {0};		//MS single {for 5 magnitude intervals}
	int N_msms_sq[5] = {0};		//MS-MS binary with q < 0.5 (small)
	int N_msms_lq[5] = {0};		// q > 0.5 (large)
	int N_msc[5] = {0};		// MS - compact
	int **N_tot;			// Total number of objects considered in terms of mag bins and for certain type of objects. 
					//
					//	________|_4:5_|_5:6_|_6:7_|_7:8_|_8:9_|___
					//	ms	|
					//	msms_sq	|	4 x 5 = 20
					//	msms_lq	|	Matrix N_tot[4][5]
					//	msc	|
	if(NULL == (N_tot = calloc(4,sizeof(int*)))){perror("calloc\n");exit(1);}
	for(i=0;i!=4;i++){
		if(NULL == (N_tot[i] = calloc(5,sizeof(int)))){perror("calloc\n");exit(1);}
	}

	FILE *out;
	
	char filename[60];
	sprintf(filename,"bnr_new_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("bnr fopen failure");exit(1);}
//MASS-RADIUS
	double **mass_fun;
	if(NULL == (mass_fun = malloc((STEPS)*sizeof(double*)))){perror("malloc\n");exit(1);}
	for(i=0;i!=STEPS;i++){
		if(NULL == (mass_fun[i] = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	}

	double *radi,*mass;
	if(NULL == (radi=malloc(STEPS*sizeof(double)))){perror("RadialObsBinFracFunction malloc fail");exit(1);}
	if(NULL == (mass=malloc(STEPS*sizeof(double)))){perror("RadialObsBinFracFunction malloc fail");exit(1);}

	double r_L50;

	r_max=MaxRadius(snap);

	FILE *fmassfun;
	fmassfun=fopen("fmassfun.dat","w");

	mass_fun = returnRadialMassFunction(snap,STEPS);
	for(i=0;i!=STEPS;i++){
		radi[i]=mass_fun[i][0];
		mass[i]=mass_fun[i][1];
		fprintf(fmassfun,"%lf\t%lf\n",radi[i],mass[i]);
	}
	fclose(fmassfun);
	
	//Half-mass radius
	r_L50 = LinearInterpolation(mass,radi,STEPS,0.5);
	printf("\trL50=%lf\n",r_L50*r_max);
/*
	for(i=0;i!=5;i++){

	if( 0 ==i )	{fprintf(out,"0-10\t");}
	else if( 1 == i){fprintf(out,"10-30\t");}
	else if( 2 == i){fprintf(out,"30-50\t");}
	else if( 3 == i){fprintf(out,"50-70\t");}
	else if( 4 == i){fprintf(out,"70-100\t");}

	if( 0 ==i )	{printf("\t0-10\t");}
	else if( 1 == i){printf("\t10-30\t");}
	else if( 2 == i){printf("\t30-50\t");}
	else if( 3 == i){printf("\t50-70\t");}
	else if( 4 == i){printf("\t70-100\t");}
*/
	fprintf(out,"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n");

	for(i=0;i!=steps;i++){
		for(r=0;r!=ROW;r++){

//			if( 0 == i ){N_tot++;}//Total number of objects 	

			r_tmp = snap[r][rad]/r_max; //Radius of certain object
		
			if( (double)i/(double)steps < r_tmp && (double)(i+1)/(double)steps >= r_tmp ){
				for(m=4;m!=9;m++){
					if( (double)m <= snap[r][mi] && (double)(m+1) > snap[r][mi] ){

			if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){N_ms[m-4]++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){N_msms_sq[m-4]++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim <= snap[r][sm2]/snap[r][sm1] ){N_msms_lq[m-4]++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){N_msc[m-4]++;}
			if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){N_msc[m-4]++;}

					}
				}
			}

		r_av = (double)(2*i+1)/(double)(2*steps); //Avarage radius of interval [r_i,r_{i+1}]

		}
//Here we devide r_av over r_L50% to normalize it to half-mass radius

		fprintf(out,"%lf %i %i %i %i %i ",r_av/r_L50,N_ms[0],N_ms[1],N_ms[2],N_ms[3],N_ms[4]);
		fprintf(out,"%i %i %i %i %i ",N_msms_sq[0],N_msms_sq[1],N_msms_sq[2],N_msms_sq[3],N_msms_sq[4]);
		fprintf(out,"%i %i %i %i %i ",N_msms_lq[0],N_msms_lq[1],N_msms_lq[2],N_msms_lq[3],N_msms_lq[4]);
		fprintf(out,"%i %i %i %i %i\n",N_msc[0],N_msc[1],N_msc[2],N_msc[3],N_msc[4]);
	}

	for( m = 0 ; m != 5 ; m++ ){ N_tot[0][m] = N_ms[m]; }
	for( m = 0 ; m != 5 ; m++ ){ N_tot[1][m] = N_msms_sq[m]; }
	for( m = 0 ; m != 5 ; m++ ){ N_tot[2][m] = N_msms_lq[m]; }
	for( m = 0 ; m != 5 ; m++ ){ N_tot[3][m] = N_msc[m]; }

	fclose(out);

	return N_tot;

}

void SysNumberRadius(double **snap, char *model){
	
	int i,r;
	int ms;
	int msms0;
	int msms1;
	int mswd;
	int tot = 0;

	double **mass_fun;
	if(NULL == (mass_fun = malloc((STEPS)*sizeof(double*)))){perror("malloc\n");exit(1);}
	for(i=0;i!=STEPS;i++){
		if(NULL == (mass_fun[i] = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	}

	double *radi,*mass;
	if(NULL == (radi=malloc(STEPS*sizeof(double)))){perror("RadialObsBinFracFunction malloc fail");exit(1);}
	if(NULL == (mass=malloc(STEPS*sizeof(double)))){perror("RadialObsBinFracFunction malloc fail");exit(1);}

	double r_MAX;
	double r_TMP;
	double L[6];
/*	L[0] = 0
 *	L[1] = L10%
 *	L[2] = L30%
 *	L[3] = L50%
 *	L[4] = L70%
 *	L[5] = L100%
 */
	r_MAX=MaxRadius(snap);

	FILE *fmassfun;
	fmassfun=fopen("fmassfun.dat","w");

	mass_fun = returnRadialMassFunction(snap,STEPS);
	for(i=0;i!=STEPS;i++){
		radi[i]=mass_fun[i][0];
		mass[i]=mass_fun[i][1];
		fprintf(fmassfun,"%lf\t%lf\n",radi[i],mass[i]);
	}
	fclose(fmassfun);

	FILE *out;

	char filename[60];
	sprintf(filename,"snr_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("snr fopen failure");exit(1);}

	L[0]=0;						printf("\trL0=%lf\n",L[0]);
	L[1]=LinearInterpolation(mass,radi,STEPS,0.1);	printf("\trL10=%lf\n",L[1]);
	L[2]=LinearInterpolation(mass,radi,STEPS,0.3);	printf("\trL30=%lf\n",L[2]);
	L[3]=LinearInterpolation(mass,radi,STEPS,0.5);	printf("\trL50=%lf\n",L[3]);
	L[4]=LinearInterpolation(mass,radi,STEPS,0.7);	printf("\trL70=%lf\n",L[4]);
	L[5]=1;						printf("\trL100=%lf\n",L[5]);

	printf("\tReal number of systems in globular cluster\n");
	printf("\tInterval[r_L]\tMS\t\tMS-MS(q<0.5)\tMS-MS(q>0.5)\tMS-WD\n");

	for(i=0;i!=5;i++){

	if( 0 ==i )	{fprintf(out,"0-10\t");}
	else if( 1 == i){fprintf(out,"10-30\t");}
	else if( 2 == i){fprintf(out,"30-50\t");}
	else if( 3 == i){fprintf(out,"50-70\t");}
	else if( 4 == i){fprintf(out,"70-100\t");}

	if( 0 ==i )	{printf("\t0-10\t");}
	else if( 1 == i){printf("\t10-30\t");}
	else if( 2 == i){printf("\t30-50\t");}
	else if( 3 == i){printf("\t50-70\t");}
	else if( 4 == i){printf("\t70-100\t");}

	ms = 0;
	msms0 = 0;
	msms1 = 0;
	mswd = 0;

		for(r=0;r!=ROW;r++){

			r_TMP=snap[r][rad]/r_MAX; //as a part of max radius

			if(r_TMP>L[i] && r_TMP<=L[i+1]){				
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){ms++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){msms0++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim <= snap[r][sm2]/snap[r][sm1] ){msms1++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){mswd++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){mswd++;tot++;}
			}
		}		

//Normalization by total number of systems
	double nms = (double)ms/tot; 
	double nmsms0 = (double)msms0/tot;
	double nmsms1 = (double)msms1/tot;
	double nmswd = (double)mswd/tot;

	printf("\t%lf\t%lf\t%lf\t%lf\n",nms,nmsms0,nmsms1,nmswd);
	fprintf(out,"%lf\t%lf\t%lf\t%lf\n",nms,nmsms0,nmsms1,nmswd);

	}
	fclose(out);	
}

void SysNumberMagnitude(double **snap, char *model){
	
	printf("\n\t11 intervals in magnitude range [4:13] mag\n");
	printf("Real number of systems in globular cluster\n");
	printf("\tInterval[mag]\tMS\t\tMS-MS(q<0.5)\tMS-MS(q>0.5)\tMS-WD\n");

	int i,r;
	int ms;
	int msms0;
	int msms1;
	int mswd;
	int tot = 0;
	
	FILE *out;

	char filename[60];
	sprintf(filename,"snm_%s.dat",model);

	if( NULL == (out=fopen(filename,"w"))){perror("snm fopen failure");exit(1);}

	for(i=4;i!=13;i++){

	ms = 0;
	msms0 = 0;
	msms1 = 0;
	mswd = 0;

		for(r=0;r!=ROW;r++){
			if( (double)i <= snap[r][mi] &&  (double)i+1 > snap[r][mi] ){
				if( 0 == snap[r][ikb] && 2 > snap[r][ik1] ){ms++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim >= snap[r][sm2]/snap[r][sm1] ){msms0++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 2 > snap[r][ik2] && q_lim <= snap[r][sm2]/snap[r][sm1] ){msms1++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik1] && 12 >= snap[r][ik2] && 10 <= snap[r][ik2] ){mswd++;tot++;}
				if( 0 != snap[r][ikb] && 2 > snap[r][ik2] && 12 >= snap[r][ik1] && 10 <= snap[r][ik1] ){mswd++;tot++;}
			}
		}		

//Normalization by total number of systems
	double nms = (double)ms/tot; 
	double nmsms0 = (double)msms0/tot;
	double nmsms1 = (double)msms1/tot;
	double nmswd = (double)mswd/tot;

	printf("\t[%i,%i]\t\t%lf\t%lf\t%lf\t%lf\n",i,i+1,nms,nmsms0,nmsms1,nmswd);
	fprintf(out,"[%i,%i]\t%lf\t%lf\t%lf\t%lf\n",i,i+1,nms,nmsms0,nmsms1,nmswd);

	}
	fclose(out);	
}

void SysNumberDistribution(double **snap, char *model){

	char string[4];
systems_more:
	printf("Number of systems(magnitude)\t-> snm\n");
	printf("Number of systems(radius)\t-> snr\n");

	printf("\nReturn\t-> R\n");

	scanf("%s",string);

	if( 0 == strcmp(string,"R") ){system("clear");}
	else if( 0 == strcmp(string,"snm") ){system("clear");SysNumberMagnitude(snap,model);goto systems_more;}
	else if( 0 == strcmp(string,"snr") ){system("clear");SysNumberRadius(snap,model);goto systems_more;}
	else{system("clear");printf("I don't get it. :c\n");goto systems_more;}
}

void Plot_SysNrRadialDistribution_New(double **snap, char *model){

	int **N_tot;

	if(NULL == (N_tot = calloc(4,sizeof(int*)))){perror("calloc\n");exit(1);}
	for(int i=0;i!=4;i++){
		if(NULL == (N_tot[i] = calloc(5,sizeof(int)))){perror("calloc\n");exit(1);}
	}

	N_tot = SysNrRadialDistribution_New(snap,model);


	FILE *gnuplot;

	char filename[60];
	char path[90];

	sprintf(path,"bnr_new_%s.dat",model);
//_________________________________________________MS stars
	sprintf(filename,"ms_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'ms_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of MS stars - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set xrange [0:%i]\n",XLIM);
		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:($2/%lf) w l t '[4:5] mag','bnr_new_%s.dat' u 1:($3/%lf) w l t '[5:6] mag','bnr_new_%s.dat' u 1:($4/%lf) w l t '[6:7] mag','bnr_new_%s.dat' u 1:($5/%lf) w l t '[7:8] mag','bnr_new_%s.dat' u 1:($6/%lf) w l t '[8:9] mag'\n",model,(double)N_tot[0][0],model,(double)N_tot[0][1],model,(double)N_tot[0][2],model,(double)N_tot[0][3],model,(double)N_tot[0][4]);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: ms_new_%s.gpi\n\n",model);
//_________________________________________________MS-MS, q < 0.5
	sprintf(filename,"msms_sq_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'msms_sq_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of MS-MS for q<0.5 - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set xrange [0:%i]\n",XLIM);
		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:($7/%lf) w l t '[4:5] mag','bnr_new_%s.dat' u 1:($8/%lf) w l t '[5:6] mag','bnr_new_%s.dat' u 1:($9/%lf) w l t '[6:7] mag','bnr_new_%s.dat' u 1:($10/%lf) w l t '[7:8] mag','bnr_new_%s.dat' u 1:($11/%lf) w l t '[8:9] mag'\n",model,(double)N_tot[1][0],model,(double)N_tot[1][1],model,(double)N_tot[1][2],model,(double)N_tot[1][3],model,(double)N_tot[1][4]);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: msms_sq_new_%s.gpi\n\n",model);
//_________________________________________________MS-MS, q > 0.5
	sprintf(filename,"msms_lq_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'msms_lq_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of MS-MS for q>0.5 - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set xrange [0:%i]\n",XLIM);
		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:($12/%lf) w l t '[4:5] mag','bnr_new_%s.dat' u 1:($13/%lf) w l t '[5:6] mag','bnr_new_%s.dat' u 1:($14/%lf) w l t '[6:7] mag','bnr_new_%s.dat' u 1:($15/%lf) w l t '[7:8] mag','bnr_new_%s.dat' u 1:($16/%lf) w l t '[8:9] mag'\n",model,(double)N_tot[2][0],model,(double)N_tot[2][1],model,(double)N_tot[2][2],model,(double)N_tot[2][3],model,(double)N_tot[2][4]);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: msms_lq_new_%s.gpi\n\n",model);

//_________________________________________________MS-compact
	sprintf(filename,"msc_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'msc_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of MS-compact stars - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set xrange [0:%i]\n",XLIM);
		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:($17/%lf) w l t '[4:5] mag','bnr_new_%s.dat' u 1:($18/%lf) w l t '[5:6] mag','bnr_new_%s.dat' u 1:($19/%lf) w l t '[6:7] mag','bnr_new_%s.dat' u 1:($20/%lf) w l t '[7:8] mag','bnr_new_%s.dat' u 1:($21/%lf) w l t '[8:9] mag'\n",model,(double)N_tot[3][0],model,(double)N_tot[3][1],model,(double)N_tot[3][2],model,(double)N_tot[3][3],model,(double)N_tot[3][4]);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: msc_new_%s.gpi\n\n",model);
//__________________________________________________Cumulative distribution (TOTAL)

int N_tot_ms = 0;
int N_tot_msms_sq = 0;
int N_tot_msms_lq = 0;
int N_tot_msc = 0;

	for(int i = 0 ; i != 4 ; i++){
		if(0 == i) for(int j = 0 ; j != 5 ; j++){N_tot_ms += N_tot[i][j];}
		if(1 == i) for(int j = 0 ; j != 5 ; j++){N_tot_msms_sq += N_tot[i][j];}
		if(2 == i) for(int j = 0 ; j != 5 ; j++){N_tot_msms_lq += N_tot[i][j];}
		if(3 == i) for(int j = 0 ; j != 5 ; j++){N_tot_msc += N_tot[i][j];}
	}

	sprintf(filename,"bnr_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'bnr_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Dystrybuanty rozkładów liczby poszczególnych obiektów w zależności od promienia - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'Względna liczba układów N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Promień połowy masy'\n");

		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set xrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:(($2+$3+$4+$5+$6)/%lf) w l t 'MS','bnr_new_%s.dat' u 1:(($7+$8+$9+$10+$11)/%lf) w l t 'MS-MS, q<0.5','bnr_new_%s.dat' u 1:(($12+$13+$14+$15+$16)/%lf) w l t 'MS-MS, q>0.5','bnr_new_%s.dat' u 1:(($17+$18+$19+$20+$21)/%lf) w l t 'MS-compact',\n",model,(double)N_tot_ms,model,(double)N_tot_msms_sq,model,(double)N_tot_msms_lq,model,(double)N_tot_msc);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: n45_new_%s.gpi\n\n",model);
//__________________________________________________[4:5] mag
	sprintf(filename,"n45_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'n45_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of objects for [4:5] mag - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:2 w l t 'MS','bnr_new_%s.dat' u 1:7 w l t 'MS-MS, q<0.5','bnr_new_%s.dat' u 1:12 w l t 'MS-MS, q>0.5','bnr_new_%s.dat' u 1:17 w l t 'MS-compact(without mag limits)',\n",model,model,model,model);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: n45_new_%s.gpi\n\n",model);
//__________________________________________________[5:6] mag
	sprintf(filename,"n56_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'n56_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of objects for [5:6] mag - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:3 w l t 'MS','bnr_new_%s.dat' u 1:8 w l t 'MS-MS, q<0.5','bnr_new_%s.dat' u 1:13 w l t 'MS-MS, q>0.5','bnr_new_%s.dat' u 1:17 w l t 'MS-compact(without mag limits)',\n",model,model,model,model);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: n56_new_%s.gpi\n\n",model);
//__________________________________________________[6:7] mag
	sprintf(filename,"n67_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'n67_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of objects for [6:7] mag - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:4 w l t 'MS','bnr_new_%s.dat' u 1:9 w l t 'MS-MS, q<0.5','bnr_new_%s.dat' u 1:14 w l t 'MS-MS, q>0.5','bnr_new_%s.dat' u 1:17 w l t 'MS-compact(without mag limits)',\n",model,model,model,model);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: n67_new_%s.gpi\n\n",model);
//__________________________________________________[7:8] mag
	sprintf(filename,"n78_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'n78_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of objects for [7:8] mag - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:5 w l t 'MS','bnr_new_%s.dat' u 1:10 w l t 'MS-MS, q<0.5','bnr_new_%s.dat' u 1:15 w l t 'MS-MS, q>0.5','bnr_new_%s.dat' u 1:17 w l t 'MS-compact(without mag limits)',\n",model,model,model,model);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: n78_new_%s.gpi\n\n",model);
//__________________________________________________[8:9] mag
	sprintf(filename,"n89_new_%s.gpi",model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'n89_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Cumulative distribution of objects for [8:9] mag - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N/N_{total}'\n");
		fprintf(gnuplot,"set xlabel 'Half-Mass Radius'\n");

		fprintf(gnuplot,"set yrange [0:]\n");
		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bnr_new_%s.dat' u 1:6 w l t 'MS','bnr_new_%s.dat' u 1:11 w l t 'MS-MS, q<0.5','bnr_new_%s.dat' u 1:16 w l t 'MS-MS, q>0.5','bnr_new_%s.dat' u 1:17 w l t 'MS-compact(without mag limits)',\n",model,model,model,model);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: n89_new_%s.gpi\n\n",model);

}

void Plot_RadialBinFracFunction_New(double **snap, char *model){

	FILE *gnuplot;

	char filename[60];
	char path[90];

	sprintf(path,"bfr_new_%s.dat",model);

	sprintf(filename,"bfr_new_%s.gpi",model);

	RadialBinFracFunction_New(snap,model);

	//double bin_frac_total = RealBinaryFraction(snap);
	double bin_frac_true = TrueBinaryFraction(snap);

	double bin_frac_o_47;	//binary fraction observationally for 4:7 magnitudes
	bin_frac_o_47 = BinaryFractionDistribution_New(snap, model);

	if(NULL == (gnuplot = fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'bfr_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Rozkład radialny udziału układów podwójnych - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'Udział układów podwójnych'\n");
		fprintf(gnuplot,"set xlabel 'Przedziały procentu promienia Lagranga'\n");

		fprintf(gnuplot,"set style data histogram\n");
		fprintf(gnuplot,"set style histogram clustered gap 1\n");
		fprintf(gnuplot,"set style fill pattern 2 border -1\n");

		fprintf(gnuplot,"set yrange [0:]\n");

		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"f(x) = %.7lf\n",bin_frac_true);
		fprintf(gnuplot,"g(x) = %.7lf\n",bin_frac_o_47);

		fprintf(gnuplot,"plot 'bfr_new_%s.dat' u 2:xtic(1) t 'układy podwójne, podejście obserwacyjne w przedziale jasności [4:9] mag', 'bfr_new_%s.dat' u 3:xtic(1) t 'układy MS-MS w przedziale jasności [4:9] mag', 'bfr_new_%s.dat' u 4:xtic(1) t 'Układy podwójne, podejście rzeczywiste, bez ograniczeń w jasności',f(x) t 'Udział układów podwójnych w ujęciu globalnym'\n",model,model,model);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: bfr_new_%s.gpi\n\n",model);

}


void Plot_BinaryFractionDistribution_New(double **snap, char *model){

	FILE *gnuplot;

	char filename[60];
	char path[90];
	
	double bin_frac_o_47;	//binary fraction observationally for 4:7 magnitudes

	sprintf(path,"bfm_new_%s.dat",model);

	sprintf(filename,"bfm_new_%s.gpi",model);

	bin_frac_o_47 = BinaryFractionDistribution_New(snap, model);

	//double bin_frac_total = RealBinaryFraction(snap);
	double bin_frac_true = TrueBinaryFraction(snap);
	
	delta_bin_frac[g] = (bin_frac_true - bin_frac_o_47)/bin_frac_true;
	g++;

	if(NULL == (gnuplot=fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}
	
		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'bfm_new_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Rozkład f_{bin} w przedziałach jasności w filtrze I  - %s'\n",model);
		
		fprintf(gnuplot,"set yrange [0:]\n");

		fprintf(gnuplot,"set style data histogram\n");
		fprintf(gnuplot,"set style histogram clustered gap 1\n");
		fprintf(gnuplot,"set style fill pattern 8 border -1\n");

		fprintf(gnuplot,"set xlabel 'Przedział jasności, filtr I [mag]'\n");
		fprintf(gnuplot,"set ylabel 'Udział układów podwójnych'\n");

		fprintf(gnuplot,"set key font ',10'\n");
		fprintf(gnuplot,"set key outside right tmargin\n");
		
		fprintf(gnuplot,"f(x) = %.7lf\n",bin_frac_true);
		fprintf(gnuplot,"g(x) = %.7lf\n",bin_frac_o_47);
		
		fprintf(gnuplot,"plot '%s' u 2:xtic(1) t 'MS-MS, q>0.5','%s' u 3:xtic(1) t 'układy MS-MS ', '%s' u 4:xtic(1) t 'MS-MS, MS-comp', f(x) t 'f_{bin} globalne', g(x) t 'MS-MS z q>0.5 dla jasności [4,7] mag'\n",path,path,path);
		
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: bfm_new_%s.gpi\n\n",model);

}

void Plot_MagnitudeDependance(){

	char path[90];
	char model[60];
	
COL = 30;

ROW = 559221;
	sprintf(model,"diogo-ibp-modified-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 242457;
	sprintf(model,"diogo-ibp-modified-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 12691;
	sprintf(model,"diogo-ibp-modified-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 1565455;
	sprintf(model,"diogo-ibp-modified-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 683393;
	sprintf(model,"diogo-ibp-modified-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 193538;
	sprintf(model,"diogo-ibp-modified-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 552773;
	sprintf(model,"diogo-ibp-original-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 242971;
	sprintf(model,"diogo-ibp-original-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 11956;
	sprintf(model,"diogo-ibp-original-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 1547011;
	sprintf(model,"diogo-ibp-original-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 710530;
	sprintf(model,"diogo-ibp-original-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);
ROW = 212316;
	sprintf(model,"diogo-ibp-original-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	ObservationalMagnitudeDependance(LoadFile(ROW,COL,model),model);

	FILE *gnuplot;
	if(NULL == (gnuplot=fopen("magdep.gp","w"))){perror("gnuplot fopen failure");exit(1);}
	
		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'results/magdep.pdf'\n");

		fprintf(gnuplot,"set xlabel 'Przedział jasności, filtr I [mag]'\n");
		fprintf(gnuplot,"set ylabel 'f_{obs}^{[i,i+1]}/f_{tot}'\n");

		fprintf(gnuplot,"unset key\n");
		
		fprintf(gnuplot,"plot 'magdep_diogo-ibp-original-1.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-original-2.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-original-3.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-original-4.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-original-5.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-original-6.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-modified-1.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-modified-2.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-modified-3.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-modified-4.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-modified-5.dat' u 2:xtic(1) w p pt 7 lc 'black','magdep_diogo-ibp-modified-6.dat' u 2:xtic(1) w p pt 7 lc 'black'");
		
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: magdep.gp\n\n");

}

void Plot_GlobalParameters_New(){
	//First modified, next original. km/s
	//double vel_disp[12] = {33.8, 6.0, 4.7, 66.9, 9.6, 10.8, 34.5, 6.1, 4.5, 45.7, 10.1, 11.0};
	//First modified, next original. 10^5
	//double N_total[12] = {6.24, 2.99, 0.16, 16.98, 8.22, 2.25, 6.2, 3.0, 0.2, 16.9, 8.5, 2.5};
	double R_tidal[12] = {49.4, 39.5, 8.9, 50.6, 40.3, 14.5, 49.4, 39.2, 8.9, 50.6, 40.8, 14.9};
	double R_half[12] = {5.5, 8.1, 1.3, 4.9, 8.5, 2.3, 5.6, 8.1, 1.4, 4.9, 8.6, 2.4};
	double R_tidal_initial[12] = {60, 60, 30, 60, 60, 30, 60, 60, 30, 60, 60, 30 };
	double R_half_initial[12] = {2.4, 8.8, 4.4, 2.4, 8.8, 4.4, 2.4, 8.8, 4.4, 2.4, 8.8, 4.4};
	
	FILE *gnuplot;
	FILE *out;
	
	char filename[60];
	char path[90];
	
	sprintf(path,"global.dat");
	sprintf(filename,"global.gpi");
	
	if(NULL == (out=fopen(path,"w"))){perror("global fopen failure");exit(1);}
		for(int i = 0; i != 12; i++){
			fprintf(out,"%lf %lf %lf\n",delta_bin_frac[i],R_tidal_initial[i]/R_half_initial[i],R_tidal[i]/R_half[i]);
		}
	fclose(out);
	
	if(NULL == (gnuplot=fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}
	
		fprintf(gnuplot,"set terminal pdf enhanced\n");
		fprintf(gnuplot,"set out 'global.pdf'\n");

		//
		fprintf(gnuplot,"set title 'Funkcja stosunku promienia pływowego do promienia połowy masy(Warunki początkowe)'\n");
		
		fprintf(gnuplot,"set xlabel 'Stosunek promienia pływowego do promienia połowy masy'\n");
		fprintf(gnuplot,"set ylabel '{/Symbol D}f=(f_{true}-f_{obs})/f_{true}'\n");

		fprintf(gnuplot,"set key font ',8'\n");
		fprintf(gnuplot,"set key right top\n");
		
		fprintf(gnuplot,"plot '%s' u 2:1 w p pt 7 lc 'blue' t 'R_{tidal}/R_{half-mass}' \n",path);
		//
		fprintf(gnuplot,"set title 'Funkcja stosunku promienia pływowego do promienia połowy masy(Po 12 mld lat)'\n");
		
		fprintf(gnuplot,"set xlabel 'Stosunek promienia pływowego do promienia połowy masy'\n");
		fprintf(gnuplot,"set ylabel '{/Symbol D}f'\n");

		fprintf(gnuplot,"set key font ',8'\n");
		fprintf(gnuplot,"set key right top\n");
		
		fprintf(gnuplot,"plot '%s' u 3:1 w p pt 7 lc 'red' t 'R_{tidal}/R_{half-mass}'\n",path);
	fclose(gnuplot);
}


//___________________________________________________________
void PlotSysNumberMagnitude(double **snap, char *model){
	printf("I dunno what to do :c\n");
}

void PlotSysNumberRadius(double **snap, char *model){
	printf("I dunno what to do :c\n");
}

void PlotCMD(double **snap, char *model,char *path){

	FILE *gnuplot;
	char filename[60];

	sprintf(filename,"cmd_binaries_%s.gpi",model);	

	printf("%s\n",filename);

	if(NULL == (gnuplot=fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'cmd_binaries_%s.pdf'\n",model);
		
		fprintf(gnuplot,"set title '%s'\n",model);
		
		fprintf(gnuplot,"set xlabel 'V-I[mag]'\n");
		fprintf(gnuplot,"set ylabel 'I[mag]'\n");
		
		fprintf(gnuplot,"set yrange [11:1]\n");
		fprintf(gnuplot,"set xrange [-0.5:2.5]\n");
		
		fprintf(gnuplot,"set key left bottom\n");
	
		fprintf(gnuplot,"plot '%s' u ($22-$24):( $7==1 && $8<2 && $9<2 ? $24 : 1/0) w d lc 'grey' t 'grey: MS-MS', '%s' u ($22-$24):( $7==0 && $8<2 ? $24 : 1/0) w d lc 'black' t 'black: MS single', '%s' u ($22-$24):( ($7==1 && $8<2 && $9>=10 && $9<=12) || ($7==1 && $9<2 && $8>=10 && $8 >=12) ? $24 : 1/0) w d lc 'red' t 'red: MS-WD', '%s' u ($22-$24):( ($7==1 && $8<2 && $9>=2 && $9<=9) || ($7==1 && $9<2 && $8>=2 && $8<=9) ? $24 : 1/0) w d lc 'green' t 'green: MS-HG,GB,He', '%s' u ($22-$24):( ($7==1 && $8>=10 && $8<=12 && $9>=2 && $9<=9) || ($7==1 && $9>=10 && $9<=12 && $8>=2 && $8<=9) ? $24 : 1/0) w d lc 'blue' t 'blue: WD-HG,GB,He'",path,path,path,path,path);

	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRIPT: cmd_binaries_%s.gpi\n\n",model);
}

void PlotBinFracMagnitude(double **snap, char *model){

	FILE *gnuplot;

	char filename[60];
	char path[90];

	sprintf(path,"bfm_%s.dat",model);

	sprintf(filename,"bfm_%s.gpi",model);

	BinaryFractionDistribution(snap, model);

	if(NULL == (gnuplot=fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}
	
		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'bfm_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Binary fraction - magnitude distribution  - %s'\n",model);
		
		fprintf(gnuplot,"set yrange [0:]\n");

		fprintf(gnuplot,"set style data histogram\n");
		fprintf(gnuplot,"set style histogram clustered gap 1\n");
		fprintf(gnuplot,"set style fill pattern 2 border -1\n");

		fprintf(gnuplot,"set xlabel 'I magnitude interval [mag]'\n");
		fprintf(gnuplot,"set ylabel 'Binary fraction'\n");

		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot '%s' u 2:xtic(1) t 'observational binary fraction', '%s' u 3:xtic(1) t 'real binary fraction'\n",path,path);
		
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: bfm_%s.gpi\n\n",model);
}

void PlotBinFracMassRatio(double **snap, char *model){

	FILE *gnuplot;

	char filename[60];
	char path[90];

	sprintf(path,"bfq_%s.dat",model);

	sprintf(filename,"bfq_%s.gpi",model);

	double bin_frac_tot = MassRatioDistribution(snap, model);

	if(NULL == (gnuplot=fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}
	
		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'bfq_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Binary fraction - mass ratio distribution - %s, f_{TOT}=%lf'\n",model,bin_frac_tot);

		fprintf(gnuplot,"max = %lf\n",bin_frac_tot);

		fprintf(gnuplot,"set yrange [0:]\n");

		fprintf(gnuplot,"set style data histogram\n");
		fprintf(gnuplot,"set style histogram rowstacked\n");
		fprintf(gnuplot,"set style fill pattern 2 border -1\n");

		fprintf(gnuplot,"set xlabel 'Mass ratio interval'\n");
		fprintf(gnuplot,"set ylabel 'Binary fraction'\n");

		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot '%s' u 2:xtic(1) t 'real binary fraction'\n",path);
		
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: bfq_%s.gpi\n\n",model);
}

void PlotBinNumberRadius(double **snap, char *model){

	FILE *gnuplot;

	char filename[60];
	char path[90];

	sprintf(path,"bnr_%s.dat",model);

	sprintf(filename,"bnr_%s.gpi",model);

	double limit_r = RadialRealBinFracFunction(snap,model);
	double limit_o = RadialObsBinFracFunction(snap,model);

	if(NULL == (gnuplot=fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'bnr_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Radial distribution of binaries - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N_{i,bin}'\n");
		fprintf(gnuplot,"set xlabel 'Percent of the Lagrange Radius interval'\n");

		fprintf(gnuplot,"set style data histogram\n");
		fprintf(gnuplot,"set style histogram clustered gap 1\n");
		fprintf(gnuplot,"set style fill pattern 2 border -1\n");

		fprintf(gnuplot,"set yrange [0:]\n");

		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bfr_r_%s.dat' u 3:xtic(1) t 'real binary fraction(%.1lf mag limit )', 'bfr_o_%s.dat' u 3:xtic(1) t 'observational binary fraction (%.1lf mag limit)'\n",model,limit_r,model,limit_o);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: bnr_%s.gpi\n\n",model);

}

void PlotBinFracRadius(double **snap, char *model){

	FILE *gnuplot;

	char filename[60];
	char path[90];

	sprintf(path,"bfr_%s.dat",model);

	sprintf(filename,"bfr_%s.gpi",model);

	double limit_r = RadialRealBinFracFunction(snap,model);
	double limit_o = RadialObsBinFracFunction(snap,model);

	if(NULL == (gnuplot=fopen(filename,"w"))){perror("gnuplot fopen failure");exit(1);}

		fprintf(gnuplot,"set terminal pdf\n");
		fprintf(gnuplot,"set out 'bfr_%s.pdf'\n",model);

		fprintf(gnuplot,"set title 'Radial distribution of binary fraction - %s'\n",model);

		fprintf(gnuplot,"set ylabel 'N_{i,bin}/(N_{i,bin}+N_{i,sin})'\n");
		fprintf(gnuplot,"set xlabel 'Percent of the Lagrange Radius interval'\n");

		fprintf(gnuplot,"set style data histogram\n");
		fprintf(gnuplot,"set style histogram clustered gap 1\n");
		fprintf(gnuplot,"set style fill pattern 2 border -1\n");

		fprintf(gnuplot,"set yrange [0:]\n");

		fprintf(gnuplot,"set key right top\n");

		fprintf(gnuplot,"plot 'bfr_r_%s.dat' u 2:xtic(1) t 'real binary fraction(%.1lf mag limit )', 'bfr_o_%s.dat' u 2:xtic(1) t 'observational binary fraction (%.1lf mag limit)'\n",model,limit_r,model,limit_o);
	
	fclose(gnuplot);

	printf("\n\t---> OUTPUT SCRPIT: bfr_%s.gpi\n\n",model);

}

void PlotRadialMassDistribution(double **snap){

	printf("There is nothing interesting. Go away.\n");

}

void PlotAll(double **snap, char *model, char *path){

	PlotCMD(snap,model,path);
	PlotBinFracMagnitude(snap,model);
	PlotBinFracMassRatio(snap,model);
	PlotBinFracRadius(snap,model);
	PlotBinNumberRadius(snap,model);
}

//___________________________________________________________Latter Bachelor analysis
void OverallAnalysis(double **snap, char *model){

/* In bins of magnitude: 
 * 1) f_bin_r_i - fraction of all binaries (for 4-9 mag)
 * 2) f_bin_o_i - fraction of MS-MS_q<.5 binaries (for 4-9 mag)
 * 3) f_MS-MS_all - fraction of all MS-MS systems (for 4-9 mag)
 *
 * In bins of Lagrange radius:
 * 4) 1) but without mag range
 * 5) 2)
 * 6) 3)
 */
	Plot_BinaryFractionDistribution_New(snap,model);
	Plot_RadialBinFracFunction_New(snap,model);
	Plot_SysNrRadialDistribution_New(snap,model);

}

void Plotting(double **snap, char *model, char *path){

	char string[4];
plot_more:
	printf("What do you want to plot?\n");
	printf("binaries on CMD\t\t\t-> cmd\n");
	printf("binary fraction(magnitude)\t-> bfm\n");
	printf("binary fraction(mass-ratio)\t-> bfq\n");
	printf("binary fraction(radius)\t\t-> bfr\n");
	printf("number of binaries(radius)\t-> bnr\n");
	printf("radial mass distribution\t-> mr\n");
	printf("number of systems(magnitude)\t-> snm\n");
	printf("number of systems(radius)\t-> snr\n");
	printf("full service\t\t\t-> full\n");
	printf("new overall analysis\t\t-> fuln\n");
	
	printf("\nReturn\t-> R\n");
	
	scanf("%s",string);

	if( 0 ==strcmp(string,"cmd") ){system("clear");PlotCMD(snap,model,path);goto plot_more;}
	else if( 0 == strcmp(string,"bfm") ){system("clear");PlotBinFracMagnitude(snap,model);goto plot_more;}
	else if( 0 == strcmp(string,"bfq") ){system("clear");PlotBinFracMassRatio(snap,model);goto plot_more;}
	else if( 0 == strcmp(string,"bfr") ){system("clear");PlotBinFracRadius(snap,model);goto plot_more;}
	else if( 0 == strcmp(string,"bnr") ){system("clear");PlotBinNumberRadius(snap,model);goto plot_more;}
	else if( 0 == strcmp(string,"snm") ){system("clear");PlotSysNumberMagnitude(snap,model);goto plot_more;}
	else if( 0 == strcmp(string,"snr") ){system("clear");PlotSysNumberRadius(snap,model);goto plot_more;}
	else if( 0 == strcmp(string,"mr") ){system("clear");PlotRadialMassDistribution(snap);goto plot_more;}
	else if( 0 == strcmp(string,"full") ){system("clear");PlotAll(snap,model,path);goto plot_more;} 
	else if( 0 == strcmp(string,"fuln") ){system("clear");OverallAnalysis(snap,model);goto plot_more;} 
	else if( 0 == strcmp(string,"R") ){system("clear");}
	else{system("clear");printf("I don't get it. :c\n");goto plot_more;}

}

void ShowModel(char *filename){
	char message[90];
	sprintf(message,"--->\tyour model:%s\n\n",filename);
	printf(message);
}

double **LoadFile(int ROW, int COL, char *model){

	int c,r;
	FILE *fsnap;
	double **snap;
	char filename[60];

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

return snap;

}

double *Fcomparison(double **snap, char *model){
	
	double f_obs = BinaryFractionDistribution_New(snap, model);
	double f_tot = TrueBinaryFraction(snap);
	
	double *f;
	if(NULL == (f = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	
	f[0] = f_obs;
	f[1] = f_tot;
	
	return f;
}

void Ultimate_Fcomparison(void){
	
	char path[90];
	char model[60];
	
	double *f;
	if(NULL == (f = malloc(2*sizeof(double)))){perror("malloc\n");exit(1);}
	
	FILE *fcomp;
	fcomp = fopen("fcomparison.dat","w+");
	
COL = 30;

ROW = 559221;
	sprintf(model,"diogo-ibp-modified-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"M1 %lf %lf\n",f[0],f[1]);
ROW = 242457;
	sprintf(model,"diogo-ibp-modified-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"M2 %lf %lf\n",f[0],f[1]);
ROW = 12691;
	sprintf(model,"diogo-ibp-modified-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"M3 %lf %lf\n",f[0],f[1]);
ROW = 1565455;
	sprintf(model,"diogo-ibp-modified-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"M4 %lf %lf\n",f[0],f[1]);
ROW = 683393;
	sprintf(model,"diogo-ibp-modified-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"M5 %lf %lf\n",f[0],f[1]);
ROW = 193538;
	sprintf(model,"diogo-ibp-modified-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"M6 %lf %lf\n",f[0],f[1]);
ROW = 552773;
	sprintf(model,"diogo-ibp-original-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"O1 %lf %lf\n",f[0],f[1]);
ROW = 242971;
	sprintf(model,"diogo-ibp-original-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"O2 %lf %lf\n",f[0],f[1]);
ROW = 11956;
	sprintf(model,"diogo-ibp-original-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"O3 %lf %lf\n",f[0],f[1]);
ROW = 1547011;
	sprintf(model,"diogo-ibp-original-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"O4 %lf %lf\n",f[0],f[1]);
ROW = 710530;
	sprintf(model,"diogo-ibp-original-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"O5 %lf %lf\n",f[0],f[1]);
ROW = 212316;
	sprintf(model,"diogo-ibp-original-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	f = Fcomparison(LoadFile(ROW,COL,model),model);
	fprintf(fcomp,"O6 %lf %lf\n",f[0],f[1]);
	
	fclose(fcomp);
}
void UltimatePlot_New(void){
	
	char path[90];
	char model[60];

COL = 30;

ROW = 559221;
	sprintf(model,"diogo-ibp-modified-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 242457;
	sprintf(model,"diogo-ibp-modified-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 12691;
	sprintf(model,"diogo-ibp-modified-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 1565455;
	sprintf(model,"diogo-ibp-modified-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 683393;
	sprintf(model,"diogo-ibp-modified-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 193538;
	sprintf(model,"diogo-ibp-modified-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 552773;
	sprintf(model,"diogo-ibp-original-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 242971;
	sprintf(model,"diogo-ibp-original-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 11956;
	sprintf(model,"diogo-ibp-original-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 1547011;
	sprintf(model,"diogo-ibp-original-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 710530;
	sprintf(model,"diogo-ibp-original-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);
ROW = 212316;
	sprintf(model,"diogo-ibp-original-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	OverallAnalysis(LoadFile(ROW,COL,model),model);

	//Global parameters
	Plot_GlobalParameters_New();
}

void UltimatePlot(void){
	
	char path[90];
	char model[60];

COL = 30;

ROW = 559221;
	sprintf(model,"diogo-ibp-modified-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 242457;
	sprintf(model,"diogo-ibp-modified-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 12691;
	sprintf(model,"diogo-ibp-modified-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 1565455;
	sprintf(model,"diogo-ibp-modified-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 683393;
	sprintf(model,"diogo-ibp-modified-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 193538;
	sprintf(model,"diogo-ibp-modified-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 552773;
	sprintf(model,"diogo-ibp-original-1");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 242971;
	sprintf(model,"diogo-ibp-original-2");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 11956;
	sprintf(model,"diogo-ibp-original-3");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 1547011;
	sprintf(model,"diogo-ibp-original-4");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 710530;
	sprintf(model,"diogo-ibp-original-5");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);
ROW = 212316;
	sprintf(model,"diogo-ibp-original-6");
	sprintf(path,"data/%s/new_snap.dat",model);
	PlotAll(LoadFile(ROW,COL,model),model,path);

}

