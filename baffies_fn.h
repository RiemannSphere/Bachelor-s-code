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
double LinearInterpolation(double *x,double *y,int n,double x_p);
double CubicSpline(double *x_tab_in, double *y_tab_in, int n, double p);
void BinariesStatistics(double **snap);
void ObservationalBinaryFraction(double **snap);
double RealBinaryFraction(double **snap);
double TrueBinaryFraction(double **snap);//new new new
void BinaryFraction(double **snap);
double BinaryFractionDistribution_New(double **snap, char *model);//new
void ObservationalMagnitudeDependance(double **snap, char *model);//neeeeeeeeeeew
void BinaryFractionDistribution(double **snap, char *model);
double MassRatioDistribution(double **snap, char *model);
double MaxRadius(double **snap);
double TotalMass(double **snap, double from, double to,double r_MAX);
double **returnRadialMassFunction(double **snap,int INTERVALS);
void RadialMassFunction(double **snap);
void RadialBinFracFunction_New(double **snap, char *model);//new
double RadialRealBinFracFunction(double **snap, char *model);
double RadialObsBinFracFunction(double **snap, char *model);
void RadialBinFracFunction(double **snap, char *model);
void RadialDistribution(double **snap, char *model);
int **SysNrRadialDistribution_New(double **snap, char *model);//new
void SysNumberRadius(double **snap, char *model);
void SysNumberMagnitude(double **snap, char *model);
void SysNumberDistribution(double **snap, char *model);
void Plot_SysNrRadialDistribution_New(double **snap, char *model);//new
void Plot_RadialBinFracFunction_New(double **snap, char *model);//new
void Plot_BinaryFractionDistribution_New(double **snap, char *model);//new
void Plot_MagnitudeDependance(void);//neeeeeeeeeeew
void Plot_GlobalParameters_New();//new new
void PlotSysNumberMagnitude(double **snap, char *model);
void PlotSysNumberRadius(double **snap, char *model);
void PlotCMD(double **snap, char *model,char *path);
void PlotBinFracMagnitude(double **snap, char *model);
void PlotBinFracMassRatio(double **snap, char *model);
void PlotBinNumberRadius(double **snap, char *model);
void PlotBinFracRadius(double **snap, char *model);
void PlotRadialMassDistribution(double **snap);
void PlotAll(double **snap, char *model, char *path);
void OverallAnalysis(double **snap, char *model);
void Plotting(double **snap, char *model, char *path);
void ShowModel(char *filename);
double **LoadFile(int ROW, int COL, char *model);
double *Fcomparison(double **snap, char *model);// same as Ultimate_Fcomparison
void Ultimate_Fcomparison(void); //new, compare fobs and ftot if fits 2*fobs=ftot
void UltimatePlot_New(void);
void UltimatePlot(void);
