#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <nlopt.h>
#include "matncs5.h"

int M,N,m,n;
double **fmatR, **fmatI;
double scale;
double F[1000000];
double sumF;
int ntrial =0;
int logfc ;

///double rfactor(double *,double ,double *, double ** ,double ** , int, int );
double rfactor(int , double *,double *, void *);

int main(int argc, char *argv[])
{
 
  FILE *fpi, *fhkl, *fpo;
  char s[120],*sxyz;
  char *satom="ATOM",*smtr="MTRI",*scry="CRYS";
  int i,iii,j,iatom,incs,isym;
  double x[100000],y[100000],z[100000],occ[100000];
  double fx[100000],fy[100000],fz[100000];
  char pdbl[100000][90];
  char dumm55[55],dtail60a[30],wline[90];
  int ihkl;
  char sh[100];
  double h[1000000], k[1000000], l[1000000], sF[1000000], FC[1000000];
  double PI,TWOPI;
  double in01,in02,in03,in04,in05,in06,in07,in08,in09,in10,in11,in12;
  double hmfx,hmfy,hmfz,kmfx,kmfy,kmfz,lmfx,lmfy,lmfz;
  double aaa=244.87,iaaa=0.004083799567;

  double rf,minf;

  int ires;



  double fx2,fy2,fz2,fx3,fy3,fz3;
  double *grad;
  void *f_data;

  PI = 4.0*atan(1.0);
  TWOPI = 2.00*PI;


  fpi=fopen(argv[1],"r");
  if (fpi == NULL){
    printf("cannot open file %s\n",argv[1]);
    return -1;
      }
  iatom=0;
  while (fgets(s,100,fpi) != NULL) {
    if (strncmp(s,satom,4) == 0){

      sxyz=s+30;
      sscanf(sxyz,"%8lf%8lf%8lf %5lf ",&x[iatom],&y[iatom],&z[iatom],&occ[iatom]);
      fx[iatom]=x[iatom]*iaaa;
      fy[iatom]=y[iatom]*iaaa;
      fz[iatom]=z[iatom]*iaaa;
      strncpy(pdbl[iatom],s,90);
      iatom++;

    } else if (strncmp(s,smtr,4) == 0){
    } else if (strncmp(s,scry,4) == 0){      
    }
  }
  fclose(fpi);
  N=iatom;

  fhkl=fopen(argv[2],"r");
  if (fpi == NULL){
    printf("cannot open file %s\n",argv[2]);
    return -1;
      }
  ihkl=0;
  while (fgets(sh,100,fpi) != NULL) {
    sscanf(sh," %lf %lf %lf %lf %lf ",&h[ihkl],&k[ihkl],&l[ihkl],&F[ihkl],&sF[ihkl]);
    sumF=sumF+F[ihkl];
      ihkl++;
  }
  fclose(fhkl);
  M=ihkl;

  /// fmat[M][N]  
  fmatR = malloc(sizeof(double *) * M);
  for (i=0;i<M;i++) fmatR[i]= malloc(sizeof(double) * N );
  fmatI = malloc(sizeof(double *) * M);
  for (i=0;i<M;i++) fmatI[i]= malloc(sizeof(double) * N );


  for (m=0;m<M;m++){
  for (n=0;n<N;n++){
    fmatR[m][n]=0.;    fmatI[m][n]=0.;
    for (incs=1;incs<6;incs++){
      fx2 = ncsmat[incs][1][1]*fx[n] + ncsmat[incs][1][2]*fy[n] + ncsmat[incs][1][3]*fz[n] ;
      fy2 = ncsmat[incs][2][1]*fx[n] + ncsmat[incs][2][2]*fy[n] + ncsmat[incs][2][3]*fz[n] ;
      fz2 = ncsmat[incs][3][1]*fx[n] + ncsmat[incs][3][2]*fy[n] + ncsmat[incs][3][3]*fz[n] ;


    /// For P23  symmetry    Unrolling for fast calculation
      /*
        fx3 = sym195[isym][1][1]*fx2                                                   ;
        fy3 =                          sym195[isym][2][2]*fy2                          ;
        fz3 =                                                   sym195[isym][3][3]*fz2 ;
      */

       hmfx=h[m]*fx2;
       hmfy=h[m]*fy2;
       hmfz=h[m]*fz2;

       kmfx=k[m]*fx2;
       kmfy=k[m]*fy2;
       kmfz=k[m]*fz2;

       lmfx=l[m]*fx2;
       lmfy=l[m]*fy2;
       lmfz=l[m]*fz2;



      /*        
        fx3 =  fx2;
        fy3 =  fy2;
        fz3 =  fz2;
      */
        in01 =  hmfx + kmfy + lmfz ;
	in02 =  hmfz + kmfx + lmfy ;
        in03 =  hmfy + kmfz + lmfx ;
	/*
        fx3 = -fx2;
        fy3 = -fy2;
        fz3 =  fz2;
	*/
        in04 = -hmfx - kmfy + lmfz ;
	in05 =  hmfz - kmfx - lmfy ;
        in06 = -hmfy + kmfz - lmfx ;
        /*
        fx3 = -fx2;
        fy3 =  fy2;
        fz3 = -fz2;
        */
        in07 = -hmfx + kmfy - lmfz ;
	in08 = -hmfz - kmfx + lmfy ;
        in09 =  hmfy - kmfz - lmfx ;
        /*
        fx3 =  fx2;
        fy3 = -fy2;
        fz3 = -fz2;
        */
        in10 =  hmfx - kmfy - lmfz ;
	in11 = -hmfz + kmfx - lmfy ;
        in12 = -hmfy - kmfz + lmfx ;


        fmatR[m][n] =  fmatR[m][n] + 
                       cos(TWOPI*in01) +
                       cos(TWOPI*in02) + 
	               cos(TWOPI*in03) +
                       cos(TWOPI*in04) +
                       cos(TWOPI*in05) + 
	               cos(TWOPI*in06) +
                       cos(TWOPI*in07) +
                       cos(TWOPI*in08) + 
	               cos(TWOPI*in09) +
                       cos(TWOPI*in10) +
                       cos(TWOPI*in11) + 
	               cos(TWOPI*in12) ;



        fmatI[m][n] =  fmatI[m][n] + 
                       sin(TWOPI*in01) +
                       sin(TWOPI*in02) +
                       sin(TWOPI*in03) +
                       sin(TWOPI*in04) +
                       sin(TWOPI*in05) +
                       sin(TWOPI*in06) +
                       sin(TWOPI*in07) +
                       sin(TWOPI*in08) +
                       sin(TWOPI*in09) +
                       sin(TWOPI*in10) +
                       sin(TWOPI*in11) +
	               sin(TWOPI*in12) ;

   }

    ///    if(fabs(fmatR[m][n])<0.001)fmatR[m][n] == 0.;
    ///if(fabs(fmatI[m][n])<0.001)fmatI[m][n] == 0.;

  }
  ///   if(m/100 == 0)printf("  %d  ",m);
  }
  ///  printf("\n");

  ///    for (i=0;i<10;i++) for (j=0;j<10;j++)  printf( " %d %d %f %f \n ", i,j,fmatI[i][j],fmatR[i][j]  );



  scale = 30.0;
  if (strcmp(argv[4],"test") != 0   ) {
  
    ///  scale=1.05-1.0/200.*2;
    double rfmin = 1000, minscale;

    
  for (iii=0;iii<200;iii++){
    scale=scale+0.1;
  rf=rfactor(N,occ,grad,f_data);
  if (rf < rfmin ){
    rfmin = rf;
    minscale = scale;
  }
  }
  scale=minscale;

  logfc = 1;
  rf=rfactor(N,occ,grad,f_data);
  printf("original rfac %f  at scale = %f \n ",rfmin,scale );  
  logfc = 0;
  

  




  nlopt_opt opt1;
    opt1=nlopt_create(NLOPT_LN_NELDERMEAD,N);
  ///  opt1=nlopt_create(NLOPT_LD_LBFGS,N);
  /// opt1=nlopt_create(NLOPT_LD_TNEWTON,N);
  double lb = -0.50000;
  double ub = 1.00000;
  nlopt_set_lower_bounds1(opt1,lb);  
  nlopt_set_upper_bounds1(opt1,ub);
  nlopt_set_min_objective(opt1, rfactor, f_data) ; 
  nlopt_set_xtol_rel(opt1,1e-5);

  double flb[1],fub[1];
  nlopt_algorithm algo=nlopt_get_algorithm(opt1);
  nlopt_get_lower_bounds(opt1,flb);  
  nlopt_get_upper_bounds(opt1,fub);  
  printf("limit %f <  %f  ? \n",flb[0],fub[0]);
  printf("algorithm %s \n ",nlopt_algorithm_name( algo ) );  
  printf("dimetion %d  \n",nlopt_get_dimension(opt1) );

  ires=nlopt_optimize(opt1,occ,&minf);

  printf("result %d \n",ires);

  nlopt_destroy(opt1);
  

  }

  logfc = 2;
  rf=rfactor(N,occ,grad,f_data);
  logfc = 0;

  printf("new rfac %f  minf %f \n ",rf,minf );  
 


  for (i=0;i<M;i++)free(fmatR[i]);
  free(fmatR);
  for (i=0;i<M;i++)free(fmatI[i]);
  free(fmatI);



  fpo=fopen(argv[3],"w");
  for(n=0;n<N;n++){
    strncpy(dumm55,pdbl[n],55);
    strncpy(dtail60a,pdbl[n]+60,20);
    sprintf(wline,"%s% 1.2f%s",dumm55,occ[n],dtail60a);

    fputs(wline,fpo);

  }
  fclose(fpo);

  return 0;




}

double rfactor(int N, double *occ, double *grad, void *f_data){

  int m,n;
  double FCR[M],FCI[M],FC[M];
  double diff,sumdiff;
  FILE *fplog;

     if ( logfc == 1 )       fplog = fopen("fcdata","w");
     else if ( logfc == 2)   fplog = fopen("fcdata2","w");

  sumdiff=0.;
  for (m=0;m<M;m++){
    FCR[m] = 0.; FCI[m] = 0.;
  for (n=0;n<N;n++){
    FCR[m] = FCR[m] + fmatR[m][n] * occ[n];
    FCI[m] = FCI[m] + fmatI[m][n] * occ[n];
  }
    FC[m] = sqrt( FCR[m]*FCR[m] + FCI[m]*FCI[m] );
    diff = fabs(F[m] - scale * FC[m]);
    sumdiff = sumdiff + diff;

     if ( logfc > 0 )  fprintf(fplog,"Fo = %f  FC = %f  FCR, %f FCI %f \n", F[m],FC[m],FCR[m],FCI[m] );
  }

  if ( logfc > 0 )  fclose(fplog);


  ntrial++;
  if ( ntrial%1 == 0) printf("inside rf = %f  at trial %d  scale %f \n",sumdiff/sumF,ntrial,scale);
  return sumdiff/sumF;
}
