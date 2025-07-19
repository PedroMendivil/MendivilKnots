//Libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Alex_pol.h"

//Defines

#define M_PI 3.14159265358979323846

//Structs

struct sim_params
{
  float T;
  float eta;
  int N;
  float f_pull;
  char IC_type;
  float meas_freq;
  float sim_time;
};

struct rand_var_props
{
  double m_1;
  double m_2;
  double var;
  double sem;
  int n_term;
  int term_end;
  int n_blocks;
};

typedef struct sim_params sim_params;
typedef struct rand_var_props rand_var_props;

//Functions

void read_parameters( sim_params *sp, FILE *f)
{
  fscanf(f,"T\t%f\n",&(sp->T));
  fscanf(f,"eta\t%f\n",&(sp->eta));
  fscanf(f,"N\t%d\n",&(sp->N));
  fscanf(f,"f_pull\t%f\n",&(sp->f_pull));
  fscanf(f,"IC_type\t%c\n",&(sp->IC_type));
  fscanf(f,"meas_freq\t%f\n",&(sp->meas_freq));
  fscanf(f,"sim_time\t%f\n",&(sp->sim_time));
}

int open_trr_file( char *sim_dir, int f_idx, FILE **f)
{
  char filename[256];
  snprintf(filename,sizeof(filename),"%s/trajectory-file-%d.trr",sim_dir,f_idx);
  *f=fopen(filename,"rb");
  if( *f==NULL)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

int read_trr_frame( int N, float *r, float *t, FILE *f)
{
  //header
  int magickvalue;
  if( fread(&magickvalue,sizeof(int),1,f)==0){ return 0;}
  if( magickvalue!=1993){ printf("Wrong format.\n"); exit(-1);}
  char trrversion[13];
  int len_s_a;
  int len_s_b;
  fread(&len_s_a,sizeof(int),1,f);
  fread(&len_s_b,sizeof(int),1,f);
  fread(trrversion,sizeof(char),sizeof(trrversion)-1,f);
  int zero;
  for( int i=0; i<7; i++)
  {
    fread(&zero,sizeof(int),1,f);
  }
  int x_size;
  fread(&x_size,sizeof(int),1,f);
  if( x_size!=3*N*sizeof(float)){ printf("Error reading trajectory.\n"); exit(-1);}
  int v_size;
  fread(&v_size,sizeof(int),1,f);
  if( v_size!=0){ printf("Error reading trajectory.\n"); exit(-1);}
  int f_size;
  fread(&f_size,sizeof(int),1,f);
  if( f_size!=0){ printf("Error reading trajectory.\n"); exit(-1);}
  int natoms;
  fread(&natoms,sizeof(int),1,f);
  if( natoms!=N){ printf("Error reading trajectory.\n"); exit(-1);}
  int step;
  fread(&step,sizeof(int),1,f);
  int time;
  fread(&zero,sizeof(int),1,f);
  fread(&time,sizeof(float),1,f);
  fread(&zero,sizeof(int),1,f);
  *t=time;
  //coordinates
  fread(r,sizeof(float),3*N,f);
  return 1;
}

void calc_cm( int N, float *r, double *r_cm)
{
  r_cm[0]=r_cm[1]=r_cm[2]=0.0;
  for( int i_p=0; i_p<N; i_p++)
  {
    r_cm[0] += r[3*i_p+0];
    r_cm[1] += r[3*i_p+1];
    r_cm[2] += r[3*i_p+2];
  }
  r_cm[0] /= N;
  r_cm[1] /= N;
  r_cm[2] /= N;
}

void calc_Rg2( int N, float *r, double *r_cm, double *Rg2)
{
  *Rg2 = 0.0;
  for( int i_p=0; i_p<N; i_p++)
  {
    *Rg2 += (r[3*i_p+0]-r_cm[0])*(r[3*i_p+0]-r_cm[0]);
    *Rg2 += (r[3*i_p+1]-r_cm[1])*(r[3*i_p+1]-r_cm[1]);
    *Rg2 += (r[3*i_p+2]-r_cm[2])*(r[3*i_p+2]-r_cm[2]);
  }
  *Rg2 = *Rg2/N;
}

void calc_nce( int N, float *r, double *nce, double theta, double varphi)
{
  double f_dir[3];
  f_dir[0] = sin(theta)*cos(varphi);
  f_dir[1] = sin(theta)*sin(varphi);
  f_dir[2] = cos(theta);
  *nce = 0.0;
  for( int i_c=0; i_c<3; i_c++)
  {
    *nce += (r[3*(N-1)+i_c]-r[3*0+i_c])*f_dir[i_c];
  }
  *nce /= (N-1.0);
}

void block_data( int *n_data, double *x)
{
  int i;
  for( i=0; (2*i)<(*n_data-1); i++)
  {
    x[i] = 0.5*(x[2*i]+x[2*i+1]);
  }
  *n_data = i;
}

void calc_moments( int n_data, double *x, double *m_1, double *m_2)
{
  *m_1 = 0.0;
  *m_2 = 0.0;
  for( int i=0; i<n_data; i++)
  {
    *m_1 += x[i];
    *m_2 += x[i]*x[i];
  }
  *m_1 /= n_data;
  *m_2 /= n_data;
}

int calc_sig_mean( int n_data, double *x, double *sem)
{
  int n_data_b = n_data;
  double m_1, m_2, var_m_1, uplim_var_m_1;
  calc_moments(n_data_b,x,&m_1,&m_2);
  var_m_1 = (m_2-m_1*m_1)/(n_data_b-1.0);
  uplim_var_m_1 = var_m_1*(1.0+sqrt(2.0/(n_data_b-1.0)));
  while( n_data_b>3)
  {
    block_data(&n_data_b,x);
    calc_moments(n_data_b,x,&m_1,&m_2);
    var_m_1 = (m_2-m_1*m_1)/(n_data_b-1.0);
    if( var_m_1>uplim_var_m_1)
    {
      uplim_var_m_1 = var_m_1*(1.0+sqrt(2.0/(n_data_b-1.0)));
    }
    else
    {
      *sem = sqrt(var_m_1);
      return n_data_b;
    }
  }
  *sem = sqrt(var_m_1);
  return n_data_b; 
}

int estimate_term( int n_data, double *x, int *opt_n_term)
{
  int i_term;
  int n_term;
  double m_1, m_2, smer;
  double smer_min = INFINITY;
  for( i_term=0; i_term<50; i_term++)
  {
    n_term = i_term*n_data/100;
    calc_moments(n_data-n_term,&x[n_term],&m_1,&m_2);
    smer = (m_2-m_1*m_1)/(n_data-n_term);
    if( smer<smer_min)
    {
      *opt_n_term = n_term;
      smer_min = smer;
    }
  }
  if( *opt_n_term==n_term)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

int main( int argc, char const *argv[])
{
  if( argc!=2){
    if( argc<2){ printf("You forgot the input.\n"); exit(-1);}
    else{ printf("Too many arguments.\n"); exit(-1);}
  }

  if( sizeof(argv[1])>128){ printf("Directory name too long.\n"); exit(-1);}
  char sim_dir[128];
  snprintf(sim_dir,sizeof(sim_dir),"%s",argv[1]);

  FILE *file_in;
  FILE *file_o1;
  FILE *file_o2;

  char filename[256];

  //Simulation parameters and variables

  sim_params sp;

  snprintf(filename,sizeof(filename),"%s/parameters.dat",sim_dir);
  file_in=fopen(filename,"rt");
  if( file_in==NULL){ printf("Error opening parameters file.\n"); exit(-1);}
  read_parameters(&sp,file_in);
  fclose(file_in);

  float *r;
  r=(float*)malloc(3*sp.N*sizeof(float));

  //Statistical variables calculation

  float t;

  int f_idx=0;
  int i_f=0;

  snprintf(filename,sizeof(filename),"%s/time-series.dat",sim_dir);
  file_o1=fopen(filename,"wt");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  snprintf(filename,sizeof(filename),"%s/knot-type.dat",sim_dir);
  file_o2=fopen(filename,"wt");
  if( file_o2==NULL){ printf("Error opening file.\n"); exit(-1);}

  while( open_trr_file(sim_dir,f_idx,&file_in))
  {
    double r_cm[3];
    double Rg2;
    double nce;
    while( read_trr_frame(sp.N,r,&t,file_in)) 
    {
      calc_cm(sp.N,r,r_cm);
      calc_Rg2(sp.N,r,r_cm,&Rg2);
      calc_nce(sp.N,r,&nce,M_PI/2.0,0.0);
      fprintf(file_o1,"%15.6f%15.6lf%15.9lf\n",t,Rg2,nce);
      i_f++;
    }
    fclose(file_in);
    char knot_name[128];
    calc_invariant(sp.N,r,knot_name);
    fprintf(file_o2,"%s\n",knot_name);
    f_idx++;
  }

  fclose(file_o1);

  fclose(file_o2);

  free(r);

  int n_files=f_idx;
  if( n_files==0){ printf("No trajectory files.\n"); exit(-1);}

  int n_data=i_f;
  if( n_data==0){ printf("No data.\n"); exit(-1);}

  //Time series reading

  double *Rg2;
  double *nce;
  Rg2=(double*)malloc(n_data*sizeof(double));
  nce=(double*)malloc(n_data*sizeof(double));

  snprintf(filename,sizeof(filename),"%s/time-series.dat",sim_dir);
  file_in=fopen(filename,"rt");
  if( file_in==NULL){ printf("Error opening file.\n"); exit(-1);}

  for( i_f=0; i_f<n_data; i_f++)
  {
    fscanf(file_in,"%15f%15lf%15lf\n",&t,&Rg2[i_f],&nce[i_f]);
  }

  fclose(file_in);

  //Analysis of the statisical variables

  rand_var_props Rg2_p;
  rand_var_props nce_p;

  //Termalisation time estimation (marginal standard error rule)

  Rg2_p.term_end=estimate_term(n_data,Rg2,&Rg2_p.n_term);
  nce_p.term_end=estimate_term(n_data,nce,&nce_p.n_term);

  //Estimation of the first two moments

  calc_moments(n_data-Rg2_p.n_term,&Rg2[Rg2_p.n_term],&Rg2_p.m_1,&Rg2_p.m_2);
  Rg2_p.var = (Rg2_p.m_2-Rg2_p.m_1*Rg2_p.m_1);
  Rg2_p.var *= (n_data-Rg2_p.n_term)/(n_data-Rg2_p.n_term-1.0);
  calc_moments(n_data-nce_p.n_term,&nce[nce_p.n_term],&nce_p.m_1,&nce_p.m_2);
  nce_p.var = (nce_p.m_2-nce_p.m_1*nce_p.m_1);
  nce_p.var *= (n_data-nce_p.n_term)/(n_data-nce_p.n_term-1.0);

  //Estimation of the standard deviation of the sample mean

  Rg2_p.n_blocks=calc_sig_mean(n_data-Rg2_p.n_term,&Rg2[Rg2_p.n_term],&Rg2_p.sem);
  nce_p.n_blocks=calc_sig_mean(n_data-nce_p.n_term,&nce[nce_p.n_term],&nce_p.sem);

  //Results

  snprintf(filename,sizeof(filename),"%s/statistics.dat",sim_dir);
  file_o1=fopen(filename,"wt");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  fprintf(file_o1,"#Analysis of %d files from %s.\n",n_files,sim_dir);
  fprintf(file_o1,"#%9s%10s%10s%10s\n","N","T","eta","f_pull");
  fprintf(file_o1,"%10d%10.4f%10.4f%10.4f\n",sp.N,sp.T,sp.eta,sp.f_pull);

  if( Rg2_p.n_blocks<32 || !Rg2_p.term_end)
  {
    fprintf(file_o1,"#Warning: Not enough data to analyse Rg2.\n");
  }
  else
  {
    fprintf(file_o1,"\n");
  }
  fprintf(file_o1,"#Rg2: n_term = %d, n_blocks = %d\n",Rg2_p.n_term,Rg2_p.n_blocks);
  fprintf(file_o1,"#%14s%15s%15s\n","<Rg2>","SEM(Rg2)","Var(Rg2)");
  fprintf(file_o1,"%15.6lf%15.6lf%15.6lf\n",Rg2_p.m_1,Rg2_p.sem,Rg2_p.var);

  if( nce_p.n_blocks<32 || !nce_p.term_end)
  {
    fprintf(file_o1,"#Warning: Not enough data to analyse nce.\n");
  }
  else
  {
    fprintf(file_o1,"\n");
  }
  fprintf(file_o1,"#nce: n_term = %d, n_blocks = %d\n",nce_p.n_term,nce_p.n_blocks);
  fprintf(file_o1,"#%14s%15s%15s\n","<nce>","SEM(nce)","Var(nce)");
  fprintf(file_o1,"%15.9lf%15.9lf%15.9lf\n",nce_p.m_1,nce_p.sem,nce_p.var);

  fclose(file_o1);

  printf("Analysed %d files from %s.\n",n_files,sim_dir);

  return 0;
}
