///////////////////////////////////////////////////////////////////////
// PROGRAMA TFG DE MARCO MENDIVIL                                    // 
// CALCULA LA EVOLUCION DE CADENAS DE POLIMEROS                      //
// MODIFICACION PARA EL ANALISIS DE NUDOS PEDRO MENDIVIL Y CHAT GPT  //
// VERSION V1.X de 09 de noviembre de 2024                           // 
///////////////////////////////////////////////////////////////////////

//Libraries  /////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <curand_kernel.h>
#include <string.h>

// Defines  ///////////////////////////////////////////////////////////

#define M_PI 3.14159265358979323846

#define h 0.001     // h = incremento de tiempo (delta t)
#define m 1.0       // m = masa
#define k_e 100.0   // k_e = constante elástica
#define b_pol 1.0   // b_pol = distancia a la próxima partícula 
#define k_b 0.0     // k_b = constante de doblado
#define sig 1.5     // sig = sigma (de potencial Lennar-Jones)
#define eps 3.0     // eps = epsilon (de potencial Lennar-Jones)
#define cutoff pow(2,1.0/6.0)*sig // = cutoff del potencial Lennar-Jones

// Variables globales  ///////////////////////////////////////////////

// npartic = número de partículas. Para random usamos global_N.
int npartic; // se modifica con cada xx.gro leido.
// file_copy[100] array global, copia de file_to_copy de void read_parameters()
char file_copy[100]; // cadena que contiene el nomre del fichero de nudo
// Todas las variables de parameters a global variables
float global_T;
float global_eta;
int global_N;
int global_N_random;
float global_f_pull;
char global_IC_type;
float global_meas_freq;
float global_sim_time;

// Structs   ///////////////////////////////////////////////////////

struct sim_params
{
  float T;           //          T   = Temperature
  float eta;         //         eta  = Viscosity
  int N;             //          N   = Number of particules (only for random)
  float f_pull;      //       f_pull = Pull force
  char IC_type;      //      IC_type = Initial condition type (f!!,r o l (load)) 
  // f do not use(checkpoint), r => chain random, l => load file XX.gro (knot)
  float meas_freq;   //    meas_freq = Frecuency of measurement (frames)
  float sim_time;    //     sim_time = Time of simulation
  char file_to_load[100]; // file_to_load = File format .gro to be read
};

typedef struct sim_params sim_params;
typedef struct curandStatePhilox4_32_10 PRNGstate;

//Functions   ///////////////////////////////////////////////////////////////

/////////////////// read_parameters & write global_vars  ////////////////////

void read_parameters( struct sim_params *sp, FILE *f)
{
  fscanf(f,"T\t%f\n",&(sp->T));
  fscanf(f,"eta\t%f\n",&(sp->eta));
  fscanf(f,"N\t%d\n",&(sp->N));
  fscanf(f,"f_pull\t%f\n",&(sp->f_pull));
  fscanf(f,"IC_type\t%c\n",&(sp->IC_type));
  fscanf(f,"meas_freq\t%f\n",&(sp->meas_freq));
  fscanf(f,"sim_time\t%f\n",&(sp->sim_time));
   
  if (fscanf(f, "file_to_load\t%99s\n", sp->file_to_load) != 1) {
      fprintf(stderr, "Error al leer file_to_load\n");
      exit(EXIT_FAILURE);
  }
  // Copiamos el contenido de file_to_load en file_copy que es global
  strcpy(file_copy, sp->file_to_load); // Copia el contenido de file_to_load en file_copy
  
  // Copiamos los demás valores de la estructura a las variables globales
  global_T = sp->T;
  global_eta = sp->eta;
  global_N = sp->N; // Redefinir inmediatamente con make_global_N
  global_f_pull = sp->f_pull;
  global_IC_type = sp->IC_type;
  global_meas_freq = sp->meas_freq;
  global_sim_time = sp->sim_time;
} 
//////////////////////////////// make_global_N ///////////////////////

int make_global_N() 
{
  FILE *knot = fopen(file_copy, "r");
  if (knot == NULL) {
     fprintf(stderr, "Error en la apertura del archivo de nudo\n");
     exit(EXIT_FAILURE);
  }
  char line[256];
  fgets(line, sizeof(line), knot); // Primera linea del archivo (MD simul etc)
  fgets(line, sizeof(line), knot); // Segunda linea del archivo (num partic. )
  int npartic = atoi(line); // Número de partículas en la segunda linea 
  fclose(knot);
    if (global_IC_type == 'k') { global_N = npartic; }  
  return global_N; // Devuelve global_N para uso posterior...
} 
////////////////////////// print parameters ////////////////////////////

void print_parameters(struct sim_params *sp, FILE *f)
{  
  printf("\n\n\n");
  printf("         PARAMETERS:\n");
  printf("\n%10s", "");
  printf(" T = %06.3f\n         eta = %06.3f\n           N = %04d\n      f_pull = %06.3f\n",sp->T,sp->eta,sp->N,sp->f_pull);
  printf("%4s", "");
  printf(" IC_type = %c\n   meas_freq = %06.3f\n    sim_time = %08.1f\nfile_to_load = %s\n",sp->IC_type,sp->meas_freq,sp->sim_time,sp->file_to_load);
} 
////////////////////////// write_gro_frame //////////////////////////////

void write_gro_frame( int global_N, float *r, FILE *f)
{
  fprintf(f,"MD simulation of a knot of polymer, t = 0.0\n");
  fprintf(f,"%5d\n",global_N);
  for( int i_p=0; i_p<global_N; i_p++)
  {
    fprintf(f,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",i_p+1,"X","X",i_p+1,r[3*i_p+0],r[3*i_p+1],r[3*i_p+2]);
  }
  fprintf(f,"%10.5f%10.5f%10.5f\n",0.0,0.0,0.0);
}
////////////////////////// write_trr_frame //////////////////////////////

void write_trr_frame( int global_N, float *r, int i_f, float t, FILE *f)
{
  //header
  int magickvalue=1993;
  fwrite(&magickvalue,sizeof(int),1,f);
  char trrversion[]="GMX_trn_file";
  int len_s_a=sizeof(trrversion);
  int len_s_b=sizeof(trrversion)-1;
  fwrite(&len_s_a,sizeof(int),1,f);
  fwrite(&len_s_b,sizeof(int),1,f);
  fwrite(trrversion,sizeof(char),sizeof(trrversion)-1,f);
  int zero=0;
  for( int i=0; i<7; i++)
  {
    fwrite(&zero,sizeof(int),1,f);
  }
  int x_size=3*global_N*sizeof(float);
  fwrite(&x_size,sizeof(int),1,f);
  int v_size=0;
  fwrite(&v_size,sizeof(int),1,f);
  int f_size=0;
  fwrite(&f_size,sizeof(int),1,f);
  int natoms=global_N;
  fwrite(&natoms,sizeof(int),1,f);
  int step=i_f;
  fwrite(&step,sizeof(int),1,f);
  int time=t;
  fwrite(&zero,sizeof(int),1,f);
  fwrite(&time,sizeof(float),1,f);
  fwrite(&zero,sizeof(int),1,f);
  //coordinates
  fwrite(r,sizeof(float),3*global_N,f);
}
////////////////////////// write_checkpoint //////////////////////////////

void write_checkpoint( int global_N, float *r, float *v, float t, int n_threads, PRNGstate *state, int f_idx, FILE *f)
{
  int natoms=global_N;
  fwrite(&natoms,sizeof(int),1,f);
  int time=t;
  fwrite(&time,sizeof(float),1,f);
  int index=f_idx;
  fwrite(&index,sizeof(int),1,f);
  fwrite(r,sizeof(float),3*global_N,f);
  fwrite(v,sizeof(float),3*global_N,f);
  fwrite(state,sizeof(PRNGstate),n_threads,f);
}

void read_checkpoint( int global_N, float *r, float *v, float *t, int n_threads, PRNGstate *state, int *f_idx, FILE *f)
{
  int natoms;
  fread(&natoms,sizeof(int),1,f);
  if( natoms!=global_N){ printf("Error reading checkpoint.\n"); exit(-1);}
  int time;
  fread(&time,sizeof(float),1,f);
  int index;
  fread(&index,sizeof(int),1,f);
  fread(r,sizeof(float),3*global_N,f);
  fread(v,sizeof(float),3*global_N,f);
  fread(state,sizeof(PRNGstate),n_threads,f);
  *t=time;
  *f_idx=index+1;
}
////////////////////////// set_random_IC //////////////////////////////

void set_random_IC( int global_N, float T, float *r, float *v)
{
  curandGenerator_t gen;
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(NULL));
  float random, theta, varphi, bondlen, bondangle;
  float dir_old[3], dir_new[3], perpdir[3], perpdirnorm;
  curandGenerateUniform(gen,&random,1); theta = acos(1.0-2.0*random); 
  curandGenerateUniform(gen,&random,1); varphi = 2.0*M_PI*random;
  dir_old[0]=sin(theta)*cos(varphi);
  dir_old[1]=sin(theta)*sin(varphi);
  dir_old[2]=cos(theta);
  r[0]=r[1]=r[2]=0.0;
  for( int i_p=1; i_p<global_N; i_p++)
  {
    curandGenerateUniform(gen,&random,1); theta = acos(1.0-2.0*random); 
    curandGenerateUniform(gen,&random,1); varphi = 2.0*M_PI*random;
    perpdir[0] = dir_old[1]*cos(theta)-dir_old[2]*sin(theta)*sin(varphi);
    perpdir[1] = dir_old[2]*sin(theta)*cos(varphi)-dir_old[0]*cos(theta);
    perpdir[2] = dir_old[0]*sin(theta)*sin(varphi)-dir_old[1]*sin(theta)*cos(varphi);
    perpdirnorm = sqrt(perpdir[0]*perpdir[0]+perpdir[1]*perpdir[1]+perpdir[2]*perpdir[2]);
    perpdir[0] /= perpdirnorm; perpdir[1] /= perpdirnorm; perpdir[2] /= perpdirnorm;
    curandGenerateUniform(gen,&random,1);
    if( k_b<__FLT_MIN__){ bondangle = acos(1.0-2.0*random);}
    else{ bondangle = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b/T)))*random)/(k_b/T));}
    dir_new[0] = dir_old[0]*cos(bondangle)+perpdir[0]*sin(bondangle);
    dir_new[1] = dir_old[1]*cos(bondangle)+perpdir[1]*sin(bondangle);
    dir_new[2] = dir_old[2]*cos(bondangle)+perpdir[2]*sin(bondangle);
    curandGenerateUniform(gen,&random,1);
    bondlen = b_pol+sqrt(2.0)*sqrt(T/k_e)*erfinv(2.0*random-1.0);
    r[3*i_p+0] = bondlen*dir_new[0]+r[3*(i_p-1)+0];
    r[3*i_p+1] = bondlen*dir_new[1]+r[3*(i_p-1)+1];
    r[3*i_p+2] = bondlen*dir_new[2]+r[3*(i_p-1)+2];
    int overlap = 0;
    for(int j_p=0; j_p<(i_p-1); j_p++)
    {
      float dist = 0.0;
      dist += (r[3*j_p+0]-r[3*i_p+0])*(r[3*j_p+0]-r[3*i_p+0]);
      dist += (r[3*j_p+1]-r[3*i_p+1])*(r[3*j_p+1]-r[3*i_p+1]);
      dist += (r[3*j_p+2]-r[3*i_p+2])*(r[3*j_p+2]-r[3*i_p+2]);
      dist = sqrt(dist);
      if( dist<cutoff){ overlap=1; break;}
    }
    if( overlap==0)
    {
      dir_old[0] = dir_new[0];
      dir_old[1] = dir_new[1];
      dir_old[2] = dir_new[2];
    }
    else{ i_p--;}
  }
  for( int i_p=0; i_p<global_N; i_p++)
  {
    for( int i_c=0; i_c<3; i_c++)
    {
      curandGenerateUniform(gen,&random,1);
      v[3*i_p+i_c] = sqrt(2.0)*sqrt(T/m)*erfinv(2.0*random-1.0);
    }
  }
  curandDestroyGenerator(gen);
}
////////////////////////// set_knot_IC //////////////////////////////

void set_knot_IC( int global_N, float *r, float *v)
{
  FILE *knot = fopen(file_copy, "r"); 
  if (knot == NULL) {
      fprintf(stderr, "Error en la apertura del archivo de nudo\n");
      exit(EXIT_FAILURE); // Salir en caso de error de apertura
  }
  float x=0;// float x, y, z; // Variables de las posiciones de los atomos
  float y=0;
  float z=0;
  char line[256]; // Los archivos .gro tendran lineas de menos de 256 caracteres
  
  // Leemos las dos primeras lineas del fichero de polimero que no nos interesan
  fgets(line, sizeof(line), knot);
  printf(" \n%s",line);
  fgets(line, sizeof(line), knot);
  npartic = atoi(line); //Pasa el valor de line a entero y lo asigna a npartic
  printf("Number of particles = %d\n", npartic);
 
  for( int i_p=0; i_p<npartic; i_p++) 
  {
    // Leer las npartic-2 siguientes lineas con las posiciones de los atomos
    fgets(line, sizeof(line), knot);
    // Usamos sscanf para extraer los 3 últimos valores de la línea
    sscanf(line, "%*s %*s %*s %f %f %f", &x, &y, &z);
   
    r[3*i_p+0] = x; // La coordenada x leida la metemos en el array r[ ]
    r[3*i_p+1] = y; // La coordenada y leida la metemos en el array r[ ]
    r[3*i_p+2] = z; // La coordenada z leida la metemos en el array r[ ]
  }
  
  for( int i_p=0; i_p<npartic; i_p++)
  {
    for( int i_c=0; i_c<3; i_c++)
    {
      v[3*i_p+i_c] = 0.0;
    }
  }
  fclose(knot); // Cerrar el archivo correctamente
} 


//Kernels: se ejecutan en la GPU ///////////////////////////////////////

__global__
void setup_PRNG( int seed, PRNGstate *state)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(seed, i_p, 0, &state[i_p]);
}

__global__
void call_PRNG( float c_rn, float *nrn, PRNGstate *state)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  for( int i_c = 0; i_c < 3; i_c++)
  {
    nrn[3*i_p+i_c] = c_rn*curand_normal(&state[i_p]);
  }
}

__global__
void calc_extern_f( int global_N, float eta, float *v, float *f_c, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<global_N)
  {
    for( int i_c = 0; i_c < 3; i_c++)
    {
      f[3*i_p+i_c] = f_c[3*i_p+i_c];

      f[3*i_p+i_c] += -eta*v[3*i_p+i_c];
    }
  }
}

__global__ 
void calc_bonds( int global_N, float *r, float *b, float *invlen)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<global_N-1)
  {
    invlen[i_p+2] = 0.0;
    for( int i_c = 0; i_c < 3; i_c++)
    {
      b[3*(i_p+2)+i_c] = r[3*(i_p+1)+i_c]-r[3*i_p+i_c];
      invlen[i_p+2] += b[3*(i_p+2)+i_c]*b[3*(i_p+2)+i_c];
    }
    invlen[i_p+2] = 1.0/sqrt(invlen[i_p+2]);
  }
}

__global__
void calc_cosines( int global_N, float *b, float *invlen, float *cosine)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<global_N-2)
  {
    cosine[i_p+3] = 0.0;
    for( int i_c = 0; i_c < 3; i_c++)
    {
      cosine[i_p+3] += b[3*(i_p+3)+i_c]*b[3*(i_p+2)+i_c];
    }
    cosine[i_p+3] *= invlen[i_p+3]*invlen[i_p+2];
  }
}

__global__
void calc_intern_f( int global_N, float *b, float *invlen, float *cosine, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<global_N)
  {
    for( int i_c = 0; i_c < 3; i_c++)
    {
      f[3*i_p+i_c] += k_e*(1.0-b_pol*invlen[i_p+1])*(-b[3*(i_p+1)+i_c]);

      f[3*i_p+i_c] += k_e*(1.0-b_pol*invlen[i_p+2])*(+b[3*(i_p+2)+i_c]);

      f[3*i_p+i_c] += k_b*(+b[3*(i_p+0)+i_c])*invlen[i_p+1]*invlen[i_p+0];

      f[3*i_p+i_c] += k_b*(+b[3*(i_p+2)+i_c]-b[3*(i_p+1)+i_c])*invlen[i_p+2]*invlen[i_p+1];

      f[3*i_p+i_c] += k_b*(-b[3*(i_p+3)+i_c])*invlen[i_p+3]*invlen[i_p+2];

      f[3*i_p+i_c] += k_b*(-cosine[i_p+2]-cosine[i_p+1])*b[3*(i_p+1)+i_c]*invlen[i_p+1]*invlen[i_p+1];

      f[3*i_p+i_c] += k_b*(+cosine[i_p+3]+cosine[i_p+2])*b[3*(i_p+2)+i_c]*invlen[i_p+2]*invlen[i_p+2];
    }
  }
}

__global__
void calc_exclvol_f( int global_N, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<global_N)
  {
    int skip;
    float d2;
    float k_LJ;
    for( int j_p=i_p-2; j_p>=0; j_p-=1+skip)
    {
      d2 = 0.0;
      for( int i_c=0; i_c<3; i_c++)
      {
        d2 += (r[3*i_p+i_c]-r[3*j_p+i_c])*(r[3*i_p+i_c]-r[3*j_p+i_c]);
      }
      if( d2<pow(cutoff,2))
      {
        k_LJ = 4.0*eps*(12.0*pow(sig,12)/(d2*d2*d2*d2*d2*d2*d2)-6.0*pow(sig,6)/(d2*d2*d2*d2));
        for( int i_c = 0; i_c < 3; i_c++)
        {
          f[3*i_p+i_c] += k_LJ*(r[3*i_p+i_c]-r[3*j_p+i_c]);
        }
        skip=0;
      }
      else
      {
        skip=((sqrt(d2)-cutoff)/(1.5*b_pol));
      }
    }
    for( int j_p=i_p+2; j_p<global_N; j_p+=1+skip)
    {
      d2 = 0.0;
      for( int i_c=0; i_c<3; i_c++)
      {
        d2 += (r[3*i_p+i_c]-r[3*j_p+i_c])*(r[3*i_p+i_c]-r[3*j_p+i_c]);
      }
      if( d2<pow(cutoff,2))
      {
        k_LJ = 4.0*eps*(12.0*pow(sig,12)/(d2*d2*d2*d2*d2*d2*d2)-6.0*pow(sig,6)/(d2*d2*d2*d2));
        for( int i_c = 0; i_c < 3; i_c++)
        {
          f[3*i_p+i_c] += k_LJ*(r[3*i_p+i_c]-r[3*j_p+i_c]);
        }
        skip=0;
      }
      else
      {
        skip=((sqrt(d2)-cutoff)/(1.5*b_pol));
      }
    }
  }
}

__global__
void RK_stage_1( int global_N, float eta, float *r_1, float *r_2, float *v_1, float *v_2, float *f_1, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<global_N)
  {
    for( int i_c = 0; i_c < 3; i_c++)
    {
      r_2[3*i_p+i_c] = r_1[3*i_p+i_c]+v_1[3*i_p+i_c]*h;

      v_2[3*i_p+i_c] = v_1[3*i_p+i_c]+f_1[3*i_p+i_c]*h/m+nrn[3*i_p+i_c]/m;
    }
  }
}

__global__
void RK_stage_2( int global_N, float eta, float *r_1, float *v_1, float *v_2, float *f_1, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<global_N)
  {
    for( int i_c = 0; i_c < 3; i_c++)
    {
      r_1[3*i_p+i_c] = r_1[3*i_p+i_c]+0.5*(v_1[3*i_p+i_c]+v_2[3*i_p+i_c])*h;

      v_1[3*i_p+i_c] = v_1[3*i_p+i_c]+0.5*(f_1[3*i_p+i_c]+f_2[3*i_p+i_c])*h/m+nrn[3*i_p+i_c]/m;
    }
  }
}

//////////////////////////    MAIN     //////////////////////////////

int main( int argc, char const *argv[])
{
///// argv[0]: nombre del programa y argv[1]: directorio de simulación (sim_dir) 

  if( argc!=2){  //si el número de argumentos es distinto de 2
    if( argc<2){ printf("You forgot the input.\n"); exit(-1);}//si es menor de 2
    else{ printf("Too many arguments.\n"); exit(-1);} // si mayor de 2 => exit
  }
  if( sizeof(argv[1])>128){ printf("Directory name too long.\n"); exit(-1);}
  char sim_dir[128];
  snprintf(sim_dir,sizeof(sim_dir),"%s",argv[1]);
 
  FILE *file_in;      // Puntero a archivo de entrada
  FILE *file_out;     // Puntero a archivo de salida
  char filename[256]; // Cadena para nombre de archivo
 
  // Simulation parameters and variables ////////////////////////////////////

  struct sim_params sp;
    FILE *param_file; // Declaramos el puntero a archivo
    param_file = fopen("/home/marcomc/Documentos/Program/Simulations/Test4/parameters.dat", "r");
    // Abrimos el archivo en modo lectura ("r") y asignamos el puntero
    if (param_file == NULL) {
        perror(" Error al abrir parameters.dat");
        return 1;
    }
  // Leer el archivo DEL DISCO: parameters.dat y meterlos en param_file
  read_parameters(&sp, param_file); // se leen los parámetros 
  print_parameters(&sp, stdout); // se imprimen los parámetros en pantalla
  fclose(param_file); // cerramos el archivo. 
 
  make_global_N(); // make_global_N hacia linea 90
  printf("global_N = number of particles = %d\n", global_N);
 
  float *r_1;
  float *r_2;

  float *v_1;
  float *v_2;

  float *f_1;
  float *f_2;

  float *nrn;
  PRNGstate *state;

  float c_rn = sqrt(2*sp.eta*sp.T*h);

  float *b;
  float *invlen;
  float *cosine;

  float *f_c;

  size_t threads_block = 256;
  size_t n_blocks = (global_N+threads_block-1)/threads_block;
  int n_threads = n_blocks*threads_block;

  // Memory allocation   //////////////////////////////////////////////

  cudaMallocManaged( &r_1, 3*global_N*sizeof(float));
  cudaMallocManaged( &r_2, 3*global_N*sizeof(float));

  cudaMallocManaged( &v_1, 3*global_N*sizeof(float));
  cudaMallocManaged( &v_2, 3*global_N*sizeof(float));

  cudaMallocManaged( &f_1, 3*global_N*sizeof(float));
  cudaMallocManaged( &f_2, 3*global_N*sizeof(float));

  cudaMallocManaged( &nrn, 3*n_threads*sizeof(float));
  cudaMallocManaged( &state, n_threads*sizeof(PRNGstate));

  cudaMallocManaged( &b, 3*(global_N+3)*sizeof(float));
  cudaMallocManaged( &invlen, (global_N+3)*sizeof(float));
  cudaMallocManaged( &cosine, (global_N+4)*sizeof(float));

  cudaMallocManaged( &f_c, 3*global_sim_time*sizeof(float));

  // Exceptions for the polymer ends ////////////////////////////////////

  b[0]=b[1]=b[2]=0.0;
  b[3]=b[4]=b[5]=0.0;
  b[3*(global_N+2)+0]=b[3*(global_N+2)+1]=b[3*(global_N+2)+2]=0.0;
  b[3*(global_N+1)+0]=b[3*(global_N+1)+1]=b[3*(global_N+1)+2]=0.0;

  invlen[0]=invlen[1]=0.0;
  invlen[global_N+2]=invlen[global_N+1]=0.0;

  cosine[0]=cosine[1]=cosine[2]=0.0;
  cosine[sp.N+3]=cosine[sp.N+2]=cosine[sp.N+1]=0.0;

  // Constant force  ////////////////////////////////////////////7

  for( int i_p=0; i_p<global_N; i_p++)
  {
    for( int i_c=0; i_c<3; i_c++)
    {
      f_c[3*i_p+i_c]=0.0;
    }
  }
  f_c[3*0+0] = -sp.f_pull;
  f_c[3*(global_N-1)+0] = sp.f_pull;

  // PRNG initialization  ////////////////////////////////////////

  setup_PRNG<<<n_blocks,threads_block>>>(time(NULL),state);
  cudaDeviceSynchronize();

  // Initial condition  ///////////////////////////////////////////

  float t = 0.0;

  int f_idx = 0;

  snprintf(filename,sizeof(filename),"%s/checkpoint.bin",sim_dir);
  file_in=fopen(filename,"rb");

  if( file_in==NULL)
  {
    if( sp.IC_type=='f')
    {
      snprintf(filename,sizeof(filename),"%s/initial-condition.bin",sim_dir);
      file_in=fopen(filename,"rb");
      if( file_in==NULL){ printf("Error opening file.\n"); exit(-1);}
      read_checkpoint(global_N,r_1,v_1,&t,n_threads,state,&f_idx,file_in);
      fclose(file_in);
      t = 0.0;
      f_idx = 0;
    }
    else if( sp.IC_type=='r')
    {
      set_random_IC(global_N,sp.T,r_1,v_1);
    }
    else if( sp.IC_type=='k')
    {
      set_knot_IC(npartic,r_1,v_1);
    }
    else
    {
      printf("Unknown IC.\n");
      exit(-1);
    }
  }
  else
  {
    read_checkpoint(global_N,r_1,v_1,&t,n_threads,state,&f_idx,file_in);
    fclose(file_in);
  }

  // Simulation  ////////////////////////////////////////////////////
  
  if( f_idx==0)
  {
    snprintf(filename,sizeof(filename),"%s/initial-condition.gro",sim_dir);
    file_out=fopen(filename,"wt");
    if( file_out==NULL){ printf("Error opening file.\n"); exit(-1);}
    write_gro_frame(global_N,r_1,file_out);
    fclose(file_out);
  }

  if( f_idx < 10)  
  {   
    snprintf(filename, sizeof(filename), "%s/initial-condition%d.gro", sim_dir, f_idx);
    file_out = fopen(filename, "wt");
    if (file_out == NULL) {
    printf("Error opening file %s.\n", filename);
    exit(-1);
    }
    write_gro_frame(global_N, r_1, file_out);
    fclose(file_out);
  }
    
  snprintf(filename,sizeof(filename),"%s/trajectory-file-%d.trr",sim_dir,f_idx);
  file_out=fopen(filename,"wb");
  if( file_out==NULL){ printf("Error opening file.\n"); exit(-1);}

  int i_s;
  int n_steps = round((1.0/sp.meas_freq)/h);
  int i_f;
  int n_frames = round(sp.sim_time*sp.meas_freq);

  for( i_f = 0; i_f < n_frames; i_f++)
  {
    printf("   Progress:%05.1lf%%\r",(100.0*i_f)/(1.0*n_frames));
    fflush(stdout);

    for( i_s = 0; i_s < n_steps; i_s++)
    {
      call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

      calc_extern_f<<<n_blocks,threads_block>>>(global_N,sp.eta,v_1,f_c,f_1);
      calc_bonds<<<n_blocks,threads_block>>>(global_N,r_1,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(global_N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(global_N,b,invlen,cosine,f_1);
      calc_exclvol_f<<<n_blocks,threads_block>>>(global_N,r_1,f_1);

      RK_stage_1<<<n_blocks,threads_block>>>(global_N,sp.eta,r_1,r_2,v_1,v_2,f_1,nrn);

      calc_extern_f<<<n_blocks,threads_block>>>(global_N,sp.eta,v_2,f_c,f_2);
      calc_bonds<<<n_blocks,threads_block>>>(global_N,r_2,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(global_N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(global_N,b,invlen,cosine,f_2);
      calc_exclvol_f<<<n_blocks,threads_block>>>(global_N,r_2,f_2);

      RK_stage_2<<<n_blocks,threads_block>>>(global_N,sp.eta,r_1,v_1,v_2,f_1,f_2,nrn);

      cudaDeviceSynchronize();
    }

    t+=n_steps*h;

    write_trr_frame(global_N,r_1,i_f,t,file_out);
  }

  fclose(file_out);

  snprintf(filename,sizeof(filename),"%s/checkpoint.bin",sim_dir);
  file_out=fopen(filename,"wb");
  if( file_out==NULL){ printf("Error opening file.\n"); exit(-1);}
  write_checkpoint(global_N,r_1,v_1,t,n_threads,state,f_idx,file_out);
  fclose(file_out);

  // Memory deallocation  ///////////////////////////////////////////

  cudaFree(r_1);
  cudaFree(r_2);

  cudaFree(v_1);
  cudaFree(v_2);

  cudaFree(f_1);
  cudaFree(f_2);

  cudaFree(nrn);
  cudaFree(state);

  cudaFree(b);
  cudaFree(invlen);
  cudaFree(cosine);

  cudaFree(f_c);
  
  return 0;
}
