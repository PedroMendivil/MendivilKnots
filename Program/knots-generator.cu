/////////////////////////////////////////////////////////////////////
// PROGRAMA GENERADOR DE NUDOS DE PEDRO MENDIVIL Y CHAT GPT        //
// CALCULA VECTORES DE POSICION INICIALES DE CADENAS DE POLIMEROS  //
// SIRVE PARA ALIMENTAR AL PROGRAMA simulate-polimer.cu            //
// VERSION V1.1    noviembre 2024                               //
/////////////////////////////////////////////////////////////////////

//Libraries     
#include <stdio.h>
#include <stdlib.h> // Necesario para system()
#include <math.h>
#include <time.h>
#include <curand_kernel.h>
#include <string.h>

// Global variables
int subir = 1;
int bajar = 0;
int n_atom = 0;
int input_in_x = 1;
int input_in_y = 0;
int input_in_z = 0;
int input_in_n_atom = 0;
int atomo_validado =0;
float x=0;
float y=0;
float z=0;
int tecla = 0;
int tecla_ant = 0;

////// AQUI SE COPIAN PARTES DE TODA UNA SECCION DE simulate-polymer.cu ///////

// Variables globales  ///////////////////////////////////////////////
int no_atomo; // Antes N lo usaremos para ir rellenando el fichero de nudo
// npartic = número de partículas. Para random usamos global_N.
int npartic; // se modifica con cada xx.gro leido.
// file_copy[100] array global, copia de file_to_copy de void read_parameters()
char file_copy[100]; // cadena que contiene el nombre del fichero de nudo
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

// Función para limpiar la pantalla
void limpiar_pantalla() {
    system("clear"); // Ejecuta el comando "clear" en sistemas basados en Unix
}

// Funciones
void cursor_a_1_1() { printf("\033[1;1H"); }     // Cursor a fila 1, columna 1
void cursor_a_16_20() { printf("\033[16;20H"); } // Cursor a fila 16, columna 20
void cursor_a_17_20() { printf("\033[17;20H"); } // Cursor a (fila, columna)
void cursor_a_16_40() { printf("\033[16;40H"); } // Cursor a (fila, columna)
void cursor_a_17_40() { printf("\033[17;40H"); } //     ""
void cursor_a_16_60() { printf("\033[16;60H"); } //     ""
void cursor_a_17_60() { printf("\033[17;60H"); } //     ""
void cursor_a_16_80() { printf("\033[16;80H"); } //     ""
void cursor_a_17_80() { printf("\033[17;80H"); } //     ""
void cursor_a_17_98() { printf("\033[20;40H"); } // Cursor a fila 17, columna 98
void cursor_a_20_40() { printf("\033[20;40H"); } // Cursor a fila 20, columna 40
void cursor_a_22_20() { printf("\033[22;20H"); } // Cursor a fila 22, columna 20

void estado_input_x() {
    if (tecla == '1') {
      input_in_x = 1;
      input_in_y = 0;
      input_in_z = 0;
      input_in_n_atom = 0;
    }
}    
void estado_input_y() {
    if (tecla == '2') {
      input_in_x = 0;
      input_in_y = 1;
      input_in_z = 0;
      input_in_n_atom = 0;
    }
}
void estado_input_z() {
    if (tecla == '3') {
      input_in_x = 0;
      input_in_y = 0;
      input_in_z = 1;
      input_in_n_atom = 0;
    }
}
void estado_input_n_atom() {
    if (tecla == '4') {
      input_in_x = 0;
      input_in_y = 0;
      input_in_z = 0;
      input_in_n_atom = 1;
    }
}
////////////////////////////  PRESENTACION DE INICIO  //////////////////////////////////

void presentacion_inicio() 
{
    printf("///////////////////////////////////////////////////////////////////////////////////////////\n");
    printf("//                                                                                       //\n");
    printf("//                                 PARA QUE QUEDE CLARO:                                 //\n");
    printf("//                                                                                       //\n");
    printf("//        ESTE PROGRAMA GENERA UN FICHERO DE POSICIONES INICIALES DE LOS MONOMEROS       //\n");
    printf("//        CON LOS QUE LUEGO ALIMENTAREMOS AL PROGRAMA simulate-polimer.bin               //\n");   
    printf("//        INSTRUCCIONES: moveremos la posición de cada partícula mediante el teclado     //\n");
    printf("//        VER INSTRUCCIONES PULSANDO LA TECLA ENTER                                      //\n");
    printf("//                                                                                       //\n");
    printf("///////////////////////////////////////////////////////////////////////////////////////////\n\n");
 
    printf("Este texto de encima será borrado al pulsar enter.\n");
   
    getchar();
    // Limpia la pantalla y muestra las posiciones iniciales
    printf("\033[2J"); // Limpiar pantalla
    cursor_a_1_1(); // IR A ORIGEN
}
void instrucciones() 
{
    printf("             /////////////////////////////////////////////////////////////////////////////////////////////////////\n");
    printf("             //                                                                                                 //\n");
    printf("             //                                 INSTRUCCIONES:                                                  //\n");
    printf("             //     PARA MODIFICAR X PULSAREMOS 1 Y DESPUES ENTER, DESPUES LA TECLA + O LA - Y DESPUES ENTER    //\n");
    printf("             //     PARA MODIFICAR Y PULSAREMOS 2 Y DESPUES ENTER, DESPUES LA TECLA + O LA - Y DESPUES ENTER    //\n");
    printf("             //     PARA MODIFICAR Z PULSAREMOS 3 Y DESPUES ENTER, DESPUES LA TECLA + O LA - Y DESPUES ENTER    //\n");
    printf("             //     PARA MODIFICAR n PULSAREMOS 4 Y DESPUES ENTER, DESPUES LA TECLA + O LA - Y DESPUES ENTER    //\n");
    printf("             //        n es el número de la partícula que vamos a almacenar en el fichero de nudo               //\n");
    printf("             //        PARA VALIDAR LA POSICION DE LA PARTICULA PULSAR 9 Y DESPUES ENTER                        //\n");
    printf("             //        PARA SALIR DEL PROGRAMA PULSAR q Y DESPUES ENTER                                         //\n");
    printf("             //  NOTA: para generar completo el archivo de nudo hay que meter antes en el directorio            //\n");
    printf("             //  de knots-generator.bin un fichero de texto (vacio) con el nombre del nudo. Este                //\n");
    printf("             //  nombre hay que cambiarlo antes de lanzar el generador de nudos en el fichero parameters.dat    //\n");
    printf("             /////////////////////////////////////////////////////////////////////////////////////////////////////\n");
}
void linea_estados_de_entrada_de_datos() {
      // Indicar la posición del dato a meter dinámicamente
      cursor_a_16_20();
      if (input_in_x == 1) { printf("↑ con +/↓ con -                                                                 "); }
      cursor_a_16_20();
      if (input_in_y == 1) { printf("                    ↑ con +/↓ con -                                             "); }
      cursor_a_16_20();
      if (input_in_z == 1) { printf("                                        ↑ con +/↓ con -                          "); }
      cursor_a_16_20();
      if (input_in_n_atom == 1) { printf("                                                           ↑ con +/↓ con -"); }
      
      cursor_a_20_40();  // Después de imprimir cursor a posición fija.
      }

void linea_con_xyz_y_n_atom() {
      // Mostrar vector posición y número de átomo
      cursor_a_17_20();
      printf(" Pos X = %.3f", x);
      cursor_a_17_40();
      printf(" Pos Y = %.3f", y);
      cursor_a_17_60();
      printf(" Pos Z = %.3f", z);
      cursor_a_17_80();
      printf(" no atom = %d  ", n_atom);
      
      cursor_a_20_40();  // Después de imprimir cursor a posición fija.
}
void analiza_tecla()  {
      // Comprobar si la tecla es '+', '-',
      if (tecla == '+') { bajar = 0; subir = 1;}
      if (tecla == '-') { bajar = 1; subir = 0;}
      // Comprobar si la tecla es '9'
      if (tecla == '9') { atomo_validado = 1;}
}
void gosub_actualizar_x() {
      if (input_in_x == 1 && subir == 1) { x = x + 1; }
      if (input_in_x == 1 && bajar == 1) { x = x - 1; }
}
void gosub_actualizar_y() {
      if (input_in_y == 1 && subir == 1) { y = y + 1; }
      if (input_in_y == 1 && bajar == 1) { y = y - 1; }
}
void gosub_actualizar_z() {
      if (input_in_z == 1 && subir == 1) { z = z + 1; }
      if (input_in_z == 1 && bajar == 1) { z = z - 1; }
}
void gosub_actualizar_n_atom() {
      if (input_in_n_atom == 1 && subir == 1) { n_atom = n_atom + 1; }
      if (input_in_n_atom == 1 && bajar == 1) { n_atom = n_atom - 1; }
}
///////////////////  ESCRITURA EN FICHERO DE LAS POSICIONES (GPT) /////////
void gosub_write_particle_to_file(const char *file_copy, int n_atom, float x, float y, float z) { 

    FILE *file = fopen(file_copy, "a+");
    if (file == NULL) {
        perror("Error al abrir el archivo");
        return;
    }
    // Escribir la partícula en el archivo con el formato especificado
    fprintf(file, "    %dX        X    %d  %7.3f  %7.3f  %7.3f\n",
            n_atom, n_atom, x, y, z);

    fclose(file);
}

//////////////////////////    MAIN     //////////////////////////////     

////////// OTRA SECCION MAS COPIADA DE MAIN DE simulate-polimer.cu /////////// 

int main( int argc, char const *argv[])
{
    // ================== DIRECTORIO =====================
    /// argv[0]: nombre del programa y argv[1]: directorio de simulación (sim_dir) 
    if( argc!=2){  //si el número de argumentos es distinto de 2
      if( argc<2){ printf("You forgot the input.\n"); exit(-1);}//si es menor de 2
      else{ printf("Too many arguments.\n"); exit(-1);} // si mayor de 2 => exit
    }
    if( sizeof(argv[1])>128){ printf("Directory name too long.\n"); exit(-1);}
    char sim_dir[128];
    snprintf(sim_dir,sizeof(sim_dir),"%s",argv[1]);

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
    /////////////// FIN COPIA de simulate polimer.cu  /////////////

    presentacion_inicio(); // Gran pantallazo paquetenteres!

    instrucciones();  // Otro gran pantallazo paquetenteres!
      
    linea_estados_de_entrada_de_datos();  // Mostrar la posición del dato a meter
    
    linea_con_xyz_y_n_atom();  // Mostrar vector posición y número de átomo
    
    while (1) {
        
        tecla = getchar();     // Leer la tecla presionada

        // Si el carácter actual es '\n', el anterior está en tecla_ant (bucle)      
        if (tecla == '\n') { tecla_ant = tecla; }
                   
        estado_input_x();      // bit de estado: genera input_in_x con tecla=1
        estado_input_y();      // bit de estado: genera input_in_y
        estado_input_z();      // bit de estado: genera input_in_z
        estado_input_n_atom(); // bit de estado: genera input_in_n_atom

        analiza_tecla();  // Controla si pulsamos + ó - y cambia subir/bajar
             
        fflush(stdout); // Asegurarse de que todo se imprima correctamente
                                                      
        if (subir == 1 || bajar == 1) {
        gosub_actualizar_x(); // Adivina: actualiza x 
        gosub_actualizar_y(); // Actualiza y ademas de 
        gosub_actualizar_z(); // Actualiza z ademas de 
        gosub_actualizar_n_atom(); // Actualiza n_atom
        //                                                          linea_estados_de_entrada_de_datos();
        //                                                          linea_con_xyz_y_n_atom();
        subir = 0;
        bajar = 0;
        cursor_a_20_40();
        }
        linea_estados_de_entrada_de_datos();
        linea_con_xyz_y_n_atom();
        
        if (atomo_validado == 1) {
        gosub_write_particle_to_file(file_copy, n_atom, x, y, z);
        cursor_a_17_98(); 
        printf("atomo_validado = %d", n_atom);
        fflush(stdout); // Asegurarse de que todo se imprima correctamente
        atomo_validado = 0;
        cursor_a_20_40();
        }

        fflush(stdout); // Asegurarse de que el cambio se vea de inmediato
        cursor_a_22_20();
        printf("Recuerde: x: tecla 1, y: tecla 2, z:tecla 3, n atomo: 4 y validar atomo con tecla 9");
        fflush(stdout); 
        cursor_a_20_40();  // Cursor: te he dicho que a (13,40).

        // Salir si se presiona 'q'
        if (tecla == 'q') {
            printf("Saliendo del programa...");
            break;
        }
        // Limpiamos el buffer de teclado.
        while (getchar() != '\n'); // Limpiar cualquier entrada sobrante
        cursor_a_20_40();  // Cursor: te he dicho que a (13,40).
    }

    return 0;
}

