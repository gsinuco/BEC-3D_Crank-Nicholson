#ifndef CONSTANTES_H
#define CONSTANTES_H

double length_scale = 1.0; 

//--------- PHYSICAL CONSTANTS ------------------------------------

double e        = 1.602176462E-19, epsilon_cero = 8.854187817E-12, mu_cero = 12.566370614E-7;
double mu       = 9.2740154E-24,   epsilon      = 12.5;
double hbar     = 1.054E-34;
double mass_at  = 1.325E-25,       a_s = 5.4E-9; //^{87}Rb
double mu_B     = 9.2740154E-24;
double k_B      = 1.38E-23;

// ------------ PARAMETERS OF 2DEG -----------------------------------

double sigma    = 3.224287E15;
double mobility = 0.5;
double mass_e   = 0.067*9.1094E-31;
double h        = 529.183E-10;
double cte      = mass_e/(M_PI*hbar*hbar);

double channels;
double current_2deg = 0.0E-4;
double applied_voltage = 5.0E-4;
double Number_wires     = 11;
double wires_separation = 1E-6;


// ------------ GRID PARAMETERS FOR RANDOM FLUCTUATIONS --------------806071750

double a  = 2.0E-6, b  = 2.0E-6, thick =0.0;
double a2 = 2.0E-6, b2 = 2.0E-6, at =  2.0E-6;
int    N_d = int(sigma*a2*b2), N_x = 100, N_y = 100,   N_z = 0;
int    N = (N_x+1)*(N_y+1);
//double deltax=a2/N_x,deltay=b2/N_y;


//------------- GRID PARAMETERS FOR MAGNETIC FIELD

double a3  = 200.0E-6, b3  = 20.0E-6, c   = 20.0E-6, z_0 = 0.0E-6; // para la trampa
int    M_x = 399,   M_y = 39,  M_z =  39;

int    NN = (M_x+1)*(M_y+1)*(M_z+1);
double deltax_M=a3/(M_x+1),deltay_M=b3/(M_y+1), deltaz_M=c/(M_z+1);

// ------------ PARAMETERS FOR TRAPING ------------

double delta_x=50E-6,delta_y=0.5E-6,delta_z=0.5E-6; // needed for frequency definition

double N_0     =    1.0E5;
double B_bias  =    0.0E-4;
double B_off   =    30.0E-4;

double w_x = 2.0*M_PI*15.0;
double w_y = 2.0*M_PI*200.0;
double w_z = 2.0*M_PI*200.0;

int numero_de_cajas   = 1;
int M_y_2             = int((M_y)/numero_de_cajas);
int NN_2              = int((M_x+1)*(M_y_2+1)*(M_z+1));

int  num_procs    = 4; 
int  N_J          = int((N_x+1)*(N_y+1)*(N_z+1));     
int  N_J_p        = int(N_J/num_procs);
int  N_p          = int((M_x+1)*(M_y+1)*(M_z+1));
int  N_p_p        = int(N_p/num_procs);

//------------- 2DEG 


#endif 


