#include<stdio.h>
#include<math.h>
#include <vector>
#include <string>

using namespace std;
typedef vector<vector<double> > MyMatrixDataType_2D;

#include "constantes.h" 

void write_densidad(const double *,FILE *);


int main(){

  FILE *disco0,*disco1,*disco2;

  int numprocs = 8;
  int n_p_p,myproc;

  int  iterations = 1;
  int  frames = 1;
  int  D      = 3;
  int  N_x    = M_x+1;
  int  N_y    = M_y+1;
  int  N_z    = M_z+1;
  int  NN     = N_x*N_y*N_z;

  double BOX_SIZE_X = a3;
  double BOX_SIZE_Y = b3;
  double BOX_SIZE_Z = c;
  double deltax,deltay,deltaz;

  deltax = BOX_SIZE_X/N_x;
  deltay = BOX_SIZE_Y/N_y;
  deltaz = BOX_SIZE_Z/N_z;

  int i,j,k,l,index;
  float x_aux;
  char root[50]="test_";
  char input_file[50];
  
  double x[NN];	

  disco0 = fopen("gs3D_test_paralelo.dat","w");
  disco2 = fopen("gs1D_test_paralelo.dat","w");


  n_p_p = int(NN/numprocs);
        
   for(myproc=0;myproc<numprocs;myproc++){
     sprintf(input_file,"%s%i",root,myproc);
     disco1 = fopen(input_file,"r");      
     for(i=0;i<n_p_p;i++){
       fscanf(disco1,"%E",&x_aux);
       index = i +  myproc*n_p_p;
       x[index] = x_aux;	
//       printf("%6i,%15.6E \n",index,x_aux);
     }
   }
      
   for(i=0;i<NN;i++){
     fprintf(disco0,"%15.6E\n",x[i]);
    } 

    write_densidad(x,disco2);
    

  return 0;
}


void write_densidad(const double *phi,FILE *disco5){

  int i,j,k,l;
  double x_j,y_k,z_l;

  double atoms, range, densidad[(M_x+1)];

   MyMatrixDataType_2D::size_type x_dim = (M_x+1);
   MyMatrixDataType_2D::size_type y_dim = (M_y+1);

   MyMatrixDataType_2D densidad2D(x_dim,vector<double>(y_dim,0.0));


  for(j=0;j<=M_x;j++){
    densidad[j] = 0.0;
  }

  for(k=0;k<=M_y;k++){
      for(j=0;j<=M_x;j++){
	densidad2D[j][k] =0.0;;
      }
   }
  
  for(l=0;l<=M_z;l++){
    for(k=0;k<=M_y_2;k++){
      for(j=0;j<=M_x;j++){
	i     = j + k*(M_x+1) + l*(M_x+1)*(M_y+1);    
	densidad[j] += phi[i]*phi[i];
      }
    }
  }    
  
  for(j=0;j<=M_x;j++){
    x_j = -0.5*a3 + deltax_M*j;
    densidad[j] = densidad[j];
    fprintf(disco5,"%15.6E %15.6E\n",1E6*x_j,N_0*densidad[j]*deltay_M*deltaz_M);
  }

  for(l=0;l<=M_z;l++){
    for(k=0;k<=M_y;k++){
      for(j=0;j<=M_x;j++){
	i     = j + k*(M_x+1) + l*(M_x+1)*(M_y+1);    
	densidad2D[j][k] += phi[i]*phi[i];
      }
    }
  }    
  
  fprintf(disco5,"\n\n");

  for(k=0;k<=M_y;k++){
  y_k = -0.5*b3 + deltay_M*k;
  for(j=0;j<=M_x;j++){
    x_j = -0.5*a3 + deltax_M*j;
    fprintf(disco5,"%15.6E %15.6E %15.6E\n",1E6*x_j,1E6*y_k,N_0*densidad2D[j][k]*deltaz_M);
  }
  fprintf(disco5,"\n");
  }
  
}
