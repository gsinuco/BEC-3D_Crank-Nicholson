#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <vector>
#include <string>

using namespace std;
typedef complex<double> dcmplx;
dcmplx I(0,1);

typedef vector<vector<double> > MyMatrixDataType_2D;

#include "constantes.h" 

void write_densidad(const dcmplx *,FILE *);


int main(){

  FILE *disco0,*disco1,*disco2;

  int numprocs = 8;
  int n_p_p,myproc;

  int  iterations = 10000;
  int  frames = 10;
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

  int i,config;
  int j,k,l,index;
  int channels,r;
  float x_aux,x_aux_;

  char root[50]="test_te_grid_";
  char input_file[50];
  char comment[100];

  dcmplx *x_;
  dcmplx **x;
  dcmplx x_r,x_i;

  x_= new dcmplx [NN];
  x = new dcmplx* [frames + 2];

  for(j=0;j<frames+2;j++){
      x[j] = new dcmplx [NN];
  }

  /*dcmplx x[frames+2][NN],x_[NN],x_r,x_i;*/
  double x_j,density_1D[N_x];
  
  disco0 = fopen("te3D_test.dat","w+");
  disco2 = fopen("te1D_test.dat","w+");

  
  n_p_p = int(NN/numprocs);

  for(j=0;j<N_x;j++){
    density_1D[j] = 0.0;
  }
  
  
    for(myproc=0;myproc<numprocs;myproc++){
      sprintf(input_file,"%s%i",root,myproc);
      disco1 = fopen(input_file,"r");      
      for(channels=0;channels<1;channels++){
      for(r=0;r<frames+2;r++){
      config = r + channels*2;
      for(i=0;i<n_p_p;i++){
	fscanf(disco1,"%E%E",&x_aux,&x_aux_);
	index = i +  myproc*n_p_p;
        x_r = x_aux;
	x_i = x_aux_;
	x[config][index] = x_r + x_i*I;	
        //fprintf(disco0,"%i %i %i %i %i %15.6E\n",myproc,channels,r,config,index,x[config][index]);
      }
      }
      }
    }
    

         
    for(config=0;config<frames+2;config++){
       for(i=0;i<NN;i++){
	  x_[i] = x[config][i];
	  fprintf(disco0,"%15.6E %15.6E\n",real(x_[i]),imag(x_[i]));
       } 
       write_densidad(x_,disco2);
    }
    
  

  return 0;
}

/*
void write_densidad(const dcmplx *phi,FILE *disco5){

  int i,j,k,l;
  double x_j,y_k,z_l;

  double atoms, range, densidad[(M_x+1)];

  for(j=0;j<=M_x;j++){
    i     = j + k*(M_x+1) + l*(M_x+1)*(M_y+1);    
    densidad[j] = 0.0;
  }
  
  for(l=0;l<=M_z;l++){
    for(k=0;k<=M_y_2;k++){
      for(j=0;j<=M_x;j++){
	i     = j + k*(M_x+1) + l*(M_x+1)*(M_y+1);    
	densidad[j] += abs(phi[i])*abs(phi[i]);
      }
    }
  }    
  
  for(j=0;j<=M_x;j++){
    x_j = -0.5*a3 + deltax_M*j;
    densidad[j] = densidad[j]*deltay_M*deltaz_M;
    fprintf(disco5,"%15.6E %15.6E\n",1E6*x_j,N_0*densidad[j]/1E6);
  }

//   atoms = 0.0;
//   range = int((0.5*a3 + 1E-6)/deltax_M);
//   for(j=1990;j<=2010;j++){
//     x_j = -0.5*a3 + deltax_M*j;
//     atoms  += densidad[j]*deltax_M;
//   } 
//   fprintf(disco5,"#Atomos en 2 mu-m: %15.6E\n",N_0*atoms);
//   
  
  fprintf(disco5,"\n\n"); 
  
}
*/

void write_densidad(const dcmplx *phi,FILE *disco5){

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
	densidad[j] += abs(phi[i])*abs(phi[i]);
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
	densidad2D[j][k] += abs(phi[i])*abs(phi[i]);
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
  fprintf(disco5,"\n");
}
