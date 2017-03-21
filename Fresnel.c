#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<malloc.h>
#include<complex.h>

// Function Prototypes
int *VEC_INT(int dim);
double *VEC_DOUBLE(int dim);
char *VEC_CHAR(int dim);
double complex *VEC_CDOUBLE(int dim);
void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C);
void TransferMatrix(double thetaI, double k0, double complex *rind, double *d,
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21);

// Global variables
int Nlayer;
int polflag;

int main(int argc, char* argv[]) {
  // integer variables
  int i, j;

  // complex double precision variables
  double complex m11, m21, r, t, cosL;
  double complex *rind;

  // real double precision variables
  double refl, *d, thetaI, lambda, k0, alpha, beta, pi;

  // File pointer(s)
  FILE *fp;

  // Character string(s)
  char *write;

  write = VEC_CHAR(1000);

  //  Did we pass a filename to the program?
  if (argc==1) {
    printf("\n\n  Please pass a file name to write the reflection to as an argument! \n\n");
    printf("  That is, please type './Fresnel.exe Reflection.txt' to compute the reflection\n");
    printf("  and store it in a file called 'Reflection.txt'\n\n");
    exit(0);
  }

  // If program makes it out alive, tell us the filename so we know where to 
  // look!
  printf("\n\n  GOING TO WRITE REFLECTION TO FILE %s\n",argv[1]);

  strcpy(write,argv[1]);

  // Open the file for writing!
  fp = fopen(write,"w");

  // Go on to some more important items for actually solving the Transfer Matrix Equations!
 
  // Polarization convention (polflag=1 is p-polarized, polflag=2 is s-polarized
  polflag=1;

  // Total number of layers in the structure, including terminal layers (which are taken to have infinite thickness)
  Nlayer=4;

  pi = 3.14159265359;

  // Vector to store the thicknesses in micrometers of each layer
  d = VEC_DOUBLE(Nlayer);

  d[0] = 0.;
  d[1] = 0.299999999999999;
  d[2] = 0.265;
  d[3] = 0.;

  //  Vector to store the complex refractive index of each layer
  rind = VEC_CDOUBLE(Nlayer);
 
  //  Multilayer structure with some absorbing layers
/*  rind[0] = 3.700 + 0.009400*I;
  rind[1] = 1.453 + 0.000000*I;
  rind[2] = 4.180 + 2.390e-3*I;
  rind[3] = 3.500 + 0.000000*I;
*/

  //  4-layer structure that mimics a 2-layer structure
  rind[0] = 4.33 + 0.*I;
  rind[1] = 2.25 + 0.*I;
  rind[2] = 2.25 + 0.*I;
  rind[3] = 2.25 + 0.*I;


  // Wavelength in nanometers
  lambda = 500.;

  // Wavenumber in inverse micrometers
  k0 = 1000./lambda;

  //  Loop over incident angle and/or wavelength

  fprintf(fp,"# Incident Angle Wavelength  Reflectance \n");
  for (i=0; i<360; i++) {

     // Increment incident angle
     thetaI=(0.25*i)*3.14159/180.;

     // Solve transfer matrix equations for incident angle thatI, wavenumber k0,
     // and structure defined by layers with refractive indices stored in the vector rind
     // and geometries stored in the vector d
     TransferMatrix(thetaI, k0, rind, d, &cosL, &beta, &alpha, &m11, &m21);

     // Fresnel reflection coefficient (which is complex if there are absorbing layers)
     r = m21/m11; 
  
     // Fresnel transmission coefficient (also complex if there are absorbing layers)
     t = 1./m11;

     // Reflectance, which is a real quantity between 0 and 1
     refl = creal(r*conj(r));

     //  Print the reflectance, incident angle, and wavelength to the output file
     fprintf(fp,"  %12.10f %12.10f  %12.10f \n",thetaI*180./pi,lambda,refl);
  }

fclose(fp);
return 0;

}


// Functions
int *VEC_INT(int dim){
  int *v,i;
  v = (int *)malloc(dim*sizeof(int));
  if (v==NULL) {
     printf("\n\nVEC_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0;
  return v;
}
double *VEC_DOUBLE(int dim){
  int i;
  double *v;
  v = (double *)malloc(dim*sizeof(double));
  if (v==NULL) {
     printf("\n\nVEC_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0;
  return v;
}

double complex *VEC_CDOUBLE(int dim) {
  int i;
  double complex *v;
  v = (double complex *)malloc(dim*sizeof(double complex));
  if (v==NULL) {
     printf("\n\nVEC_CDOUBLE:  Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0 + I*0.;
  return v;
}

char *VEC_CHAR(int dim){
  char *v;
  v = (char *)malloc(dim*sizeof(char));
  if (v==NULL) {
     printf("\n\nVEC_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  return v;
}

void TransferMatrix(double thetaI, double k0, double complex *rind, double *d, 
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21) { 
  
  int i, j, k, indx;
  double complex *kz, *phiL, *D, *Dinv, *Pl, *phil;
  double complex *EM, ctheta, tmp, *tmp2, *tmp3, c0, c1, ci, kx;

  kz   = VEC_CDOUBLE(Nlayer);
  phil = VEC_CDOUBLE(Nlayer);
  D    = VEC_CDOUBLE(4*Nlayer);
  Dinv = VEC_CDOUBLE(4*Nlayer);
  Pl   = VEC_CDOUBLE(4*Nlayer);
  EM   = VEC_CDOUBLE(4);
  tmp2 = VEC_CDOUBLE(4); 
  tmp3 = VEC_CDOUBLE(4);

  c0 = 0. + I*0.;
  c1 = 1. + I*0.;
  ci = 0. + I*1.;

  //  x-component of incident wavevector...
  //  should be in dielectric material, so the imaginary 
  //  component should be 0.
  kx = k0*rind[0]*sin(thetaI);

  //  Now get the z-components of the wavevector in each layer
  for (i=0; i<Nlayer; i++) {
     kz[i] = (rind[i]*k0)*(rind[i]*k0) - kx*kx;
     kz[i] = csqrt(kz[i]);
     // Want to make sure the square root returns the positive imaginary branch
     if (cimag(kz[i])<0.)  {
        kz[i] = creal(kz[i]) - cimag(kz[i]);
     }
   }



   //  Calculate the P matrix
   for (i=1; i<Nlayer-1; i++) {
     phil[i]=kz[i]*d[i];

     //  Upper left (diagonal 1)
     Pl[i*4] = cexp(-ci*phil[i]);  
     //  upper right (off diagonal 1)
     Pl[i*4+1] = c0;
     //  lower left (off diagonal 2)
     Pl[i*4+2] = c0;
     //  lower right (diagonal 2)
     Pl[i*4+3] = cexp(ci*phil[i]);

   }

 
   //  Calculate the D and Dinv matrices
   for (i=0; i<Nlayer; i++) {
     ctheta = kz[i]/(rind[i]*k0);
     //  p-polarized incident waves
     if (polflag==1) {  

       //  Upper left (diagonal 1)
       D[i*4] = ctheta;
       // upper right
       D[i*4+1] = ctheta;
       // lower left
       D[i*4+2] = rind[i];
       // lower right
       D[i*4+3] = -rind[i];

     } 
     //  s-polarized incident waves
     if (polflag==2) {

       // upper left
       D[i*4] = 1;
       // upper right
       D[i*4+1] = 1;
       // lower left
       D[i*4+2] = rind[i]*ctheta;
       // lower right
       D[i*4+3] = -1*rind[i]*ctheta;

     }
     //  Now compute inverse
     //  Compute determinant of each D matrix
     tmp = D[i*4]*D[i*4+3]-D[i*4+1]*D[i*4+2];
     tmp = 1./tmp;

     //printf("  tmp is %12.10f  %12.10f\n",creal(tmp),cimag(tmp));    
     Dinv[i*4]=tmp*D[i*4+3];
     Dinv[i*4+1]=-1*tmp*D[i*4+1];
     Dinv[i*4+2]=-1*tmp*D[i*4+2];
     Dinv[i*4+3]=tmp*D[i*4];
 
   }


   // Initial EM matrix
   EM[0] = c1;
   EM[1] = c0;
   EM[2] = c0;
   EM[3] = c1;
   for (i=Nlayer-2; i>0; i--) {
     CMatMult2x2(i, Pl  , i,  Dinv, 0, tmp2);
     CMatMult2x2(i, D   , 0, tmp2,  0, tmp3);
     CMatMult2x2(0, tmp3, 0, EM  ,  0, tmp2); 

     for (j=0; j<2; j++) {
       for (k=0; k<2; k++) {
          EM[2*j+k] = tmp2[2*j+k];
       }
     }
   }
   CMatMult2x2(0, EM  , Nlayer-1, D   , 0, tmp2);
   CMatMult2x2(0, Dinv, 0, tmp2, 0, EM); 


   //  Finally, collect all the quantities we wish 
   //  to have available after this function is called
   *m11 = EM[0*2+0];  //  
   *m21 = EM[1*2+0];
   *beta = creal(kx);
   *alpha = cimag(kx);
   *cosL = ctheta;

   free(kz);  
   free(phil);
   free(D);
   free(Dinv);
   free(Pl);
   free(EM);
   free(tmp2);
   free(tmp3);

}

 

void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C) {
     int i, j, k, m, n;

     double complex sum;

     for (k=0; k<2; k++) {
       for (i=0; i<2; i++) {

          sum = 0. + 0*I;
          for (j=0; j<2; j++) {

            m = 2*i + j;
            n = 2*j + k;           
            sum += A[Aidx*4+m]*B[Bidx*4+n];

          }
          
          C[Cidx*4 + (2*i+k)] = sum;
        }
      }

}
