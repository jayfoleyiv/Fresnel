#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include</usr/include/malloc/malloc.h>
#include</usr/include/complex.h>

// Function Prototypes
int *VEC_INT(int dim);
double *VEC_DOUBLE(int dim);
char *VEC_CHAR(int dim);
double complex *VEC_CDOUBLE(int dim);
void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C);
void TransferMatrix(double thetaI, double k0, double complex *rind, double *d,
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21);
double SpectralEfficiency(double *emissivity, int N, double *lambda, double lambdabg, double Temperature, double *P);
void Bruggenman(double f, double epsD, double complex epsM, double *eta, double *kappa);
void MaxwellGarnett(double f, double epsD, double complex epsM, double *eta, double *kappa);
//void Lorentz(double we, double de, double w, double *epsr, double *epsi);
void Lorentz(int num, double *params, double w, double *epsreal, double *epsimag);
int ReadDielectric(char *file, double *lambda, double complex *epsM);
int IsDominated(int idx, int LENGTH, double *O1,double *O2);
void ReadBRRind(int numBR, double lambda, double *BRlambda, double complex *BRind, double *n, double *k); 
// Global variables
int Nlayer;
int polflag;
double c = 299792458;
double pi=3.141592653589793;
double hbarev = 6.582119e-16;
int main(int argc, char* argv[]) {
  // integer variables
  int i, j, k, F1, F2, TK, NI, NV;
  // complex double precision variables
  double complex m11, m21, r, t, st, cosL;
  double complex *rind, eps_metal;
  double eps_real, eps_imag, nlow, nhi;
  double PU;

  // Drude + 2L parameters
  double *d2l;
  double h=6.626e-34;
  double kb = 1.38064852e-23;
  double rho;



  // real double precision variables
  double R, T, A, *d, thetaI, lambda, k0, alpha, beta;
  double sti, n1, n2, thetaT, rp, Rp, Tangle;
  double eta, kappa;
  double we, de, w;
  double dalloy, d1, d2, d3, d4, vf1, vf2, epsbg, fac1, fac2;

  // Lists for spectral efficiency
  double *LamList, *Emiss, *clam;
  // Variables for Spectral Efficiency
  double Temp, lbg;
  // This is intentionally larger than number of wavelength values in data files
  int NumLam=4000;
  // Number of variations to try
  int numVars;


  //  Allocate arrays for spectral efficiency
  LamList = VEC_DOUBLE(NumLam);
  Emiss   = VEC_DOUBLE(NumLam);
  clam    = VEC_DOUBLE(NumLam);

  FILE *fp, *pfp;

  // Character string(s)
  char *write, *line;

  write   = VEC_CHAR(1000);
  line    = VEC_CHAR(1000);

  //  Did we pass a filename to the program?
  if (argc==1) {
    exit(0);
  }

  strcpy(write,argv[1]);

  // initialize variables to be read
  nlow=0.;
  nhi=0.;
  epsbg=3.097;
  d1 = 0.180;
  d2 = 0.120;
  d3=0.01;
  d4=0.9;
  lbg=0.0;
  Temp=0.0;
  // Open the file for writing!
  fp = fopen(write,"r");
 

  fscanf(fp,"%s",line);  // Number_of_variations
  fscanf(fp,"%i",&numVars);  
  fscanf(fp,"%s",line);   // dalloy
  fscanf(fp,"%lf",&dalloy);
  fscanf(fp,"%s",line);  //  nlow
  fscanf(fp,"%lf",&nlow);
  fscanf(fp,"%s",line);  // nhi
  fscanf(fp,"%lf",&nhi);
  fscanf(fp,"%s",line);  // epsbg
  fscanf(fp,"%lf",&epsbg);
  fscanf(fp,"%s",line);  //Lambda bg
  fscanf(fp,"%lf",&lbg);
    
  polflag=1;

  if (numVars>7) {

    printf("  currently numVars is limited to 7, reducing to 7!\n");
    numVars=7;

  }

  int numVf, numNlayers, numFac, *NLa, *PF, numT;
  double *VFa, *SEA, *SFAC, *SDA, *Tem;
  
  numNlayers=numVars;
  numFac=numVars;
  numVf = numVars;
  numT  = numVars;

  // d2l parameters should have (9*number_of_variations) elements in general
  d2l = (double *)malloc((numVars*9)*sizeof(double));
  NLa = (int*)malloc((numNlayers*sizeof(int)));
  SEA = (double*)malloc((numT*numFac*numVf*numNlayers*numFac*sizeof(double)));
  SDA = (double*)malloc((numT*numFac*numVf*numNlayers*numFac*sizeof(double)));
  PF  = (int*)malloc((numT*numFac*numVf*numNlayers*numFac*sizeof(int)));
  SFAC= (double*)malloc((numFac*sizeof(double)));
  VFa = (double*)malloc((numVf*sizeof(double)));
  Tem = (double*)malloc(numT*sizeof(double));

  pfp = fopen("PARAMS.txt","r");
  double pval;
  for (i=0; i<numVars; i++) {

    for (j=0; j<9; j++) {

      fscanf(pfp,"%lf",&pval);
      d2l[i*9+j] = pval;
    }
  }
        



  d = VEC_DOUBLE(1000);
  rind = VEC_CDOUBLE(1000);
  char *type;
  type = (char *)malloc(1000*sizeof(char));

  printf("  Type  epsinf       ampD          gammaD        ampL1         omL1          gammaL1       ampL2         omL2          gammaL2       NL vf            d1            d2            Temp         SE                SD\n");

  // ML inrements through the various materials
  for (int ML=0; ML<numVars; ML++) {
 
    if (ML==0) strcpy(type,"Cr");
    else if (ML==1) strcpy(type,"Pd");
    else if (ML==2) strcpy(type,"Pt");
    else if (ML==3) strcpy(type,"Rh");
    else if (ML==4) strcpy(type,"Ta");
    else if (ML==5) strcpy(type,"V");
    else if (ML==6) strcpy(type,"W");
    else {

      printf("  not a valid type!\n");
      exit(0);
    }
    for (int TI=0; TI<numT; TI++) {

      Temp = 1200+(1000/17.)*TI;
      //Temp = 1207.123+(1320/17.)*TI;
      Tem[TI] = Temp;

      for (TK=0; TK<numVf; TK++) {

        //0.10122+(1./20)*TK;
        vf1 = 0. + (1./20)*TK;
        VFa[TK] = vf1;
  
        // Loop over different factors to multiply d1 by
        for (F1=0; F1<numFac; F1++) {

          fac1 = 0.6 + (1./20)*F1;   
          //fac1 = 0.55322+(1./20)*F1;
          SFAC[F1] = fac1;

          // Loop over different factors to multiply d2 by
          for (F2=0; F2<numFac; F2++) {

            fac2 = 0.6 + (1./20)*F2;
            //fac2 = 0.71123+(1./20)*F2;
            //  Loop over different number of layers
       
            for (NI=0; NI<numNlayers; NI++) { 

              Nlayer = 7 + NI;
              //Nlayer = 9;

              // NLa vector stores the value of Nlayer_i for layter use
              NLa[NI] = Nlayer;

              // Vector to store the thicknesses in micrometers of each layer
              d[0] = 0.;

              d[1] = dalloy;

              // Refractive index of air
              rind[0] = 1.00 + 0.*I;

              rind[1] = 1.00 + 0.*I;

              // Now start the Bragg Reflector
              for (i=2; i<Nlayer-2; i++) {

                if (i%2==0) {
                  d[i] = d1*fac1;
                  rind[i] = nlow + 0.*I;
                }
                else {
                  d[i] = d2*fac2;
                  rind[i] = nhi + 0.*I;
                }
              }

              d[Nlayer-3] = 0.01;
              rind[Nlayer-3] = sqrt(epsbg) + 0.*I;
              // W layer that is the substrate for the Bragg Reflector
              d[Nlayer-2] = 0.9;
              // Temporary - will replace with Tungsten!
              rind[Nlayer-2] = 1.0 + 0.*I;
 
              // Air underneath
              d[Nlayer-1] = 0.;
              rind[Nlayer-1] = 1.0 + 0.*I;
 
 
             //  Top/Bottom layer RI for Transmission calculation
             n1 = creal(rind[0]);
             n2 = creal(rind[Nlayer-1]);
   
             // Normal incidence
             thetaI = 0;
             for (i=0; i<NumLam; i++) {
    
               lambda = 400e-9 + i*2e-9; 
               LamList[i] = lambda;    // Lambda in meters
               k0 = 2*pi*1e-6/lambda;  // k0 in inverse microns - verified
               w=2*pi*c/lambda;        // angular frequency 
  
           
               //w_ald[i]*w_ald[i];
               // Alloy superstrate Layer (Layer 1 in the structure [Layer 0 is air!])
               Lorentz(ML*numVars, d2l, w*hbarev, &eps_real, &eps_imag);
               eps_metal = eps_real + I*eps_imag;
  
               MaxwellGarnett(vf1, epsbg, eps_metal, &eta, &kappa);
  
               //Bruggenman(vf1,epsbg, epsald, &eta, &kappa);
               rind[1] = eta + I*kappa; 
  
  
               // Alumina layer
               rind[Nlayer-3] = sqrt(epsbg);
  
               // W substrate layer (Layer N-2 in the structure [layer N-1 is air!])
               MaxwellGarnett(1.0, epsbg, eps_metal, &eta, &kappa);
               rind[Nlayer-2] = eta + I*kappa;
   
               // Solve the Transfer Matrix Equations
               TransferMatrix(thetaI, k0, rind, d, &cosL, &beta, &alpha, &m11, &m21);
 
               rho = (2*pi*h*c*c/pow(lambda,5))*(1/(exp(h*c/(lambda*kb*Temp))-1));
 
               // Fresnel reflection coefficient (which is complex if there are absorbing layers)
               r = m21/m11; 
 
               // Stored energy 
               st = (r + 1);
               st = st*conj(st);
               // Fresnel transmission coefficient (also complex if there are absorbing layers)
               t = 1./m11;
 
               // Reflectance, which is a real quantity between 0 and 1
               R = creal(r*conj(r));
               Tangle =  n2*creal(cosL)/(n1*cos(thetaI));
               T = creal(t*conj(t))*Tangle;
               A = 1 - R - T;
 
               // Store absorbance/emissivity in array Emiss
               Emiss[i] = A;
 
 
             }
         
             double SE = SpectralEfficiency(Emiss, NumLam, LamList, lbg, Temp, &PU);
             //printf("HI %8.6f    %8.6f   %8.6f    %8.6e    %8.6e    %8.6f    %i     %12.10f  %12.10e\n",
             //           vf1,     nhi,    nlow,    d1*fac1, d2*fac2, Temp,  Nlayer,  SE,      PU);

printf("  %s   %8.6e  %8.6e  %8.6e  %8.6e  %8.6e  %8.6e  %8.6e  %8.6e  %8.6e  %i  %8.6e  %8.6e  %8.6e  %8.6e  %8.6e  %8.6e\n",
         type,d2l[ML*9+0],d2l[ML*9+1],d2l[ML*9+2],d2l[ML*9+3],d2l[ML*9+4],d2l[ML*9+5],d2l[ML*9+6],d2l[ML*9+7],d2l[ML*9+8],
         Nlayer, vf1, d1*fac1, d2*fac2, Temp, SE, PU);
//printf("  Type  epsinf  ampD  gammaD  ampL1  omL1  gammaL1   ampL2   omL2   gammaL2   NL   vf    d1        d2        Temp         SE                SD\n");
//             printf("  %8.6f  %i  %8.6f  %8.6f  %8.6f  %12.10e  %12.10e\n",vf1,Nlayer, d1*fac1, d2*fac2, Temp, SE, PU);
//             SEA[TI*numT*numVf*numFac*numFac+TK*numVf*numFac*numFac+F1*numFac*numFac+F2*numFac+NI] = SE;
//             SDA[TI*numT*numVf*numFac*numFac+TK*numVf*numFac*numFac+F1*numFac*numFac+F2*numFac+NI] = PU;

           }
         }
       }
     }
   }
 }
 return 0;

}
 
 /*
  FILE *pf;
  pf = fopen("Pareto_18_BruggenmanAlloy_Vendor_BR_Temp.txt","w");
  int id;

  for (int TI=0; TI<numT; TI++) {

  for (TK=0; TK<numVf; TK++) {

  for (F1=0; F1<numFac; F1++) {

  for (F2=0; F2<numFac; F2++) {
  
  for (NI=0; NI<numNlayers; NI++) {

    i = TI*numT*numVf*numFac*numFac+TK*numVf*numFac*numFac+F1*numFac*numFac+F2*numFac+NI;
    // id is 1 if member i is dominated by at least one other member j!=i
    // a member is pareto optimal only if it is NOT dominated 
    id = IsDominated(i, numT*numVf*numFac*numNlayers*numFac, SEA, SDA);
    if (id) PF[i] = 0;
    else {
      PF[i] = 1;
      fprintf(pf,"  %f  %f     %f       %f   %i     %12.10f  %12.10e\n",
                  VFa[TK],SFAC[F1]*d1,SFAC[F2]*d2,Tem[TI],NLa[NI],SEA[i],SDA[i]);
    }


  }
  }
  }
  }
  }
  

return 0;

}
*/

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

double SpectralEfficiency(double *emissivity, int N, double *lambda, double lbg, double T, double *P){
    int i;
    double dlambda, sumD, sumN;
    double l, em;
    double h = 6.626e-34;
    double kb = 1.38064852e-23;
    double rho;

    sumD = 0;
    sumN = 0;

    for (i=1; i<N-1; i++) {

      em = emissivity[i];
      l = lambda[i];
      rho = (2.*pi*h*c*c/pow(l,5))*(1/(exp(h*c/(l*kb*T))-1));

      dlambda = fabs((lambda[i+1]- lambda[i-1])/(2));
      sumD += em*rho*dlambda;

      if (l<=lbg) {

        sumN += (l/lbg)*em*rho*dlambda;

      }
    }

    *P = sumN;
 
    return sumN/sumD;

}


void Bruggenman(double f, double epsD, double complex epsM, double *eta, double *kappa) {
  // medium 1 is surrounding medium (dielectric)
  // medium 2 is inclusion (W) - f passed to function is volume fraction of inclusion
  double f1, f2;
  double complex b, eps1, eps2, epsBG;
  eps1 = epsD + 0.*I;
  eps2 = epsM;


  f1 = (1 - f);
  f2 = f;
  b = (2*f1 - f2)*eps1 + (2*f2 - f1)*eps2;

  epsBG = (b + csqrt(8.*eps1*eps2 + b*b))/4.;

  // test to see that epsBG satisfy Bruggenman condition
  double complex test;
   *eta   = creal(csqrt(epsBG));
   *kappa = cimag(csqrt(epsBG));

}


void MaxwellGarnett(double f, double epsD, double complex epsM, double *eta, double *kappa) {
   double complex num, denom;

   num   = epsD*(2*f*(epsM - epsD) + epsM + 2*epsD);
   denom = 2*epsD + epsM + f*(epsD-epsM); 

   *eta   = creal(csqrt(num/denom));
   *kappa = cimag(csqrt(num/denom));

}

//  Evaluates real and imaginary part of refractive index from 
//  the Lorent oscillator model given omega_0, gamma_0, and omega
void Lorentz(int num, double *params, double w, double *epsreal, double *epsimag) {

  double epsinf, ampD, gammaD, ampL1, ampL2, omL1, omL2, gammaL1, gammaL2;
  epsinf  = params[num+0];
  ampD    = params[num+1];
  gammaD  = params[num+2];
  ampL1   = params[num+3];
  omL1    = params[num+4];
  gammaL1 = params[num+5];
  ampL2   = params[num+6];
  omL2    = params[num+7];
  gammaL2 = params[num+8];  

  double complex epsilon;

  // Drude
  epsilon =  epsinf -  ampD/(w*w + I*gammaD*w);
  // Lorentz1
  epsilon += ampL1/((omL1*omL1-w*w) - I*gammaL1*w);
  // Lorentz2
  epsilon += ampL2/((omL2*omL2-w*w) - I*gammaL2*w);

  //printf("  w:  %12.10e  we:  %12.10f  de:  %12.10f  epsr:  %12.10f  epsi:  %12.10f\n",w,we,de,creal(epsilon),cimag(epsilon));

  *epsreal = creal(epsilon);
  *epsimag = cimag(epsilon);

}

int ReadDielectric(char *file, double *lambda, double complex *epsM) {
   int i;
   FILE *fp;
   double lam, epsr, epsi;

   fp = fopen(file,"r");

   i=0;
   while(!feof(fp)) {

     fscanf(fp, "%lf",&lam);
     fscanf(fp, "%lf",&epsr);
     fscanf(fp, "%lf",&epsi);

     lambda[i] = lam;
     epsM[i]   = epsr + I*epsi;

     i++;
   }

   printf("#  There are %i elements in file %s\n",i,file);
   fflush(stdout);
   return i;
   fclose(fp);
}

void ReadBRRind(int numBR, double lambda, double *BRlambda, double complex *BRind, double *n, double *k) {
  int i, fdx, bdx, die;
  double temp, eta, kappa;

  // The wavelength we are interested in is smaller than any in the range of data
  if (lambda<BRlambda[0]) {

    *n = creal(BRind[0]) + (lambda - BRlambda[0])*((creal(BRind[1]) - creal(BRind[0]))/(BRlambda[1] - BRlambda[0]));
    *k = cimag(BRind[0]) + (lambda - BRlambda[0])*((cimag(BRind[1]) - cimag(BRind[0]))/(BRlambda[1] - BRlambda[0]));


  }
  // The wavelength we are interested in is larger than any in the range of data
  else if (lambda>BRlambda[numBR-2]) {

    *n = creal(BRind[numBR-2]) +(lambda - BRlambda[numBR-2])*((creal(BRind[numBR-2]) - creal(BRind[numBR-3]))/(BRlambda[numBR-2] - BRlambda[numBR-3]));
    *k = cimag(BRind[numBR-2]) +(lambda - BRlambda[numBR-2])*((cimag(BRind[numBR-2]) - cimag(BRind[numBR-3]))/(BRlambda[numBR-2] - BRlambda[numBR-3]));


  }
  // We need to scan the data to find the BRlambda for two lambdas that straddle the lambda of interest
  else {

    i=0; 
    die=1;
    do {

      temp = BRlambda[i];
      if (temp>lambda) {
      
        die=0;
        fdx = i;
        bdx = i-1; 

      }
      else i++; 

    }while(die);

    *n = creal(BRind[bdx]) + (lambda - BRlambda[fdx])*((creal(BRind[fdx]) - creal(BRind[bdx]))/(BRlambda[fdx] - BRlambda[bdx]));
    *k = cimag(BRind[bdx]) + (lambda - BRlambda[fdx])*((cimag(BRind[fdx]) - cimag(BRind[bdx]))/(BRlambda[fdx] - BRlambda[bdx]));
  
  }

}


int IsDominated(int idx, int LENGTH, double *O1,double *O2) {
  int i, is, rval;
  double Val1, Val2;

  Val1 = O1[idx];
  Val2 = O2[idx];

  // start by assuming solution is NOT dominated
  rval = 0;
  for (i=0; i<LENGTH; i++)

      if (i!=idx) {

        // Trying to maximize the function, xi dominates xidx if 
        // fj(xi) >= fj(xidx) for all j and fj(xi) < fj(xidx) for at least one j
        if ((O1[i]>=Val1 && O2[i]>=Val2) && (O1[i]>Val1 || O2[i]>Val2)) {

          //printf("  x%i is dominated by x%i\n",idx,i);
          //printf("  f1(%i):  %12.10f  f1(%i): %12.10f  f2(%i): %12.10f  f2(%i):  %12.10f\n",
          //idx,Val1,i,O1[i],idx,Val2,i,O2[i]);
          //printf("  terminating early!  i is %i out of %i\n",i,LENGTH);
          i=LENGTH;
          rval = 1;

        }
      }
  return rval;
}

