- Contains c-code that varies over geometric and material parameters and computes figures of merits
- Can be used to create both test and training data sets

- To compile, type

`g++ -o Train.exe Train.c -O3`

- To run and save output to a textfile called "Test.txt", type:

`./Train.exe input.txt > Test.txt &`

- Features that can be modified in input.txt:

-- number_of_variations:  Basically how many different values of a particular feature will be tried.  
   Note that currently the code will go sequentially through a list of Drude-Lorentz parameters fit to 
   the permittivity of different materials, and will try number_of_variations different materials.  
   For example, if number_of_variations is 2, then Cr and Pd will be tried.  

-- n1 and n2 - These will define the refractive indices of the materials that make up the Bragg Reflector

- Additional features *can* be modified, but one should think carefully before doing so!

- Currently the file "Test.txt" contains the output from running the program with number_of_variations = 3, so Cr, Pd, and Pt are considered at
  along with all possible combinations of 3 different values of d1, d2, Temp, volume fraction, and number of layers in the Bragg reflector 

- About the Drude + 2 Lorentz Permittivity model:

c Drude :

       cc=epsinf - ampD/(om2 + ri   gammaD om)
       
c Lorentz 1 :

        cc = cc + ampL/( (omL^2 - om2) - ri gammaL om )
        
c Lorentz 2 :

        cc = cc + ampL/( (omL^2 - om2) - ri gammaL om )
        

PARAMS:
Cr
0.221474E+01  0.261386E+02  0.957129E-01
0.392081D+01  0.585496D+00  0.314290D+00
0.173123D+03  0.186294D+01  0.282379D+01
Pd
0.276486E+01  0.333670E+02  0.194411E-01
0.400000D+02  0.189432D+01  0.208485D+01
0.409110D+02  0.537494D+00  0.944586D+00
Pt
0.667131E+01  0.308192E+02  0.950157E-01
0.400000D+02  0.201085D+01  0.157733D+01
0.440199D+02  0.890291D+00  0.831543D+00
Rh
0.807728E+00  0.632151E+02  0.620410E-01
0.816751D+01  0.120125D+01  0.529252D+00
0.251214D+03  0.548138D+00  0.490489D+01
Ta
0.669500E+01  0.568556E+02  0.658926E-01
0.400000D+02  0.135609D+01  0.596268D+01
0.638460D+02  0.323698D+01  0.184396D+01
V
0.265683E+01  0.274292E+02  0.657368E-01
0.930849D+00  0.590744D+00  0.326088D+00
0.746247D+02  0.220001D+01  0.258214D+01
W
0.657569E+01  0.356285E+02  0.523842E-01
0.146288D+01  0.962884D+00  0.185501D+00
0.400000D+03  0.336333D+01  0.877130D+01
