- Contains c-code that varies over geometric and material parameters and computes figures of merits
- Can be used to create both test and training data sets

- To compile, type

`g++ -o Train.exe Train.c -O3`

- To run and save output to a textfile called "Test.txt", type:

`./Train.exe input.txt > Test.txt &`

- input.txt need not be modified unless you are interested in changing lambda_bg for the PV

c Drude :

       cc=epsinf - ampD/(om2 + ri   gammaD om)
       
c Lorentz 1 :

        cc = cc + ampL/( (omL^2 - om2) - ri gammaL om )
        
c Lorentz 2 :

        cc = cc + ampL/( (omL^2 - om2) - ri gammaL om )
        
TO DO:
- Code in way to read D2L parameters from file
- Validate the D2L function is working (reproduces the fit)
- Post results to github
