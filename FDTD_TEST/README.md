- These files produced with the commercial software Lumerical.

Ref_Trans_Air_n_1.5_0.1i_Air_d_400nm_lambda_500nm.txt
Ref_Trans_Air_n_3.5_n_2_d_400nm_lambda_500nm.txt
Ref_Trans_Air_n_1.5_Air_d_400nm_lambda_500nm.txt
Ref_Trans_n_2_n_3.5_Air_d_400nm_lambda_500nm.txt
Ref_Trans_Air_n_1.5_n_2_d_400nm_lambda_500nm.txt

- Lumerical is presumably well tested for calculating reflection and transmission through multi-layer planar stacks

- Each calculation was done with 3 layers:
-- the first being semi-infinite (d[0]=0.)
-- the second being 400 nm (d[1]=0.4 microns),
-- the third being semi-infinite (d[0]=0.)

- Column 1 is the angle of incidence in degrees
- Column 2 is the reflection
- Column 3 is the transmission

- 500 nm light is used for all simulations

Several Variations of the refractive index are given for layer 0, 1, and 2 as follows:

Ref_Trans_Air_n_1.5_0.1i_Air_d_400nm_lambda_500nm.txt ->  n[0] = 1.0, n[1] = 1.5 + 0.1*I, n[2] = 1.0
Ref_Trans_Air_n_3.5_n_2_d_400nm_lambda_500nm.txt      ->  n[0] = 1.0, n[1] = 3.5, n[2] = 2.0
Ref_Trans_Air_n_1.5_Air_d_400nm_lambda_500nm.txt      ->  n[0] = 1.0, n[1] = 1.5, n[2] = 1.0
Ref_Trans_n_2_n_3.5_Air_d_400nm_lambda_500nm.txt      ->  n[0] = 2.0, n[1] = 3.5, n[2] = 1.0
Ref_Trans_Air_n_1.5_n_2_d_400nm_lambda_500nm.txt      ->  n[0] = 1.0, n[1] = 1.5, n[2] = 2.0
