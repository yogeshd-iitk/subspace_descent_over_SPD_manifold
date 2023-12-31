# Low complexity subspace descent over symmeric positive definite manifold (Yogesh Darmwal, Ketan Rajawat)

--------------------------------------------------------------------------------------------------------------------------------------------------------------

Files

---------------------------------------------------------------------------------------------------------------------------------------------------------------

1) RSD_extension.m
 
	It is the main code file which implements Riemannian subspace descent algorithm proposed in the paper (https://arxiv.org/pdf/2305.02041.pdf)

2) user_input.m

	It contains user provided data and the symbolic expression of the objective function.

3) number_of_comp.m

	It calculates the number of terms of each type in the objective function g.

4) user_fun.m

	It is the MATLAB function corresponding to the symbolic expression in the file 'user_input.m', generated using 'matlabFunction'.

5) user_fun_grad.m

	It is the MATLAB function generated with 'matlabFunction' to compute the gradient of 'user_fun.m'.

6) user_fun_hess.m

	It is the MATLAB function generated with 'matlabFunction' to compute the Hessian of 'user_fun.m'.

7) inti_var.m

	It is the function used to calculate intermediate matrices M1p(X0), M2q(X0), M41r(X0), M42r(X0), M51s(X0), M52s(X0), M61m(X0), M62m(X0), and the function value at the initialization X0. Additionally, it computes trace(M41r(X0)*M42r(X0)), trace(M51s(X0)*M52s(X0)), and trace(M61m(X0)*transpose(M62m(X0)), serving as initializations in the unidirectional case, for the recursive computation of trace(M41r(X)*M42r(X)), trace(M51s(X)*M52s(X)), and trace(M61m(X)*transpose(M62m(X))) at any point X with O(n) complexity. This is in contrast to the direct calculation's O(n^2) cost.

8) direct_deriv.m

	It computes the directional derivatives of the functions h1p(X), h2q(X), h3(X), h4r(X), h5s(X), and h6m(X) at the point X.

9) subspace_select_second_order_greedy.m

	It generates a subspace using a greedy heuristic method.

10) subspace_select_second_order_random.m

	It generates a subspace using a randomized method.

11) mod_entry.m

	In the randomized unidirectional case, this function is utilized to store values that undergo modifications in the current iteration. These stored values are then used to compute trace(M41r(X)*M42r(X)), trace(M51s(X)*M52s(X)), and trace(M61m(X)*transpose(M62m(X))) at any point X with O(n) complexity.

12) Bup_multiplcatn.m

	It post-multiplies any matrix A by the update matrix 'Bup,' denoted as A*Bup.

13) Bupintrns_multiplcatn.m

	It pre-multiplies any matrix A by the inverse of the update matrix 'Bup,' denoted as inverse(Bup)*A.

14) intrmdt_var.m

	It computes intermediate variables utilized in the functions 'trace_hessV_V_fun.m' and 'trace_V_grad_V_Xin_fun.m'.

15) trace_grad_V_fun.m

	It computes 'trace(gradf *V)'.

16) trace_hessV_V_fun.m

	It computes 'trace(hessf[V]*V)'.

17) trace_V_grad_V_Xin_fun.m

	It computes 'trace(V*gradf*V*inverse(X))'.

18) fun_min_con.m 

	It implements the alternate method of solving positive definite matrix-constrained problems using the MATLAB function 'fmincon'.

19) user_fun_con.m

	It performs function evaluation with respect to the variable X using the 'user_fun.m' function in the 'fun_min_con.m' code.

20) lambda_star_fun.m

	It calculates the optimal values of lambdas for each basis vector in the chosen subspace.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

Data files generated by the code

---------------------------------

1)  user_data.mat
				
   	It stores user-provided data, including constant matrices Cp, Dq, Ar, Hr, Fs, Gs, Pm, Qm, as well as initializations such as X0, B0, inverse(B0), etc.

2)  function_value_second_order_W_iter_scale_step_subspace_dimensionRRSD.mat
				
   	It stores the results of the RRSD algorithm with randomized subspace selection method.

3)  function_value_second_order_W_iter_scale_step_subspace_dimensionRGSD.mat
				
   	It stores the results of the RGSD algorithm with greedy subspace selection method.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

Options with the code 'RSD_extension.m'

-------------------------------------------

1) stepsz_method: The stepsize can be selected adaptively (stepsz_method='adaptive') or fixed to a constant value (stepsz_method='constant'). In the case where the stepsize is constant, the user must provide the specific constant value to the variable 'stepsize'.

2) subspace_dimension: It indicates the cardinality of the chosen subspace and can be set either to 'multi' for a multidimensional subspace or to 'uni' for a unidirectional subspace. However, the 'uni' option does not work for the greedy case.

3) scale_step : In cases where the function grows faster than its quadratic approximation, the scaling factor 'scale_step' should be set to an appropriate value to reduce the step sizes selected by the adaptive method.

4) Algo : Subspaces can be selected using either the greedy method (algo='greedy') or the randomized method (algo='randomized').

---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Running the code

-------------------------------------------

1) Specify the constant matrices Cp, Dq, Ar, Hr, Fs, Gs, Pm, Qm in the file 'user_input.m'. Additionally, provide the functional form of the function 'g,' as defined in the paper, in terms of the following function blocks: trace(Cp * inverse(X)), trace(Dq * X), logdetX, trace(Ar * X * Hr * X), trace(Fs * inverse(X) * Gs * inverse(X)), and trace(Pm * X * Qm *inverse(X)).
By default, these matrices are set to randomly generated positive (semi)definite matrices, and the functional form of 'g' is defined as f_user = trace(Cp * inverse(X)) + trace(Dq * X) - logdetX + trace(Ar * X * Hr * X) + trace(Fs * inverse(X) * Gs * inverse(X)) + trace(Pm * X * Qm * inverse(X)).

2) In the main file 'RSD_extension.m', set the following variables :'stepsz_method' (stepsize selection method), 'scale_step' (scaling factor of stepsize), 'subspace_dimension' (direction of subspace chosen), 'algo' (greedy or randomized algorithm).
Default settings for these variables are: stepsz_method='adaptive', scale_step=1, subspace_dimension='multi', algo='greedy'.

3) Run the code file 'RSD_extension.m' and provide values for 'W' (size of matrices), 'iter'(no of iterations) when prompted.


-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Performance of the algorithm

----------------------------------------------------------------------------------------------------------------------------------------------------------------
For comparison purposes, an alternative method using the MATLAB function 'fmincon' is implemented in the file 'fun_min_con.m'. The positive definiteness constraint is enforced using the change of variable X = L * L' + eps * I, where eps = 2.220446049250313e-16, and 'I' is the identity matrix. Here, 'X' is the variable of the original problem, and 'L' is the new variable used in the code file 'fun_min_con.m'.  To ensure comparable performance, we used 3000 iterations for the proposed RGSD algorithm, while 6000 iterations were used for 'fun_min_con.m'. Here, 'W' denotes the matrix size, 'f' denotes the function value at the last iteration, and 'time' represents the time taken by the algorithm in seconds. Subscripts 'RRSD' and 'fmincon' are used to distinguish the results of the algorithms. The results obtained are given below. To reproduce these results, first run the code 'RSD_extension.m' with default setting and set W to 10, 20, or 50, and set iter to 3000 when prompted. This will produce the results of the RRSD algorithm and store user data, i.e., constant matrices Cp, Dq, Ar, Hr, Fs, Gs, Pm, Qm, in the file 'user_data.mat'. Next, run the code file 'fun_min_con.m'. It will use the 'user_data.mat' file and automatically detect the size of matrices from the file.

---------------------------------------------------------------------------------------------------------------------------------------------------------------
W=10
---------------------------------------------------------------------------------------------------------------------------------------------------------------
f_RGSD=9.750333900089568e+02;

time_RGSD=10.120367699999996;

f_fmincon=9.750334873443385e+02;

time_fmincon=10.836257;

f_fmincon-f_RGSD =9.733538172440603e-05
--------
W=20
--------
f_RGSD   = 9.201035544563601e+03;

time_RGSD   = 12.666394699999998;

f_fmincon= 9.359874987444849e+03;

time_fmincon= 20.420358

f_fmincon-f_RGSD=1.588394428812480e+02
--------
W=50
--------
f_RGSD=1.282001654616017e+05;

time_RGSD=23.138362700000012;

f_fmincon=2.440208588165893e+05;

time_fmincon=37.116243;

f_fmincon-f_RGSD=1.158206933549876e+05
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


