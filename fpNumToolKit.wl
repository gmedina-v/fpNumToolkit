(* ::Package:: *)

(* Automatically call the module 'flowEqs' with the functions to initialise the flow equation *)
(* Comment out this line if the module is not installed. *)
(* <<"flowEqs`" *)
		
betaFunctions::usage =
"betaFunctions[eq,boundaryConditions,approxN]
Computes the symbolic beta functions from a given flow equation with the specified boundary conditions at the given order.

Parameters
----------------------------------------------------
eq: Flow equation. Must be a list containing the coefficients of the Taylor expansion of the flow equation in small curvature.
boundaryConditions: Boundary conditions of couplings to neglect. A list of 2 elements. Can be left as symbolic variables.
approxN: Approximation order, i.e., number of beta functions to compute. Must be smaller than or equal to the number of elements in eq.

Returns
----------------------------------------------------
A list containing the general beta functions

Notes
----------------------------------------------------
This function should only be used for illustrative purposes at low approximation orders.
Using this for large approximation orders is not recommended.
"
betaFunctions[eq_,boundaryConditions_,approxN_]:=
	Block[{u,v,matV,matU,bcValues,betas},
		u[i_]:=-(eq[[i+1]]/.D[\[Lambda][_][t],t]->0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_]:=-Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		bcValues={\[Lambda][approxN][t]->boundaryConditions[[1]],\[Lambda][approxN+1][t]->boundaryConditions[[2]]};
		matV=Table[v[i,j]Boole[i!=j],{i,0,approxN-1},{j,0,approxN-1}];
		matU=Table[u[i],{i,0,approxN-1}];
		betas=Inverse[IdentityMatrix[approxN]-(matV/.bcValues)].(matU/.bcValues);
		Return[betas];
		];

valFP::usage = 
"valFP[{approxN,flowEq,fp0,boundary_:{0,0},options_:{maxIter_:100,wp_:100,fvalTol_:10^-60,jumpTol_:0.03,ignoreLastN_:2}}]
Find a solution of the flow equation at a given order N.
Returns: A list of (3) elements with: solution, residual, test pass indicator

Parameters
----------------------------------------------------
approxN: The approximation order N
flowEq: A list of (approxN) elements. The coefficients of the flow equation as a polynomial of the scalar curvature R. Beta functions must be set to zero
fp0: A list of (approxN) elements. Initial guess to start looking for the solution
boundary: A list of (2) elements. The values of the two highest couplings (boundary conditions)
The options is a list of five elements in the following order:
maxIter: Maximum number of iterations for the root search algorithm
wp: Working Precision
fvalTol: Tolerance of the residuals
jumpTol: Tolerance for the relative difference between the solution found and the initial guess
ignoreLastN: Number of couplings to ignore in the jump test from highest to lowest

Returns
----------------------------------------------------
solution: A list of (approxN) elements. The numeric solution found by the root search algorithm within the limitations set by maxIter
fval: A list of (approxN) elements. The equations evaluated at the values of the solution, i.e. the residuals
test pass indicator: A boolean stating if the found solution passed the consistency tests (see function 'runTests')"

valFP[{approxN_,flowEq_,fp0_,boundary_:{0,0},options_:{100,100,10^-60,0.03,2}}]:=
(* valFP[{approxN_,flowEq_,fp0_,boundary_:{0,0},{maxIter_,wp_,fvalTol_,jumpTol_,ignoreLastN_}_:{100,100,10^-60,0.03,2}}]:= *)
	Block[{maxIter,wp,fvalTol,jumpTol,ignoreLastN,eq,initialGuess,boundaryConditions,root,fval,passTests,testResults},
		(* Read options *)
		{maxIter,wp,fvalTol,jumpTol,ignoreLastN} = options;
		(* Define boundary conditions and initial guess for the root search *)
		boundaryConditions = {\[Lambda][approxN][t] -> boundary[[1]], \[Lambda][approxN + 1][t] -> boundary[[2]]};
		initialGuess = Table[{\[Lambda][j][t],fp0[[j+1]]},{j,0,approxN-1}];
		(* eq = (#==0)&/@(flowEq/.boundaryConditions); *)
		eq = flowEq/.boundaryConditions;
		Quiet[root = FindRoot[eq,initialGuess,WorkingPrecision -> wp,MaxIterations -> maxIter][[;;,2]];];
		(* Run tests to verify if the solution is good. Tests returns True for PASSED and False for FAILED *)
		{fval,testResults} = runTests[approxN,eq,root,fp0,fvalTol,jumpTol,ignoreLastN];
	 	passTests = AllTrue[testResults,TrueQ];
	 	(* If[AllTrue[testResults,TrueQ],
			passTests = True;,
			passTests = False;
			]; *)
		Return[{root,fval,passTests}];
		];

runTests::usage = 
"runTests[approxN,eq,root,fp0,fvalTol,jumpTol,ignoreLastN]
Run consistency tests given a proposed solution at order N and a reference solution at order N-1.
Returns: A list of (2) elements with: fval and {testFval,testGaussian,testJump}

Parameters
----------------------------------------------------
approxN: The approximation order N
eq: A list of (approxN) elements. The coefficients of the flow equation as a polynomial of the scalar curvature R. Beta functions must be set to zero
root: A list of (approxN) elements. The numeric solution to test
fp0: A list of (approxN-1) elements. The reference solution to test against
fvalTol: Tolerance of the residuals
jumpTol: Tolerance for the relative difference between the solution found and the initial guess
ignoreLastN: Number of couplings to ignore in the jump test from highest to lowest

Returns
----------------------------------------------------
fval: A list of (approxN) elements. The equations evaluated at the values of the solution, i.e. the residuals
testFval: A boolean for the result of the residuals test. Returns True for passed
testGaussian: A boolean for the result of the Gaussian test. Returns True for passed
testJump: A boolean for the result of the jump test. Returns True for passed

Notes
----------------------------------------------------
There are three consistency tests:
1. Verify that the sum of the absolute value of the residuals of the proposed solution are smaller than a given tolerance
2. Verify that the proposed solution has not vanishing value for the leading two couplings
3. Verify that the relative difference between the coupling values of the proposed solution against the reference solution is not larger than a tolerance, 
	excluding a number of higher order couplings. The geometric mean is used to average accross the difference of each coupling."

runTests[approxN_,eq_,root_,fp0_,fvalTol_,jumpTol_,ignoreLastN_]:=
	Block[{fval,uptoLambda,prevFP,newFP,testFval,testGaussian,testJump},
		(* The tests return True for PASSED and False for FAILED *)
		(* 1. Is it really is a solution of the objective functions? *)
		(* fval = eq[[;;,1]]/.Table[\[Lambda][j][t] -> root[[j+1]],{j,0,approxN-1}]; *)
		fval = eq/.Table[\[Lambda][j][t] -> root[[j+1]],{j,0,approxN-1}];
		testFval = Not[Total[Abs[fval]] > fvalTol];
		(* 2. Is it Gaussian in any of the first two couplings? *)
		testGaussian = Not[MatchQ[Chop[root],{a_,b_,___}/;a==0||b==0]];
		(* 3. Is the jump from the testSol too large? *)
		(* Be careful that some of the couplings may be practically zero. So the relative change goes to infinity *)
		uptoLambda = -Min[2+ignoreLastN,approxN];
		If[MemberQ[Abs[Chop[fp0[[;;uptoLambda]]]],0],
			prevFP = Delete[fp0[[;;uptoLambda]],Position[Abs[Chop[fp0[[;;uptoLambda]]]],0]];
			newFP = Delete[root[[;;uptoLambda]],Position[Abs[Chop[fp0[[;;uptoLambda]]]],0]];,
			(* else *)
			prevFP = fp0[[;;uptoLambda]];
			newFP = root[[;;uptoLambda]];
			];
		testJump = Not[GeometricMean[Abs[newFP/prevFP-1]] > jumpTol];
		Return[{fval,{testFval,testGaussian,testJump}}];
		];

randomSearch::usage = 
"randomSearch[approxN,flowEq,randomFlag,numRandomPoints,valOptions,seedFileName_:'',randFactor_:1,boundaryConditions_:{0,0},filterFlag_:True,fp0_:{}]
Search for a solution of the flow equation at order N using either random points, a given solution, or randomizing given solutions.
Returns: A list of (3) elements with: fp, fval, passTests

Parameters
----------------------------------------------------
approxN: The approximation order N
flowEq: A list of (approxN) elements. The coefficients of the flow equation as a polynomial of the scalar curvature R. Beta functions must be set to zero
randomFlag: A numeric index indicating how to compute the starting points 1: random numbers, 0: read solutions from a file, -1: Randomize solutions from a file, -2. Use solutions passed in the arguments
numRandomPoints: Number of random points or random variations to use
valOptions: A list with (5) elements with the optional parameters for 'valFP': {maxIter,wp,fvalTol,jumpTol,ignoreLastN}
seedFileName: A string with the absolute filename where to read solutions from (only used if randomFlag == 0,-1)
randFactor: A number which specifies the multiplicative factor of the random point
boundaryConditions: A list of (2) elements. The values of the two highest couplings (boundary conditions)
filterFlag: A boolean specifying whether to filter the solutions by coupling values or not
fp0: A list of (approxN) elements with a proposed solution to use

Returns
----------------------------------------------------
fp: A list of (approxN) elements. The numeric solution found by the root search algorithm within the limitations set by maxIter
fval: A list of (approxN) elements. The equations evaluated at the values of the solution, i.e. the residuals
passTests: A boolean stating if the found solution passed the consistency tests (see function 'runTests')

Notes
----------------------------------------------------
The file in seedFileName is expected to have dimensions (z x fpN x N), where z can be any number but the first element of this dimension must contain the solutions to try. 
	This element should be a (fpN x N) list, where fpN is the number of distinct solutions to try and N is the number of orders of each solution. 
	Of all the orders N, only the last one is used. It is expected to be a list of (approxN) elements. If it's not, it is truncated or padded with zeros as necessary.
For randomFlag == 1 the random points are generated from a pseudo random uniform distribution between -1 and 1 and multiplied by randFactor.
For randomFlag == -1 a relative variation is applied to the solution, where the relative difference is given by a pseudo random number from the uniform distribution between -1 and 1.
A small imaginary part is added to all solutions to enable search in the complex plane.
The filter by coupling values is implemented calling 'filterCouplings', which can be adjusted to filter by all tests, or only by the residuals test."

randomSearch[approxN_,flowEq_,randomFlag_,numRandomPoints_,valOptions_,seedFileName_:"",randFactor_:1,boundaryConditions_:{0,0},filterFlag_:True,fp0_:{}] :=
	Block[{randomFP0,fp,fval,passTests,fvalTol},
		(* Choose which way *)
		Switch[randomFlag,
			1,  (* Random points *)
				randomFP0 = randFactor*RandomReal[{-1,1},{numRandomPoints,approxN}];,
			0,  (* Use given solutions *)
				randomFP0 = Import[seedFileName][[1,;;,-1]];
				randomFP0 = PadRight[randomFP0,{Length[randomFP0],approxN},0];,
			-1, (* Randomize given solutions *)
				randomFP0 = Import[seedFileName][[1,;;,-1]];
				randomFP0 = PadRight[randomFP0,{Length[randomFP0],approxN},0];
				randomFP0 = Flatten[(ConstantArray[#,numRandomPoints]&/@randomFP0)*RandomReal[{-1,1},{Length[randomFP0],numRandomPoints,approxN}],1];,
			-2, (* Use a solution passed in the inputs *)
				randomFP0 = PadRight[fp0,{Length[fp0],approxN},0];
			];
		randomFP0 = randomFP0 + 10^-20 I;

		(* Call valFP *)
		{fp,fval,passTests} = Transpose[Map[valFP[{approxN,flowEq,#,boundaryConditions,valOptions}]&,randomFP0]];
		
		(* Filter by couplings? *)
		If[TrueQ[filterFlag],
			fvalTol = valOptions[[3]];
			{fp,fval,passTests} = filterCouplings[fp,fval,passTests,fvalTol];
			];

		(* Reshape arrays to dimenstions (fpN x 1 approxN) *)
		fp = ArrayReshape[fp,{Length[fp],1,approxN}];
		fval = ArrayReshape[fval,{Length[fp],1,approxN}];
		passTests = ArrayReshape[passTests,{Length[fp],1}];
		Return[{fp,fval,passTests}];
		];

filterCouplings::usage = 
"filterCouplings[fp,fval,passTests,fvalTol_:10^-60,fvalFlag_:True,testFlag_:True]
Apply filters to a list of solutions of the flow equaiton, including at the very least deleting duplicates.
Returns: A (3 x m) list with: fp, fval, passTests for each of the m solutions which passed the test.

Parameters
----------------------------------------------------
fp: A list of (n) elements with the solutions to filter
fval: A list of (n) elements with the residuals for each solution
passTests: A list of (n) elements with the results of 'valFP' for each solution
fvalTol: Tolerance of the residuals
fvalFlag: A boolean specifying whether to filter by residuals
testFlag: A boolean specifying whether to filter by the result of passTests

Returns
----------------------------------------------------
{fp, fval, passTests} as stated above but with the (m) elements that passed the tests

Notes
----------------------------------------------------
Solutions are counted as duplicates if they match in the first 50 digits."

filterCouplings[fp_,fval_,passTests_,fvalTol_:10^-60,fvalFlag_:True,testFlag_:True] := 
	Block[{indOk,indFval,indPass},
		indOk = Flatten[Map[Position[N[Chop[fp],50],#,\[Infinity],1]&,DeleteDuplicates[N[Chop[fp],50]]]];
		If[fvalFlag,
			indFval = Flatten[Position[Not/@(#>fvalTol&/@(Total/@Abs[fval])),True]];
			indOk = Intersection[indOk,indFval];
			];
		If[testFlag,
			indPass = Flatten[Position[passTests,True]];
			indOk = Intersection[indOk,indPass];
			];
		Return[Table[Extract[m,n],{m,{fp,fval,passTests}},{n,indOk}]];
		];

eigSystem::usage = 
"eigSystem[eq,fp,boundaryConditions,minN,maxN,parallel_:False]
Compute eigenvalues and eigenvectors of a solution in a range of approximations.
Returns: A (fp x N x 2) list with: eigenvalues and eigenectors for each solution at each order N.

Parameters
----------------------------------------------------
eq: A list of (maxN) elements. The coefficients of the flow equation as a polynomial of the scalar curvature R. Must include beta functions
fp: A list of (fpN x maxN-minN+1, approxN_i) elements with the solutions. Each element is a list with the corresponding number of elements of (approxN)
boundaryConditions: A list of (fpN x maxN-minN+1, 2) elements with the boundary conditions. Each element is a list of (2) elements with the values of highest order couplings
minN: The lowest order at which to compute the eigensystem
maxN: The highest order at which to compute the eigensystem
parallel: A boolean specifying whether to perform parallel evaluation

Returns
----------------------------------------------------
eig: A list of (fp x maxN-minN+1 x 2, approxN_i) elements. The first element of the third dimension contains a list with the eigenvalues, and the second one, the eigenvectors.
	Note the order does matter, so if sorting by eigenvalues, the sorting should also be considered for the eigenvectors.

Notes
----------------------------------------------------
The dimenstions of the inputs fp and boundaryConditions are mandatory, even if the computation is being done for a single solution at a single order.
In that case, the inputs can have dimensions (1 x 1 x approxN) and (1 x 1 x 2), respectively, where approxN is the order at which the eigensystem will be computed."

eigSystem[eq_,fp_,boundaryConditions_,minN_,maxN_,parallel_:False] :=
	Block[{dims,fpN,numN,nsteps,u,v,matV,matdU,fpValues,bcValues,eig,temp1,temp2,tempValues},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)
		dims = Dimensions/@{fp,boundaryConditions};
		{fpN,numN} = Transpose[dims[[;;,;;2]]];

		If[Not[And @@ (Equal @@ # & /@ {fpN,numN})],
			Print["Inputs are not the same dimensions. Check inputs."];
			Abort[];,
			{fpN,numN} = First/@{fpN,numN};
			];

		If[numN != (maxN-minN+1),
			Print["Input does not have the correct number of orders requested N, it should be (maxN - minN + 1)."];
			Abort[];
			];

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		Print["Create matV: "<>ToString[AbsoluteTiming[matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Create matdU: "<>ToString[AbsoluteTiming[matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Make matV a function: "<>ToString[AbsoluteTiming[matV = Function@@{matV}/.Dispatch@Thread[Array[\[Lambda],maxN+2,{0,maxN+1}]->Array[Slot,maxN+2,{1,maxN+2}]];]]];
		Print["Make matdU a function: "<>ToString[AbsoluteTiming[matdU = Function@@{matdU}/.Dispatch@Thread[Array[\[Lambda],maxN+2,{0,maxN+1}]->Array[Slot,maxN+2,{1,maxN+2}]];]]];
		tempValues = MapThread[PadRight[Join[#1,#2],maxN+2]&,{fp,boundaryConditions},2];

		If[TrueQ[parallel],
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = ParallelMap[matV@@#&,tempValues,{2}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{2}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = ParallelMap[matdU@@#&,tempValues,{2}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = ParallelTable[-Eigensystem[(temp1[[m,i]]).(temp2[[m,i]])],{m,fpN},{i,numN}];]]];
			(* eig: {fpN, totalN, 2} *)
			(* eig[[fpN, totalN, 1]]: {minN:maxN} *)
			(* eig[[fpN, totalN, 2]]: {minN:maxN, minN:maxN} *)
			,
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = Map[matV@@#&,tempValues,{2}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{2}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = Map[matdU@@#&,tempValues,{2}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = Table[-Eigensystem[(temp1[[m,i]]).(temp2[[m,i]])],{m,fpN},{i,numN}];]]];
			];
		Return[eig];
		];

eigSystemExtraBC::usage = 
"
eigSystemExtraBC[eq,fp,boundaryConditions,minN,maxN,parallel:False]
Same as eigSystem but allows for a total of up to 6 boundary conditions. Suitable when the flow equation is expanded around non-zero curvature.
"
eigSystemExtraBC[eq_,fp_,boundaryConditions_,minN_,maxN_,parallel_:False] :=
	Block[{dims,fpN,numN,nsteps,u,v,matV,matdU,fpValues,bcValues,eig,temp1,temp2,tempValues},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)
		dims = Dimensions/@{fp,boundaryConditions};
		{fpN,numN} = Transpose[dims[[;;,;;2]]];

		If[Not[And @@ (Equal @@ # & /@ {fpN,numN})],
			Print["Inputs are not the same dimensions. Check inputs."];
			Abort[];,
			{fpN,numN} = First/@{fpN,numN};
			];

		If[numN != (maxN-minN+1),
			Print["Input does not have the correct number of orders requested N, it should be (maxN - minN + 1)."];
			Abort[];
			];

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		Print["Create matV: "<>ToString[AbsoluteTiming[matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Create matdU: "<>ToString[AbsoluteTiming[matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Make matV a function: "<>ToString[AbsoluteTiming[matV = Function@@{matV}/.Dispatch@Thread[Array[\[Lambda],maxN+6,{0,maxN+5}]->Array[Slot,maxN+6,{1,maxN+6}]];]]];
		Print["Make matdU a function: "<>ToString[AbsoluteTiming[matdU = Function@@{matdU}/.Dispatch@Thread[Array[\[Lambda],maxN+6,{0,maxN+5}]->Array[Slot,maxN+6,{1,maxN+6}]];]]];
		tempValues = MapThread[PadRight[Join[#1,#2],maxN+6]&,{fp,boundaryConditions},2];

		If[TrueQ[parallel],
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = ParallelMap[matV@@#&,tempValues,{2}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{2}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = ParallelMap[matdU@@#&,tempValues,{2}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = ParallelTable[-Eigensystem[(temp1[[m,i]]).(temp2[[m,i]])],{m,fpN},{i,numN}];]]];
			(* eig: {fpN, totalN, 2} *)
			(* eig[[fpN, totalN, 1]]: {minN:maxN} *)
			(* eig[[fpN, totalN, 2]]: {minN:maxN, minN:maxN} *)
			,
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = Map[matV@@#&,tempValues,{2}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{2}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = Map[matdU@@#&,tempValues,{2}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1]],{1,i},{1,i}],{m,fpN},{i,minN,maxN}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = Table[-Eigensystem[(temp1[[m,i]]).(temp2[[m,i]])],{m,fpN},{i,numN}];]]];
			];
		Return[eig];
		];

eigFilter::usage = 
"eigFilter[maxN,fp,fval,passTests,eig,evNCompare_:2,evRankCompare_:2,maxEvVar_:1,largeEvTol_:3]
Filters solutions based on the eigenspectrum according to consistency tests
Returns: A (4 x fpN) list with: {fp,fval,passTests,eig} that pass the filters.

Parameters
----------------------------------------------------
maxN: The highest available approximation order of the solutions to filter
fp: A list of (fpN x N x approxN_i) elements with the coupling values to filter
fval: A list of (fpN x N x approxN_i) elements with the residuals to filter
passTests: A list of (fpN x N) elements with the results of coupling consistency tests to filter
eig: A list of (fpN x Max[2,maxN] x 2 x approxN_i ) elements with the eigensystem to which the tests will be applied
evNCompare: The number of approximation orders N to compare
evRankCompare: The number of leading eigenvalues to compare
maxEvVar: Maximum allowed relative variation in the eigenvalues between two consecutive orders
largeEvTol: Factor controlling the tolerance for large numerical eigenvalues

Returns
----------------------------------------------------
values: A list of (4 x fpN) elements with the filtered inputs {fp,fval,passTests,eig}.

Notes
----------------------------------------------------
This function performs two tests to the eigensystem:
1. Large eigenvalue tests. Solutions with large eigenvalues are removed. The tolerance for a large eigenvalue is: largeEvTol x maxN. e.g. at order N = 5 with largeEvTol = 3 the largest
	absolute value of anny eigenvalue that is allowed is 5 x 3 = 15. Solutions with eigenvalues larger than that will be dropped.
2. Relative difference test. The relative difference between the leading eigenvalues across different orders must be small. The number of eigenvalues to compare is controlled
	directly by evRankCompare. The number of orders to compare, by evNCompare. The absolute value and the mean of the relative difference is taken to determine the outcome.
	The maximum tolerance is given by maxEvVar. e.g. maxEvVar = 1 indicates that a relative difference of up to 100% on average is allowed."

eigFilter[maxN_,fp_,fval_,passTests_,eig_,evNCompare_:2,evRankCompare_:2,maxEvVar_:1,largeEvTol_:3] :=
	Block[{temp1,temp2,indLarge,indVar,values},
		temp1 = Sort[#,Re[#1]>Re[#2]&]&/@eig[[;;,-1,1]];
		If[Not[evNCompare>=maxN || evRankCompare-1>=maxN] && evNCompare>1,
			temp2 = Re[(Map[Sort[#,Re[#1]>Re[#2]&]&,eig[[;;,;;,1]],2])[[;;,Max[evRankCompare-1,maxN-evNCompare];;,;;evRankCompare]]];
			indLarge = Position[AnyTrue[#,#>largeEvTol*maxN&]&/@Abs[temp1],True];
			indVar = Position[AnyTrue[Abs[#],#>maxEvVar&]&/@Map[Mean,Transpose[Ratios[#]&/@temp2-1,{1,3,2}],{2}],True];
			values = Delete[#,Union[indLarge,indVar]]&/@{fp,fval,passTests,eig};
			log[StringJoin[ToString[Length[indLarge]]," solutions removed due to large eigenvalues, ",ToString[Length[indVar]],
				" removed due to high variation (possible doublecounting). ",ToString[Length[values[[1]]]]," solutions remain."]];,
			log["Not enough orders to filter by eigenvalues. Change parameters. Abort."];
			Abort[];
			];
		(* values = {fp,fval,passTests,eig} (clean) *)
		Return[values];
		];
		
verticalSearchNFP::usage = 
"verticalSearchNFP[minN,maxN,flowEq,fp1,bc_:{0,0},valOptions_:{100,100,10^-60,0.03,2},imagPart_10^-20 I]
This function performs a vertical search from minN to maxN using as starting point fp1.

Parameters
----------------------------------------------------
minN: Starting approximation order
maxN: Ending approximation order
flowEq: A list containing the terms of the flow equation with beta functions set to zero.
fp1: The initial solution at order minN (minN)
bc: Boundary conditions

Returns
----------------------------------------------------
A list of (maxN-minN+1 x 3) elements"

verticalSearchNFP[minN_,maxN_,flowEq_,fp1_,bc_:{0,0},valOptions_:{100,100,10^-60,0.03,2},imagPart_:10^-20 I] := 
	Block[{values,fp0,tempFlowEq,fp,fval,passTests},
		If[Not[Length[flowEq] == maxN],
			Print["Flow equation does not have the correct number of elements."];
			Abort[];
			];
		If[Not[Length[fp1] == minN],
			Print["Initial fixed point does not have the correct number of elements."];
			Abort[];
			];

		fp0 = fp1 + imagPart;
		values = Reap[Do[
			tempFlowEq = flowEq[[;;approxN]];
			{fp,fval,passTests} = valFP[{approxN,tempFlowEq,fp0,bc,valOptions}];
			fp0 = Join[fp,{0}] + imagPart;
			Sow[{fp,fval,passTests}];
			,{approxN,minN,maxN}];
		][[2,1]];
		(* values = Partition[#,maxN-minN+1]&/@Transpose[values]; *)
		Return[values];
		];

fpSolveVBC::usage = 
"fpSolveVBC[approxN,flowEq,fp1,bc,p1,p2,nsteps,valOptions_:{100,100,10^-60,0.1,2},history_:True,spacing_:'linear',logPrecision_:3,imagPart_:10^-20 I]
This function performs a smooth deformation of the solution between two boundary coordinates p1 -> p2

Parameters
----------------------------------------------------
flowEq must be a list of approxN elements with beta functions set to 0
fp1 must be a list of approxN elements
bc must be a list of two elements
p1 and p2 must each be a list of two elements
nsteps must be a positive integer larger than 1

Returns
----------------------------------------------------
Returns a list of (3) elements if history is set to False (default)
Returns a list of (nsteps x 3) elements if history is set to True

Notes
----------------------------------------------------
spacing: {'linear','logSparseToDense','logDenseToSparse'}
logSparsetoDense: The initial point is exact, the end point is approximate 
logDenseToSparse: The initial point is approximate, the end point is exact
"
fpSolveVBC[approxN_,flowEq_,fp1_,bc_,p1_,p2_,nsteps_,valOptions_:{100,100,10^-60,0.1,2},history_:True,spacing_:"linear",logPrecision_:3,imagPart_:10^-20 I] :=
	Block[{warning,logSpacing,r,tempBC,fp0,fp,fval,passTests,values},
		If[Not[Length[flowEq] == approxN],
			Print["Flow equation does not have the correct number of elements."];
			Abort[];
			];
		If[Not[Length[fp1] == approxN],
			Print["Initial fixed point does not have the correct number of elements."];
			Abort[];
			];

		warning = 0;
		(* Straight line *)
		(* lineEq = LeastSquares[{{#1[[1]],1},{#2[[1]],1}},{#1[[2]],#2[[2]]}].{x,1}&@@{p1,p2}; *)
		(* In logSpacing, c controls the precision. 
			c = 1: lower precision at the dense point and more evenly distributed points; 
			c = high: higher precision at the dense point and less evenly distributed points; *)
		logSpacing[a_,b_,c_] := (10^Subdivide[Log10[10^-c],Log10[1],Max[nsteps-2,2]])(b-a)+a;
		(* Choose spacing between p1 and p2 *)
		Switch[spacing,
			"linear",
				(* Straight line *)
				points = Subdivide[p1,p2,Max[nsteps-1,2]];,
			"logSparseToDense",
				(* Log spacing sparse to dense *)
				(* The initial point is exact, the end point is approximate *)
				points = Reverse[Prepend[MapThread[{#1,#2}&,{logSpacing[p2[[1]],p1[[1]],logPrecision],logSpacing[p2[[2]],p1[[2]],logPrecision]}],p2]];,
			"logDenseToSparse",
				(* Log spacing dense to spare *)
				(* The initial point is approximate, the end point is exact *)
				points = Prepend[MapThread[{#1,#2}&,{logSpacing[p1[[1]],p2[[1]],logPrecision],logSpacing[p1[[2]],p2[[2]],logPrecision]}],p1];
				,
			_, Print["Spacing option not valid. Aborting."];
			Abort[];
			];
		tempBC = bc*#&/@points;

		fp0 = fp1;
		values = Reap[Do[
			{fp,fval,passTests} = valFP[{approxN,flowEq,fp0,tempBC[[n]],valOptions}];
			(* Print a warning if the solution was real and is now complex *)
			If[AnyTrue[Im[Chop[fp0]],#==0&] && AnyTrue[Im[Chop[fp]],#!=0&],
				Print["Warning: solution has drifted into the complex plane"];
				warning = 1;
				];
			(* Print a warning if the solution was complex and is now real, and set the test indicator to FAIL *)
			If[AnyTrue[Im[Chop[fp0]],#!=0&] && AnyTrue[Im[Chop[fp]],#==0&],
				Print["Critical warning: solution has emerged from the complex plane into the real axis. This solution cannot be trusted anymore."];
				warning = 2;
				passTests = False;
				];
			fp0 = fp + imagPart;
			Sow[{fp,fval,passTests,warning}];
			,{n,nsteps}]][[2,1]];
		values = MapThread[Join[#1, {#2}] &, {values, points}];
		If[history,
			(* Return the whole history *)
			Return[values];,
			(* Return only the last point *)
			Return[values[[-1]]];
			];
		];
		
histSearch::usage = "
histSearch[minN,maxN,flowEq,fp0,valOptions,nsteps:100]
This function performs a history search calling fpSolveVBC, using a high N solution to look for low N solutions

Parameters
----------------------------------------------------
minN and maxN are the search orders
flowEq must be a list of maxN elements with beta functions set to 0
fp0 must be a list of (fpN x maxN) elements

Returns
----------------------------------------------------
Returns a list of (3 x numFP x maxN-minN+1) elements
"

histSearch[minN_,maxN_,flowEq_,fp0_,valOptions_,nsteps_:100] :=
	Block[{numFP,values,optFP,optBC,fp,fval,passTests},
		numFP = Length[fp0];
		values = Reap[Do[
			tempFlowEq = flowEq[[;;approxN]];
			optFP = fp0[[fpN,;;approxN]];
			optBC = fp0[[fpN,approxN+1;;approxN+2]];
			{fp,fval,passTests} = fpSolveVBC[approxN,tempFlowEq,optFP,optBC,{1,1},{0,0},nsteps,valOptions];
			Sow[{fp,fval,passTests}];
			,{fpN,numFP}
			,{approxN,minN,maxN}]][[2,1]];
		values = Partition[#,maxN-minN+1]&/@Transpose[values];
		Return[values];
		];

varMatter::usage = 
"
varMatter[approxN,flowEq,fp1,matterContent,p1,p2,nsteps,bc_:{0,0},valOptions_:{100,100,10^-60,0.1,2},spacing_:'linear',logPrecision_:3,imagPart_:10^-20 I]
Looks for a fixed point variating the matter content by tracing a straight line from p1 to p2
Returns a list of three elements with fp, fval and passTests
Flow equation must be in terms of ns, nm and nd
"
varMatter[approxN_,flowEq_,fp1_,matterContent_,p1_,p2_,nsteps_,bc_:{0,0},valOptions_:{100,100,10^-60,0.1,2},spacing_:"linear",logPrecision_:3,imagPart_:10^-20 I] := 
	Block[{logSpacing,points,tempMat,tempFlowEq,fp0,values,fp,fval,passTests},
		If[Not[Length[flowEq] == approxN],
			Print["Flow equation does not have the correct number of elements."];
			Abort[];
			];
		If[Not[Length[fp1] == approxN],
			Print["Initial fixed point does not have the correct number of elements."];
			Abort[];
			];

		(* Straight line *)
		(* lineEq = LeastSquares[{{#1[[1]],1},{#2[[1]],1}},{#1[[2]],#2[[2]]}].{x,1}&@@{p1,p2}; *)
		(* In logSpacing, c controls the precision. 
			c = 1: lower precision at the dense point and more evenly distributed points; 
			c = high: higher precision at the dense point and less evenly distributed points; *)
		logSpacing[a_,b_,c_] := (10^Subdivide[Log10[10^-c],Log10[1],Max[nsteps-2,2]])(b-a)+a;
		Switch[spacing,
			"linear",
				(* Straight line *)
				points = Subdivide[p1,p2,Max[nsteps-1,2]];,
			"logSparseToDense",
				(* Log spacing sparse to dense *)
				(* The initial point is exact, the end point is approximate *)
				points = Reverse[Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p2[[1]],p1[[1]],logPrecision],logSpacing[p2[[2]],p1[[2]],logPrecision],logSpacing[p2[[3]],p1[[3]],logPrecision]}],p2]];,
			"logDenseToSparse",
				(* Log spacing dense to spare *)
				(* The initial point is approximate, the end point is exact *)
				points = Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p1[[1]],p2[[1]],logPrecision],logSpacing[p1[[2]],p2[[2]],logPrecision],logSpacing[p1[[3]],p2[[3]],logPrecision]}],p1];
				,
			_, Print["Spacing option not valid. Aborting."];
			Abort[];
			];
		tempMat = matterContent*#&/@points;

		fp0 = fp1;
		values = Reap[Do[
			tempFlowEq = flowEq/.{ns->tempMat[[n,1]],nm->tempMat[[n,2]],nd->tempMat[[n,3]]};
			{fp,fval,passTests} = valFP[{approxN,tempFlowEq,fp0,bc,valOptions}];
			fp0 = fp + imagPart;
			Sow[{fp,fval,passTests}];
			,{n,nsteps}]][[2,1]];
		values = MapThread[Join[#1, {#2}] &, {values, points}];
		Return[values];
		]

eigSysVarMat::usage = 
"eigSysVarMat[eq,fp,boundaryConditions,matter,minN,maxN,parallel_:False]
Compute eigenvalues and eigenvectors of a solution in a range of approximations.
Returns: A (fp x N x matter x 2) list with: eigenvalues and eigenectors for each solution at each order N and matter configuration.

Parameters
----------------------------------------------------
eq: A list of (maxN) elements. The coefficients of the flow equation as a polynomial of the scalar curvature R. Must include beta functions
fp: A list of (fpN x maxN-minN+1, approxN_i) elements with the solutions. Each element is a list with the corresponding number of elements of (approxN)
boundaryConditions: A list of (fpN x maxN-minN+1, 2) elements with the boundary conditions. Each element is a list of (2) elements with the values of highest order couplings
matter: A list of (m x 3) elements with the m matter configurations to try, given by (ns, nm, nd) in that order.
minN: The lowest order at which to compute the eigensystem
maxN: The highest order at which to compute the eigensystem
parallel: A boolean specifying whether to perform parallel evaluation

Returns
----------------------------------------------------
eig: A list of (fp x maxN-minN+1 x # of matter configurations x 2, approxN_i) elements. The first element of the third dimension contains a list with the eigenvalues, and the second one, the eigenvectors.
	Note the order does matter, so if sorting by eigenvalues, the sorting should also be considered for the eigenvectors.

Notes
----------------------------------------------------
The dimenstions of the inputs fp and boundaryConditions are mandatory, even if the computation is being done for a single solution at a single order.
In that case, the inputs can have dimensions (1 x 1 x approxN) and (1 x 1 x 2), respectively, where approxN is the order at which the eigensystem will be computed. 
This function assumes that all fixed point branches share the same list of matter configurations at every order N"

eigSysVarMat[eq_,fp_,boundaryConditions_,matter_,minN_,maxN_,parallel_:False] :=
	Block[{dims,fpN,numN,nsteps,u,v,matV,matdU,temp1,temp2,tempValues,eig},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)
		dims = Dimensions/@{fp,boundaryConditions,matter};
		{fpN,numN,nsteps} = Transpose[dims[[;;,;;3]]];

		If[Not[And @@ (Equal @@ # & /@ {fpN,numN,nsteps})],
			Print["Inputs are not the same dimensions. Check inputs."];
			Abort[];,
			{fpN,numN,nsteps} = First/@{fpN,numN,nsteps};
			];

		If[numN != (maxN-minN+1),
			Print["Input does not have the correct number of orders requested N, it should be (maxN - minN + 1)."];
			Abort[];
			];

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		(* fpValues: {fpN, totalN, minN:maxN} *)
		(* bcValues: {fpN, totalN, 2} *)

		Print["Create matV: "<>ToString[AbsoluteTiming[matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Create matdU: "<>ToString[AbsoluteTiming[matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Make matV a function: "<>ToString[AbsoluteTiming[matV = Function@@{matV}/.Dispatch@Join[Thread[{ns,nm,nd}->Array[Slot,3]],Thread[Array[\[Lambda],maxN+2,{0,maxN+1}]->Array[Slot,maxN+2,{4,maxN+5}]]];]]];
		Print["Make matdU a function: "<>ToString[AbsoluteTiming[matdU = Function@@{matdU}/.Dispatch@Join[Thread[{ns,nm,nd}->Array[Slot,3]],Thread[Array[\[Lambda],maxN+2,{0,maxN+1}]->Array[Slot,maxN+2,{4,maxN+5}]]];]]];
		tempValues = MapThread[PadRight[Join[#1,#2,#3],maxN+5]&,{matter,fp,boundaryConditions},3];

		If[TrueQ[parallel],
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = ParallelMap[matV@@#&,tempValues,{3}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{3}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = ParallelMap[matdU@@#&,tempValues,{3}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = ParallelTable[-Eigensystem[(temp1[[m,i,q]]).(temp2[[m,i,q]])],{m,fpN},{i,numN},{q,nsteps}];]]];
			,
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = Map[matV@@#&,tempValues,{3}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{3}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = Map[matdU@@#&,tempValues,{3}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = Table[-Eigensystem[(temp1[[m,i,q]]).(temp2[[m,i,q]])],{m,fpN},{i,numN},{q,nsteps}];]]];
			];
		Return[eig];
		];

varGrav::usage = 
"
varGrav[approxN,flowEq,fp1,gravContent,p1,p2,nsteps,bc_:{0,0},valOptions_:{100,100,10^-60,0.1,2},spacing_:'linear',logPrecision_:3,imagPart_:10^-20 I]
Looks for a fixed point variating the gravitational content by tracing a straight line from p1 to p2
Returns a list of three elements with fp, fval and passTests
Flow equation must be in terms of a, b and c
"
varGrav[approxN_,flowEq_,fp1_,gravContent_,p1_,p2_,nsteps_,bc_:{0,0},valOptions_:{100,100,10^-60,0.1,2},spacing_:"linear",logPrecision_:3,imagPart_:10^-20 I] := 
	Block[{logSpacing,points,tempGrav,tempFlowEq,fp0,values,fp,fval,passTests},
		If[Not[Length[flowEq] == approxN],
			Print["Flow equation does not have the correct number of elements."];
			Abort[];
			];
		If[Not[Length[fp1] == approxN],
			Print["Initial fixed point does not have the correct number of elements."];
			Abort[];
			];
		(* In logSpacing, c controls the precision. 
			c = 1: lower precision at the dense point and more evenly distributed points; 
			c = high: higher precision at the dense point and less evenly distributed points; *)
		logSpacing[a_,b_,c_] := (10^Subdivide[Log10[10^-c],Log10[1],Max[nsteps-2,2]])(b-a)+a;
		(* Choose spacing between p1 and p2 *)
		Switch[spacing,
			"linear",
				(* Straight line *)
				points = Subdivide[p1,p2,Max[nsteps-1,2]];,
			"logSparseToDense",
				(* Log spacing sparse to dense *)
				(* The initial point is exact, the end point is approximate *)
				points = Reverse[Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p2[[1]],p1[[1]],logPrecision],logSpacing[p2[[2]],p1[[2]],logPrecision],logSpacing[p2[[3]],p1[[3]],logPrecision]}],p2]];,
			"logDenseToSparse",
				(* Log spacing dense to spare *)
				(* The initial point is approximate, the end point is exact *)
				points = Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p1[[1]],p2[[1]],logPrecision],logSpacing[p1[[2]],p2[[2]],logPrecision],logSpacing[p1[[3]],p2[[3]],logPrecision]}],p1];
			];
		tempGrav = gravContent*#&/@points;

		fp0 = fp1;
		values = Reap[Do[
			tempFlowEq = flowEq/.{a->tempGrav[[n,1]],b->tempGrav[[n,2]],c->tempGrav[[n,3]]};
			{fp,fval,passTests} = valFP[{approxN,tempFlowEq,fp0,bc,valOptions}];
			fp0 = fp + imagPart;
			Sow[{fp,fval,passTests}];
			,{n,nsteps}]][[2,1]];
		values = MapThread[Join[#1, {#2}] &, {values, points}];
		Return[values];
		];

eigSysVarGrav::usage = 
"eigSysVarGrav[eq,fp,boundaryConditions,gravity,minN,maxN,parallel_:False]
Compute eigenvalues and eigenvectors of a solution in a range of approximations.
Returns: A (fp x N x gravity x 2) list with: eigenvalues and eigenectors for each solution at each order N and gravitational configuration.

Parameters
----------------------------------------------------
eq: A list of (maxN) elements. The coefficients of the flow equation as a polynomial of the scalar curvature R. Must include beta functions
fp: A list of (fpN x maxN-minN+1 x m(grav), approxN_i) elements with the solutions. Each element is a list with the corresponding number of elements of (approxN)
boundaryConditions: A list of (fpN x maxN-minN+1 x m(grav), 2) elements with the boundary conditions. Each element is a list of (2) elements with the values of highest order couplings
gravity: A list of (m x 3) elements with the m gravitational configurations to try, given by (a, b, c) in that order.
minN: The lowest order at which to compute the eigensystem
maxN: The highest order at which to compute the eigensystem
parallel: A boolean specifying whether to perform parallel evaluation

Returns
----------------------------------------------------
eig: A list of (fp x maxN-minN+1 x # of gravitational configurations x 2, approxN_i) elements. The first element of the third dimension contains a list with the eigenvalues, and the second one, the eigenvectors.
	Note the order does matter, so if sorting by eigenvalues, the sorting should also be considered for the eigenvectors.

Notes
----------------------------------------------------
The dimenstions of the inputs fp and boundaryConditions are mandatory, even if the computation is being done for a single solution at a single order.
In that case, the inputs can have dimensions (1 x 1 x approxN) and (1 x 1 x 2), respectively, where approxN is the order at which the eigensystem will be computed. 
This function assumes that all fixed point branches share the same list of matter configurations at every order N"

eigSysVarGrav[eq_,fp_,boundaryConditions_,gravity_,minN_,maxN_,parallel_:False] :=
	Block[{dims,fpN,numN,nsteps,u,v,matV,matdU,temp1,temp2,tempValues,eig},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)
		dims = Dimensions/@{fp,boundaryConditions,gravity};
		{fpN,numN,nsteps} = Transpose[dims[[;;,;;3]]];

		If[Not[And @@ (Equal @@ # & /@ {fpN,numN,nsteps})],
			Print["Inputs are not the same dimensions. Check inputs."];
			Abort[];,
			{fpN,numN,nsteps} = First/@{fpN,numN,nsteps};
			];

		If[numN != (maxN-minN+1),
			Print["Input does not have the correct number of orders requested N, it should be (maxN - minN + 1)."];
			Abort[];
			];

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		Print["Create matV: "<>ToString[AbsoluteTiming[matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Create matdU: "<>ToString[AbsoluteTiming[matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
		Print["Make matV a function: "<>ToString[AbsoluteTiming[matV = Function@@{matV}/.Dispatch@Join[Thread[{a,b,c}->Array[Slot,3]],Thread[Array[\[Lambda],maxN+2,{0,maxN+1}]->Array[Slot,maxN+2,{4,maxN+5}]]];]]];
		Print["Make matdU a function: "<>ToString[AbsoluteTiming[matdU = Function@@{matdU}/.Dispatch@Join[Thread[{a,b,c}->Array[Slot,3]],Thread[Array[\[Lambda],maxN+2,{0,maxN+1}]->Array[Slot,maxN+2,{4,maxN+5}]]];]]];
		tempValues = MapThread[PadRight[Join[#1,#2,#3],maxN+5]&,{gravity,fp,boundaryConditions},3];

		If[TrueQ[parallel],
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = ParallelMap[matV@@#&,tempValues,{3}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{3}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = ParallelMap[matdU@@#&,tempValues,{3}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = ParallelTable[-Eigensystem[(temp1[[m,i,q]]).(temp2[[m,i,q]])],{m,fpN},{i,numN},{q,nsteps}];]]];
			,
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = Map[matV@@#&,tempValues,{3}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{3}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = Map[matdU@@#&,tempValues,{3}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-minN+1,q]],{1,i},{1,i}],{m,fpN},{i,minN,maxN},{q,nsteps}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[eig = Table[-Eigensystem[(temp1[[m,i,q]]).(temp2[[m,i,q]])],{m,fpN},{i,numN},{q,nsteps}];]]];
			];
		Return[eig];
		];


(* Testing functions *)
(***************************************************************************)
(***************************************************************************)
(* Same as varGrav but prints the iteration number *)
varGrav2[approxN_,flowEq_,fp1_,gravContent_,p1_,p2_,nsteps_,bc_:{0,0},valOptions_:{100,100,10^-60,0.1,2},spacing_:"linear",logPrecision_:3,imagPart_:10^-20 I] := 
	Block[{logSpacing,points,tempGrav,tempFlowEq,fp0,values,fp,fval,passTests},
		If[Not[Length[flowEq] == approxN],
			Print["Flow equation does not have the correct number of elements."];
			Abort[];
			];
		If[Not[Length[fp1] == approxN],
			Print["Initial fixed point does not have the correct number of elements."];
			Abort[];
			];

		(* In logSpacing, c controls the precision. 
			c = 1: lower precision at the dense point and more evenly distributed points; 
			c = high: higher precision at the dense point and less evenly distributed points; *)
		logSpacing[a_,b_,c_] := (10^Subdivide[Log10[10^-c],Log10[1],Max[nsteps-2,2]])(b-a)+a;
		(* Choose spacing between p1 and p2 *)
		Switch[spacing,
			"linear",
				(* Straight line *)
				points = Subdivide[p1,p2,Max[nsteps-1,2]];,
			"logSparseToDense",
				(* Log spacing sparse to dense *)
				(* The initial point is exact, the end point is approximate *)
				points = Reverse[Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p2[[1]],p1[[1]],logPrecision],logSpacing[p2[[2]],p1[[2]],logPrecision],logSpacing[p2[[3]],p1[[3]],logPrecision]}],p2]];,
			"logDenseToSparse",
				(* Log spacing dense to spare *)
				(* The initial point is approximate, the end point is exact *)
				points = Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p1[[1]],p2[[1]],logPrecision],logSpacing[p1[[2]],p2[[2]],logPrecision],logSpacing[p1[[3]],p2[[3]],logPrecision]}],p1];
			];
		tempGrav = gravContent*#&/@points;

		fp0 = fp1;
		values = Reap[Do[
		log["Iteration: "<>ToString[n]<>" of "<>ToString[nsteps]];
			tempFlowEq = flowEq/.{a->tempGrav[[n,1]],b->tempGrav[[n,2]],c->tempGrav[[n,3]]};
			{fp,fval,passTests} = valFP[{approxN,tempFlowEq,fp0,bc,valOptions}];
			fp0 = fp + imagPart;
			Sow[{fp,fval,passTests}];
			,{n,nsteps}]][[2,1]];
		values = MapThread[Join[#1, {#2}] &, {values, points}];
		Return[values];
		];

(* Same as varMatter2 but prints the iteration number *)
varMatter2[approxN_,flowEq_,fp1_,matterContent_,p1_,p2_,nsteps_,bc_:{0,0},valOptions_:{100,100,10^-60,0.1,2},spacing_:"linear",logPrecision_:3,imagPart_:10^-20 I] := 
	Block[{logSpacing,points,tempMat,tempFlowEq,fp0,values,fp,fval,passTests},
		If[Not[Length[flowEq] == approxN],
			Print["Flow equation does not have the correct number of elements."];
			Abort[];
			];
		If[Not[Length[fp1] == approxN],
			Print["Initial fixed point does not have the correct number of elements."];
			Abort[];
			];

		(* Straight line *)
		(* lineEq = LeastSquares[{{#1[[1]],1},{#2[[1]],1}},{#1[[2]],#2[[2]]}].{x,1}&@@{p1,p2}; *)
		(* In logSpacing, c controls the precision. 
			c = 1: lower precision at the dense point and more evenly distributed points; 
			c = high: higher precision at the dense point and less evenly distributed points; *)
		logSpacing[a_,b_,c_] := (10^Subdivide[Log10[10^-c],Log10[1],Max[nsteps-2,2]])(b-a)+a;
		Switch[spacing,
			"linear",
				(* Straight line *)
				points = Subdivide[p1,p2,Max[nsteps-1,2]];,
			"logSparseToDense",
				(* Log spacing sparse to dense *)
				(* The initial point is exact, the end point is approximate *)
				points = Reverse[Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p2[[1]],p1[[1]],logPrecision],logSpacing[p2[[2]],p1[[2]],logPrecision],logSpacing[p2[[3]],p1[[3]],logPrecision]}],p2]];,
			"logDenseToSparse",
				(* Log spacing dense to spare *)
				(* The initial point is approximate, the end point is exact *)
				points = Prepend[MapThread[{#1,#2,#3}&,{logSpacing[p1[[1]],p2[[1]],logPrecision],logSpacing[p1[[2]],p2[[2]],logPrecision],logSpacing[p1[[3]],p2[[3]],logPrecision]}],p1];
				,
			_, Print["Spacing option not valid. Aborting."];
			Abort[];
			];
		tempMat = matterContent*#&/@points;

		fp0 = fp1;
		values = Reap[Do[
			log["Iteration: "<>ToString[n]<>" of "<>ToString[nsteps]];
			tempFlowEq = flowEq/.{ns->tempMat[[n,1]],nm->tempMat[[n,2]],nd->tempMat[[n,3]]};
			{fp,fval,passTests} = valFP[{approxN,tempFlowEq,fp0,bc,valOptions}];
			fp0 = fp + imagPart;
			Sow[{fp,fval,passTests}];
			,{n,nsteps}]][[2,1]];
		values = MapThread[Join[#1, {#2}] &, {values, points}];
		Return[values];
		]

eigSystemTest[eq_,fp_,boundaryConditions_,minN_,maxN_,parallel_:False] :=
	Block[{u,v,matV,matdU,fpValues,bcValues,eig,temp1,temp2,tempValues},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		fpValues = Table[\[Lambda][k] -> fp[[m,n-minN+1,k+1]],{m,Length[fp]},{n,minN,maxN},{k,0,n-1}];
		bcValues = Table[{\[Lambda][k] -> boundaryConditions[[m,k-minN+1,1]],\[Lambda][k+1] -> boundaryConditions[[m,k-minN+1,2]]},{m,Length[fp]},{k,minN,maxN}];
		(* fpValues: {fpN, totalN, minN:maxN} *)
		(* bcValues: {fpN, totalN, 2} *)

		If[TrueQ[parallel],
			Print[AbsoluteTiming[matV = ParallelTable[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]];
			Print[AbsoluteTiming[matdU = ParallelTable[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]];
			(* matV: {maxN, maxN} *)
			(* matdU: {maxN, maxN} *)
			
			Print[AbsoluteTiming[eig = ParallelTable[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]])]
				,{m,Length[fp]},{i,minN,maxN}];]];
			(* eig: {fpN, totalN, 2} *)
			(* eig[[fpN, totalN, 1]]: {minN:maxN} *)
			(* eig[[fpN, totalN, 2]]: {minN:maxN, minN:maxN} *)
			,
			Print["Create matV: "<>ToString[AbsoluteTiming[matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
			Print["Create matdU: "<>ToString[AbsoluteTiming[matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
			Return[{matdU,matV}];
			
			Print["Make matV a function: "<>ToString[AbsoluteTiming[matV = Function@@{matV/.{\[Lambda][n_]->\[Lambda][n+1]}}/.Dispatch@Thread[Array[\[Lambda],maxN+2]->Array[Slot,maxN+2]];]]];
			Print["Make matdU a function: "<>ToString[AbsoluteTiming[matdU = Function@@{matdU/.{\[Lambda][n_]->\[Lambda][n+1]}}/.Dispatch@Thread[Array[\[Lambda],maxN+2]->Array[Slot,maxN+2]];]]];
			tempValues = MapThread[PadRight[Join[#1,#2],{Length[#1],maxN+2}]&,{fp,boundaryConditions}];
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = Map[matV@@#&,tempValues,{2}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[IdentityMatrix[i]-Take[temp1[[m,i-1]],{1,i},{1,i}],{m,Length[fp]},{i,minN,maxN}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{2}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = Map[matdU@@#&,tempValues,{2}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-1]],{1,i},{1,i}],{m,Length[fp]},{i,minN,maxN}];]]];
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[Monitor[eig = Table[-Eigensystem[(temp1[[m,i]]).(temp2[[m,i]])],{m,Length[fp]},{i,maxN-minN+1}];,{m,i}]]]];
			];
		Return[eig];
		];

eigSystemTest2[eq_,fp_,boundaryConditions_,minN_,maxN_,parallel_:False] :=
	Block[{u,v,matV,matdU,fpValues,bcValues,eig,temp1,temp2,tempValues},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)

		(* Functions for building the stability matrix *)
		u[i_,j_] := D[-(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1],\[Lambda][j][t]];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		fpValues = Table[\[Lambda][k] -> fp[[m,n-minN+1,k+1]],{m,Length[fp]},{n,minN,maxN},{k,0,n-1}];
		bcValues = Table[{\[Lambda][k] -> boundaryConditions[[m,k-minN+1,1]],\[Lambda][k+1] -> boundaryConditions[[m,k-minN+1,2]]},{m,Length[fp]},{k,minN,maxN}];
		(* fpValues: {fpN, totalN, minN:maxN} *)
		(* bcValues: {fpN, totalN, 2} *)

		If[TrueQ[parallel],
			Print[AbsoluteTiming[matV = ParallelTable[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]];
			Print[AbsoluteTiming[matdU = ParallelTable[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]];
			(* matV: {maxN, maxN} *)
			(* matdU: {maxN, maxN} *)
			
			Print[AbsoluteTiming[eig = ParallelTable[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]])]
				,{m,Length[fp]},{i,minN,maxN}];]];
			(* eig: {fpN, totalN, 2} *)
			(* eig[[fpN, totalN, 1]]: {minN:maxN} *)
			(* eig[[fpN, totalN, 2]]: {minN:maxN, minN:maxN} *)
			,
			Print["Create matV: "<>ToString[AbsoluteTiming[matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];
			Print["Create matdU: "<>ToString[AbsoluteTiming[matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];]]];

			numU = Function@@{u[i,j]/.{\[Lambda][n_][t]->\[Lambda][n+1]}}/.Dispatch@Thread[Array[\[Lambda],maxN+2]->Array[Slot,maxN+2]];

			Print["Make matV a function: "<>ToString[AbsoluteTiming[matV = Function@@{matV/.{\[Lambda][n_]->\[Lambda][n+1]}}/.Dispatch@Thread[Array[\[Lambda],maxN+2]->Array[Slot,maxN+2]];]]];
			Print["Make matdU a function: "<>ToString[AbsoluteTiming[matdU = Function@@{matdU/.{\[Lambda][n_]->\[Lambda][n+1]}}/.Dispatch@Thread[Array[\[Lambda],maxN+2]->Array[Slot,maxN+2]];]]];
			tempValues = MapThread[PadRight[Join[#1,#2],{Length[#1],maxN+2}]&,{fp,boundaryConditions}];
			Print["Evaluate matV: "<>ToString[AbsoluteTiming[temp1 = Map[matV@@#&,tempValues,{2}];]]];
			Print["Split matV into submatrices: "<>ToString[AbsoluteTiming[temp1 = Table[Identity[i]-Take[temp1[[m,i-1]],{1,i},{1,i}],{m,Length[fp]},{i,minN,maxN}];]]];
			Print["Invert matV: "<>ToString[AbsoluteTiming[temp1 = Map[Inverse,temp1,{2}];]]];
			Print["Evaluate matdU: "<>ToString[AbsoluteTiming[temp2 = Map[matdU@@#&,tempValues,{2}];]]];
			Print["Split matdU into submatrices: "<>ToString[AbsoluteTiming[temp2 = Table[Take[temp2[[m,i-1]],{1,i},{1,i}],{m,Length[fp]},{i,minN,maxN}];]]];
			(* Print["ok"]; *)
			(* Print[ToString[Dimensions/@{temp1,temp2}]]; *)
			Print["Compute eigensystem: "<>ToString[AbsoluteTiming[Monitor[eig = Table[-Eigensystem[(temp1[[m,i]]).(temp2[[m,i]])],{m,Length[fp]},{i,maxN-minN+1}];,{m,i}]]]];
			];
		Return[eig];
		];


(* Legacy functions *)
(***************************************************************************)
(***************************************************************************)
(* Same as eigSystem but using replacement rules instead of anonymous functions *)
eigSystemL[eq_,fp_,boundaryConditions_,minN_,maxN_,parallel_:False] :=
	Block[{u,v,matV,matdU,fpValues,bcValues,eig},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		fpValues = Table[\[Lambda][k] -> fp[[m,n-minN+1,k+1]],{m,Length[fp]},{n,minN,maxN},{k,0,n-1}];
		bcValues = Table[{\[Lambda][k] -> boundaryConditions[[m,k-minN+1,1]],\[Lambda][k+1] -> boundaryConditions[[m,k-minN+1,2]]},{m,Length[fp]},{k,minN,maxN}];
		(* fpValues: {fpN, totalN, minN:maxN} *)
		(* bcValues: {fpN, totalN, 2} *)

		If[TrueQ[parallel],
			matV = ParallelTable[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			matdU = ParallelTable[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			(* matV: {maxN, maxN} *)
			(* matdU: {maxN, maxN} *)
			
			eig = ParallelTable[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]])]
				,{m,Length[fp]},{i,minN,maxN}];
			(* eig: {fpN, totalN, 2} *)
			(* eig[[fpN, totalN, 1]]: {minN:maxN} *)
			(* eig[[fpN, totalN, 2]]: {minN:maxN, minN:maxN} *)
			,
			matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			eig = Table[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]])]
				,{m,Length[fp]},{i,minN,maxN}];
			];
		Return[eig];
		];

(* Same as eigSystemExtraBC but using replacement rules instead of anonymous functions *)
eigSystemLExtraBC[eq_,fp_,boundaryConditions_,minN_,maxN_,parallel_:False] :=
	Block[{u,v,matV,matdU,fpValues,bcValues,eig},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		fpValues = Table[\[Lambda][k] -> fp[[m,n-minN+1,k+1]],{m,Length[fp]},{n,minN,maxN},{k,0,n-1}];
		bcValues = Table[\[Lambda][j] -> boundaryConditions[[m,k-minN+1,j-k+1]],{m,Length[fp]},{k,minN,maxN},{j,k,k+5}];
		(* fpValues: {fpN, totalN, minN:maxN} *)
		(* bcValues: {fpN, totalN, 2} *)

		If[TrueQ[parallel],
			matV = ParallelTable[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			matdU = ParallelTable[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			(* matV: {maxN, maxN} *)
			(* matdU: {maxN, maxN} *)
			
			eig = ParallelTable[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]])]
				,{m,Length[fp]},{i,minN,maxN}];
			(* eig: {fpN, totalN, 2} *)
			(* eig[[fpN, totalN, 1]]: {minN:maxN} *)
			(* eig[[fpN, totalN, 2]]: {minN:maxN, minN:maxN} *)
			,
			matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			eig = Table[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1]]/.fpValues[[m,i-minN+1]])]
				,{m,Length[fp]},{i,minN,maxN}];
			];
		Return[eig];
		];

(* Same as eigsysVarGrav but using replacement rules instead of anonymous functions *)
eigSysVarGravOld[eq_,fp_,boundaryConditions_,gravity_,minN_,maxN_,parallel_:False] :=
	Block[{u,v,matV,matdU,fpValues,bcValues,gravValues,eig},
		(* Dimensions: *)
		(* fp: {fpN, totalN, maxN} *)
		(* boundaryConditions: {fpN, totalN, 2} *)
		(* totalN = maxN - minN + 1 *)

		(* Functions for building the stability matrix *)
		u[i_] := -(eq[[i+1]]/.D[\[Lambda][_][t],t] -> 0)/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];
		v[i_,j_] := -Coefficient[eq[[i+1]],D[\[Lambda][j][t],t],1]/Coefficient[eq[[i+1]],D[\[Lambda][i][t],t],1];

		fpValues = Table[\[Lambda][k] -> fp[[m,n-minN+1,q,k+1]],{m,Length[fp]},{n,minN,maxN},{q,Length[gravity]},{k,0,n-1}];
		bcValues = Table[{\[Lambda][k] -> boundaryConditions[[m,k-minN+1,q,1]],\[Lambda][k+1] -> boundaryConditions[[m,k-minN+1,q,2]]},{m,Length[fp]},{k,minN,maxN},{q,Length[gravity]}];
		gravValues = Table[{a->gravity[[q,1]],b->gravity[[q,2]],c->gravity[[q,3]]},{q,Length[gravity]}];
		(* fpValues: {fpN, totalN, minN:maxN} *)
		(* bcValues: {fpN, totalN, 2} *)

		If[TrueQ[parallel],
			matV = ParallelTable[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			matdU = ParallelTable[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			(* matV: {maxN, maxN} *)
			(* matdU: {maxN, maxN} *)
			
			eig = ParallelTable[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1,q]]/.fpValues[[m,i-minN+1,q]]/.gravValues[[q]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1,q]]/.fpValues[[m,i-minN+1,q]]/.gravValues[[q]])]
				,{m,Length[fp]},{i,minN,maxN},{q,Length[gravity]}];
			(* eig: {fpN, totalN, # of gravitational configurations, 2} *)
			(* eig[[fpN, totalN, # of gravitational configurations, 1]]: {minN:maxN} *)
			(* eig[[fpN, totalN, # of gravitational configurations, 2]]: {minN:maxN, minN:maxN} *)
			,
			matV = Table[v[i,j]Boole[i!=j]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			matdU = Table[D[u[i],\[Lambda][j][t]]/.n_[t] -> n,{i,0,maxN-1},{j,0,maxN-1}];
			eig = Table[-Eigensystem[Quiet[Inverse[(IdentityMatrix[i]-matV[[;;i,;;i]])/.bcValues[[m,i-minN+1,q]]/.fpValues[[m,i-minN+1,q]]/.gravValues[[q]]]]
				.(matdU[[;;i,;;i]]/.bcValues[[m,i-minN+1,q]]/.fpValues[[m,i-minN+1,q]]/.gravValues[[q]])]
				,{m,Length[fp]},{i,minN,maxN},{q,Length[gravity]}];
			];
		Return[eig];
		];

