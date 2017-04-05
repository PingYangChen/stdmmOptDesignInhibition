#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

typedef struct {
  rowvec GBest;
  double fGBest;
  rowvec fGBestHist;
  mat PBest;
	vec fPBest;		
} PSO_Result, *Ptr_PSO_Result;

void PSO_MAIN(Ptr_PSO_Result Ptr_PSO_Result, const int &LOOPINFO, const mat &OPTIONMAT, const vec &DESIGNINFO, const rowvec &FIXEDVALUE);
void initialRandSwarm(mat &swarm, mat &vStep, const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO);
void GetBoundary(mat &vmax, mat &subd, mat &slbd, const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO);
void updatePos(mat &swarm, const mat &vStep, const mat &subd, const mat &slbd, const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO);
void fEval(vec &fSwarm, const mat &swarm, const mat &OPTIONMAT, const vec &DESIGNINFO, const rowvec &FIXEDVALUE, const int &LOOPINFO);
mat inforMat(rowvec DESIGN, rowvec para, const vec &DESIGNINFO);
colvec fDev(rowvec x, rowvec p, const vec &DESIGNINFO);
rowvec Sol_locallyOptimal(rowvec &para, const vec &DESIGNINFO);
double Calc_DISP(const rowvec &x, const rowvec &DESIGN, const rowvec &ASSIST, const vec &DESIGNINFO, const double &fGBest);
void SWARM_INIT_TRICK(mat &swarm, const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO);
void ASSISTANT_SEARCH_CORNER(int &nAst, rowvec &ASSIST, vec &BestF, const rowvec &DESIGN, const double &fBest, const mat &OPTIONMAT, const vec &DESIGNINFO);

//[[Rcpp::export]]
List NESTEDPSO(const int LOOPINFO, NumericMatrix OPTIONMAT_INPUT, NumericVector DESIGNINFO_INPUT, NumericVector FIXEDVALUE_INPUT)
{
	arma_rng::set_seed_random();
	int NCPU = omp_get_max_threads();
	omp_set_num_threads(NCPU - 1);
	
	arma::mat OPTIONMAT(OPTIONMAT_INPUT.begin(), OPTIONMAT_INPUT.nrow(), OPTIONMAT_INPUT.ncol(), false);
  arma::vec DESIGNINFO(DESIGNINFO_INPUT.begin(), DESIGNINFO_INPUT.size(), false);
	arma::rowvec FIXEDVALUE(FIXEDVALUE_INPUT.begin(), FIXEDVALUE_INPUT.size(), false);
	
	PSO_Result Result = {};
	//mexPrintf("Calling Cpp PSO Kernel... \n");
	
  PSO_MAIN(&Result, LOOPINFO - 1, OPTIONMAT, DESIGNINFO, FIXEDVALUE);
		
	return List::create(Named("GBest") = wrap(Result.GBest),
											Named("fGBest") = wrap(Result.fGBest),
											Named("fGBestHist") = wrap(Result.fGBestHist),
											Named("PBest") = wrap(Result.PBest),
											Named("fPBest") = wrap(Result.fPBest));
}

//[[Rcpp::export]]
List MUSEARCH(double doublefVal, NumericVector DESIGN_INPUT, NumericMatrix OPTIONMAT_INPUT, NumericVector DESIGNINFO_INPUT)
{
	arma_rng::set_seed_random();
	int NCPU = omp_get_max_threads();
	omp_set_num_threads(NCPU - 1);

	arma::rowvec DESIGN(DESIGN_INPUT.begin(), DESIGN_INPUT.size(), false);
	arma::mat OPTIONMAT(OPTIONMAT_INPUT.begin(), OPTIONMAT_INPUT.nrow(), OPTIONMAT_INPUT.ncol(), false);
  arma::vec DESIGNINFO(DESIGNINFO_INPUT.begin(), DESIGNINFO_INPUT.size(), false);
	
	//mexPrintf("Calling Cpp ASSISTANT_SEARCH Kernel... \n");
	int nAst;
	arma::rowvec ASSIST;
	arma::vec BestF;
	
	ASSISTANT_SEARCH_CORNER(nAst, ASSIST, BestF, DESIGN, doublefVal, OPTIONMAT, DESIGNINFO);
	
	arma::rowvec ASSIST_WEIGHT;
	if (nAst > 1) {
		//mexPrintf("nAst = %d \nPSO Search Weight..\n", nAst);
		int dPara = (int)DESIGNINFO[3];
		OPTIONMAT(3, 3) = (int)(ASSIST.n_cols/dPara - 1);
		int ASTLOOP = 3;
		uword dLength = (uword)DESIGN.n_elem;
		uword aLength = (uword)ASSIST.n_elem;
		arma::rowvec FIXEDVALUE(1 + dLength + aLength);
		FIXEDVALUE(0) = doublefVal;
		FIXEDVALUE.subvec(1, dLength) = DESIGN;
		FIXEDVALUE.subvec(1 + dLength, dLength + aLength) = ASSIST;
		
		PSO_Result ASTWT_Result = {};
		PSO_MAIN(&ASTWT_Result, ASTLOOP, OPTIONMAT, DESIGNINFO, FIXEDVALUE);
		ASSIST_WEIGHT = ASTWT_Result.GBest;
	} else {
		//mexPrintf("nAst = 1 \n");
		ASSIST_WEIGHT = arma::ones(1);
	}

	return List::create(Named("Mu") = ASSIST,
											Named("Mu_w") = ASSIST_WEIGHT,
											Named("BestF") = BestF);	
}

//[[Rcpp::export]]
List QUICKDISPERSION(double fGBest, NumericVector DESIGN_INPUT, NumericVector ASSIST_INPUT, NumericVector DESIGNINFO_INPUT)
{
	int NCPU = omp_get_max_threads();
	omp_set_num_threads(NCPU - 1);
		
	arma::rowvec ASSIST(ASSIST_INPUT.begin(), ASSIST_INPUT.size(), false);
	arma::rowvec DESIGN(DESIGN_INPUT.begin(), DESIGN_INPUT.size(), false);
  arma::vec DESIGNINFO(DESIGNINFO_INPUT.begin(), DESIGNINFO_INPUT.size(), false);	

	/* REVISE HERE */
	int nGrid = 120;
	int dSupp = (int)DESIGNINFO(1);
	arma::vec ub = conv_to<vec>::from(DESIGNINFO.subvec(4, 3 + dSupp));
	arma::vec lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + dSupp, 3 + 2*dSupp));	
	arma::vec sLine = linspace<vec>(lb(0), ub(0), nGrid);
	arma::vec iLine = linspace<vec>(lb(1), ub(1), nGrid);
	
	arma::mat sMesh(nGrid, nGrid), iMesh(nGrid, nGrid), dispVal(nGrid, nGrid);
	sMesh.zeros(); iMesh.zeros(); dispVal.zeros();
	int iDraw;
	double maxVal;
	#pragma omp parallel private(iDraw) 
	{
		#pragma omp for
		for (iDraw = 0; iDraw < nGrid; iDraw++) {
			for (int j = 0; j < nGrid; j++) {
				sMesh(iDraw, j) = sLine(iDraw);
				iMesh(iDraw, j) = iLine(j);
				arma::rowvec x; 
				x << sLine(iDraw) << iLine(j) << endr;
				dispVal(iDraw, j) = Calc_DISP(x, DESIGN, ASSIST, DESIGNINFO, fGBest);
			}
		}
	}
	maxVal = max(max(dispVal));
	/* END REVISION */
	
	/* OUTPUT */	
	return List::create(Named("sGrid") = wrap(sLine),
											Named("iGrid") = wrap(iLine),
											Named("sMesh") = wrap(sMesh),
											Named("iMesh") = wrap(iMesh),
											Named("dispVal") = wrap(dispVal),
											Named("maxVal") =wrap(maxVal));	
}

// PSO options
// ( 0) 1: Outer Loop; 2: Inner Loop; 3: Local Loop
// ( 1) 1: Maximize; 0: Minimize
// ( 2) 1: Terminate by MaxIter; 0: Terminate by Conv
// ( 3) swarm dimension
// ( 4) swarm size
// ( 5) maximum iteration
// ( 6) Tolerence for Terminate by Conv
// ( 7) c1
// ( 8) c2
// ( 9) Vmax coefficient
// (10) w start
// (11) w end
// (12) w var
// (13) Chi
// (14) Number of neighborhood
void PSO_MAIN(Ptr_PSO_Result Ptr_PSO_Result, const int &LOOPINFO, const mat &OPTIONMAT, 
							const vec &DESIGNINFO, const rowvec &FIXEDVALUE) 
{
	// SET PSO PARAMETERS
	arma::rowvec OPTIONS = conv_to<rowvec>::from(OPTIONMAT.row(LOOPINFO));
	//mexPrintf("IN PSO LOOP: %d \n", LOOPINFO);
	int dSwarm = (int)OPTIONS(3);
	int nSwarm = (int)OPTIONS(4);
	int maxIter		= (int)OPTIONS(5);
	int maximize 	= (int)OPTIONS(1);
	double c1 = OPTIONS(7);
	double c2 = OPTIONS(8);	
	double w_end       = OPTIONS(11);         // Final inertia weight value
	double w_start     = OPTIONS(10);       // Initial inertia weight value
	int w_varyfor   = (int)(OPTIONS(12)*maxIter);
	double w_now       = w_start;
	double inertdec    = (w_start - w_end)/w_varyfor;      // Inertia weight change per iteration step

	arma::mat swarm(nSwarm, dSwarm), vStep(nSwarm, dSwarm), PBest(nSwarm, dSwarm);
	arma::mat subd(nSwarm, dSwarm), slbd(nSwarm, dSwarm), vmax(nSwarm, dSwarm), GBmat(nSwarm, dSwarm);
	arma::rowvec GBest(dSwarm);
	arma::vec fSwarm(nSwarm), fPBest(nSwarm);
	double fGBest;
	arma::uword GBestIdx;
	arma::rowvec fGBestHist(maxIter);
	
	////// INITIALIZE
	//if (LOOPINFO == 0) { mexPrintf("Outer Loop: Initialize .. "); mexEvalString("drawnow;"); }
	
	// GENERATE THE VMAX AND SWARM BOUNDARY MATRIX
	GetBoundary(vmax, subd, slbd, OPTIONS, DESIGNINFO, LOOPINFO); 
	// INITIALIZE RANDOM SWARM
	initialRandSwarm(swarm, vStep, OPTIONS, DESIGNINFO, LOOPINFO);
	SWARM_INIT_TRICK(swarm, OPTIONS, DESIGNINFO, LOOPINFO);	
	// --- Initialization of function evaluation
	fEval(fSwarm, swarm, OPTIONMAT, DESIGNINFO, FIXEDVALUE, LOOPINFO); 
	// Initialize LOCAL BEST
	fPBest = fSwarm;	PBest = swarm;
	// Initialize GLOBAL BEST
	if (maximize == 0) {
		fGBest = fPBest.min(GBestIdx);	
		GBest = PBest.row(GBestIdx);		
	} else {
		fGBest = fPBest.max(GBestIdx);	
		GBest = PBest.row(GBestIdx);	
	}
	//if (LOOPINFO == 0) mexPrintf("DONE \n"); 
	////// PSO LOOP
	int t; // SET ITERATION COUNTER
	for (t = 0; t < maxIter; t++) {
		/*if (LOOPINFO == 0) {
			if (t == 0) { mexPrintf("Outer Loop: Updating ..    "); }
			mexPrintf("\b\b\b%2.0f%%", (double)((t+1)*100/maxIter)); mexEvalString("drawnow;");
			if (t == (maxIter - 1)) { mexPrintf("\n"); }
		}*/
		// UPDATING THE INERTIA WEIGHT
		if (t <= w_varyfor) { w_now = w_now - inertdec; }
		// UPDATING VELOCITY
		GBmat = repmat(GBest, nSwarm, 1); 
		vStep = w_now*vStep + c1*randu(nSwarm, dSwarm) % (PBest - swarm) + c2*randu(nSwarm, dSwarm) % (GBmat - swarm);
		vStep = min(vStep, vmax);
		vStep = max(vStep, (-1)*vmax);
		// UPDATING POSITION
		updatePos(swarm, vStep, subd, slbd, OPTIONS, DESIGNINFO, LOOPINFO);	
		// UPDATING OBJECTIVE FUNCTION VALUES
		fEval(fSwarm, swarm, OPTIONMAT, DESIGNINFO, FIXEDVALUE, LOOPINFO); 
		// UPDATING THE LOCAL BEST
		// UPDATING THE GLOBAL BEST
		if (maximize == 0) {
			if (any(fSwarm < fPBest)) {
				uvec RowChange = find(fSwarm < fPBest);
				fPBest.elem(RowChange) = fSwarm.elem(RowChange);
				PBest.rows(RowChange) = swarm.rows(RowChange);
			}
			if (min(fPBest) < fGBest) {
				fGBest = fPBest.min(GBestIdx);	
				GBest = PBest.row(GBestIdx);	
			}		
		} else {
			if (any(fSwarm > fPBest)) {
				uvec RowChange = find(fSwarm > fPBest);
				fPBest.elem(RowChange) = fSwarm.elem(RowChange);
				PBest.rows(RowChange) = swarm.rows(RowChange);
			}
			if (max(fPBest) > fGBest) {
				fGBest = fPBest.max(GBestIdx);	
				GBest = PBest.row(GBestIdx);	
			}		
		}
		// RECORDING THE CURRENT GLOBAL BEST VALUE
		fGBestHist(t) = fGBest; 
	}
	// OUTPUT
	Ptr_PSO_Result->GBest = GBest;
	Ptr_PSO_Result->fGBest = fGBest;
	Ptr_PSO_Result->fGBestHist = fGBestHist;
	Ptr_PSO_Result->PBest = PBest;
	Ptr_PSO_Result->fPBest = fPBest;
}

void initialRandSwarm(mat &swarm, mat &vStep, const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO) 
{
	int dSwarm = (int)OPTIONS(3);
	int nSwarm = (int)OPTIONS(4);
	
	int dSupp = (int)DESIGNINFO(1);
	int nSupp = (int)DESIGNINFO(2);
	int dPara = (int)DESIGNINFO(3);
	//mexPrintf("INIT  LOOP: %d\n", LOOPINFO);
	vec ub, lb;
	mat bddTmp;
	// INITIALIZE POSITION
	switch (LOOPINFO) {
		// OUTER LOOP: DESIGNS
		case 0: 
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4, 3 + dSupp));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + dSupp, 3 + 2*dSupp));
			// SUPPORT POINTS
			for (int i = 0; i < dSupp; i++) {
				bddTmp = (ub(i)-lb(i))*ones(nSwarm, nSupp);
				swarm.cols(nSupp*i, nSupp*(i+1) - 1) = randu(nSwarm, nSupp) % bddTmp + lb(i);
			}
			// WEIGHTS
			swarm.col(dSupp*nSupp) = randu(nSwarm, 1)*0.5;
			for (int i = 1; i < (nSupp - 1); i++) {
				swarm.col(dSupp*nSupp + i) = randu(nSwarm, 1) % (1 - sum(swarm.cols(dSupp*nSupp, dSupp*nSupp + i - 1), 1));
			}
		break;
		// INNER LOOP: MODEL PARAMETERS
		case 1:
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp, 3 + 2*dSupp + dPara));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp + dPara, 3 + 2*dSupp + 2*dPara));
			swarm = randu(nSwarm, dPara) % repmat((ub.t() - lb.t()), nSwarm, 1) + repmat(lb.t(), nSwarm, 1);
		break;
		// LOCAL LOOP: DESIGNS
		case 2: 
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4, 3 + dSupp));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + dSupp, 3 + 2*dSupp));
			// SUPPORT POINTS
			for (int i = 0; i < dSupp; i++) {
				bddTmp = (ub(i)-lb(i))*ones(nSwarm, nSupp);
				swarm.cols(nSupp*i, nSupp*(i+1) - 1) = randu(nSwarm, nSupp) % bddTmp + lb(i);
			}
			// WEIGHTS
			swarm.col(dSupp*nSupp) = randu(nSwarm, 1)*0.5;
			for (int i = 1; i < (nSupp - 1); i++) {
				swarm.col(dSupp*nSupp + i) = randu(nSwarm, 1) % (1 - sum(swarm.cols(dSupp*nSupp, dSupp*nSupp + i - 1), 1));
			}
		break;
		// ASSISTANT WEIGHT LOOP: WEIGHTS OF THE ASSISTANT DESIGN
		case 3:
			swarm.col(0) = randu(nSwarm, 1)*0.5;
			for (int i = 1; i < dSwarm; i++) {
				swarm.col(i) = randu(nSwarm, 1) % (1 - sum(swarm.cols(0, i - 1), 1));
			}
		break;
	}
	// INITIALIZE VELOCITY
	vStep = randu(nSwarm, dSwarm);
}

void GetBoundary(mat &vmax, mat &subd, mat &slbd, const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO) 
{
	int dSwarm = (int)OPTIONS(3);
	//int nSwarm = (int)OPTIONS(4);
	int dSupp = (int)DESIGNINFO(1);
	int nSupp = (int)DESIGNINFO(2);
	int dPara = (int)DESIGNINFO(3);
	double vk = OPTIONS(9);
	//mexPrintf("GET VMAX  LOOP: %d\n", LOOPINFO);
	vec ub, lb, vmaxVal;
	// Fixing velocity if it absolute value over Vmax
	subd.zeros();
	slbd.zeros();
	switch (LOOPINFO) {
		// OUTER LOOP: DESIGNS
		case 0: 
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4, 3 + dSupp));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + dSupp, 3 + 2*dSupp));
			vmaxVal = (ub - lb)/vk;
			// VMAX AND SWARM BOUNDARY FOR SUPPORT POINTS
			for (int i = 0; i < dSupp; i++) {
				subd.cols(nSupp*i, nSupp*(i+1) - 1).fill(ub(i));
				slbd.cols(nSupp*i, nSupp*(i+1) - 1).fill(lb(i));
				vmax.cols(nSupp*i, nSupp*(i+1) - 1).fill(vmaxVal(i));
			}
			// VMAX FOR WEIGHTS
			// BOUNDARY FOR WEIGHTS SHELL NOT HAVE TO BE ASSIGNED, THE WEIGHT WILL BE UPDATED CELL BY CELL
			vmax.cols(nSupp*dSupp, dSwarm - 1).fill(1.0/vk);
		break;
		// INNER LOOP: MODEL PARAMETERS
		case 1:
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp, 3 + 2*dSupp + dPara));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp + dPara, 3 + 2*dSupp + 2*dPara));
			vmaxVal = (ub - lb)/vk;
			for (int i = 0; i < dPara; i++) {
				subd.col(i).fill(ub(i));
				slbd.col(i).fill(lb(i));
				vmax.col(i).fill(vmaxVal(i));
			}
		break;
		// LOCAL LOOP: DESIGNS
		case 2: 
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4, 3 + dSupp));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + dSupp, 3 + 2*dSupp));
			vmaxVal = (ub - lb)/vk;
			// VMAX AND SWARM BOUNDARY FOR SUPPORT POINTS
			for (int i = 0; i < dSupp; i++) {
				subd.cols(nSupp*i, nSupp*(i+1) - 1).fill(ub(i));
				slbd.cols(nSupp*i, nSupp*(i+1) - 1).fill(lb(i));
				vmax.cols(nSupp*i, nSupp*(i+1) - 1).fill(vmaxVal(i));
			}
			// VMAX FOR WEIGHTS
			// BOUNDARY FOR WEIGHTS SHELL NOT HAVE TO BE ASSIGNED, THE WEIGHT WILL BE UPDATED CELL BY CELL
			vmax.cols(nSupp*dSupp, dSwarm - 1).fill(1.0/vk);
		break;
		// ASSISTANT WEIGHT LOOP: WEIGHTS OF THE ASSISTANT DESIGN
		case 3:
			vmax.fill(1.0/vk);
		break;
	}
}

void updatePos(mat &swarm, const mat &vStep, const mat &subd, const mat &slbd, \
								const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO) 
{
	int dSwarm = (int)OPTIONS(3);
	double chi = OPTIONS(13);
	int dSupp = (int)DESIGNINFO(1);
	int nSupp = (int)DESIGNINFO(2);
	//mexPrintf("P LOOP: %d\n", LOOPINFO);
	swarm += chi*vStep;
	
	vec ub, lb;
	mat swarmTmp1;
	vec swarmTmp2;
	vec WEIGHTSUM;
	uvec RowChange;
	// UPDATE POSITION
	switch (LOOPINFO) {
		// OUTER LOOP: DESIGNS
		case 0: 
			// SUPPORT POINTS
			swarmTmp1 = swarm.cols(0, dSupp*nSupp - 1);
			swarmTmp1 = min(swarmTmp1, subd.cols(0, dSupp*nSupp - 1));
			swarmTmp1 = max(swarmTmp1, slbd.cols(0, dSupp*nSupp - 1));
			swarm.cols(0, dSupp*nSupp - 1) = swarmTmp1;
			// WEIGHTS
			swarmTmp2 = swarm.col(dSupp*nSupp);
			RowChange = find(swarmTmp2 > 1 || swarmTmp2 < 0);
			swarmTmp2.elem(RowChange) = randu(RowChange.n_rows, 1)*0.5;	
			swarm.col(dSupp*nSupp) = swarmTmp2;
			for (int i = 1; i < (nSupp - 1); i++) {
				WEIGHTSUM = sum(swarm.cols(dSupp*nSupp, dSupp*nSupp + i - 1), 1);
				swarmTmp2 = swarm.col(dSupp*nSupp + i);
				RowChange = find(swarmTmp2 > (1 - WEIGHTSUM) || swarmTmp2 < 0);
				swarmTmp2.elem(RowChange) = randu(RowChange.n_rows, 1) % (1 - WEIGHTSUM(RowChange));	
				swarm.col(dSupp*nSupp + i) = swarmTmp2;
			}
		break;
		// INNER LOOP: MODEL PARAMETERS
		case 1:
			swarm = min(swarm, subd);
			swarm = max(swarm, slbd);
		break;
		// LOCAL LOOP: DESIGNS
		case 2: 
			// SUPPORT POINTS
			swarmTmp1 = swarm.cols(0, dSupp*nSupp - 1);
			swarmTmp1 = min(swarmTmp1, subd.cols(0, dSupp*nSupp - 1));
			swarmTmp1 = max(swarmTmp1, slbd.cols(0, dSupp*nSupp - 1));
			swarm.cols(0, dSupp*nSupp - 1) = swarmTmp1;
			// WEIGHTS
			swarmTmp2 = swarm.col(dSupp*nSupp);
			RowChange = find(swarmTmp2 > 1 || swarmTmp2 < 0);
			swarmTmp2.elem(RowChange) = randu(RowChange.n_rows, 1)*0.5;	
			swarm.col(dSupp*nSupp) = swarmTmp2;
			for (int i = 1; i < (nSupp - 1); i++) {
				WEIGHTSUM = sum(swarm.cols(dSupp*nSupp, dSupp*nSupp + i - 1), 1);
				swarmTmp2 = swarm.col(dSupp*nSupp + i);
				RowChange = find(swarmTmp2 > (1 - WEIGHTSUM) || swarmTmp2 < 0);
				swarmTmp2.elem(RowChange) = randu(RowChange.n_rows, 1) % (1 - WEIGHTSUM(RowChange));	
				swarm.col(dSupp*nSupp + i) = swarmTmp2;
			}
		break;
		// ASSISTANT WEIGHT LOOP: WEIGHTS OF THE ASSISTANT DESIGN
		case 3:
			swarmTmp2 = swarm.col(0);
			RowChange = find(swarmTmp2 > 1 || swarmTmp2 < 0);
			swarmTmp2.elem(RowChange) = randu(RowChange.n_rows, 1)*0.5;	
			swarm.col(0) = swarmTmp2;
			for (int i = 1; i < dSwarm; i++) {
				WEIGHTSUM = sum(swarm.cols(0, i - 1), 1);
				swarmTmp2 = swarm.col(i);
				RowChange = find(swarmTmp2 > (1 - WEIGHTSUM) || swarmTmp2 < 0);
				swarmTmp2.elem(RowChange) = randu(RowChange.n_rows, 1) % (1 - WEIGHTSUM(RowChange));	
				swarm.col(i) = swarmTmp2;
			}
		break;
	}
}

void ASSISTANT_SEARCH_CORNER(int &nAst, rowvec &ASSIST, vec &BestF, const rowvec &DESIGN, const double &fBest, const mat &OPTIONMAT, const vec &DESIGNINFO) 
{
	int dSupp = (int)DESIGNINFO(1);
	int dPara = (int)DESIGNINFO(3);
	arma::vec pub = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp, 3 + 2*dSupp + dPara));
	arma::vec plb = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp + dPara, 3 + 2*dSupp + 2*dPara));
	arma::uvec pFixed = find((pub - plb) == 0); 
	arma::uvec pSearch = find((pub - plb) != 0);
	
	/* (Abandoned) Modify the search criterion to fit the def. of $\mathcal{N}(\xi^*)$ 
	   in Dette, H., Haines, L. M., and Imhof, L. A. (2007). */
	int sdPara = (int)pSearch.n_rows;
	int nGrid = std::pow(2, sdPara);
	arma::mat s = zeros(nGrid, dPara); 
	arma::mat d = zeros(nGrid, dPara); 
	
	s.cols(pFixed) = repmat(pub.elem(pFixed), 1, nGrid).t();
	d.cols(pFixed).zeros();

	for (int i = 0; i < sdPara; i++) {
		s.col(pSearch(i)) = repmat( \
			join_vert(ones(std::pow(2, (sdPara-1-i)))*plb(pSearch(i)), 
								ones(std::pow(2, (sdPara-1-i)))*pub(pSearch(i))
			), std::pow(2, i), 1);
		d.col(pSearch(i)) = repmat( \
			join_vert(ones(std::pow(2, (sdPara-1-i)))*(pub(pSearch(i))-plb(pSearch(i)))/99, 
								ones(std::pow(2, (sdPara-1-i)))*(plb(pSearch(i))-pub(pSearch(i)))/99
			), std::pow(2, i), 1);
	}
	arma::mat eyeMat(dPara, dPara); eyeMat.eye();
	arma::vec fvals(nGrid), fmins(nGrid); fmins.zeros();
	arma::vec nbVals(dPara); 
	arma::rowvec FIXEDVALUE;
	int F_LOOPINDEX = 1;
	
	int iAstPt; int flag = 0;
	while ((all(fmins == 0)) & (flag < 100)) {
		flag++;
		#pragma omp parallel private(iAstPt) 
		{
			#pragma omp for
			for (iAstPt = 0; iAstPt < nGrid; iAstPt++) {
				arma::vec MIN_EFF(1), COMP_EFF(sdPara);
				arma::mat paraTmp1 = conv_to<mat>::from(s.row(iAstPt));
				fEval(MIN_EFF, paraTmp1, OPTIONMAT, DESIGNINFO, DESIGN, F_LOOPINDEX);  
				arma::mat paraTmp2 = repmat(s.row(iAstPt), sdPara, 1) + (repmat(d.row(iAstPt), sdPara, 1) % eyeMat.rows(pSearch));
				fEval(COMP_EFF, paraTmp2, OPTIONMAT, DESIGNINFO, DESIGN, F_LOOPINDEX);	
				//if (all(COMP_EFF >= MIN_EFF(0))) { fmins(iAstPt) = 1; }
				if (std::abs(fBest - MIN_EFF(0)) < 1e-5) { fmins(iAstPt) = 1; }
				fvals(iAstPt) = (double)MIN_EFF(0);
			}
		}
	}
	arma::mat BestS;
	if (all(fmins == 0)) {
		//mexPrintf("Cannot find the assistant design. \nReturn the parameter point with highest min(eff.) \n");
		BestF(0) = fvals.min();
		ASSIST = vectorise(s.rows(find(fvals == BestF(0)))).t();
		nAst = ASSIST.n_cols/dPara;
	} else {
		BestS = s.rows(find(fmins == 1));
		BestF = fvals.elem(find(fmins == 1));
		nAst = (int)BestF.n_rows;
		ASSIST = vectorise(BestS).t();
	}
}

void fEval(vec &fSwarm, const mat &swarm, const mat &OPTIONMAT, const vec &DESIGNINFO, \
						const rowvec &FIXEDVALUE, const int &LOOPINFO) 
{	
	// SET PSO PARAMETERS
	int nSwarm = (int)swarm.n_rows;
	int dSupp = (int)DESIGNINFO[1];
	int nSupp = (int)DESIGNINFO[2];
	int dPara = (int)DESIGNINFO[3];
	arma::mat m; 
	int F_LOOPINDEX;
	
	switch (LOOPINFO) {
		case 0:
			// OUTER LOOP
			F_LOOPINDEX = 1;
			int iOut;
			#pragma omp parallel private(iOut) 
			{
				#pragma omp for
				for (iOut = 0; iOut < nSwarm; iOut++) {
					PSO_Result InnerResult = {};
					rowvec F_DESIGN = conv_to<rowvec>::from(swarm.row(iOut));
					PSO_MAIN(&InnerResult, F_LOOPINDEX, OPTIONMAT, DESIGNINFO, F_DESIGN);
					if (InnerResult.fGBest > 1e120) { InnerResult.fGBest = -1*datum::inf; }
					fSwarm(iOut) = InnerResult.fGBest;
				}
			}
		break;
		
		case 1:
			// INNER LOOP: MINIMUM EFFICIENCY
			F_LOOPINDEX = 2;	
			for (int i = 0; i < nSwarm; i++) {
				rowvec F_PARA = conv_to<rowvec>::from(swarm.row(i));
				rowvec L_DESIGN;
				mat L_m; 
				double L_VAL;
				if (nSupp == 3 + (int)(DESIGNINFO(0)/4)) {
					L_DESIGN = Sol_locallyOptimal(F_PARA, DESIGNINFO);
					L_m = inforMat(L_DESIGN, F_PARA, DESIGNINFO);
					L_VAL = std::log(det(L_m)); // Maximin 
				} else {
					PSO_Result LocalResult = {};
					PSO_MAIN(&LocalResult, F_LOOPINDEX, OPTIONMAT, DESIGNINFO, F_PARA);
					L_VAL = LocalResult.fGBest;
				}
				m = inforMat(FIXEDVALUE, F_PARA, DESIGNINFO);
				if (rcond(m) < (dPara*datum::eps)) {
					fSwarm(i) = datum::inf;
				} else {
					fSwarm(i) = std::log(det(m)) - L_VAL;// + std::log(100.00); // Maximin 
				}
			}
		break;
		
		case 2:
			// LOCALLY OPTIMAL DESIGN
			for (int i = 0; i < nSwarm; i++) {
				rowvec F_DESIGN = conv_to<rowvec>::from(swarm.row(i));
				m = inforMat(F_DESIGN, FIXEDVALUE, DESIGNINFO);
				if (rcond(m) < (dPara*datum::eps)) {
					fSwarm[i] = (-1)*datum::inf;
				} else {
					fSwarm[i] = std::log(det(m));
				}
			}
		break;
		
		case 3:
			// WEIGHTS OF THE ASSISTANT DESIGN
			double fGBest = (double)FIXEDVALUE(0);
			rowvec DESIGN = conv_to<rowvec>::from(FIXEDVALUE.subvec(1, dSupp*nSupp + nSupp - 1));
			rowvec AstPt = conv_to<rowvec>::from(FIXEDVALUE.subvec(dSupp*nSupp + nSupp, FIXEDVALUE.n_cols - 1));
			rowvec ASSIST;
			double vTmp;
			rowvec x(dSupp);
			for (int i = 0; i < nSwarm; i++) {
				fSwarm(i) = 0;
				ASSIST = join_horiz(AstPt, swarm.row(i));
				for (int j = 0; j < nSupp; j++) {
					for (int k = 0; k < dSupp; k++) {
						x(k) = DESIGN(j + k*nSupp);
					}
					vTmp = Calc_DISP(x, DESIGN, ASSIST, DESIGNINFO, fGBest);
					fSwarm(i) += vTmp*vTmp + (double)(1.0/nSupp);
				}
			}
		break;
	}
}

mat inforMat(rowvec DESIGN, rowvec para, const vec &DESIGNINFO) 
{
	int dSupp = (int)DESIGNINFO[1];
	int nSupp = (int)DESIGNINFO[2];
	int dPara = (int)DESIGNINFO[3];
	vec w = zeros(nSupp);
	double wTmp = 0.0;
	for (int i = 0; i < nSupp; i++) {
		if (i < (nSupp - 1)) {
			w(i) = DESIGN(dSupp*nSupp + i);
			wTmp += w(i);
		} else {
			w(i) = 1 - wTmp;
		}
	}
	// CALCULATE THE FISHER'S INFORMATION MATRIX
	mat m = zeros(dPara, dPara);
	vec d;
	rowvec x(dSupp);
	for (int i = 0; i < nSupp; i++) {
		for (int j = 0; j < dSupp; j++) {
			x(j) = DESIGN(i + j*nSupp);
		}
		d = fDev(x, para, DESIGNINFO);
		m += w(i) * d * d.t();
	}
	return(m);
}

colvec fDev(rowvec x, rowvec p, const vec &DESIGNINFO) 
{
	int MODEL_ID = (int)DESIGNINFO(0);
	colvec d(p.n_cols);
	switch (MODEL_ID) {
		// COMPETITIVE MODEL
		case 1:
			d(0) = x(0)/(p(1)*(1+x(1)/p(2))+x(0));
			d(1) = (-1)*p(0)*x(0)*(1+x(1)/p(2))/((p(1)*(1+x(1)/p(2))+x(0))*(p(1)*(1+x(1)/p(2))+x(0)));
			d(2) = p(0)*x(0)*p(1)*x(1)/(p(2)*p(2))/((p(1)*(1+x(1)/p(2))+x(0))*(p(1)*(1+x(1)/p(2))+x(0)));					
		break;					
		// NONCOMPETITIVE MODEL						
		case 2:
			d(0) = x(0)/((p(1)+x(0))*(1+x(1)/p(2)));
			d(1) = (-1)*p(0)*x(0)/(1+x(1)/p(2))/((p(1)+x(0))*(p(1)+x(0)));
			d(2) = p(0)*x(0)*x(1)/(p(1)+x(0))/((1+x(1)/p(2))*(1+x(1)/p(2)))/(p(2)*p(2));
		break;					
		// UNCOMPETITIVE MODEL
		case 3:
			d(0) = x(0)/(p(1)+(1+x(1)/p(1))*x(0));
			d(1) = (-1)*p(0)*x(0)/((p(1)+x(0)*(1+x(1)/p(2)))*(p(1)+x(0)*(1+x(1)/p(2))));
			d(2) = p(0)*x(0)*x(0)*x(1)/(p(2)*p(2))/((p(1)+x(0)+x(0)*x(1)/p(2))*(p(1)+x(0)+x(0)*x(1)/p(2)));
		break;		
		// MIXED MODEL
		case 4:
			d(0) = x(0)/(p(1)*(1+x(1)/p(2))+x(0)*(1+x(1)/p(3)));
			d(1) = -p(0)*x(0)*(1+x(1)/p(2))/((p(1)*(1+x(1)/p(2))+x(0)*(1+x(1)/p(3)))*(p(1)*(1+x(1)/p(2))+x(0)*(1+x(1)/p(3))));
			d(2) = p(0)*x(0)*p(1)*x(1)/((p(1)*(1+x(1)/p(2))+x(0)*(1+x(1)/p(3)))*(p(1)*(1+x(1)/p(2))+x(0)*(1+x(1)/p(3))))/(p(2)*p(2));
			d(3) = p(0)*x(0)*x(0)*x(1)/((p(1)*(1+x(1)/p(2))+x(0)*(1+x(1)/p(3)))*(p(1)*(1+x(1)/p(2))+x(0)*(1+x(1)/p(3))))/(p(3)*p(3));
		break;
	}
	return(d);
}

double Calc_DISP(const rowvec &x, const rowvec &DESIGN, const rowvec &ASSIST, 
									const vec &DESIGNINFO, const double &fGBest) 
{
	uword dPara = (uword)DESIGNINFO(3);
	double DELTA = (double)dPara;
	int nAst;
	if (ASSIST.n_cols == dPara + 1) {
		nAst = 1;
	}	else {
		int a = (int)ASSIST.n_cols;
		nAst = (int)((a + 1)/(dPara + 1)); 
	}

	vec astWt = zeros(nAst);
	if (nAst == 1) {
		astWt(0) = 1.0;
	} else {
		double wTmp = 0;
		for (int i = 0; i < nAst; i++) {
			if (i < (nAst - 1)) {
				astWt(i) = ASSIST(dPara*nAst + i);
				wTmp += ASSIST(dPara*nAst + i);
			} else {
				astWt(i) = 1.0 - wTmp;
			}
		}	
	}
	
	double DISPVAL = 0.0;
	vec DISPVALTMP(nAst); DISPVALTMP.zeros();
	rowvec para(dPara); para.zeros();
	mat m;
	colvec d;
	for (int i = 0; i < nAst; i++) {
		for (uword j = 0; j < dPara; j++) {
			para(j) = ASSIST(i + j*nAst);
		}
		m = inforMat(DESIGN, para, DESIGNINFO);
		if (rcond(m) < (dPara*datum::eps)) {
			DISPVALTMP(i) = 1e8;
			//mexPrintf("WARNING: SINGULAR INFORMATION MATRIX. \n");
		} else {
			d = fDev(x, para, DESIGNINFO);
			DISPVALTMP(i) = as_scalar((d.t()*m.i()*d))*astWt(i);
		}
	}
	DISPVAL = accu(DISPVALTMP);
	DISPVAL = DISPVAL - DELTA;	
	return(DISPVAL);
}

rowvec Sol_locallyOptimal(rowvec &para, const vec &DESIGNINFO) 
{
	int MODEL_ID = (int)DESIGNINFO(0);
	double ubs = DESIGNINFO(4); 
	double ubi = DESIGNINFO(5);
	double lbs = DESIGNINFO(6); 
	double lbi = DESIGNINFO(7);
	double km, kic, kiu; //v = para(0);
	rowvec S, I, P;
	
	rowvec DESIGN;
	switch (MODEL_ID) {
		// COMPETITIVE MODEL
		case 1:
			km = para(1); kic = para(2);
			S << ubs << \
					 std::max(lbs, ubs*km*(kic + lbi)/(2*km*kic + 2*km*lbi + ubs*kic)) << \
					 std::max(lbs, std::min(km*(kic + ubi)/kic, ubs)) << endr;
			I << lbi << \
					 lbi << \
					 std::min((2*km*lbi + ubs*kic + km*kic)/km, ubi) << endr;
			P << 1.0/3.0 << 1.0/3.0 << endr;
		break;					
		// NONCOMPETITIVE MODEL						
		case 2:
			km = para(1); kic = para(2);
			S << ubs << \
					 std::max(lbs, ubs*km/(ubs + 2*km)) << \
					 ubs << endr;
			I << lbi << \
					 lbi << \
					 std::min(kic + 2*lbi, ubi) << endr;
			P << 1.0/3.0 << 1.0/3.0 << endr;
		break;					
		// UNCOMPETITIVE MODEL
		case 3:
			km = para(1); kiu = para(2);
			S << ubs << \
					 std::max(lbs, ubs*km*kiu/(ubs*lbi + ubs*kiu + 2*km*kiu)) << \
					 ubs << endr;
			I << lbi << \
					 lbi << \
					 std::min((2*ubs*lbi + ubs*kiu + km*kiu)/ubs, ubi) << endr;
			P << 1.0/3.0 << 1.0/3.0 << endr;
		break;		
		// MIXED MODEL
		case 4:
			km = para(1); kic = para(2); kiu = para(3);
			I << lbi << \
					 lbi << \
					 std::min((2*ubs*kic*lbi + ubs*kic*kiu + 2*km*kiu*lbi + km*kiu*kic)/(km*kiu + ubs*kic), ubi) << \
					 std::min(lbi + sqrtf((kic + lbi)*(km*kiu*kic + km*kiu*lbi + ubs*kic*kiu + ubs*kic*lbi)/(km*kiu + ubs*kic)), ubi) << endr;
			S << ubs << \
					 std::max(lbs, ubs*km*kiu*(kic + lbi)/(ubs*kic*lbi + ubs*kic*kiu + 2*km*kiu*lbi + 2*km*kiu*kic)) << \
					 ubs << \
					 std::max((-1)*km*kiu*(kic + 2*lbi - I(3))/(kic*(kiu + 2*lbi - I(3))), lbs) << endr;
			P << 1.0/4.0 << 1.0/4.0 << 1.0/4.0 << endr;
		break;
	}
	DESIGN = arma::join_horiz(S, arma::join_horiz(I, P)); 
	return(DESIGN);
}

void SWARM_INIT_TRICK(mat &swarm, const rowvec &OPTIONS, const vec &DESIGNINFO, const int &LOOPINFO) 
{
	int MODEL_ID = (int)DESIGNINFO(0);
	int dSwarm = (int)OPTIONS(3);
	int nSwarm = (int)OPTIONS(4);
	int dSupp = (int)DESIGNINFO(1);
	int nSupp = (int)DESIGNINFO(2);
	int dPara = (int)DESIGNINFO(3);
	arma::vec ub, lb;
	arma::rowvec A;
	double ini_w;
	// VARIABLES USED IN THE INNER LOOP
	arma::uvec pFixed, pSearch;
	int sdPara, nGrid;
	arma::mat s;	
	// INITIALIZE POSITION
	switch (LOOPINFO) {
		// OUTER LOOP: DESIGNS
		case 0: 
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4, 3 + dSupp));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + dSupp, 3 + 2*dSupp));
			ini_w = (double)(1.0/nSupp);			
			// SUPPORT POINTS
			if (MODEL_ID < 4) {
				if (nSupp == 3) {
					swarm.col(0).fill(ub(0));
					swarm.col(nSupp).fill(lb(1));
					swarm.col(nSupp + 1).fill(lb(1));
					swarm.cols(dSupp*nSupp, dSwarm - 1).fill(ini_w);
				}
				if (nSupp == 4) {
					swarm.cols(0,1).fill(lb(0));
					swarm.cols(2,3).fill(ub(0));
					swarm.col(nSupp + 1).fill(lb(1));
					swarm.col(nSupp + 3).fill(lb(1));
					swarm.col(nSupp + 2).fill(ub(1));		
				}
			} else {
				if (nSupp == 4) {
					swarm.col(0).fill(ub(0));
					swarm.col(2).fill(ub(0));
					swarm.col(nSupp).fill(lb(1));
					swarm.col(nSupp + 1).fill(lb(1));	
					swarm.cols(dSupp*nSupp, dSwarm - 1).fill(ini_w);
				}
			}
			// WEIGHTS: EQUAL WEIGHTS
			swarm.submat(0, dSupp*nSupp, 0, dSwarm - 1).fill(ini_w);
		break;
		// INNER LOOP: MODEL PARAMETERS
		case 1:
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp, 3 + 2*dSupp + dPara));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + 2*dSupp + dPara, 3 + 2*dSupp + 2*dPara));
			pFixed = find((ub - lb) == 0); 
			pSearch = find((ub - lb) != 0);
			sdPara = (int)pSearch.n_rows;
			nGrid = std::pow(2, sdPara);
			s = zeros(nGrid, dPara);
			if (nSwarm >= nGrid) {
				s.cols(pFixed) = repmat(ub.elem(pFixed), 1, nGrid).t();
				for (int i = 0; i < sdPara; i++) {
					s.col(pSearch(i)) = repmat( \
						arma::join_vert(ones(std::pow(2, (sdPara-1-i)))*lb(pSearch(i)), 
														ones(std::pow(2, (sdPara-1-i)))*ub(pSearch(i))
						), std::pow(2, i), 1);
				}			
				swarm.submat(0, 0, (nGrid-1), (dPara-1)) = s;
			}
		break;
		// LOCAL LOOP: DESIGNS
		case 2: 
			ub = conv_to<vec>::from(DESIGNINFO.subvec(4, 3 + dSupp));
			lb = conv_to<vec>::from(DESIGNINFO.subvec(4 + dSupp, 3 + 2*dSupp));
			ini_w = (double)(1.0/nSupp);
			// SUPPORT POINTS
			if (MODEL_ID < 4) {
				if (nSupp == 3) {
					swarm.col(0).fill(ub(0));
					swarm.col(nSupp).fill(lb(1));
					swarm.col(nSupp + 1).fill(lb(1));
					swarm.cols(dSupp*nSupp, dSwarm - 1).fill(ini_w);
				}
				if (nSupp == 4) {
					swarm.cols(0,1).fill(lb(0));
					swarm.cols(2,3).fill(ub(0));
					swarm.col(nSupp + 1).fill(lb(1));
					swarm.col(nSupp + 3).fill(lb(1));
					swarm.col(nSupp + 2).fill(ub(1));						
				}				
			} else {
				if (nSupp == 4) {
					swarm.col(0).fill(ub(0));
					swarm.col(2).fill(ub(0));
					swarm.col(nSupp).fill(lb(1));
					swarm.col(nSupp + 1).fill(lb(1));	
					swarm.cols(dSupp*nSupp, dSwarm - 1).fill(ini_w);
				}
			}
			// WEIGHTS: EQUAL WEIGHTS
			swarm.submat(0, dSupp*nSupp, 0, dSwarm - 1).fill(ini_w);
		break;
		// ASSISTANT WEIGHT LOOP: WEIGHTS OF THE ASSISTANT DESIGN
		case 3:
			ini_w = (double)(1.0/(dSwarm + 1));
			swarm.row(0).fill(ini_w);
		break;
	}
}


