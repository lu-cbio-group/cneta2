#include "optimization.hpp"

// using namespace std;


/********************** Functions used in GSL optimization ************************/
double my_f(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    enew[i].length = exp(gsl_vector_get(v, i));
  }
  evo_tree new_tree(tree->nleaf, enew);
  if(lnl_type.model == MK){
      new_tree.mu = tree->mu;
  }else{
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}

double my_f_mu(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    enew[i].length = exp(gsl_vector_get(v, i));
  }
  evo_tree new_tree(tree->nleaf, enew);
  if(lnl_type.model == MK){
      new_tree.mu = exp(gsl_vector_get(v,tree->edges.size()-1));  // 0 to nedge-2 are epars, nedge - 1  is mu
  }else{
      new_tree.dup_rate = exp(gsl_vector_get(v,tree->edges.size()-1));
      new_tree.del_rate = exp(gsl_vector_get(v,tree->edges.size()));
      new_tree.chr_gain_rate = exp(gsl_vector_get(v,tree->edges.size() + 1));
      new_tree.chr_loss_rate = exp(gsl_vector_get(v,tree->edges.size()+2));
      new_tree.wgd_rate = exp(gsl_vector_get(v,tree->edges.size()+3));
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}


double my_f_cons(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }

  // parameters coming in are internal branch lengths followed by total time
  int count = 0;
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    if(enew[i].end > tree->nleaf - 1){
      enew[i].length = exp(gsl_vector_get(v, count));
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp(gsl_vector_get(v, count)));
  if(lnl_type.model == MK){
      new_tree.mu = tree->mu;
  }else{
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}

double my_f_cons_mu(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }

  // parameters coming in are internal branch lengths followed by mutation rate
  int count = 0;
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    if(enew[i].end > tree->nleaf - 1){
      enew[i].length = exp(gsl_vector_get(v, count));
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp(gsl_vector_get(v, count)));
  if(lnl_type.model == MK){
      new_tree.mu = exp(gsl_vector_get(v, count + 1));
  }
  else{
      new_tree.dup_rate = exp(gsl_vector_get(v, count + 1));
      new_tree.del_rate = exp(gsl_vector_get(v, count + 2));
      new_tree.chr_gain_rate = exp(gsl_vector_get(v, count + 3));
      new_tree.chr_loss_rate = exp(gsl_vector_get(v, count + 4));
      new_tree.wgd_rate = exp(gsl_vector_get(v, count + 5));
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}


// given a tree, maximise the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
// use GSL simplex optimization, less used than BFGS version
void max_likelihood(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double& min_nlnl, const double ssize){
  if(opt_type.bsr_mode > 0){
    cerr << "Error: simplex optimizer (optim=0) does not support bsr_mode>0. Use BFGS (optim=1) instead." << endl;
    exit(EXIT_FAILURE);
  }
  int debug = 0;
  vector<double> tobs = opt_type.tobs;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  /* save original handler, install new handler */
  gsl_error_handler_t  *old_handler = gsl_set_error_handler(&my_err_handler);

  int model = lnl_type.model;
  int cn_max = lnl_type.cn_max;
  int cn_type = lnl_type.cn_type;
  int is_total = lnl_type.is_total;

  int cons = lnl_type.cons;
  int estmu = opt_type.estmu;
  double tolerance = opt_type.tolerance;
  int miter = opt_type.miter;

  int nparams_est;
  int npar;
  gsl_vector *x;

  int nedge = 2 * rtree.nleaf - 2;
  int nintedge = rtree.nleaf - 2;

  if(!cons){
    nparams_est = nedge - 1 ;
    if(!estmu){
      npar = nparams_est;
    }else{  // estimate mutation rates
      if(model == MK){
          npar = nparams_est + 1;
      }
      if(model == BOUNDT){
          npar = nparams_est + 5;
      }
    }

    // initialise the best guess branch length and mu if required
    x = gsl_vector_alloc(npar);
    for(int i = 0; i < nparams_est; ++i){
      gsl_vector_set(x, i, log(rtree.edges[i].length));
    }
    if(estmu){
      if(model == MK){
          gsl_vector_set(x, nparams_est, log(rtree.mu));
      }
      if(model == BOUNDT){
          gsl_vector_set(x, nparams_est, log(rtree.dup_rate));
          gsl_vector_set(x, nparams_est + 1, log(rtree.del_rate));
          gsl_vector_set(x, nparams_est + 2, log(rtree.chr_gain_rate));
          gsl_vector_set(x, nparams_est + 3, log(rtree.chr_loss_rate));
          gsl_vector_set(x, nparams_est + 4, log(rtree.wgd_rate));
      }
      minex_func.f = my_f_mu;
    }else{
      minex_func.f = my_f;
    }
  }else{
    nparams_est = nintedge + 1;
    if(!estmu){
      npar = nparams_est;
    }else{
        if(model == MK){
            npar = nparams_est + 1;
        }
        if(model == BOUNDT){
            npar = nparams_est + 5;
        }
    }

    x = gsl_vector_alloc(npar);
    // initialise with internal edges
    vector<edge*> intedges = rtree.get_internal_edges();
    for(int i = 0; i < nintedge; ++i){
      gsl_vector_set(x, i, log(intedges[i]->length));
    }

    // initialise with current total tree time
    gsl_vector_set(x, nintedge, log(get_total_time(rtree.get_node_times(), lnl_type.max_tobs)));
    if(estmu){
      if(model == MK){
          gsl_vector_set(x, nparams_est, log(rtree.mu));
      }
      if(model == BOUNDT){
          gsl_vector_set(x, nparams_est, log(rtree.dup_rate));
          gsl_vector_set(x, nparams_est + 1, log(rtree.del_rate));
          gsl_vector_set(x, nparams_est + 2, log(rtree.chr_gain_rate));
          gsl_vector_set(x, nparams_est + 3, log(rtree.chr_loss_rate));
          gsl_vector_set(x, nparams_est + 4, log(rtree.wgd_rate));
      }
      minex_func.f = my_f_cons_mu;
    }else{
      minex_func.f = my_f_cons;
    }
  }

  // Set initial step sizes to 1
  //cout << "numbers: nedge / nintedge / npar / x->size: " << nedge << "\t" << nintedge << "\t" << npar << "\t" << x->size << endl;
  gsl_vector* ss = gsl_vector_alloc(npar);
  gsl_vector_set_all(ss, ssize);

  // Initialize method and iterate
  GSL_PARAM gsl_param = {rtree, vobs, lnl_type};
  GSL_PARAM* p = &gsl_param;
  void* pv = p;
  minex_func.n = npar;
  minex_func.params = pv;

  s = gsl_multimin_fminimizer_alloc(T, npar);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  size_t iter = 0;
  int status;
  double size;
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if(status) break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, tolerance);
    // gsl_set_error_handler_off();

    if(debug){
      printf("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
	      iter,
	      exp(gsl_vector_get(s->x, 0)),
	      exp(gsl_vector_get(s->x, s->x->size-1)),
	      s->fval, size);
    }

    if(status == GSL_SUCCESS){
      if(debug){
        	printf("converged to minimum at\n");

        	printf("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
        		iter,
        		exp(gsl_vector_get(s->x, 0)),
        		exp(gsl_vector_get(s->x, s->x->size-1)),
        		s->fval, size);
        }
    }
  }
  while(status == GSL_CONTINUE && iter < miter);

  if(status == GSL_CONTINUE ){
    cout << "### WARNING: maximum likelihood did not converge" << endl;
  }

  if(!cons){
    for(int i = 0; i < nparams_est; ++i){
      rtree.edges[i].length = exp(gsl_vector_get(s->x, i));
    }
    if(estmu){
        if(model == MK){
            rtree.mu = exp(gsl_vector_get(s->x, nparams_est));
        }
        if(model == BOUNDT){
            rtree.dup_rate = exp(gsl_vector_get(s->x, nparams_est));
            rtree.del_rate = exp(gsl_vector_get(s->x, nparams_est + 1));
            rtree.chr_gain_rate = exp(gsl_vector_get(s->x, nparams_est+2));
            rtree.chr_loss_rate = exp(gsl_vector_get(s->x, nparams_est+3));
            rtree.wgd_rate = exp(gsl_vector_get(s->x, nparams_est+4));
            rtree.mu = 0;
        }
    }

    min_nlnl = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    /* restore original handler */
    gsl_set_error_handler(old_handler);

  }else{
    int count = 0;
    // update internal lengths
    for(int i = 0; i < nedge - 1 ; ++i){
      if(rtree.edges[i].end > rtree.nleaf - 1){
        	rtree.edges[i].length = exp(gsl_vector_get(s->x, count));
        	count++;
      }else{
	        rtree.edges[i].length = 0;
      }
    }

    // Fill external edge lengths by looping over nodes
    double total_time = exp(gsl_vector_get(s->x, count));
    for(int i = 0; i < rtree.nleaf - 1; ++i){
      vector<int> es = rtree.get_ancestral_edges(rtree.nodes[i].id);
      reverse(es.begin(), es.end());

      rtree.edges[es.back()].length = total_time + tobs[ rtree.nodes[i].id ];
      for(int j = 0; j < es.size()-1; ++j){
        rtree.edges[es.back()].length -= rtree.edges[es[j]].length;
      }
    }

    if(estmu){
        if(model == MK){
            rtree.mu = exp(gsl_vector_get(s->x, nparams_est));
        }
        if(model == BOUNDT){
            rtree.dup_rate = exp(gsl_vector_get(s->x, nparams_est));
            rtree.del_rate = exp(gsl_vector_get(s->x, nparams_est + 1));
            rtree.chr_gain_rate = exp(gsl_vector_get(s->x, nparams_est+2));
            rtree.chr_loss_rate = exp(gsl_vector_get(s->x, nparams_est+3));
            rtree.wgd_rate = exp(gsl_vector_get(s->x, nparams_est+4));
            rtree.mu = 0;
        }
    }

    min_nlnl = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    /* restore original handler */
    gsl_set_error_handler(old_handler);

  }

}

/********************** End of functions used in GSL optimization ************************/


/*****************************************************
    One dimensional optimization with Brent method
*****************************************************/

/** 
 * @brief Compute the likelihood of the tree with one parameter (one branch or one mutation rate)
 * Used in Brent optimization
 * return negative likelihood function for minimalization
 * @param rtree the evolutionary tree
 * @param vobs the observed data
 * @param vobs_change the observed changes
 * @param obs_decomp the observed decompositions
 * @param comps the set of components
 * @param lnl_type the likelihood type
 * @param value the parameter value to be updated
 * @param type the type of parameter to be updated: -1: branch length; 0: mu; 1: dup_rate; 2: del_rate; 3: chr_gain_rate; 4: chr_loss_rate; 5: wgd_rate
 *  @return the negative likelihood value
 * TODO: to be extended to arm-level events
 */
double computeFunction(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double value, int type){
    // Update the branch
    switch(type){
        case -1:{
            rtree.current_it->length = value;
            rtree.current_it_back->length = value;
            // Need to update node times and ages?
            rtree.edges[rtree.current_eid].length = value;
            break;
        }
        case 0:{
            rtree.mu = value;
            break;
        }
        case 1:{
            rtree.dup_rate = value;
            break;
        }
        case 2:{
            rtree.del_rate = value;
            break;
        }
        case 3:{
            rtree.chr_gain_rate = value;
            break;
        }
        case 4:{
            rtree.chr_loss_rate = value;
            break;
        }
        case 5:{
            rtree.wgd_rate = value;
            break;
        }
        default:
            break;
    }

    double nlnl = 0.0;
    if(lnl_type.model == DECOMP){
        // nlnl = -get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        nlnl = -get_likelihood_change(rtree, vobs_change, obs_decomp, lnl_type, 0);  // debug = 0
    }else{
        nlnl = -get_likelihood_revised(rtree, vobs, lnl_type);
    }
    // cout << "\n result of computeFunction " << nlnl << endl;

    return nlnl;
}


#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d)(a)= (b);(b)= (c);(c)= (d);
#define SIGN(a,b)((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Brents method in one dimension */
double brent_opt(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double ax, double bx, double cx, double tol, double *foptx, double *f2optx, double fax, double fbx, double fcx){
	int iter;
	double a,b,d = 0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double xw,wv,vx;
	double e = 0.0;

	a= (ax < cx ? ax : cx);
	b= (ax > cx ? ax : cx);
	x=bx;
	fx=fbx;
	if(fax < fcx){
		w=ax;
		fw=fax;
		v=cx;
		fv=fcx;
	} else{
		w=cx;
		fw=fcx;
		v=ax;
		fv=fax;
	}

	for(iter=1;iter<=ITMAX;iter++){
		xm = 0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if(fabs(x-xm) <= (tol2 - 0.5*(b-a))){
			*foptx = fx;
			xw = x-w;
			wv = w-v;
			vx = v-x;
			*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
			         (v*v*xw + x*x*wv + w*w*vx);
			return x;
		}

		if(fabs(e) > tol1){
			r= (x-w)*(fx-fv);
			q= (x-v)*(fx-fw);
			p= (x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if(q > 0.0)
				p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e= (x >= xm ? a-x : b-x));
			else{
				d=p/q;
				u=x+d;
				if(u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else{
			d=CGOLD*(e= (x >= xm ? a-x : b-x));
		}

		u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu = computeFunction(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, u, type);
		if(fu <= fx){
			if(u >= x)
				a=x;
			else
				b=x;

			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else{
			if(u < x)
				a=u;
			else
				b=u;
			if(fu <= fw || w == x){
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else
				if(fu <= fv || v == x || v == w){
					v=u;
					fv=fu;
				}
		}
	}

	*foptx = fx;
	xw = x-w;
	wv = w-v;
	vx = v-x;
	*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/(v*v*xw + x*x*wv + w*w*vx);

	return x;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN
#undef GOLD
#undef GLIMIT
#undef TINY


double minimizeOneDimen(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double xmin, double xguess, double xmax, double tolerance, double *fx, double *f2x){
	double eps, optx, ax, bx, cx, fa, fb, fc;
	//int    converged;	/* not converged error flag */

	/* first attempt to bracketize minimum */
	if(xguess < xmin) xguess = xmin;
	if(xguess > xmax) xguess = xmax;
	eps = xguess*tolerance*50.0;
	ax = xguess - eps;
	if(ax < xmin) ax = xmin;
	bx = xguess;
	cx = xguess + eps;
	if(cx > xmax) cx = xmax;

	/* check if this works */
    // compute fb first to save some computation, if any
	fb = computeFunction(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, bx, type);
	fa = computeFunction(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, ax, type);
	fc = computeFunction(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, cx, type);

	/* if it works use these borders else be conservative */
	if((fa < fb) ||(fc < fb)){
		if(ax != xmin) fa = computeFunction(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, xmin, type);
		if(cx != xmax) fc = computeFunction(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, xmax, type);
		ax = xmin;
		cx = xmax;
	}
	/*
	const int MAX_ROUND = 10;
	for(i = 0;((fa < fb-tolerance) ||(fc < fb-tolerance)) &&(i < MAX_ROUND); i++){
		// brent method require that fb is smaller than both fa and fc
		// find some random values until fb achieve this
			bx = (((double)rand()) / RAND_MAX)*(cx-ax) + ax;
			fb = computeFunction(bx);
	}*/

/*
	if((fa < fb) ||(fc < fb)){
		if(fa < fc){ bx = ax; fb = fa; } else{ bx = cx; fb = fc; }
		//cout << "WARNING: Initial value for Brent method is set at bound " << bx << endl;
	}*/
	//	optx = brent_opt(xmin, xguess, xmax, tolerance, fx, f2x, fa, fb, fc);
	//} else
	optx = brent_opt(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, type, ax, bx, cx, tolerance, fx, f2x, fa, fb, fc);
  if(*fx > fb) // if worse, return initial value
  {
      *fx = computeFunction(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, bx, type);
      return bx;
  }

	return optx; /* return optimal x */
}




// return log likelihood
double optimize_mutation_rates(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double tolerance){
    int debug = 0;

    if(debug){
        cout << "\tUsing Brent method to optimize the likelihood of mutation rate" << endl;
    }

    double negative_lh = MAX_NLNL;
    double ferror, optx;

    // an anonymous function (lambda) to optimize one rate
    auto optimize_rate = [&](int idx, double& rate, const std::string& label) {
        optx = minimizeOneDimen(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, idx, RATE_MIN, rate, RATE_MAX, tolerance, &negative_lh, &ferror);

        rate = optx;

        if (debug) {
            cout << "\tmax Brent logl: " << negative_lh << " optimized " << label << " " << optx << endl;
        }
    };

    if(lnl_type.model == MK){
        optx = minimizeOneDimen(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, 0, RATE_MIN, rtree.mu, RATE_MAX, tolerance, &negative_lh, &ferror);
    }else{ // for other models
        switch (lnl_type.cn_type) {
        case ALL:{
            optimize_rate(1, rtree.dup_rate,      "duplication rate");
            optimize_rate(2, rtree.del_rate,      "deletion rate");
            optimize_rate(3, rtree.chr_gain_rate, "chromosome gain rate");
            optimize_rate(4, rtree.chr_loss_rate, "chromosome loss rate");
            optimize_rate(5, rtree.wgd_rate,      "WGD rate");
            break;
        }
        case ONLY_SEG:{
            optimize_rate(1, rtree.dup_rate, "duplication rate");
            optimize_rate(2, rtree.del_rate, "deletion rate");
            break;
        }
        case EXCLUDE_SEG:{
            optimize_rate(3, rtree.chr_gain_rate, "chromosome gain rate");
            optimize_rate(4, rtree.chr_loss_rate, "chromosome loss rate");
            optimize_rate(5, rtree.wgd_rate,      "WGD rate");
            break;
        }
        case EXCLUDE_CHR:{
            optimize_rate(1, rtree.dup_rate, "duplication rate");
            optimize_rate(2, rtree.del_rate, "deletion rate");
            optimize_rate(5, rtree.wgd_rate, "WGD rate");
            break;
        }
        case EXCLUDE_WGD:{
            optimize_rate(1, rtree.dup_rate,      "duplication rate");
            optimize_rate(2, rtree.del_rate,      "deletion rate");
            optimize_rate(3, rtree.chr_gain_rate, "chromosome gain rate");
            optimize_rate(4, rtree.chr_loss_rate, "chromosome loss rate");
            break;
        }
        default:
            cerr << "### ERROR: Unknown copy number alteration type for optimization!" << endl;
            exit(EXIT_FAILURE);
        }
    } 
    return -negative_lh;
}



// return log likelihood
double optimize_one_branch(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps,
LNL_TYPE& lnl_type, double tolerance, Node* node1, Node* node2){
    int debug = 0;
    if(debug){
        cout << "\tOptimizing the branch " << node1->id + 1 << ", " << node2->id + 1 << endl;
    }

    assert(!((node1->id == rtree.nleaf && node2->id == rtree.nleaf - 1) || (node2->id == rtree.nleaf && node1->id == rtree.nleaf - 1)));

    rtree.current_it = (Neighbor*) node1->findNeighbor(node2);
    assert(rtree.current_it);
    rtree.current_it_back = (Neighbor*) node2->findNeighbor(node1);
    assert(rtree.current_it_back);
    assert(rtree.current_it->length == rtree.current_it_back->length);

    int eid = rtree.get_edge_id(node1->id, node2->id);
    assert(eid >= 0);
    rtree.current_eid = eid;

    double current_len = rtree.edges[eid].length;
    assert(current_len >= 0.0);
    double negative_lh = MAX_NLNL;
    double ferror, optx;

    // Brent method
    optx = minimizeOneDimen(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, -1, BLEN_MIN, current_len, BLEN_MAX, tolerance, &negative_lh, &ferror);

    rtree.current_it->length = optx;
    rtree.current_it_back->length = optx;
    rtree.edges[eid].length = optx;

    if(debug){
        cout << "\tUsing Brent method to optimize the likelihood of one branch length of edge " << eid + 1 << endl;
        cout << "\tmax Brent logl: " << -negative_lh << " optimized branch length " << optx << endl;
    }

    return -negative_lh;
}


void compute_best_traversal(evo_tree& rtree, NodeVector &nodes, NodeVector &nodes2){
    Node* farleaf = rtree.find_farthest_leaf();

    // double call to farthest leaf to find the longest path on the tree
    rtree.find_farthest_leaf(farleaf);

    rtree.get_preorder_branches(nodes, nodes2, farleaf);
}


double optimize_all_branches(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps,
LNL_TYPE& lnl_type, int my_iterations, double tolerance){
    int debug = 0;

    NodeVector nodes, nodes2;
    compute_best_traversal(rtree, nodes, nodes2);

    if(debug){
        cout << "\nOptimizing branch lengths (max " << my_iterations << " loops)..." << endl;
        cout << "nodes in best traversal: ";
        for(auto n : nodes){
            cout << "\t" << n->id + 1;
        }
        cout << endl;
        cout << "nodes2 in best traversal: ";
        for(auto n : nodes2){
            cout << "\t" << n->id + 1;
        }
        cout << endl;
    }

    int model = lnl_type.model;
    double tree_lh = rtree.score;

    // if(model == DECOMP){
    //     tree_lh = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
    // }else{
    //     tree_lh = get_likelihood_revised(rtree, vobs, lnl_type);
    // }
    if(debug){
        cout << "Initial tree log-likelihood: " << tree_lh << endl;
    }

    DoubleVector lenvec;
    double new_tree_lh = 0.0;
    for(int i = 0; i < my_iterations; i++){
    	save_branch_lengths(rtree, lenvec, 0);

        for(int j = 0; j < nodes.size(); j++){
            if(debug){
                cout << " optimizing branch " << nodes[j]->id << " " << nodes2[j]->id << endl;
            }
            // skip normal branch, which will always be 0
            if(nodes[j]->id == rtree.nleaf - 1 || nodes2[j]->id == rtree.nleaf - 1)   continue;
            new_tree_lh = optimize_one_branch(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, tolerance, (Node* )nodes[j],(Node* )nodes2[j]);
        }

        if(debug){
            cout << "tree log likelihood after iteration " << i + 1 << " : " << new_tree_lh << endl;
        }

        if(new_tree_lh < tree_lh - tolerance * 0.1){
            // IN RARE CASE: tree log-likelihood decreases, revert the branch length and stop
            if(debug){
                cout << "tree log-likelihood decreases" << endl;
                cout << "NOTE: Restoring branch lengths as tree log-likelihood decreases after branch length optimization: " << tree_lh << " -> " << new_tree_lh << endl;
            }

            restore_branch_lengths(rtree, lenvec);

            double max_delta_lh = 1.0;

            if(model == DECOMP){
              // new_tree_lh = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
              new_tree_lh = get_likelihood_change(rtree, vobs_change, obs_decomp, lnl_type, 0);  // debug = 0
            }else{
              new_tree_lh = get_likelihood_revised(rtree, vobs, lnl_type);
            }

            if(fabs(new_tree_lh - tree_lh) > max_delta_lh){
              cout << endl;
              cout << "new_tree_lh: " << new_tree_lh << "   tree_lh: " << tree_lh << endl;
            }
        	assert(fabs(new_tree_lh - tree_lh) < max_delta_lh);

        	return new_tree_lh;
        }

        // only return if the new_tree_lh >= tree_lh!(in rare case that likelihood decreases, continue the loop)
        if(tree_lh <= new_tree_lh && new_tree_lh <= tree_lh + tolerance){
          if(debug) cout << "tree log-likelihood increases" << endl;
        	return new_tree_lh;
        }

        tree_lh = new_tree_lh;
    }

    // recompute node times and ages to be consistent with branch lengths
    rtree.calculate_node_times();
    rtree.calculate_age_from_time();

    if(debug) cout << "current score " << tree_lh << endl;

    return tree_lh;
}


/*****************************************************
    L-BFGS-B method
*****************************************************/
// not used so far
void update_variables(evo_tree& rtree, int model, int cons, int estmu, double *x){
    int debug = 0;

    int nedge = 2 * rtree.nleaf - 2;
    // create a new tree from current value of parameters
    vector<edge> enew;
    for(int i = 0; i < nedge; ++i){
      enew.push_back(rtree.edges[i]);
    }

    if(!cons){
      // The first element of x is not used for optimization
      // The index of x is added by 1 compared with index for simplex method
      for(int i = 0; i < nedge - 1 ; ++i){
          enew[i].length = x[i + 1];
      }
      evo_tree new_tree(rtree.nleaf, enew);
      if(estmu){
          if(model == MK){
              new_tree.mu = x[nedge];
              if(debug){
                  for(int i = 0; i < nedge + 1; i++){ cout << x[i] << '\n';}
                  cout << "mu value so far: " << new_tree.mu << endl;
              }
          }else{
              new_tree.dup_rate = x[nedge];
              new_tree.del_rate = x[nedge + 1];
              new_tree.chr_gain_rate = x[nedge + 2];
              new_tree.chr_loss_rate = x[nedge + 3];
              new_tree.wgd_rate = x[nedge + 4];
              new_tree.mu = 0;

              if(debug){
                  for(int i = 0; i < nedge + 2; i++){ cout << x[i] << '\n';}
                  cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                  cout << "del_rate value so far: " << new_tree.del_rate << endl;
                  cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                  cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                  cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
              }
          }
      }
      rtree = evo_tree(new_tree);
    }else{  // constrained optimization
      int count = 0;
      for(int i = 0; i < nedge - 1 ; ++i){
        if(enew[i].end > rtree.nleaf - 1){
          	enew[i].length = x[count + 1];
          	count++;
        }else{
  	        enew[i].length = 0;
        }
      }
      if(debug){
          cout << "total tree height so far: " << x[count + 1] << endl;
      }
      evo_tree new_tree(rtree.nleaf, enew, x[count + 1]);
      if(estmu){
        int nintedge = rtree.nleaf - 2;

        if(model == MK){
            new_tree.mu = x[nintedge + 2];
            if(debug){
                for(int i = 0; i <= nintedge + 2; i++){ cout << x[i] << '\n';}
                cout << "mu value so far: " << new_tree.mu << endl;
            }
        }else{
            new_tree.dup_rate = x[nintedge + 2];
            new_tree.del_rate = x[nintedge + 3];
            new_tree.chr_gain_rate = x[nintedge + 4];
            new_tree.chr_loss_rate = x[nintedge + 5];
            new_tree.wgd_rate = x[nintedge + 6];
            new_tree.mu = 0;

            if(debug){
                for(int i = 0; i <= nintedge + 6; i++){ cout << x[i] << '\n';}
                cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                cout << "del_rate value so far: " << new_tree.del_rate << endl;
                cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
            }
         }
      }
      rtree = evo_tree(new_tree);
    }
}

void update_tree_rates(LNL_TYPE & lnl_type, evo_tree & rtree, double *x, int nparams_est, int bsr_mode){
    // Global reference rates must stay frozen when bsr_mode>0; this assert fails loudly if a
    // future edit re-enables one of the disabled call sites elsewhere in this file.
    assert(bsr_mode == 0);
    switch (lnl_type.cn_type)
    {
        case ALL:
            rtree.dup_rate = x[nparams_est + 1];
            rtree.del_rate = x[nparams_est + 2];
            rtree.chr_gain_rate = x[nparams_est + 3];
            rtree.chr_loss_rate = x[nparams_est + 4];
            rtree.wgd_rate = x[nparams_est + 5];
            break;

        case EXCLUDE_WGD:
            rtree.dup_rate = x[nparams_est + 1];
            rtree.del_rate = x[nparams_est + 2];
            rtree.chr_gain_rate = x[nparams_est + 3];
            rtree.chr_loss_rate = x[nparams_est + 4];
            break;

        case ONLY_SEG:
            rtree.dup_rate = x[nparams_est + 1];
            rtree.del_rate = x[nparams_est + 2];
            break;

        case EXCLUDE_SEG:
            rtree.chr_gain_rate = x[nparams_est + 1];
            rtree.chr_loss_rate = x[nparams_est + 2];
            rtree.wgd_rate = x[nparams_est + 3];
            break;

        case EXCLUDE_CHR:
            rtree.dup_rate = x[nparams_est + 1];
            rtree.del_rate = x[nparams_est + 2];
            rtree.wgd_rate = x[nparams_est + 3];
            break;

        default:
            cerr << "Unknown CN type" << endl;
            exit(EXIT_FAILURE);
    }
}

// Update the tree after each iteration in the BFGS optimization
// Estimate ratios rather than branch length in order to avoid negative terminal branch lengths
// Sort node times in increasing order and take the first Ns intervals
void update_variables_transformed(evo_tree& rtree, double *x, LNL_TYPE& lnl_type, OPT_TYPE& opt_type){
    int debug = 0;
    if(debug){
        cout << "Update the tree after each iteration in the BFGS optimization" << endl;
        cout << "Tree before optimization " << rtree.make_newick() << endl;
        cout << "All the node ages so far: ";
        for(auto n : rtree.nodes){
          cout << "\t" << n.age;
        }
        cout << endl;
    }

    auto print_rate = [&](const char* name, double value) {
        cout << name << " value so far: " << value << endl;
    };


    int cn_max = lnl_type.cn_max;
    int cn_type = lnl_type.cn_type;
    int is_total = lnl_type.is_total;
    int nparams_est = 0;

    if(opt_type.opt_one_branch){
        if(debug){
          cout << "transform only one branch " << rtree.current_eid + 1 << endl;
        }
        if(!lnl_type.cons){
            nparams_est = 1;
            rtree.edges[rtree.current_eid].length = x[1];
            rtree.current_it->length = x[1];
            rtree.current_it_back->length = x[1];
        }else{
            nparams_est = 1;
            vector<double> ratios = rtree.get_ratio_from_age();

            // ratios[0] = x[1]; // need to update root edge
            int eend = rtree.edges[rtree.current_eid].end;
            if(eend <= rtree.nleaf){
                cerr << "Error: constrained single-branch optimisation requires an internal "
                     << "non-root edge (eend=" << eend << ", nleaf=" << rtree.nleaf << ")" << endl;
                exit(EXIT_FAILURE);
            }
            int nid = eend - rtree.root_node_id;
            ratios[nid] = x[1];
            // update based on all ratios together to avoid inconsistencies
            rtree.update_edges_from_ratios(ratios, lnl_type.knodes);
            // rtree.update_edge_from_ratio(x[1], rtree.current_eid);
            // get length of optimized edge
            double blen = rtree.edges[rtree.current_eid].length;
            rtree.current_it->length = blen;
            rtree.current_it_back->length = blen;

            if(debug){
              cout << "estimated x: ";
              for(int i = 0; i < nparams_est; i++){
                  double val = x[i + 1];
                  cout << "\t" << val;
              }
              cout << endl;
              cout << "tree after optimization " << rtree.make_newick() << endl;
              cout << "all the node ages so far: ";
              for(auto n : rtree.nodes){
                cout << "\t" << n.age;
              }
              cout << endl;
            } // debug
        } // else
    }else{
        if(!lnl_type.cons){
            // nparams_est = nedge - 1;
            nparams_est = 2 * rtree.nleaf - 3;
            // The first element of x is not used for optimization
            // The index of x is added by 1 compared with index for simplex method
            for(int i = 0; i < nparams_est; i++){
              rtree.edges[i].length = x[i + 1];
            }
        }else{
            // nparams_est = nintedge + 1;
            nparams_est = rtree.nleaf - 1;          // number of branches

            // store original ratios
            vector<double> ratios = rtree.get_ratio_from_age();

            double min_root = lnl_type.max_tobs + rtree.nleaf * BLEN_MIN;

            if(debug){
                int nratios = rtree.nleaf - 1;
                cout << "There are " << nratios << " ratios" << endl;
                cout << "Original values of estimated variables: " << endl;
                for(int i = 0; i < nratios; i++){
                    cout << i + 1 << "\t" << "\t" << ratios[i] << endl;
                }
                cout << "estimated x: ";
                for(int i = 0; i < nparams_est; i++){
                    double val = x[i + 1];
                    cout << "\t" << val;
                }
                cout << endl;
                cout << "min age of root allowed " << min_root << endl;
            } // debug

            // The estimated value may be nan
            for(int i = 0; i < nparams_est; i++){
                double val = x[i + 1];
                if(std::isnan(val)){  // return previous values
                    if(debug) cout << "nan returned in BFGS!" << endl;
                    val = ratios[i];
                }
                ratios[i] = val;
            }

            // If the optimized root is close to the boundary, revert to original values to avoid very small branches
            if(fabs(x[1] - min_root) > SMALL_VAL){
              // cout << "\n\nupdate all branches" << endl;
              rtree.update_edges_from_ratios(ratios, lnl_type.knodes);
            }

            if(!is_age_time_consistent(rtree.get_node_times(), rtree.get_node_ages())){
              cout << "Updated tree has inconsistent time/ages after updating the optimized tree!" << endl;
              rtree.print();
              string newick = rtree.make_newick(8);
              cout << newick << endl;
            }

            if(debug){
                cout << "Current values of estimated ratio variables: " << endl;
                for(int i = 0; i < ratios.size(); i++){
                    cout << i + 1 << "\t" << "\t" << ratios[i] << endl;
                }
            } // debug
          } // else

        if(debug){
          cout << "New branch lengths: \n";
          for(int i = 0; i < rtree.edges.size(); i++){
              cout << i + 1 << "\t" << "\t" << rtree.edges[i].length << endl;
          }
        } // debug
    }

    if(lnl_type.cons && !is_tip_age_valid(rtree.get_node_ages(), opt_type.tobs)){
      rtree.print();
      cout << rtree.make_newick() << endl;
      cout << "Tip timings inconsistent with observed data after updating the optimized tree!" << endl;
      exit(EXIT_FAILURE);
    }

    if(opt_type.bsr_mode == 1){
        // edge_rates[eid] = global_rates * m_i (one shared multiplier per branch)
        vector<int> active_eids;

        // cout << "\n[BSR1 DEBUG] before get_active_bsr_eids" << endl;
        // cout << "global rates before active_eids:" << endl;
        // cout << "  dup_rate      = " << rtree.dup_rate << endl;
        // cout << "  del_rate      = " << rtree.del_rate << endl;
        // cout << "  chr_gain_rate = " << rtree.chr_gain_rate << endl;
        // cout << "  chr_loss_rate = " << rtree.chr_loss_rate << endl;
        // cout << "  wgd_rate      = " << rtree.wgd_rate << endl;

        // if(!opt_type.opt_one_branch)
        active_eids = get_active_bsr_eids(rtree);

        // cout << "\n[BSR1 DEBUG] after get_active_bsr_eids" << endl;
        // cout << "active_eids.size() = " << active_eids.size() << endl;

        // int nprint = min(5, (int)active_eids.size());
        // cout << "first active eids: ";
        // for(int k = 0; k < nprint; k++){
        //     cout << active_eids[k] << " ";
        // }
        // cout << endl;     

        // [2026-07-14 disabled] read-back counterpart of the set_tree_rates removal
        // above: this let each BFGS evaluation overwrite the frozen global reference
        // rates with whatever was in the per-edge multiplier slots. See
        // docs/flowcharts.md "Two-step calibration workflow".
        // // update global rate for all active edges
        // update_tree_rates(lnl_type, rtree, x, nparams_est);

        RateSet global_rates(0, rtree.dup_rate, rtree.del_rate,
                             rtree.chr_gain_rate, rtree.chr_loss_rate, rtree.wgd_rate);

        for(int k = 0; k < (int)active_eids.size(); k++){
            int eid = active_eids[k];
            double m = x[nparams_est + k + 1];
            rtree.edge_rates[eid] = global_rates * m;
        }
 
        // cout << "\n[BSR1 DEBUG] after edge_rates update" << endl;
        // for(int k = 0; k < nprint; k++){
        //     int eid = active_eids[k];
        //     int idx = nparams_est + k + 1;
        //     cout << "  eid=" << eid
        //         << " m=x[" << idx << "]=" << x[idx]
        //         << " edge dup=" << rtree.edge_rates[eid].dup
        //         << " edge del=" << rtree.edge_rates[eid].del
        //         << " edge chr_gain=" << rtree.edge_rates[eid].chr_gain
        //         << " edge chr_loss=" << rtree.edge_rates[eid].chr_loss
        //         << " edge wgd=" << rtree.edge_rates[eid].wgd
        //         << endl;
        // }        
    } else if(opt_type.bsr_mode == 3){
        // RLC: root-to-tip propagation with inheritance; shift edges get parent_rate * m_k
        // Keep existing edge_rates frozen; do not attempt to read multipliers from x[].
        // [2026-07-14 disabled] temporary diagnostic prints used to confirm the bug below:
        // fired on every single BFGS objective/gradient evaluation (thousands of times per
        // optimization, across many parallel candidates), so left uncommented they flood
        // the log and add real cout-lock contention under OpenMP. Commented out (not
        // deleted) rather than removed, so a later reader puzzled by huge log files from
        // this period of the project knows what produced them and why they were disabled.
        // cout << "[DEBUG-CHECK] K=" << opt_type.rlc_shift_eids.size()
        //      << " dup_rate BEFORE=" << rtree.dup_rate
        //      << " x[nparams_est+1]=" << x[nparams_est + 1] << endl;
        // [2026-07-14 disabled] this call used to overwrite the frozen global dup_rate/
        // del_rate/... with whatever was in this shift edge's multiplier slots (confirmed
        // 2026-07-14 via the temporary debug print above: dup_rate jumped from 0.00143 to
        // 1 as soon as K>=1). See docs/flowcharts.md "Two-step calibration workflow".
        // update_tree_rates(lnl_type, rtree, x, nparams_est);
        // cout << "[DEBUG-CHECK] dup_rate AFTER=" << rtree.dup_rate << endl;

        update_edge_rates_rlc_from_x(rtree, opt_type.rlc_shift_eids, x, nparams_est, opt_type.bsr_slots_cache);

    } else if(opt_type.bsr_mode == 2){
        // edge_rates[eid].field_t = global_rate_t * m_t  (independent multiplier per type per branch)
        // idx(k, t) = nparams_est + k * n_types_per_edge + t + 1
        // [2026-07-14 disabled] same bug as bsr_mode==1/3 above: this let the frozen
        // global reference rates be overwritten by the first edge's per-type multiplier
        // values. See docs/flowcharts.md "Two-step calibration workflow".
        // update_tree_rates(lnl_type, rtree, x, nparams_est);
        // [2026-07-14 disabled] recomputed bsr_slots on every call; it's invariant for the
        // whole optimization, so max_likelihood_BFGS now computes it once into
        // opt_type.bsr_slots_cache instead. See docs/flowcharts.md.
        // auto bsr_slots = get_bsr_rate_slots(rtree, cn_type);
        const vector<BsrRateSlot>& bsr_slots = opt_type.bsr_slots_cache;
        int n_types = (int)bsr_slots.size();
        vector<int> active_eids;
        if(!opt_type.opt_one_branch)
            active_eids = get_active_bsr_eids(rtree);
        for(int k = 0; k < (int)active_eids.size(); k++){
            int eid = active_eids[k];
            for(int t = 0; t < n_types; t++){
                // [2026-07-14 disabled] see bsr_var_index in optimization.hpp
                // double m = x[nparams_est + k * n_types + t + 1];
                double m = x[bsr_var_index(nparams_est, k, n_types, t)];
                rtree.edge_rates[eid].*(bsr_slots[t].edge_field) = rtree.*(bsr_slots[t].tree_rate) * m;
            }
        }
    } else if(opt_type.estmu){
        if(lnl_type.model == MK){
            // nintedge + 2 for constrained branches
            rtree.mu = x[nparams_est + 1];
            if(debug){
                for(int i = 0; i <= nparams_est + 1; i++){ cout << x[i] << '\n'; }
                cout << "mu value so far: " << rtree.mu << endl;
            }
        }else{ // for other models
          // only update those estimated rates
          update_tree_rates(lnl_type, rtree, x, nparams_est, opt_type.bsr_mode);

          if(debug){
              switch(cn_type){
                case ALL:{
                    print_rate("dup_rate",      rtree.dup_rate);
                    print_rate("del_rate",      rtree.del_rate);
                    print_rate("chr_gain_rate", rtree.chr_gain_rate);
                    print_rate("chr_loss_rate", rtree.chr_loss_rate);
                    print_rate("wgd_rate",      rtree.wgd_rate);
                    break;
                }
                case ONLY_SEG:{
                    print_rate("dup_rate", rtree.dup_rate);
                    print_rate("del_rate", rtree.del_rate);
                    break;
                }
                case EXCLUDE_SEG:{
                    print_rate("chr_gain_rate", rtree.chr_gain_rate);
                    print_rate("chr_loss_rate", rtree.chr_loss_rate);
                    print_rate("wgd_rate",      rtree.wgd_rate);
                    break;
                }
                case EXCLUDE_CHR:{
                    print_rate("dup_rate", rtree.dup_rate);
                    print_rate("del_rate", rtree.del_rate);
                    print_rate("wgd_rate", rtree.wgd_rate);
                    break;
                }
                case EXCLUDE_WGD:{
                    print_rate("dup_rate",      rtree.dup_rate);
                    print_rate("del_rate",      rtree.del_rate);
                    print_rate("chr_gain_rate", rtree.chr_gain_rate);   
                    print_rate("chr_loss_rate", rtree.chr_loss_rate);
                    break;
                }
                default:
                    cerr << "Unknown cn_type while printing rates" << endl;
                    exit(EXIT_FAILURE);
              }
          } // debug
        }
    }
}



/**
    the target function which needs to be optimized (negative log likelihood)
    @param x the input vector x
    @return the function value at x (negative log likelihood)
*/
double targetFunk(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double x[]){
  int debug = 0;
  update_variables_transformed(rtree, x, lnl_type, opt_type);

  double nlnl;
  if(lnl_type.model == DECOMP){
    // nlnl = -1.0 * get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
    if(debug) cout << "checking targetFunk" << endl;
    if(opt_type.bsr_mode > 0){  // variable rate: all modes use per-edge rates in edge_rates
      nlnl = -1.0 * get_likelihood_change_variable_rate(rtree, vobs_change, obs_decomp, lnl_type, 0);
    } else {  // constant rate: single global rate shared across all branches
      nlnl = -1.0 * get_likelihood_change(rtree, vobs_change, obs_decomp, lnl_type, 0);
    }
  }else{
    nlnl = -1.0 * get_likelihood_revised(rtree, vobs, lnl_type);
  }

  if(!isfinite(nlnl)){
    cout << "[DEBUG] non-finite nlnl value check in targetFunk, nlnl = " << nlnl << endl;
    if(opt_type.bsr_mode > 0){
      cout << "        bsr_mode=" << opt_type.bsr_mode
           << " global dup_rate=" << rtree.dup_rate
           << " del_rate=" << rtree.del_rate
           << " chr_gain_rate=" << rtree.chr_gain_rate
           << " chr_loss_rate=" << rtree.chr_loss_rate
           << " wgd_rate=" << rtree.wgd_rate << endl;
      if(!rtree.edge_rates.empty()){
        cout << "        edge_rates[0]: dup=" << rtree.edge_rates[0].dup
             << " del=" << rtree.edge_rates[0].del
             << " chr_gain=" << rtree.edge_rates[0].chr_gain
             << " chr_loss=" << rtree.edge_rates[0].chr_loss
             << " wgd=" << rtree.edge_rates[0].wgd << endl;
      }
    } else {
      cout << "        dup_rate=" << rtree.dup_rate
           << " del_rate=" << rtree.del_rate
           << " chr_gain_rate=" << rtree.chr_gain_rate
           << " chr_loss_rate=" << rtree.chr_loss_rate
           << " wgd_rate=" << rtree.wgd_rate << endl;
    }
    int max_to_print = 10;   // print first 10 variables of x[i], which are
    cout << "        first "<< max_to_print <<" x: ";
    for (int i = 1; i <= max_to_print; ++i) {
      cout << " x[" << i << "]=" << x[i];
    }
    cout << endl;
    exit(EXIT_FAILURE); //nlnl = MAX_NLNL;// assign a large value, so that optimizer can continue, DEBUG ONLY
  }
  return nlnl;
}



/**
	the approximated derivative function
	@param x variables to estimate 
	@param dfx the derivative at x
	@return the function value at x
*/
double derivativeFunk(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int ndim, double x[], double dfx[]){
    int debug = 0;

    double *h = new double[ndim + 1];
    double temp;
    int dim;

    double fx = targetFunk(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, x);

	for(dim = 1; dim <= ndim; dim++){
        temp = x[dim];
        h[dim] = ERROR_X * fabs(temp);
        if(h[dim] == 0.0) h[dim] = ERROR_X;
        x[dim] = temp + h[dim];
        h[dim] = x[dim] - temp;

        dfx[dim] = (targetFunk(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, x));
        x[dim] = temp;
	}

	for(dim = 1; dim <= ndim; dim++){
        dfx[dim] = (dfx[dim] - fx) / h[dim];
        if(debug){
            cout << "dfx[dim] " << dfx[dim] << endl;
        }
    }

    delete [] h;

    // restore rtree (including edge_rates for bsr_mode>0) to the state at x,
    // so callers see a consistent tree after gradient computation
    targetFunk(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, x);

	return fx;
}


double optimFunc(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *vars){
    return targetFunk(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, vars-1);
}

double optimGradient(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *x, double *dfx){
    return derivativeFunk(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, nvar, x - 1, dfx - 1);
}


// internal function to interface with L-BFGS-B
void lbfgsb(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT){
	char task[60];
	double f, *g, dsave[29], *wa;
	int tr = -1, iter = 0, *iwa, isave[44], lsave[4];

	/* shut up gcc -Wall in 4.6.x */

	for(int i = 0; i < 4; i++) lsave[i] = 0;

	if(n == 0){ /* not handled in setulb */
		*fncount = 1;
		*grcount = 0;

    *Fmin = optimFunc(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, n, u);
		strcpy(msg, "NOTHING TO DO");
		*fail = 0;
		return;
	}

	if(nREPORT <= 0){
		cerr << "REPORT must be > 0(method = \"L-BFGS-B\")" << endl;
		exit(EXIT_FAILURE);
	}

	switch(trace){
        case 2: tr = 0; break;
        case 3: tr = nREPORT; break;
        case 4: tr = 99; break;
        case 5: tr = 100; break;
        case 6: tr = 101; break;
        default: tr = -1; break;
	}

	*fail = 0;
	g = (double*) malloc(n * sizeof(double));
	/* this needs to be zeroed for snd in mainlb to be zeroed */
	wa = (double *) malloc((2*m*n+4*n + 11*m*m+8*m) * sizeof(double));
	iwa = (int *) malloc(3*n * sizeof(int));
	strcpy(task, "START");

	while(1){
		/* Main workhorse setulb() from ../appl/lbfgsb.c : */
		setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task, tr, lsave, isave, dsave);
		/*    Rprintf("in lbfgsb - %s\n", task);*/
		if(strncmp(task, "FG", 2) == 0){
            f = optimGradient(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, n, x, g);
			if(!isfinite(f)){
				cerr << "L-BFGS-B needs finite values of 'fn'" << endl;
				exit(EXIT_FAILURE);
			}
		}else if(strncmp(task, "NEW_X", 5) == 0){
			iter++;
			if(trace == 1 &&(iter % nREPORT == 0)){
				cout << "iter " << iter << " value " << f << endl;
			}
			if(iter > maxit){
				*fail = 1;
				break;
			}
		}else if(strncmp(task, "WARN", 4) == 0){
			*fail = 51;
			break;
		}else if(strncmp(task, "CONV", 4) == 0){
			break;
		}else if(strncmp(task, "ERROR", 5) == 0){
			*fail = 52;
			break;
		}else{ /* some other condition that is not supposed to happen */
			*fail = 52;
			break;
		}
	}

	*Fmin = f;
	*fncount = *grcount = isave[33];
	if(trace){
		cout << "final value " << *Fmin << endl;
		if(iter < maxit && *fail == 0)
			cout << "converged" << endl;
		else
			cout << "stopped after " << iter << " iterations\n";
	}
	strcpy(msg, task);

	free(g);
	free(wa);
	free(iwa);
}


/**
 Function to access the L-BFGS-B function, taken from IQ-TREE package which is further taken from HAL_HAS software package
 1. int nvar or n : The number of the variables
 2. double* vars or x : initial values of the variables
 3. double* lower or l: lower bounds of the variables
 4. double* upper or u : upper bounds of the variables
 5. double pgtol: gradient tolerance
 5. int maxit : max # of iterations
 @return minimized function value
 After the function is invoked, the values of x will be updated
*/
double L_BFGS_B(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, double* x, double* l, double* u){
  int debug = 0;

  int i;
  double Fmin;
  int fail;
  int fncount;
  int grcount;
  char msg[100];

  int m = 10;          // number of BFGS updates retained in the "L-BFGS-B" method. It defaults to 5.

  int *nbd;           // 0: unbounded; 1: lower bounded; 2: both lower & upper; 3: upper bounded
  nbd = new int[n];
  for(i = 0; i < n; i++) nbd[i] = 2;

  double factr = 1e+7; // control the convergence of the "L-BFGS-B" method.
  // Convergence occurs when the reduction in the object is within this factor
  // of the machine tolerance.
  // Default is 1e7, that is a tolerance of about 1e-8

  double pgtol = opt_type.tolerance;
  double maxit = opt_type.miter;

  //	double pgtol = 0;   // helps control the convergence of the "L-BFGS-B" method.
  //    pgtol = 0.0;
  // It is a tolerance on the projected gradient in the current search direction.
  // Default is zero, when the check is suppressed

  int trace = 0;      // non-negative integer.
  if(debug)
      trace = 1;
  // If positive, tracing information on the progress of the optimization is produced.
  // Higher values may produce more tracing information.

  int nREPORT = 10;   // The frequency of reports for the "L-BFGS-B" methods if "trace" is positive.
  // Defaults to every 10 iterations.

  /*#ifdef USE_OLD_PARAM
  lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr1, &fail, ex,
  	factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
  #else*/

  // cout << "initial values in bfgs: ";
  // for(int mi = 0; mi < 3; mi++){
  //     cout << "\t" << x[mi];
  // }
  // cout << "\n";
  //
  // cout << "upper bound in bfgs: ";
  // for(int mi = 0; mi < 3; mi++){
  //     cout << "\t" << u[mi];
  // }
  // cout << "\n";

  lbfgsb(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, n, m, x, l, u, nbd, &Fmin, &fail, factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
  //#endif

  if(fail == 51 || fail == 52){
      cout << msg << endl;
  }

  delete[] nbd;

  return Fmin;
}

void set_param(int offset, double value, int nparams_est, double* variables, double* lower_bound, double* upper_bound){
    int idx = nparams_est + offset;
    variables[idx] = value;
    lower_bound[idx] = MIN_MRATE;
    upper_bound[idx] = MAX_MRATE;
}

void set_tree_rates(int cn_type, evo_tree &rtree, int nparams_est, double* variables, double* lower_bound, double* upper_bound, int bsr_mode){
    // Global reference rates must stay frozen when bsr_mode>0, never seeded as free BFGS
    // variables; this assert fails loudly if a future edit re-enables one of the disabled
    // call sites elsewhere in this file.
    assert(bsr_mode == 0);
    switch (cn_type)
    {
    case ALL:
    {
        set_param(1, rtree.dup_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(2, rtree.del_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(3, rtree.chr_gain_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(4, rtree.chr_loss_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(5, rtree.wgd_rate, nparams_est, variables, lower_bound, upper_bound);
        break;
    }
    case ONLY_SEG:
    {
        set_param(1, rtree.dup_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(2, rtree.del_rate, nparams_est, variables, lower_bound, upper_bound);
        break;
    }
    case EXCLUDE_SEG:
    {
        set_param(1, rtree.chr_gain_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(2, rtree.chr_loss_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(3, rtree.wgd_rate, nparams_est, variables, lower_bound, upper_bound);
        break;
    }
    case EXCLUDE_CHR:
    {
        set_param(1, rtree.dup_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(2, rtree.del_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(3, rtree.wgd_rate, nparams_est, variables, lower_bound, upper_bound);
        break;
    }
    case EXCLUDE_WGD:
    {
        set_param(1, rtree.dup_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(2, rtree.del_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(3, rtree.chr_gain_rate, nparams_est, variables, lower_bound, upper_bound);
        set_param(4, rtree.chr_loss_rate, nparams_est, variables, lower_bound, upper_bound);
        break;
    }
    default:
        cerr << "Error: unknown cn_type when optimizing mutation rates!" << endl;
        exit(EXIT_FAILURE);
    }
}

void max_likelihood_BFGS(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double &min_nlnl, int debug){
    // debug = 1;
    // initialize variables for L-BFGS-B optimization
    int model = lnl_type.model;
    int cn_max = lnl_type.cn_max;
    int cn_type = lnl_type.cn_type;
    int is_total = lnl_type.is_total;
    int patient_age = lnl_type.patient_age;
    int cons = lnl_type.cons;

    int estmu = opt_type.estmu;
    int opt_one_branch = opt_type.opt_one_branch;

    // bsr_mode>0 + estmu=1 is not supported: global rates serve as fixed reference values;
    // per-branch multipliers m_i are estimated instead (actual rate on branch i = global_rate * m_i).
    // Estimating both global rates and multipliers simultaneously is non-identifiable.
    if(opt_type.bsr_mode > 0 && estmu){
        cerr << "Error: estmu=1 is not compatible with bsr_mode>0 (global rates are fixed as reference; per-branch multipliers are estimated instead)" << endl;
        exit(EXIT_FAILURE);
    }

    // Set variables
    int nparams_est = 0;    // number of parameters to estimate

    if(opt_one_branch){
      if(debug){
        cout << "Updating only one branch " << rtree.current_eid + 1 << endl;
      }
      nparams_est = 1;
    }else{
      if(!cons){  // only estimate internal branches
          nparams_est = 2 * rtree.nleaf - 3;
      }else{    // estimate internal branches and mutation rate
          nparams_est = rtree.nleaf - 1;
      }
    }

    // build active edge list for variable rate (excludes the normal-sample edge)
    // empty when opt_one_branch=1, so no multipliers enter the optimisation in that case
    vector<int> active_eids;
    if(opt_type.bsr_mode > 0 && !opt_one_branch){
        if(opt_type.bsr_mode == 3)
            active_eids = opt_type.rlc_shift_eids;   // only the δ=1 shift edges
        else
            active_eids = get_active_bsr_eids(rtree); // all optimizable edges
    }
    int n_bsr_edges = (int)active_eids.size();  // correct count regardless of constraints

    // Active rate slots for BSR: which rate types participate (depends on cn_type).
    // bsr_mode=1: shared multiplier — n_types_per_edge=1 regardless of slot count.
    // bsr_mode=2: independent per-type multipliers — n_types_per_edge=slots.size().
    vector<BsrRateSlot> bsr_slots;
    if(opt_type.bsr_mode > 0){
        bsr_slots = get_bsr_rate_slots(rtree, cn_type);
        if(bsr_slots.empty())
            cerr << "Warning: bsr_mode=" << opt_type.bsr_mode
                 << " but all reference rates for cn_type=" << cn_type
                 << " are zero — BSR multipliers will have no effect. "
                 << "Check that non-zero rates are passed when using bsr_mode>0." << endl;
        // [2026-07-14 added] cache so update_variables_transformed/update_edge_rates_rlc_from_x
        // don't each recompute this on every BFGS objective/gradient evaluation — it's
        // invariant for the whole optimization (global reference rates are frozen).
        opt_type.bsr_slots_cache = bsr_slots;
    }

    int n_types_per_edge = 0;
    if(opt_type.bsr_mode == 1)      n_types_per_edge = 1;
    else if(opt_type.bsr_mode == 2) n_types_per_edge = (int)bsr_slots.size();
    else if(opt_type.bsr_mode == 3) n_types_per_edge = (int)bsr_slots.size();  

    // [2026-07-14 disabled] nrates (5th arg) was always bsr_slots.size(), reserving ndim
    // space for the 5 global reference rates even when bsr_mode>0 (where estmu is forced
    // 0 and nothing writes/reads those slots anymore) — wasted, always-
    // [0,0]-bounded optimizer dimensions. nrates is only meaningful when estmu=1 (global
    // rates genuinely being estimated, which only happens for bsr_mode=0). See
    // docs/flowcharts.md.
    // int ndim = get_ndim(estmu, nparams_est, model, cn_type, bsr_slots.size(), n_types_per_edge, n_bsr_edges);
    int ndim = get_ndim(estmu, nparams_est, model, cn_type, estmu ? bsr_slots.size() : 0, n_types_per_edge, n_bsr_edges);

    if(debug){
      cout << "\nThere are " << ndim << " parameters excluding global mutation rates to optimise " << endl;
    }

    double* variables = new double[ndim + 1];
    memset(variables, 0.0, (ndim + 1) * sizeof(double));
    double* upper_bound = new double[ndim + 1];
    memset(upper_bound, 0.0, (ndim + 1) * sizeof(double));
    double* lower_bound = new double[ndim + 1];
    memset(lower_bound, 0.0, (ndim + 1) * sizeof(double));

    if(cons){    // edges converted to ratio to incorporate time constraints
        // check tip validity
        if(!is_tip_age_valid(rtree.get_node_ages(), opt_type.tobs)){
          cout << "Tip timings inconsistent with observed data when doing BFGS!" << endl;
          rtree.print();
          cout << rtree.make_newick() << endl;
          exit(EXIT_FAILURE);
        }

        vector<double> ratios = rtree.get_ratio_from_age();

        if(debug){
            cout << "\nInitializing variables related to branch length" << endl;
            cout << "time ratios obtained from node ages: ";
            for(int i = 0; i < rtree.nleaf - 1; i++){
                cout << "\t" << ratios[i];
            }
            cout << endl;
        }

        if(opt_one_branch){   // only optimize one branch
            int nid = rtree.edges[rtree.current_eid].end - rtree.root_node_id;
            int idx = 1;
            variables[idx] = ratios[nid];
            lower_bound[idx] = MIN_RATIO;
            upper_bound[idx] = MAX_RATIO;
        }else{
            // age of root
            variables[1] = ratios[0];
            // make minimal age of root much larger than maximum observed time to avoid very small branch lengths and failing in optimization
            lower_bound[1] = lnl_type.max_tobs * opt_type.scale_tobs + rtree.nleaf * BLEN_MIN;
            // age at 1st sample, so need to add time until last sample
            upper_bound[1] = lnl_type.max_tobs + patient_age;

            for(int i = 1; i < nparams_est; ++i){
              variables[i + 1] = ratios[i];
              lower_bound[i + 1] = MIN_RATIO;
              upper_bound[i + 1] = MAX_RATIO;
            }
            
        }
    }else{
        if(opt_one_branch){
            variables[1] = rtree.edges[rtree.current_eid].length;
            lower_bound[1] = BLEN_MIN;
            upper_bound[1] = patient_age;
        }else{
            for(int i = 0; i < nparams_est; ++i){
              variables[i + 1] = rtree.edges[i].length;
              lower_bound[i + 1] = BLEN_MIN;
              upper_bound[i + 1] = patient_age;
            }
        }
    }

    // [2026-07-14 disabled] this pushed the fixed global reference rates (dup_rate etc.)
    // into the BFGS optimization variables for any bsr_mode>0, so they got estimated/
    // corrupted instead of staying frozen (they overlap with the per-edge multiplier
    // slots set up below). See docs/flowcharts.md "Two-step calibration workflow" for
    // the full writeup. Kept here commented out, not deleted, per team request.
    // // reinitialise edge_rates if missing or stale (size mismatch after tree changes)
    // if(opt_type.bsr_mode > 0){ //&& (int)rtree.edge_rates.size() != (int)rtree.edges.size()
    //     set_tree_rates(cn_type, rtree, nparams_est, variables, lower_bound, upper_bound);
    // }

    if(opt_type.bsr_mode == 1){
        // Shared bounds derived from the tightest constraint across all active slots.
        // double m_lower = 0.0, m_upper = std::numeric_limits<double>::infinity();
        // for(auto& slot : bsr_slots){
        //     m_lower = std::max(m_lower, MIN_MRATE / rtree.*(slot.tree_rate));
        //     m_upper = std::min(m_upper, MAX_MRATE / rtree.*(slot.tree_rate));
        // }
        // if(m_upper == std::numeric_limits<double>::infinity()) m_upper = MAX_MRATE;
        // if(m_lower <= 0.0) m_lower = MIN_MRATE;
        // if(m_lower > m_upper){
        //     cerr << "Warning: inconsistent multiplier bounds for cn_type=" << cn_type
        //          << "; resetting to [MIN_MRATE, MAX_MRATE]" << endl;
        //     m_lower = MIN_MRATE; m_upper = MAX_MRATE;
        // }

        for(int k = 0; k < n_bsr_edges; k++){
            int eid = active_eids[k];
            // warm-start: extract shared m from the first available slot
            double m = (!bsr_slots.empty())
                       ? rtree.edge_rates[eid].*(bsr_slots[0].edge_field) / rtree.*(bsr_slots[0].tree_rate)
                       : 1.0;
            int idx = nparams_est + k + 1;
            variables[idx]   = m;
            lower_bound[idx] = M_MIN;
            upper_bound[idx] = M_MAX;
        }
    } else if(opt_type.bsr_mode == 2 || opt_type.bsr_mode == 3){
        // Independent per-type multiplier per branch.
        // idx(k, t) = nparams_est + k * n_types_per_edge + t + 1
        for(int k = 0; k < n_bsr_edges; k++){
            int eid = active_eids[k];
            for(int t = 0; t < (int)bsr_slots.size(); t++){
                auto& slot = bsr_slots[t];
                double m = rtree.edge_rates[eid].*(slot.edge_field) / rtree.*(slot.tree_rate);
                // [2026-07-14 disabled] see bsr_var_index in optimization.hpp
                // int idx = nparams_est + k * n_types_per_edge + t + 1;
                int idx = bsr_var_index(nparams_est, k, n_types_per_edge, t);
                variables[idx]   = m;
                lower_bound[idx] = M_MIN;
                upper_bound[idx] = M_MAX;
            }
        }
    }
    // else if(opt_type.bsr_mode == 3){
    //     // RLC: one multiplier per shift edge. Use RLC_MULT bounds (not absolute rate bounds).
    //     // Warm-start at m=1.0 (no shift); stepwise outer loop provides a good starting point.
    //     for(int k = 0; k < n_bsr_edges; k++){
    //         int idx = nparams_est + k + 1;
    //         variables[idx]   = 1.0;
    //         lower_bound[idx] = RLC_MULT_MIN;
    //         upper_bound[idx] = RLC_MULT_MAX;
    //     }
    // } 
    else if(estmu){    // estimate mutation rates, bsr_mode=0
        if(model == MK){
            int i = nparams_est;
            variables[i + 1] = rtree.mu;
            lower_bound[i + 1] = MIN_MRATE;
            upper_bound[i + 1] = MAX_MRATE;
        }else{ // other models
            set_tree_rates(cn_type, rtree, nparams_est, variables, lower_bound, upper_bound, opt_type.bsr_mode);
       }
    }

    if(debug){
        if(opt_one_branch){
            cout << "only optimize one branch" << endl;
        }
        cout << "Variables to estimate: " << endl;
        for(int i = 0; i < ndim + 1; i++){
            cout << "\t" << variables[i] << "\tLB:" << lower_bound[i] << "\tUB:" << upper_bound[i] << endl;
        }
    }

    // variables contains the parameters to estimate(branch length, mutation rate)
    min_nlnl = L_BFGS_B(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, ndim, variables + 1, lower_bound + 1, upper_bound + 1);

    // ensure rtree (branch lengths, edge_rates) reflects the final optimised variables,
    // since the last internal gradient call may have left rtree at a perturbed point
    update_variables_transformed(rtree, variables, lnl_type, opt_type);

    // Check the validity of the tree
    if(cons && !is_tree_valid(rtree, lnl_type.max_tobs, patient_age, cons)){
      cout << "The optimized tree after BFGS " << rtree.make_newick() << " is not valid!" << endl;
      exit(EXIT_FAILURE);
    }

    if (debug)
    {
        cout << "log likelihood of current ML tree " << rtree.make_newick() << " is " << -min_nlnl << endl;
    }

    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;
}



// Optimizing the branch length of(node1, node2) with BFGS to incorporate constraints imposed by patient age and tip timings.
// When any branch length is updated, the neighbour length needs to be updated
// Optimize mutation rates if necessary
double optimize_one_branch_BFGS(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, Node* node1, Node* node2){
    int debug = 0;
    if(debug){
        cout << "\tOptimizing the branch " << node1->id + 1 << ", " << node2->id + 1 << endl;
    }

    // does not optimize virtual branch from root
    assert(!((node1->id == rtree.nleaf && node2->id == rtree.nleaf - 1) || (node2->id == rtree.nleaf && node1->id == rtree.nleaf - 1)));

    int estmu = opt_type.estmu;    // record orignal estmu

    rtree.current_it = (Neighbor*) node1->findNeighbor(node2);
    assert(rtree.current_it);
    rtree.current_it_back = (Neighbor*) node2->findNeighbor(node1);
    assert(rtree.current_it_back);
    assert(rtree.current_it->length == rtree.current_it_back->length);

    int eid = rtree.get_edge_id(node1->id, node2->id);
    assert(eid >= 0);
    rtree.current_eid = eid;  // used to track the branch to optimize

    if(debug){
        cout << "tree before optimization by BFGS " << rtree.make_newick() << endl;
        rtree.print();
        cout << "\tUsing BFGS method to optimize the likelihood of one branch length of edge " << eid + 1 << endl;
        // cout << "\taddress of rtree before optimization " << &rtree << endl;
        // cout << "\taddress of node1 before optimization " << &(rtree.node1) << endl;
        // cout << "\taddress of current_it before optimization " << rtree.current_it << endl;
    }

    opt_type.opt_one_branch = 1;
    opt_type.estmu = 0;
    // bsr_mode is intentionally NOT set to 0: when opt_one_branch=1, active_eids is empty
    // so no multipliers enter the optimisation, but the variable-rate likelihood is still used
    // with the current (frozen) edge_rates. This is the correct behaviour.
    double negative_lh = MAX_NLNL;
    // optimize ratio based on NNI branch(branch length, node times and ages have been updated during optimization)
    max_likelihood_BFGS(rtree, vobs, vobs_change, obs_decomp, comps, lnl_type, opt_type, negative_lh);

    // cout << "\taddress of rtree after optimization " << &rtree << endl;
    // cout << "\taddress of current_it after optimization " << rtree.current_it << endl;

    opt_type.opt_one_branch = 0;
    opt_type.estmu = estmu;

    if(debug){
        cout << "\tmax logl: " << -negative_lh << " optimized branch length " << rtree.edges[eid].length << endl;
        cout << "tree after optimization by BFGS " << rtree.make_newick() << endl;
        rtree.print();
    }

    return -negative_lh;
}


/****** bsr_mode=3 ML-RLC implementation ******/

void update_edge_rates_rlc(
    evo_tree& rtree,
    const vector<int>& shift_eids,
    const vector<RateSet>& multipliers)
{
    assert(shift_eids.size() == multipliers.size());

    // resize edge_rates if needed
    if((int)rtree.edge_rates.size() != (int)rtree.edges.size())
        rtree.edge_rates.resize(rtree.edges.size());

    // build eid -> multiplier lookup
    unordered_map<int,RateSet> shift_m;
    for(int k = 0; k < (int)shift_eids.size(); ++k)
        shift_m[shift_eids[k]] = multipliers[k];

    // per-node current local RateSet; initialise all to global rates
    RateSet global_rates(0, rtree.dup_rate, rtree.del_rate,
                         rtree.chr_gain_rate, rtree.chr_loss_rate, rtree.wgd_rate);
    vector<RateSet> node_rate(rtree.nodes.size(), global_rates);

    // preorder traversal from root
    stack<int> stk;
    stk.push(rtree.root_node_id);
    while(!stk.empty()){
        int nid = stk.top(); stk.pop();
        for(int child : rtree.nodes[nid].daughters){
            int eid = rtree.nodes[child].e_in;
            RateSet parent_rate = node_rate[nid];
            RateSet child_rate;
            auto it = shift_m.find(eid);
            if(it != shift_m.end()){
                const RateSet& m = it->second;
                child_rate = parent_rate * m;
            } else {
                child_rate = parent_rate;
            }
            rtree.edge_rates[eid] = child_rate;
            node_rate[child]       = child_rate;
            stk.push(child);
        }
    }
}


void update_edge_rates_rlc_from_x(
    evo_tree& rtree,
    const vector<int>& shift_eids,
    double* x,
    int nparams_est,
    const vector<BsrRateSlot>& bsr_slots)
{
    // [2026-07-14 disabled] recomputed bsr_slots on every call (fires on every BFGS
    // objective/gradient evaluation for bsr_mode=3); it's invariant for the whole
    // optimization, so the caller now passes max_likelihood_BFGS's cached
    // opt_type.bsr_slots_cache instead. See docs/flowcharts.md.
    // vector<BsrRateSlot> bsr_slots = get_bsr_rate_slots(rtree, cn_type);
    int n_types_per_edge = (int)bsr_slots.size();

    // neutral (1.0) multiplier for every field; only the active bsr_slots get overwritten below
    vector<RateSet> multipliers(shift_eids.size(), RateSet(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
    for(int k = 0; k < (int)shift_eids.size(); ++k){
        for(int t = 0; t < n_types_per_edge; ++t){
            // [2026-07-14 disabled] see bsr_var_index in optimization.hpp
            // int idx = nparams_est + k * n_types_per_edge + t + 1;
            int idx = bsr_var_index(nparams_est, k, n_types_per_edge, t);
            multipliers[k].*(bsr_slots[t].edge_field) = x[idx];
        }
    }
    update_edge_rates_rlc(rtree, shift_eids, multipliers);
}


void stepwise_search_shift_edges(
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug)
{
    // debug = 1;
    assert(opt_type.bsr_mode == 3);
    assert(opt_type.estmu == 0);

    const int criterion       = opt_type.rlc_criterion;  // 0=AIC, 1=BIC (see compute_rlc_ic)
    const int n_sites         = lnl_type.n_sites_for_ic;
    const double eps          = 1e-6;  // small improvement threshold
    // [2026-07-14 added] number of independent parameters a shift edge unlocks (5 for
    // cn_type=ALL); the penalty must scale by this, not just by the number of edges.
    // See compute_rlc_ic in optimization.hpp.
    const int n_types_per_edge = (int)get_bsr_rate_slots(rtree, lnl_type.cn_type).size();

    // full candidate pool = all optimizable edges (excludes normal-sample edge)
    vector<int> all_candidates = get_active_bsr_eids(rtree);

    // --- K=0 baseline ---
    OPT_TYPE cur_opt = opt_type;
    cur_opt.rlc_shift_eids = {};
    evo_tree cur_tree = rtree;

    double cur_nlnl = MAX_NLNL;
    max_likelihood_BFGS(cur_tree, vobs, vobs_change, obs_decomp, comps,
                        lnl_type, cur_opt, cur_nlnl, debug);
    RlcIC cur_ic = compute_rlc_ic(cur_nlnl, 0, n_types_per_edge, n_sites, criterion);
    double cur_score = cur_ic.score;

    if(debug)
        cout << "[RLC] K=0 baseline: nlnl=" << cur_nlnl
             << " raw_logL=" << cur_ic.raw_logL << " AIC=" << cur_ic.aic << " BIC=" << cur_ic.bic
             << " score=" << cur_score << endl;

    // ----- start here: original code ----- //
    // // --- forward stepwise ---
    // while(true){
    //     double best_score = cur_score;
    //     double best_nlnl  = cur_nlnl;
    //     int    best_edge  = -1;
    //     evo_tree best_tree = cur_tree;
    //     OPT_TYPE best_opt  = cur_opt;

    //     for(int eid : all_candidates){
    //         // skip already-selected shift edges
    //         bool already = false;
    //         for(int s : cur_opt.rlc_shift_eids)
    //             if(s == eid){ already = true; break; }
    //         if(already) continue;

    //         // build candidate shift set
    //         vector<int> cand_shifts = cur_opt.rlc_shift_eids;
    //         cand_shifts.push_back(eid);

    //         // start from the current accepted tree (copy)
    //         evo_tree cand_tree = cur_tree;
    //         OPT_TYPE cand_opt  = cur_opt;
    //         cand_opt.rlc_shift_eids = cand_shifts;

    //         // warm-start: new shift edge gets m=1.0 (set via init)
    //         // update_edge_rates_rlc with uniform multipliers=1 keeps rates consistent
    //         vector<double> ones(cand_shifts.size(), 1.0);
    //         update_edge_rates_rlc(cand_tree, cand_shifts, ones);

    //         double cand_nlnl = MAX_NLNL;
    //         max_likelihood_BFGS(cand_tree, vobs, vobs_change, obs_decomp, comps,
    //                             lnl_type, cand_opt, cand_nlnl, debug);

    //         double cand_score = -cand_nlnl - lambda * (double)cand_shifts.size();

    //         if(debug)
    //             cout << "[RLC] candidate eid=" << eid
    //                  << " K=" << cand_shifts.size()
    //                  << " nlnl=" << cand_nlnl
    //                  << " score=" << cand_score << endl;

    //         if(cand_score > best_score + eps){
    //             best_score = cand_score;
    //             best_nlnl  = cand_nlnl;
    //             best_edge  = eid;
    //             best_tree  = cand_tree;
    //             best_opt   = cand_opt;
    //         }
    //     }

    //     if(best_edge == -1) break;  // no improvement found

    //     cur_score = best_score;
    //     cur_nlnl  = best_nlnl;
    //     cur_tree  = best_tree;
    //     cur_opt   = best_opt;

    //     cout << "[RLC] accepted shift edge " << best_edge
    //          << " K=" << cur_opt.rlc_shift_eids.size()
    //          << " nlnl=" << cur_nlnl
    //          << " penalized_score=" << cur_score << endl;

    //     // --- optional backward cleanup ---
    //     // bool removed = true;
    //     // while(removed){
    //     //     removed = false;
    //     //     for(int k = 0; k < (int)cur_opt.rlc_shift_eids.size(); ++k){
    //     //         vector<int> try_shifts = cur_opt.rlc_shift_eids;
    //     //         try_shifts.erase(try_shifts.begin() + k);

    //     //         evo_tree try_tree = cur_tree;
    //     //         OPT_TYPE try_opt  = cur_opt;
    //     //         try_opt.rlc_shift_eids = try_shifts;

    //     //         vector<double> ones(try_shifts.size(), 1.0);
    //     //         update_edge_rates_rlc(try_tree, try_shifts, ones);

    //     //         double try_nlnl = MAX_NLNL;
    //     //         max_likelihood_BFGS(try_tree, vobs, vobs_change, obs_decomp, comps,
    //     //                             lnl_type, try_opt, try_nlnl, debug);

    //     //         double try_score = -try_nlnl - lambda * (double)try_shifts.size();
    //     //         if(try_score > cur_score + eps){
    //     //             int removed_eid = cur_opt.rlc_shift_eids[k]; // save before cur_opt is replaced
    //     //             cur_score = try_score;
    //     //             cur_nlnl  = try_nlnl;
    //     //             cur_tree  = try_tree;
    //     //             cur_opt   = try_opt;
    //     //             removed   = true;
    //     //             cout << "[RLC] backward: removed shift edge " << removed_eid
    //     //                  << " K=" << cur_opt.rlc_shift_eids.size()
    //     //                  << " score=" << cur_score << endl;
    //     //             break;
    //     //         }
    //     //     }
    //     // }
    
    // }
    // ----- end here: original code ----- //
    
    
    // ----- start here: if using OpenMP, parallelize candidate calculation only ----- //
    while(true){
        int ncand = (int)all_candidates.size();

        // Store each candidate result at its corresponding index.
        vector<char> evaluated(ncand, 0);
        vector<double> candidate_scores(ncand, -MAX_NLNL);
        vector<double> candidate_nlnls(ncand, MAX_NLNL);
        vector<double> candidate_aics(ncand, MAX_NLNL);
        vector<double> candidate_bics(ncand, MAX_NLNL);
        vector<evo_tree> candidate_trees(ncand, cur_tree);
        vector<OPT_TYPE> candidate_opts(ncand, cur_opt);

        // ============================================================
        // FORWARD: parallel candidate evaluation
        // ============================================================
        // Parallel section: only evaluate candidates; do not update the global best result here.
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < ncand; ++i){
            int eid = all_candidates[i];

            // Skip edges that have already been selected as shift edges.
            bool already = false;
            for(int s : cur_opt.rlc_shift_eids){
                if(s == eid){
                    already = true;
                    break;
                }
            }
            if(already) continue;

            // Build the candidate shift-edge set.
            vector<int> cand_shifts = cur_opt.rlc_shift_eids;
            cand_shifts.push_back(eid);

            // Each candidate uses its own copies.
            evo_tree cand_tree = cur_tree;
            OPT_TYPE cand_opt  = cur_opt;
            LNL_TYPE cand_lnl_type = lnl_type;

            cand_opt.rlc_shift_eids = cand_shifts;

            // Warm start: initialize the new shift multiplier at 1.0 (all rate types).
            vector<RateSet> ones(cand_shifts.size(), RateSet(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
            update_edge_rates_rlc(cand_tree, cand_shifts, ones);

            double cand_nlnl = MAX_NLNL;

            // Keep debug off inside the parallel region to avoid interleaved output.
            max_likelihood_BFGS(cand_tree, vobs, vobs_change, obs_decomp, comps,
                                cand_lnl_type, cand_opt, cand_nlnl, 0);

            RlcIC cand_ic = compute_rlc_ic(cand_nlnl, (int)cand_shifts.size(), n_types_per_edge, n_sites, criterion);
            double cand_score = cand_ic.score;

            // Each thread writes only to its own index.
            evaluated[i]        = 1;
            candidate_scores[i] = cand_score;
            candidate_nlnls[i]  = cand_nlnl;
            candidate_aics[i]   = cand_ic.aic;
            candidate_bics[i]   = cand_ic.bic;
            candidate_trees[i]  = std::move(cand_tree);
            candidate_opts[i]   = std::move(cand_opt);
        }

        // ============================================================
        // FORWARD: serial selection
        // ============================================================
        // Serial section: compare all evaluated candidates and select the best one.
        // Print per-candidate results here (outside the parallel region) so output
        // isn't interleaved across threads. Always on (independent of the general
        // debug flag) since this is the main RLC search diagnostic.
        for(int i = 0; i < ncand; ++i){
            if(!evaluated[i]) continue;
            cout << "[RLC] candidate eid=" << all_candidates[i]
                 << " K=" << cur_opt.rlc_shift_eids.size() + 1
                 << " nlnl=" << candidate_nlnls[i]
                 << " AIC=" << candidate_aics[i]
                 << " BIC=" << candidate_bics[i]
                 << " score=" << candidate_scores[i] << endl;
        }

        // [2026-07-14 disabled] see pick_best_improving in optimization.hpp — this ~25-line
        // loop was hand-copied here and in the backward removal selection below.
        // int best_i = -1;
        // for(int i = 0; i < ncand; ++i){
        //     if(!evaluated[i]) continue;
        //     double cand_score = candidate_scores[i];
        //     int eid = all_candidates[i];
        //     // A candidate must improve over the current accepted model.
        //     if(cand_score <= cur_score + eps) continue;
        //     if(best_i == -1){
        //         best_i = i;
        //         continue;
        //     }
        //     double best_candidate_score = candidate_scores[best_i];
        //     int best_candidate_eid = all_candidates[best_i];
        //     bool better =
        //         cand_score > best_candidate_score + eps;
        //     bool tie_but_smaller =
        //         fabs(cand_score - best_candidate_score) <= eps &&
        //         eid < best_candidate_eid;
        //     if(better || tie_but_smaller){
        //         best_i = i;
        //     }
        // }
        int best_i = pick_best_improving(ncand, evaluated, candidate_scores, all_candidates, cur_score, eps);

        // Stop the stepwise search if no candidate improves the current model.
        if(best_i == -1) break;

        int best_edge = all_candidates[best_i];

        cur_score = candidate_scores[best_i];
        cur_nlnl  = candidate_nlnls[best_i];
        cur_tree  = std::move(candidate_trees[best_i]);
        cur_opt   = std::move(candidate_opts[best_i]);

        cout << "[RLC] accepted shift edge " << best_edge
            << " K=" << cur_opt.rlc_shift_eids.size()
            << " nlnl=" << cur_nlnl
            << " penalized_score=" << cur_score << endl;

        // ============================================================
        // BACKWARD CLEANUP
        // After accepting one forward edge, repeatedly test whether
        // removing one currently selected edge improves the score.
        // ============================================================
        bool removed = true;
        while(removed){
            removed = false;
            int nshift = (int)cur_opt.rlc_shift_eids.size();
            if(nshift == 0) break;
            // save current shifted-edge IDs because cur_opt will be modified during the loop
            vector<int> current_shifts = cur_opt.rlc_shift_eids;

            //store one result for remoing each current shift edge
            vector<char> remove_evaluated(nshift, 0);
            vector<double> remove_scores(nshift, -MAX_NLNL);
            vector<double> remove_nlnls(nshift, MAX_NLNL);
            vector<double> remove_aics(nshift, MAX_NLNL);
            vector<double> remove_bics(nshift, MAX_NLNL);
            vector<evo_tree> remove_trees(nshift, cur_tree);
            vector<OPT_TYPE> remove_opts(nshift, cur_opt);

        // ========================================================
        // BACKWARD: parallel candidate evaluation
        // ========================================================
        #pragma omp parallel for schedule(dynamic)
        for(int k = 0; k < nshift; ++k){
            vector<int> try_shifts = current_shifts;

            // Candidate k removes current_shifts[k].
            try_shifts.erase(try_shifts.begin() + k);

            evo_tree try_tree = cur_tree;
            OPT_TYPE try_opt = cur_opt;
            LNL_TYPE try_lnl_type = lnl_type;

            try_opt.rlc_shift_eids = try_shifts;

            vector<RateSet> ones(try_shifts.size(), RateSet(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
            update_edge_rates_rlc(try_tree, try_shifts, ones);

            double try_nlnl = MAX_NLNL;

            // Keep debug off inside the parallel region to avoid interleaved output.
            max_likelihood_BFGS(try_tree, vobs, vobs_change, obs_decomp, comps,
                                try_lnl_type, try_opt, try_nlnl, 0);

            RlcIC try_ic = compute_rlc_ic(try_nlnl, (int)try_shifts.size(), n_types_per_edge, n_sites, criterion);
            double try_score = try_ic.score;

            remove_evaluated[k] = 1;
            remove_scores[k]    = try_score;
            remove_nlnls[k]     = try_nlnl;
            remove_aics[k]      = try_ic.aic;
            remove_bics[k]      = try_ic.bic;
            remove_trees[k]     = std::move(try_tree);
            remove_opts[k]      = std::move(try_opt);
        }

        // ========================================================
        // BACKWARD: serial selection
        // ========================================================
        // Print per-candidate removal results here (outside the parallel region) so
        // output isn't interleaved across threads. Always on (independent of the
        // general debug flag) since this is the main RLC search diagnostic.
        for(int k = 0; k < nshift; ++k){
            if(!remove_evaluated[k]) continue;
            cout << "[RLC] backward candidate remove_eid=" << current_shifts[k]
                 << " K=" << nshift - 1
                 << " nlnl=" << remove_nlnls[k]
                 << " AIC=" << remove_aics[k]
                 << " BIC=" << remove_bics[k]
                 << " score=" << remove_scores[k] << endl;
        }

        // [2026-07-14 disabled] see pick_best_improving in optimization.hpp — same
        // selection logic as the forward candidate selection above.
        // int best_remove_k = -1;
        // for(int k = 0; k < nshift; ++k){
        //     if(!remove_evaluated[k]) continue;
        //     double try_score = remove_scores[k];
        //     int removed_eid = current_shifts[k];
        //     // Removal must improve the current accepted model.
        //     if(try_score <= cur_score + eps) continue;
        //     if(best_remove_k == -1){
        //         best_remove_k = k;
        //         continue;
        //     }
        //     double best_remove_score = remove_scores[best_remove_k];
        //     int best_removed_eid = current_shifts[best_remove_k];
        //     bool better = try_score > best_remove_score + eps;
        //     bool tie_but_smaller = fabs(try_score - best_remove_score) <= eps && removed_eid < best_removed_eid;
        //     if(better || tie_but_smaller){
        //         best_remove_k = k;
        //     }
        // }
        int best_remove_k = pick_best_improving(nshift, remove_evaluated, remove_scores, current_shifts, cur_score, eps);

        // No backward deletion improves the model.
        if(best_remove_k == -1) break;

        int removed_eid = current_shifts[best_remove_k];

        cur_score = remove_scores[best_remove_k];
        cur_nlnl = remove_nlnls[best_remove_k];
        cur_tree = std::move(remove_trees[best_remove_k]);
        cur_opt = std::move(remove_opts[best_remove_k]);

        removed = true;

        cout << "[RLC] backward: removed shift edge "
             << removed_eid
             << " K=" << cur_opt.rlc_shift_eids.size()
             << " nlnl=" << cur_nlnl
             << " penalized_score=" << cur_score
             << endl;
    }
        
}
        
        


// ----- end here: if using OpenMP, parallelize candidate calculation only ----- //




    // write back final results
    // For bsr_mode=3, tree comparisons use the AIC/BIC score from compute_rlc_ic (per
    // opt_type.rlc_criterion). rtree.score holds that score (larger = better).
    // min_nlnl holds -(score) so minimizing min_nlnl = maximizing score.
    rtree                   = cur_tree;
    opt_type.rlc_shift_eids = cur_opt.rlc_shift_eids;

    int K = (int)opt_type.rlc_shift_eids.size();
    RlcIC final_ic = compute_rlc_ic(cur_nlnl, K, n_types_per_edge, n_sites, criterion);
    double raw_logL        = final_ic.raw_logL;
    double penalized_score = final_ic.score;

    opt_type.rlc_raw_logL = raw_logL;
    opt_type.rlc_penalized_score = penalized_score;
    rtree.rlc_shift_eids = opt_type.rlc_shift_eids;
    rtree.rlc_raw_logL = raw_logL;
    rtree.rlc_penalized_score = penalized_score;
    rtree.score = penalized_score;
    min_nlnl    = -penalized_score;   // negative penalized score, minimized by caller

    cout << "[RLC] final: K=" << K << " shift_eids=[ ";
    for(int s : opt_type.rlc_shift_eids) cout << s << " ";
    cout << "] raw_logL=" << raw_logL
         << " AIC=" << final_ic.aic << " BIC=" << final_ic.bic
         << " score=" << penalized_score << endl;
}


// One candidate solution: which candidate edges (indices into all_candidates,
// not raw eids) are shift edges. vector<char> rather than vector<bool> because
// population evaluation will run under #pragma omp parallel for (like
// stepwise_search_shift_edges), and vector<bool>'s bit-packed storage is not
// safe to write concurrently at different indices.
struct GaIndividual {
    vector<char> genes;              // genes[i]=1 means all_candidates[i] is a shift edge
    double nlnl  = MAX_NLNL;
    double score = -MAX_NLNL;
    double aic   = MAX_NLNL;  // [2026-07-21 added] reported alongside score/bic, not used for selection
    double bic   = MAX_NLNL;  // [2026-07-21 added] reported alongside score/aic, not used for selection
};

// Genetic-algorithm search for shift edges under bsr_mode=3 (ML-RLC). Chromosome =
// one bit per candidate edge (1 = shift edge). Roulette-wheel selection with elitism,
// uniform crossover, bit-flip mutation; terminates on max generations or a no-improvement
// streak. Not yet validated against stepwise_search_shift_edges on real data.
void genetic_search_shift_edges(
    gsl_rng* r,
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug)
{
    vector<int> all_candidates = get_active_bsr_eids(rtree);
    int ncand = (int)all_candidates.size();

    const int criterion = opt_type.rlc_criterion;  // 0=AIC, 1=BIC (see compute_rlc_ic)
    const int n_sites = lnl_type.n_sites_for_ic;
    const int n_types_per_edge = (int)get_bsr_rate_slots(rtree, lnl_type.cn_type).size();

    // Cache keyed by the raw gene pattern (e.g. "0010100..."), shared across the whole
    // search: the same shift-edge combination can easily reappear across generations
    // (elitism keeps it around, crossover/mutation can regenerate it by chance), and each
    // miss costs a full BFGS run. Guarded by #pragma omp critical since evaluate() below
    // runs inside #pragma omp parallel for.
    struct GaCacheEntry { double nlnl, score, aic, bic; };
    unordered_map<string, GaCacheEntry> score_cache;  // gene pattern -> cached result

    // Fitness of one individual: decode its genes into a shift-edge set, warm-start
    // that set's multipliers at 1.0 (same convention as stepwise_search_shift_edges),
    // run BFGS from a fresh copy of the input rtree (not an evolving "current best"
    // tree — unlike the stepwise search, each individual is an independent candidate
    // solution, not an incremental extension of a shared one), and score it.
    auto evaluate = [&](GaIndividual& indiv){
        string key(indiv.genes.begin(), indiv.genes.end());

        bool cache_hit = false;
        #pragma omp critical(ga_score_cache)
        {
            auto it = score_cache.find(key);
            if(it != score_cache.end()){
                indiv.nlnl  = it->second.nlnl;
                indiv.score = it->second.score;
                indiv.aic   = it->second.aic;
                indiv.bic   = it->second.bic;
                cache_hit = true;
            }
        }
        if(cache_hit) return;

        vector<int> shift_eids;
        for(int i = 0; i < ncand; ++i){
            if(indiv.genes[i]) shift_eids.push_back(all_candidates[i]);
        }

        evo_tree cand_tree = rtree;
        OPT_TYPE cand_opt  = opt_type;
        LNL_TYPE cand_lnl_type = lnl_type;
        cand_opt.rlc_shift_eids = shift_eids;

        vector<RateSet> ones(shift_eids.size(), RateSet(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
        update_edge_rates_rlc(cand_tree, shift_eids, ones);

        double cand_nlnl = MAX_NLNL;
        // debug off inside the parallel region below, same reasoning as stepwise_search_shift_edges
        max_likelihood_BFGS(cand_tree, vobs, vobs_change, obs_decomp, comps,
                            cand_lnl_type, cand_opt, cand_nlnl, 0);

        RlcIC cand_ic = compute_rlc_ic(cand_nlnl, (int)shift_eids.size(), n_types_per_edge, n_sites, criterion);
        indiv.nlnl  = cand_nlnl;
        indiv.score = cand_ic.score;
        indiv.aic   = cand_ic.aic;
        indiv.bic   = cand_ic.bic;

        #pragma omp critical(ga_score_cache)
        {
            score_cache[key] = {indiv.nlnl, indiv.score, indiv.aic, indiv.bic};
        }
    };

    const int pop_size = 20;
    const double p_init = 0.1;

    vector<GaIndividual> population(pop_size);
    for(auto& indiv : population){
        indiv.genes.assign(ncand, 0);
        for(int i = 0; i < ncand; ++i){
            if(gsl_rng_uniform(r) < p_init) indiv.genes[i] = 1;
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for(int p = 0; p < pop_size; ++p){
        evaluate(population[p]);
    }

    if(debug){
        for(int p = 0; p < pop_size; ++p){
            int k = 0;
            cout << "[GA init] indiv " << p << ": shift_eids=[ ";
            for(int i = 0; i < ncand; ++i){
                if(population[p].genes[i]){ cout << all_candidates[i] << " "; k++; }
            }
            cout << "] K=" << k
                 << " nlnl=" << population[p].nlnl
                 << " score=" << population[p].score << endl;
        }
    }

    // Elites: the nelites best-scoring individuals bypass crossover/mutation entirely
    // and are copied unchanged into the next generation, so the best solution found so
    // far can never be lost to the stochastic steps below.
    const int nelites = std::min(2, pop_size);
    // Expected ~1 gene flips per individual, a standard GA default; guarded against
    // ncand==0 (no candidate edges) to avoid a division by zero.
    const double mutation_rate = (ncand > 0) ? (1.0 / ncand) : 0.0;
    const double mate_prob     = 0.5;  // uniform crossover mixing probability
    const int max_generation   = 30;
    const int max_no_improve   = 10;
    const double eps           = 1e-6;

    auto by_score_desc = [](const GaIndividual& a, const GaIndividual& b){
        return a.score > b.score;
    };

    // Roulette-wheel selection: builds a pop_size-long pool of parents to draw pairs from
    // for crossover. Scores are shifted so the worst individual in pop gets a weight just
    // above 0 (same rescaling idea as ga_roulette_wheel_selection in ga.c, which shifts by
    // the minimum fitness when fitness can be negative — compute_rlc_ic's score here is
    // always deeply negative, being the negated AIC/BIC value). The pool itself is not the
    // next generation: elites (above) bypass it entirely.
    auto select_parent_pool = [&](const vector<GaIndividual>& pop) -> vector<GaIndividual> {
        int n = (int)pop.size();
        double min_score = pop[0].score;
        for(const auto& indiv : pop) min_score = std::min(min_score, indiv.score);

        vector<double> weights(n);
        double total_weight = 0;
        for(int i = 0; i < n; ++i){
            weights[i] = (pop[i].score - min_score) + 1e-6;
            total_weight += weights[i];
        }

        vector<GaIndividual> pool;
        pool.reserve(n);
        for(int k = 0; k < n; ++k){
            double target = gsl_rng_uniform(r) * total_weight;
            double acc = 0;
            int chosen = n - 1;
            for(int i = 0; i < n; ++i){
                acc += weights[i];
                if(acc >= target){ chosen = i; break; }
            }
            pool.push_back(pop[chosen]);
        }
        return pool;
    };

    // Uniform crossover: child starts as a copy of parent 1, then each gene independently
    // has probability mate_prob of being overwritten by parent 2's corresponding gene.
    auto crossover = [&](const GaIndividual& p1, const GaIndividual& p2) -> GaIndividual {
        GaIndividual child;
        child.genes = p1.genes;
        for(int i = 0; i < ncand; ++i){
            if(gsl_rng_uniform(r) < mate_prob) child.genes[i] = p2.genes[i];
        }
        return child;
    };

    // Bit-flip mutation: each gene independently has probability mutation_rate of flipping.
    auto mutate = [&](GaIndividual& indiv){
        for(int i = 0; i < ncand; ++i){
            if(gsl_rng_uniform(r) < mutation_rate) indiv.genes[i] = 1 - indiv.genes[i];
        }
    };

    double best_score_so_far = -MAX_NLNL;
    int n_no_improve = 0;

    for(int gen = 0; gen < max_generation; ++gen){
        std::sort(population.begin(), population.end(), by_score_desc);

        double gen_best = population.front().score;
        if(gen_best > best_score_so_far + eps){
            best_score_so_far = gen_best;
            n_no_improve = 0;
        } else {
            n_no_improve++;
        }

        cout << "[GA] generation " << gen
             << " best_nlnl=" << population.front().nlnl
             << " best_score=" << gen_best
             << " n_no_improve=" << n_no_improve << endl;

        // Stop once the population hasn't improved for max_no_improve generations, or
        // once max_generation is reached — same two termination conditions as
        // ga_default_termination in ga.c (max generations / no-improvement streak).
        if(n_no_improve >= max_no_improve) break;
        if(gen == max_generation - 1) break;

        vector<GaIndividual> parents = select_parent_pool(population);

        // Elites carry over unchanged; the remaining slots are filled by crossover+mutation.
        vector<GaIndividual> next_population(population.begin(), population.begin() + nelites);
        while((int)next_population.size() < pop_size){
            int i1 = gsl_rng_uniform_int(r, pop_size);
            int i2 = gsl_rng_uniform_int(r, pop_size);
            GaIndividual child = crossover(parents[i1], parents[i2]);
            mutate(child);
            next_population.push_back(child);
        }

        // Only the newly created offspring need evaluating — elites already carry a
        // valid nlnl/score from the previous generation (or from cache, if this exact
        // gene pattern was seen before).
        #pragma omp parallel for schedule(dynamic)
        for(int p = nelites; p < pop_size; ++p){
            evaluate(next_population[p]);
        }

        population = std::move(next_population);
    }

    // ------------------------------------------------------------------
    // Write back final results, same convention as stepwise_search_shift_edges:
    // rtree.score holds the penalized score (larger = better); min_nlnl holds
    // -(penalized_score) so minimizing min_nlnl = maximizing penalized score.
    // Unlike stepwise (which keeps its best tree copy, cur_tree, updated throughout the
    // search), GaIndividual doesn't carry an evo_tree — keeping a full tree copy per
    // individual per generation would be wasteful, and isn't needed for anything but the
    // eventual winner. So the winning gene pattern is decoded once more here and its tree
    // rematerialized with one extra BFGS run (negligible next to the pop_size*generations
    // runs already spent finding it).
    // ------------------------------------------------------------------
    std::sort(population.begin(), population.end(), by_score_desc);
    const GaIndividual& best = population.front();

    vector<int> best_shift_eids;
    for(int i = 0; i < ncand; ++i){
        if(best.genes[i]) best_shift_eids.push_back(all_candidates[i]);
    }

    opt_type.rlc_shift_eids = best_shift_eids;

    vector<RateSet> ones(best_shift_eids.size(), RateSet(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
    update_edge_rates_rlc(rtree, best_shift_eids, ones);

    double final_nlnl = MAX_NLNL;
    max_likelihood_BFGS(rtree, vobs, vobs_change, obs_decomp, comps,
                        lnl_type, opt_type, final_nlnl, debug);

    int K = (int)best_shift_eids.size();
    RlcIC final_ic = compute_rlc_ic(final_nlnl, K, n_types_per_edge, n_sites, criterion);
    double raw_logL        = final_ic.raw_logL;
    double penalized_score = final_ic.score;

    opt_type.rlc_raw_logL        = raw_logL;
    opt_type.rlc_penalized_score = penalized_score;
    rtree.rlc_shift_eids         = opt_type.rlc_shift_eids;
    rtree.rlc_raw_logL           = raw_logL;
    rtree.rlc_penalized_score    = penalized_score;
    rtree.score                  = penalized_score;
    min_nlnl                     = -penalized_score;

    cout << "[GA] final: K=" << K << " shift_eids=[ ";
    for(int s : opt_type.rlc_shift_eids) cout << s << " ";
    cout << "] raw_logL=" << raw_logL
         << " AIC=" << final_ic.aic << " BIC=" << final_ic.bic
         << " score=" << penalized_score << endl;
}


// [2026-07-21 added] Brute-force baseline for bsr_mode=3 ML-RLC: evaluates every subset of
// get_active_bsr_eids's candidate edges (2^ncand subsets) instead of searching greedily
// (stepwise) or heuristically (GA). Only feasible for small ncand (9 candidates on a
// 5-sample tree -> 512 subsets). Because every subset is evaluated, the AIC-best and
// BIC-best subsets are both known from one run at no extra BFGS cost — reported alongside
// each other regardless of opt_type.rlc_criterion, which only decides which one is written
// back into rtree/opt_type.
void exhaustive_search_shift_edges(
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug)
{
    assert(opt_type.bsr_mode == 3);
    assert(opt_type.estmu == 0);

    const int criterion = opt_type.rlc_criterion;  // 0=AIC, 1=BIC (see compute_rlc_ic)
    const int n_sites = lnl_type.n_sites_for_ic;
    const int n_types_per_edge = (int)get_bsr_rate_slots(rtree, lnl_type.cn_type).size();

    vector<int> all_candidates = get_active_bsr_eids(rtree);
    int ncand = (int)all_candidates.size();
    long long n_subsets = 1LL << ncand;

    if(ncand > 20){
        cerr << "[RLC][exhaustive] warning: ncand=" << ncand << " candidate edges -> "
             << n_subsets << " subsets to evaluate; this brute-force search is meant as a "
                "small-tree validation baseline, not routine use on larger trees." << endl;
    }

    vector<double> subset_nlnls(n_subsets, MAX_NLNL);
    vector<double> subset_aics(n_subsets, MAX_NLNL);
    vector<double> subset_bics(n_subsets, MAX_NLNL);
    vector<double> subset_scores(n_subsets, -MAX_NLNL);

    // [2026-07-22 added] progress counter so a long exhaustive run has some visible heartbeat instead of total silence until the whole parallel loop finishes. Prints roughly every 10% of subsets done (~10 lines total, regardless of n_subsets), guarded by a critical section since cout isn't thread-safe against interleaving from concurrent threads.
    long long completed = 0;
    long long print_every = std::max(1LL, n_subsets / 10);

    #pragma omp parallel for schedule(dynamic)
    for(long long mask = 0; mask < n_subsets; ++mask){
        vector<int> shift_eids;
        for(int i = 0; i < ncand; ++i){
            if(mask & (1LL << i)) shift_eids.push_back(all_candidates[i]);
        }

        evo_tree cand_tree = rtree;
        OPT_TYPE cand_opt  = opt_type;
        LNL_TYPE cand_lnl_type = lnl_type;
        cand_opt.rlc_shift_eids = shift_eids;

        vector<RateSet> ones(shift_eids.size(), RateSet(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
        update_edge_rates_rlc(cand_tree, shift_eids, ones);

        double cand_nlnl = MAX_NLNL;
        // Keep debug off inside the parallel region to avoid interleaved output.
        max_likelihood_BFGS(cand_tree, vobs, vobs_change, obs_decomp, comps,
                            cand_lnl_type, cand_opt, cand_nlnl, 0);

        RlcIC ic = compute_rlc_ic(cand_nlnl, (int)shift_eids.size(), n_types_per_edge, n_sites, criterion);
        subset_nlnls[mask]  = cand_nlnl;
        subset_aics[mask]   = ic.aic;
        subset_bics[mask]   = ic.bic;
        subset_scores[mask] = ic.score;

        long long done;
        #pragma omp atomic capture
        done = ++completed;
        if(done % print_every == 0 || done == n_subsets){
            #pragma omp critical(rlc_exhaustive_progress)
            {
                cout << "[RLC][exhaustive] progress: " << done << "/" << n_subsets << " subsets done" << endl;
            }
        }
    }

    // Serial scan: best-by-active-criterion (drives the write-back below), plus
    // best-by-AIC and best-by-BIC purely for reporting — both already known from the same
    // pass, unlike stepwise/GA which only explore one criterion-dependent path.
    const double eps = 1e-6;
    long long best_mask = 0, best_aic_mask = 0, best_bic_mask = 0;
    double best_score = subset_scores[0];
    double best_aic = subset_aics[0], best_bic = subset_bics[0];
    for(long long mask = 1; mask < n_subsets; ++mask){
        if(subset_scores[mask] > best_score + eps){ best_score = subset_scores[mask]; best_mask = mask; }
        if(subset_aics[mask] < best_aic - eps){ best_aic = subset_aics[mask]; best_aic_mask = mask; }
        if(subset_bics[mask] < best_bic - eps){ best_bic = subset_bics[mask]; best_bic_mask = mask; }
    }

    auto mask_to_shifts = [&](long long mask){
        vector<int> shifts;
        for(int i = 0; i < ncand; ++i) if(mask & (1LL << i)) shifts.push_back(all_candidates[i]);
        return shifts;
    };

    auto print_combo = [&](const char* label, long long mask){
        vector<int> shifts = mask_to_shifts(mask);
        cout << "[RLC][exhaustive] " << label << ": K=" << shifts.size() << " shift_eids=[ ";
        for(int s : shifts) cout << s << " ";
        cout << "] nlnl=" << subset_nlnls[mask]
             << " AIC=" << subset_aics[mask] << " BIC=" << subset_bics[mask] << endl;
    };

    // Always-on diagnostic (independent of the general debug flag), same convention as
    // stepwise_search_shift_edges's per-candidate prints — full dump of all subsets is
    // debug-gated since it can be hundreds of lines.
    if(debug){
        for(long long mask = 0; mask < n_subsets; ++mask) print_combo("candidate", mask);
    }
    print_combo("AIC-best", best_aic_mask);
    print_combo("BIC-best", best_bic_mask);

    // --- write back final results using the subset selected by opt_type.rlc_criterion ---
    vector<int> best_shift_eids = mask_to_shifts(best_mask);
    opt_type.rlc_shift_eids = best_shift_eids;

    vector<RateSet> ones(best_shift_eids.size(), RateSet(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
    update_edge_rates_rlc(rtree, best_shift_eids, ones);

    double final_nlnl = MAX_NLNL;
    max_likelihood_BFGS(rtree, vobs, vobs_change, obs_decomp, comps,
                        lnl_type, opt_type, final_nlnl, debug);

    int K = (int)best_shift_eids.size();
    RlcIC final_ic = compute_rlc_ic(final_nlnl, K, n_types_per_edge, n_sites, criterion);
    double raw_logL        = final_ic.raw_logL;
    double penalized_score = final_ic.score;

    opt_type.rlc_raw_logL = raw_logL;
    opt_type.rlc_penalized_score = penalized_score;
    rtree.rlc_shift_eids = opt_type.rlc_shift_eids;
    rtree.rlc_raw_logL = raw_logL;
    rtree.rlc_penalized_score = penalized_score;
    rtree.score = penalized_score;
    min_nlnl = -penalized_score;

    cout << "[RLC][exhaustive] final (criterion=" << (criterion == 0 ? "AIC" : "BIC")
         << "): K=" << K << " shift_eids=[ ";
    for(int s : opt_type.rlc_shift_eids) cout << s << " ";
    cout << "] raw_logL=" << raw_logL
         << " AIC=" << final_ic.aic << " BIC=" << final_ic.bic
         << " score=" << penalized_score
         << " (evaluated " << n_subsets << " subsets)" << endl;
}


void optimize_tree_by_bsr_mode(
    gsl_rng* r,
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug)
{
    // [2026-07-14 added] auto-calibrate the 5 global reference rates once (bsr_mode=0,
    // estmu=1) before running the real bsr_mode>0 search, so callers don't have to
    // manually run the program twice (once to calibrate, once to search). Jointly
    // estimating global rates and per-edge multipliers in one run is non-identifiable
    // (scaling rates up and multipliers down by the same factor leaves every edge rate,
    // and the likelihood, unchanged), so calibration must happen as a separate,
    // isolated sub-optimization first.
    // [2026-07-14 disabled] ran the calibration sub-optimization on every call, including
    // retries from maximize_tree_likelihood's outer `while(!(nlnl < MAX_NLNL))` loop
    // (cnetml.cpp) — if calibration itself had already succeeded and only the later
    // shift-edge search failed, a retry redid the whole calibration for nothing, with no
    // record of whether the redo was actually necessary.
    // if(opt_type.bsr_mode > 0){
    if(opt_type.bsr_mode > 0 && !opt_type.bsr_calibrated){
        OPT_TYPE calib_opt = opt_type;
        calib_opt.bsr_mode = 0;
        calib_opt.estmu = 1;
        calib_opt.rlc_shift_eids = {};

        evo_tree calib_tree = rtree;
        double calib_nlnl = MAX_NLNL;
        max_likelihood_BFGS(calib_tree, vobs, vobs_change, obs_decomp, comps,
                            lnl_type, calib_opt, calib_nlnl, debug);

        // [2026-07-14 added] only accept and remember the calibration if it actually
        // succeeded (calib_nlnl finite); on failure leave rtree's rates untouched and
        // bsr_calibrated false, so a retry will genuinely retry calibration too instead
        // of skipping it based on a failed attempt.
        if(calib_nlnl < MAX_NLNL){
        // only carry the 5 calibrated rate values over, not branch lengths — matches
        // what a manual second run would see (branch lengths come fresh from rtree)
        rtree.dup_rate      = calib_tree.dup_rate;
        rtree.del_rate      = calib_tree.del_rate;
        rtree.chr_gain_rate = calib_tree.chr_gain_rate;
        rtree.chr_loss_rate = calib_tree.chr_loss_rate;
        rtree.wgd_rate      = calib_tree.wgd_rate;

        RateSet global_rates(0, rtree.dup_rate, rtree.del_rate,
                             rtree.chr_gain_rate, rtree.chr_loss_rate, rtree.wgd_rate);
        rtree.init_edge_rates(global_rates);

        opt_type.bsr_calibrated = true;

        if(debug)
            cout << "[CALIBRATE] estimated reference rates: dup=" << rtree.dup_rate
                 << " del=" << rtree.del_rate
                 << " chr_gain=" << rtree.chr_gain_rate
                 << " chr_loss=" << rtree.chr_loss_rate
                 << " wgd=" << rtree.wgd_rate << endl;
        }
    }

    if(opt_type.bsr_mode == 3){
        // 0=exhaustive, 1=stepwise (default), 2=GA — exhaustive is the baseline method;
        // stepwise remains the practical default via cnetml.cpp's rlc_search_method
        // default_value(1).
        if(opt_type.rlc_search_method == 0){
            exhaustive_search_shift_edges(rtree, vobs, vobs_change, obs_decomp,
                                          comps, lnl_type, opt_type, min_nlnl, debug);
        } else if(opt_type.rlc_search_method == 1){
            stepwise_search_shift_edges(rtree, vobs, vobs_change, obs_decomp,
                                        comps, lnl_type, opt_type, min_nlnl, debug);
        } else if(opt_type.rlc_search_method == 2){
            genetic_search_shift_edges(r, rtree, vobs, vobs_change, obs_decomp,
                                       comps, lnl_type, opt_type, min_nlnl, debug);
        } else {
            cerr << "Error: unknown opt_type.rlc_search_method=" << opt_type.rlc_search_method << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        max_likelihood_BFGS(rtree, vobs, vobs_change, obs_decomp,
                            comps, lnl_type, opt_type, min_nlnl, debug);
    }
}
