#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <sys/time.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "dispersion.h"
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define GW 2.0
#define PI 3.14149
#define REL 1.0e-3
#define ABS 1e-7
#define DEBUG_MODE false
#define T_TRESHOLD 1e-12



//Bose-distributions
double nB(double p, double T){
    double ans = 0.0;
    ans = 1.0/(exp(p/T)-1.0);
    return ans;
}



//Fermi-Dirac distribution
    double nF(double p, double T, double M){
        double ans = 0.0;
        ans =1.0/(exp(sqrt(p*p+M*M)/T)+1.0);
        return ans;
    }

//function which return -1 for positive nunber
int whatsthatsign(double x){
	if (x > 0) return 1;
	if (x < 0) return -1;
return 0;
}




double LpB(double p, double M, double k0, double k1){
    double KK = k0*k0-k1*k1;
    
    double top1 = KK-M*M+2*p*(k0+k1);
    double bot1 = KK-M*M+2*p*(k0-k1);
    
    double top2 = KK-M*M-2*p*(k0-k1);
    double bot2 = KK-M*M-2*p*(k0+k1);
    double arg1 = top1/bot1;
    double arg2 = top2/bot2;
    double arg3 = (top1*top2)/(bot1*bot2);


    
//    if(arg1<0){std::cout<<"ERROR LpB: arg1 <0"<<std::endl;};
//    if(arg1==0){std::cout<<"ERROR LpB: arg1 =0"<<std::endl;};
//    if(arg2<0){std::cout<<"ERROR LpB: arg2 <0"<<std::endl;};
//    if(arg2==0){std::cout<<"ERROR LpB: arg2 =0"<<std::endl;};
    

    double arg4= gsl_sf_log_abs(top1)-gsl_sf_log_abs(bot1)+gsl_sf_log_abs(top2)-gsl_sf_log_abs(bot2);
    
//   if(arg3==0|arg3<0){std::cout<<" LMB "<<" k1 "<<k1<<" k0 "<<k0<<" k0-k1 "<<k0-k1<<" p "<<p<<" top1 "<<top1<<" top2 "<<top2<<" bot1 "<<bot1<<" bot2 "<<bot2<<" arg3 "<<arg3<<" arg4 "<<arg4<<std::endl;};


    return arg4;

   

}





double LmB(double p, double M, double k0, double k1){
    
    double KK = k0*k0-k1*k1;
    
    double top1 = KK-M*M+2*p*(k0+k1);
    double bot1 = KK-M*M+2*p*(k0-k1);
    
    double top2 = KK-M*M-2*p*(k0-k1);
    double bot2 = KK-M*M-2*p*(k0+k1);
    double arg1 = top1/bot1;
    double arg2 = top2/bot2;
    double arg3 = (top1*bot2)/(bot1*top2);
    //if(arg3==0||arg3<0){std::cout<<"ERROR LmB: arg3 =0"<<" k1 "<<k1<<" k0 "<<k0<<" k0-k1 "<<k0-k1<<" p "<<p<<" M "<<M<<std::endl;};

    
    
//    if(arg1<0){std::cout<<"ERROR LmB: arg1 <0"<<std::endl;};
//    if(arg1==0){std::cout<<"ERROR LmB: arg1 =0"<<std::endl;};
//    if(arg2<0){std::cout<<"ERROR LmB: arg2 <0"<<std::endl;};
//    if(arg2==0){std::cout<<"ERROR LmB: arg2 =0"<<std::endl;};
//
    return gsl_sf_log_abs(top1)-gsl_sf_log_abs(bot1)-gsl_sf_log_abs(top2)+gsl_sf_log_abs(bot2);

    //return gsl_sf_log_abs(top1*bot2/(bot1*top2));
   
//    if(top1==0||bot1==0){std::cout<<"ERROR LmB: top1: "<<top1<<" bot1: "<<bot1<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
//    if(top2==0||bot2==0){std::cout<<"ERROR LmB: top2: "<<top2<<" bot2: "<<bot2<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
//    

}


double LpF(double p, double M, double k0, double k1){
    
    double KK = k0*k0-k1*k1;
    double E = sqrt(p*p+M*M);
    
    double top1 = KK+M*M+2*E*k0+2*p*k1;
    double bot1 = KK+M*M+2*E*k0-2*p*k1;
    
    double top2 = KK+M*M-2*E*k0+2*p*k1;
    double bot2 = KK+M*M-2*E*k0-2*p*k1;
    double arg1 = top1/bot1;
    double arg2 = top2/bot2;
    double arg3 = (top1*top2)/(bot1*bot2);
    return gsl_sf_log_abs(top1)-gsl_sf_log_abs(bot1)+gsl_sf_log_abs(top2)-gsl_sf_log_abs(bot2);

    
}


double LmF(double p, double M, double k0, double k1){
    
    double KK = k0*k0-k1*k1;
    double E = sqrt(p*p+M*M);
    
    double top1 = KK+M*M+2*E*k0+2*p*k1;
    double bot1 = KK+M*M+2*E*k0-2*p*k1;
    
    double top2 = KK+M*M-2*E*k0+2*p*k1;
    double bot2 = KK+M*M-2*E*k0-2*p*k1;
    double arg1 = top1/bot1;
    double arg2 = top2/bot2;
    double arg3 = (top1*bot2)/(bot1*top2);

    
    return gsl_sf_log_abs(top1)-gsl_sf_log_abs(bot1)-gsl_sf_log_abs(top2)+gsl_sf_log_abs(bot2);
    
 
}



//integrand Tr(k.Sigma)
double TRKS_integrand(double p, void *d){
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;
    double f1=(4+((M*M-k1*k1+k0*k0)/(2*p*k1))*LpB(p,M,k0,k1))*nB(p,T)/p;
    
    double f2=(4-((M*M-k1*k1+k0*k0)/(2*p*k1))*LpF(p,M,k0,k1))*nF(p,T,M)/sqrt(p*p+M*M);
    
    
    return (1/(8*PI*PI))*p*p*(f1+f2);
    
}




// integrating Tr(k.Sigma)
double TRKS(double T, double M,double k0, double k1)
{
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (50000);
    struct dispParams params = {T,M,k0,k1};
    double result,error;
    gsl_function F;
    /*
     *Base is the generic starting value for what I base the integral size off.
     * Then this increases temperatyre in steps of "T" until I have reached a value below that the predefined threshold, T_TRESH
     * It then integrates this from 0 to value of T where cross the T_Treshold.
     */
    double base = TRKS_integrand(T,&params);
    double iT = 2;
    
    for(iT=2;fabs(TRKS_integrand(iT*T,&params)/base)>=T_TRESHOLD;iT=iT+1)
    {
        if(DEBUG_MODE){std::cout<<std::setprecision(9)<<" # TRKS integral base "<<base<<" iT: "<<iT<<" rat:  "<<fabs(TRKS_integrand(iT*T,&params)/base)<<std::endl;}
    }
    
    
    F.function = &TRKS_integrand;
    F.params = &params;
    
    gsl_integration_qags(&F,0.0, iT*T, ABS, REL, 50000, w, &result, &error);
    gsl_integration_workspace_free (w);
    
    
    return result;
    
    
}


//integrand Tr(u.Sigma)
double TRUS_integrand(double p, void *d){
    struct dispParams * params = (struct dispParams *)d;
    double M = params->M;
    double k1 = params->k1;
    double k0 = params->k0;
    double T = params->T;
    double f1 = (LmB(p,M,k0,k1)+(k0/p)*LpB(p,M,k0,k1))*nB(p,T);
    double f2 = nF(p,T,M)*LmF(p,M,k0,k1);
    
    return 1/(8*PI*PI)*(p/k1)*(f1+f2);
}


// integrating Tr(u.Sigma)
double TRUS(double T, double M,double k0, double k1)
{
    
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (50000);
    struct dispParams params = {T,M,k0,k1};
    double result,error;
    gsl_function F;
    F.function = &TRUS_integrand;
    F.params = &params;
    
    double base = TRUS_integrand(T,&params);
    double iT = 2;
    
    for(iT=2;fabs(TRUS_integrand(iT*T,&params)/base)>=T_TRESHOLD;iT=iT+1)
    {
        if(DEBUG_MODE){	std::cout<<std::setprecision(9)<<"# TRKU integral base "<<base<<" iT: "<<iT<<" rat:  "<<fabs(TRUS_integrand(iT*T,&params)/base)<<std::endl;}
    }
    //std::cout << "now: " << iT << std::endl;
    gsl_integration_qags(&F,0.0, iT*T, ABS, REL, 50000, w, &result, &error);
    
    
    gsl_integration_workspace_free (w);
    
    return result;
    
}
double TRS_integrand(double p, void*d){
	
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;


	return -1.0/(4.0*PI*PI)*p/k1*(1.0/k0*nF(p,T,M)*LpF(p,M,k0,k1)-1.0/p*nB(p,T)*LpB(p,M,k0,k1));


}

double TRS(double T,double M,double k0,double k1)
{

    
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (50000);
	struct dispParams params = {T,M,k0,k1};
	double result,error;
	gsl_function F;
	F.function = &TRS_integrand;
	F.params = &params;
	
	double base = TRS_integrand(T,&params);
	double iT = 2;

	for(iT=2;fabs(TRS_integrand(iT*T,&params)/base)>=T_TRESHOLD;iT=iT+1)
	  {
	 // if(DEBUG_MODE){	std::cout<<std::setprecision(9)<<"# TRKU integral base "<<base<<" iT: "<<iT<<" rat:  "<<fabs(TRUS_integrand(iT*T,&params)/base)<<std::endl;}
		}
    //std::cout << "now: " << iT << std::endl;
	gsl_integration_qags(&F,0.0, iT*T, ABS, REL, 50000, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;

}



double sigmaA(double T, double M, double k0, double k1){
    
    
    	return (1/(k1*k1))*(TRKS(T,M,k0,k1)-k0*TRUS(T,M,k0,k1));
}


double sigmaB(double T, double M, double k0, double k1){
	return (pow(k0/k1,2)-1)*TRUS(T,M,k0,k1)-(k0/pow(k1,2))*TRKS(T,M,k0,k1);
   }

double sigmaC(double T, double M, double k0, double k1){
	return TRS(T,M,k0,k1);
}

//opposite now, put in k0 get k1
double dispEquation(double k0, void * d){
	struct dispParams2 * params = (struct dispParams2 *)d;
	double T = params->T;
	double M = params->M;
    double k1 = params->k1;
    
    double ans_mass= (1+sigmaA(T,M,k0,k1))*k0+sigmaB(T,M,k0,k1)-k1*(1+sigmaA(T,M,k0,k1));
    return ans_mass;
  
}

double dispSolved(double T, double M, double k1){

	int iter = 0, max_iter = 250;
   	int status;
	double r = 0;
	// INitialise gsl root solving algorithm, going to use brent as its a brackeing algorithm that as long as we get a point on either sie, is guarrenteed to
	// find the contained root.
	
	const gsl_root_fsolver_type * S   = gsl_root_fsolver_brent;
	gsl_root_fsolver * s     = gsl_root_fsolver_alloc (S);
    struct dispParams2 p = {T,M,k1};

	//This for loop will start a bit to the left and right of input k0 and increas:e the bracket until we get opposite signs
	//Will only fail if its all negative or positive, which can happen at very low k/T (presumably due to not complete expressions)
	
	double x_lo=k1;
	double x_hi=k1;
	

	double m = 1e-8;
    for(m = 1e-5; whatsthatsign(dispEquation(k1*pow(10,m),&p))==whatsthatsign(dispEquation(k1*pow(10,-m),&p)); m=m+0.005)
   	{
        if(DEBUG_MODE){std::cout<<m<<" "<<pow(10,m)<<" "<<whatsthatsign(dispEquation(k1*pow(10,m),&p))<<" "<<whatsthatsign(dispEquation(k1*pow(10,-m),&p))<<std::endl;}
        if(DEBUG_MODE){std::cout<<m<<" "<<pow(10,m)<<" "<<dispEquation(k1*pow(10,m),&p)<<" "<<-m<<" "<<pow(10,-m)<<" "<<dispEquation(k1*pow(10,-m),&p)<<std::endl;}
		if(m>4)
		{
			//std::cout<<"# ERROR: dispSolved@'dispersion.c' rootfinder: Wandered from k0, didnt bracket the root. "<<std::endl;
		}
	
	}
		//std::cout<<m*multiplier<<" "<<dispEquation(w0*(m*multiplier),&p)<<" "<<dispEquation(w0/(m*multiplier),&p)<<std::endl;
		x_lo=std::min(k1*pow(10,-m),k1*pow(10,m));
		x_hi=std::max(k1*pow(10,-m),k1*pow(10,m));


	//Same as with integration, we need to set it up in terms of a gsl_function F .
	gsl_function F;
	F.function =&dispEquation;
	F.params =&p;

	gsl_root_fsolver_set(s,&F,x_lo,x_hi);
	

	// In general dont really need much iterations, very efficient algorithm
	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_x_lower(s);
		x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo,x_hi,0,1e-8);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
    
   // if(iter>max_iter-2){std::cout<<" Hitting max iter?"<<std::endl;}

	// Always free your fsolver for memory reasons
	gsl_root_fsolver_free (s); 


    if(r<0){std::cout<<"#ERROR: k0 is negative energy state"<<std::endl;}
    
	return r;

}

