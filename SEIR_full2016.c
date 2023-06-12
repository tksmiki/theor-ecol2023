//this is a time-scale-approximation version of SIER of KHV model
//from 20090227, revised at 090316

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "./util/nrutil.h"
#include "./util/forode.h"
#include "./util/nrutil.c"
#include "./util/rk4.c"

#define DELTAZERO 5.0e-4

#define EPSS 1.0e-6
#define ZERO 1.0e-20
#define PPAI 3.14159265

double larger(double a, double b);

void differential(double time, double in[], double out[]);

double dS_Cdt(double time, double vr[]);
double dE_Cdt(double time, double vr[]);
double dI_Cdt(double time, double vr[]);
double dR_Cdt(double time, double vr[]);
double dS_Pdt(double time, double vr[]);
double dE_Pdt(double time, double vr[]);
double dI_Pdt(double time, double vr[]);
double dR_Pdt(double time, double vr[]);

double alfa(double wt);    //temperature dependent resistance
double bet(double wt);  //temperature dependent re-infection
double gam(double wt);  //temperature dependent infection
double resident(double wt);  //resident_time within patch

int jS_C, jE_C, jI_C, jR_C;
int jS_P, jE_P, jI_P, jR_P;

double r, K, a, d0, dF, dV, f, c;
double alfa0, bet0, gam0;
double p;

double wt_C, wt_P;  //water temperature of patch 1 & 2.
double minwtKHV, maxwtKHV;
double WT_min, WT_max;  //maximum and minimum temperature of patch 1
double  WT_opt; //optimal temperature for KOI with tendancy to staying
double WT_diff; //difference in temperature between patches; wt2 = wt1 - WT_def

double *temp0, *ampli, *phase;
int *pwer;

double eff_R0; //effective R0

int main(int arg, char *argv[])
{	
  double t, tzero, end_time, deltat;
  double *v;
  double *dfdt;
  
  long n_variable;
  long i, j, k, m;
  long write_index;
  
  long ran_seed;
  unsigned long ran_seedmt;
  
  double ave_start_time, ave_end_time, cv_start_time, cv_end_time;
  
  double average_total, average_infective, average_infected, average_prevalence;
  double cumulative_mortality;

  double *average;
  double *cv;
  
  tzero = 0.0;
  end_time = 365*2000.0; //10000days
  //end_time = 365*100.0;
  deltat = DELTAZERO;
  
  n_variable = 8;
  v = dvector(1, n_variable);
  dfdt = dvector(1, n_variable);
  average = dvector(1, n_variable);
  cv = dvector(1, n_variable); 
  temp0 = dvector(1,2);
  ampli = dvector(1,2);
  phase = dvector(1,2);
  pwer = ivector(1,2);
  
  jS_C = 1;
  jE_C = 2;
  jI_C = 3;
  jR_C = 4;
  jS_P = 5;
  jE_P = 6;
  jI_P = 7;
  jR_P = 8;
  
  write_index = 1;
  
  r = log(2.0)/(4.4*365.0); //doubling time = 4.4 years
  //a = -log(0.68)/(15.0/1440.0); //15 minの暴露で32%の個体が感染 from Perelberg et al. 2005
  //a /= 100.0;
  //a = -log(0.05)/(40.0/1440.0);
  
  d0 = 1.0/(365.0*10.0);// average lifetime 10 years 
  dF = 0.0;//dF = 0.0
  dV = 0.1;
  a = dV/0.119;//from  (Perelberg et al. 2005) and estimate from ODE model
  //a *= 10.0;   //discount to 10%
  
  f = 0.0; //f = 100
  //c = 1.0/365.0;//1年間は免疫が続く
  c = 0.0/365.0;
  
  bet0 = 0.1; //仮定！！ (alpha_EtoS)
  alfa0 = (0.42/0.40)*bet0;  //by solving ODE Romen et al. 2003, alpha_EtoR
  gam0 = (1.0/7.0);  //from qualitatively, 7 days (Perelberg et al. 2005) alpha_EtoI
  
  //sensitivity analysis
  //gam0 *= 1.1;
  gam0 *= 0.95;
  
  WT_min = 2.2;
  WT_max = 34.5;
  WT_opt = 30.0;
  WT_diff = 2.0; 
  
  //minwtKHV = 18.0; //15 or 18 from Gilad 2003  []should be modified for Ecological Research]
  //maxwtKHV = 28.0;//28 or 30 from Gilad 2003
  
  minwtKHV = 15.0; //15 or 18 from Gilad 2003  []should be modified for Ecological Research]
  maxwtKHV = 30.0;//28 or 30 from Gilad 2003
  
  p = 0.5; //fraction of coastal region

  //for(k = 1; k <= 19; k++) {
    
    //p = 0.05*k;//fraction of Coastal area
  
  //we should start with E,I,R > 0 for rapid conversion to periodic solution
  v[jS_C] = 0.3;  //0.3
  v[jE_C] = 0.1;  //0.1
  v[jI_C] = 0.1;  //0.1
  v[jR_C] = 0.1;  //0.1
  v[jS_P] = 0.3;
  v[jE_P] = 0.1;
  v[jI_P] = 0.1;
  v[jR_P] = 0.1;
  
  for(j = 1;j <= 8; j++) average[j] = 0.0;
  for(t = tzero; t <= end_time +ZERO;) {
    
    //wt_C = 0.5*(WT_min + WT_max) + 0.5*(WT_max - WT_min)*sin(2.0*PPAI*t/365.0 - 0.5* PPAI); 
    //wt_P = wt_C - WT_diff;
    //wt_P = wt_C + 3.0*sin(2.0*PPAI*t/1.0);
    
    //from 'day_ave.csv' and '内湖温度ver2.xls'cd n	cd
    wt_C = 0.5*(31.6 + 2.5) + 0.5*(31.6 - 2.5)*sin(2.0*PPAI*(floor(t)/365.0 - 0.25)) + (3.5/2.0)*sin(2.0*PPAI*t/1.0);
    wt_P = 0.5*(25.4 + 4.0) + 0.5*(25.4 - 4.0)*sin(2.0*PPAI*(floor(t)/365.0 - 0.25)) + (1.2/2.0)*sin(2.0*PPAI*t/1.0);
    
    //wt_C = 25.0;
    //wt_P = 25.0;
    
    //wt_P = wt_C;
    //wt1 = temp0[1] + ampli[1]*power(sin(PPAI*t/1.0 - phase[1]), pwer[1]);
    //wt2 = temp0[2] + ampli[2]*power(sin(PPAI*t/1.0 - phase[2]), pwer[2]);
    
    
    //if (t > 365.0*10.0 && t < 365.0*10.0 + 0.1) {
    //v[jI_C] += 0.001*deltat;
      //v[jI_P] += 0.01*deltat;
      //}

    differential(t, v, dfdt);
    rk4(v, dfdt, n_variable, t, deltat, v, differential);
    t += deltat;
    
    //calculating average and cv//
    if(t >= 1949.0*365.0 && t <= 1999.0*365.00) {
       for(j = 1; j <= 8;j ++) average[j] += deltat*v[j]/(1999.0*365.0 - 1949.0*365.0);
    }
    
#ifdef TIME
    average_total = p*(v[jS_C] + v[jE_C] + v[jI_C] + v[jR_C]) + (1.0 - p)*(v[jS_P] + v[jE_P] + v[jI_P] + v[jR_P]);
    average_infective = p*v[jI_C] + (1.0 - p)*v[jI_P];
    average_infected = p*(v[jE_C] + v[jI_C]) + (1.0 - p)*(v[jE_P] + v[jI_P]);
    average_prevalence = average_infected/average_total;
    
    if(fmod(t, 365) > 0.001) cumulative_mortality += dV*average_infective*deltat;
    else cumulative_mortality = 0.0;
    
    if(t > 0.0*365.0 && write_index < 1.0/deltat)  write_index++;
    
    else if(t > 0.0*365.0){
      printf("%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t/365.0, v[jS_C], v[jE_C], v[jI_C], v[jR_C], v[jS_P], v[jE_P], v[jI_P], v[jR_P], wt_C, wt_P, average_total, average_infective, average_prevalence, cumulative_mortality);
      write_index = 1;
    }		
#endif  
  }//end of for t
#ifdef STAT
  printf("bet0_increase\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", p, average[1], average[2], average[3],average[4], average[5],average[6],average[7],average[8]);
#endif   
  //}//end of for k 
  
  return 0;
}  


void differential(double time, double in[], double out[])
{
  out[jS_C] = dS_Cdt(time, in);
  out[jE_C] = dE_Cdt(time, in);
  out[jI_C] = dI_Cdt(time, in);
  out[jR_C] = dR_Cdt(time, in);
  out[jS_P] = dS_Pdt(time, in);
  out[jE_P] = dE_Pdt(time, in);
  out[jI_P] = dI_Pdt(time, in);
  out[jR_P] = dR_Pdt(time, in);
}

double larger(double a, double b)
{
  if(a >= b) return a;
  else return b;
}

double alfa(double wt)
{
  double  temp1, temp2, alfa1, alfa2;
  
  //temp1 = 30.0;
  //temp2 = 22.0;
  //alfa1 = 0.023;
  //alfa2 = 0.013;
  if(wt >= maxwtKHV) return alfa0;
  else return 0.0;
  //return alfa1;
}		

double bet(double wt)
{
  if(wt >=  maxwtKHV) return bet0;
  else return 0.0;
  //return bet0;
}

double gam(double wt)
{
  if(wt >= minwtKHV && wt <= maxwtKHV) return gam0;
  else return 0.0;
  //return gam0;
}

double resident(double wt)
{
  double a16, b16, c16, a17, b17, c17, aH, bH, cH,d;
  
  a16 = 84.9430/1440.0; 
  b16 = 26.44;
  c16= 11.3107;
  
  a17 = 110.4103/1440.0; 
  b17 = 28.0842;
  c17= 2.5006;
  
  aH = a16;  //hypothetical example
  bH = b16 + 4.0;  //***_bp4 (plus 4.0 oC)
  cH = c16;
  
  d = 20.0/1440.0; //minimum value = 20 min
  
  return (a17*exp(-(wt-b17)*(wt-b17)/c17) + d);
  //return (a16*exp(-(wt-b16)*(wt-b16)/c16) + d);
  //return (aH*exp(-(wt-bH)*(wt-bH)/cH) + d);
  //return (a17 + d); //random behavior tends not to move
  //return (a16 + d);   //random behaviro tends not to move
  
 
}

double dS_Cdt(double time, double vr[])
{
  double temp;
  temp = r*(1.0 - vr[jS_C] - vr[jR_C])*(vr[jS_C] + vr[jR_C]) - a*vr[jS_C]*vr[jI_C] - (d0 + dF)*vr[jS_C]  + bet(wt_C)*vr[jE_C] + c*vr[jR_C] - (1.0 - p)*(1.0/resident(wt_C))*vr[jS_C] + (1.0 - p)*(1.0/resident(wt_P))*vr[jS_P];
  return temp;
}

double dE_Cdt(double time, double vr[])
{
  double temp = a*vr[jS_C]*vr[jI_C] - (alfa(wt_C) + bet(wt_C) + gam(wt_C))*vr[jE_C] - (d0 + dF)*vr[jE_C] - (1.0 - p)*(1.0/resident(wt_C))*vr[jE_C] +  (1.0 - p)*(1.0/resident(wt_P))*vr[jE_P];
  return temp;
}

double dI_Cdt(double time, double vr[])
{
  double temp;
  temp = gam(wt_C)*vr[jE_C] - (d0 + dF + dV)*vr[jI_C] - (1.0 - p)*(1.0/resident(wt_C))*vr[jI_C] + (1.0 - p)*(1.0/resident(wt_P))*vr[jI_P];
  return temp;	
}

double dR_Cdt(double time, double vr[])
{
  double temp;
  temp = alfa(wt_C)*vr[jE_C]  - (d0 + dF + c)*vr[jR_C] - (1.0 - p)*(1.0/resident(wt_C))*vr[jR_C] + (1.0 - p)*(1.0/resident(wt_P))*vr[jR_P];
  return temp;
}

double dS_Pdt(double time, double vr[])
{
  double temp;
  temp = r*(1.0 - vr[jS_P] - vr[jR_P])*(vr[jS_P] + vr[jR_P]) - a*vr[jS_P]*vr[jI_P] - (d0 + dF)*vr[jS_P] + c*vr[jR_P]  + bet(wt_P)*vr[jE_P] + p*(1.0/resident(wt_C))*vr[jS_C] -  p*(1.0/resident(wt_P))*vr[jS_P];
  //temp = 0.0;
  return temp;
}

double dE_Pdt(double time, double vr[])
{
  double temp = a*vr[jS_P]*vr[jI_P] - (alfa(wt_P) + bet(wt_P) + gam(wt_P))*vr[jE_P] - (d0 + dF)*vr[jE_P] + p*(1.0/resident(wt_C))*vr[jE_C] -  p*(1.0/resident(wt_P))*vr[jE_P];
  //temp = 0.0;
  return temp;
}

double dI_Pdt(double time, double vr[])
{
  double temp;
  temp = gam(wt_P)*vr[jE_P] - (d0 + dF + dV)*vr[jI_P] + p*(1.0/resident(wt_C))*vr[jI_C] - p*(1.0/resident(wt_P))*vr[jI_P];
  //temp = 0.0;
  return temp;	
}

double dR_Pdt(double time, double vr[])
{
  double temp;
  temp = alfa(wt_P)*vr[jE_P]  - (d0 + dF + c)*vr[jR_P] + p*(1.0/resident(wt_C))*vr[jR_C] - p*(1.0/resident(wt_P))*vr[jR_P];
  //temp = 0.0;
  return temp;
}


