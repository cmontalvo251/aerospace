#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <complex>

// %function uvw = Dryden(x,N,wupper,wlower,h)
// % function uvw = Dryden(x,N,wupper,wlower,h)
// % x is a 3xP vector or points to sample
// % wupper is a vector containing the upper limits of the waves
// % wlower is a vector containing the lower limits of the waves
// % Ni is the number of intervals such that the waves sampled are:
// % wi = linspace(wlower,wupper,Ni)
// % h is the height above ground
// %
// % The Dryden model uses the energy function
// %
// % E(W) = 8*o^2/(pi)*(1/(1+(L*W)^2)^3)
// %
// % This is placed in the Spectral Density Function
// %
// % Sij(w) = E(W)/(4*pi*W^4)*(W^2*Kdij-wi*wj)
// %
// % Note W = norm(w) where w = [w1,....wm] and Kdij is the Kronecker
// % delta function
// %
// % This is then factorized into H(w) where
// %
// % S(w) = H*conj(H)^T
// %
// % Finally the system is realized for a given set of random numbers
// %
// % Using Shinozuka we arrive at the equation for uvw
// %
// % uvw = [f1(x),f2(x),...,fm(x)];
// %
// % N=(N1,...,Nm)
// %
// % fj(x) = sum(m=1,j)sum(l=1,N)norm(Hjm(wl))sqrt(2*dw)*cos(wlp'*x+Ojm(wl)+PHIml)
// %
// % wlp = wl + delwl
// %
// % delwl = random number between -dwl/2,dwl/2
// %
// % wl = wlower + (l-1/2)*dwl
// %
// % dw = dw1*...*dwm
// %
// % Ojm(wl) = atan(Im(Hjm(wl))/Re(Hjm(wl)))
// %
// % PHIml = random number between 0,2pi

#define PI 3.14159265358979
//3/14 1:59:26

double randUniform()
{
  double out = 0.0;
  out = ((double)(rand() % 100)/(double)100);
  return out;
}

int main()
{

  //%%Example problem
  double xlimit = 500; //%%%%How big is our grid
  double res = 1; //%%%and at what resolution
  double wupper[3],wlower[3],dw[3],dwp[3];
  int N = 70;
  int ii,jj,kk;
  for (ii = 0;ii < 3; ii++)
    {
      wupper[ii] = 20; //again we may as well set this big
      wlower[ii] = -20;
      dw[ii] = (wupper[ii]-wlower[ii])/N;
      dwp[ii] = 0.1*dw[ii];
    }
  double h = (double)5000/(double)3.28;
  //double h = 298;
  double z = 200; //%%let z be a constant since it does not change that much

  //%%1D
  int P = xlimit/res;

  //%%Generate Length Scales and Variance
  double sigu = 1;
  double Lu = h;

  //%%Generate PHI
  double PHI1[N][N][N],PHI2[N][N][N],PHI3[N][N][N];
  double dwp_rand;

  for (ii = 0;ii<N;ii++)
    {
      for (jj = 0;jj<N;jj++)
	{
	  for (kk = 0;kk<N;kk++)
	    {
	      PHI1[ii][jj][kk] = 2*PI*randUniform();
	      PHI2[ii][jj][kk] = 2*PI*randUniform();
	      PHI3[ii][jj][kk] = 2*PI*randUniform();
	      //dwp_rand[ii][jj][kk] = randUniform();
	    }
	}
    }
  printf("Random Numbers Generated\n");

  //%%Generate f
  double alfa = 4*(sigu*sigu)*(pow(Lu,5))/(PI*PI);
  double factor = sqrt(2*dw[0]*dw[1]*dw[2]);
  int idx,jdx,k1,k2,k3;
  double x,y,wp1,wp2,wp3,W,e,f,u,v,w,S11,S12,S22,S31,S32,S33;
  double w1k1,w2k2,w3k3,delw1,delw2,delw3,W2;
  double H11,H21,H22,H31,H32,wpstate,f1,f2,f3,H33,L33;
  FILE *ufile,*vfile,*wfile;

  ufile = fopen("Uturb.txt","w");
  vfile = fopen("Vturb.txt","w");
  wfile = fopen("Wturb.txt","w");

  for (idx = 0;idx<P;idx++)
    {
      x = ((double)idx)*res; 
      for (jdx = 0;jdx<P;jdx++)
	{
	  y = ((double)jdx)*res;
	  printf("X = %f Y = %f \n",x,y);
	  u = 0;v=0;w=0;
	  for (k1 = 0;k1<N;k1++)
	    {
	      w1k1 = wlower[0] + ((double)k1+0.5)*dw[0];
	      //printf("%f \n",w1k1);
	      for (k2 = 0;k2<N;k2++)
		{
		  w2k2 = wlower[1] + ((double)k2+0.5)*dw[1];
		  for (k3 = 0;k3<N;k3++)
		    {
		      w3k3 = wlower[2] + ((double)k3+0.5)*dw[2];
		      dwp_rand = randUniform();
		      delw1 = -dwp[0]/2 + dwp[0]*dwp_rand;
		      delw2 = -dwp[1]/2 + dwp[1]*dwp_rand;
		      delw3 = -dwp[2]/2 + dwp[2]*dwp_rand;
		      wp1 = w1k1+delw1;
		      wp2 = w2k2+delw2;
		      wp3 = w3k3+delw3;
		      W = sqrt(wp1*wp1+wp2*wp2+wp3*wp3);
		      e = pow(1+(Lu*W*W*Lu),-3);
		      f = alfa*e;
		      //%S
		      W2 = W*W;
		      S11 = f*(W2-wp1*wp1);
		      S12 = f*-wp1*wp2;
		      S22 = f*(W2-wp2*wp2);
		      S31 = f*-wp1*wp3;
		      S32 = f*-wp2*wp3;
		      S33 = f*(W2-wp3*wp3);
		      //%H
		      H11 = sqrt(S11);
		      H21 = S12/H11;
		      H22 = sqrt(S22-H21*H21);
		      H31 = S31/H11;
		      H32 = (S32-H31*H21)/H22;
		      L33 = S33-H31*H31-H32*H32;
		      if (L33 < 0) L33 = -L33;
		      H33 = sqrt(L33);
		      wpstate = wp1*x + wp2*y + wp3*z;
		      f1 = cos(wpstate+PHI1[k1][k2][k3]);
		      f2 = cos(wpstate+PHI2[k1][k2][k3]);
		      f3 = cos(wpstate+PHI3[k1][k2][k3]);
		      u = u + fabs(H11)*f1;
		      v = v + fabs(H21)*f1 + fabs(H22)*f2;
		      w = w + fabs(H31)*f1 + fabs(H32)*f2 + fabs(H33)*f3;
		    }
		}
	    }
	  u = factor*u;
	  v = factor*v;
	  w = factor*w;
	  printf("u = %f v = %f w = %f \n",u,v,w);
	  fprintf(ufile,"%f ",u);
	  fprintf(vfile,"%f ",v);
	  fprintf(wfile,"%f ",w);
	} //end ycoordinate
      fprintf(ufile,"\n");
      fprintf(vfile,"\n");
      fprintf(wfile,"\n");
    } //end xcoordinate


}//end main
