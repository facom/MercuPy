#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <SpiceUsr.h>

#define RAD2DEG 180/M_PI
#define DEG2RAD M_PI/180

int main(int argc,char *argv[])
{
  FILE *flog;
  int i;
  double to,mu;
  SpiceDouble elements[8];
  SpiceDouble state[6];

  //IF NUMBER OF PARAMTER IS UNAPPROPIATE RETURN ERROR
  if(argc<8){
    fprintf(stdout,"-1 -1 -1 -1 -1 -1 -1 -1\n");
    exit(0);
  }

  //GET ELEMENTS
  for(i=0;i<6;i++)
    state[i]=atof(argv[i+1]);
  //GET EPHEMERIS TIME
  to=atof(argv[i+1]);i++;
  mu=atof(argv[i+1]);

  //CONVERT ELEMENTS TO STATE
  if((state[0]*state[0]+state[1]*state[1]+state[2]*state[2])>0
     and
     (state[3]*state[3]+state[4]*state[4]+state[5]*state[5])>0){
    oscelt_c(state,to,mu,elements);
  }
  else
    elements[0]=elements[1]=elements[2]=elements[3]=
      elements[5]=elements[5]=elements[6]=elements[7]=-0.0;

  //*
  flog=fopen("state2elem.log","w");
  fprintf(flog,"State:x=%e y=%e z=%e vx=%e vy=%e vz=%e to=%e mu=%e\n",
	  state[0],state[1],state[2],
	  state[3],state[4],state[5],
	  to,mu);
  fprintf(flog,"Elements:q=%e e=%e i=%e n=%e g=%e M=%e to=%e mu=%e\n",
	  elements[0],elements[1],elements[2],elements[3],
	  elements[4],elements[5],elements[6],mu);
  fprintf(flog,"Ephemeris time:to=%e\n",to);
  fclose(flog);
  //*/

  //PRINT STATE
  for(i=0;i<8;i++){
    if(i>=2 && i<=5)
      elements[i]*=RAD2DEG;
    fprintf(stdout,"%e ",elements[i]);
  }
  fprintf(stdout,"\n");
}
