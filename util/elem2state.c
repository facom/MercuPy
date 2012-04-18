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
  double to;
  SpiceDouble elements[8];
  SpiceDouble state[6];

  //IF NUMBER OF PARAMTER IS UNAPPROPIATE RETURN ERROR
  if(argc<10){
    fprintf(stdout,"-1 -1 -1 -1 -1 -1\n");
    exit(0);
  }
  for(i=0;i<8;i++){
    elements[i]=atof(argv[i+1]);
    //ANGLES CONVERTED TO RADIANS
    if(i>=2 && i<=5)
      elements[i]*=DEG2RAD;
  }
  to=atof(argv[i+1]);

  //*
  flog=fopen("elem2state.log","w");
  fprintf(flog,"Elements:q=%e e=%e i=%e n=%e g=%e M=%e to=%e mu=%e\n",
	  elements[0],elements[1],elements[2],elements[3],
	  elements[4],elements[5],elements[6],elements[7]);
  fprintf(flog,"Ephemeris time:to=%e\n",to);
  fclose(flog);
  //*/

  //CONVERT ELEMENTS TO STATE
  if(elements[0]>0){
    conics_c(elements,to,state);
  }else{
    state[0]=state[1]=state[2]=state[3]=state[4]=state[5]=-0.0;
  }

  //PRINT STATE
  for(i=0;i<6;i++)
    fprintf(stdout,"%e ",state[i]);
  fprintf(stdout,"\n");
}
