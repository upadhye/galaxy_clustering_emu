#include<stdio.h>
#include<stdlib.h>
#include"cpgplot.h"

int main(int argc, char**argv)
{

  cpgopen("?");
  cpgenv(0.,1.,0.,1.,0,0);
  cpgend();
  return(0);
}
