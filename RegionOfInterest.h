#ifndef RegionOfInterest_h
#define RegionOfInterest_h

double ROI_func(int region, double eget){
  double p[5];

  if(region == 1){
  // fit for median
  p[0] = 0.00130283;
  p[1] = -0.317612;
  p[2] = -0.942851;
  p[3] = -0.107378;
  p[4] =  1.27837;
  }

  if(region == 2){
  // fit for median
  p[0] = -0.00154114;
  p[1] = -0.401004;
  p[2] = -0.43574;
  p[3] = 0.243757;
  p[4] = 0.828763;
  }

  if(region == 3){
  // fit for median
  p[0] = 0.00643157;
  p[1] = -0.23502;
  p[2] = -0.467785;
  p[3] = 0.100745;
  p[4] = 0.855804;
  }

  if(region == 4){
  // fit for median
  p[0] = 0.00029069;
  p[1] = -0.163572;
  p[2] = -0.46938;
  p[3] = 0.197002;
  p[4] = 1.08246;
  }

  if(region == 5){
  // fit for median
  p[0] = 0.00612733;
  p[1] = -0.129135;
  p[2] = -0.13285;
  p[3] = 0.174302;
  p[4] = 0.0624058;
  }

  return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
}

#endif 
