#include <iostream>
#include <math.h>       /* pow */
#include <vector>
 
double mean(std::vector<double> const &a)
{
  double aux_out = 0;
  for(size_t i=0; i<a.size(); i++)
  {
   aux_out += a[i];
  }

  return aux_out / a.size();
}

double std_dev(std::vector<double> const &a)
{
  if (a.size() < 2)
    return -1;
  double aux_mean = mean(a);
  
  double aux_out = 0;
  for(size_t i=0; i<a.size(); i++)
  {
   aux_out += (a[i]*a[i]);
  }
  
  aux_out -= (a.size()*aux_mean*aux_mean);
  aux_out /= (a.size()-1);
  return sqrt(aux_out);
}
