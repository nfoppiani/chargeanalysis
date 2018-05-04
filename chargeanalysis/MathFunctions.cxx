#include <iostream>
#include <math.h> /* pow */
#include <vector>
#include <algorithm>

void remove_nan(std::vector<double> &a)
{
  a.erase(std::remove_if(a.begin(), 
                         a.end(), 
                         [](double x){return std::isnan(x);}), 
                         a.end());
}

double mean(std::vector<double> a)
{
  remove_nan(a);
  double aux_out = 0;
  for (size_t i = 0; i < a.size(); i++)
  {
    aux_out += a[i];
  }

  return aux_out / a.size();
}

double std_dev(std::vector<double> a)
{
  remove_nan(a);
  if (a.size() < 2)
    return -1;
  double aux_mean = mean(a);

  double aux_out = 0;
  for (size_t i = 0; i < a.size(); i++)
  {
    aux_out += (a[i] * a[i]);
  }

  aux_out -= (a.size() * aux_mean * aux_mean);
  aux_out /= (a.size() - 1);
  return sqrt(aux_out);
}

double median(std::vector<double> a)
{
  remove_nan(a);
  if (a.size() > 0)
  {
    std::nth_element(a.begin(), a.begin() + a.size() / 2, a.end());
    return a[a.size() / 2];
  }
  else
  {
    return -1;
  }
}