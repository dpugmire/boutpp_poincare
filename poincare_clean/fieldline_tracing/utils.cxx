#include "utils.h"
#include <fstream>
#include <iomanip>
#include <vector>


void writeArray1DToFile(std::vector<double>& array, const std::string& fname)
{
  auto fname2 = fname + ".c.txt";
  std::ofstream out(fname2, std::ofstream::out);
  auto nx = array.size();
  out << "(" << nx << ")" << std::endl;

  out << std::scientific << std::setprecision(10);
  int cnt = 0;
  for (size_t i = 0; i < nx; i++)
  {
    if (cnt > 5000)
      break;
    out << i << ", " << array[i] << std::endl;
    cnt++;
  }

  out.close();
}

void writeArray2DToFile(std::vector<std::vector<double>>& array, const std::string& fname)
{
  auto fname2 = fname + ".c.txt";
  std::ofstream out(fname2, std::ofstream::out);
  auto nx = array.size();
  auto ny = array[0].size();
  out << "(" << nx << ", " << ny << ")" << std::endl;

  out << std::scientific << std::setprecision(10);
  size_t x0 = 0, x1 = nx, y0 = 0, y1 = ny;
  int cnt = 0;
  int maxCnt = 5000;
  maxCnt = -1;

  //x1 /= 2;
  //y1 /= 2;

  for (size_t i = x0; i < x1; i++)
    for (size_t j = y0; j < y1; j++)
    {
      if (maxCnt > 0 && cnt > maxCnt)
        break;
      out << i << ", " << j << ", " << array[i][j] << std::endl;
      cnt++;
    }
  out.close();
}

void writeArray3DToFile(const std::vector<std::vector<std::vector<double>>>& array, const std::string& fname)
{
  auto val0 = array[123][79][101];
  auto val1 = array[192][47][200];

  auto fname2 = fname + ".c.txt";
  std::ofstream out(fname2, std::ofstream::out);
  auto nx = array.size();
  auto ny = array[0].size();
  auto nz = array[0][0].size();
  out << "(" << nx << ", " << ny << ", " << nz << ")" << std::endl;

  out << std::scientific << std::setprecision(10);
  size_t x0 = 0, x1 = nx, y0 = 0, y1 = ny, z0 = 0, z1 = nz;

  x1 /= 2;
  y1 /= 2;
  z1 /= 2;

  int cnt = 0;
  for (size_t i = x0; i < x1; i++)
    for (size_t j = y0; j < y1; j++)
      for (size_t k = z0; k < z1; k++)
      {
        if (cnt > 5000)
          break;
        out << i << ", " << j << ", " << k << ", " << array[i][j][k] << std::endl;
        cnt++;
      }
  out.close();
}
