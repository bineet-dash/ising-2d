#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
// #include "../lib.h"

using namespace std;
ofstream outfiledata("data.txt");
ofstream outfiledebug("debug.txt");

float J; int progress_percent;

int periodic(int base, int addendum, int limit) //this is addition for Z_n group starting from 0.
{
  int result=base+addendum;
  if(result>limit){result= result%(limit+1);}
  else if(result < 0){
    while(result < 0)
    {
      result+=limit+1;
    }
    result= result%(limit+1);
  }
  return result;
}

float calc_delE(int** arr, int x, int y, int L) //L=lattice_size
{
  float del_E=2*J*arr[x][y]; //below take L-1 as the L=lattice_size starting from 1
  float sum_nn= arr[periodic(x,-1,L-1)][y]+ arr[periodic(x,1,L-1)][y]+arr[x][periodic(y,-1,L-1)]+arr[x][periodic(y,1,L-1)];
  del_E*=sum_nn;
  return del_E;
}

float energy(int** arr, int L)
{
  int energy= 0;
  for(int x=0; x<L; x++)
  {
    for(int y=0;y<L; y++)
    {
      energy += arr[x][y]*arr[periodic(x,-1,L-1)][y]+arr[x][y]*arr[periodic(x,1,L-1)][y]+ arr[x][y]*arr[x][periodic(y,-1,L-1)]+arr[x][y]*arr[x][periodic(y,1,L-1)];
    }
  }
  float result= -1*float(energy)*J/(L*L);
  return result;
}

float magnetization(int** arr, int L)
{
  int mag=0;
  for(int x=0; x<L; x++)
  {
    for(int y=0; y<L; y++)
      mag+= arr[x][y];
  }
  return abs(mag)/(L*L);
}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double ran0(long *idum)
{
   long  k;
   double ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

int main()
{
  int lattice_size,suggested_x, suggested_y, number_cycles, spin;
  float original_energy, suggested_energy,jump_probability, beta, del_E, total_energy,debug, initial_temperature, final_temperature;
  int** spinarray;
  cout << "Enter the size of the lattice and J: ";
  cin >> lattice_size >> J;
  cout << "Enter the no of cycles: ";
  cin >> number_cycles;
  cout << "Enter the initial and final temperature: ";
  cin >> initial_temperature >> final_temperature;
  long idum = -1;

  spinarray = new int*[lattice_size];
  for(int i=0;i<lattice_size;i++)
  {
    spinarray[i]= new int[lattice_size];
  }

  for(int i=0; i<lattice_size;i++)
  {
    for(int j=0;j<lattice_size;j++)
    {
      if(ran0(&idum)<=0.5) {spinarray[i][j]=1;}
      else {spinarray[i][j]=-1;}
    }
  }

  float temperature_step = 0.05;

  for(float temperature=final_temperature; temperature >= initial_temperature; temperature-=temperature_step)
  {
    beta= 1/temperature;

    for(int cycle=1; cycle<=number_cycles/2; cycle++)
    {
      for(int x=0; x<lattice_size; x++)
      {
        for(int y=0; y<lattice_size; y++)
        {
          suggested_x=floor(lattice_size*ran0(&idum));
          suggested_y=floor(lattice_size*ran0(&idum));
          del_E= calc_delE(spinarray, suggested_x, suggested_y, lattice_size);
          jump_probability= exp(-beta*del_E);
          if(ran0(&idum)<=jump_probability)
          {
            spinarray[suggested_x][suggested_y]*=-1;
          }
        }
       }
    }

    double total_M  = 0.0;
    double total_M2 = 0.0;
    double total_E  = 0.0;
    double total_E2 = 0.0;

    for(int cycle=number_cycles/2; cycle<number_cycles; cycle++)
    {
      for(int x=0; x<lattice_size; x++)
      {
        for(int y=0; y<lattice_size; y++)
        {
          suggested_x=floor(lattice_size*ran0(&idum));
          suggested_y=floor(lattice_size*ran0(&idum));
          del_E= calc_delE(spinarray, suggested_x, suggested_y, lattice_size);
          jump_probability= exp(-beta*del_E);
          if(ran0(&idum)<=jump_probability)
          {
            spinarray[suggested_x][suggested_y]*=-1;
          }
        }
       }
      total_M  += magnetization(spinarray, lattice_size);
      total_M2 += pow(magnetization(spinarray, lattice_size),2);
      total_E  += energy(spinarray, lattice_size);
      total_E2 += pow(energy(spinarray, lattice_size),2);
    }

    int no_msmt = number_cycles/2;

    double avg_M  = total_M/no_msmt;
    double avg_M2 = total_M2/no_msmt;
    double avg_E  = total_E/no_msmt;
    double avg_E2 = total_E2/no_msmt;

    double err_E = sqrt(avg_E2 - pow(avg_E,2))/sqrt(no_msmt);
    double err_M = sqrt(avg_M2 - pow(avg_M,2))/sqrt(no_msmt);

    double cv = (avg_E2 - pow(avg_E,2)) /pow(temperature,2); 
    double xi = (avg_M2 - pow(avg_M,2)) /temperature; 

    outfiledata << setw(5) << temperature << " " << setw(10) <<  avg_M <<  " " << setw(10) <<  err_M << " " << setw(10)
                 <<  avg_E << " " << setw(10) <<  err_E << " " << setw(10) <<  cv << " " << setw(10) <<  xi << endl;

    progress_percent=int((final_temperature-temperature)*100/(final_temperature-initial_temperature));
    cout.flush();
    cout << "\r [ "<< progress_percent+1 << "% ] ";
    for(int i=0; i<progress_percent; i++)
    {
      cout << "#";
    }
  }
  cout << endl;
}