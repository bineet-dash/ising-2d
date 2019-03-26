#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

double J; int progress_percent;


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

double calc_delE(int** arr, int x, int y, int L) //L=L
{
  double del_E=2*J*arr[x][y]; //below take L-1 as the L=L starting from 1
  double sum_nn= arr[periodic(x,-1,L-1)][y]+ arr[periodic(x,1,L-1)][y]+arr[x][periodic(y,-1,L-1)]+arr[x][periodic(y,1,L-1)];
  del_E*=sum_nn;
  return del_E;
}

double calc_energy(int** arr, int L)
{
  double energy= 0;
  for(int x=0; x<L; x++)
  {
    for(int y=0;y<L; y++)
    {
      energy += arr[x][y]*arr[periodic(x,-1,L-1)][y]+arr[x][y]*arr[periodic(x,1,L-1)][y]+ arr[x][y]*arr[x][periodic(y,-1,L-1)]+arr[x][y]*arr[x][periodic(y,1,L-1)];
    }
  }
  double result= -1*double(energy)*J;
  return result/double(L*L);
}

double calc_magnetization(int** arr, int L)
{
  double mag=0;
  for(int x=0; x<L; x++)
  {
    for(int y=0; y<L; y++)
      mag+= arr[x][y];
  }
  return abs(mag)/double(L*L);
}

int main(int argc, char* argv[])
{
  if(argc!=3)
  {
    std::cerr << "Enter (1) lattice size, (2) number of MC sweeps. " << endl;
    exit(1);
  }

   J = 1.0;

  int L = std::atoi(argv[1]);
  int no_sweeps = std::atoi(argv[2]);
  int no_therm = no_sweeps/2;
  double no_msmt = no_sweeps-no_therm;
  long idum = time(NULL);
  double initial_temperature=0.01, final_temperature=5.0;

  // double initial_temperature, final_temperature;
  // cout << "Enter the size of the lattice and J: ";
  // cin >> L >> J;
  // cout << "Enter the no of cycles: ";
  // cin >> number_cycles;
  // int no_therm = number_cycles/2;
  // cout << "Enter the initial and final temperature: ";
  // cin >> initial_temperature >> final_temperature;
  // long idum = -1;
  

  ofstream outfiledata("results.dat");
  int** spinarray;
  spinarray = new int*[L];
  for(int i=0;i<L;i++)
  {
    spinarray[i]= new int[L];
  }

  for(int i=0; i<L;i++)
  {
    for(int j=0;j<L;j++)
    {
      if(ran0(&idum)<=0.5) {spinarray[i][j]=1;}
      else {spinarray[i][j]=-1;}
    }
  }

  double temperature_step = 0.05;

  for(double temperature=final_temperature; temperature >= initial_temperature; temperature-=temperature_step)
  {
    double total_energy = 0;
    double total_magnetization = 0;
    double total_energy2 = 0;
    double total_magnetization2 = 0;

    for(int cycle=1; cycle<=no_therm; cycle++)
    {
      for(int x=0; x<L; x++)
      {
        for(int y=0; y<L; y++)
        {
          int suggested_x=floor(L*ran0(&idum));
          int suggested_y=floor(L*ran0(&idum));
          double del_E= calc_delE(spinarray, suggested_x, suggested_y, L);
          if(ran0(&idum)<=exp(-del_E/temperature))
          {
            spinarray[suggested_x][suggested_y]*=-1;
          }
        }
      }
    }

    for(int cycle=no_therm; cycle<no_sweeps; cycle++)
    {
      for(int x=0; x<L; x++)
      {
        for(int y=0; y<L; y++)
        {
          int suggested_x=floor(L*ran0(&idum));
          int suggested_y=floor(L*ran0(&idum));
          double del_E= calc_delE(spinarray, suggested_x, suggested_y, L);
          if(ran0(&idum)<=exp(-del_E/temperature))
          {
            spinarray[suggested_x][suggested_y]*=-1;
          }
        }
       }
      double energy = calc_energy(spinarray, L);
      double magnetization = calc_magnetization(spinarray, L);

      total_energy += energy;
      total_magnetization += magnetization;
      total_energy2 += pow(energy,2);
      total_magnetization2 += pow(magnetization,2); 
    }

    double avg_E =  total_energy/no_msmt;
    double avg_E2 = total_energy2/no_msmt;
    double avg_M =  total_magnetization/no_msmt;
    double avg_M2 = total_magnetization2/no_msmt;

    double err_E = sqrt(avg_E2 - pow(avg_E,2));
    double err_M = sqrt(avg_M2 - pow(avg_M,2));

    double cv = (avg_E2 - pow(avg_E,2))/pow(temperature,2);
    double xi = (avg_M2 - pow(avg_M,2))/temperature;


    outfiledata << setw(5) << temperature << " " << setw(10) <<  avg_M <<  " " << setw(10) <<  err_M << " " << setw(10)
                 <<  avg_E << " " << setw(10) <<  err_M << " " << setw(10) <<  cv << " " << setw(10) <<  xi << endl;

    cout << temperature << " done. \r"; cout.flush();
  }
  cout << endl;
}



/* progress_percent=int((final_temperature-temperature)*100/(final_temperature-initial_temperature));
cout.flush();
cout << "\r [ "<< progress_percent+1 << "% ] ";
for(int i=0; i<progress_percent; i++)
{
  cout << "#";
} */