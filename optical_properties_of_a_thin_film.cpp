#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#define _USE_MATH_DEFINES
using namespace std;

const int N = 13;
void indice_s (vector<double> Tsustrato, vector<double> & s);
void indice_n1 (vector<double> TM, vector<double> Tm, vector<double> s, vector<double> & n);
void grosor1 (vector<int> lambda, vector<double> n, vector<double> & d1);
void stand_deviation (vector<double> & data);
void number_m (vector<double> n, vector<int> lambda, vector<double> d1, vector<double> & m);
void grosor2 (vector<double> n1, vector<int> lambda, vector<double> m, vector<double> & d2);
void stand_deviation2 (vector<double> & data);
void indice_n2 (vector<double> d2, vector<int> lambda, vector<double> m, vector<double> & n2);
void Ti_Talf (vector<double> TM, vector<double> Tm, vector<double> & Ti, vector<double> & Talfa);
void funcion_x (vector<double> s, vector<double> n2, vector<double> TM, vector<double> Tm, vector<double> Ti, vector<double> Talfa, vector<double> & xTM, vector<double> & xTm, vector<double> & xTi, vector<double> & xTalfa);
void alpha (vector<double> d2, vector<double> xTM, vector<double> xTm, vector<double> xTi, vector<double> xTalfa, vector<double> & alfaTM, vector<double> & alfaTm, vector<double> & alfaTi, vector<double> & alfaTalfa);
void extincion_k (vector<int> lambda, vector<double> alfxTM, vector<double> alfxTm, vector<double> alfxTi, vector<double> alfxTalfa, vector<double> & kTM, vector<double> & kTm, vector<double> & kTi, vector<double> & kTalfa);
void print (vector<int> lambda, vector<double> Trasmitancia, vector<double> Ts, vector<double> TM, vector<double> Tm, vector<double> Ti, vector<double> Talfa, vector<double> n1, vector<double> n2, vector<double> d1, vector<double> d2, vector<double> m, vector<double> alfxTM, vector<double> alfxTm, vector<double> alfxTi, vector<double> alfxTalfa, vector<double> kTM, vector<double> kTm, vector<double> kTi, vector<double> kTalfa);

int main()
{
  std::cout.precision(5);
  std::cout.setf(std::ios::scientific);
 
  vector<int> lambda = {571, 588, 611, 633, 660, 690, 723, 761, 804, 853, 910, 975, 1055};
  
  vector<double> Transmitancia = {0.61877, 0.43551, 0.83093, 0.49233, 0.90196, 0.51235, 0.91135, 0.52226, 0.91424, 0.52989, 0.91483, 0.53637, 0.91531};

  vector<double> TM = {0.61877, 0.73522, 0.83093, 0.87813, 0.90196, 0.91135, 0.91353, 0.91400, 0.91424, 0.91451, 0.91483, 0.91509, 0.91531};

  vector<double> Tm = {0.38720, 0.43551, 0.47442, 0.49233, 0.50461, 0.51235, 0.51783, 0.52226, 0.52623, 0.52989, 0.53322, 0.53637, 0.54032};
  
  vector<double> Tsustrato = {0.91112, 0.91152, 0.91200, 0.91239, 0.91281, 0.91320, 0.91356, 0.91391, 0.91424, 0.91455, 0.91483, 0.91510, 0.91535};

  vector<double> s (N, 0.0);
  vector<double> n1 (N, 0.0);
  vector<double> d1 (N, 0.0);
  vector<double> m (N, 0.0);
  vector<double> d2 (N+2, 0.0);
  vector<double> n2 (N, 0.0);
  vector<double> Ti (N, 0.0);
  vector<double> Talfa (N, 0.0);
  vector<double> xTM (N, 0.0);
  vector<double> xTm (N, 0.0);
  vector<double> xTi (N, 0.0);
  vector<double> xTalfa (N, 0.0);
  vector<double> alfxTM (N, 0.0);
  vector<double> alfxTm (N, 0.0);
  vector<double> alfxTi (N, 0.0);
  vector<double> alfxTalfa (N, 0.0);
  vector<double> kTM (N, 0.0);
  vector<double> kTm (N, 0.0);
  vector<double> kTi (N, 0.0);
  vector<double> kTalfa (N, 0.0);
  
  indice_s (Tsustrato, s);
  indice_n1 (TM, Tm, s, n1);
  Ti_Talf (TM, Tm, Ti, Talfa);
  grosor1 (lambda, n1, d1);
  number_m (n1, lambda, d1, m);
  m = {12, 11.5, 11, 10.5, 10, 9.5, 9, 8.5, 8, 7.5, 7, 6.5, 6};
  grosor2 (n1, lambda, m, d2);
  indice_n2 (d2, lambda, m, n2);
  funcion_x (s, n2, TM, Tm, Ti, Talfa, xTM, xTm, xTi, xTalfa);
  alpha (d2, xTM, xTm, xTi, xTalfa, alfxTM, alfxTm, alfxTi, alfxTalfa);
  extincion_k (lambda, alfxTM, alfxTm, alfxTi, alfxTalfa, kTM, kTm, kTi, kTalfa);
  
  print (lambda, Transmitancia, Tsustrato, TM, Tm, Ti, Talfa, n1, n2, d1, d2, m, alfxTM, alfxTm, alfxTi, alfxTalfa, kTM, kTm, kTi, kTalfa);
  
  return 0;
}

void indice_s (vector<double> Tsustrato, vector<double> & s)
{
  double Ts = 0.0;
  
  for(int ii=0; ii<N; ii++)
    {
      Ts = Tsustrato[ii];
      s[ii] = 1.0/Ts + sqrt(1.0/(pow(Ts, 2)) - 1);
    }
}

void indice_n1 (vector<double> TM, vector<double> Tm, vector<double> s, vector<double> & n)
{
  for(int ii=0; ii<N; ii++)
    {
      double T_M = TM[ii], T_m = Tm[ii], s_ind = s[ii];

      if (ii < 8) //Región debil/media
	{
	  double N = 2*s_ind*((T_M-T_m)/(T_M*T_m))+(s_ind*s_ind+1)/2;
	  n[ii] = sqrt(N + sqrt(N*N-s_ind*s_ind));
	}

      else{ //Region transparente
	double M = (2*s_ind)/T_m - (s_ind*s_ind+1)/2;
	n[ii] = sqrt(M + sqrt(M*M-s_ind*s_ind));
      }
    }
}

void grosor1 (vector<int> lambda, vector<double> n, vector<double> & d1)
{
  for (int ii=0; ii<(N-2); ii++)
    {
      double l1 = lambda[ii], l2 = lambda[ii+2], n1 = n[ii], n2 = n[ii+2];
      d1[ii] = (l1*l2)/(2*(l2*n1-l1*n2));
    }
  stand_deviation (d1);
}

void stand_deviation (vector<double> & data)
{
  double sum = 0.0;
  double E = 0.0;
  int size = data.size()-2;
  
  for(int ii=0; ii<size; ii++){
    sum+=data[ii];
  }
  
  double mean = sum/size;
  data[N-2] = mean;
  
  for(int jj=0; jj<size; jj++){
    E += (data[jj] - mean)*(data[jj] - mean);
  }

  data[N-1] = sqrt((1.0/size)*E);
}

void number_m (vector<double> n, vector<int> lambda, vector<double> d1, vector<double> & m)
{
  for (int ii=0; ii<N; ii++)
    {
      m[ii] = (2*n[ii]*d1[N-2])/(lambda[ii]);
    }
}

void grosor2 (vector<double> n1, vector<int> lambda, vector<double> m, vector<double> & d2)
{
  for (int ii=0; ii<N; ii++)
    {
      d2[ii] = (m[ii]*lambda[ii])/(2*n1[ii]);
    }
  stand_deviation2 (d2);
}

void stand_deviation2 (vector<double> & data)
{
  double sum = 0.0;
  double E = 0.0;
  int size = data.size()-2;
  
  for(int ii=0; ii<size; ii++){
    sum+=data[ii];
  }
  
  double mean = sum/size;
  data[N] = mean;
  
  for(int jj=0; jj<size; jj++){
    E += (data[jj] - mean)*(data[jj] - mean);
  }

  data[N+1] = sqrt((1.0/size)*E);
}

void indice_n2 (vector<double> d2, vector<int> lambda, vector<double> m, vector<double> & n2)
{
  double d = d2[N];
  
  for (int ii=0; ii<N; ii++)
    {
      n2[ii] = (m[ii]*lambda[ii])/(2*d);
    }
}

void Ti_Talf (vector<double> TM, vector<double> Tm, vector<double> & Ti, vector<double> & Talfa)
{
  for (int ii=0; ii<N; ii++)
    {
      Ti[ii] = (2*TM[ii]*Tm[ii])/(TM[ii]+Tm[ii]);
    }

  for (int jj=0; jj<N; jj++)
    {
      Talfa[jj] = sqrt(TM[jj]*Tm[jj]);
    }
}

void funcion_x (vector<double> s, vector<double> n2, vector<double> TM, vector<double> Tm, vector<double> Ti, vector<double> Talfa, vector<double> & xTM, vector<double> & xTm, vector<double> & xTi, vector<double> & xTalfa)
{
  for (int ii=0; ii<N; ii++)
    {
      double T_M = TM[ii], n= n2[ii], s_ind = s[ii], EM = (8*n*n*s_ind)/(T_M)+(n*n-1)*(n*n-s_ind*s_ind);
      xTM[ii] = (EM-sqrt(EM*EM-pow(n*n-1, 3)*(n*n-pow(s_ind, 4))))/(pow(n-1, 3)*(n-s_ind*s_ind));
      //xTM[ii] = EM;
    }

  for (int ii=0; ii<N; ii++)
    {
      double T_m = Tm[ii], n= n2[ii], s_ind = s[ii], Em = (8*n*n*s_ind)/(T_m)-(n*n-1)*(n*n-s_ind*s_ind);
      xTm[ii] = (Em-sqrt(Em*Em-pow(n*n-1, 3)*(n*n-pow(s_ind, 4))))/(pow(n-1, 3)*(n-s_ind*s_ind));
      //xTm[ii] = Em;
    }

  for (int ii=0; ii<N; ii++)
    {
      double T_i = Ti[ii], n= n2[ii], s_ind = s[ii], F = (8*n*n*s_ind)/(T_i);
      xTi[ii] = (F-sqrt(F*F-pow(n*n-1, 3)*(n*n-pow(s_ind, 4))))/(pow(n-1, 3)*(n-s_ind*s_ind));
    }
  
  for (int ii=0; ii<N; ii++)
    {
      double T_alfa = Talfa[ii], n= n2[ii], s_ind = s[ii], G = (128*pow(n, 4)*s_ind*s_ind)/(T_alfa*T_alfa)+n*n*pow(n*n-1, 2)*pow(s_ind*s_ind-1, 2)+pow(n*n-1, 2)*pow(n*n-s_ind*s_ind, 2);
      xTalfa[ii] = sqrt(G - sqrt(G*G-pow(n*n-1, 6)*pow(n*n-pow(s_ind, 4), 2)))/(pow(n-1, 3)*(n-s_ind*s_ind));
    }
}

void alpha (vector<double> d2, vector<double> xTM, vector<double> xTm, vector<double> xTi, vector<double> xTalfa, vector<double> & alfaTM, vector<double> & alfaTm, vector<double> & alfaTi, vector<double> & alfaTalfa)
{
  double d = d2[N];

  for (int ii=0; ii<N; ii++)
    {
      alfaTM[ii] = -(log(xTM[ii])/d)*pow(10, 4); //nm^(-1)*10^(7)cm^(-1)/nm^(-1)= 10^4*10^3cm^(-1)------> unidades 10^3cm^(-1)
      alfaTm[ii] = -(log(xTm[ii])/d)*pow(10, 4);
      alfaTi[ii] = -(log(xTi[ii])/d)*pow(10, 4);
      alfaTalfa[ii] = -(log(xTalfa[ii])/d)*pow(10, 4);
    }

  
}

void extincion_k (vector<int> lambda, vector<double> alfxTM, vector<double> alfxTm, vector<double> alfxTi, vector<double> alfxTalfa, vector<double> & kTM, vector<double> & kTm, vector<double> & kTi, vector<double> & kTalfa)
{
  for (int ii=0; ii<N; ii++)
    {
      kTM[ii] = ((lambda[ii]*pow(10, -7))*alfxTM[ii])/(4*M_PI); //longitud de onda en cm
      kTm[ii] = ((lambda[ii]*pow(10, -7))*alfxTm[ii])/(4*M_PI);
      kTi[ii] = ((lambda[ii]*pow(10, -7))*alfxTi[ii])/(4*M_PI);
      kTalfa[ii] = ((lambda[ii]*pow(10, -7))*alfxTalfa[ii])/(4*M_PI);
    }
}

void print (vector<int> lambda, vector<double> Trasmitancia, vector<double> Ts, vector<double> TM, vector<double> Tm, vector<double> Ti, vector<double> Talfa, vector<double> n1, vector<double> n2, vector<double> d1, vector<double> d2, vector<double> m, vector<double> alfxTM, vector<double> alfxTm, vector<double> alfxTi, vector<double> alfxTalfa, vector<double> kTM, vector<double> kTm, vector<double> kTi, vector<double> kTalfa)
{
  cout<<"lambda"<<"\t"<<"\t"<<"Trans"<<"\t"<<"\t"<<"Ts"<<"\t"<<"\t"<<"TM"<<"\t"<<"\t"<<"Tm"<<"\t"<<"\t"<<"Ti"<<"\t"<<"\t"<<"Talfa"<<endl;
  
  for (int ii=0; ii<N; ii++)
    {
      cout<<lambda[ii]<<"\t"<<Trasmitancia[ii]<<"\t"<<Ts[ii]<<"\t"<<TM[ii]<<"\t"<<Tm[ii]<<"\t"<<Ti[ii]<<"\t"<<Talfa[ii]<<endl;
    }
  
  cout<<"lambda"<<"\t"<<"\t"<<"TM"<<"\t"<<"\t"<<"Tm"<<"\t"<<"\t"<<"n1"<<"\t"<<"\t"<<"d1"<<"\t"<<"\t"<<"m"<<"\t"<<"\t"<<"d2"<<"\t"<<"\t"<<"n2"<<endl;
  
  for (int ii=0; ii<N; ii++)
    {
      cout<<lambda[ii]<<"\t"<<TM[ii]<<"\t"<<Tm[ii]<<"\t"<<n1[ii]<<"\t"<<d1[ii]<<"\t"<<m[ii]<<"\t"<<d2[ii]<<"\t"<<n2[ii]<<endl;
    }
  
  cout<<"mean d2: "<<d2[N]<<"\t"<<"std dev d2: "<<d2[N+1]<<endl;

  cout<<"lambda"<<"\t"<<"alfa_TM"<<"\t"<<"\t"<<"alfa_Tm"<<"\t"<<"\t"<<"alfa_Ti"<<"\t"<<"\t"<<"alfa_Talfa"<<"\t"<<"kTM"<<"\t"<<"\t"<<"kTm"<<"\t"<<"\t"<<"kTi"<<"\t"<<"\t"<<"kTalfa"<<endl;

  for (int ii=0; ii<N; ii++)
    {
      cout<<lambda[ii]<<"\t"<<alfxTM[ii]<<"\t"<<alfxTm[ii]<<"\t"<<alfxTi[ii]<<"\t"<<alfxTalfa[ii]<<"\t"<<kTM[ii]<<"\t"<<kTm[ii]<<"\t"<<kTi[ii]<<"\t"<<kTalfa[ii]<<endl;
    }
}
