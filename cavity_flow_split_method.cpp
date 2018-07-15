/* Author: Juan Carlos Cuevas Bautista
*  Date:Dec 20,2014.
*  
*  Program to compute the numerical solution of the Navier Stokes Equations 
*  for the 2D cavity flow in a square grid using the splitting method
*  x:   u_t+uu_x+vu_y=-p_x+1/Re_l*(u_{xx}+u_{yy}),  0<x<2
*  y:   v_t+uv_x+vv_y=-p_y+1/Re_l*(v_{xx}+v_{yy}),  0<y<2
*  Boundary Conditions square function for the top free-stream velocity
*  u=0 0<x<dx and 2-dx<x<dx @y=2  
*  u=1 dx<=x<=2-dx @y=2  
*  u=v=0 other boundaries
*  Boundary conditions pressure correction
*  p=0 @ y=2 and p_x=0 at x=0,2 and p_y=0 @ y=0.
*  with initial condition u=v=p=0 everywhere.
*/

#include<iostream>  
#include<string>    
#include<math.h>
#include<fstream>  
#include<sstream>
#include<iomanip>
#include <vector>
#define PI 3.14159265 

using namespace std;
/*Prototype functions to write an output file (output)
fill the velocity vector and the pressure field with the initial conditions 
(initial_conditions), apply the boundary conditions (boundary_conditions) 
after to calculate the solution in each step and for calling the differents schemes, 
(abmultistep). The abmultistep function allows to the user select the desired method, in
this case splitting method abstep=1 is available. Finally the pressure function which
computes the pressure correction term for each step in the splitting method.*/
void output(int, int, int, vector< vector<double> > &, int);
void initial_conditions(int, int, int, double, double, vector< vector<double> > &, 
                            vector< vector<double> > &,vector< vector<double> > &);
void boundary_conditions(vector< vector<double> > &, int, int, double, int);
void abmultistep(int, int, int, int, double, double , double, double, 
      double, double, double, vector< vector<double> > &,
      vector< vector<double> > &, vector< vector<double> > &);
void splitting(double, double, double, int, int, int, double, double, double, double, 
      vector< vector<double> > &, vector< vector<double> > &, vector< vector<double> > &);
void pressure(double, double, double, int, int, int, vector< vector<double> > &,
          vector< vector<double> > &, vector< vector<double> > &);

int main()
{
/*Define the initial parameters for the numerical solution, the grid
xl<x<xr, where xl=0 and xr=2, yl<y<yr, where yl=0 and yr=2. The time domain 
for the solution, where to is the initial time and tf is the final time. nx and 
ny are the number of grid points for the square grid in the y an x direction 
respectively. N is the number of time steps and re is the Reynolds number defined 
by rho*U*Lx/mu where Lx=Ly.*/
int xl=0;
double xr=2.0;
double yl=0.0;
double yr=2.0;
int to=0;
double tf=0.1;
double N=100;
double nx=50; 
double ny=50;
double re=10;
/*delta_x,y,t are the size of the space steps in the grid and
the size of the time step respectively. Alpha, beta coefficients 
are the ratio between the time step and the space steps. sigma_x,y
are the ratios between the time step and the square space steps times
the reynolds number. Select the method to compute the numerical 
solution, i.e. abstep =1 splitting method. Define the u,v and p vectors
that will be used to storage the solution of the Navier Stokes Equations.*/
double delta_x=(xr-xl)/(nx);
double delta_y=(yr-yl)/(ny);
double delta_t=(tf-to)/N;
double alpha=delta_t/delta_x;
double beta=delta_t/delta_y;
double sigma_x=delta_t/(delta_x*delta_x*re);
double sigma_y=delta_t/(delta_y*delta_y*re);
cout<<"delta_x= "<<delta_x<<endl;
cout<<"delta_t= "<<delta_t<<endl;
int xteps=int(nx);
int yteps=int(ny);
int tsteps=int(N);
int abstep=1;
vector< vector<double> > u(xteps+1,vector<double> (yteps+1));
vector< vector<double> > v(xteps+1,vector<double> (yteps+1));
vector< vector<double> > p(xteps+1,vector<double> (yteps+1));
abmultistep(abstep, xteps, yteps, tsteps, delta_x, delta_y, 
            delta_t, alpha, beta, sigma_x, sigma_y, u, v, p);
return 0;
}
/*Function to call the splitting method and control different multisteps methods
if they are available.*/
void abmultistep(int abstep, int xteps, int yteps, int tsteps, double delta_x, double delta_y, 
                    double delta_t, double alpha, double beta, double sigma_x, double sigma_y,
         vector< vector<double> > &u,vector< vector<double> > &v, vector< vector<double> > &p)
{
initial_conditions(xteps, yteps, tsteps, delta_x, delta_y, u, v, p);
switch (abstep)
  {
   case 1:
      splitting(delta_x, delta_y, delta_t, tsteps, xteps, yteps, alpha, beta, sigma_x, sigma_y, u, v, p);
      break;
  }

}
/* Function to compute the initial conditions at time t=0,
the initial conditions u=v=p=0 everywhere. After initialize the
u, v an p vectors, the boundary conditions (B.C) are applied and the results 
printed by the output function.
*/
void initial_conditions(int xteps, int yteps, int tsteps, double delta_x, 
        double delta_y, vector< vector<double> > &u, vector< vector<double> > &v, 
                                                    vector< vector<double> > &p)
{
int t=0;
for (int i=0; i<=xteps; i++)
  {
   for (int j=0; j<=yteps; j++)
     {
      u[i][j]=0;
      v[i][j]=0; 
      p[i][j]=0;
     }
  }
boundary_conditions(u, yteps, xteps, delta_x, 1);
boundary_conditions(v, yteps, xteps, delta_x, 2);
boundary_conditions(p, yteps, xteps, delta_x, 3);
output(xteps, yteps, t, u, 1);
output(xteps, yteps, t, v, 2);
output(xteps, yteps, t, p, 3);
}
/*The Boundary conditions function apply Dirichlet homogeneus boundary conditions
in each time step (delta_t) of the solution for the u and v vectors except in the
top for u. Dirichlet and Neumann boundary condtions are applied to the pressure field
vector. A special field number (fieldv) allows to the function to know which boundary
conditions apply. For instance u=1, v=2 and p=3, with this number the field vector
fv is filled with the respective boundary conditions.
Velocity boundary conditions:
square wave at the top 
u=0 if 0<x<dx and 2-dx<x<dx @y=2  
u=2 if dx<=x<=2-dx @y=2  
u=v=0 other boundaries
Pressure boundary conditions
p=0 @ y=2 and p_x=0 at x=0,2 and p_y=0 @ y=0.
*/
void boundary_conditions(vector< vector<double> > &fv, int yteps, int xteps, double dx, int fieldv)
{
for (int i=0; i<=xteps; i++)
  {
   for (int j=0; j<=yteps; j++)
     {
      switch(fieldv)
       {
        case 1:
            if ((j==yteps)&&(dx<=i*dx && i*dx<=(xteps-1)*dx))
               {
                fv[i][j]=2;
               }
            else if ((j==yteps)&&((0<=i*dx && i*dx<dx)||((xteps-1)*dx<i*dx && i*dx<=xteps)))
               {
                fv[i][j]=0;
               }
            else if ((i==0 || i==xteps) || j==0)
               {
                fv[i][j]=0;
               }
            break;
        case 2:
            if ((i==0 || i==xteps)||(j==0 || j==yteps))
               {
                fv[i][j]=0;
               }
            else
              continue;
            break;
        case 3:
             if (j==0)
               {
                fv[i][j]=fv[i][j+1];
               }
             else if (j==yteps)
               {
                fv[i][j]=0;
               } 
             else if (i==0)
               {
                fv[i][j]=fv[i+1][j];
               }
             else if (i==xteps)
               {
                fv[i][j]=fv[i-1][j];
               }
            else
              continue;
            break;
       }
     }  
  }   
}
/*Compute the solution for the cavity flow using the splitting method in
four steps, an example of the method in the x direction
uhat=u[n]-(u[n]*u_x[n]+v[n]*u_y[n])                  \\predictor step-1 AB one step in time
nabla2p=div.uhat                                     \\compute pressure corrector term step-2
u2hat=uhat[n]-p_x[n]*dt                              \\add pressure correction   step-3
u[n+1]=u2hat+1/Re_l*(u_{xx}[n]+u_{yy}[n])*dt         \\final velocity  step-4 AB one step in time
To increase the rate of convergence in the convective term a predictor corrector
method is used, Euler Trapezoidal, one time step in the x direction is
//-----step 1/2
utilde[n+1][p]=u[n]+dt*f(u[n],t[n])                     \\Predictor
uhat=u[n]+0.5*dt*(f(utilde[n+1][p],t[n+1])+f(u[n],t[n])) \\corrector
where f(u[n],t[n])=u[n]*u_x[n]+v[n]*u_y[n]
The central difference in space scheme is used for the first and second derivatives
for the u, v and p fields.
*/
void splitting(double dx, double dy, double dt, int tsteps, int xteps, int yteps, double alpha, 
              double beta, double sigma_x, double sigma_y, vector< vector<double> > &un, 
               vector< vector<double> > &vn, vector< vector<double> > &pn)
{
vector< vector<double> > utilde(xteps+1,vector<double> (yteps+1));
vector< vector<double> > uhat(xteps+1,vector<double> (yteps+1));
vector< vector<double> > u2hat(xteps+1,vector<double> (yteps+1));
vector< vector<double> > vtilde(xteps+1,vector<double> (yteps+1));
vector< vector<double> > vhat(xteps+1,vector<double> (yteps+1));
vector< vector<double> > v2hat(xteps+1,vector<double> (yteps+1));

 for (int t=0; t<tsteps; t++)
   {
//----------step 1/2
    for (int i=1; i<xteps; i++)
      {
       for (int j=1; j<xteps; j++)
          {
           utilde[i][j]=un[i][j]-0.5*alpha*(un[i][j]*(un[i+1][j]-un[i-1][j])+
                                      vn[i][j]*(un[i][j+1]-un[i][j-1]));
           vtilde[i][j]=vn[i][j]-0.5*alpha*(un[i][j]*(vn[i+1][j]-vn[i-1][j])+
                                      vn[i][j]*(vn[i][j+1]-vn[i][j-1]));
          }
      } 
    boundary_conditions(utilde, yteps, xteps, dx, 1);
    boundary_conditions(vtilde, yteps, xteps, dx, 2);
//------------step 1
    for (int i=1; i<xteps; i++)
      {
       for (int j=1; j<xteps; j++)
          {
           uhat[i][j]=un[i][j]-0.25*alpha*(utilde[i][j]*(utilde[i+1][j]-utilde[i-1][j])+
                                      vtilde[i][j]*(utilde[i][j+1]-utilde[i][j-1])+ 
                                           un[i][j]*(un[i+1][j]-un[i-1][j])+ 
                                           vn[i][j]*(un[i][j+1]-un[i][j-1]));

           vhat[i][j]=vn[i][j]-0.25*alpha*(utilde[i][j]*(vtilde[i+1][j]-vtilde[i-1][j])+
                                      vtilde[i][j]*(vtilde[i][j+1]-vtilde[i][j-1])+
                                           un[i][j]*(vn[i+1][j]-vn[i-1][j])+ 
                                           vn[i][j]*(vn[i][j+1]-vn[i][j-1]));
          }
      } 
//----------intermediate step update pressure(step 2)
    pressure(dx, dy, dt, t, xteps, yteps, uhat, vhat, pn);
//----------step 3 and 4
    for (int i=1; i<xteps; i++)
      {
       for (int j=1; j<xteps; j++)
          {
           u2hat[i][j]=uhat[i][j]-(0.5*alpha)*(pn[i+1][j]-pn[i-1][j]);
           un[i][j]=u2hat[i][j]+sigma_x*(un[i+1][j]-2*un[i][j]+un[i-1][j])+
                                sigma_y*(un[i][j+1]-2*un[i][j]+un[i][j-1]);
           v2hat[i][j]=vhat[i][j]-(0.5*beta)*(pn[i][j+1]-pn[i][j-1]);
           vn[i][j]=v2hat[i][j]+sigma_x*(vn[i+1][j]-2*vn[i][j]+vn[i-1][j])+
                                sigma_y*(vn[i][j+1]-2*vn[i][j]+vn[i][j-1]);
          }
      }
//-----------end time step
    boundary_conditions(un, yteps, xteps, dx, 1);
    boundary_conditions(vn, yteps, xteps, dx, 2);
    output(xteps, yteps, t+1, un, 1);
    output(xteps, yteps, t+1, vn, 2);
   }
} 
/*
Pressure function compute the pressure correction using Successive Overrelaxation 
(SOR) method
pn[i][j] =(1-w)*p_old[i][j]+w*0.25*
              (p_old[i+1][j]+pn[i-1][j]+p_old[i][j+1]+pn[i][j-1]+div.uhat/dt)
where pn is the new solution of the pressure an w is the over relaxation parameter.
diff parameter storages the error in the SOR method.
*/
void pressure(double delta_x, double delta_y, double delta_t, int t, 
                  int xteps, int yteps, vector< vector<double> > &un, 
            vector< vector<double> > &vn, vector< vector<double> > &pn)
{
double dx=delta_x;
double dy=delta_y;
double dx2=dx*dx;
double dy2=dy*dy;
double dt=delta_t;
t++;
double diff=1;
double w=1.68;  
int iterations = 0;
vector< vector<double> > p_temp(xteps+1,vector<double> (yteps+1));
//-----Succesive Over Relaxation Method SOR
while( (iterations <= 60) && ( diff > 0.00001) )
  {
   p_temp = pn; 
   diff = 0.;  
   for (int i=1; i<xteps; i++)
    {
     for(int j=1; j<yteps; j++)
       {
        pn[i][j] =(1-w)*p_temp[i][j]+w*0.25*(p_temp[i+1][j]+pn[i-1][j]+p_temp[i][j+1]+pn[i][j-1]
         -(0.5*dx/dt)*((un[i+1][j]-un[i-1][j])+(vn[i][j+1]-vn[i][j-1])));
        diff=diff+fabs(p_temp[i][j]-pn[i][j]);
       }
    }
   iterations++;
   diff /= pow((xteps+1.0),2.0);
  } 
boundary_conditions(pn, yteps, xteps, dx, 3);
output(xteps, yteps, t, pn, 3);
}
/*Function to export the results to an output file,
the results are saved in three different folders
depending of the field variable U, V and P. 
*/
void output(int xteps, int yteps, int t, vector< vector<double> > &un, int fieldv)
{
int time=t;
time++;
string s, dir, ext, filename;
ext=".data";
if(fieldv==1)
   {
    dir="U/2Dlinear_";
   }
else if(fieldv==2)
   {
    dir="V/2Dlinear_";
   }
else if(fieldv==3)
   {
    dir="P/2Dlinear_";
   }
else
   {
    cout<<"error in directory"<<endl;
   }
stringstream convert; 
convert<<time;
s = convert.str();
filename=dir+s+ext;
ofstream myfile (filename.c_str());
 if (myfile.is_open())
  {
   for (int i=0; i<=yteps; i++)
       {
        for (int j=0; j<=xteps; j++)
           {
             myfile<<setprecision(10)<<un[i][j]<<' ';
           }
           myfile <<"\n"; 
       }
  }    
 else cout <<"Unable to open file";
myfile.close();
}
