%--Author: Juan Carlos Cuevas Bautista
%--Date: 11/20/2014
%--Program to compute the numerical solution to Burger's Equation
%--u_t+u*u_x=nu*u_xx
%--Using the Chebysev-gauss-lobato points and the maxima 
%--cardinal function with domain x E (-1,1)
%--Dirichlet Boundary conditions u(1)=u(-1)=0
%--Initial condition u(x,0)=-sin(pi*x)
clear all
close all
graphics_toolkit ("gnuplot")
%--Define the size of the mesh; M number of gridpoints,
%--N time steps, ua and ub boundary conditions,
%--Chebysev-Gauss Lobato points j=0,1,....,M, 
%--the time interval tf-to and diffusion coefficient nu.
%--I the identity matrix, xj is the interpolation points
%--and dt the size of time step.
M=140;
N=8;
j=(0:M);
ua=0;
ub=0;
nu=0.01;
to=0;
tf=0.2;
I=diag(ones(M-1,1),0);  
delta_t=(tf-to)/N;
xj=cos(j.*pi/M)';
%--Compute the initial condition uo and reshape
%--the matrix u(xj,t) to compute just in the inner gridpoints 2:M
%--star timing
tic,
uo=-sin(pi*xj);
u=uo(2:M,1);
%--To compute the numerical solution for the Burger's equation
%--is necessary to use two indexes, for the grid-points and the
%--time respectively, since this requires the use of loops and since 
%--Matlab is not a compiled language then some loops are avoided
%--by vectorize the operations.
%--function diff1 compute the first chebysev derivative
%--for the maxima cardinal function
%--psi_k=(-1)^(k+1)*(1-x^2)T'_N/(c_k*N^2*(x-x_k))
%--where psi_k(x_j)=0 if j~=k and 1 if j=k
%--function taken from "Spectral Methods in Matlab" Author:Trefethen.
function [D,xj]=diff1(space_steps,grid)
M=space_steps;
xj=grid;
if M==0
 D=0;
 xj=0;
end
c=[2; ones(M-1,1);2].*(-1).^(0:M)';
X=repmat(xj,1,M+1);
dX=X-X';
D=(c*(1./c)')./(dX+(eye(M+1))); 
D= D-diag(sum(D')); 
end
%--Compute the first chebyshev differentiation matrix d1
%--reshape it to include just the inner points DM
%--do the same for the second differentiation matrix D2M.
%--Burger's equation for the Chebysev Gauss-Lobatto points
%--u_(M,t)+u_M*DM*u_M=nu*D2M*u_M
%--solve the system Ax=b => u=A\b
[d1,xpoints]=diff1(M,xj);
DM=d1(2:M,2:M);
D2=d1*d1;
D2M=D2(2:M,2:M);
%--First time step t=1 of Adam Bashforth 
%--and Crank-Nicolson is computed here.
A=I-0.5*nu*delta_t*D2M;
B=I+0.5*nu*delta_t*D2M;
u(:,2)=A\(B*u(:,1)-delta_t*(u(:,1).*(DM*u(:,1))));
%--n+1 time steps for Adams-Basforth and Crank-Nicolson are computed here.
%--Solve the system Ax=b => u=A\b, then the matriz solution usol 
%--is constructed along with the two boundary conditions matrices 
%--ua and ub respectively.
for t=2:N
  u(:,t+1)=A\(B*u(:,t)-0.5*delta_t*...
             (3*(u(:,t).*(DM*u(:,t)))-(u(:,t-1).*(DM*u(:,t-1)))));
end
ua=ua*ones(1,N+1);
ub=ub*ones(1,N+1);
usol=[ua;u;ub];
%--Stop timing
toc,
%%%%%%%-----Post Processing
%--2D plot
h1=figure(1);
plot(xj,usol);
xlabel('x_j','fontsize',12);
ylabel('u(x,t)','fontsize',12);
saveas(h1,'gauss_lobato.epsc2');
