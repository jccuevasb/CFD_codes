%--Author: Juan Carlos Cuevas Bautista
%--Date: 11/07/2014
%--Program to compute the numerical solution to Burger's Equation
%--u_t+u*u_x=nu*u_xx
%--By Galerking Fourier Method
%--With domain x E (-1,1)
%--Boundary conditions u(1)=u(-1)=0
%--Burger's equation is re-scaled using the transformation x'=pi*x
%--u_t+pi*u*u_x'=nu*pi^2*u_xx'.
clear all
close all
graphics_toolkit ("gnuplot")
%--Define the size of the mesh; M gridpoints
%--N time steps, xl and xr domain of the function,
%--the time interval tf-to and diffusion coefficient nu
M=8192;  
N=8; 
xl=-pi;
xr=pi;
to=0;
tf=0.2;
nu=0.01;
%--Compute the time vector for plotting and
%--the fourier wavenumber k E 0,1,...,M/2,-M/2+1, -M/2+2,...,-1.
%--Matlab storage positive and negative wavenumbers.
%--Define the size of the time steps dt.
%--The grid xj=[-pi,dx,pi], where xj=pi*x and x E [-1,1]
k=2*pi/(xr-xl)*[0:M/2-1,0,-M/2+1:-1];
delta_x=(xr-xl)/M;
delta_t=(tf-to)/N;
time=[0:delta_t:tf];
xj=xl+(0:M-1)'*delta_x;
%xj=xl+(0:M)'*delta_x;

%--star timing
tic,
%--Compute the initial condition 
%--u(x,0)=-sin(xj) 
%--Transform uo to fourier domain u_hat using Fast Fourier Transform fft
%--build-it function in matlab.
uo=-sin(xj);
u_hat(:,1)=fft(uo);
%--The convolution in the fourier domain of the convection term F(uu_x)
%--is computed by w_hat. The arguments are uhat and the wavenumber k. 
%--First compute uhat and the derivative u_x=vhat
%--in the fourier domain, then transform to physical domain using ifft
%--buil-it Matlab function and compute the simple product u*u_x in the
%--physical domain, finally return the result of the product in the
%--fourier domain w_hat.  
function w_hat=uux_hat(u,wavenumber)
k=wavenumber;
u_hat=u;
v_hat=i*k'.*u_hat;
wj=ifft(u_hat).*ifft(v_hat);
w_hat=fft(wj);
end
%--To solve the burger's equation, Crank-Nicolson method for linear
%--term nu*u_xx and second order Adam-Bashforth method for nonlinear
%--term u*u_x are used. First time step t=1 of Adam Bashforth 
%--and Crank-Nicolson is computed here. Burger's equation in the 
%--Fourier domain is 
%--uhat_t+pi*F(u*u_x')+pi^2*k^2*nu*u_hat.
%--First compute u_hat(n+1) in the step n+1, then update w_hat(n+1). 
w_hat(:,1)=uux_hat(fft(uo),k);
for t=1:1
  for m=1:length(k)
     u_hat(m,t+1)=(u_hat(m,t)*(1-0.5*nu*pi^2*k(m)^2*delta_t)-...
           pi*delta_t*w_hat(m,t))/(1+0.5*nu*pi^2*k(m)^2*delta_t);
  end
  w_hat(:,t+1)=uux_hat(u_hat(:,t+1),k);
end
%--n+1 time steps for Adams-Basforth and Crank-Nicolson are computed here.
%--First compute u_hat(n+1) in the step n+1, then update w_hat(n+1).
for t=2:N
  for m=1:length(k)
     u_hat(m,t+1)=(u_hat(m,t)*(1-0.5*nu*pi^2*k(m)^2*delta_t)-...
           0.5*pi*delta_t*(3*w_hat(m,t)-w_hat(m,t-1)))/(1+0.5*nu*pi^2*k(m)^2*delta_t);
  end
  w_hat(:,t+1)=uux_hat(u_hat(:,t+1),k);
end  
%--Transform the spectral solution u_hat into the physical space uj and stop
%--time counting for execution program.
uj=ifft(u_hat);
toc,
%%%%%%%-----Post Processing
%--2D plot
%--3D plot
%vj=ifft(v_hat);
h1=figure(1);
alp=1/pi;  
plot(xj*alp,uj)
xlabel('x_j','fontsize',12);
ylabel('u(x,t)','fontsize',12);
saveas(h1,'fourier_galerkin.epsc2');
%hold on
%plot(xj,vj)
%plot(xj,wj,'rc')
%plot(xj,sin(xj).*cos(xj),'r--')
%hold off
%mesh(xj,time,real(uj'))
