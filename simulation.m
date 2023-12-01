tic

% flow parameters
k=16;
m=5;
eta=1e-4;
mu=0;
nu=1e-4;
sigma=0.05;
beta=5;
B0=1e-2;

%Number of Fourier components and grid points
Nx=256;
Ny=256;

%Grid points in x and y
x_arr=linspace(0,2*pi/k,Nx+1);
y_arr=linspace(0,2*pi/m,Ny+1);
x_arr=x_arr(1:end-1);
y_arr=y_arr(1:end-1);
dx=x_arr(2)-x_arr(1);
dy=y_arr(2)-y_arr(1);

[X,Y]=meshgrid(x_arr,y_arr);

%Wavenumbers in x and y
k1_arr=0:Nx-1;
k2_arr=0:Ny-1;
[K1,K2]=meshgrid(k1_arr,k2_arr);

%Assign the wavenumbers for Fourier transforms used for second-order derivatives with de-aliasing
K1_new=K1;
K2_new=K2;

if rem(Nx,2)==0
   K1_new(:,Nx/2+2:end)=K1(:,Nx/2+2:end)-Nx;
else
   K1_new(:,(Nx+1)/2+1:end)=K1(:,(Nx+1)/2+1:end)-Nx;
end

if rem(Ny,2)==0
   K2_new(Ny/2+2:end,:)=K2(Ny/2+2:end,:)-Ny;
else
   K2_new((Ny+1)/2+1:end,:)=K2((Ny+1)/2+1:end,:)-Ny;
end

%The time step and the number of steps
dt=0.005;
N_t=250000/5;

%Generate the real and imaginary part of white noise; 
%W represents a Broanian motion
%Remove this part if one loads the sample white noise from the data files
obj = bm(0, 1);
[W, T] = obj.simByEuler(N_t, 'DeltaTime', dt);
e=ones(N_t+1,1);
M1=spdiags([-e e]/(dt),[0,1],N_t+1,N_t+1);
noise_r=M1*W;

obj = bm(0, 1);
[W, T] = obj.simByEuler(N_t, 'DeltaTime', dt);
e=ones(N_t+1,1);
M1=spdiags([-e e]/(dt),[0,1],N_t+1,N_t+1);
noise_i=M1*W;

%The initial wave field; b1 and b2 are the magnetic field in x and y
psi=real(0.0001*exp(1i*m*Y));
b1=0*Y;
b2=0*Y;

%Computing the gradients in x and y
%fft_df_2d is a function file that computes the derivative in x and y using
%Fourier transform
[psi_x,  psi_y]=fft_df_2d(Nx, Ny, k, m, psi);

%The initial mean velocity and mean field
U_mean=zeros(length(T),Ny);
B_mean=zeros(length(T),Ny);
U_mean_rms=zeros(1,length(T));
B_mean_rms=zeros(1,length(T));

U_mean(1,:)=-k/2/pi*sum(psi_y,2)*dx;
U_mean_rms(1)=sqrt(m/2/pi*sum((U_mean(1,:)).^2)*dy);

%sample time steps to collect data
n_skip=1000;
T_save=T(1:n_skip:end);
zeta_save=zeros(length(T_save),Ny,Nx);
j_save=zeros(length(T_save),Ny,Nx);
ns=2;

for n=2:length(T)
%The vorticity field.
%Laplacian_fft is a function file that computes the Laplacian using Fourier
%transforms
zeta=Laplacian_fft(k,m,Nx,Ny,psi);

%The forcing with the complex white noise
F=2*real(sigma*ones(Ny,1)*exp(1i*k*x_arr)*(noise_r(n)+1i*noise_i(n))/sqrt(2));

%Computing the gradients in x and y
%fft_df_2d is a function file that computes the derivative in x and y using
%Fourier transform
[psi_x,  psi_y]=fft_df_2d(Nx, Ny, k, m, psi);
[psi_xx,  ~]=fft_df_2d(Nx, Ny, k, m, psi_x);
[psi_xy, psi_yy]=fft_df_2d(Nx, Ny, k, m, psi_y);

[zeta_x, zeta_y]=fft_df_2d(Nx, Ny, k, m,zeta);
[b1_x, b1_y]=fft_df_2d(Nx, Ny, k, m, b1);
[b2_x, b2_y]=fft_df_2d(Nx, Ny, k, m, b2);
j=b2_x-b1_y;
[j_x, j_y]=fft_df_2d(Nx, Ny, k, m, j);

%The Fourier transforms
%f_NT represents the Fourier transform of the nonlinear terms
f_zeta=fft2(zeta);
f_NT=fft2( psi_y.*zeta_x-psi_x.*zeta_y-beta*psi_x+F +j_x.*(b1+B0)+j_y.*b2 );

%Computing the Fourier transform of vorticity at the next step using
%Crank-Nicolson scheme. 
f_zeta=f_zeta.*(1-mu/2*dt+nu/2*( -(k*K1_new).^2-(m*K2_new).^2)*dt)./(1+mu/2*dt-nu/2*(  -(k*K1_new).^2-(m*K2_new).^2)*dt)+dt*f_NT./(1+mu/2*dt+nu/2*(  (k*K1_new).^2+(m*K2_new).^2)*dt);

%Performing the inverse Fourier transform to find the vorticity
zeta=ifft2(f_zeta);

%Performing the inverse Laplacian to find the stream function
psi=Laplacian_fft_inv(k,m,Nx,Ny,zeta);

%Similar computation method for the magnetic field
f_b1=fft2(b1);
f_NT=fft2(-psi_yy.*b2+psi_y.*b1_x-psi_xy.*(B0+b1)-psi_x.*b1_y);
f_b1=f_b1.*(1+eta/2*( -(k*K1_new).^2-(m*K2_new).^2)*dt)./(1-eta/2*(  -(k*K1_new).^2-(m*K2_new).^2)*dt)+dt*f_NT./(1+eta/2*(  (k*K1_new).^2+(m*K2_new).^2)*dt);
b1=ifft2(f_b1);

f_b2=fft2(b2);
f_NT=fft2(psi_xy.*b2+psi_y.*b2_x+psi_xx.*(B0+b1)-psi_x.*b2_y);
f_b2=f_b2.*(1+eta/2*( -(k*K1_new).^2-(m*K2_new).^2)*dt)./(1-eta/2*(  -(k*K1_new).^2-(m*K2_new).^2)*dt)+dt*f_NT./(1+eta/2*(  (k*K1_new).^2+(m*K2_new).^2)*dt);
b2=ifft2(f_b2);

%The mean velocity and the mean magnetic field
U_mean(n,:)=-k/2/pi*sum(psi_y,2)*dx;
U_mean_rms(n)=sqrt(m/2/pi*sum((U_mean(n,:)).^2)*dy);

B_mean(n,:)=k/2/pi*sum(b1,2)*dx;
B_mean_rms(n)=sqrt(m/2/pi*sum((B_mean(n,:)).^2)*dy);

%In case of numerical instability which makes the mean velocity unreasonably large, let the computation stop
if abs(U_mean_rms(n))>1
figure(1)
plot(T(1:n-5),transpose(U_mean_rms(1:n-5)))
drawnow    
   break 
end

%Show sample results during computation
if ceil(n/1500)==n/1500

figure (2)
imagesc(x_arr,y_arr,real(zeta))
title(['$\zeta$ at $t=$',num2str(T(n))],'interpreter','latex')
colorbar
drawnow

figure(3)
plot(T(1:n),U_mean_rms(1:n))
drawnow
set(gca, 'YScale', 'log')
xlabel('$t$','interpreter','latex')
ylabel('$U_\mathrm{rms}$','Interpreter','latex')

figure(4)
plot(T(1:n),B_mean_rms(1:n))
drawnow
set(gca, 'YScale', 'log')
xlabel('$t$','interpreter','latex')
ylabel('$B_\mathrm{rms}$','Interpreter','latex')
ylim([1e-10,0.01])

figure(5)
imagesc(T(1:n),y_arr, U_mean(1:n,:)' )
xlabel('$t$','interpreter','latex')
ylabel('$y$','interpreter','latex')
text(100,-0.05,'$U(y,t)$','interpreter','latex')
text(-27,-0.05,'$(c)$','interpreter','latex')
set(gcf,'Position',[272 380 490 181])
colorbar

figure(6)
imagesc(T(1:n),y_arr, B_mean(1:n,:)' )
xlabel('$t$','interpreter','latex')
ylabel('$y$','interpreter','latex')
text(100,-0.05,'$B(y,t)$','interpreter','latex')
set(gcf,'Position', [271 100 492 202])
text(-27,-0.05,'$(e)$','interpreter','latex')
colorbar

figure(7)
imagesc(x_arr,y_arr,real(j))
title(['$j$ at $t=$',num2str(T(n))],'interpreter','latex')
colorbar
drawnow

end

if n==(ns-1)*n_skip+1  
   zeta_save(ns,:,:)=zeta;
   j_save(ns,:,:)=j;
   ns=ns+1;
end

end

toc
