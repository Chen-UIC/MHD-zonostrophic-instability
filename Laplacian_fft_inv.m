function Lu_fft=Laplacian_fft_inv(kx,ky,N1,N2,Lu)

L1=2*pi/kx;
L2=2*pi/ky;

fu=fft2(Lu);

k1_arr=0:N1-1;
k2_arr=0:N2-1;
[K1,K2]=meshgrid(k1_arr,k2_arr);

K1_new=K1;
K2_new=K2;

if rem(N1,2)==0
K1_new(:,N1/2+2:end)=K1(:,N1/2+2:end)-N1;
else
K1_new(:,(N1+1)/2+1:end)=K1(:,(N1+1)/2+1:end)-N1;
end

if rem(N2,2)==0
K2_new(N2/2+2:end,:)=K2(N2/2+2:end,:)-N2;
else
K2_new((N2+1)/2+1:end,:)=K2((N2+1)/2+1:end,:)-N2;
end

U=-fu./( (2*pi/L1*K1_new).^2+(2*pi/L2*K2_new).^2 );
U(1,1)=0;

Lu_fft=ifft2(U);
end