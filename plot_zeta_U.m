%%
clim1=[0 0.02];
clim2=[0 0.02];
clim3=[0 0.02];
clim4=[-0.9 0.9];
clim5=[-0.9 0.9];
clim6=[-0.9 0.9];

figure (10)
set(gcf,'Position',[703 200 986 489])

subplot('position',[0.09 0.55 0.2 0.4])
T_sample=60;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
xticks([])
%caxis(clim1)
ylabel('$y$','interpreter','latex')
text(0.1,-0.05,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
text(-0.2,-0.05,'$\zeta$','interpreter','latex')

subplot('position',[0.32 0.55 0.2 0.4])
T_sample=80;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
xticks([])
%caxis(clim2)
yticks([])
text(0.1,-0.05,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
 %cl.Position=[0.925 0.73 0.0342 0.2424];

subplot('position',[0.55 0.55 0.2 0.4])
T_sample=100;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
%caxis(clim3)
xticks([])
yticks([])
text(0.1,-0.05,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
% cl.Position=[0.47 0.42 0.0342 0.2424];


subplot('position',[0.78 0.55 0.15 0.4])
T_sample=100;
[~,index]=min(abs(T_sample-T));
plot(U_mean(index,:),flip(y_arr))
yticks([])
xticks([])
ylim([0 2*pi/m])
xlim([-0.04 0.02])
text(-0.05,1.3,'$U$','Interpreter','latex')
text(-0.005,1.3,['$t=$',num2str(T_sample)],'Interpreter','latex')

subplot('position',[0.09 0.09 0.2 0.4])
T_sample=120;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
%caxis(clim4)
text(0.1,-0.05,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
ylabel('$y$','interpreter','latex')
 %cl.Position=[0.925 0.42 0.0342 0.2424];
%text(-2,-1,'$\zeta$','Interpreter','latex')

subplot('position',[0.32 0.09 0.2 0.4])
T_sample=140;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
%caxis(clim5)
yticks([])
xlabel('$x$','interpreter','latex')
text(0.1,-0.05,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
 %cl.Position=[0.47 0.1015 0.0342 0.2424];


subplot('position',[0.55 0.09 0.2 0.4])
T_sample=160;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
cl=colorbar('eastoutside')
%caxis(clim6)
text(0.1,-0.05,['$t=$',num2str(T_sample)],'interpreter','latex')
 %cl.Position=[0.925 0.1015 0.0342 0.2424];

subplot('position',[0.78 0.09 0.15 0.4])
T_sample=160;
[~,index]=min(abs(T_sample-T));
plot(U_mean(index,:),flip(y_arr))
yticks([])
ylim([0 2*pi/m])
text(-0.005,1.3,['$t=$',num2str(T_sample)],'Interpreter','latex')
xlabel('$U$','interpreter','latex')
