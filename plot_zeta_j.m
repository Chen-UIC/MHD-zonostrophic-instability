figure (11)
set(gcf,'Position',[703 200 986 489])

subplot('position',[0.09 0.55 0.22 0.4])
T_sample=60;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
xticks([])
ylabel('$y$','interpreter','latex')
text(0.05,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
text(-0.2,1.3,'$\zeta$','interpreter','latex')
set(gca,'YDir','normal')


text(-0.45, 1.3,'$(a)$','interpreter','latex')

subplot('position',[0.32 0.55 0.22 0.4])
T_sample=80;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
xticks([])
yticks([])
text(0.05,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
set(gca,'YDir','normal')

subplot('position',[0.55 0.55 0.22 0.4])
T_sample=120;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
xticks([])
yticks([])
text(0.05,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
set(gca,'YDir','normal')

subplot('position',[0.78 0.55 0.22 0.4])
T_sample=160;
[~,index]=min(abs(T_sample-T_save));
zeta_sample=zeros(Ny,Nx);
zeta_sample(:,:)=zeta_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[zeta_sample,zeta_sample])
xticks([])
yticks([])
text(0.05,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
set(gca,'YDir','normal')

subplot('position',[0.09 0.09 0.22 0.4])

T_sample=60;
[~,index]=min(abs(T_sample-T_save));
j_sample=zeros(Ny,Nx);
j_sample(:,:)=j_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[j_sample,j_sample])
text(-0.02,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
ylabel('$y$','interpreter','latex')
xlabel('$x$','interpreter','latex')
text(-0.2,1.3,'$j$','interpreter','latex') 
set(gca,'YDir','normal')
text(-0.45, 1.3,'$(b)$','interpreter','latex')

subplot('position',[0.32 0.09 0.22 0.4])
T_sample=80;
[~,index]=min(abs(T_sample-T_save));
j_sample=zeros(Ny,Nx);
j_sample(:,:)=j_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[j_sample,j_sample])
yticks([])
xlabel('$x$','interpreter','latex')
text(0.05,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
cl=colorbar('eastoutside')
set(gca,'YDir','normal')

subplot('position',[0.55 0.09 0.22 0.4])
T_sample=120;
[~,index]=min(abs(T_sample-T_save));
j_sample=zeros(Ny,Nx);
j_sample(:,:)=j_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[j_sample,j_sample])
cl=colorbar('eastoutside')
text(0.05,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
xlabel('$x$','interpreter','latex')
yticks([])
set(gca,'YDir','normal')

subplot('position',[0.78 0.09 0.22 0.4])
T_sample=160;
[~,index]=min(abs(T_sample-T_save));
j_sample=zeros(Ny,Nx);
j_sample(:,:)=j_save(index,:,:);
imagesc([x_arr,2*x_arr],y_arr,[j_sample,j_sample])
cl=colorbar('eastoutside')
text(0.05,1.3,['$t=$',num2str(T_sample)],'interpreter','latex')
xlabel('$x$','interpreter','latex')
yticks([])
set(gca,'YDir','normal')