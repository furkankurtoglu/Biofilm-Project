close all
clear
clc



cd C:\Users\Furkan\Documents\GitHub\Biofilm-Project\PhysiCell_Folder
cd output

% 
% cd C:\Users\Furkan\Desktop\output_100_um_tumor_1e-6_degradation

%%
% Obtaining names of mat files
s = what;
MatFiles = s.mat;
MEMatFiles = MatFiles(contains(MatFiles,'micro'));
MEMatFiles(1) = [];



%%
for i = 1:length(MEMatFiles)
    
    load(MEMatFiles{i});
    O2 = multiscale_microenvironment(5,:);
    Glucose = multiscale_microenvironment(7,:);
    ECM = multiscale_microenvironment(6,:);
%     Autoinducer = multiscale_microenvironment(8,:);
    O2 = reshape(O2, [25 25]);
    Glucose = reshape(Glucose, [25 25]);
%     Autoinducer = reshape(Autoinducer, [25 25]);
    ECM = reshape(ECM, [25 25]);
    
    XPos = multiscale_microenvironment(1,:);
    YPos = multiscale_microenvironment(2,:);
    XPos = reshape(XPos, [25 25]);
    YPos = reshape(YPos, [25 25]);
    
    
    ax2 = subplot(2,2,1);
    contourf(XPos,YPos,O2,'linecolor','none');
    caxis(ax2,[0 20]);
    title('Oxygen');
    colorbar('eastoutside');
    axis image;
    
    ax3 = subplot(2,2,2);
    contourf(XPos,YPos,Glucose,'linecolor','none');
   % caxis(ax3,[-0.001 1]);
    colorbar('eastoutside');
    set(gca,'color',[0.2422,0.1504,0.6603]);
    title('Glucose');
    axis image;
    
    
    ax1 = subplot(2,2,3);
    contourf(XPos,YPos,ECM,'linecolor','none');
 %   caxis(ax3,[-0.001 1]);
    colorbar('eastoutside');
    set(gca,'color',[0.2422,0.1504,0.6603]);
    title('ECM');
    axis image;    
    
%     ax4 = subplot(2,2,4);
%     contourf(XPos,YPos,Autoinducer,'linecolor','none');
%  %   caxis(ax3,[-0.001 1]);
%     colorbar('eastoutside');
%     set(gca,'color',[0.2422,0.1504,0.6603]);
%     title('Autoinducer');
%     axis image;    
    
    
    BigTitle = num2str(i/24);
    suptitle(['Time = ',BigTitle,' days']);
    pause(0.05);
end

%%
cd ..
