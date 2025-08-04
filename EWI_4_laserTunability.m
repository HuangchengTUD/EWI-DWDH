clear
clc
close all

mm = 1e-3; um = 1e-6; nm = 1e-9;
set(0,'defaultAxesFontName', 'times new roman','defaultAxesFontSize',16);
set(0,'defaultTextFontName', 'times new roman','defaultTextFontSize',16);
%%
lambda1 = 632.992*nm;   lambda2 = 633.069*nm;
dx = 3.45e-6; dy = 3.45e-6; % pixel size
datapath = 'Experiments\SteppedSamples\';

fieldRef = {'field_mirror_hn.mat','field_mirror_ld.mat'};
if exist ([datapath,fieldRef{1}],"file") == 2
    load([datapath,fieldRef{1}]); load([datapath,fieldRef{2}]);
    disp('Reference field are exist and loaded.')
else
    disp('Reference field are not exist!')
    [holo_mir,M,N,X,Y] = funcs.holo_read([datapath,'mirror.bmp']);
    [field_mirror] = funcs.manual_get_plus1(holo_mir,M,N,X,Y,dx,dy,2); % select 2nd quadrant

    field_mirror_hn = field_mirror{1};    field_mirror_ld = field_mirror{2};
end
%%
Files = dir(fullfile(datapath,'*grooves_T085_A87_001*.bmp'));
loopnum = length(Files);
bar = waitbar(0,'Loading your data');

for kk = 1:loopnum
    tic;
    raw_img = Files(kk).name;
    [holo,M,N,X,Y] = funcs.holo_read([datapath,raw_img]);

    if kk == 1
        disp(['The raw hologram is : ', raw_img]);
        [field_pls1,filter] = funcs.manual_get_plus1(holo,M,N,X,Y,dx,dy,2);

        slopePhase = exp(1i*2*pi*((+4000)*X*dx + (-8000)*Y*dx)); % for grooves

        kz1 = 2*pi*sqrt((1/lambda1).^2-(X/M/dx).^2-(Y/N/dy).^2);
        kz2 = 2*pi*sqrt((1/lambda2).^2-(X/M/dx).^2-(Y/N/dy).^2);
        ForwardPropagate = @(field,dis,kz) ifft2(fftshift(exp(1i*kz*dis)).*fft2(field));

        phiHN = zeros(N,M,loopnum); phiLD = zeros(N,M,loopnum);
        %%%%%%%%%%%  corp a small area for saving memory   %%%%%%%%%%%
%         phiHN = zeros(600,1500,loopnum); phiLD = zeros(600,1500,loopnum);
    else
        SPE = funcs.FT(holo);
        field_pls1{1} = fftshift(ifft2(fftshift(filter{1}.*SPE)));
        field_pls1{2} = fftshift(ifft2(fftshift(filter{2}.*SPE)));
    end

    % 1. compensate the fields on camera plane
    field_camPlane_HN = field_pls1{1}./exp(1i*angle(field_mirror_hn)).*slopePhase;
    field_camPlane_LD = field_pls1{2}./exp(1i*angle(field_mirror_ld)).*slopePhase;

    % 2. Back propagate
    dis = 45.50*mm; % steps: 40.50*mm grooves: 45.5
    field_onfocus_HN = ForwardPropagate(field_camPlane_HN,-dis,kz1);
    field_onfocus_LD = ForwardPropagate(field_camPlane_LD,-dis,kz2);
    % 3. Field registration
    field_onfocus_LD = circshift(field_onfocus_LD,[0 +1]);
    
    % 4. Save the phases for two sources
    phiHN(:,:,kk) = angle(field_onfocus_HN); 
    phiLD(:,:,kk) = angle(field_onfocus_LD);

    currentprogress = roundn((kk/loopnum)*100,-1);
    remainingtime = roundn((loopnum-kk)*toc/60,-1);
    barString = sprintf('Current Progress: %.2f%%, Remainging Time: %.2f min :', ...
        currentprogress,remainingtime);
    waitbar(kk/loopnum,bar,barString);
end
I_HN = abs(field_onfocus_HN); I_HN = funcs.nmlz(I_HN);
I_LD = abs(field_onfocus_LD); I_LD = funcs.nmlz(I_LD);


figure(1);colormap jet;
subplot 222;imagesc(angle(field_onfocus_HN./field_onfocus_LD));axis image;colorbar;title('PhiB after Prop.');drawnow;
subplot 221;imagesc(I_LD);axis image;colorbar;title('Amp. for LD');colormap(gca,"gray");drawnow;
subplot 223;imagesc(phiHN(:,:,kk));axis image;colorbar;title('Phase for He-Ne');drawnow;
subplot 224;imagesc(phiLD(:,:,kk));axis image;colorbar;title('Phase for LD');drawnow;

return
%% compensate for slope phase for PHIB
% residual_slope_Beat = exp(1i*2*pi*((+5)*X*dx + (+0)*Y*dx)).*exp(1i*(+0.0));
phiB_prop_0 = exp(1i*(phiHN(:,:,301)))./exp(1i*(phiLD(:,:,301)));
[coeff,z_fit] = funcs.aberfittingbydraw(phiB_prop_0,1,X*2/M,Y*2/N);
residual_slope_Beat = exp(1i*z_fit);

phiB_prop = angle(phiB_prop_0./residual_slope_Beat);
% phiB_prop(I_HN < 0.1 | I_LD < 0.1) = nan;

% two line profiles on phiB
lh_idx = [701 1300 2700 1350]; lv_idx = [1500 601 1550 2000];
figure(1);
subplot 223;imagesc(phiB_prop);axis image;colormap(gca,"hsv");colorbar;title('Beat phase (rad)');drawnow;
    hold on; 
    plot([lh_idx(1),lh_idx(3)],[lh_idx(2),lh_idx(4)],'--r',LineWidth=2);
    plot([lv_idx(1),lv_idx(3)],[lv_idx(2),lv_idx(4)],'--k',LineWidth=2);hold off;

lineh = exp(1i*phiB_prop(lh_idx(2):lh_idx(4),lh_idx(1):lh_idx(3))); 
lineh = angle(mean(lineh,1,"omitnan"));

linev = exp(1i*phiB_prop(lv_idx(2):lv_idx(4),lv_idx(1):lv_idx(3)));
linev = angle(mean(linev,2,"omitnan"));
subplot 224;plot(lineh,'-r','LineWidth',2); axis tight;grid on;grid minor;ylim([-4 4]);
hold on;plot(linev,'-k','LineWidth',2); hold off; 
legend('Horizontal','Vertical','Location','southwest');

return
%% regions selection
[~,~,dim3] = size(phiHN);
for kk = 1 : loopnum
    phiB_s = angle(exp(1i*phiHN(:,:,kk))./exp(1i*phiLD(:,:,kk))./residual_slope_Beat);

    if kk == 1
        [phases,idxs,stds] = funcs.drawSeveral(phiB_s);
    else
        for jj = 1:length(idxs)
            y1 = idxs{jj}(2); y2 = idxs{jj}(2) + idxs{jj}(4); 
            x1 = idxs{jj}(1); x2 = idxs{jj}(1) + idxs{jj}(3);
            values = exp(1i*phiB_s(y1:y2-1,x1:x2-1));
            phases(kk,jj) = angle(mean(values,"all","omitnan"));
            stds(kk,jj) = std(values(:),"omitnan");
        end
    end
end
phaseDiff = wrapTo2Pi(angle(exp(1i*phases)./exp(1i*phases(:,1))));

key = 2;    % 1: height quantification; 2: wavelength calibration
switch key
    case 1
        lambda2 = 633.000*nm;
        Lambda = lambda1*lambda2/abs(lambda2-lambda1)
        stepsheight = phaseDiff/4/pi*Lambda/mm;

        fprintf(['The average height is: \n',num2str(mean(stepsheight)),'\n']);
    case 2
        h = 993.554*um;
        Lambda = 4*pi*h./phaseDiff(:,2);
        lambda2 = Lambda*lambda1./(Lambda-lambda1)/nm;

        fprintf('The average of Lambda is %.2f + %.3f(mm); lambda2 is %.3f + %.1f(pm);\n', ...
            mean(Lambda)/mm,std(Lambda)/mm,mean(lambda2),std(lambda2)*1e3);
end

%% Laser tunability in a range of current (81mA - 95mA)
groupSize = 30;
groupNum = length(lambda2)/groupSize;

% initialize the variables for 'mean' & 'std' 
l2Means = zeros(groupNum, 1);
l2Std = zeros(groupNum, 1);

% plot the raw phase diff
figure(7); set(gcf,'Position',[400 500 550 330]);

% Through each groups
for kk = 1:groupNum
    % take out the ith group data
    groupidx1 = (kk-1)*groupSize+1; groupidx2 = kk*groupSize;
    x_ = groupidx1:groupidx2;

    groupData = lambda2(groupidx1:groupidx2);
    plot(x_,groupData, '.', 'MarkerSize', 12);hold on;

    % Average & Sigma calculation
    l2Means(kk) = mean(groupData(:));     l2Std(kk) = std(groupData(:));
end
grid on;axis tight;
ylim([632.990 633.15]);

% 绘制每组的均值线
for kk = 1:groupNum
    % 在每个分组范围内绘制均值线
    x_ = (kk-1)*groupSize + 1:kk*groupSize;  % 每组对应的 x 坐标范围
    y_ = l2Means(kk) * ones(size(x_));  % 创建与分组对应的均值线
    plot(x_, y_, '-r', 'LineWidth', 2);  % 绘制均值线，颜色为红色虚线
end

% legend & title & label
xticks(30:30:loopnum);
xlabel('Frame index (#)'); 
ylabel('Wavelength (nm)'); 
hold off;

% set(gca,'XTickLabel',drivecurrent);
driveCurrent = 81:1:95;
figure(8); set(gcf,'Position',[400 100 530 300]);
errorbar(driveCurrent,l2Means,l2Std,'.','MarkerSize',18,'MarkerEdgeColor','auto','LineWidth',2.0);
grid on; axis padded;
% plot wavelength of He-Ne as reference
hold on; plot([81 95],[632.990 632.990],'--r','LineWidth',2.5); hold off

set(gca,'GridAlpha',0.3,'MinorGridAlpha',0.2, ...
    'FontSize',18,'LineWidth',1.0,'XTick',driveCurrent,'YLim',[632.980 633.13],'YTick',632.990:0.03:633.150);
legend(['T = 18.0',char(176),'C'],'Location','northwest','Interpreter','tex');












