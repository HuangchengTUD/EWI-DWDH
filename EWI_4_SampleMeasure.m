clear
clc
close all

mm = 1e-3; um = 1e-6; nm = 1e-9;
set(0,'defaultAxesFontName', 'times new roman','defaultAxesFontSize',16);
set(0,'defaultTextFontName', 'times new roman','defaultTextFontSize',16);
%%
lambda1 = 632.992*nm;   lambda2 = 633.036*nm;
dx = 3.45e-6; dy = 3.45e-6; % pixel size
datapath = 'Experiments\SteppedSamples\';

fieldRef = {'field_mirror_hn.mat','field_mirror_ld.mat'};
if exist ([datapath,fieldRef{1}],"file") == 2
    load([datapath,fieldRef{1}]); load([datapath,fieldRef{2}]);
    disp('Reference field are exist and loaded.')
else
    disp('Reference field are not exist!')
    [holo_mir,M,N,X,Y] = funcs.holo_read([datapath,'mirror_T085_A87.bmp']);
    [field_mirror] = funcs.manual_get_plus1(holo_mir,M,N,X,Y,dx,dy,2); % select 2nd quadrant

    field_mirror_hn = field_mirror{1};    field_mirror_ld = field_mirror{2};
end

%%
Files = dir(fullfile(datapath,'*steps_T085_A87_001*.bmp'));
loopnum = length(Files);
bar = waitbar(0,'Loading your data');

for kk = 1:loopnum
    tic;
    raw_img = Files(kk).name;
    [holo,M,N,X,Y] = funcs.holo_read([datapath,raw_img]);

    if kk == 1
        disp(['The raw hologram is : ', raw_img]);
        [field_pls1,filter] = funcs.manual_get_plus1(holo,M,N,X,Y,dx,dy,2);

        slopePhase = exp(1i*2*pi*((+5000)*X*dx + (-3500)*Y*dx)); % for brassSteps

        kz1 = 2*pi*sqrt((1/lambda1).^2-(X/M/dx).^2-(Y/N/dy).^2);
        kz2 = 2*pi*sqrt((1/lambda2).^2-(X/M/dx).^2-(Y/N/dy).^2);
        ForwardPropagate = @(field,dis,kz) ifft2(fftshift(exp(1i*kz*dis)).*fft2(field));

        phiHN = zeros(N,M,loopnum); phiLD = zeros(N,M,loopnum);

    else
        SPE = funcs.FT(holo);
        field_pls1{1} = fftshift(ifft2(fftshift(filter{1}.*SPE)));
        field_pls1{2} = fftshift(ifft2(fftshift(filter{2}.*SPE)));
    end

    % 1. compensate the fields on camera plane
    field_camPlane_HN = field_pls1{1}./exp(1i*angle(field_mirror_hn)).*slopePhase;
    field_camPlane_LD = field_pls1{2}./exp(1i*angle(field_mirror_ld)).*slopePhase;

    % 2. Back propagate
    dis = 56.50*mm; % steps: 56.50*mm grooves: 45.5
    field_onfocus_HN = ForwardPropagate(field_camPlane_HN,-dis,kz1);
    field_onfocus_LD = ForwardPropagate(field_camPlane_LD,-dis,kz2);
    % 3. Field registration
    field_onfocus_LD = circshift(field_onfocus_LD,[+3 +0]);
    
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

figure(1);colormap hsv;
subplot 222;imagesc(angle(field_onfocus_HN./field_onfocus_LD));axis image;colorbar;title('PhiB after Prop.');drawnow;
subplot 221;imagesc(I_LD);axis image;colorbar;title('Amp. for LD');colormap(gca,"gray");drawnow;
subplot 223;imagesc(phiHN(:,:,kk));axis image;colorbar;title('Phase for He-Ne');drawnow;
subplot 224;imagesc(phiLD(:,:,kk));axis image;colorbar;title('Phase for LD');drawnow;

return
%% Strategy 2: Averaging afterwards
residual_slope_Beat = exp(1i*2*pi*((+0)*X*dx + (+0)*Y*dx)).*exp(1i*(+0.0));

phiB_prop = angle(exp(1i*phiHN(:,:,1))./exp(1i*phiLD(:,:,1)).*residual_slope_Beat);
% phiB_prop(I_HN < 0.1 | I_LD < 0.1) = nan;

% two line profiles on phiB
lh_idx = [1001 1300 3000 1350]; lv_idx = [2000 501 2050 2000];
figure(1);
subplot 223;imagesc(phiB_prop);axis image;colorbar;title('Beat phase (rad)');drawnow;
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

%% Strategy 2: Averaging afterwards
[~,~,dim3] = size(phiHN);
for kk = 1 : dim3
    phiB_s = angle(exp(1i*phiHN(:,:,kk))./exp(1i*phiLD(:,:,kk)).*residual_slope_Beat);
    if kk == 1
        % determine the areas
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

phaseDiff = round(angle(exp(1i*phases)./exp(1i*phases(:,1))),4);
phaseDiff = wrapTo2Pi(phaseDiff);

key = 1;    % 1: height quantification; 2: wavelength calibration
switch key
    case 1
        lambda2 = 633.036*nm;
        Lambda = lambda1*lambda2/abs(lambda2-lambda1)
        stepsheight = phaseDiff/4/pi*Lambda/mm;

        fprintf(['The average height is: \n',num2str(mean(stepsheight,1)),'\n']);
    case 2
        h = 993.554*um;
        Lambda = 4*pi*h./phaseDiff(:,2);
        lambda2 = Lambda*lambda1./(Lambda-lambda1)/nm;

        fprintf('The average of Lambda is %.2f + %.3f(mm); lambda2 is %.3f + %.1f(pm);\n', ...
            mean(Lambda)/mm,std(Lambda)/mm,mean(lambda2),std(lambda2)*1e3);
end

%% Reference measurements Comparison ------1
stepResults_1 = load('Experiments\20241024\stepResults_1.mat');
nomi_h = 0:0.5:3.0;
measured_h = mean(stepsheight,1);
sigmah = std(stepsheight,[],1);

h_wli = [0.0528 0.5358 1.0242 1.5119 2.0 2.4873 2.9739] - 0.0528; % Measurement of step height from WLI
std_wli = [0.0042  0.0064  0.0067  0.0066   0.0066  0.0061  0.0066];

figure(9); set(gcf,'Position',[400 100 480 300]);

errorbar(nomi_h,stepResults_1.measured_h,stepResults_1.sigmah,'-o', ...
    'MarkerSize',6,'MarkerFaceColor',[0.8 0.8 0.8],'CapSize',12,'LineWidth',1.5); % Lambda = 4.84;
hold on; errorbar(nomi_h,measured_h,sigmah,'-o', ...
    'MarkerSize',6,'MarkerFaceColor',[0.8 0.8 0.8],'CapSize',12,'LineWidth',1.5); % Lambda = 9.17;
errorbar(nomi_h,h_wli,std_wli,'-o', ...
    'MarkerSize',6,'MarkerFaceColor',[0.8 0.8 0.8],'CapSize',12,'LineWidth',1.5);
hold off;

grid on; axis padded;
set(gca,'FontWeight','bold','FontSize',16,'XTick',nomi_h,'YTick',0:0.5:3.0);

xlabel('Nominal height (mm)'); ylabel('Experiment (mm)');
legend('\Lambda = 4.84 mm','\Lambda = 9.17 mm','WLI','Location','northwest');

% save([datapath,'stepResults_1.mat'],'measured_h',"sigmah");
% save([datapath,'stepResults_2.mat'],'measured_h',"sigmah");
%% Reference measurements Comparison ------  2
stepResults_1 = load('Experiments\20241024\stepResults_1.mat');
stepResults_2 = load('Experiments\20241203\stepResults_2.mat');

h_wli = [0.0528 0.5358 1.0242 1.5119 2.0 2.4873 2.9739] - 0.0528; % Measurement of step height from WLI
std_wli = [0.0042  0.0064  0.0067  0.0066   0.0066  0.0061  0.0066];

figure(10); set(gcf,'Position',[400 100 480 300]);
plot(h_wli,h_wli,'--k','LineWidth',1.5);
hold on;
errorbar(h_wli,stepResults_1.measured_h,stepResults_1.sigmah,'o','Color',"#0072BD", ...
    'MarkerSize',6,'MarkerFaceColor',[0.8 0.8 0.8],'CapSize',12,'LineWidth',1.5); % Lambda = 4.84;
errorbar(h_wli,stepResults_2.measured_h,stepResults_2.sigmah,'o','Color',"#D95319", ...
    'MarkerSize',6,'MarkerFaceColor',[0.8 0.8 0.8],'CapSize',12,'LineWidth',1.5); % Lambda = 9.17;
hold off;

grid on; axis padded;
set(gca,'FontWeight','bold','FontSize',16,'XTick',round(h_wli,2),'YTick',0:0.5:3.0);
xlabel('Height reference of WLI (mm)'); ylabel('Experiment (mm)');
legend('','\Lambda = 5.04 mm','\Lambda = 9.17 mm','Location','northwest');

%%  Ploting  ------ 1. Field of individual
E_hn = sqrt(I_HN); E_hn(E_hn > 0.8) = nan;
Field1 = E_hn.*exp(1i*angle(field_onfocus_HN));
figure(91);set(gcf,'Position',[400 600 400 300]);
imagesc(X(1,:)*dx/mm,Y(:,1)*dy/mm,field2pic(Field1));axis image off;
% xlabel({'{\itx} (mm)'});ylabel({'{\ity} (mm)'}); title('Amp. for HN');drawnow;

E_ld = sqrt(I_LD); E_ld(E_ld > 0.7) = 0.7;
Field2 = E_ld.*exp(1i*angle(field_onfocus_LD));
figure(92);set(gcf,'Position',[400 100 400 300]);
imagesc(X(1,:)*dx/mm,Y(:,1)*dy/mm,field2pic(Field2));axis image off;
% xlabel({'{\itx} (mm)'});ylabel({'{\ity} (mm)'}); title('Phase for LD');drawnow;

%% Ploting  ------ 2. Beat Phase/ Height & cross profile
phiB_prop = angle(field_onfocus_HN./field_onfocus_LD.*exp(1i*(+2.102)));
phiB_prop(I_HN < 0.03 | I_LD < 0.03) = nan;

height_map = (phiB_prop)/4/pi*Lambda/mm;
% % median filter
height_map = medfilt2(height_map, [7 7]);

figure(95);set(gcf,'Position',[400 400 400 300]);
imagesc(X(1,:)*dx/mm,Y(:,1)*dy/mm,height_map,'AlphaData',~isnan(height_map));axis image off;
colormap(gca,"hsv");
hold on;     plot([-2.5 5.5],[0,0],'-.k',LineWidth = 2.0);    hold off;
%%
lh_idx = [400 1300 4000 1350];
lineh = exp(1i*phiB_prop(lh_idx(2):lh_idx(4),:)); lineh = angle(mean(lineh,1,"omitnan"));
lineh_height = lineh/4/pi*Lambda/mm;
% median filter
lineh_height = medfilt1(lineh_height, 7);

figure(97);set(gcf,'Position',[400 100 400 180]);
plot(X(1,:)*dx/mm,lineh_height,'-k','LineWidth',1.0,'MarkerSize',6); axis padded;grid on;
set(gca,'FontWeight','Bold','FontSize',16,'XLim',[-2.5 5.5]);
xlabel({'{\itx} (mm)'},'fontsize',18); ylabel({'{\ith} (mm)'},'fontsize',18); 

%% Ploting  ------ 3. Height map 3D 
ridx = [1401 501 3900 2000]; % [x1 y1 x2 y2]
sub_map = height_map;  sub_map(I_HN < 0.1 | I_LD < 0.1) = nan;
sub_map = sub_map(ridx(2):ridx(4),ridx(1):ridx(3));

% median filter
submap_filtered = medfilt2(sub_map, [11 11]);

figure(98);colormap hsv;set(gcf,'Position',[200 400 500 360]);
mesh(X(1,ridx(1):ridx(3))*dx/mm,Y(ridx(2):ridx(4),1)*dy/mm,submap_filtered,'EdgeAlpha',0.5);box on;
zlim([-3.5 5.0]);xlim([-3 5]);ylim([-2 2]);view([12 24]);
set(gca,'YDir','reverse','FontWeight','bold','FontSize',16);

%% Data From WLI
steps_wli = load('C:\Users\hshangguan\OneDrive - Delft University of Technology\BackUp\ExperimentalData\Peiyu_doc\Copper_upper_side.mat');
steps_wli = steps_wli.rawDataMatrix;

h_wli = -flipud(steps_wli)*um;
px_wli = 3.895*um;
X_wli = -(2827-1)/2:1:(2827-1)/2;

figure(97);hold on;
% plot((X_wli*px_wli/mm)-1.3,h_wli(:,7)-1.0,'.r','MarkerSize',8);grid on;
scatter((X_wli*px_wli/mm)-1.3, h_wli(:,7)-1.0, 6, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
hold off

%% Decorrelation noise
sq = 5*um;
sigma_h = Lambda/4/pi*(7/4)*(1-exp(-8*(pi*sq/Lambda).^2)).^(2/5);

fprintf('The precision in height measurement is %f (um)\n',sigma_h/um);


