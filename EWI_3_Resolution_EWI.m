clear
clc
close all

mm = 1e-3; um = 1e-6; nm = 1e-9;
set(0,'defaultAxesFontName', 'times new roman','defaultAxesFontSize',16);
set(0,'defaultTextFontName', 'times new roman','defaultTextFontSize',16);
%%  read hologram & get the off-axis terms
% Variables setting
lambda1 = 633e-9;
lambda2 = 633.09e-9;
dx = 3.45e-6; dy = 3.45e-6; % pixel size
datapath = 'Experiments\ResolutionTarget\';

raw_img = 'dw_ew2_usaf.tif'; disp(['The raw hologram is : ', raw_img]);
[holo,M,N,X,Y] = funcs.holo_read([datapath,raw_img]);
[field_plus1,dig_ref,field_centered] = funcs.get_DW_offaxis_term(holo,M,N,X,Y,dx,dy,lambda1);

figure;
subplot 121;imagesc(angle(field_centered{2}));axis image;title('Off-axis term of HN');drawnow;
subplot 122;imagesc(angle(field_centered{1}));axis image;title('Off-axis term of LD');drawnow;
return
%% phase aberration correction --- HN part
kz1 = 2*pi*sqrt((1/lambda1).^2-(X/M/dx).^2-(Y/N/dy).^2);
kz2 = 2*pi*sqrt((1/lambda2).^2-(X/M/dx).^2-(Y/N/dy).^2);
ForwardPropagate = @(field,dis,kz) ifft2(fftshift(exp(1i*kz*dis)).*fft2(field));

%%%%%%%%%%%%%%%% Background aberration fitting (1. He-Ne) %%%%%%%%%%%%%%%%%
% 1. Back propagate
field_onfocus_HN = ForwardPropagate(field_centered{2},-139.00*mm,kz1);
% 2. poly fitting (iterative)
input_0 = field_onfocus_HN;  coeff_k = zeros(15,1); totalFit = zeros(N,M);

iterNum = 3;
for kk = 1:iterNum
    [coeff,z_fit] = funcs.aberfittingbydraw(input_0,4,X*2/M,Y*2/N);

    coeff_k = coeff_k + coeff;
    totalFit = totalFit + z_fit;
    input_0 = field_onfocus_HN./exp(1i*totalFit);
end

figure;set(gcf,'Position',[600 300 465 270]);
bar(0:14,-coeff_k); axis auto; grid on;
xlabel('Term Index'); ylabel('\it p_{\alpha,\beta}'); title('Polynomials coefficients');
drawnow;

% 3. The total digital Reference beam
totalRefPhase_hn = dig_ref{2}./exp(1i*totalFit);
%%
%%%%%%%%%%%%%%%% Background aberration fitting (2. LD) %%%%%%%%%%%%%%%%%%%%
% 1. Back propagate
field_onfocus_LD = ForwardPropagate(conj(field_centered{1}),-0.00*mm,kz2);
% 2. poly fitting (iterative)
input_0 = field_onfocus_LD; coeff_k = zeros(15,1); totalFit = zeros(N,M);

iterNum = 3;
for kk = 1:iterNum
    [coeff,z_fit] = funcs.aberfittingbydraw(input_0,4,X*2/M,Y*2/N);

    coeff_k = coeff_k + coeff;
    totalFit = totalFit + z_fit;
    input_0 = field_onfocus_LD./exp(1i*totalFit);
end

figure;set(gcf,'Position',[600 300 465 270]);
bar(0:14,-coeff_k); axis auto; grid on;
xlabel('Term Index'); ylabel('\it p_{\alpha,\beta}'); title('Polynomials coefficients');
drawnow;

% 3. The total digital Reference beam
totalRefPhase_ld = dig_ref{1}.*exp(1i*totalFit);

%% Re-process the 'field_centered'
[field_plus1,~,~] = funcs.manual_get_plus1(holo,M,N,X,Y,dx,dy,2);

field_onfocus_HN = ForwardPropagate(field_plus1{2}.*totalRefPhase_hn,-138.30*mm,kz1);
field_onfocus_LD = ForwardPropagate(conj(field_plus1{1}.*totalRefPhase_ld),-138.92*mm,kz2);
% 1. field registration
field_onfocus_LD = circshift(field_onfocus_LD,[-17,10]);

% 2. determing magnification
subarea = abs(field_onfocus_LD(1780:1807,1744:2600)).^2;

yline = mean(subarea,1);    yline = yline./max(yline);
figure(13);set(gcf,'Position',[1000 400 900 200]);
plot(yline,'LineWidth',2.0);axis tight;grid on;

bbb = imbinarize(yline,0.3); 
hold on;plot(bbb,'-r','Linewidth',1.5);hold off;
set(gca,'FontSize',18);
xlabel('Pixels');ylabel('Intensity');legend('Data','Binarized','Location','bestoutside');

locs = find(diff(bbb) ~= 0);
widths = locs(2:2:end) - locs(1:2:end-1);
widths = reshape(widths,3,[]);
my_lw = mean(widths)*dx/um;

nomial_lw = [17.54 19.69 22.10 24.80 27.84];
ma_ = my_lw./nomial_lw;
ma = round(mean(ma_),3);
fprintf('The magnification is around: %.3f.\n',ma);

%% 3. The beat amplitude
U_EWI = field_onfocus_HN.*conj(field_onfocus_LD);
V = funcs.nmlz(abs(U_EWI));
V_new = min(max(V, 0), 0.2); % restrain V in [0 0.2]
U_EWI = V_new.*exp(1i*angle(U_EWI));

% 4. plotting
figure;colormap gray;set(gcf,'Position',[200 200 420 300]);
imagesc(X(1,:)*dx./ma/mm,Y(:,1)*dy./ma/mm,abs(U_EWI));axis image off;clim([0 0.2]);
xlim([0 0.25]);ylim([-0.25 0]);

line1 = funcs.nmlz(U_EWI(1300,2236:2266)); % group 7-1 = 3.91 um
figure(2),set(gcf,'Position',[200 200 400 150]);
plot((1:size(line1,2))*dx./ma/um,line1,'-r','LineWidth',2);axis tight;grid on;

line2 = funcs.nmlz(U_EWI(1310,2050:2110)); % group 6-1 = 7.81 um
figure(3),set(gcf,'Position',[700 200 400 150]);
plot((1:size(line2,2))*dx./ma/um,line2,'-r','LineWidth',2);axis tight;grid on;

return
%% searching distance
field_0 = conj(field_plus1{1}.*totalRefPhase_ld);
U_temp = funcs.rangePropagation(field_0,M,N,X,Y,dx,dy,lambda2,-138.70*mm,-139.20*mm,26, ...
    2000,2350,1250,1550);
%     dis = 150.9*mm;

%% Image registration
clear corrI
U1 = abs(field_onfocus_HN(1260:1380,2140:2300));
U2 = abs(field_onfocus_LD(1260:1380,2140:2300));

yshift = -5:1:5; xshift = -5:1:5; 

for ii = 1:length(yshift)
    for jj = 1:length(xshift)
        I2_shift = circshift(U2,[yshift(ii) xshift(jj)]);
        corrI(ii,jj) = corr2(U1,I2_shift);

        imagesc(U1 + I2_shift);axis image; title([num2str(ii),'; ',num2str(jj)]);drawnow
%         pause(0.2)
    end
end

%%
dis = 138.06*mm;
zs = dis/ma; zd = zs/(ma-1);
z_obj = abs(dis)/(ma-1);

fprintf('The zs, zd, zo are: %.3f; %.3f; %.3f;\n',zs/mm,zd/mm,z_obj/mm);

return
phi_s_obj = exp(1i*(pi/lambda1/z_obj)*( (X*dx).^2 + (Y*dy).^2));

aaa = funcs.Fresnelprop(corrected_image_v2.*phi_s_obj,lambda1,-zs,M,N,X,Y,dx,dy);
figure(2),colormap jet;
imagesc(abs(aaa)),axis image;colorbar




