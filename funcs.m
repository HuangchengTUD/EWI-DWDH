classdef funcs
    methods(Static)

        function [holo,M,N,X,Y] = holo_read(varargin)
            filename = varargin{1};
            holo = double(imread(filename));
            holo = holo(:,:,1);
            if nargin > 1
                disp('Zero padding the hologram.')
                [a,b] = size(holo);
                holo = padarray(holo,[a/2,b/2]);
            end
            [N,M] = size(holo); % N rows; M columns
            [X,Y] = meshgrid(-M/2:M/2-1,-N/2:N/2-1); 
        end

        function F = FT(f)
            F = fftshift(fft2(fftshift(f)));
        end

        function BW = threshold_FT(FT_holo, M, N)
            I = sqrt(abs(FT_holo));
            px = 30; I(N/2-px:N/2+px, M/2-px:M/2+px) = 0;
            I = I./max(I(:));

            T = graythresh(I);
            BW = imbinarize(I, T);
%             BW = imbinarize(I, 'adaptive', 'Sensitivity', 0.4);
        end

        function [plus_coor,m,n,p,q] = get_plus1(BW, M, N)
            cc = bwconncomp(BW,4);
            numPixels = cellfun(@numel,cc.PixelIdxList);
            idx_list = zeros(2,1);
            for i = 1:2
                [~,idx] = max(numPixels);
                idx_list(i) = idx;
                numPixels(idx) = 0;
            end

            for i = 1:length(cc.PixelIdxList)
                if ~(i == idx_list(1) || i == idx_list(2) )
                    BW(cc.PixelIdxList{i}) = 0;
                end
            end
            props = regionprops(BW); % 'Area' 'Centroid' 'BoundingBox'

            box_size = props(1).BoundingBox;
            plus_coor = fix(props(1).Centroid);
            dc_coor = [M/2, N/2]+1;
            p_and_q = plus_coor - dc_coor;

            figure, imshow(BW);set(gca,'FontSize',16);
            axis on;
            hold on;
            plot(plus_coor(1), plus_coor(2), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
            plot(dc_coor(1), dc_coor(2), 'r*', 'MarkerSize', 10, 'LineWidth', 2);

            line([plus_coor(1),dc_coor(1)],[plus_coor(2),dc_coor(2)],'Color','green','LineStyle','--','linewidth',2)
            for i = 1:numel(props)
                rectangle('Position', props(i).BoundingBox, ...
                    'LineWidth', 3, 'EdgeColor', 'r', 'LineStyle', '--');
            end

            % values for p,q,m and n (min and kemper paper)
            m = box_size(3);
            n = box_size(4);
            p = p_and_q(1);
            q = p_and_q(2);
            disp('p: '+string(p)+' q: '+string(q));
            disp('m: '+string(m)+' n: '+string(n));
        end

%         function Img_field = filter_center_plus1(FT_holo,plus_coor,M,N,m,n,p,q,X,Y,dx,dy)
%             Filter = zeros(N,M);
%             a = min([q-20 n]); b = min([p-20 m]);
%             Filter((plus_coor(2)-a):(plus_coor(2)+a), (plus_coor(1)-b):(plus_coor(1)+ b)) = 1;
%             FT_holo_filter = FT_holo.*Filter;
%             holo_filtered = fftshift(ifft2(fftshift(FT_holo_filter)));
% 
%             Reference = exp(1i*2*pi*(p/(M*dx)*X*dx + q/(N*dy)*Y*dy));
%             Img_field = holo_filtered.*Reference;
%         end

        function [field_obj,slope_phase,field_plus1,C0] = get_offaxis_term(holo,M,N,X,Y,dx,dy,lambda)
            FT_holo = funcs.FT(holo);
            BW = funcs.threshold_FT(FT_holo, M, N);
            [plus_coor,m,n,p,q] = funcs.get_plus1(BW, M, N);
            C0 = mean([(M*dx).^2/lambda/m, (N*dy).^2/lambda/n]);
            disp('Initial curvature: '+string(C0));

            % create filter window
            Filter = zeros(N,M);
            a = fix(min([abs(q)-20, 1.3*n])); b = fix(min([abs(p)-20, 1.3*m]));
            x1 = plus_coor(1)-b; x2 = plus_coor(1)+ b; y1 = plus_coor(2)-a; y2 = plus_coor(2)+a;
            Filter( y1:y2, x1:x2 ) = 1;
                figure;colormap jet;
                imagesc(log10(abs(FT_holo)));axis image;colorbar
                hold on;
                rectangle('Position',[x1,y1,2*b,2*a], 'LineWidth',3, 'EdgeColor','k', 'LineStyle','-.');
                hold off;

            field_plus1 = fftshift(ifft2(fftshift(FT_holo.*Filter)));
            slope_phase = exp(-1i*2*pi*(p/(M*dx)*X*dx + q/(N*dy)*Y*dy));
            field_obj = field_plus1.*slope_phase;
        end

        function [g,h] = get_g_and_h(field_obj,M,N)
            bw = imbinarize(real(field_obj));

            cc = bwconncomp(bw,4);
            numPixels = cellfun(@numel, cc.PixelIdxList);
            num_of_term = 3;
            idx_list = zeros(num_of_term,1);
            for i = 1:num_of_term
                [~,idx] = max(numPixels);
                idx_list(i) = idx;
                numPixels(idx) = 0;
            end

            for i = 1:length(cc.PixelIdxList)
                if ~(i == idx_list(1) || i == idx_list(2) || i == idx_list(3))
                    bw(cc.PixelIdxList{i}) = 0;
                end
            end

            props = regionprops(bw,'Centroid','BoundingBox','Eccentricity');
            [~,best_term_idx] = min([props.Eccentricity]);
            term_center = props(best_term_idx).Centroid;
            dc_center = [M/2, N/2]+1;
            g_and_h = (term_center - dc_center);

            figure, imshow(bw); set(gca,'FontSize',16);
            axis on;
            hold on;
            plot(dc_center(1), dc_center(2), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
            plot(term_center(1), term_center(2), 'b*', 'MarkerSize', 10, 'LineWidth', 2);
            line([dc_center(1),term_center(1)],[dc_center(2),term_center(2)],'Color','yellow','LineStyle','--')
            rectangle('Position', props(best_term_idx).BoundingBox, 'Linewidth', 3, 'EdgeColor', 'r', 'LineStyle', '--');

            % values for g and h
            g = g_and_h(1);
            h = g_and_h(2);
            disp('g: '+string(g)+' h: '+string(h));
        end

        function sign = determine_SP_sign(func, C0, field_obj)
            phi_spherical = func(C0);

            phase_mask_p = exp(1i * phi_spherical);
            J1 = angle(field_obj .* phase_mask_p);
            phase_mask_n = exp(-1i * phi_spherical);
            J2 = angle(field_obj .* phase_mask_n);

            [~, idx] = min([std(J1(:)), std(J2(:))]);
            if idx == 1
                sign = 1;
                disp('The sign of spherical phase is +');
            else
                sign = 0;
                disp('The sign of spherical phase is -');
            end
        end

        function [J] = std_cost_func_sphetical(fun,curv,field_obj)
            phi_spherical = fun(curv);
            corrected_image = field_obj .* phi_spherical;
            phase = angle(corrected_image);
            J = std(phase,[],'all');
            if curv == 0
                J = 0.5;
            end
        end

        function [J] = bin_CF_noTele_BAR_1d(fun,curv,field_obj,M,N,sign)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Title: CF_noTele_BAR0531                                                     %
            %                                                                              %                                                                        
            % Authors: Raul Castaneda and Ana Doblas                                       %
            % Department of Electrical and Computer Engineering, The University of Memphis,% 
            % Memphis, TN 38152, USA.                                                      %   
            %                                                                              %
            % Email: adoblas@memphis.edu                                                   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            phi_spherical = fun(curv);
            if sign
                phase_mask = exp((1j)*phi_spherical);
            else
                phase_mask = exp((-1j)*phi_spherical);
            end
            corrected_image = field_obj .* phase_mask;
            phase = angle(corrected_image);
            phase = phase + pi; 
            ib = imbinarize(phase, 0.5);
            J = M*N - sum(ib(:));
          end 

        function [corrected_image,spherical_phase] = automatic_method(field_obj,M,N,X,Y,lambda,dx,dy,C0,algo,cost)
            
            % get the center of the remaining spherical phase 'g' and 'h'
            [g,h] = funcs.get_g_and_h(field_obj,M,N);

            phi_spherical_creater = @(C) exp(1i*(pi/lambda/C)*(( (X-g)*dx ).^2 + ( (Y-h)*dy ).^2));

            % determin the sign of the spherical phase + or -
            sign = funcs.determine_SP_sign(phi_spherical_creater,C0,field_obj);

            %
            if cost == 1
                minfunc = @(C) funcs.std_cost_func_sphetical(phi_spherical_creater,C,field_obj);
            elseif cost == 2
                minfunc = @(C) funcs.bin_CF_noTele_BAR_1d(phi_spherical_creater,C,field_obj,M,N,sign);
            end

            cost_fun = {'STD cost function','BIN cost function'};
            alg_array = ['GA','GA+PS','FMC','FMU','FSO','SA','PTS',];
            alg = alg_array(algo);

            cur0 = C0; lb = min([0.5*cur0, 1.5*cur0]); ub = max([0.5*cur0, 1.5*cur0]);
            disp(['initial curvature: ', num2str(cur0)]);

            if alg == 'GA'
                disp(['Running the GA algorithm with the ', cost_fun{1}])
                options = optimoptions('ga','Display','final', ...
                    'InitialPopulationMatrix',cur0, ...
                    'InitialPopulationRange',[lb;ub], ...
                    'PopulationSize',15, ... % smaller one speeds up convergence
                    'SelectionFcn','selectionremainder', ...
                    'PlotFcn',@gaplotbestf);
                [out,fval] = ga(minfunc,1,[],[],[],[],[],[],[],options);

            elseif alg == 'GA+PS'
                disp(['Running the hybrid (GA + PS) algorithm with the ', cost_fun{1}])
                hybridopts = optimoptions('patternsearch','Display','final');
                options = optimoptions('ga', 'Display', 'final', ...
                    'HybridFcn', {@patternsearch,hybridopts},...
                    'InitialPopulationRange', [lb;ub],...
                    'PopulationSize',15, ... % smaller one speeds up convergence but not needed
                    'SelectionFcn','selectionremainder', ...
                    'PlotFcn',{@gaplotbestf, @gaplotstopping});
                [out,fval] = ga(minfunc,1,[],[],[],[],[],[],[],options);

            elseif alg == 'FMC'
                disp(['Running the fmincon algorithm with the ', cost_fun{1}])
                options = optimoptions('fmincon','Display','final');
                [out,fval] = fmincon(minfunc,cur0,[],[],[],[],lb,ub,[],options);

            elseif alg == 'FMU'
                disp(['Running the fminunc algorithm with the ', cost_fun{1}])
                options = optimoptions('fminunc','Display','final');
                [out,fval] = fminunc(minfunc,cur0,options);
            else
                disp('No proper optimization method found')
            end
            
            disp(['Optimized curvature: ', num2str(out)]);
%             disp(['value of the costfnc: ', num2str(fval)]);

            spherical_phase = phi_spherical_creater(out);
            corrected_image = field_obj .* spherical_phase;
            figure();colormap gray;
            imagesc(angle(corrected_image));axis image; colorbar; title('Corrected image');
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%     slopephase correction     %%%%%%%%%%%%%%%  
        function J = std_costfunc_slope(func,a,uncorrected_img,M,N)
            phase_mask = func(a);
            corrected_image = uncorrected_img .* phase_mask;
            phase = angle(corrected_image);
            J = std(phase,[],'all');
        end

        function [corrected_image,phi_slope] = slope_corrector(uncorrected_img,M,N,X,Y,dx,dy)
            phi_slope_creater = @(a) exp(-1i*2*pi*(a(1)/(M*dx)*X*dx + a(2)/(N*dy)*Y*dy));
            
            minfunc = @(a) funcs.std_costfunc_slope(phi_slope_creater,a,uncorrected_img,M,N);
            initial_guess = [0, 0];
            lb = [-3,-3];
            ub = [3, 3];
            options = optimoptions('ga','Display','final', ...
                'InitialPopulationMatrix',initial_guess, ...
                'PopulationSize',15,...
                'SelectionFcn','selectionremainder', ...
                'PlotFcn',@gaplotbestf);
            a_opt = ga(minfunc,2,[],[],[],[],lb,ub,[],options);

            phi_slope = phi_slope_creater(a_opt);
            corrected_image = uncorrected_img .* phi_slope;

            disp(['Optimized slope frequencies: ', num2str(a_opt)]);
            figure,colormap gray;
            imagesc(angle(corrected_image));axis image;colorbar;title('After correction')

        end
%%%%%%%%%%%%%%%     slopephase correction     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [field_plus1,dig_ref,field_centered] = get_DW_offaxis_term(holo,M,N,X,Y,dx,dy,lambda)
            FT_holo = funcs.FT(holo);
            BW = funcs.threshold_FT(FT_holo, M, N);

            cc = bwconncomp(BW,4); % 4-connected area
            numPixels = cellfun(@numel,cc.PixelIdxList);
            idx_list = zeros(4,1);
            for i = 1:length(idx_list)
                [~,idx] = max(numPixels);
                idx_list(i) = idx;
                numPixels(idx) = 0;
            end

            for i = 1:length(cc.PixelIdxList)
                if ~ismember(i, idx_list)
                    BW(cc.PixelIdxList{i}) = 0;
                end
            end

            props = regionprops(BW); % 'Area' 'Centroid' 'BoundingBox'
            dc_coor = [M/2, N/2]+1;

            fig1 = figure(10); imshow(BW);set(gca,'FontSize',16);
            hold on;
            plot(dc_coor(1), dc_coor(2), 'r*', 'MarkerSize', 10, 'LineWidth', 2);

            field_centered = cell(2,1); dig_ref = cell(2,1);   field_plus1 = cell(2,1);
            for kk = 1:2
                % get +1 coordinate
                box_size = props(kk).BoundingBox;
                plus_coor = fix(props(kk).Centroid);
                p_and_q = plus_coor - dc_coor;

                figure(fig1);
                plot(plus_coor(1), plus_coor(2), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
                line([plus_coor(1),dc_coor(1)],[plus_coor(2),dc_coor(2)],'Color',[2-kk kk-1 0],'LineStyle','--','linewidth',2);
                rectangle('Position', props(kk).BoundingBox, ...
                    'LineWidth', 2, 'EdgeColor', [2-kk kk-1 0], 'LineStyle', '--');

                % values for p, q and m, n (中心点和宽�?�?
                m = box_size(3); n = box_size(4);
                p = p_and_q(1);  q = p_and_q(2);
                text(box_size(1),box_size(2)-100,['m: '+string(m)+'  n: '+string(n) ...
                    'p: '+string(p)+'  q: '+string(q)],'fontsize',16,'Color','yellow');

                % Define the filter windows
                Filter = zeros(N,M);
                a = fix(min([abs(q)-40, 1.2*n])); b = fix(min([abs(p)-40, 1.2*m]));
                x1 = plus_coor(1)-b; x2 = plus_coor(1)+ b; y1 = plus_coor(2)-a; y2 = plus_coor(2)+a;
                Filter( y1:y2, x1:x2 ) = 1;

                % 1. get the +1 term
                field_plus1{kk} = fftshift(ifft2(fftshift(FT_holo.*Filter)));
                % 2. get slope phase from 'p & q'
                slope_phase = exp(-1i*2*pi*(p/(M*dx)*X*dx + q/(N*dy)*Y*dy));
                field_shift = field_plus1{kk}.*slope_phase;
                % 3. initial curvature of the parabolic phase
                c0 = mean([(M*dx).^2/lambda/m, (N*dy).^2/lambda/n]);
                [g,h] = funcs.get_g_and_h(field_shift,M,N);
                % 4. determine the sign 
                phase_mask_p = exp(1i*(pi/lambda/c0)*(((X-g)*dx).^2 + ((Y-h)*dy).^2));
                J1 = field_shift .* phase_mask_p;
                [dpx, dpy] = gradient(angle(J1));
                jump_count1 = sum((abs(dpx(:)) > 0.8*pi) | (abs(dpy(:)) > 0.8*pi));

                phase_mask_n = 1./phase_mask_p;
                J2 = field_shift .* phase_mask_n;
                [dpx, dpy] = gradient(angle(J2));
                jump_count2 = sum((abs(dpx(:)) > 0.8*pi) | (abs(dpy(:)) > 0.8*pi));

                if jump_count1 < jump_count2
                    field_centered{kk} = J1;
                    dig_ref{kk} = slope_phase.*phase_mask_p;
                else
                    field_centered{kk} = J2;
                    dig_ref{kk} = slope_phase.*phase_mask_n;
                end

            end
        end

        function [cleaned] = cpx_filloutliers(cpxinput)
            amp = abs(cpxinput); phi = angle(cpxinput);
            amp_cleaned = filloutliers(amp,'nearest', 'quartiles', 'ThresholdFactor', 2.0);
            cleaned = amp_cleaned.*exp(1i*phi);

            num = nnz(amp ~= amp_cleaned);
            fprintf('The number of outliers: %f \n', num);
        end


        %{
            ##################################################################
            ################Functions for manual determination################
            ##################################################################
        %}

        function varargout = manual_get_plus1(holo,M,N,X,Y,dx,dy,holo_type)
            if nargin < 8
                holo_type = 1;
            end 
            
            % Initialize variables
            SPE = funcs.FT(holo)/sqrt(M*N);   PSD_log = log10(abs(SPE));
            field_plus1 = cell(holo_type,1);
            filter = cell(holo_type,1);
            dig_ref = cell(holo_type,1);
            
            % Display log(PSD) and select ROI
            fig1 = figure(1); fig1.WindowState = 'maximized'; imagesc(PSD_log); colorbar; title('Select the ROI');daspect([1 1 1])

            for k = 1: holo_type
                % select ROI of +1order interactively
                figure(fig1);
                r = drawrectangle('Label',['+1 order for lambda',num2str(k)],'Color',[2-k k-1 0]);
                ROI_idx = r.Position;

                % Region index & creat filter mask
                y1 = fix(ROI_idx(2)); y2 = y1 + fix(ROI_idx(4));
                x1 = fix(ROI_idx(1)); x2 = x1 + fix(ROI_idx(3));

                % Creat filter window -- 'cosineTaperWindow'
                filter{k} = zeros(N,M);

                row = tukeywin(fix(ROI_idx(4)),0.4);
                col = tukeywin(fix(ROI_idx(3)),0.4);
                cosineWin = row * col';
                filter{k}(y1:y2-1 , x1:x2-1) = cosineWin; % sharpEdgeWindow: 'filter{k}(y1:y2,x1:x2) = 1;'

                % Filter our the field(s) of +1
                field_plus1{k} = fftshift(ifft2(fftshift(filter{k}.*SPE))).*sqrt(M*N);

                if nargout >= 3
                    % 1. extract ROI 
                    roi_mass = abs(SPE(y1:y2-1 , x1:x2-1)).^2;
                    % 2. find the maximun
                    [~, max_idx] = max(roi_mass(:));
                    [max_row, max_col] = ind2sub(size(roi_mass), max_idx);
                    % 3. define calculation window
                    window_size = 20; % half size
                    row_min = max(1, max_row - window_size);
                    row_max = min(size(roi_mass,1), max_row + window_size);
                    col_min = max(1, max_col - window_size);
                    col_max = min(size(roi_mass,2), max_col + window_size);
                    % 4. extract subregion's coordinate and weight
                    fx_roi = X(y1+row_min-1:y1+row_max-1, x1+col_min-1:x1+col_max-1)/M/dx;
                    fy_roi = Y(y1+row_min-1:y1+row_max-1, x1+col_min-1:x1+col_max-1)/N/dy;
                    roi_mass_sub = roi_mass(row_min:row_max, col_min:col_max);
                    % 5. calculate centroid
                    total_mass = sum(roi_mass_sub(:));
                    p = sum(sum(fx_roi.*roi_mass_sub)) / total_mass;
                    q = sum(sum(fy_roi.*roi_mass_sub)) / total_mass; % get by pick [p,q] = ginput(1);   pause(.1)

                    dig_ref{k} = exp(-1i*2*pi*(p*X*dx + q*Y*dy));
                end
            end
            % Set varargout based on the number of requested output arguments
            varargout{1} = field_plus1;
            if nargout >= 2
                varargout{2} = filter;
            end
            if nargout >= 3
                varargout{3} = dig_ref;
            end
        end


        %{
            ##################################################################
            ################        Auxiliary Functions       ################
            ##################################################################
        %}
        function I1 = nmlz(I0,varargin)
            mn = 0; mx = 1; % default min & max
            Imax = max(I0(:)); Imin = min(I0(:));

            for idx_var = 1:2:length(varargin)
                key = varargin{idx_var};
                switch key
                    case 'min'
                        mn = varargin{idx_var+1};
                    case 'max'
                        mx = varargin{idx_var+1};
                    otherwise
                        error(['undefined key : ',key])
                end
            end

            I1 = (I0-Imin)./(Imax-Imin).*(mx-mn) + mn;
        end

        function Uout = rangePropagation(Uin,M,N,X,Y,dx,dy,lambda,dis1,dis2,S,x1,x2,y1,y2)
            kz = 2*pi*sqrt((1/lambda).^2-(X/M/dx).^2-(Y/N/dy).^2);
            ForwardPropagate = @(field,dis,kz) ifft2(fftshift(exp(1i*kz*dis)).*fft2(field));

            apodized_mask = funcs.CosineTaperWindow(M,N,0.1);
            dis = linspace(dis1,dis2,S);
            
            for kk = 1:length(dis)

                Uout = ForwardPropagate(Uin.*apodized_mask,dis(kk),kz);
                subarea = Uout(y1:y2, x1:x2);

                figure(30);colormap gray
                subplot 121;imagesc(abs(subarea));axis image;colorbar;title(['Amplitude; dis = ',num2str(dis(kk)*1000)]);

                subplot 122;imagesc(angle(subarea));axis image;colorbar;title('Phase');drawnow;
%                 pause(0.2)
            end
        end


        function S = LaguerreGaussian(X,Y,dx,dy,p,l,w0s)
            syms  variable
            Laguerre = 0;

            for k = 0:p
                C = factorial(p+abs(l))/factorial(abs(l)+k)/factorial(k)/factorial(p-k);
                Laguerre = Laguerre + C*(-variable).^k;
            end
            Laguerre = eval(['@(variable)',vectorize(Laguerre)]);


            S= (sqrt(2).*(X*dx+sign(l)*1i*Y*dy)./w0s).^(abs(l)).*Laguerre(2*((X*dx).^2+(Y*dy).^2)/w0s.^2).*exp(-(sqrt((X*dx).^2+(Y*dy).^2)/w0s).^2);
        end

        
        %%%%%%%%%%%%%%%%%%%%%   Propagation   %%%%%%%%%%%%%%%%%%%%%
        function [U1,Dx,Dy] = Fresnelprop(U0,lambda,dis,M,N,X,Y,dx,dy)
            Dx = lambda*dis/M/dx; Dy = lambda*dis/N/dy;

            A0 = exp(1i*pi/lambda/dis*((X*dx).^2 + (Y*dy).^2));
            A1 = exp(1i*pi/lambda/dis*((X*Dx).^2 + (Y*Dy).^2));

            U1 = exp(1i*2*pi*dis/lambda)/(1i*lambda*dis).*A1.*funcs.FT(U0.*A0)*dx*dy;
        end

        function padedfield = doublepad (infield)
            [a,b] = size(infield);
            padedfield = padarray(infield,[a/2,b/2]);
        end

        function [outfield] = BLASprop(infield,M,N,X,Y,dx,dy,lambda,dis,methd)
            if nargin ~= 10
                methd = 1;
            end

            if methd == 2
                infield = funcs.doublepad(infield);
                M = 2*M;  N = 2*N;
                [X,Y] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
            end

            X_c = M*dx./sqrt(1+(2*dis./(M*dx)).^2)/lambda;  X_c = fix(X_c);
            Y_c = N*dy./sqrt(1+(2*dis./(N*dy)).^2)/lambda;  Y_c = fix(Y_c);
            kz = 2*pi*sqrt((1/lambda).^2-(X/M/dx).^2-(Y/N/dy).^2);
            kz( abs(X) > X_c | abs(Y) > Y_c ) = 0;

            mask = ones(N,M);
            mask( abs(X) > X_c | abs(Y) > Y_c ) = 0;

            outfield = ifft2(fftshift(mask.*exp(1i*kz*dis)).*fft2(infield));

            switch methd
                case 1
                    % No padding
                case 2
                    outfield = outfield(N/4+1 : N/4*3, M/4+1 : M/4*3);
            end

        end









        function meanvalue = drawandmean(input)
            if ~isreal(input)
                input = angle(input);
            end

            fig1 = figure('WindowState','maximized');colormap jet;
            imagesc(input);axis image;colorbar

            figure(fig1);
            r = drawrectangle('Label','Region of Interest','Color','white');
            
            % Region index 
            ridx = r.Position;
            y1 = fix(ridx(2)); y2 = y1 + fix(ridx(4)); x1 = fix(ridx(1)); x2 = x1 + fix(ridx(3));

            values = exp(1i*input(y1:y2-1,x1:x2-1));
            meanvalue = angle(mean(values,'all','omitnan'));
            fprintf('The mean value of this region is: (%f)\n', meanvalue);
        end

        function [meanvalue,idx,stdvalue] = drawSeveral(input)
            if ~isreal(input)
                input = angle(input);
            end
            figure('WindowState','maximized');colormap jet;
            imagesc(input);axis image;colorbar

            meanvalue = []; idx = {}; stdvalue = [];

            while true
                h = imrect;
                pos = wait(h);

                if isempty(pos)
                    break
                end

                ridx = fix(pos);
                y1 = ridx(2); y2 = ridx(2) + ridx(4); x1 = ridx(1); x2 = ridx(1) + ridx(3);

                rectangle('Position', ridx, 'EdgeColor', 'k', 'LineWidth', 1.5);
                text(x1-40,y1-60,['ROI:',num2str(length(meanvalue)+1)],'Color','k','FontSize',16,'FontWeight','bold'); % showing the region

                values = exp(1i*input(y1:y2-1,x1:x2-1));
                meanvalue(end+1) = angle(mean(values,'all','omitnan'));
                stdvalue(end+1) = std(values(:),'omitnan');
                idx{end+1} = ridx;
                fprintf('The mean value of this region is: (%f)\n', meanvalue(end));
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   background fitting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%$$$$$$$$$$$$$$$ 1.fitting the input $$$$$$$$$$$$$$$$$$
        function [coeff,phase_rgst] = aberfitting(phasemap,order)
            % phasemap: real/complex phase term
            % X,Y: normalized coordinate [-1 1]
            % Output: phase_rgst -- real phase value
            
            if ~isreal(phasemap)
                phasemap = angle(phasemap);
            end
            [N,M] = size(phasemap);  x_ = linspace(-1,1,M); y_ = linspace(-1,1,N);
            [X, Y] = meshgrid(x_,y_);

            % creating vectors
            validIdx = find(~isnan(phasemap));
            Fxy = phasemap(validIdx);
            x_idx = X(validIdx); y_idx = Y(validIdx);

            % creating A matrix
            nt = (order+1)*(order+2)/2; % number of polynomial terms
            A = zeros([length(Fxy),nt]);
            h = 1;
            for k = 0 : order
                for beta = 0 : k
                    alpha = k - beta;
                    A(:,h) = (x_idx.^alpha).*(y_idx).^beta;
                    h = h + 1;
                end
            end

            % Fitting polynomial coefficients by least-square method
            AT = A.';
            coeff = (AT*A)\(AT*Fxy);

            % validate the coefficients and calculate errors
            phase_fit = 0;
            h = 1;
            for k = 0 : order
                for beta = 0:k
                    alpha = k - beta;
                    phase_fit = phase_fit + coeff(h) .* X.^alpha .* Y.^beta;
                    h = h + 1;
                end
            end

            phase_rgst = angle(exp(1i*phasemap)./exp(1i*phase_fit));
            error = mean(abs(phase_rgst),'all');

%             fprintf('Average error: (%.4e)\n', error);
        end

        %%%%%%%%%%$$$$$$$$$$$$$$$ 2.fitting by drawing $$$$$$$$$$$$$$$$$$
        function [coeff,z_fit] = aberfittingbydraw(phasemap,order,X,Y)
            if ~isreal(phasemap)
                phasemap = angle(phasemap);
            end
            fig1 = figure('WindowState','maximized');colormap jet;
            imagesc(phasemap);axis image;colorbar; title('Draw ROIs and double-click to confirm; Press Esc to exit');drawnow;

            % select ROI and get index
            rectangles = {};
            while true
                h = imrect;
                pos = wait(h);

                if isempty(pos)
                    break; %
                end

                rectangles{end+1} = pos;
                rectangle('Position', pos, 'EdgeColor', 'w', 'LineWidth', 1.5);
            end

            % get indices within each 'rect' and concatenate
            indices = cell(size(rectangles));

            for k = 1:length(rectangles)
                pos = rectangles{k};
                x1 = round(pos(1));
                y1 = round(pos(2));
                x2 = round(pos(1) + pos(3));
                y2 = round(pos(2) + pos(4));

                [subX, subY] = meshgrid(x1:x2, y1:y2);
                indices{k} = sub2ind(size(phasemap), subY(:), subX(:));
            end

            totalidx = vertcat(indices{:}); % vertical concatenation
            idx = false(size(phasemap)); idx(totalidx) = true;

            %
            Fxy = phasemap(idx);
            x_idx = X(idx); y_idx = Y(idx);

            % creating A matrix
            nt = (order+1)*(order+2)/2; % number of polynomial terms
            A = zeros([length(Fxy),nt]);
            h = 1;
            for k = 0 : order
                for beta = 0 : k
                    alpha = k - beta;
                    A(:,h) = (x_idx.^alpha).*(y_idx).^beta;
                    h = h + 1;
                end
            end

            % Fitting polynomial coefficients by least-square method
            AT = A.';
            coeff = (AT*A)\(AT*Fxy);

            % validate the coefficients and calculate errors
            z_fit = 0;
            h = 1;
            for k = 0 : order
                for beta = 0:k
                    alpha = k - beta;
                    z_fit = z_fit + coeff(h) .* X.^alpha .* Y.^beta;
                    h = h + 1;
                end
            end

            error = z_fit - phasemap;
            fprintf('The mean and std of the error: (%.4e, %.4e)\n', std(error(:),'omitnan'), mean(error(:),'omitnan'));

            figure();colormap hsv;
            subplot 131;imagesc(z_fit);axis image;colorbar;title('Fitted');drawnow;
            subplot 132;imagesc(angle(exp(1i*phasemap)./exp(1i*z_fit)));axis image;colorbar;drawnow;
            subplot 133;bar(0:nt-1,coeff); axis square; grid on;
            xlabel('Standard poly coefficients');drawnow;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   background fitting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Visibility,Mask,AC,DC] = fringeVis(hologram,M,N,X,Y,dx,dy)
            [Field,Mask] = funcs.manual_get_plus1(hologram,M,N,X,Y,dx,dy,2);
            AC = abs(Field{1}); DC = abs(Field{2});
            Visibility = 2*AC./DC; 

            figure();subplot 121;imagesc(Visibility);colorbar;axis image;title('Fringe visibility');
            subplot 122;histogram(Visibility,'EdgeColor','none');
        end

        function filterWindow = CosineTaperWindow(M,N,alpha)
            % alpha - percentage of taper width
            if nargin < 3
                alpha = 0.5;
            end

            row = tukeywin(N,alpha);
            col = tukeywin(M,alpha);

            filterWindow = row * col';
        end









    end
end

