% ========================================================================
% Introduction
% ========================================================================
% This code provides a simple demonstration of the extended ptychographical
% iterative engine (ePIE) using simulation data. This code is mainly
% based on the following papers:
%   - A. M. Maiden and J. M. Rodenburg, "An improved ptychographical phase
%     retrieval algorithm for diffractive imaging," Ultramicroscopy 109,
%     1256-1262 (2009).
%   - F. Zhang, I. Peterson, J. Vila-Comamala, A. Diaz, F. Berenguer, 
%     R. Bean, B. Chen, A. Menzel, I. K. Robinson, and J. M. Rodenburg, 
%     "Translation position determination in ptychographic coherent 
%     diffraction imaging," Optics Express 21, 13592-13606 (2013).
%
% Author: Yunhui Gao (gyh21@mails.tsinghua.edu.cn)
% =========================================================================
%%
% =========================================================================
% Data generation
% =========================================================================
clear;clc
close all

% load functions
addpath(genpath('./utils'))

% simulation settings
N1 = 512;   % image dimension (height)
N2 = 512;   % image dimension (width)

% physical parameters
params.pxsize = 2.740e-3;           % pixel size (mm)
params.wavlen = 0.532e-3;           % wavelength (mm)
params.dist_1 = 2;                  % object-to-diffuser distance (mm)
params.dist_2 = 10;                 % diffuser-to-sensor distance (mm)

% objecct settings
feature_size = 8;
obj_amp = rand(round(N1/feature_size),round(N2/feature_size));
obj_amp(obj_amp < 0.5) = 0;
obj_amp(obj_amp >= 0.5) = 1;
obj_amp = imresize(obj_amp,[N1,N2],'nearest');
obj_pha = zeros(size(obj_amp));
obj = obj_amp.*exp(1i*obj_pha);

bias = 0.02;  figw = 0.50;  figh = 0.40;
figure,set(gcf,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,2,[bias bias],[bias bias+0.04],[bias bias]);
ax = subplot(1,2,1);
ax.Position = pos{1};
imshow(abs(obj),[],'border','tight')
title('Amplitude of the object')
ax = subplot(1,2,2);
ax.Position = pos{2};
imshow(angle(obj),[],'border','tight');
title('Phase of the object')

% probe settings
radius = 100;
probe = propagate(aperture(N1,N2,N1/2,N2/2,radius),params.dist_1,params.pxsize,params.wavlen);
bias = 0.02;    figw = 0.50;    figh = 0.40;
figure,set(gcf,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,2,[bias bias],[bias bias+0.04],[bias bias]);
ax = subplot(1,2,1);
ax.Position = pos{1};
imshow(abs(probe),[],'border','tight')
title('Amplitude of the probe')
ax = subplot(1,2,2);
ax.Position = pos{2};
imshow(angle(probe),[],'border','tight');
title('Phase of the probe')

% probe positions
K1 = 6;             % number of positions (along x-axis)
K2 = 6;             % number of positions (along y-axis)
overlap = 0.8;      % overlapping ratio (between 0 and 1)
step = radius*params.pxsize*(1-overlap)*2;
[shifts_1,shifts_2] = meshgrid(linspace(-step*(K1-1)/2,step*(K1-1)/2,K1),linspace(-step*(K2-1)/2,step*(K2-1)/2,K2));
shifts_1 = shifts_1 + 5e-3*randn(size(shifts_1));   % add random offsets to avoid grid-like artifact
shifts_2 = shifts_2 + 5e-3*randn(size(shifts_2));   % add random offsets to avoid grid-like artifact

% calculate diffraction patterns
K = K1*K2;          % total number of measurements
snr_val = inf;      % signal-to-noise ratio
y = zeros(N1,N2,K); 
bias = 0.02;  figw = 0.30;  figh = 0.40;
figure,set(gcf,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,1,[bias bias],[bias bias+0.04],[bias bias]);

for k = 1:K
    exit_wave = probe.*imshift(obj,shifts_1(k)/params.pxsize, shifts_2(k)/params.pxsize);
    y(:,:,k) = abs(propagate(exit_wave,params.dist_2,params.pxsize,params.wavlen)).^2;
    y(:,:,k) = max(awgn(y(:,:,k),snr_val),0);
    im = imshow(y(:,:,k),[],'border','tight');
    im.Parent.Position = pos;
    title(['Image No.', num2str(k)])
    drawnow
end

%%
% =========================================================================
% Estimation
% =========================================================================

pos_error = 0.03;   % estimation error for the probe positions (mm)

% add perturbations to the probe positions
perturbs_1 = pos_error*rand(size(shifts_1)) - pos_error/2;
perturbs_2 = pos_error*rand(size(shifts_2)) - pos_error/2;
shifts_1_est = shifts_1 + perturbs_1;
shifts_2_est = shifts_2 + perturbs_2;

figure
set(gcf,'color','w','unit','normalized','position',[0.3,0.3,0.4,0.4])
plot([shifts_1(:),shifts_1_est(:)]',[shifts_2(:),shifts_2_est(:)]','-','linewidth',1,'color','k')
hold on
p1 = plot(shifts_1(:),shifts_2(:),'o','color','b','markerfacecolor','b');
hold on
p2 = plot(shifts_1_est(:),shifts_2_est(:),'o','color','r','markerfacecolor','r');
set(gca,'XDir','reverse')
xlabel('Position (mm)')
ylabel('Position (mm)')

legend([p1,p2],'True position','Estimated position','location','eastoutside')

%%
% =========================================================================
% Reconstruction
% =========================================================================

% initialization
filter = fspecial('gaussian',[50,50],50);    % Gaussian filtering to smooth the probe boundary
probe_est = imfilter(aperture(N1,N2,N1/2,N2/2,radius+20),filter);
obj_est  = ones(N1,N2);        % initial estimate of the object

% display settings
bias = 0.02;    figw = 0.70;    figh = 0.40;
figure,set(gcf,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,3,[bias bias],[bias bias],[bias bias]);

n_iters = 20;                   % number of iterations
pos_shift_dir = zeros(K,2);     % used to store the update direction of each probe position
beta = 10*ones(K,1);            % algorithm parameter (for probe position update)
alpha_1 = 1;                    % algorithm parameter (for object function update in ePIE)
alpha_2 = 0.2;                  % algorithm parameter (for probe function update in ePIE)

iter_probe = 1;                 % iteration to start probe function update
iter_position = 1;              % iteration to start probe position update
shift_max = 0.5;                % maximum allowable position shift in one iteration (pixel)
calib_channel = 'amplitude';    % channel (amplitude, phase, or both) for position calibration
disp_channel  = 'amplitude';    % display channel (amplitude or phase)

% main loop
for iter = 1:n_iters

    % traverse all positions in a random order
    for k = randperm(K)
        
        exit_wave = probe_est.*imshift(obj_est,shifts_1_est(k)/params.pxsize, shifts_2_est(k)/params.pxsize);
        u_est = propagate(exit_wave,params.dist_2,params.pxsize,params.wavlen);
            u_est = sqrt(y(:,:,k)).*exp(1i*angle(u_est));
        exit_wave_new = propagate(u_est,-params.dist_2,params.pxsize,params.wavlen);
        obj_est_new = obj_est + alpha_1 * conj(imshift(probe_est,-shifts_1_est(k)/params.pxsize,-shifts_2_est(k)/params.pxsize))./max(abs(probe_est(:)).^2) .* imshift(exit_wave_new - exit_wave,-shifts_1_est(k)/params.pxsize,-shifts_2_est(k)/params.pxsize);
        
        % probe function update
        if iter >= iter_probe
            probe_est_new = probe_est + alpha_2 * conj(imshift(obj_est,shifts_1_est(k)/params.pxsize, shifts_2_est(k)/params.pxsize))./max(abs(obj_est(:)).^2) .* (exit_wave_new - exit_wave);
            probe_est = probe_est_new;
        end
        
        % probe position update
        if iter >= iter_position

            window_calib = radius*2;    % window size for positional calibration
            
            obj_est_calib = obj_est(round(N1/2-window_calib/2 - shifts_2_est(k)/params.pxsize):round(N1/2+window_calib/2 - shifts_2_est(k)/params.pxsize)-1,...
                round(N2/2-window_calib/2 - shifts_1_est(k)/params.pxsize):round(N2/2+window_calib/2 - shifts_1_est(k)/params.pxsize)-1);
            obj_est_new_calib = obj_est_new(round(N1/2-window_calib/2 - shifts_2_est(k)/params.pxsize):round(N1/2+window_calib/2 - shifts_2_est(k)/params.pxsize)-1,...
                round(N2/2-window_calib/2 - shifts_1_est(k)/params.pxsize):round(N2/2+window_calib/2 - shifts_1_est(k)/params.pxsize)-1);
            
            % image registration
            if strcmpi(calib_channel,'phase')
                img_calib = angle(obj_est_calib);
                img_new_calib = angle(obj_est_new_calib); 
            elseif strcmpi(calib_channel,'amplitude')
                img_calib = abs(obj_est_calib);
                img_new_calib = abs(obj_est_new_calib);
            else
                img_calib = obj_est_calib;
                img_new_calib = obj_est_new_calib;
            end
            [shift,~] = dftregistration(fft2(img_calib),fft2(img_new_calib),1000);
            shift = max_shift(fliplr(shift(3:4)), shift_max);
            shifts_1_est(k) = shifts_1_est(k) - beta(k)*shift(1)*params.pxsize;
            shifts_2_est(k) = shifts_2_est(k) - beta(k)*shift(2)*params.pxsize;
            
            % accelerate convergence
            if iter >= iter_position
                pos_dir_new = [shift(1),shift(2)];
                if dot(pos_dir_new, pos_shift_dir(k,:)) > 0.3 * dot(pos_dir_new,pos_dir_new)
                    beta(k) = beta(k)*1.1;
                elseif dot(pos_dir_new, pos_shift_dir(k,:)) < -0.3 * dot(pos_dir_new,pos_dir_new)
                    beta(k) = beta(k)*0.9;
                end
            end
            pos_shift_dir(k,:) = [shift(1), shift(2)];

        end
        
        fprintf('Iteration %3d / %3d | Position error: %7.5e | beta: %5.2f | shift: %5.2e \n', ...
            iter, n_iters, std(shifts_2_est(:) - shifts_2(:)) + std(shifts_1_est(:) - shifts_1(:)), ...
            beta(k), sqrt((beta(k)*shift(1)*params.pxsize)^2 + (beta(k)*shift(2)*params.pxsize)^2));

        obj_est = obj_est_new;
        
        % display
        ax = subplot(1,3,1);
        ax.Position = pos{1};
        imshow(abs(probe_est),[]);
        ax = subplot(1,3,2);
        ax.Position = pos{2};
        if strcmpi(disp_channel,'phase')
            imshow(angle(obj_est),[]);
        else
            imshow(abs(obj_est),[]);
        end
        ax = subplot(1,3,3);
        ax.OuterPosition = pos{3};
        cla
        plot([shifts_1(:),shifts_1_est(:)]',[shifts_2(:),shifts_2_est(:)]','-','linewidth',1,'color','k')
        hold on,plot(shifts_1(:),shifts_2(:),'o','color','b','markerfacecolor','b')
        hold on,plot(shifts_1_est(:),shifts_2_est(:),'o','color','r','markerfacecolor','r')
        set(gca,'XDir','reverse')
        
        drawnow;
    end
end

%%
shifts_1_est = shifts_1 + perturbs_1;
shifts_2_est = shifts_2 + perturbs_2;

% display
        ax = subplot(1,3,1);
        ax.Position = pos{1};
        imshow(abs(imfilter(aperture(N1,N2,N1/2,N2/2,radius+20),filter)),[]);
        ax = subplot(1,3,2);
        ax.Position = pos{2};
        if strcmpi(disp_channel,'phase')
            imshow(angle(obj_est),[]);
        else
            imshow(abs(ones(size(obj_est))),[0,2]);
        end
        ax = subplot(1,3,3);
        ax.OuterPosition = pos{3};
        cla
        plot([shifts_1(:),shifts_1_est(:)]',[shifts_2(:),shifts_2_est(:)]','-','linewidth',1,'color','k')
        hold on,plot(shifts_1(:),shifts_2(:),'o','color','b','markerfacecolor','b')
        hold on,plot(shifts_1_est(:),shifts_2_est(:),'o','color','r','markerfacecolor','r')
        set(gca,'XDir','reverse')
        
        drawnow;
%%
% =========================================================================
% Auxiliary functions
% =========================================================================

function shift_new = max_shift(shift, val)
% =========================================================================
% This function ensures that the estimated shift does not exceed a given
% value.
% -------------------------------------------------------------------------
% Input:    - shift     : Estimated shift.
%           - val       : Threshold value (in pixel).
% Output:   - shift_new : Updated shift.
% =========================================================================
deno = abs(shift);
index_1 = (deno <= val);
index_2 = (deno >  val);
deno(index_1) = 1;
deno(index_2) = deno(index_2) / val;
shift_new = shift./deno;

end