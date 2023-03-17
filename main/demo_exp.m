% ========================================================================
% Introduction
% ========================================================================
% This code provides a simple demonstration of the extended ptychographical
% iterative engine (ePIE) using experimental data. This code is mainly
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
% Initialization
% =========================================================================
clear;clc
close all

% load functions
addpath(genpath('./utils'))

%%
% =========================================================================
% Auto-focusing
% =========================================================================

% load reference image for auto-focusing
prefix = 'VimbaImage_';
img_ref = im2double(imread(['../data/experiment/','ref.bmp']));

% select area of calculation
disp('Please select the area for diffraction calculation ...')
figure
[temp,rect_cal] = imcrop(img_ref);
if rem(size(temp,1),2) == 1
    rect_cal(4) = rect_cal(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect_cal(3) = rect_cal(3) - 1;
end
close
disp('Selected.')
img_obj_crop = imcrop(img_ref,rect_cal);

% select area of display
disp('Please select the area for display ...')
figure
[temp,rect_dis] = imcrop(img_obj_crop);
if rem(size(temp,1),2) == 1
    rect_dis(4) = rect_dis(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect_dis(3) = rect_dis(3) - 1;
end
close
disp('Selected.')

% physical parameters
params.pxsize = 2.740e-3;           % pixel size (mm)
params.wavlen = 0.532e-3;           % wavelength (mm)

%% determine the optimal object-to-sensor distance
bias = 0.02;    figw = 0.50;    figh = 0.40;
figure,set(gcf,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,2,[bias bias],[bias bias+0.04],[bias bias]);
for d = 3.69
    wavefront = propagate(sqrt(img_obj_crop),-d, params.pxsize, params.wavlen);
    ax = subplot(1,2,1);
    ax.Position = pos{1};
    imshow(imcrop(abs(wavefront), rect_dis),[],'border','tight');
    title(['Amplitude (d = ',num2str(d),' mm)'])
    ax = subplot(1,2,2);
    ax.Position = pos{2};
    imshow(imcrop(angle(wavefront), rect_dis),[],'border','tight');
    title(['Phase (d = ',num2str(d),' mm)'])
    drawnow;
    pause(1)
end

%% set the distance
params.dist = 3.69;     % object-to-sensor distance (mm)

%% load the entire dataset

K = 49; % total number of images

img_obj = im2double(imread(['../data/experiment/',prefix,'1.bmp']));

% select area of position calibration
disp('Please select the area for position calibration ...')
figure
[temp,rect_pos] = imcrop(img_obj);
if rem(size(temp,1),2) == 1
    rect_pos(4) = rect_pos(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect_pos(3) = rect_pos(3) - 1;
end
close
disp('Selected.')
img_obj = imcrop(img_obj,rect_pos);

shifts_relative = zeros(K-1,2);     % relative lateral shifts between adjacent images

% display & calculate
bias = 0.02;    figw = 0.50;    figh = 0.40;
figure,set(gcf,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,2,[bias bias],[bias bias+0.04],[bias bias]);
for k = 2:K
    img_obj_shift = im2double(imread(['../data/experiment/',prefix,num2str(k),'.bmp']));
    img_obj_shift = imcrop(img_obj_shift, rect_pos);
    [shift,~] = dftregistration(fft2(img_obj),fft2(img_obj_shift),100);
    shifts_relative(k-1,:) = fliplr(shift(3:4));
    ax = subplot(1,2,1);
    ax.Position = pos{1};
    imshow(img_obj,[],'border','tight');
    title('Previous image')
    ax = subplot(1,2,2);
    ax.Position = pos{2};
    imshow(img_obj_shift,[],'border','tight')
    title('Next image')
    drawnow;
    
    fprintf('Image No. %3d -> %3d | Relative shift = (Left: %7.2f px, Up: %7.2f px) \n', ...
            k-1, k, shift(4), shift(3));

    img_obj = img_obj_shift;

end

% calculate absolute shifts
shifts_absolute = zeros(K,2);
shifts_absolute(1,:) = [0,0];
for k = 2:K
    shifts_absolute(k,:) = shifts_absolute(k-1,:) + shifts_relative(k-1,:);
end

shifts_absolute(:,1) = shifts_absolute(:,1) - min(shifts_absolute(:,1));
shifts_absolute(:,2) = shifts_absolute(:,2) - min(shifts_absolute(:,2));
% figure,plot(shifts_absolute(:,1),shifts_absolute(:,2),'o','color','r','markerfacecolor','r')
labels = cell(K,1);
for k = 1:K
    labels{k} = num2str(k);
end
% text(shifts_absolute(:,1),shifts_absolute(:,2),labels,'color','r','VerticalAlignment','bottom','HorizontalAlignment','left')

shifts_normalized = shifts_absolute;
shifts_normalized(:,1) = shifts_normalized(:,1) - (max(shifts_normalized(:,1)) + min(shifts_normalized(:,1)))/2;
shifts_normalized(:,2) = shifts_normalized(:,2) - (max(shifts_normalized(:,2)) + min(shifts_normalized(:,2)))/2;

figure,plot(shifts_normalized(:,1),shifts_normalized(:,2),'o','color','r','markerfacecolor','r')
text(shifts_normalized(:,1),shifts_normalized(:,2),labels,'color','r','VerticalAlignment','bottom','HorizontalAlignment','left')

%% load experiment data

img_obj = im2double(imread(['../data/experiment/',prefix,'1.bmp']));

% select probe area
disp('Please select the image area ...')
figure
[temp,rect_aoi] = imcrop(img_obj);
if rem(size(temp,1),2) == 1
    rect_aoi(4) = rect_aoi(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect_aoi(3) = rect_aoi(3) - 1;
end
close
disp('Selected.')
temp = imcrop(img_obj,rect_aoi);
[M1,M2] = size(temp);

margin = 20;    % add additional spaces
N1 = ceil((max(shifts_normalized(:,1)) - min(shifts_normalized(:,1))) + M1 + margin);
N2 = ceil((max(shifts_normalized(:,2)) - min(shifts_normalized(:,2))) + M2 + margin);
if rem(N1,2) == 1
    N1 = N1+1;
end
if rem(N2,2) == 1
    N2 = N2+1;
end

% load images
y = zeros(N1,N2,K);
for k = 1:K
    img_obj = im2double(imread(['../data/experiment/',prefix,num2str(k),'.bmp']));
    y(N1/2-M1/2:N1/2+M1/2-1,N2/2-M2/2:N2/2+M2/2-1,k) = imcrop(img_obj,rect_aoi);
end

%% 
disp('Please select the probe area. Press Enter to confirm ...')
figure,imshow(y(:,:,1),[],'border','tight');
h = drawellipse('Center',[N2/2,N1/2],'SemiAxes',[M2/2,M1/2], 'RotationAngle',0,'StripeColor','m');
% wait for return/enter key press
flag = 0;
while flag ~= 1
    pause;
    key = get(gcf,'CurrentKey');
    if strcmpi(key,'return')
        flag = 1;
    else
        flag = 0;
    end
end
disp('Selected.')
probe_radius_est = mean(h.SemiAxes);    % initial estimate of the probe radius (pixel)
probe_center_est = h.Center;
close
% probe_radius_est = 500;
% probe_center_est = [N2/2, N1/2];

%% run the ePIE algorithm
                     
filter = fspecial('gaussian',[200,200],300);    % Gaussian filtering to smooth the probe boundary
probe_est = imfilter(aperture(N1,N2,probe_center_est(2),probe_center_est(1),probe_radius_est),filter);
obj_est = 0.5*ones(N1,N2);                     % initial estimate of the object

shifts_1_est = -shifts_normalized(:,1)*params.pxsize;
shifts_2_est = -shifts_normalized(:,2)*params.pxsize;
shifts_1_est_init = shifts_1_est;
shifts_2_est_init = shifts_2_est;

% display settings
bias = 0.02;    figw = 0.70;    figh = 0.40;
fig = figure;
set(fig,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,3,[bias bias],[bias bias],[bias bias]);

n_iters = 20;                   % number of iterations
pos_shift_dir = zeros(K,2);     % used to store the update direction of each probe position
beta = 10*ones(K,1);            % algorithm parameter (for probe position update)
alpha_1 = 1;                    % algorithm parameter (for object function update in ePIE)
alpha_2 = 0.2;                  % algorithm parameter (for probe function update in ePIE)

iter_probe = 1;                 % iteration to start probe function update
iter_position = 5;              % iteration to start probe position update
shift_max = 0.5;                % maximum allowable position shift in one iteration (pixel)
calib_channel = 'amplitude';    % channel (amplitude, phase, or both) for position calibration
disp_channel  = 'amplitude';    % display channel (amplitude or phase)

errs = zeros(n_iters,1);
disp('Start iteration ...')
% main loop
for iter = 1:n_iters

    % traverse all positions in a random order
    for k = randperm(K)
        
        exit_wave = probe_est.*imshift(obj_est,shifts_1_est(k)/params.pxsize, shifts_2_est(k)/params.pxsize);
        u_est = propagate(exit_wave,params.dist,params.pxsize,params.wavlen);
        errs(iter) = errs(iter) + 1/K*norm2(abs(u_est) - sqrt(y(:,:,k)))^2;
        u_est = sqrt(y(:,:,k)).*exp(1i*angle(u_est));
        exit_wave_new = propagate(u_est,-params.dist,params.pxsize,params.wavlen);
        obj_est_new = obj_est + alpha_1 * conj(imshift(probe_est,-shifts_1_est(k)/params.pxsize,-shifts_2_est(k)/params.pxsize))./max(abs(probe_est(:)).^2) .* imshift(exit_wave_new - exit_wave,-shifts_1_est(k)/params.pxsize,-shifts_2_est(k)/params.pxsize);
        
        % probe function update
        if iter >= iter_probe
            probe_est_new = probe_est + alpha_2 * conj(imshift(obj_est,shifts_1_est(k)/params.pxsize, shifts_2_est(k)/params.pxsize))./max(abs(obj_est(:)).^2) .* (exit_wave_new - exit_wave);
            probe_est = probe_est_new;
        end
        
        % probe position update
        if iter >= iter_position

            window_calib = 200;    % window size for positional calibration
            
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
        plot([shifts_1_est_init(:),shifts_1_est(:)]',[shifts_2_est_init(:),shifts_2_est(:)]','-','linewidth',1,'color','k')
        hold on,plot(shifts_1_est_init(:),shifts_2_est_init(:),'o','color','b','markerfacecolor','b')
        hold on,plot(shifts_1_est(:),shifts_2_est(:),'o','color','r','markerfacecolor','r')
        set(gca,'XDir','reverse')
            
        drawnow;

    end

    fprintf('Iteration %3d / %3d | Error: %7.2f \n', iter, n_iters, errs(iter));

end

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


function n = norm2(x)
n = norm(x(:),2);
end