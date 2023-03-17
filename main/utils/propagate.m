function [w_o] = propagate(w_i, dist, pxsize, wavlen)

% dimension of the wavefield
[ny,nx] = size(w_i);

% sampling in the frequency domain
kx = pi/pxsize*(-1:2/nx:1-2/nx);
ky = pi/pxsize*(-1:2/ny:1-2/ny);
[KX,KY] = meshgrid(kx,ky);

% wave number
k = 2*pi/wavlen;

% circular convoluion via ffts
inputFT = fftshift(fft2(w_i));
H = exp(1i*dist*sqrt(k^2-KX.^2-KY.^2));
w_o = ifft2(fftshift(inputFT.*H));

end

