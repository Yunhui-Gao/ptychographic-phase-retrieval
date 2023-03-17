function img = imshift(img,s1,s2)

[n2,n1] = size(img);

f1 = -n1/2:1:n1/2-1;
f2 = -n2/2:1:n2/2-1;
[u1,u2] = meshgrid(f1,f2);

img = ifft2(fftshift(fftshift(fft2(img)).*exp(-1i*2*pi*(s1*u1/n1 + s2*u2/n2))));

end

