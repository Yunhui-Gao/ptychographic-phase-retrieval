function x = aperture(n1, n2, c1, c2, radius)

x = zeros(n1,n2);
for i = 1:n1
    for j = 1:n2
        if (i-c1)^2 + (j-c2)^2 <= radius^2
            x(i,j) = 1;
        end
    end
end

end

