function out = shiftPSFAxial(in,n)

c = (n+1)/2;
kl = squeeze(in(c(1),c(2),:));


[Y, I] = max(kl');

dz = I-c(3);

out = zeros(size(in));

if dz > 0
    for i = 1:n(3)-dz
        out(:,:,i) = in(:,:,i+dz);
    end
elseif dz < 0
        for i = abs(dz)+1:n(3)
            out(:,:,i) = in(:,:,i+dz);
        end
else
    out = in;
end
             

end






