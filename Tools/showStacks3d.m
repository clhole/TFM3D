function showStacks3d(in)

for i = 1:size(in,3)
    figure(1)
    imagesc(squeeze(in(:,:,i)))
    axis equal
    pause
end


for j = 1:max(size(in))
    figure(2)
    imagesc(squeeze(in(j,:,:)))
    axis equal
    pause
end


% for k = 1:max(size(in))
%     figure(4)
%     imagesc(squeeze(in(:,k,:)))
%     axis equal
%     pause
% end

end