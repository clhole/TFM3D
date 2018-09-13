function [ du ] = GradientAt( u, delta, n , ind )
%GRADIENTAT   calculates the gradient averaged over a n*n*n-window
%             located at ind=[i,j,k]
%Tobias Schoch, 2012

% adapt inconvenient n
l = floor(n/2);
N = size(u);
N = N(1:3);

ld = ones(2,3)*l;

% % check for border and eventually reduce window size
for d=1:3
    a=ind(d)-l;
    if a<=l
        ld(1,d)=a;
    end
    a=N(d)-ind(d);
    if a<=l
        ld(2,d)=a;
    end
end
% for d=1:3
%     a=l-ind(d)+1;
%     if a>0
%         ld(1,d)=l-a;
%         ld(2,d)=l+a;
%     end
%     a=l+ind(d)-N(d);
%     if a>0
%         ld(1,d)=l+a;
%         ld(2,d)=l-a;
%     end
% end

n = 1+ld(1,:)+ld(2,:);

% create window around location
rangeX = [-ld(1,1):ld(2,1)];
rangeY = [-ld(1,2):ld(2,2)];
rangeZ = [-ld(1,3):ld(2,3)];
win = u(ind(1)+rangeX,ind(2)+rangeY,ind(3)+rangeZ,:);

% calculate all common values
M = zeros(3,1);
M(1) = (delta(1)^2*n(2)*n(3))*sum(rangeX.^2);
M(2) = (delta(2)^2*n(1)*n(3))*sum(rangeY.^2);
M(3) = (delta(3)^2*n(2)*n(1))*sum(rangeZ.^2);

m1 = permute(sum(sum(win,3),2),[4 1 2 3]);
su = sum(m1,2);
xu = delta(1).*sum((ones(3,1)*rangeX).*m1,2);
yu = delta(2).*sum((ones(3,1)*rangeY).*permute(sum(sum(win,3),1),[4 2 1 3]),2);
zu = delta(3).*sum((ones(3,1)*rangeZ).*permute(sum(sum(win,2),1),[4 3 1 2]),2);

% initialize and fill diagonals
A = diag([M;prod(n)]);
b = ones(4,1);

% fill simple sums
tmp = transpose(delta).*[n(2)*n(3)*sum(rangeX),n(1)*n(3)*sum(rangeY),n(1)*n(2)*sum(rangeZ)];
A(4,1:3)=tmp;
A(1:3,4)=transpose(tmp);

% fill other elements
xy=delta(1).*delta(2).*n(3).*sum(rangeX)*sum(rangeY);
yz=n(1).*delta(2).*delta(3).*sum(rangeY)*sum(rangeZ);
xz=delta(1).*n(2).*delta(3).*sum(rangeX)*sum(rangeZ);

A(1,2)=xy;
A(2,1)=xy;
A(1,3)=xz;
A(3,1)=xz;
A(2,3)=yz;
A(3,2)=yz;

% preallocate 
du = zeros(3);

% gradient for d th component
for i=1:3
    % fill component dependent b
    b(1)=xu(i);
    b(2)=yu(i);
    b(3)=zu(i);
    b(4)=su(i);

    % solve problem
    opts.SYM = true;
    x=linsolve(A,b,opts);
    du(:,i)=x(1:3);
end


end

