function [ u, exitflags, out, m, fit ] = GaussFit( m1, m2, hw)
%GAUSSFIT   determines the shift vector <u> of two 3D bead images <m1,m2>
%           by fitting a 3D gaussian <fit> to the correlation values <m>
%           with <hw> as the initial guess for the half-width of the gaussian.
%           Additional output: - <out> error messages of lsqfit
%                              - <exitflags> error flags of lsqfit
%Tobias Schoch, 2012
    
    warning off

    % Set step tolerance for lsqfit
    FitTol = 1e-4;

    winsize = transpose(size(m1));

    % FT of subvolumes in image 1
    Ff = fftn(m1,winsize.*2); % double padding
    % FT of subvolumes in image 2
    Fg = fftn(m2,winsize.*2); % double padding
    m = ifftn(conj(Ff).*Fg); 

    winpdsize = transpose(size(m));

    % Calculate mesh for padding
    ml = -ceil((winpdsize-1)./2);
    pl = floor((winpdsize-1)./2);
    x_win_ran = ml(1):pl(1);
    y_win_ran = ml(2):pl(2);
    z_win_ran = ml(3):pl(3);
    [x_win_pd,y_win_pd,z_win_pd] = ndgrid(x_win_ran,y_win_ran,z_win_ran); 
    % win_mesh = cat(4,x_win_pd,y_win_pd,z_win_pd);
    
    % Shift to center
    m = fftshift(m);
    
    % Prepare fitting function with 0 as zero line
    pfun = @(a,x)Gauss3D([a,0],x);

    % Prepare initial values for fitting
    a0 = zeros(1,7);
    a0(1:3) = 0; % displacement x,y,z
    a0(4:6) = hw; % half-width
    % a0(8) = min(reshape(m,numel(m),1)); % zero line

    [a0(7),I] = max(reshape(m,numel(m),1));
    [i,j,k] = ind2sub(size(m),I);
    a0(1:3) = [x_win_pd(i,j,k),y_win_pd(i,j,k),z_win_pd(i,j,k)]; % displacement is location of max.
    
    % Bounds for fit
    lb(1:3) = -ceil(winpdsize./2); %displacement
    ub(1:3) = ceil(winpdsize./2);
    
    lb(4:6) = 0.1; % half-width
    ub(4:6) = winpdsize./2;
    
    lb(7) = 0;  % amplitude
    up(7) = Inf;
    
    % Focus on peak to reduce computation time for fitting
    winsize = max(round(4*hw),[5,5,5]); % minimal size is 10px  
    % Take window of <winsize> around maximum peak
    m_red = GetWindowAt(m,winsize,[i,j,k]);
    l = [-floor(winsize'/2),floor(winsize'/2)-(1-mod(winsize',2))];
    i_ran = l(1,1):l(1,2);
    j_ran = l(2,1):l(2,2);
    k_ran = l(3,1):l(3,2);
    
    % Calculate corresponding ranges
    [x_win_pd,y_win_pd,z_win_pd] = ndgrid(i_ran+x_win_ran(i),j_ran+y_win_ran(j),k_ran+z_win_ran(k));  
    win_mesh_red = cat(4,x_win_pd,y_win_pd,z_win_pd);

    % Fitting
    opts = optimset('Display','off','MaxFunEvals',2000,'MaxIter',3000,'TolFun',FitTol,'TolX',FitTol); %,'Algorithm','levenberg-marquardt');
    [a,resnorm,residual,exitflags,out] = lsqcurvefit(pfun,a0,win_mesh_red,m_red,lb,ub,opts);
    % [a,resnorm,residual,exitflags,out]=lsqcurvefit(pfun,a0,win_mesh,m,lb,ub,opts);
    u = a(1:3);
    fit = pfun(a,win_mesh_red);
    
end

function [ y ] = Gauss3D( a, x )
%GAUSS3D caclulates the 3D-gauss function

	y = a(8)+a(7).*exp(-((x(:,:,:,1)-a(1)).^2)/(2*a(4))^2-((x(:,:,:,2)-a(2)).^2)/(2*a(5))^2-((x(:,:,:,3)-a(3)).^2)/(2*a(6))^2);

end