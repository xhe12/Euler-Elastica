% -----------------------------------------------------------
% -----------------------------------------------------------
% f------The Observed Picture(M*N);
% phi------Signed Distance for Segamentation(M*N);
% p------The introduced Auxillary Variable, p = (phi_x, phi_y), (M*N,M*N);
% n------The introduced Auxillary Variable, n = (p_x, p_y)/|p|, (M*N,M*N);
% c1-----Fitting grey intensity inside;
% c2-----Fitting grey intensity outside;
% a-----parameter;
% b-----parameter;
% r1------penalty parameter 1;
% r2------penalty parameter 2;
% r3------penalty parameter 3;
% lambda1-----The Lagrangian Multiplier 1,
%              lambda1 = (lambda11, lambda12), (M*N,M*N);
% lambda2-----The Lagrangian Multiplier 2,
%              lambda2 = (lambda21, lambda22), (M*N,M*N);
% lambda3-----The Lagrangian Multiplier 3, (M*N);
% h------Laplacian Operator Kernel;
% FL-----Fourier Transformation of Laplacian, M*N;
% ------------------------------------------------------------
% ------------------------------------------------------------

im1 = imread('pizza kid.jpg');
im1 = rgb2gray(im1);
%im1 = im1(50:end-100,100:end-100);
f = im2double(im1);
f = f*255;
N = size(f,2); % # of columns, x direction;
M = size(f,1); % # of rows, y direction;
area = M*N;
a =.001;
b =1;
r1 = .5;
r2 = .4;
r3 = .1;
r4 = 8;
eps = 1;
epsilon = 10^(-6);
maxit = 300;

iteration = 0;



[FDx,FDy,BDx,BDy] = fun;

% Calculate the Fourier Transform of Laplacian;
L = zeros(M,N);
L(2,1) = 1;
L(1,2) = 1;
L(M,1) = 1;
L(1,N) = 1;
L(1,1) = -4;
FL = fft2(L);

% % Initialize phi;
 [~,u,~,~] = CV(f,M,N,FL);
index = zeros(M,N);
index(u<0.99) = 1;

phi = bwdist(index);

index = zeros(M,N);
index(u>=0.99) = 1;
phi = bwdist(index)*(-1)+phi;


deltat = 0.01;


for k = 1:100
    oldphi = phi;
    aa = oldphi(2:end-1,2:end-1)-oldphi(2:end-1,1:end-2);
    bb = oldphi(2:end-1,3:end)-oldphi(2:end-1,2:end-1);
    cc = oldphi(2:end-1,2:end-1)-oldphi(1:end-2,2:end-1);
    dd = oldphi(3:end,2:end-1)-oldphi(2:end-1,2:end-1);

    splus = oldphi(2:end-1,2:end-1)./sqrt(oldphi(2:end-1,2:end-1).^2 ...
          + (1-sqrt(bb.^2+dd.^2)).^2);
    
    sminus = oldphi(2:end-1,2:end-1)./sqrt(oldphi(2:end-1,2:end-1).^2 ...
          + (1-sqrt(aa.^2+cc.^2)).^2);

    ap = max(aa,0);
    am = min(aa,0);
    bp = max(bb,0);
    bm = min(bb,0);
    cp = max(cc,0);
    cm = min(cc,0);
    dp = max(dd,0);
    dm = min(dd,0);

    phi(2:end-1,2:end-1) = oldphi(2:end-1,2:end-1) ...
                         - deltat*splus.*(sqrt(ap.^2+bm.^2+cp.^2+dm.^2)-1) ...
                         - deltat*sminus.*(sqrt(am.^2+bp.^2+cm.^2+dp.^2)-1);
end 


f = f/255;

figure(1)
imshow(f)
hold on
contour(phi,[-.1,.1],'r')
% Initialize p,n,q,lambdas;
p1 = zeros(M,N);
p2 = zeros(M,N);
absp = zeros(M,N);
n1 = zeros(M,N);
n2 = zeros(M,N);
m1 = zeros(M,N);
m2 = zeros(M,N);
q = zeros(M,N);
% p1 = BDx(phi);
% p2 = BDy(phi);
% m1 = p1;
% m2 = p2;
% absp = sqrt(p1.^2+p2.^2)+epsilon;
% n1 = p1./absp;
% n2 = p2./absp;
% q = BDx(n1)+BDy(n2);

lambda11 = zeros(M,N);
lambda12 = zeros(M,N);
lambda21 = zeros(M,N);
lambda22 = zeros(M,N);
lambda31 = zeros(M,N);
lambda32= zeros(M,N);
lambda4 = zeros(M,N);

% mesh(phi)

% Construct g1, g2, g3, g4;
g1 = zeros(M,N);
g2 = g1;
g3 = g1;
g4 = g1;

g1(1,2) = 1;
g1(1,N) = 1;
g1(1,1) = -2;
Fg1 = fft2(g1);

g2(1,N) = 1;
g2(2,1) = 1;
g2(2,N) = -1;
g2(1,1) = -1;
Fg2 = fft2(g2);

g3(M,1) = 1;
g3(1,2) = 1;
g3(M,2) = -1;
g3(1,1) = -1;
Fg3 = fft2(g3);

g4(M,1) = 1;
g4(2,1) = 1;
g4(1,1) = -2;
Fg4 = fft2(g4);

delta = .01;

%%
for k = 1:maxit

%%  
    k
    % Update c0 and c1;
    dem = .5*(1 + 2/pi*atan(phi/eps));% Denominator
    dem = dem(:);
    num = f(:).*dem;% Numerator
    c0 = sum(num)/sum(dem);

    dem = .5*(1 - 2/pi*atan(phi/eps));% Denominator
    dem = dem(:);
    num = f(:).*dem;% Numerator
    c1 = sum(num)/sum(dem);


    % Update phi;

    F = ((f-c0).^2-(f-c1).^2)/pi;
    G = (a + b*abs(q)).*absp/pi;
    
     oldphi = phi;
%     delta = .01;
    l = 0;
    tol = 1;
     while tol > 10^(-3) && l<=50
         oldphi1 = phi;
        H = eps*F./(eps^2+phi.^2) - 2*eps*phi.*G./(eps^2+phi.^2).^2 ...
            + r1*(BDx(p1)+BDy(p2)) + (BDx(lambda11)+ BDy(lambda12))-delta*phi;
        Fphi = fft2(H)./(r1*FL-delta);
        phi = real(ifft2(Fphi));
        tol = sum(abs((phi(:)-oldphi1(:))))./sum(abs(oldphi1(:)));% Relative Error
        l = l+1;
     end

      
%         dif = eps*F./(eps^2+phi.^2) - 2*eps*phi.*G./(eps^2+phi.^2).^2 ...
%                + r1*(BDx(p1)+BDy(p2)) + BDx(lambda11)+BDy(lambda12)...
%                -r1*(BDx(FDx(phi))+BDy(FDy(phi)));
%         dif_max = sum(abs(dif(:)))/area
%%
     %% Update q;

    v = r4*(BDx(n1)+BDy(n2))-lambda4;
    w = b*eps./((eps^2+phi.^2)*pi).*absp;
    q = max((abs(v)-w)./abs(v),0).*v/r4;


    % Update p;
    abs_n = sqrt(n1.^2+n2.^2);
    eta = eps*(a+b*abs(q))./(pi*(eps^2+phi.^2)) + lambda21.*n1 ...
          + lambda22.*n2;
    mu = r1 + r2*abs_n.^2 + r2 + r3;
    a1 =  (r1*FDx(phi)-lambda11+lambda21+r3*m1-lambda31)./mu;
    a2 =  (r1*FDy(phi)-lambda12+lambda22+r3*m2-lambda32)./mu;
    abs_a = sqrt(a1.^2+a2.^2);
    nu1 = - r2*n1;
    nu2 = - r2*n2;
    abs_nu = sqrt(nu1.^2+nu2.^2);
    
    epsilon = 10^(-6);
    
    % eta>=mu.*abs_a;
    index = (eta>=mu.*abs_a);
    p1(index) = 0;
    p2(index) = 0;
    
    % eta<mu.*abs_a & a==0 & nu==0 & eta>=0 
    index = 1-index;
    index1 = (abs_a<epsilon).*(abs_nu<epsilon);
    index2 = (eta>=0);
    cond = index.* index1.* index2;
    cond = logical(cond);
    p1(cond) = 0;
    p2(cond) = 0;
    
    % eta<mu.*abs_a & a==0 & nu==0 & eta<0 
    index2 = 1 - index2;
    cond = index.* index1.* index2;
    cond = logical(cond);
    p1(cond) = -sqrt(2)/2*eta(cond)./mu(cond);
    p2(cond) = -sqrt(2)/2*eta(cond)./mu(cond);
    
    % eta<mu.*abs_a & a~=0 & nu==0
    index1 = (abs_a>=epsilon).*(abs_nu<epsilon);
    cond = index.* index1;
    cond = logical(cond);
    coef = 1 - eta(cond)./(mu(cond).*abs_a(cond));
    p1(cond) = coef.*a1(cond);
    p2(cond) = coef.*a2(cond);
    
    % eta<mu.*abs_a & a==0 & nu~=0
    index1 = (abs_a<epsilon).*(abs_nu>=epsilon);
    cond = index.* index1;
    cond = logical(cond);
    coef = eta(cond)./((mu(cond)-2*abs_nu(cond)).*abs_nu(cond));
    p1(cond) = coef.*nu1(cond);
    p2(cond) = coef.*nu2(cond);
    
    % eta<mu.*abs_a & a~=0 & nu~=0
    index1 = (abs_a>=epsilon).*(abs_nu>=epsilon);
    cond = index.* index1;
    cond = logical(cond);
    nu1_cond = nu1(cond);
    nu2_cond = nu2(cond);
    a1_cond = a1(cond);
    a2_cond = a2(cond);
    abs_a_cond = abs_a(cond);
    abs_nu_cond = abs_nu(cond);
    mu_cond = mu(cond);
    eta_cond = eta(cond);
    
    cos_alpha = (a1_cond.*nu1_cond+a2_cond.*nu2_cond)...
                  ./(abs_a_cond.*abs_nu_cond);
              
    index = (cos_alpha>1);
    cos_alpha(index) = 1;
    index = (cos_alpha<-1);
    cos_alpha(index) = -1;
    
    alpha = acos(cos_alpha);
    lt = mu_cond.*abs_a_cond.*abs_nu_cond.*sin(alpha);% Last term of
                                                      % the equation;
%      theta = acos(a1_cond./abs_a_cond);%Initialize theta;
% %    theta = acos(nu1_cond./abs_nu_cond);%Initialize theta;
    
 
%         y = (mu_cond(i)*abs_a_cond(i)*cos(-pi:.01:pi)-eta_cond(i)).^2 ...
%             ./(mu_cond(i)+2*abs_nu_cond(i)*cos((-pi:.01:pi)+alpha(i)));
    theta = zeros(size(alpha,1),1);
        
    
    tol = 1;
    j = 1;
    
    while tol>10^(-3) & j<=10
        oldtheta = theta;
        
        Func = mu_cond.*sin(theta).*abs_a_cond.*(mu_cond...
               + abs_nu_cond.*cos(theta+alpha))...
               + eta_cond.*abs_nu_cond.*sin(theta+alpha)...
               - lt;
           
        DFunc = mu_cond.*abs_a_cond.*(mu_cond.*cos(theta)...
               + abs_nu_cond.*cos(2*theta+alpha))...
               + eta_cond.*abs_nu_cond.*cos(theta+alpha);
        
        theta = theta - Func./DFunc;
        
        rel =  sum(abs((oldtheta-theta)))./sum(abs(theta)); %Relative Error;
        tol = max(rel);

        j = j+1;
    end
    
%         Func = mu_cond.*sin(theta).*abs_a_cond.*(mu_cond...
%                + abs_nu_cond.*cos(theta+alpha))...
%                + eta_cond.*abs_nu_cond.*sin(theta+alpha)...
%                - lt;
%     
%      dif_max = max(abs(Func))
%     
    Det = nu1_cond.*a2_cond - nu2_cond.*a1_cond;
    
    theta_tilda = (Det>=0).* theta + (Det<0).*(-theta);% May need to check;
    
    b1 = (cos(theta_tilda).*a1_cond - sin(theta_tilda).*a2_cond)...
         ./abs_a_cond;
    b2 = (sin(theta_tilda).*a1_cond + cos(theta_tilda).*a2_cond)...
         ./abs_a_cond;
    
    coef = (mu_cond.* (b1.* a1_cond + b2.* a2_cond) - eta_cond)...
           ./(mu_cond + 2* (nu1_cond.* b1 + nu2_cond.* b2));
    
    p1(cond) = coef.* b1;
    p2(cond) = coef.* b2;
    

    %update m;
    m1 = p1+lambda31/r3;
    m2 = p2+lambda32/r3;
    abs_m = sqrt(m1.^2+m2.^2)+epsilon;
    m1(abs_m> .01)=m1(abs_m>.01)./abs_m(abs_m>.01);
    m2(abs_m>.01)=m2(abs_m>.01)./abs_m(abs_m>.01);

%     m1 = m1./abs_m;
%     m2 = m2./abs_m;
%%
    % Update n;
%     AP1 = [p1(:,1)+p1(:,N), p1(:,2:N)+p1(:,1:N-1)];
%     AP1 = AP1/2;
%     AP2 = [p2(1,:)+p2(M,:);p2(2:M,:)+p2(1:M-1,:)];
%     AP2 = AP2/2;
%     A_absp = sqrt(AP1.^2 + AP2.^2);
   absp = sqrt(p1.^2+p2.^2);
    
    D = max(r2*absp(:).^2);
    
    coef1 = -r2*absp.*p1 + absp.*lambda21 + r4*FDx(q) + FDx(lambda4);
           
    coef2 = -r2*absp.*p2 + absp.*lambda22 + r4*FDy(q) + FDy(lambda4);
    
    n1old = n1;
    n2old = n2;
%     tol = 1;
%     i = 1;
%     while tol>10^(-6) & i<=100

    
        a11 = r4 * Fg1 - D;
        a12 = r4 * Fg2;
        f1 = r2*absp.^2.*n1 + coef1 - D*n1;
        Ff1 = fft2(f1);

        a21 = r4 * Fg3;
        a22 = r4 * Fg4 - D;
        f2 = r2*absp.^2.*n2 + coef2 - D*n2;
        Ff2 = fft2(f2);

        Det = a11.* a22 - a12.* a21;
        Fn1 = (a22.*Ff1 - a12.*Ff2)./Det;
        Fn2 = (-a21.*Ff1 + a11.*Ff2)./Det;
        n1 = real(ifft2(Fn1));
        n2 = real(ifft2(Fn2));

%         tol = sum(sqrt((n1(:) - n1old(:)).^2 + (n2(:) - n2old(:)).^2))...
%             /area;
%             
%         i = i+1;
% 
%     end


%     dif1 = -FDx(coef1.*(BDx(n1)+BDy(n2))) + (r2*(absp.*n1-p1)+lambda21).*absp;
%     dif2 = -FDy(coef2.*(BDx(n1)+BDy(n2))) + (r2*(absp.*n2-p2)+lambda22).*absp;
%     dif1max = max(abs(dif1(:)))
%     dif2max = max(abs(dif2(:)))
%     dif1max = sum(abs(dif1(:)))/area
%     dif2max = sum(abs(dif2(:)))/area
   



%%
    % Update lambdas;
    lambda11old = lambda11;
    lambda12old = lambda12;
    lambda21old = lambda21;
    lambda22old = lambda22;
    lambda31old = lambda31;
    lambda32old = lambda32;
    lambda4old = lambda4;
    
    lambda11 = lambda11 + r1*(p1 - FDx(phi));
    lambda12 = lambda12 + r1*(p2 - FDy(phi));
    lambda21 = lambda21 + r2*(absp.*n1 - p1);
    lambda22 = lambda22 + r2*(absp.*n2 - p2);
    lambda31 = lambda31 + r3*(p1 - m1);
    lambda32 = lambda32 + r3*(p2 - m2);
    lambda4 = lambda4 + r4*(q - (BDx(n1)+BDy(n2)));
    

    
    %----------------------------------------------------------------------
    %-------------------- Relative Errors 1 through 4 ---------------------
    L1(k) = sum(sqrt((lambda11(:)-lambda11old(:)).^2 ...
            + (lambda12(:)-lambda12old(:)).^2))...
            /sum(sqrt(lambda11old(:).^2+lambda12old(:).^2));
    
    L2(k) = sum(sqrt((lambda21(:)-lambda21old(:)).^2 ...
            + (lambda22(:)-lambda22old(:)).^2))...
            /sum(sqrt(lambda21old(:).^2+lambda22old(:).^2));
    
    L3(k) = sum(sqrt((lambda31(:)-lambda31old(:)).^2 ...
            + (lambda32(:)-lambda32old(:)).^2))...
            /sum(sqrt(lambda31old(:).^2+lambda32old(:).^2));
    L4(k) = sum(sqrt((lambda4(:)-lambda4old(:)).^2)) ...
        /sum(sqrt(lambda4old(:).^2));
    
    %----------------------------------------------------------------------
    %----------------------- Residues 1 through 4 -------------------------

    R = sqrt((p1 - FDx(phi)).^2 +(p2 - FDy(phi)).^2);
    R1(k) = sum(R(:))/area;
    R = sqrt((absp.*n1 - p1).^2 + (absp.*n1 - p1).^2);
    R2(k) = sum(R(:))/area;
    R = sqrt((p1 - m1).^2 +(p2 - m2).^2);
    R3(k) = sum(R(:))/area;
    R = abs(q-(BDx(n1)+BDy(n2)));
    R4(k) = sum(R(:))/area;
%     dif_R11 = p1-FDx(phi);
%     dif_R12 = p2-FDy(phi);
%     figure(1)
%     mesh(dif_R11)
%     figure(2)
%     mesh(dif_R12)
    
    %----------------------------------------------------------------------
    %------------------------------ Energy --------------------------------
    Hphi = .5*(1 + 2/pi*atan(phi/eps));
    Fit = (f - c0).^2.*Hphi + (f - c1).^2.*(1-Hphi);%Fitting Term
    abs_deltaphi = sqrt(FDx(phi).^2 + FDy(phi).^2);
    Ra1 = FDx(phi)./abs_deltaphi;
    Ra2 = FDy(phi)./abs_deltaphi;
    R = (a + abs(BDx(Ra1)+BDy(Ra2))*b).*eps./(pi*(eps^2+phi.^2)).*absp;
%     abs_Hphi = sqrt(FDx(Hphi).^2+FDy(Hphi).^2);
%     R = (a + (BDx(Ra1)+BDy(Ra2)).^2*b).*abs_Hphi;
    E(k) = sum(R(:)+Fit(:));
    disp = E(k)
    
    %----------------------------------------------------------------------
    %-------------------------- Relative Error ----------------------------
    RE(k) = sum(abs((phi(:)-oldphi(:)))) ...
        /sum(abs(oldphi(:)));
%%
iteration = iteration + 1;
end



figure(2)
mesh(Hphi)
x = 1:maxit;

figure(3)
mesh(phi)

figure(4)
plot(x,log(L1(1:maxit)),'r', x,log(L2(1:maxit)),'g',...
    x,log(L3(1:maxit)),'b', x, log(L4(1:maxit)),'c')
legend('L1','L2','L3','L4')

figure(5)
plot(x,log(R1(1:maxit)),'r', x,log(R2(1:maxit)),'g',...
        x,log(abs(R3(1:maxit))),'b', x, log(abs(R4(1:maxit))),'c')
legend('Residue1','Residue2','Residue3','Residue4')

figure(6)
plot(x,log(E(1:maxit)))
legend('Energy')

figure(7)
plot(x,log(RE(1:maxit)))
legend('Relative error in phi')

figure(8)
imshow(f)
hold on
contour(phi,[-.01,.01],'r')
title(['b = ',num2str(b),' r1= ',num2str(r1),' r2 = ',num2str(r2),...
    ' r3 = ', num2str(r3),...
    ' r4 = ', num2str(r4),...
    ' epsilon = ', num2str(eps),...
    ' iteration = ', num2str(iteration)])
hold off

% dh = 1;
% % eps = 1.5
% 
% absxplus = sqrt((phi(2:end-1,3:end)-phi(2:end-1,2:end-1)).^2 ...
%     + 1/4*((phi(3:end,2:end-1)-phi(1:end-2,2:end-1))/2 ...
%           +(phi(3:end,3:end)-phi(1:end-2,3:end))/2).^2);
%           
% Dxplus = (phi(2:end-1,3:end)-phi(2:end-1,2:end-1))./absxplus;
%           
% 
% absxminus = sqrt((phi(2:end-1,2:end-1)-phi(2:end-1,1:end-2)).^2 ...
%     + 1/4*((phi(3:end,1:end-2)-phi(1:end-2,1:end-2))/2 ...
%           +(phi(3:end,2:end-1)-phi(1:end-2,2:end-1))/2).^2);
%           
% Dxminus = (phi(2:end-1,2:end-1)-phi(2:end-1,1:end-2))./absxminus;
%           
% 
% absyplus = sqrt((phi(3:end,2:end-1)-phi(2:end-1,2:end-1)).^2 ...
%     + 1/4*((phi(3:end,3:end)-phi(3:end,1:end-2))/2 ...
%           +(phi(2:end-1,3:end)-phi(2:end-1,1:end-2))/2).^2);
%           
% Dyplus = (phi(3:end,2:end-1)-phi(2:end-1,2:end-1))./absyplus;
%     
% 
% absyminus = sqrt((phi(2:end-1,2:end-1)-phi(1:end-2,2:end-1)).^2 ...
%     + 1/4*((phi(2:end-1,3:end)-phi(2:end-1,1:end-2))/2 ...
%           +(phi(1:end-2,3:end)-phi(1:end-2,1:end-2))/2).^2);
%           
% Dyminus = (phi(2:end-1,2:end-1)-phi(1:end-2,2:end-1))./absyminus;
%           
% 
% kappa = (Dxplus-Dxminus + Dyplus-Dyminus)/dh;
% abs_delta = 1/4*(absxplus+absxminus+absyplus+absyminus)/dh;
% 
% H = 1/pi*atan(phi/eps);
% 
% DHx = (H(2:end-1,3:end)-H(2:end-1,1:end-2))/2;
% DHy = (H(3:end,2:end-1)-H(1:end-2,2:end-1))/2;
% abs_DH = sqrt(DHx.^2+DHy.^2);
% 
% bd_curvature = abs(kappa).*abs_DH;
% 
% E = sum(bd_curvature(:))