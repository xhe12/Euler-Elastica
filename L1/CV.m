% -----------------------------------------------------------
% -----------------------------------------------------------
% f------The Observed Picture(M*N);
% u------Segmentation Function(M*N);
% v------The introduced Auxillary Variable, v = u, (M*N);
% mu-----parameter;
% r1------penalty parameter 1;
% r2------penalty parameter 2;
% lambda1-----The Lagrangian Multiplier 1,
%              lambda1 = (lambda11, lambda12), (M*N,M*N);
% h------Laplacian Operator Kernel;
% FL-----Fourier Transformation of Laplacian, M*N;
% ------------------------------------------------------------
% ------------------------------------------------------------

  function [C,u,c0,c1] = CV(f,M,N,FL)
%%

% im1 = imread('ex4.png');
% im1 = rgb2gray(im1);
% f = im2double(im1);
% f = f*255;
% N = size(f,2); % # of columns, x direction;
% M = size(f,1); % # of rows, y direction;


mu = 100;
r1 = 100;
r2 = 100;
maxit = 20;

zero = zeros(M,N);

[FDx,FDy,BDx,BDy] = fun;

% L = zeros(M,N);
% L(2,1) = 1;
% L(1,2) = 1;
% L(M,1) = 1;
% L(1,N) = 1;
% L(1,1) = -4;
% FL = fft2(L);


% Initialize u, v, c0, c1, lambda1, lambda2;
u = (f-min(f(:)))/(max(f(:))-min(f(:)));
v = u;
c0 = sum(f(:).*v(:))/sum(v(:));
c1 = sum(f(:).*(1-v(:)))/sum(1-v(:));
p1 = zero;
p2 = zero;
% p1(:,1:(N-1)) = u(:,2:N) - u(:,1:(N-1));
% p1(:,N) = u(:,1) - u(:,N);
% p2(1:(M-1),:) = u(2:M,:) - u(1:(M-1),:);
% p2(M,:) = u(1,:) - u(M,:);

p1 = FDx(u);
p2 = FDy(u);

lambda11 = zero;
lambda12 = zero;
lambda2 = zero;



k =1;
%%
while  k < maxit


%%
    oldlambda11 = lambda11;
    oldlambda12 = lambda12;
    oldlambda2 = lambda2;



    % Calculate v;
    F = (f-c0).^2 - (f-c1).^2 - lambda2 - r2*u;
    F = -F/r2;
    v(0<=F<=1) = F;
    v(F>1) = 1;
    v(F<0) = 0;
    
    
    % Calculate the Fourier Transform of u by equation (25);
    
%     dlambda11dx(:,2:N) = lambda11(:,2:N) - lambda11(:,1:N-1);
%     dlambda11dx(:,1) = lambda11(:,1) - lambda11(:,N);
% 
%     dlambda12dy(2:M,:) = lambda12(2:M,:) - lambda12(1:M-1,:);
%     dlambda12dy(1,:) = lambda12(1,:) - lambda12(M,:);
%     
%     dp1dx(:,2:N) = p1(:,2:N) - p1(:,1:N-1);
%     dp1dx(:,1) = p1(:,1) - p1(:,N);
% 
%     dp2dy(2:M,:) = p2(2:M,:) - p2(1:M-1,:);
%     dp2dy(1,:) = p2(1,:) - p2(M,:);

    dlambda11dx = BDx(lambda11);
    dlambda12dy = BDy(lambda12);
    dp1dx = BDx(p1);
    dp2dy = BDy(p2);
    
    Fr = fft2( - r1*(dp1dx+dp2dy) - (dlambda11dx+dlambda12dy) ...
        + r2*v - lambda2 );
    Fr = Fr./( r2 - r1*FL );
           
        
    % Calculate u by inverse Fourier Transformation;
        u = ifft2(Fr);        
        u = real(u);

    
        %Calculate dudx and dudy;
%         dudx(:,1:N-1) = u(:,2:N) - u(:,1:N-1);
%         dudx(:,N) = u(:,1) - u(:,N);
%         dudy(1:M-1,:) = u(2:M,:) - u(1:M-1,:);
%         dudy(M,:) = u(1,:) - u(M,:);
        dudx = FDx(u);
        dudy = FDy(u);
   
        %Calculate q = (q1,q2);
        q1 = r1*dudx - lambda11;
        q2 = r1*dudy - lambda12;
        
        
        %Calculate p;
        q_abs = sqrt(q1.^2 + q2.^2);
        c = (1 - mu./q_abs)/r1;
        p1 = max(zero,c.*q1); 
        p2 = max(zero,c.*q2);
       
        % Update lambdas;
        lambda11 = lambda11 + r1*(p1 - dudx);
        lambda12 = lambda12 + r1*(p2 - dudy);
        lambda2 = lambda2 + r2*(u - v);
        
        %Update c0, c1;
        num = f.*v;
        num = sum(num(:));
        den = sum(v(:));
        c0 = num/den;
        
        num = f.*(1-v);
        num = sum(num(:));
        den = 1-v;
        den = sum(den(:));
        c1 = num/den;
        
    
      
        % Calculate Residues;
        
        d1 = p1-dudx;
        d1 = d1.^2;
        d2 = p2-dudy;
        d2 = d2.^2;
        R1(k) = sum(sqrt(d1(:)+d2(:)));

        R2(k) = sum(abs(u(:)-v(:)));
        
        % Calculate the lambda change;
        
        R = sqrt((lambda11-oldlambda11).^2 + (lambda12-oldlambda12).^2);
        L1(k) = sum(R(:));
        R = abs(lambda2-oldlambda2);
        L2(k) = sum(R(:));
        % Energy;
        
        E = (f-c0).^2.*u + (f-c1).^2.*(1-u) + mu*sqrt(dudx.^2+dudy.^2);
        Energy(k) = sum(E(:));
        
        k = k+1;

end

%%
% figure(2)
% mesh(u)
% 
% 
% X = 1:maxit-1;
% figure(3)
% plot(X,log(R1),'r', X, log(R2), 'b')
% legend('Residues1', 'Residue2')
% 
% figure(4)
% plot(X,log(L1),'r', X, log(L2), 'b')
% legend('Lagrangian1', 'Lagrangian2')
% 
% figure(5)
% plot(X,log(Energy))
% legend('Energy')
% 
% figure(9)
% mesh(abs(p1-FDx(u)))
% figure(10)
% mesh(abs(p2-FDy(u)))

figure(1)
imshow(f/255)
hold on
C = contour(u,[.99,1.01],'r');
hold off
  end