%coded by jdy
% Load images
% psnr_DDFD=zeros(1,10);
% time_DDFD=zeros(1,10);
close all
 k=4;
    
    tic;
    % Load images
    I1=im2double(imread('images/lenag2.png'));%标准图像
    I2=im2double(imread('images/lenag1.png'));%浮动图像
    
    % Get the Key Points
    Options.upright=true;
    Options.tresh=0.0001;
    Ipts1=OpenSurf(I1,Options);
    Ipts2=OpenSurf(I2,Options);
    
    % Put the landmark descriptors in a matrix
    D1 = reshape([Ipts1.descriptor],64,[]);
    D2 = reshape([Ipts2.descriptor],64,[]);
    
    % Find the best matches
    err=zeros(1,length(Ipts1));
    cor1=1:length(Ipts1);
    cor2=zeros(1,length(Ipts1));
    for i=1:length(Ipts1),
        distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2,1);
        [err(i),cor2(i)]=min(distance);
    end
    
    % Sort matches on vector distance
    [err, ind]=sort(err);
    cor1=cor1(ind);
    cor2=cor2(ind);
    
    % Make vectors with the coordinates of the best matches
    Pos1=[[Ipts1(cor1).y]',[Ipts1(cor1).x]'];
    Pos2=[[Ipts2(cor2).y]',[Ipts2(cor2).x]'];
    Pos1=Pos1(1:30,:);
    Pos2=Pos2(1:30,:);
    
    % Show both images
    I = zeros([size(I1,1) size(I1,2)*2 size(I1,3)]);
    I(:,1:size(I1,2),:)=I1; I(:,size(I1,2)+1:size(I1,2)+size(I2,2),:)=I2;
    figure, imshow(I); hold on;
    
    % Show the best matches
    plot([Pos1(:,2) Pos2(:,2)+size(I1,2)]',[Pos1(:,1) Pos2(:,1)]','-');
    plot([Pos1(:,2) Pos2(:,2)+size(I1,2)]',[Pos1(:,1) Pos2(:,1)]','o');
    
    % Calculate affine matrix
    Pos1(:,3)=1; Pos2(:,3)=1;
    M=Pos1'/Pos2';
    
    % Add subfunctions to Matlab Search path
    functionname='OpenSurf.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/WarpFunctions'])
    
    % Warp the image
    I1_warped=affine_warp(I1,M,'bilinear');
    
    S=I1; M=I1_warped;
    % Velocity field smoothing kernel
    Hsmooth=fspecial('gaussian',[60 60],10);
    
    % The transformation fields
    Tx=zeros(size(M)); Ty=zeros(size(M));
    
    [Sy,Sx] = gradient(S);
    iteration=10;
    psnr1=zeros(1,iteration);
    alpha=2.5;
    
    for itt=1:iteration
   % Difference image between moving and static image
        Idiff=M-S;
        Relevancy=NPROD(M,S);   %相关系数  
        R=-(k^2)*log(Relevancy(1,2)); %相关度影响因子
        if R<1
            R=1;
        end
%         % Default demon force, (Thirion 1998)
%         Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
%         Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);

%         Extended demon force. With forces from the gradients from both
%         moving as static image. (Cachier 1999, He Wang 2005)
        [My,Mx] = gradient(M);
        Ux = -Idiff.*  Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)*R;
        Uy = -Idiff.*  My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)*R;
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;

        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);

        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
        
        M=movepixels(I1,Tx,Ty);        
    end
    N=M-I1;
    mi=MI(M,I1);
    Rcc=NPROD(M,I1);
    psnr=PSNR(M,I1);
    K = [0.01 0.03];
    window = fspecial('gaussian', 11, 1.5);
    L = 1;
    [mssim, ssim_map] = ssim_index(I1,M,K,window,L);
    mssim;
    toc;
    %psnr_DDFD(1,k)=psnr;
    %time_DDFD(1,k)=toc;
    

%PSNR
figure1 = figure;
plot(psnr_DDFD);

% Create xlabel
xlabel('k');

% Create ylabel
ylabel('PSNR');
%time
figure1 = figure;

% Create axes

plot(time_DDFD);

% Create xlabel
xlabel('k');

% Create ylabel
ylabel('time');