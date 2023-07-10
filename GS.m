%% Iterative Fourier Transform Algorithm%
clear; % Clear all memory
close all;
clc;
% Loading the target image
% Create grid
N = 500;                % Number of grid points
L = 1;                  % Grid size
x = linspace(-L/2, L/2, N);
y = linspace(-L/2, L/2, N);
[X, Y] = meshgrid(x, y);

% Super-Gaussian profile parametersa

m = 16;                  % Super-Gaussian exponent
A = 100;               % Super-Gaussian amplitude
w0 = L/6;               % Super-Gaussian width

% Generate super-Gaussian profile
target = A * exp(-(abs(X)/w0).^(2*m)) .* exp(-(abs(Y)/(w0/4)).^(2*m));


% Gaussian beam parameters
w0 = L/6;               % Gaussian beam waist radius

% Calculate Gaussian beam profile
R = sqrt(X.^2 + Y.^2);  % Circular shape
initial_profile = exp(-((R/w0).^2));

%m=size(target,1); % Size of the image
%scale=N/m; % Estimate the necessary scaling factor
%target=imresize(target,scale); % Resize image to the matrix
%size
%target=double(target); % Convert symbolic object to a numeric
%object
target=target./(max(max(target))); % Normalize matrix
% Defining DOE phase
DOE=2*pi*rand(N,N); % Generate a N x N matrix of random phase
%between 0 and 2p
s=1000;
% IFTA algorithm
ax = gca;
for t=1:s; %Iterate to calculate the phase value
    DOEphase=exp(1i*DOE);
    % Forward iteration
    iterf=fft2(initial_profile.*DOEphase);
    intf=abs(iterf);
    angf=angle(iterf);
    A=target.*exp(1i*angf);

    % Backward iteration
    iterb=ifft2(A);
    angb=angle(iterb);
    DOE=angb;
    error=target-intf./max(max(intf)); %Calculate error
    intf=intf./max(max(intf));
    E=sum(sum(abs(error)))/(N*N);
    differences=target-intf;
    squaredDifferences = differences.^2;
    meanSquaredDifferences = mean(squaredDifferences,"all");
    rmse = sqrt(meanSquaredDifferences);
    if E<0.000005;%0.000005
        iteration=t;
        break
    end
%     plot(target(250,:))
%     hold on
%     plot(differences(250,:))
%     hold off
%     axis auto
    imagesc(differences);
    axis image;
    colormap jet;
    colorbar
    M{t} = getframe;
end

% f = figure(2)
% f.Position(3:4) = [1000 420]
% f,ax = subplot(1,3,1)
% ax.Position = [0.04 0.1 0.3 0.7]
% imagesc(initial_profile)
% f,ax = subplot(1,3,2)
% ax.Position = [0.34 0.1 0.3 0.7]
% imagesc(intf)
% f,ax = subplot(1,3,3)
% ax.Position = [0.64 0.1 0.3 0.7]
% imagesc(target)
% figure(1)

% for i = 1:size(M,2)
%     imagesc(M{i});
%     colormap jet;
%     colorbar
% end
%movie(M, -5);

movie2gif(M)
%V=VideoWriter('myavifile.avi');
%open(V);
%writeVideo(V,M)
%close(V);
function x = conjgrad(A, b, x)
    r = b - A * x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
            break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end