%% Load images sequence
seqdir = 'DataSeq1';
seqpath = fullfile('input', seqdir);
fileimages = dir(fullfile(seqpath, '*.jpg'));

for ii = 1:numel(fileimages)
    images{ii} = imread(fullfile(seqpath, fileimages(ii).name));
    if size(images{ii}, 3) == 3; images{ii} = rgb2gray(); end;
    images{ii} = im2double(images{ii});
end

%% Lucas-Kanade Optic Flow
% compute gradient using 2x2 prewitt operator
im1 = images{1};
im2 = images{2};
Ix = conv2(im1, [-1 1; -1 1], 'same');
Iy = conv2(im1, [-1 -1; 1 1], 'same');
It = conv2(im1, ones(2), 'same') - conv2(im2, ones(2), 'same');

window_sz = 5;
sumIx2 = conv2(Ix.^2, ones(window_sz), 'valid');
sumIy2 = conv2(Iy.^2, ones(window_sz), 'valid');
sumIxy = conv2(Ix.*Iy, ones(window_sz), 'valid');
sumIxt = conv2(Ix.*It, ones(window_sz), 'valid');
sumIyt = conv2(Iy.*It, ones(window_sz), 'valid');

rows = size(sumIx2, 1); cols = size(sumIx2, 2);
u = zeros(rows,cols); v = zeros(rows, cols);
tic;for ii = 1:rows
    for jj = 1:cols
        A = [sumIx2(ii,jj) sumIxy(ii,jj); sumIxy(ii,jj) sumIy2(ii,jj)];
        b = [-sumIxt(ii,jj); -sumIyt(ii,jj)];
        x = A\b;
        u(ii,jj) = x(1); v(ii,jj) = x(2);
    end;
end; toc;

u = padarray(u, [2 2]);
v = padarray(v, [2 2]);
%% Visualize with quivers
close all;

% downsize u and v
u = u*2;
v = v*2;

u_deci = u(1:5:end, 1:5:end);
v_deci = v(1:5:end, 1:5:end);
[m, n] = size(im2);

[X,Y] = meshgrid(1:n, 1:m);
X_deci = X(1:5:end, 1:5:end);
Y_deci = Y(1:5:end, 1:5:end);

figure; imshow(im2);
figure; imshow(im1);
hold on;
% draw the velocity vectors
quiver(X_deci, Y_deci, u_deci,v_deci, 'y')

%% Gaussian Laplacian Pyramids