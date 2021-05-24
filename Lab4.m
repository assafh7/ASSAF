%% section 1 - question 2 - fourier based

close all; clear all;
 
image = imread("C:\Users\Assaf Hasky\Desktop\Image Processing Labs\Lab4DIP\Hough Images\lines.tif");
[m, n] = size(image);

BSF = ones(m, n);
BSF(m / 2 - 1 : m / 2 + 1,  :) = 0;
F = fft2(image);
F_shift = fftshift(F);
figure; imshow(log(abs(F_shift)),[]);title("FFT of the lines Image")


diagonalLines_F_domain = F_shift .* BSF;
figure; imshow(log(abs(diagonalLines_F_domain)),[]);title("FFT of the filtered lines Image ")
diagonalLines = ifft2(ifftshift(diagonalLines_F_domain));
diagonalLines = ~imbinarize(abs(diagonalLines));

verticalLines = ~(abs(image - diagonalLines));

figure;
subplot(1, 3, 1);
imshow(image);
title("original");

subplot(1, 3, 2);
imshow(verticalLines);
title("vertical lines");

subplot(1, 3, 3);
imshow(diagonalLines, []);
title("diagonal lines");

% figure;
% imshow(log(abs(fftshift(F))), []);


%% section 2 - question 2 - hough based

close all; clear all;

image = imread("C:\Users\Assaf Hasky\Desktop\Image Processing Labs\Lab4DIP\Hough Images\lines.tif");
image_neg = ~image;
% image_edge = edge(image, 'canny');
[m, n] = size(image);

diagonalLinesImage = image;
diagonalLinesImage = image;

figure;
imshow(image);
title("original");

% start_angle = 10;
% end_angle = 80;
% theta_resolution = 0.5;

%[H,T,R] = hough(b, 'Theta', start_angle:theta_resolution:end_angle);

[H, theta, rho] = hough(image_neg);

figure;
imshow(H,[], 'XData', theta, 'YData', rho,...
            'InitialMagnification', 'fit');
xlabel('\theta'), ylabel('\rho');
axis on;
axis normal;
hold on;
peaks  = houghpeaks(H, 7, 'threshold', ceil(0.3 * max(H(:))));
x = theta(peaks(:, 2)); 
y = rho(peaks(:, 1));
plot(x, y, 's', 'color', 'white');
title("HoughPeaks");

lines = houghlines(image_neg, theta, rho, peaks, "FillGap", 40, "MinLength", 100);

figure(3);
imshow(image);
title("hough diagonal lines");
hold on;
figure(4);
imshow(image);
title("hough vertical lines");
hold on;
for k = 1:size(lines, 2)
   if lines(k).theta ~= 0
        xy = [lines(k).point1; lines(k).point2];
        figure(3);
        plot(xy(:,1), xy(:,2), 'LineWidth', 0.75, 'Color', 'red');
   else
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,1), xy(:,2), 'LineWidth', 0.75, 'Color', 'green');
   end
end

for k = 1:size(lines, 2)
    if lines(k).theta == 0
        col = lines(k).point1(1);
        diagonalLinesImage(:, col - 2 : col + 2) = 1;
    end
end

verticalLinesImage = double(image) - double(diagonalLinesImage);

figure;
imshow(diagonalLinesImage);
title("only diagonal lines");
figure;
imshow(verticalLinesImage, []);
title("only vertical lines");
%% section 3 - question 3a

close all; clear all;

image = imread("C:\Users\Assaf Hasky\Desktop\Image Processing Labs\Lab4DIP\Hough Images\jeep2.jpg");
image_gray = rgb2gray(image);
image_edge = edge(image_gray, 'canny');

figure;
subplot(2, 2, 1);
imshow(image_gray);
title("Original");

subplot(2, 2, 2);
imshow(image_edge);
title("Edged"); 

[H, theta, rho] = hough(image_edge);

subplot(2, 2, 3);
imshow(H,[], 'XData', theta, 'YData', rho,...
            'InitialMagnification', 'fit');
xlabel('\theta'), ylabel('\rho');
axis on;
axis normal;
hold on;
peaks  = houghpeaks(H, 10, 'threshold', ceil(0.5 * max(H(:))));
x = theta(peaks(:, 2)); 
y = rho(peaks(:, 1));
plot(x, y, 's', 'color', 'white');
title("HoughPeaks");

lines = houghlines(image_edge, theta, rho, peaks, "FillGap", 40, "MinLength", 100);

subplot(2, 2, 4);
imshow(image_gray);
hold on;
max_len = 0;
for k = 1:size(lines, 2)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1), xy(:,2), 'LineWidth', 0.75, 'Color', 'green');
   len = norm(lines(k).point1 - lines(k).point2);
   if (len > max_len)
      max_len = len;
      xy_long = xy;
      theta = lines(k).theta;
   end
end
plot(xy_long(:,1), xy_long(:,2), 'LineWidth', 2, 'Color', 'cyan');
title("The jeep angle (cyan) = " + theta);

%% section 4 - question 3b
close all; clear all;

image = imread("C:\Users\Assaf Hasky\Desktop\Image Processing Labs\Lab4DIP\Hough Images\jeep2.jpg");
image_gray = rgb2gray(image);
imshow(image_gray);
impixelinfo();

%check with impixelinfo to assess the radius of the wheels

[centers, radii] = imfindcircles(image_gray, [10 30]);

centersStrong2 = centers(1:2,:); 
radiiStrong2 = radii(1:2);

viscircles(centersStrong2, radiiStrong2, 'EdgeColor', 'b'); hold on;

center1.row = centers(1, 2);
center1.col = centers(1, 1);
center2.row = centers(2, 2);
center2.col = centers(2, 1);

coefficients = polyfit([center1.col center2.col], [center1.row center2.row], 1);
m = coefficients (1);
b = coefficients (2);
radians = atan(m);
degrees = rad2deg(radians);
negativeDegrees = -(90 - degrees);

plot([center1.col center2.col], [center1.row center2.row],...
    'LineWidth', 2, 'Color', 'g');
title("The jeep angle: " + degrees + " / " + negativeDegrees);



