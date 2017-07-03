%% INITIALISATION
clear all;
close all;
clc;
%% INPUT SEGMENT
% The variable 'structimg' stores the input image (in .mat format) in a
% structure. As the structure can't be directly converted to a matrix
% without knowing the field names, its converted to a cell and then into a
% matrix.
strctimg=uiimport();
cellimg=struct2cell(strctimg);
matimg=cell2mat(cellimg);
siz=size(matimg);
% The variable 'start_spectra' stores the starting wevelength of the
% hyperspestral image, similarly 'end_spectra' stores the end wavelength
% of the image. The wavelength values must be given in nanometer(nm)
% scale. The user has to give a pixel value for which the spectrum is to
% be displayed.
N = input('Please enter the Julian day of the year on which picture was
taken: ');
time = input('Please enter the time in 24Hrs format [HH, MM]: ');
lati = input('Pleae enter the latitude of the site: ');
pix = input('Enter the pixel values for spectrum display: ');
nme = input('Please enter the desired name of your output RGB image: ');
% The program takes the middle value of each band (red green and blue)
% from the visivle spectrum. For future reference in visible spectrum:
% RED = 400–484 THz or 620–750 nm
% Green = 526–606 THz or 495–570 nm
% Blue = 606–668 THz or 450–495 nm
% Here we have taken:
% Red = 700nm
% Green = 546.1nm
% Blue = 435.8nm
% From the these data the band number containing the R G & B Spectrum is
% calculated and stored for RGB image display
% start_spectra= 350.505 nm;
% end_spectra = 2582.13 nm;
ext_pic = '.jpg';
ext1_mat = '_modified';
ext2_mat = '.mat';
nme_v = strcat(nme,'_visible');
nme_ve = strcat(nme,'_vegitation');
nme_S = strcat(nme,'_SWIR');
outp_pic_v = strcat(nme_v,ext_pic);
outp_pic_ve = strcat(nme_ve,ext_pic);
outp_pic_S = strcat(nme_S,ext_pic);
outp_mat = strcat(nme,ext1_mat,ext2_mat);
t = time(1) + (time(2))/60;
% red_index=39;
% green_index=24;
% blue_index=17;
red_index_visible=29;
green_index_visible=23;
blue_index_visible=16;
red_index_vegi=50;
green_index_vegi=23;
blue_index_vegi=16;
red_index_SWIR=204;
green_index_SWIR=150;
blue_index_SWIR=93;
%% RADIOMETRIC CORRECTIONS
% The raw data given has pixel values those are digital numbers recorded
% by sensors. For an image to be displayed, we need the reflactace matrix
% at the surface of earth. For this we need to followseveral steps like:
% 1. Convertion of DNs to Radiance Matrix(RM)
% 2. Conversion od RMs to Exo-atmospheric Reflectance Matrix
% 3. Application of weather correction to get Earth Reflectance matrix
%% DN TO RADIANCE MATRIX
% According to Hyperion Sensor dataset in USGS website for EO-1 satellite
% (http://eo1.usgs.gov/faq/question?id=19) it is clear that the scaling
% factor for the VNIR band (upto band 70) is 40 and for SWIR (band >70)
% it is 80.
for i=1:70
	img_radiance(:,:,i) = 0.025*(matimg(:,:,i));
end;
for i=71:siz(3)
	img_radiance(:,:,i) = 0.0125*(matimg(:,:,i));
end;
%% SR TO EXO-ATMOSPHERIC REFLECTANCE
load('JD_dist');
d = JD_dist(N);
H = 15*(t-12);
load('delta.mat');
del = delta(N);
% SZ = cosd(sz) = sind(lati)*sind(del) + cosd(lati)*cosd(del)*cosd(H);
SZ = sind(lati)*sind(del) + cosd(lati)*cosd(del)*cosd(H);
load('ESUN_pidiv.mat');
for i = 1:siz(3)
	ex_atm_ref(:,:,i) = (ESUN(i)*d*d*img_radiance(:,:,i))/SZ;
end;
%% DARK PIXEL SUBTRACTION & NOISE FILTER
% The following loop itterates over each band and subtracts the minimum
% pixel valus of the corresponding band from all piexls. This method
% removes the scattering due to the sunlight and gets rid of the negative
% valued pixels.
for i=1:(siz(3))
	ma(i)=min(min(ex_atm_ref(:,:,i)));
	img_filtered(:,:,i)=ex_atm_ref(:,:,i)-ma(i);
end
clear matimg start_spectra end_spectra band_range band_width;
%% NORMALISATION
% In normalisation we map the values of the original image to [0 1]. This
% is implemented both on raw imge and filtered image.
for i = 1:siz(3)
	matimg_filtered_norm(:,:,i) =
	img_filtered(:,:,i)/(max(max(img_filtered(:,:,i))));
end;
%% RGB MATRIX FORMATION
% The following commands stores the R G & B bands in a 3D matrix which
% can be further dispalyed. Here also canny edge detector is used for
% segmenting the image and display outlines. This is donr both for
% filtered and original images. Mind that the display is only done on the
% normalised image.
imgRf_v = matimg_filtered_norm(:,:,red_index_visible);
imgRf_v = imadjust(imgRf_v);
imgGf_v = matimg_filtered_norm(:,:,green_index_visible);
imgGf_v = imadjust(imgGf_v);
imgBf_v = matimg_filtered_norm(:,:,blue_index_visible);
imgBf_v = imadjust(imgBf_v);
img_filtered_view_v(:,:,1) = imgRf_v;
img_filtered_view_v(:,:,2) = imgGf_v;
img_filtered_view_v(:,:,3) = imgBf_v;
imgRf_ve = matimg_filtered_norm(:,:,red_index_vegi);
imgRf_ve = imadjust(imgRf_ve);
imgGf_ve = matimg_filtered_norm(:,:,green_index_vegi);
imgGf_ve = imadjust(imgGf_ve);
imgBf_ve = matimg_filtered_norm(:,:,blue_index_vegi);
imgBf_ve = imadjust(imgBf_ve);
img_filtered_view_vegi(:,:,1) = imgRf_ve;
img_filtered_view_vegi(:,:,2) = imgGf_ve;
img_filtered_view_vegi(:,:,3) = imgBf_ve;
imgRf_S = matimg_filtered_norm(:,:,red_index_SWIR);
imgRf_S = imadjust(imgRf_S);
imgGf_S = matimg_filtered_norm(:,:,green_index_SWIR);
imgGf_S = imadjust(imgGf_S);
imgBf_S = matimg_filtered_norm(:,:,blue_index_SWIR);
imgBf_S = imadjust(imgBf_S);
img_filtered_view_S(:,:,1) = imgRf_S;
img_filtered_view_S(:,:,2) = imgGf_S;
img_filtered_view_S(:,:,3) = imgBf_S;
%% SPECTRUM BUILDING
% The following loop is used to find the spectrum vector of each of the
% possible image in the program. It is applied on raw data (both filtered
% and not filtered) aswell as the normalised data. Then the plots are
% compred to see the changes agter filtering. We can extract features
% from the plot for further classification.
for i = 1:siz(3)
	spectrum_original(i) = ex_atm_ref(pix(1),pix(2),i);
end;
%% DISPLAY SECTION
% The following commands display the essential images and plots
% sequeintially.
hFig4 = figure('Toolbar','none', 'Menubar','none');
him4= imshow(img_filtered_view_v);
hSP4 = imscrollpanel(hFig4,him4);
hFig5 = figure('Toolbar','none', 'Menubar','none');
him5= imshow(img_filtered_view_vegi);
hSP5 = imscrollpanel(hFig5,him5);
hFig6 = figure('Toolbar','none', 'Menubar','none');
him6= imshow(img_filtered_view_S);
hSP6 = imscrollpanel(hFig6,him6);
figure;
plot(spectrum_original,'blue');
grid on;
xlabel('Band');
ylabel('Pixel Value');
title('Original Image Spectrum')
imwrite(img_filtered_view_v,outp_pic_v);
imwrite(img_filtered_view_vegi,outp_pic_ve);
imwrite(img_filtered_view_S,outp_pic_S);
save(outp_mat,'matimg_filtered_norm');
%% END OF CODE