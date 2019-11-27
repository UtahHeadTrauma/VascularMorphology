function [totadj, filename] = cleanimage_2(RGB)
%%
%Input:workspace image variable or file name as string
%string example: 'Image.jpg'
%
%SEARCHES FOR .TIF OR .JPG IMAGES AS INPUT IF YOU DON'T SPECIFY AN INPUT
%
%Updated to run with the Tortuosity_2 function
%
%takes an image file and cleans it up.  It will take the image and convert it
%into several different color types and then take one layer of each and
%average then all together to get a gray image that has more contrast than
%the original.
%
%Outputs: the new image and filename.
%
%writen by Shaun Evans/Kendall McMillan, University of Utah
%%

if nargin==0  %if no input into function
    %%
    D=dir('IMdata');
    D=D(3:end-1);    %get IMdata directory
    nmax=0;

    %%
    filename=uigetfile({'*.tif';'*.png';'*.jpg'},['select the image file'])  %this will let you select a file
    RGB=imread(filename);       %read file
    
elseif isa(RGB,'char')==1 %if char input
        RGB=imread(RGB);
end
%RGB = imread('5036736.tif');  %imports test image

%%
if size(RGB,3)==1
   RGB=cat(3,RGB,RGB,RGB); 

end

%%

ycbcr = rgb2ycbcr(RGB);            %changes color space to ycbcr
hsv = rgb2hsv(RGB);                %changes color space to hsv
dRGB = im2double(RGB);             %changes image to double  
dycbcr = im2double(ycbcr);         %changes ycbcr to double


cform = makecform('srgb2lab', 'AdaptedWhitePoint', whitepoint('D65'));   %creates color transformation structure 
lab = applycform(RGB,cform);                                             %applies color transformation structure 
dlab = im2double(lab);                                                   %changes to double

hsvcomp=(imcomplement(hsv(:,:,2)));
hsvcomp(hsvcomp==1)=0;
ringcomp=medfilt2(hsvcomp,[12,12]);

total1 = (dRGB(:,:,2)+dycbcr(:,:,1)+hsv(:,:,3)+dlab(:,:,1))./4;          %combines select layers of the images and averages them
total2 = (dRGB(:,:,2)+dycbcr(:,:,1)+dRGB(:,:,3)+dlab(:,:,1)+ringcomp)./5;%combines select layers of the images and averages them

total=(total2.^2+total1.^2).^.5;

totadj = adapthisteq(total,'ClipLimit',0.005);    %adjust the histogram
adj2=ringcomp.*totadj;
adj2=medfilt2(adj2,[2,2]);
totadj=wiener2(((adj2.^2+totadj.^2).^.5)); %final adjusted image

end
