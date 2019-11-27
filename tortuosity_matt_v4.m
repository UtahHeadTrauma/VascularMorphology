function [tavg tmax tmin stDev segtort] =tortuosity_matt_v4( )
%%
%
%INPUT: pre processed image from cleanimage2 function
%(if no input is provided, it will run cleanimage2)
%
%OUTPUT: tortuosity of vessels data file and associated image
%
%written by Kendall McMillan, Matt Byrne
%%
clear
global z gflag      %global variables
warning('off','all')
%%
%create a colored skeleton over the gray image. 
R=.90;  %[default=.9]
G=.0;   %[default=.0]
B=.30;  %[default=.3]

%length of unwanted spurs in skeleton. 
splth=5;    %[default=3]

%%this is the measurement cut off. 
llcutoff=12;  %[default=12]   
RH = [];
gflag=0; %flag for a button
%%
% This Section Measures Tortuosity of Vessel Segments

%load the image
    [totadj, filename]=cleanimage_2();   %run cleanimage function
    RGB=imread(filename); %read in original image
    
    %save folder for data
    sv=['IMdata/' filename(1:end-4)];

    %set up the basic gui structure
    close
    scrsz = get(groot,'ScreenSize');  %This is just setting the screen size of the gui
    scrsz(1)=10;
    scrsz(2)=50;
    scrsz(4)=scrsz(4)-130;
    scrsz(3)=scrsz(3)-20;
    
    %show original and pre-processed images
    mfig=figure('OuterPosition',scrsz,'Position',scrsz,'Color','w','Name',['Retcam image filter: ' filename],'Resize','off');  %open full figure
    subplot(1,2,1)
    imshow(RGB)
    title('original')
    subplot(1,2,2)   
    if mean(mean(totadj(totadj>0.05)))<.3        %catch for images that are very dark to lighten them up
        totadj=totadj+.3;
    end
    totadj=imcomplement(totadj); %take the reverse of the image
    imshow(totadj)      %show adjusted image now reversed
    title('preprocessed, ready to filter')
    
    %% if save original and preprocess comparison
    if exist(sv,'dir')==0
    mkdir(sv)   %create a directory if none exists
    end
    if exist([sv '/orig_' filename(1:end-4) '.jpg'])==0
    saveas(gcf, [sv '/orig_' filename(1:end-4)] ,'jpg')   %save final figure if it does not exist
    end

%%
%vessel selection
subplot(1,1,1)
disp('Using freehand, crop image around vessels you want to measure')
imshow(totadj,'InitialMagnification','fit')
title('Using freehand, crop image around vessels you want to measure')

%%
%use saved data to show previous vessel selections on this image, for
%comparison
if exist([sv '/vessimage_' filename(1:end-4) '.mat'])==0    %check if image exists
vessimage=totadj;   %if not, then create dummy variable for saving
imshow(vessimage,'InitialMagnification','fit');
title('using freehand, crop image around vessels you want to measure')
else
    vessload=load([sv '/vessimage_' filename(1:end-4)],'vessimage'); %if image exists, load it into workspace
    vessimage=vessload.vessimage; 
    imshow(vessimage,'InitialMagnification','fit');
    title('using freehand, crop image around vessels you want to measure')

    %GUI buttons initilization
guibl
guibcont
uiwait
delete(gbl);delete(gbcont);

end
t2=uicontrol('Style','text','String','Hold Left mouse and drag a circle around the vessels to be measured.',...
    'Position',[10,55,800,25],'ForegroundColor','b','FontSize',13);
hF=imfreehand();                                       %crop out what you don't need 
delete(t2);
%%
%masking and save selection image
bmask = hF.createMask(); %create black/white mask
title(' ');
edgemask=bwmorph(bmask,'remove'); %remove intorior pixels                     
vessimage=vessimage+edgemask; %black out circle around mask area
imshow(vessimage,'InitialMagnification','fit');
saveas(gcf, [sv '/vessim_' filename(1:end-4)] ,'jpg')  %save figure
crv=totadj;                                            %creates gray image
crv(~bmask) = 1;

%%
%Show new image
imshow(crv,'InitialMagnification','fit')
title('masked image')
z='y';
disp('error')
%%
%this allows user to "zoom" into selection
disp('crop to zoom (draw box around desired area and double click inside box)')   %crop out white area of image
title('draw box around desired area and double click inside box')
t2=uicontrol('Style','text','String','crop to zoom (draw box around desired area and double click inside box)',...
    'Position',[10,130,500,25],'ForegroundColor','b');
crv=imcrop(crv); %crop image to measurment area
crvo2=crv;
title('use the crop tool to zoom')
delete(t2)
disp('error 2')
% end
crvo2=crv;
title('filtering...')
crv=imcomplement(crv);                                 %make the image negative
imshow(crv,'InitialMagnification','fit')               %see the gray image before thresholding              
title('filtering...')
[a,b]=size(crv);
sdtsz=floor(sqrt(a^2+b^2)./60);                        %finds a value relative to total size
areax=floor(a*b*.008);                                 %compute a precentage of the total area for bwareaopen
areay=floor(a*b*.0045);                                %remove increasingly small areas of noise
areaz=floor(areay/2);                                  %determine the size of areas removed
%%
BW2C= cat(3, ones(size(crv)).*R, ones(size(crv)).*G, ones(size(crv)).*B);  %color matrix
crvo=crv; %reset original image.
%start over point
redo=1;
while redo==1  %allows for reseting of the image and starting over (this is a horrible way to code this)
disp('error 3')
crv=crvo;

 t3=uicontrol('Style','text','String','Skeletonization Threshold:',...
    'Position',[10,25,250,25],'ForegroundColor','b');
            s = sprintf('Changing the slider value starts the skeletonization algorithm with a different greythresh value');
            t3.TooltipString = s;
%filtering buttons
guibsl
adaptdarkb
adaptbrightb
wfiltb
histadjust
resetb
medfiltb
lowfiltb
cutout
uiwait
%get rid of buttons
delete(gbgo);delete(sliderval);delete(addarkb);delete(adbrightb);delete(wfb);delete(histadjustb);delete(gbr);
delete(medfb);delete(lowfb);delete(cutoutb);delete(t3);delete(t2);delete(gbcalc); 
redo=0;

title(' ')
%post filtering and display
BWn=bwareaopen(BWn,10);                            %remove any remaining small objects
subplot(1,2,1)                                     %display final image before crop
imshow(crv,'InitialMagnification','fit')
hold on
h=imshow(BW2C,'InitialMagnification','fit');       %plot color mask
hold off
set(h,'AlphaData',BWn);                            %create color skeleton
subplot(1,2,2)
imshow(BWn,'InitialMagnification','fit')           %plot bw skeleton

%spur removal option
title('Remove spurs?')
t2=uicontrol('Style','text','String','Remove spurs?',...
    'Position',[10,130,200,25],'ForegroundColor','b');
s = sprintf('Spurs are extra pixles that may be attached to the longer lines.');
    t2.TooltipString = s;
%GUI buttons
guibn
guiby
uiwait
delete(gbyes)
delete(gbno)
delete(t2)
q=z;

if q=='y'
    t2=uicontrol('Style','edit','String','3',...
    'Position',[100,100,200,25],'ForegroundColor','b');
    s = sprintf('type spur length here');
    t2.TooltipString = s;
guibgo                                             %go button
uiwait                                             %wait for input
sptype=str2num(get(t2,'String'));
if isnumeric(sptype)==1 %check spur type and assign to spur length
   splth=sptype; 
end
delete(t2)
delete(gbgo)
    BWn=bwmorph(BWn,'spur',splth);                 %removes extra spurs in the final skeleton. 
end

BWn=bwareaopen(BWn,6);                             %delete very small objects. 
%%
subplot(1,2,1)
imshow(crvo,'InitialMagnification','fit')          %show original image
hold on;
h=imshow(BW2C,'InitialMagnification','fit');       %make color mask
hold off
set(h,'AlphaData',BWn);                            %create color skeleton
subplot(1,2,2)
imshow(BWn,'InitialMagnification','fit')           %create bw skeleton

assignin('base', 'BWn', BWn); %assign black and white image to global space
%%
%multiple branches computation of tortuosity

%find branchpoints and remove the pixels around them
bp=bwmorph(BWn,'branchpoints');                    %create matrix of branchpoints
[k,l]=find(bwmorph(BWn,'branchpoints'));           %find branchpoints
D=bwdistgeodesic(BWn,find(bp),'quasi');            %create first geodesic matrix
[x y]=find(D<2 & D>.5);                            %find points around the branchpoints

bp2=zeros(size(bp));                               %create a matrix to preallocate space
c1=1;
[d1 d2]=size(x);                                   %get amount of pixles in x
while c1<d1                   
bp2(x(c1),y(c1))=1;                                %get rid of sqrt(2) values at branch points
c1=c1+1;
end
bp2=im2bw(bp2);
BW2=BWn-bp-bp2;                                    %create skeleton matrix without branchpoints
BW2=im2bw(BW2);

%This finds the geodesic distance of the new image and plots points
ep=bwmorph(BW2,'endpoints');                       %create maxtrix of endpoints
D2=bwdistgeodesic(BW2,find(ep),'quasi');           %create second geodesic matrix
[i,j]=find(bwmorph(BW2,'endpoints'));              %find endpoints
hold; plot(j,i,'yd'); plot(l,k,'r*');              %plot yellow diamonds on endpoints and red asterisks on branchpoints 
hold off

%measure all the seperate branches
[L, num] = bwlabel(BW2);                           %create seperate regions for each connected area (between endpoints)
%%
grnct=1;
c3=1;
segtort=0;
while c3<num+1
STATS = regionprops((L==c3), D2, 'MaxIntensity');  %find the highest value in the selected region
tcdist=STATS.MaxIntensity;      %get max intensity (value) of vessel, this is half of vessel length
if tcdist > (llcutoff/2)                           %this is the measurment cut off.
curvydist(1,c3)=(tcdist*2)+sqrt(2);                %to determine # of pixles in the line, add back the sqrt(2) that was removed earlier. 

[r c]=find(L==c3);      %this gets all the pixel locations for vessel
rc=[r c];
linedist(1,c3)=pdist2(max(rc),min(rc));            %find the straight distance to the endpoints of that line
segtort(c3)=curvydist(c3)/linedist(c3);            %calculate tortuosity, store in matrix
curvydist=curvydist(curvydist>0);

text(mean(c),mean(r),num2str(grnct),'color','green')   %this is broken, grnct needs to be iterated

%this plots the region number near the center of the vessel
text((mean(c)+sqrt(mean(c))),mean(r)+sqrt(mean(r)),num2str(segtort(c3),4),'color','yellow') %plot tortuosity on image next to vessel 
end
c3=c3+1;
end

%%this ends the redo loop
t2=uicontrol('Style','text','String','Does this look okay?  If done press yes, if not press no.',...
    'Position',[10,75,500,25],'ForegroundColor','b','FontSize',13);
%GUI buttons
guibn
guiby
uiwait
delete(gbno)
delete(gbyes)
delete(t2)

if z=='n'
    redo=1;
    gflag=0;
end

end
%if you have a tort value <1, don't include it in further calculation. 
%this can happen if you have a very small segment that you are trying to measure.
segtort=segtort(segtort>=1);        
[tmax location_max] = max(segtort)   %find the max tort value
[tmin location_min] = min(segtort)   %find the min value
tavg=mean(segtort)                   %find average value
stDev=std(segtort)                   %find standard deviation 
file2=filename(1:end-4);             %remove the file type from name

ni=1;                                %initalize indexing variables
fi=1;

%Save the figure if a workspace variable exists
while ni==1
if exist([sv '/tort_' file2 '.mat'])==0  
saveas(gcf, [sv '/fig_' file2] ,'jpg')  
ni=0;
else
    if fi>1
    file2=[file2(1:end-1) num2str(fi)];  %iterate the file number  
    else
    file2=[file2 num2str(fi)];
    end
    fi=fi+1;
end
end

%%
% This Section Measures Diameters of Vessel Segments
close all
%open full figure
mfig=figure('OuterPosition',scrsz,'Position',scrsz,'Color','w','Name',['Retcam image filter: ' filename],'Resize','off');  

c3=1;
while c3<(max(max(L))+1) %iterate through all vessels
imshow(crvo2,'InitialMagnification','fit')            %show original image
hold on;
h=imshow(BW2C,'InitialMagnification','fit');          %make color mask
hold off
set(h,'AlphaData',(L==c3));                           %create color skeleton
hold on
bpcur=bwmorph(L==c3,'endpoints');                     %get endpoints
quartlength(c3)=length(L(L==c3))./4;
geoD=bwdistgeodesic(L==c3,find(bpcur,1),'quasi');     %do geodesict distance calculation for vessel
[cp1,cp2]=find(geoD>max(max(geoD)./2),1);               %this is for placing markers on the image
[qp1,qp2]=find(geoD>max(max(geoD))./4,1);
[qp3,qp4]=find(geoD>max(max(geoD)).*.75,1);
plot(cp2,cp1,'bx','Markersize',15)
plot(qp2,qp1,'bx','Markersize',15)
plot(qp4,qp3,'bx','Markersize',15)

%start skip mini gui
t2=uicontrol('Style','text','String','skip diameter?',...
    'Position',[10,130,200,25],'ForegroundColor','b');
s = sprintf('if diameter measurement is obscured by edge, omit?');
    t2.TooltipString = s;
%GUI buttons
guibn
guiby
uiwait
delete(gbyes)
delete(gbno)
delete(t2)

if z=='n' %skip if
    
%start of diameter selection
xc=0;
     while length(xc) ~= 6                          %point selection for endpoints
      [xc, yc] = getpts;
      
      if length(xc) ~= 6
        xc=[];
        yc=[];
         warning('you did not select 6 points, do it again!')
     end
    end

 hold on 
       plot(xc,yc,'y+','MarkerSize',15)
       plot(xc,yc,'rx','MarkerSize',15)   %plot selected points
 hold off

for qq=1:2:6                                        %get distances between points
    dimss(qq)=pdist2([xc(qq),yc(qq)],[xc(qq+1),yc(qq+1)]);
end

diameters(c3)=mean(dimss(dimss~=0)); %record vessel diameters
else
   diameters(c3)=0; 
end

c3=c3+1; %increment to next vessel number
end
%%
%This Section Measures Angles Between trunks and daughter vessels

%This section is still being perfected, will probably be replaced
%by hand measurements for the thesis. 

imshow(crvo2,'InitialMagnification','fit');
[by,bx]=find(D==0)%this gets endpoints of geodesic vessel matrix
imshow(crvo2,'InitialMagnification','fit');
 hold on
for bb=1:length(by)
    disp('plot 2 points');
   plot(bx(bb),by(bb),'r*'); %plot the endpoints
    xc=[];
    yc=[];
     while length(xc) ~= 2                  %point selection for endpoints
      [xc, yc] = getpts; 
      
      if length(xc) ~= 2
        xc=[];
        yc=[];
         warning('you did not select 2 points, try again!')
     end
     end
    %vertex diameter measurement
    radss(bb)=pdist2([xc(1),yc(1)],[xc(2),yc(2)]);
    
    %plot points corresponding to diameter measurement, from the vertex. 
    tempD=D;
    tempD(tempD>-1)=1;
    tempD(isnan(tempD))=0;
	tempD=logical(tempD);
    tbp=find(bwmorph(tempD,'branchpoints'),length(by));
    tempD=bwdistgeodesic(tempD,tbp(bb),'quasi');
    tempD(tempD>radss(bb))=0;
    tempD(isnan(tempD))=0;
    ep3=bwmorph(tempD,'endpoints');
    [ey,ex]=find(ep3); 
   
    v2 = [ex(1)-bx(bb), ey(1)-by(bb)];
  v3 = [ex(2)-bx(bb), ey(2)-by(bb)];
% Normalize the two vectors:
  u2 = v2 / norm(v2);
  u3 = v3 / norm(v3);
% The angle is then theta1 = acos( dot(u2, u3) ) in radians
 theta1(bb) = acos(dot(u2, u3));
 text((abs(ex(1)-ex(2))./2+min([ex(1),ex(2)])),(abs(ey(1)-ey(2))./2+min([ey(1),ey(2)])),num2str(theta1(bb)),'FontSize',5,'color','red')
 
  v2 = [ex(3)-bx(bb), ey(3)-by(bb)];
  v3 = [ex(2)-bx(bb), ey(2)-by(bb)];
% Normalize the two vectors:
  u2 = v2 / norm(v2);
  u3 = v3 / norm(v3);
% The angle is then theta1 = acos( dot(u2, u3) ) in radians
 theta2(bb) = acos(dot(u2, u3));
  text((abs(ex(3)-ex(2))./2+min([ex(3),ex(2)])),(abs(ey(3)-ey(2))./2+min([ey(3),ey(2)])),num2str(theta2(bb)),'FontSize',5,'color','red')

  v2 = [ex(1)-bx(bb), ey(1)-by(bb)];
  v3 = [ex(3)-bx(bb), ey(3)-by(bb)];
% Normalize the two vectors:
  u2 = v2 / norm(v2);
  u3 = v3 / norm(v3);
% The angle is then theta1 = acos( dot(u2, u3) ) in radians
 theta3(bb) = acos(dot(u2, u3));
    text((abs(ex(1)-ex(3))./2+min([ex(1),ex(3)])),(abs(ey(1)-ey(3))./2+min([ey(1),ey(3)])),num2str(theta3(bb)),'FontSize',5,'color','red')

end
thetas=[theta1',theta2',theta3'];
ni=1;  %initalize indexing variables
fi=1;

while ni==1
if exist([sv '/tort_' file2 '.mat'])==0  
saveas(gcf, [sv '/fig_angs_' file2] ,'jpg')  %save final figure
ni=0;
else
    if fi>1
    file2=[file2(1:end-1) num2str(fi)];    %increment the filename
    else
    file2=[file2 num2str(fi)];
    end
    fi=fi+1;
end
end

%%
%this records the values in a .mat file
file2=filename(1:end-4);                     %remove the file type from name
if exist(sv,'dir')==0
    mkdir(sv)                                %create save directory if none exists
end

ni=1;                                        %initalize indexing variables
fi=1;
save([sv '/vessimage_' file2 '.mat'],'vessimage') 

%save vessel data and image in unique file location
while ni==1
if exist([sv '/tort_' file2 '.mat'])==0  
save([sv '/tort_' file2 '.mat'],'tmax','location_max','tmin','location_min','tavg','stDev','segtort','RH','vessimage','curvydist','L','crvo2','BW2C','diameters','bp','D','thetas') %save the variables
ni=0;
else
    if fi>1
    file2=[file2(1:end-1) num2str(fi)];    
    else
    file2=[file2 num2str(fi)];
    end
    fi=fi+1;
end
end
disp(filename)
D = dir(['IMdata']);
Num = length(D(([D.isdir])))-2;
close('all')    %end of program, return nothing

  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      button functions                         %%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function guiby     %the gui yes button
% global  gbyes z

gbyes=uicontrol('Style','pushbutton','String','yes?',...
    'Position',[10,50,100,25],'Callback',@yn);
    function yn(hObject,eventdata,handels) 
        
       z='y';  
        uiresume
    end
end
%%
function guibn        %the gui no button
% global  gbno z

gbno=uicontrol('Style','pushbutton','String','no?',...
    'Position',[10,25,100,25],'Callback',@ny);
    function ny(hObject,eventdata,handels)
        z='n';        
        uiresume
    end
end
%%
function guibl     %the gui toggle lines button
% global  gbyes z

gbl=uicontrol('Style','togglebutton','String','Toggle Lines',...
    'Position',[10,50,100,25],'Callback',@l);
s = sprintf('Press to toggle showing the previous vessel selection lines');
    gbl.TooltipString = s;
    function l(hObject,eventdata,handels) 
        button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
       imshow(totadj,'InitialMagnification','fit'); %this toggles lines for if you have done multiple calculations on the same image
       set(gbl,'ForegroundColor','red')
       set(gbl,'String','Lines Hidden')      

else
    imshow(vessimage,'InitialMagnification','fit');
    set(gbl,'ForegroundColor','black')
    set(gbl,'String','Lines Visible')  
    
end              
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adaptive darken button
function adaptdarkb
addarkb=uicontrol('Style','pushbutton','String','Adaptive darken',...
    'Position',[10,75,100,25],'Callback',@addarkbutton);
s = sprintf('Press to darken image (darker areas will be darkened more than lighter areas)');
    addarkb.TooltipString = s;
function addarkbutton(hObject,eventdata,handels)
   crv=crv.^1.3;
    crv=crv.*255; %runns the addaptive darken algorithm 
im2=(crv./10);
crv=crv+im2;
crv=crv./255;
crv(crv<0)=0;
crv(crv>1)=1;
imshow(crv,'InitialMagnification','fit')
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  wiener filter button
    function wfiltb
wfb=uicontrol('Style','pushbutton','String','wiener filter',...
    'Position',[10,150,100,25],'Callback',@wbutton);
 s = sprintf('Press to apply a wiener filter, (adaptive noise-removal)');
    wfb.TooltipString = s;
function wbutton(hObject,eventdata,handels)
crv=wiener2(crv,[sdtsz,sdtsz]); %wiener filter, neighborhood is 1/60th the diagonal dist of image.
imshow(crv,'InitialMagnification','fit')
end
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adaptive brighten button
function adaptbrightb
adbrightb=uicontrol('Style','pushbutton','String','Adaptive brighten',...
    'Position',[10,125,100,25],'Callback',@adbrightbutton);
s = sprintf('Press to apply adaptive brightening, (brighten, subtract mean, and enhance contrast)');
    adbrightb.TooltipString = s;
function adbrightbutton(hObject,eventdata,handels)
   crv=crv.*255; 
crv=(crv.^1.1);
crv=crv-(mean(mean(crv))/4);     %runs the addaptive brightnen algorithm
crv=crv./255;
lims=stretchlim(crv);
if lims(2)<.9
    lims(2)=lims(2)+.05;
end
crv=imadjust(crv,lims);
imshow(crv,'InitialMagnification','fit')
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   median filter button
    function medfiltb
medfb=uicontrol('Style','pushbutton','String','median filter',...
    'Position',[10,100,100,25],'Callback',@medbutton);
s = sprintf('Press to apply a median filter, filter matrix size is 1/60th the diagonal distance of the current image');
    medfb.TooltipString = s;
function medbutton(hObject,eventdata,handels)
crv=medfilt2(crv,[sdtsz,sdtsz]);  %median filter
imshow(crv,'InitialMagnification','fit')
end
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   edging lowpass filter button
    function lowfiltb
lowfb=uicontrol('Style','pushbutton','String','edging lowpass',...
    'Position',[10,200,100,25],'Callback',@lowbutton);
s = sprintf('Applies a filter that is a combo lowpass and edge sharpening');
    lowfb.TooltipString = s;
function lowbutton(hObject,eventdata,handels)
N=128;
wp=.2;
b=fir1(N,wp);                    %lowpass filter matrix
lp=ftrans2(b);
crv=imfilter(crv,lp);            %low pass filter
crv=abs(crv);
imshow(crv,'InitialMagnification','fit')

end
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  reset buttion function
    function resetb
gbr=uicontrol('Style','pushbutton','String','Reset',...
    'Position',[10,50,100,25],'Callback',@rbutton);
s = sprintf('Press to reset the image');
    gbr.TooltipString = s;
function rbutton(hObject,eventdata,handels)
    crv=crvo;   %set image veriable back to original
    
    imshow(crv,'InitialMagnification','fit');    
end
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   histadjust button
    function histadjust
histadjustb=uicontrol('Style','pushbutton','String','normalize histogram',...
    'Position',[10,175,100,25],'Callback',@histadjustbutton);
s = sprintf('Press to apply adaptive histogram equilization using an exponential curve');
    histadjustb.TooltipString = s;
function histadjustbutton(hObject,eventdata,handels)
crv=adapthisteq(crv,'Distribution','exponential'); %histogram equilization
imshow(crv,'InitialMagnification','fit')
end
    end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   the cutout button
    function cutout
cutoutb=uicontrol('Style','pushbutton','String','Cutout',...
    'Position',[10,225,100,25],'Callback',@cutoutf);
s = sprintf('Press to crop out areas of noise');
    cutoutb.TooltipString = s;
function cutoutf(hObject,eventdata,handels)
    t4=uicontrol('Style','text','String','draw a circle around what you want to blackout',...
    'Position',[300,30,300,25],'ForegroundColor','b');
imshow(crv)
title('circle an area to blackout')
%crop out what you don't want with imfreehand
hF2=imfreehand();                
bmask2 = hF2.createMask();   %algorithm for blacking out sections of image
crv(bmask2)=0;
imshow(crv,'InitialMagnification','fit')
delete(t4)
uiwait
end
    end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  The go button
function guibgo    %the gui go button
gbgo=uicontrol('Style','pushbutton','String','done',...
    'Position',[100,50,100,25],'Callback',@ny);
s = sprintf('Press when finished');
    gbgo.TooltipString = s;
    function ny(hObject,eventdata,handels)
       
        uiresume
    end
end
%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  The continue button
function guibcont    %the gui go button
gbcont=uicontrol('Style','pushbutton','String','Continue',...
    'Position',[100,50,100,25],'Callback',@ny);
s = sprintf('Press to continue running program');
    gbcont.TooltipString = s;
    function ny(hObject,eventdata,handels)
        
        uiresume
    end
end
%%%%%%%%%%%%%%%%%%%%%%
%%                                 
%slider
    function guibsl
sliderval=uicontrol('Style', 'slider','String','Skeleton Threshold',...
        'Min',.001,'Max',1,'Value',.4,...
        'Position', [10 10 300 20],...
        'Callback', {@thrsl});  
    s = sprintf('slide to change threshold value, this also runs the skeletonization code');
    sliderval.TooltipString = s;
        function thrsl(hObject,eventdata,handels)
            thresh=get(sliderval,'Value');
            t2=uicontrol('Style','text','String',['greythresh: ' num2str(thresh)],...
    'Position',[350,10,100,25],'ForegroundColor','b');
          %every time you move the slider the skeleton is recomputed
         BWn=im2bw(crv,thresh);       %change to a BW image using the threshold input
         BWn=bwmorph(BWn,'close');
         BWn=bwareaopen(BWn,areax);   %remove noise
         BWn=bwmorph(BWn,'thin',2);   %thin areas
         BWn=bwmorph(BWn,'spur',1);   %remove spurs
         BWn=bwareaopen(BWn,areay);   %remove more noise spots
      
         BWn=bwmorph(BWn,'thin',5);   %more thinning and spur removing
         BWn=bwmorph(BWn,'spur',1);
         BWn=bwareaopen(BWn,(areaz)); %remove whats left of the noise spots
         BWn=bwmorph(BWn,'thin',4);   %more thinning and skeletonize
         BWn=bwmorph(BWn,'skel');
         BWn=bwmorph(BWn,'thin',Inf); %thin the lines to only one pixle width
        subplot(1,2,1)
        imshow(crv,'InitialMagnification','fit')      %plot gray image
        hold on
        h=imshow(BW2C,'InitialMagnification','fit');  %plot color mask
        hold off                                      %plots the skeleton on top of the gray in color.
        set(h,'AlphaData',BWn);                       %create color skeleton
        subplot(1,2,2)                                %plot bw skeleton separately
        imshow(BWn,'InitialMagnification','fit')
        title('filtering...')
        if gflag==0
           gflag=1;
           guibgo
           guibcalc
        end
        end
  
    end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  The calculate preview button
function guibcalc    %the gui calc button

gbcalc=uicontrol('Style','pushbutton','String','Calculate',...
    'Position',[10,250,100,25],'Callback',@calc);
s = sprintf('Press to show tortuosity calculation preview');
    gbcalc.TooltipString = s;
    function calc(hObject,eventdata,handels)

        BWn=bwareaopen(BWn,10);         %remove any remaining small objects
        BWn=bwmorph(BWn,'spur',splth);  %gets rid of extra spurs in the final skeleton. 

        subplot(1,2,1)
imshow(crvo,'InitialMagnification','fit')      %show original image
hold on;
h=imshow(BW2C,'InitialMagnification','fit');   %make color mask
hold off
set(h,'AlphaData',BWn);                        %create color skeleton
subplot(1,2,2)
imshow(BWn,'InitialMagnification','fit')       %create bw skeleton
%
%find branchpoints and remove the pixels around them
bp=bwmorph(BWn,'branchpoints');                %create matrix of branchpoints
[k,l]=find(bwmorph(BWn,'branchpoints'));       %find branchpoints
D=bwdistgeodesic(BWn,find(bp),'quasi');        %create first geodesic matrix
[x y]=find(D<2 & D>.5);                        %find points around the branchpoints

bp2=zeros(size(bp));                           %create a matrix to prealocate space
c1=1;
[d1 d2]=size(x);                               %get amount of pixles in x
while c1<d1                   
bp2(x(c1),y(c1))=1;                            %get rid of sqrt(2) values at branch points
c1=c1+1;
end
bp2=im2bw(bp2);
BW2=BWn-bp-bp2;                                %create skeleton matrix without branchpoints
BW2=im2bw(BW2);

%This finds the geodesic distance of the new image and plots some points
ep=bwmorph(BW2,'endpoints');               %create maxtrix of endpoints
D2=bwdistgeodesic(BW2,find(ep),'quasi');   %create second geodesic matrix
[i,j]=find(bwmorph(BW2,'endpoints'));      %find endpoints
hold; plot(j,i,'yd'); plot(l,k,'r*');      %plot yellow diamonds on endpoints and red asterisks on branchpoints 
hold off

%measure all the seperate branches
[L, num] = bwlabel(BW2);                    %create seperate "regions" for each connected area (between endpoints)
%
grnct=1;
c3=1;
segtort=0;
while c3<num+1
STATS = regionprops((L==c3), D2, 'MaxIntensity');   %find the highest value in the selected region
tcdist=STATS.MaxIntensity;
if tcdist > (llcutoff/2)                            %this is the measurement cut off. 
curvydist(1,c3)=(tcdist*2)+sqrt(2);                 %get the # of pixles in the line, we add back the sqrt(2) that was removed earlier. 

[r c]=find(L==c3);
rc=[r c];
linedist(1,c3)=pdist2(max(rc),min(rc));             %find the straight distance to the endpoints of that line
segtort(c3)=curvydist(c3)/linedist(c3);             %calculate tortuosity, store in matrix


text(mean(c),mean(r),num2str(grnct),'color','green')   %this plots the region number near the center of the vessel
grnct=grnct+1;
text((mean(c)+sqrt(mean(c))),mean(r)+sqrt(mean(r)),num2str(segtort(c3),4),'color','yellow')   %this plots the region number near the center of the vessel

end
c3=c3+1;
end
       
    end
end
%%%%%%%%%%%%%%%%%%%%%%

end %end main function

