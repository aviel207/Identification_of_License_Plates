clear all;close all;clc;

%% Selecting an image - בחירת תמונה

[fn,pn]=uigetfile({'*.*';'*.jpg';'*.tif';'*.bmp'}, 'Select an image');

%% Image resizing and histogram equalization - שינוי גודל תמונה ושיויון היסטוגרמה

%horizontally concatenates text in arrays - משרשרת אופקית טקסט במערכים
Image = imread(strcat(pn,fn));
sizeImage = size(Image);
% resizes the image so that it has the specified number of rows and columns - משנה את גודל התמונה כך שתכלול את המספר שצוין של שורות ועמודות
if (sizeImage(2)*sizeImage(1)<800*500) 
    Image=imresize(Image,[6]);
end
% Performing histogram equalization - מבצע שיויון היסטגרמה
Image = imadjust(Image, stretchlim(Image),[0 1]); 
imshow(Image)

%% Identifying and filtering objects - זיהוי וסינון עצמים

%Identify the yellow objects in the picture - מזהה את העצמים הצהובים בתמונה 
J = (Image(:,:,1)>130)&(Image(:,:,2)>80)&(Image(:,:,3)<85); 
% Filters components that are less than 1000 pixels from the image,
% (roughly the PIX number of the license plate) - מסנן רכיבים שהם פחות מ-1000 פיקסלים מהתמונה
J = bwareaopen(J,1000);
imshow(J)
%% Finding the right object for a license plate - מציאת העצם המתאים ללוחית רישוי

% Identifies the white objects in a Black and White image - זיהוי העצמים הלבנים בתמונה מסוג שחור-לבן
bwI = bwconncomp(J);
% The number of objects detected in the image -מספר העצמים שזוהו בתמונה
numberObj = bwI.NumObjects;
% Gets the properties of the image j - j מקבל את מאפייני התמונה 
Properties = regionprops(J,'All','Image');

for i=1:numberObj
    Currentfeature = Properties(i);
    %Cut out the object from the images - Image,J חותך את העצמים מהתמונות
    License_plate = imcrop(Image, [Currentfeature.BoundingBox(1) Currentfeature.BoundingBox(2) Currentfeature.BoundingBox(3) Currentfeature.BoundingBox(4)]);  
    shapePlate = imcrop(J, [Currentfeature.BoundingBox(1) Currentfeature.BoundingBox(2) Currentfeature.BoundingBox(3) Currentfeature.BoundingBox(4)]);
    %Rotates the Image,J images using the interpolation method "bilunear" - "bilunear" מסובב את התמונות לפי שיטת האינפולציה
    License_plate = imrotate(License_plate,-Currentfeature.Orientation,'bilinear'); 
    shapePlate = imrotate(shapePlate,-Currentfeature.Orientation,'bilinear');
    % Gets the properties of the image shapePlate - shapePlate מקבל את מאפייני התמונה 
    feature = regionprops(shapePlate,'All','Image'); 
    % Gets the values of BOUNDING BOX - BOUNDING BOX -מקבל את הערכים של ה 
    Bounding = feature.BoundingBox;
    %Cut out the object from the image - plate חותך את העצמים מהתמונה 
    License_plate = imcrop(License_plate, [Bounding(1) Bounding(2) Bounding(3) Bounding(4)]);
    %The Height and Width of plate and the ratio between them - הגובה והרוחב של התמונה והיחס ביניהם
    H=Bounding(3);
    Y=Bounding(4);
    ratio = H/Y;
    % Checks if they match the ratio of the license plate% Checks if they - בדיקת תאימות ללוחית הרישוי
    if(ratio>3.9 && ratio<5.1) 
        break;
    else
        % check other object - נבדוק את העצם הבאה 
        i=i+1;
ס    end
end
imshow(License_plate)
%% Adjusting and filtering the image - התאמת וסינון התמונה
%RGB to grayscale conversion - המרה לגווני אפור
License_plate = rgb2gray(License_plate);
%Performing histogram equalization - ביצוע השוואת היסטוגרמה
License_plate = imadjust(License_plate);
% Convert to image BW above 80 grayscale - המר לתמונה BW מעל 80 גווני אפור
License_plate = (License_plate<80); 
%compare to size - השוואה לגודל 
License_plate = imresize(License_plate, [135,590]);
%Filters components that are less than 500 pixels from the image -  מסנן רכיבים שהם פחות מ-500 פיקסלים מהתמונה
License_plate = bwareaopen(License_plate,500); 
imshow(License_plate)
%% Finding the object that matches the digit of a license plate - מציאת העצם המתאים לספרה של לוחית רישוי
%Create a vector to receive the digits
PlateNumber = [];
%Gets the properties of the image shapePlate - shapePlate מקבל את מאפייני התמונה 
Nums = regionprops(License_plate,'All','Image');
%Identifying the white objects (the numbers)-  (זיהוי העצמים הלבנים(המספרים 
bwI = bwconncomp(License_plate);
% The number of digits detected in the image - מספר הספרות שזוהו בתמונה
numberObj = bwI.NumObjects;

Counter = 1;
for i = 1:numberObj
    CurrentNum = Nums(i);
    %Cut out the digits from the images - חותך את העצמים מהתמונה
    Num = imcrop(License_plate, [CurrentNum.BoundingBox(1) CurrentNum.BoundingBox(2) CurrentNum.BoundingBox(3) CurrentNum.BoundingBox(4)]);
    %The Height and Width of plate and the ratio between them - הגובה והרוחב של העצם והיחס ביניהם
    H = CurrentNum.BoundingBox(3);
    W = CurrentNum.BoundingBox(4); 
    ratio = W/H; 
    %Checking if the ratio matches the size of the license plate digits - בדיקה אם היחס מתאים לגודל של ספרות לוחית הרישוי
    if(ratio>1.2 && ratio<4) 
        % VECTOR SIZE 10
        Similarity = [0 0 0 0 0 0 0 0 0 0];
        for j = 0 : 9
            %Reading the digits for comparison - קריאת הספרות להשוואה
            Filter = imread(strcat(int2str(j),'.jpg'));
            % Convert to BW image - המרה לתמונת שחור לבן
            Filter = im2bw(Filter);
            %Resize our digit to number size - שינוי גודל הספרה שלנו לגודל מספר
            Num = imresize(Num, size(Filter));
            % Compare a number to our digit bit by bit - השוואה ביט אחר ביט
            Similar = bitand(Num,Filter); 
            %Equlaize num to digit by Cross-correlation - השוואה ע"י קרוס-קורלצייה
            Similar2cross=max(max(abs(normxcorr2(Num,Filter))));
            Similarity(j+1) = Similar2cross; 
        end

        num = 0; 
        for j = 0 : 9
            % Gets the digit that best matches an unidentified digit - מקבל את המספר הכי קרוב לספרה
            if(Similarity(j+1)==max(Similarity))
                num = j;
            end
        end
        %Inserting the number into the vector - הכנסת המספר לוקטור
        PlateNumber(Counter) = num; 
        Counter=Counter+1;
    end
end

%% Saving the license plate number in a text file           
%open new window -  פתיחת חלון חדש
figure(1);
%Showing the original image
imshow(Image);
%save number vector into text file - שמירת הווקטור בקובץ טקסט
dlmwrite('PlateNumber.txt', PlateNumber); 
% open file - פתיחת הקובץ   
winopen('PlateNumber.txt'); 