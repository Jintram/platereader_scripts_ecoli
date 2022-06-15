
% generates array of random colors of size noColors
% saves them to "myColor" 
% structure of colors = [R1 G1 B1; R2 G2 B2; ...]

% Settings ***************
% No. colors you want to generate
noColors=96
% ***************

myColor=[];
% How many different R, G and B values needed
myRGBCount = ceil(noColors^(1/3)); 

% Generate colors - simply all permutations of different combinations of 
% R = {1..myRGBCount}/(RGBCount+1). The +1 in the denominator is to avoid
% values of myColor[i] = [1, 1, 1].
for i1 = [1:myRGBCount] for i2 = [1:myRGBCount] for i3 = [1:myRGBCount]
    myColor=[myColor; [i1/(myRGBCount+1) i2/(myRGBCount+1) i3/(myRGBCount+1)]];
end; end; end;

% Mix order up to avoid similar colors next to each other.
myColor=myColor(randperm(length(myColor)),:)

% Show to user
myColor    

% Save colors
save([myRootDir myScriptDir 'myColor.mat'],'myColor')