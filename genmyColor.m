
% generates array of random colors of size noColors
% saves them to "myColor" 
% structure of colors = [R1 G1 B1; R2 G2 B2; ...]

noColors=96

myColor=[];
myRGBCount = ceil(noColors^(1/3))
for i1 = [1:myRGBCount] for i2 = [myRGBCount:-1:1] for i3 = [1:myRGBCount]
    myColor=[myColor; [i1/myRGBCount i2/myRGBCount i3/myRGBCount]];
end; end; end;
myColor=myColor(randperm(length(myColor)),:)
myColor    

save([myRootDir myScriptDir 'myColor.mat'],'myColor')