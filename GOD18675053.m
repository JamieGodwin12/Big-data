%opens dataset 
fid = fopen('car_data.txt','r');
Data =[];
while 1==1
    %gets line from file
    car = fgetl(fid);  
    if car == -1
        break;
    end    
    %converts string to array 
    a = strsplit(car,'\t');
    Data = [Data;a];        
end
i = height(Data);
for j = height(Data):-1:1
    for s = 1:9
        %remove all of the NA values 
        Na = strcmp(Data(j,s),'NA');
        if Na == 1
            Data(j,:) = [];
        end
    end
end
%this gets the size of data that has been imported
[rows,collums] = size(Data);
%sets which data to be included
q = 1:1:rows;
p = 1:1:8;
%output only double vaules. removes strings 
stringA = string(Data);
Dcar = double(stringA(q,p));
Dcar(1,:)= [];
%temporary vaules declared
mpgmeantemporary = 0;
accelerationtemporary = 0;
horsepowertemporary = 0;
weighttemporaryorary = 0;
%individual arrays
mpgarray = Dcar(:,1);
accelerationarray = Dcar(:,6);
horsearray = Dcar(:,4);
weightarray = Dcar(:,5);
%sorted arrays
mpgarray1 = insertionSort(mpgarray);
accelerationarray1 = insertionSort(accelerationarray);
horsearray1 = insertionSort(horsearray);
weightarray1 = insertionSort(weightarray);
%counter 
for ctr = 1:392
    mpgmeantemporary = mpgmeantemporary + mpgarray1(ctr);
    accelerationtemporary = accelerationtemporary + accelerationarray1(ctr);
    horsepowertemporary = horsepowertemporary + horsearray1(ctr);
    weighttemporaryorary = weighttemporaryorary + weightarray1(ctr);
end
%calculate mean
mpgmean = meanfunction(mpgarray1);
accelerationmean = meanfunction(accelerationarray1);
horsemean = meanfunction(horsearray1);
weightmean = meanfunction(weightarray1);

%median
mpgmedian = Median(mpgarray1);
accelmedian = Median(accelerationarray1);
horsemedian = Median(horsearray1);
weightmedian = uint16(Median(weightarray1));
%minimum
mpgmin = mpgarray1(1);
accelmin = accelerationarray1(1);
horsemin = horsearray1(1);
weightmin = weightarray1(1);
%maximum
mpgmax = mpgarray1(392);
accelmax = accelerationarray1(392);
horsemax = horsearray1(392);
weightmax = weightarray1(392);
%standard deviation
mpgstdeviation = stdev(mpgarray1);
acceldeviation = stdev(accelerationarray1);
horsedeviation = stdev(horsearray1);
weightdeviation = stdev(weightarray1);

%table being made for variables  
mean = [mpgmean;accelerationmean;horsemean;weightmean];
median = [mpgmedian;accelmedian;horsemedian;weightmedian];
min = [mpgmin;accelmin;horsemin;weightmin];
max = [mpgmax;accelmax;horsemax;weightmax];
standarddeviation = [mpgstdeviation;acceldeviation;horsedeviation;weightdeviation];
%table variables
variables = ["MPG";"Acceleration";"Horsepower";"Weight"];
%table heading
vars = table(variables,mean,median,min,max,standarddeviation);
disp(vars)
%ploting the data using boxplot
figure
boxplot(mpgarray)
title('MPG')
figure
boxplot(accelerationarray)
title('Acceleration')
figure
boxplot(horsearray)
title('Horsepower')
figure
boxplot(weightarray)
title('Weight')
% plotting the scatter graphs
figure
scatter(accelerationarray,mpgarray)
title('Acceleration vs MPG')
xlabel('Acceleration');
ylabel('MPG');
figure
scatter(mpgarray,horsearray)
title('Horsepower vs MPG')
xlabel('MPG');
ylabel('Horsepower');
figure
scatter(weightarray,horsearray)
title('Weight vs Horsepower')
xlabel('Weight');
ylabel('Horsepower');
%plotting the density graphs 
figure
histogram(mpgarray)
title('MPG density graph')
figure
histogram(accelerationarray)
title('Acceleration density graph')
figure
histogram(horsearray)
title('Horsepower density graph')
figure
histogram(weightarray)
title('Weight density graph')

%sets each up 2 arrays. for test and traing datasets
%splits them into 70 and 30
Dcar = shuffler(Dcar);
trainset = string.empty();
testset = string.empty();
intT = round((size(Dcar, 1)-1)*0.7);
for i = 1:size(Dcar, 2) 

    trainset(1, i) = Dcar(1, i);
    testset(1, i) = Dcar(1, i);
end
for e = 2:intT+1
    for i = 1:size(Dcar, 2)
        trainset(e, i) = Dcar(e, i);
    end
end
testR = 1;
for e = intT+2:size(Dcar, 1)
    testR = testR+1;
    for i = 1:size(Dcar, 2)
        testset(testR, i) = Dcar(e, i);
    end
end

%which collums are used in regression fuction
figure
[m , b] = Regression(trainset, 6,1);
figure
[m , b] = Regression(trainset, 4,1);
figure
[m , b] = Regression(trainset, 5,4);

% function matrix = CorrelationMatrix(Dcar)
%     matrix = double.empty(0, 4)
%     for x = l:size(Dcar,2)
%         for y = l:size(Dcar,2)
%             meanx = meanfunction(Dcar,x);
%             meany = meanfunction(Dcar,y);
%             totalnumerator= 0;
%             totaldenominatorx= 0;
%             totaldenominatory= 0;
%             for k = 2:size(Dcar,1)
%                 diffx = double(Dcar(k , x)- meanx;
%                 diffy = double(Dcar(k , x)- meany;
%                 totaldiff = diffx * diffy;
%                 totalnumerator = totalnumerator + totaldiff
%                 totaldenominatorx = totaldenominatorx +(diffx * diffy);
%                 totaldenominatory = totaldenominatory +(diffx * diffy);
%             end 
%             totaldenmoinator = totaldenmoinatorx * totaldenmoinatory;
%             matrix(x, y) = (totalnumerator / sqrt(totaldenmoinator));
%         end
%     end
% end


%insertionSort algo for sorting 
function x = insertionSort(x)
for is=2:length(x)
    b = is;
    while (b > 1 && x(b-1) > x(b))
        [x(b),x(b-1)] = deal(x(b-1),x(b));
        b = b - 1;
    end
end
end
%this is the shuffle functiuon 
function shuffle = shuffler(Dcar)
         for a = 2:size(Dcar,1)
            temp = Dcar(a,:);
            randy = randi ([2 size(Dcar, 1)]);
            Dcar(a,:)= Dcar(randy,:);
            Dcar(randy,:) = temp;
         end
         shuffle = Dcar;
end
%median function 
function med = Median(x)
    n = length(x);
    t = (n+1)/2;
    med = (x(floor(t))+x(ceil(t)))/2;
end
%standard deviation function 
function stdeviation = stdev(x)
    AvgX = sumsfunction(x)/length(x);
    stdeviation = sqrt(sumsfunction((x-AvgX).^2)/(length(x)-1));
end
%sum function
function sums = sumsfunction(x)
         sums =0;
         for i = 1 : length(x)
             sums = sums + x(i);
         end
end
%mean function
function mean = meanfunction(x)
         l = length(x);
         mean = sumsfunction(x) /l;
        
end
%regression
function [m,b] = Regression(trainset,q,p)
         %seting the sum vaules and geting x and y vaules from trainset 
         x = trainset(2:size(trainset, 1), q);
         y = trainset(2:size(trainset, 1), p);
         sumx = 0;
         sumy = 0;
         sumxy = 0;
         sumx2 = 0;
         %converts string into dobule to be used for caulations
         for i = 1 : size(x,1)
             sumx = sumx+str2double(x(i));
             sumy = sumy+str2double(y(i));
             sumxy = sumxy+str2double(x(i))*str2double(y(i));
             sumx2 = sumx2+str2double(x(i))*str2double(x(i));
         end
         %cacualtes the slope and intercept for simple L regression
         f = (size(x,1)*sumxy)-(sumx*sumy);
         h = (size(x,1)*sumx2)-(sumx*sumx);
         m = f/h; 
         b = (sumy-(m*sumx))/size(x,1);
end


