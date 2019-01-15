data = xlsread('Filtered2IronCarbonClustersEnergiesTesting.xlsx');


Q=length(data);


fid=fopen('Filtered3IronCarbonClustersEnerpiesTesting.xls', 'w');


for iQ=1:1:Q
    if data(iQ,2)==1.125 || data(iQ,3) ==1.125|| data(iQ,4)==1.125 || data(iQ,5)==1.125 || data(iQ,6)==1.125 || ...
       data(iQ,7)==1.125 ||  data(iQ,8)==1.125 ||  data(iQ,9)==1.125 ||  data(iQ,10)==1.125 
       continue;
    else
        fprintf(fid,'%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', data(iQ,1),data(iQ,2),data(iQ,3),data(iQ,4),...
            data(iQ,5),data(iQ,6),data(iQ,7),data(iQ,8),data(iQ,9),data(iQ,10),data(iQ,11),data(iQ,12),data(iQ,13));
    end
end
fclose('all');
