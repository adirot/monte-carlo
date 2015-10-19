function colorPlot(x,y)
% every row in x is an x axis to plot
% every rew in y is a y axis to plot
% this will plot all [x, y]'s in different colors
    [numOfPlots,~] = size(x);
    
  
   cmap = hsv(numOfPlots);
   
   figure; hold on;
   
   for i = 1:numOfPlots
       plot(x(i,:),y(i,:),'color',cmap(i,:));
   end
   
end
       