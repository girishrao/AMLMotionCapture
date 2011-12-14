% takes a cell from mocap data and animates the markers

function[] = animate(Markers)

[dim,frames] = size(Markers);

points = dim/3;

%find max and min dims to normalize plot
dims(1) = min(min(Markers(1:3:size(Markers,1),:)));
dims(3) = min(min(Markers(2:3:size(Markers,1),:)));
dims(5) = min(min(Markers(3:3:size(Markers,1),:)));
dims(2) = max(max(Markers(1:3:size(Markers,1),:)));
dims(4) = max(max(Markers(2:3:size(Markers,1),:)));
dims(6) = max(max(Markers(3:3:size(Markers,1),:)));

 mindim = min(dims);
 maxdim = max(dims);


for i=1:frames

    scatter3(Markers(1:3:dim,i),Markers(2:3:dim,i),Markers(3:3:dim,i));
    
    axis([mindim maxdim mindim maxdim mindim maxdim]);

    axis('square');
    drawnow;
    
end