function [ beadpos ] = RandBeadPos( range, cbeads, rbead, fvcell )
%RANDBEADPOS  returns a list of random bead positions. Beads are checked
%             for collisions with other beads and with the cell body using
%             the function inpolyhedron by Sven Holcombe.
%Input:
%  <range>    range of bead positions [xmin,xmax;ymin,ymax;zmin,zmax]
%  <cbeads>   bead concentration
%  <rbead>    bead radius
%  <fvcell>   cell model patch object
%Output:
%  <beadpos>  bead positions in the format [x,y,z]
%CL

% Seeds the random number generator based on the current time
rng('shuffle')

% Cubic boundary surrounding the cell to increase collision check efficiency
bmax = [max(fvcell.vertices(:,1))+rbead,max(fvcell.vertices(:,2))+rbead,...
    max(fvcell.vertices(:,3))+rbead];
bmin = [min(fvcell.vertices(:,1))-rbead,min(fvcell.vertices(:,2))-rbead,...
    min(fvcell.vertices(:,3))-rbead];

% Generate a basic bead with radius 1 at [0,0,0]
[xb,yb,zb] = sphere(5);

% Determine the number of beads
nbeads = round(prod(range(:,2)-range(:,1))*cbeads);

% Generate the list of bead positions
beadpos = zeros(nbeads,3);
for i=1:nbeads
    
    redo = true;
    
    while redo

        redo = false;
        
        % Generate a random bead position within the specified range
        new=rand(1,3);
        new(1) = range(1,1)+rbead+new(1).*(range(1,2)-range(1,1)-2*rbead);
        new(2) = range(2,1)+rbead+new(2).*(range(2,2)-range(2,1)-2*rbead);
        new(3) = range(3,1)+rbead+new(3).*(range(3,2)-range(3,1)-2*rbead);
        beadpos(i,:) = new;
        check = beadpos(1:i-1,:);

        % Generate a bead with radius <rbead> at the new position
        bead = zeros(size(xb(:),1),3);
        bead(:,1) = rbead*xb(:)+new(1,1);
        bead(:,2) = rbead*yb(:)+new(1,2);
        bead(:,3) = rbead*zb(:)+new(1,3);
        
        % Check if the bead is inside the cube surrounding the cell
        I = intersect(find(new>bmin),find(new<bmax));
        if size(I,2) == 3
            % Check if the bead is intersecting the cell body
            inshp = inpolyhedron(fvcell,bead);
            if any(inshp)
                redo = true;
            end
        end
            
        % Check if the bead is intersecting other beads
        if ~isempty(check)
            for d = 1:3
                I = find(abs(check(:,d)-new(d)) < (2*rbead));
                check = check(I,:);
                if isempty(I)
                    break
                end     
            end
            if ~isempty(I)
                tmp = ones(size(check,1),1)*new;
                tmp = sqrt((tmp(:,1)-check(:,1)).^2+(tmp(:,2)-check(:,2)).^2+...
                    (tmp(:,3)-check(:,3)).^2);
                I = find(tmp<=(2*rbead));
                if ~isempty(I)
                    redo = true;
                end
            end
        end
        
    end
    
end

end