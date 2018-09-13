function [warncode] = WRITE_sat(fileOUT,coordVERTICES,coordNORMALS)
% WRITE_sat  Write an ACIS SAT (v4.0) file
%==========================================================================
% FILENAME:          WRITE_sat.m
% AUTHOR:            Adam H. Aitkenhead
% INSTITUTION:       The Christie NHS Foundation Trust
% CONTACT:           adam.aitkenhead@christie.nhs.uk
% PURPOSE:           Write an ACIS SAT (v4.0) file, which can then be
%                    imported into various CAD packages as a solid model.
%
% USAGE:             WRITE_sat(fileOUT,coordVERTICES,coordNORMALS)
%                    writes the ACIS SAT file <fileOUT> based on the facet
%                    vertices <coordVERTICES> and normals <coordNORMALS>.
%
%                    <fileOUT> is a text string providing the filename of
%                    the SAT file.  eg. 'sample.sat'
%
%                    <coordVERTICES> is an Nx3x3 array defining the
%                    vertex positions for each facet, with:
%                       - 1 row for each facet
%                       - 3 cols for the x,y,z coordinates
%                       - 3 pages for the three vertices
%
%                    <coordNORMALS> is an Nx3 array defining the normal
%                    vector for each facet, with:
%                       - 1 row for each facet
%                       - 3 cols for the x,y,z components of the vector
%
% STL REQUIREMENTS:  The STL data must meet the following criteria:
%                    1. All facets are triangular.
%                    2. No duplicate or overlapping facets exist.
%                    3. The mesh is properly closed.
%                    4. The facet normals are properly defined.
%                    5. For each facet edge, there must exist an odd number
%                       of identical edges on adjacent facets.
%
%                    Although uncommon, some STL meshes may not meet
%                    requirement 5 and cannot be converted to SAT using
%                    this code.
%
% REFERENCES:        For a description of the ACIS format, refer to:
%                    http://local.wasp.uwa.edu.au/~pbourke/dataformats/sat/sat.pdf
%==========================================================================

%==========================================================================
% VERSION  USER  CHANGES
% -------  ----  -------
% 100331   AHA   Original version
% 100406   AHA   Added variable <edgelistPAIRS> to record paired edges.
%                Speeds up execution by up to 15% for large meshes.
% 100429   AHA   Added a check for duplicate facets
% 100506   AHA   Replaced calls to <ismember> with calls to <find>,
%                speeding up the code considerably.
% 100513   AHA   Improved error checking
% 100527   AHA   Added the ability to handle meshes where on odd number
%                (>1) of facets share an identical edge.
% 110822   AHA   Added error checking to ensure that each for each facet
%                edge, there exists an odd number of identical edges on
%                adjacent facets (STL requirement 5).
% 120606   AHA   If the mesh does not meet the requirements, a plot is
%                displayed to show which facets are problematic.
% 150923   CL    Added progress information
%==========================================================================

warncode = false(1);

%===========================
% Create lists of the vertices
%===========================

%Create a list of the unique vertices:
vertices    = [coordVERTICES(:,:,1);coordVERTICES(:,:,2);coordVERTICES(:,:,3)];
vertices    = unique(vertices,'rows');

vertexCOUNT = size(vertices,1);           %Total number of vertices
facetCOUNT  = size(coordVERTICES,1);      %Total number of facets

%Create a list of the vertex indices used for each facet:
facetvertexlist = zeros(facetCOUNT,3);
for loopF = 1:facetCOUNT
  for loopV = 1:3
    vlisttemp = find( vertices(:,1)==coordVERTICES(loopF,1,loopV) );
    vlisttemp = vlisttemp( find( vertices(vlisttemp,2)==coordVERTICES(loopF,2,loopV) ) );
    vlisttemp = vlisttemp( find( vertices(vlisttemp,3)==coordVERTICES(loopF,3,loopV),1,'last') );
    facetvertexlist(loopF,loopV) = vlisttemp;
  end
end

%===========================
% Check for duplicate facets
%===========================

if size(facetvertexlist,1) > size(unique(sort(facetvertexlist,2),'rows'),1)
  
  [~,uniqueIND,~] = unique(sort(facetvertexlist,2),'rows');
  duplicatefacets = ~ismember(1:size(facetvertexlist,1),uniqueIND);
  
  figure
  hold on
  [hALL] = PLOT_3D_stl_patch(coordVERTICES);
  [herr] = PLOT_3D_stl_patch(coordVERTICES(duplicatefacets,:,:));
  set(herr,'FaceColor',[1,0.5,0.5])
  hold off
      
  disp('=======================================================================')
  disp('ERROR:  STL does not meet requirement 2:  There must be no duplicate');
  disp('        facets.');
  disp('        The problem facet(s) are highlighted in the STL viewer.');
  disp('=======================================================================')
  warncode(2) = true;
  return
  
end

%===========================
% Create a list of all facet edges, labelled using the vertex indices of the two ends of the edge
%===========================

edgelistLONG = zeros(facetCOUNT,2,3);
edgelistLONG(:,:,1) = facetvertexlist(:,1:2);
edgelistLONG(:,:,2) = facetvertexlist(:,2:3);
edgelistLONG(:,:,3) = facetvertexlist(:,[3,1]);
edgelistSHORT = [edgelistLONG(:,:,1);edgelistLONG(:,:,2);edgelistLONG(:,:,3)];
edgelistSHORT = sort(edgelistSHORT,2);
edgelistSHORT = unique(edgelistSHORT,'rows');

edgeCOUNT = size(edgelistSHORT,1);    %Total number of edges

%======================
% Display file information
%======================

% disp(' ')
% disp([' Number of facets:    ',num2str(facetCOUNT)])
% disp([' Number of edges:     ',num2str(edgeCOUNT)])
% disp([' Number of vertices:  ',num2str(vertexCOUNT)])

%===========================
%Calculate the number of lines that will be required for the each section of the SAT file
%===========================

numlinesINTROSECTION         = 4;
numlinesFACETSECTION         = 6*facetCOUNT;    %Requires assumption:  All facets are triangular
numlinesEDGESECTION          = edgeCOUNT;
numlinesSTRAIGHTCURVESECTION = edgeCOUNT;

%======================
% Create a new SAT file for writing
%======================

% disp([' Exporting ',fileOUT]);
fidOUT = fopen(fileOUT,'w');

%======================
% WRITE header
%======================

headerdate = datestr(now);

fprintf(fidOUT,['400 0 1 0\n',...
                '72 adam.aitkenhead@physics.cr.man.ac.uk (The Christie NHS Foundation Trust) 8 ACIS 4.0 20 %s\n',...
                '1 9.9999999999999995e-007 1e-010\n'],...
                headerdate...
       );

%======================
% WRITE intro
%======================

fprintf(fidOUT,['-0 body $1 $2 $-1 $-1 #\n',...
                '-1 name_attrib-gen-attrib $-1 $-1 $-1 $0 keep keep_kept ignore 4 Part #\n',...
                '-2 lump $3 $-1 $4 $0 #\n',...
                '-3 rgb_color-st-attrib $-1 $-1 $-1 $2 0.80000000000000000 0.80000000000000000 0.80000000000000000 #\n',...    % 204/255 = 0.8
                '-4 shell $-1 $-1 $-1 $5 $-1 $2 #\n']...
       );

%======================
% WRITE facet
%======================

lineIND = 5;    %The line number in the SAT file to begin the first facet.

%Prepare an array to record the first use of each edge.  Each row in this
%array corresponds to a row in edgelistSHORT.  The number in the first
%column denotes the first facet to use each edge, and the number in the
%second column denotes the order of the facet edge.
edgelistPAIRS = zeros(size(edgelistSHORT,1),2);

%Prepare an array which will be used in the handling of meshes where >2
%facets share an identical edge.  Each edge must be paired with 1 other
%edge.  This array will be used to keep track of edges which have already
%been paired up, to prevent them being paired with additional edges.  The
%array is built on-the-fly as these edges are encountered, but since the
%number of edges that fall into this category is usually small, the speed
%penalty should be minimal.
matchedgelistUSED = [];

clearmsg = '';

for loopFACE = 1:facetCOUNT
    
  %Show progress
  msg = sprintf('Write SAT... %.0f%%',loopFACE/facetCOUNT*100);
  fprintf('%s',[clearmsg,msg])
  clearmsg = repmat(sprintf('\b'),1,length(msg));
    
  %Coordinates of first vertex:
  vertex1xyz  = coordVERTICES(loopFACE,1:3,1);
  vertex2xyz  = coordVERTICES(loopFACE,1:3,2);
  
  %Vector normal of the facet:
  facetNORMAL = coordNORMALS(loopFACE,:);
  
  %Vector between vertex 1 and vertex 2 on the first edge:
  vectorV1V2  = vertex2xyz - vertex1xyz;
  
  faceLINENUM = lineIND;
  if loopFACE < facetCOUNT
    facenextLINENUM = lineIND+6;
  else
    facenextLINENUM = -1;
  end
  fprintf(fidOUT,['-%d face $-1 $%d $%d $4 $-1 $%d forward single #\n',...
                  '-%d loop $-1 $-1 $%d $%d #\n',...
                  '-%d plane-surface $-1 %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f forward_v I I I I #\n'],...
                  lineIND,facenextLINENUM,lineIND+1,lineIND+2,...
                  lineIND+1,faceLINENUM+3,faceLINENUM,...
                  lineIND+2,vertex1xyz(1:3),facetNORMAL(1:3),vectorV1V2(1:3)...
         );
  lineIND = lineIND+3;
  
  %Loop around the three edges of the facet, beginning with the first edge.
  %The edges of each facet will be traced in the following order:  1->2  2->3  3->1
  currentvertexA = facetvertexlist(loopFACE,1);
  currentvertexB = facetvertexlist(loopFACE,2);
  coedgenext  = 2;
  coedgefinal = 3;
  for loopCOEDGE = 1:3
    
    %Check whether this edge is in forward or reverse:
    [currentvertexAsorted,currentvertexBsorted] = deal(min([currentvertexA,currentvertexB]),max([currentvertexA,currentvertexB]));
    if currentvertexA == currentvertexAsorted
      direction = 'forward';
    else
      direction = 'reversed';
    end
    
    %Find the line in the EDGE section which relates to the current edge:
    reftoEDGElist    = find(edgelistSHORT(:,1)==currentvertexAsorted);
    reftoEDGElist    = reftoEDGElist( find(edgelistSHORT(reftoEDGElist,2)==currentvertexBsorted,1,'last') );
    reftoEDGEsection = numlinesINTROSECTION + numlinesFACETSECTION + reftoEDGElist;
    
    %Check that this edge has not already been matched with another edge:
    if edgelistPAIRS(reftoEDGElist,1)==0
        
      %Find the COEDGEs in adjacent facets which matches the current edge:
      [rowA,colA] = find(facetvertexlist==currentvertexA);
      [rowB,colB] = find(facetvertexlist==currentvertexB);
      [matchingFACET,intercoA,intercoB] = intersect(rowA,rowB);
      colA = colA(intercoA);
      colB = colB(intercoB);
      [matchingFACET,diffcoIND] = setdiff(matchingFACET,loopFACE);
      colA = colA(diffcoIND);
      colB = colB(diffcoIND);
      colAplusB = colA+colB;
      
      % The order of the matching edge in the other SAT facet.  (Facets are looped anti-clockwise)
      COEDGEorder = zeros(size(colAplusB,1),1);
      COEDGEorder(colAplusB==3) = 1;  % 1st and 2nd vertex => Total = 3
      COEDGEorder(colAplusB==5) = 2;  % 2nd and 3rd vertex => Total = 5
      COEDGEorder(colAplusB==4) = 3;  % 3rd and 1st vertex => Total = 4

      %Record the first occurrence of this edge:
      edgelistPAIRS(reftoEDGElist,1:2) = [loopFACE,loopCOEDGE];   % row->edge# : col1->facet col2->edge# of originaloccurence
    
    elseif isempty(matchedgelistUSED) || numel(find(matchedgelistUSED(:,1)==loopFACE & matchedgelistUSED(:,2)==loopCOEDGE)==0)
      %If the edge has already been matched to ONE other then there is no
      %need to repeat the search, just look up the original result.
      matchingFACET = edgelistPAIRS(reftoEDGElist,1);
      COEDGEorder   = edgelistPAIRS(reftoEDGElist,2);
    end
    
    if isempty(matchingFACET)==1  ||  rem(numel(matchingFACET),2)==0
        
      fclose(fidOUT);
      
      figure
      hold on
      [hALL] = PLOT_3D_stl_patch(coordVERTICES);
      [herr] = PLOT_3D_stl_patch(coordVERTICES(loopFACE,:,:));
      set(herr,'FaceColor',[1,0.5,0.5])
      plot3([vertex1xyz(1),vertex2xyz(1)],[vertex1xyz(2),vertex2xyz(2)],[vertex1xyz(3),vertex2xyz(3)],'r-','LineWidth',2)
      hold off
      
      disp('\n=======================================================================')
      disp('ERROR:  The STL does not meet requirement 5:  For each facet edge,');
      disp('        there must exist an odd number of identical edges on adjacent');
      disp('        facets.');
      disp('        The problem facet and edge are highlighted in the STL viewer.');
      disp('=======================================================================')
      warncode(5) = true;
      
      return
         
    elseif numel(matchingFACET)>1
      %This part of the code deals with the situation where an edge has
      %more than coincident edge in adjacent facets.  Only one of these
      %potential matching edges can be chosen as the matched edge.  If an
      %inappropriate choice is made, the final object will not form a
      %closed solid body.  An appropriate choice is as follows:
      % 1- The matching edge must be looped around in the opposite
      %    direction to the edge in the original facet.
      % 2- The matching edge must not have been previously matched to
      %    another edge.
      
      %The current facet is given by:                         loopFACE 
      %The relevant edge on the current facet is goven by:    loopCOEDGE
      %The list of facets with a matching edge is given by:   matchingFACET
      %The relevant edge on each of these facets is given by: COEDGEorder
      
      edgecurrent = edgelistLONG(loopFACE,:,loopCOEDGE);
      
      edgeprospective = [];
      for loopa = 1:size(matchingFACET,1)
        edgeprospective = [edgeprospective;edgelistLONG(matchingFACET(loopa),:,COEDGEorder(loopa))];
      end
      
      edgeopposite = [];
      for loopa = 1:size(matchingFACET,1)
        edgeopposite = [edgeopposite;~isequal(edgecurrent,edgeprospective(loopa,:))];
      end
      
      matchingFACET = matchingFACET(logical(edgeopposite));
      COEDGEorder   = COEDGEorder(logical(edgeopposite));

      if isempty(matchedgelistUSED)==0
        temp = setdiff([matchingFACET,COEDGEorder],matchedgelistUSED,'rows');
        matchingFACET = temp(:,1);
        COEDGEorder   = temp(:,2);
      end
      
      matchingFACET = matchingFACET(1);
      COEDGEorder   = COEDGEorder(1);
      
      %Add the current edge and its matched edge to the list of paired
      %edges to make sure they are not re-used later.
      matchedgelistUSED = [matchedgelistUSED; loopFACE,loopCOEDGE; matchingFACET,COEDGEorder];
      
    end
    
    %Find the line number of the matching COEDGE:
    reftomatchingCOEDGEinFACETsection = numlinesINTROSECTION + 6*(matchingFACET-1) + 3 + COEDGEorder;
    reftonextedge  = faceLINENUM+2+coedgenext;
    reftofinaledge = faceLINENUM+2+coedgefinal;
    
    %Write the COEDGE line:
    fprintf(fidOUT,'-%d coedge $-1 $%d $%d $%d $%d %s $%d $-1 #\n',lineIND,reftonextedge,reftofinaledge,reftomatchingCOEDGEinFACETsection,reftoEDGEsection,direction,faceLINENUM+1);

    %Define the vertices in preparation for the next iteration of the loop 'loopCOEDGE':
    currentvertexA = facetvertexlist(loopFACE,coedgenext);
    currentvertexB = facetvertexlist(loopFACE,coedgefinal);
    
    %Define the linenumber shifts in preparation for the next iteration of the loop 'loopCOEDGE':
    lineIND = lineIND+1;
    coedgenext = rem(coedgenext+1,3);
    if coedgenext==0
      coedgenext = 3;
    end
    coedgefinal = rem(coedgefinal+1,3);
    if coedgefinal==0
      coedgefinal = 3;
    end
    
  end
  
end

%======================
% WRITE edge
%======================

%Calculate the line number references to the two vertices for each edge:
vertexLINES = edgelistSHORT + lineIND + numlinesEDGESECTION + numlinesSTRAIGHTCURVESECTION -1;

for loopEDGE = 1:edgeCOUNT
    
  vertexA = edgelistSHORT(loopEDGE,1);
  vertexB = edgelistSHORT(loopEDGE,2);
  
  %Find the first facet which contains this edge:
  [rowA,colA] = find(facetvertexlist==vertexA);
  [rowB,colB] = find(facetvertexlist(rowA,:)==vertexB);
  facetswiththisedgeALL   = sort(rowA(rowB));
  facetswiththisedgeFIRST = facetswiththisedgeALL(1);
  
  %Find which edge in that facet matches the edge in question:
  edgelistinthatfacet = [edgelistLONG(facetswiththisedgeFIRST,1:2,1);edgelistLONG(facetswiththisedgeFIRST,1:2,2);edgelistLONG(facetswiththisedgeFIRST,1:2,3)];
  indexofedgeinfacet = find(edgelistinthatfacet(:,1)==vertexA & edgelistinthatfacet(:,2)==vertexB,1,'last');
  if isempty(indexofedgeinfacet)
    indexofedgeinfacet = find(edgelistinthatfacet(:,1)==vertexB & edgelistinthatfacet(:,2)==vertexA,1,'last');
  end

  %Line number of the first coedge which contains the edge in question:
  coedge = numlinesINTROSECTION + 6*(facetswiththisedgeFIRST-1) + 3 + indexofedgeinfacet;
    
  STRAIGHTCURVElinenumber = lineIND + edgeCOUNT;
  fprintf(fidOUT,'-%d edge $-1 $%d $%d $%d $%d forward #\n',lineIND,vertexLINES(loopEDGE,1),vertexLINES(loopEDGE,2),coedge,STRAIGHTCURVElinenumber);
  lineIND = lineIND+1;
end

%======================
% WRITE straight-curve
%======================

for loopSTRAIGHT = 1:edgeCOUNT

  startvertexIND = edgelistSHORT(loopSTRAIGHT,1);
  endvertexIND   = edgelistSHORT(loopSTRAIGHT,2);
 
  [facetID,vertexID] = find(facetvertexlist==startvertexIND,1,'first');
  startCO            = coordVERTICES(facetID,:,vertexID);

  [facetID,vertexID] = find(facetvertexlist==endvertexIND,1,'first');
  endCO              = coordVERTICES(facetID,:,vertexID);

  lengthCO = endCO - startCO;

  fprintf(fidOUT,'-%d straight-curve $-1 %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f I I #\n',lineIND,startCO(1:3),lengthCO(1:3));
  lineIND = lineIND+1;

end

%======================
% WRITE vertex
%======================

for loopVERTEX = 1:vertexCOUNT
  [resultCOL,resultROW]      = find(vertexLINES'==lineIND,1,'first');   %Transposing vertexLINES means that the 'find' command looks along each row in turn, rather than looking down each column in turn
  EDGEtoVERTEXfirstreference = resultROW + numlinesINTROSECTION + numlinesFACETSECTION;
  POINTlinenumber            = lineIND + vertexCOUNT;
  fprintf(fidOUT,'-%d vertex $-1 $%d $%d #\n',lineIND,EDGEtoVERTEXfirstreference,POINTlinenumber);
  lineIND = lineIND+1;
end

%======================
% WRITE point
%======================

for loopPOINT = 1:vertexCOUNT
  fprintf(fidOUT,'-%d point $-1 %8.6f %8.6f %8.6f #\n',lineIND,vertices(loopPOINT,1:3));
  lineIND = lineIND+1;
end

%======================
% WRITE footer
%======================

fprintf(fidOUT,'End-of-ACIS-data');

%======================
% Close the SAT file
%======================

fclose(fidOUT);
