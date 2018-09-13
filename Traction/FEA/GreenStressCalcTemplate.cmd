/BATCH
/FILNAME,TFM3D,0

!*********************************************************************************************
!******************      FEA for Determining Green's Function u=G*T      *********************
!*********************************************************************************************
!CL

/PREP7 						!Enter preprocessor
/CWD,'***working directory***'			!Change working directory
/UNITS,uMKS					!Use MKS system of units

!Read query points from text file
!*********************************************************************************************
/INQUIRE,nlines,LINES,upos,txt			!Get number of lines
nqloc = nlines-1				!Get number of query points
*DIM,qloctable,TABLE,nqloc,3			!Define table to store query points
*TREAD,qloctable,upos,txt			!Read query point locations from file
*DIM,qloc,ARRAY,nqloc,3				!Define array for query points
*VCOL,3,3					!Number of columns to copy
*MFUN,qloc(1,1),COPY,qloctable(1,1)		!Copy query points into array

!Geometry
!*********************************************************************************************
~SATIN,'***modelname***','sat','***path***',SOLIDS,0
						!Load cell model				
BLOCK,-200,200,-200,200,-200,200		!Create 400x400x400um substrate bulk
VSBV,2,1					!Subtract cell from the substrate cube

!Element types
!*********************************************************************************************
ET,1,SOLID285					!4-node tetrahedral element
ET,2,SURF154					!Surface effect element
KEYOPT,2,4,1					!...without midside nodes
KEYOPT,2,11,2					!...and pressure vector on full area
KEYOPT,2,12,0					!...applied regardless of element orientation

!Material properties
!*********************************************************************************************
MP,EX,1,'***E***'   				!Young's modulus [Pa] of the substrate
MP,PRXY,1,'***nu***'				!Poisson's ratio of the substrate

!Meshing
!*********************************************************************************************
LSEL,S,LOC,X,-200				!Group substrate edges for mesh definition
LSEL,A,LOC,X,200	
LSEL,A,LOC,Y,-200
LSEL,A,LOC,Y,200
LSEL,A,LOC,Z,-200
LSEL,A,LOC,Z,200
CM,edges,LINE
TYPE,1						!Activate element type 1 (SOLID285)
MAT,1						!Activate material 1
LESIZE,edges,***size***				!Mesh spacing on substrate edges
ESIZE,***size***				!Define minimum mesh spacing
VMESH,3						!Mesh the substrate (volume 3)

!Apply DOF constraints (fixed substrate bottom)
!*********************************************************************************************
NSEL,S,LOC,Z,-200				!Select nodes on the substrate bottom
D,ALL,UX,0					!Fix selected nodes in x direction
D,ALL,UY,0					!Fix selected nodes in y direction
D,ALL,UZ,0					!Fix selected nodes in z direction 

!Create load steps: Apply unit tractions separately for each interior face and dimension
!*********************************************************************************************
ASEL,S,LOC,X,-150,150				!Select all areas of the interior surface
ASEL,R,LOC,Y,-150,150
ASEL,R,LOC,Z,-150,150
*GET,nareas,AREA,,COUNT				!Get number of selected areas (faces)
*GET,amin,AREA,,NUM,MIN				!Get minimum area number
*DIM,areas,ARRAY,nareas,1			!Array for area numbers
*DIM,fcent,ARRAY,nareas,3			!Array for face centroid locations
*DO,I,1,nareas					!Store area numbers
   areas(I,1) = amin				
   amin = ARNEXT(amin)				!Next higher area number in the selection
*ENDDO

!nareas = 5 !for debugging
/BCSOPTION,,INCORE

ls=0						!Set current load step number

*DO,ia,1,nareas					!Loop over all interior faces

   ASEL,S,AREA,,areas(ia,1)			!Select current face
   TYPE,2					!Activate element type 2 (SURF154)
   AMESH,ALL					!Overlay face with type 2 element(s)

   ASUM						!Get statistics of the selected face
   *GET,fcent(ia,1),AREA,,CENT,X		!Get x location of face centroid
   *GET,fcent(ia,2),AREA,,CENT,Y		!Get y location of face centroid
   *GET,fcent(ia,3),AREA,,CENT,Z		!Get z location of face centroid

   *DO,dim,1,3					!Loop over all dimensions
      
      ls=ls+1					!Update load step number
      SFEDELE,ALL,ALL,ALL			!Delete all surface loads from elements

      ASEL,S,AREA,,areas(ia,1)			!Select current face
      ESLA,S					!Select associated surface element(s)
      *IF,dim,EQ,1,THEN
         SFE,ALL,5,PRES,,1,1,0,0		!Apply unit pressure [Pa] in global x direction
      *ELSEIF,dim,EQ,2
         SFE,ALL,5,PRES,,1,0,1,0		!Apply unit pressure [Pa] in global y direction
      *ELSEIF,dim,EQ,3
         SFE,ALL,5,PRES,,1,0,0,1		!Apply unit pressure [Pa] in global z direction
      *ENDIF

      ALLSEL					!Select all entities
      LSWRITE,ls				!Save load step
      						
   *ENDDO !Dimensions

   *CFOPEN,progress,tmp				!Save progress information to file
   *CFWRITE,Face,ia,nareas			
   *CFCLOS    

*ENDDO !Faces

!Solve
!*********************************************************************************************
FINISH						!Exit from processor
/SOLU						!Enter solver
ANTYPE,STATIC,NEW				!New static analysis
ALLSEL
OUTRES,ALL,NONE
OUTRES,NSOL,LAST 
LSSOLVE,1,ls					!Read and solve all load steps

!Post-processing
!*********************************************************************************************
FINISH						!Exit from processor
/POST1						!Enter general postprocessor

INRES,NSOL					!Read all items from results file
FILE,'TFM3D','rst','.'				!Load results file

NSEL,ALL					!Select all nodes
*GET,nnodes,NODE,,COUNT				!Get number of nodes
*DIM,nodes,ARRAY,nnodes,3 			!Array for nodal coordinates
*VGET,nodes(1,1),NODE,,LOC,X 			!Get nodal x coordinates
*VGET,nodes(1,2),NODE,,LOC,Y 			!Get nodal y coordinates
*VGET,nodes(1,3),NODE,,LOC,Z 			!Get nodal z coordinates

*DIM,uxnodes,ARRAY,nnodes,3			!Array for nodal x displacements
*DIM,uynodes,ARRAY,nnodes,3			!Array for nodal y displacements
*DIM,uznodes,ARRAY,nnodes,3			!Array for nodal z displacements
*DIM,uxqloc,ARRAY,nqloc,1 			!Array for x displacements at query points
*DIM,uyqloc,ARRAY,nqloc,1 			!Array for y displacements at query points
*DIM,uzqloc,ARRAY,nqloc,1 			!Array for z displacements at query points
*DIM,green,ARRAY,3*nqloc,3*nareas		!Array for discrete Green's function

*DO,lsi,1,ls					!Loop over all load steps
   
   SET,lsi					!Defines load step to be read from results file
   
   *VGET,uxnodes,NODE,,U,X			!Get nodal x displacements
   *VGET,uynodes,NODE,,U,Y			!Get nodal y displacements
   *VGET,uznodes,NODE,,U,Z			!Get nodal z displacements

   *MOPER,uxqloc,qloc,MAP,uxnodes,nodes		!Interpolate x displacement to query points
   *MOPER,uyqloc,qloc,MAP,uynodes,nodes		!Interpolate y displacement to query points
   *MOPER,uzqloc,qloc,MAP,uznodes,nodes		!Interpolate z displacement to query points

   *DO,I,1,nqloc				!Write displacements at query points to matrix
      green(3*(I-1)+1,lsi) = uxqloc(I,1)
      green(3*(I-1)+2,lsi) = uyqloc(I,1)
      green(3*(I-1)+3,lsi) = uzqloc(I,1)
   *ENDDO

*ENDDO !Load steps

*MWRITE,green(1,1),green,txt,,JIK,3*nareas,3*nqloc,1
(***nareas-1***(E16.9',')E16.9)			!Write discrete Green's function to file

*MWRITE,fcent(1,1),fcent,txt,,JIK,3,nareas,1
(2(E16.9',')E16.9)				!Write face centroid locations to file

FINISH						!Exit from processor