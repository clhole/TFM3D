/BATCH
/FILNAME,TFM3D,0

!*************************************************************************************************
!**************************      FEA for traction calculation      *******************************
!*************************************************************************************************
!CL

!Modeling
!*************************************************************************************************
/PREP7 						!Enter preprocessor
/CWD,'***working directory***'			!Change working directory
/UNITS,uMKS					!Use MKS system of units

/INQUIRE,nlines,LINES,u,txt			!Get number of lines
nupos = nlines-1				!Get number of positions with known displacements
*DIM,utable,TABLE,nupos,6			!Define table to store displacement data
*TREAD,utable,u,txt				!Read displacement data from file to table

*DIM,upos,ARRAY,nupos,3				!Array for positions
*DIM,uxknown,ARRAY,nupos,1			!Array for known displacements in x
*DIM,uyknown,ARRAY,nupos,1			!Array for known displacements in y
*DIM,uzknown,ARRAY,nupos,1			!Array for known displacements in z

*VCOL,3,3					!Number of columns to copy
*MFUN,upos(1,1),COPY,utable(1,1)		!Copy bead positions into array
*VCOL,1,1					!Number of columns to copy
*MFUN,uxknown(1,1),COPY,utable(1,4)		!Copy bead displacements into array
*VCOL,1,1					!Number of columns to copy
*MFUN,uyknown(1,1),COPY,utable(1,5)		!Copy bead displacements into array
*VCOL,1,1					!Number of columns to copy
*MFUN,uzknown(1,1),COPY,utable(1,6)		!Copy bead displacements into array

*DIM,rtable,TABLE,3,2				!Define table for coordinate range
*TREAD,rtable,range,txt				!Read coordinate range from file to table
*DIM,rng,ARRAY,3,2				!Define array for coordinate range
*VCOL,3,3					!Number of columns to copy
*MFUN,rng(1,1),COPY,rtable(1,1)			!Copy coordinate range into array

~SATIN,'***modelname***','sat','***path***',SOLIDS,0
						!Load cell model				
BLOCK,rng(1,1),rng(1,2),rng(2,1),rng(2,2),rng(3,1),rng(3,2)
						!Create substrate bulk						
VSBV,2,1					!Subtract cell from the substrate bulk

!Group substrate edges for mesh definition
LSEL,S,LOC,X,rng(1,1)	
LSEL,A,LOC,X,rng(1,2)	
LSEL,A,LOC,Y,rng(2,1)
LSEL,A,LOC,Y,rng(2,2)
LSEL,A,LOC,Z,rng(3,1)
LSEL,A,LOC,Z,rng(3,2)
CM,edges,LINE

!Element and material definition
!*************************************************************************************************
ET,1,SOLID285					!4-node tetrahedral element
MP,EX,1,'***E***'   				!Young's modulus [Pa] of the substrate
MP,PRXY,1,'***nu***'				!Poisson's ratio of the substrate

!ET,2,SURF154					!Surface effect element
!KEYOPT,2,4,1					!...without midside nodes
!KEYOPT,2,11,2					!...and pressure vector on full area
!KEYOPT,2,12,0					!...applied regardless of element orientation

!Meshing
!*************************************************************************************************
TYPE,1						!Activate element type 1 (SOLID285)
MAT,1						!Activate material 1
LESIZE,edges,***size***				!Mesh spacing on substrate edges
ESIZE,***size***			!Define minimum mesh spacing

VMESH,3						!Mesh the substrate (volume 3)
!ASEL,S,LOC,X,-150,150				!Select all areas of the interior surface
!ASEL,R,LOC,Y,-150,150
!ASEL,R,LOC,Z,-150,150
!ASEL,S,AREA,,areas(ia,1)			!Select current face
!TYPE,2						!Activate element type 2 (SURF154)
!AMESH,ALL
!Apply loads and constraints


!*************************************************************************************************
!NSEL,ALL					!Select all nodes
ASEL,S,LOC,X,-20,20
ASEL,R,LOC,Y,-20,20
ASEL,R,LOC,Z,-20,20
!ESLA,S
NSLA,S,1
*GET,nnodes,NODE,,COUNT				!Get number of nodes
*DIM,nodes,ARRAY,nnodes,3 			!Define array for nodal coordinates
*VGET,nodes(1,1),NODE,,LOC,x 			!Nodal x coordinates
*VGET,nodes(1,2),NODE,,LOC,y 			!Nodal y coordinates
*VGET,nodes(1,3),NODE,,LOC,z 			!Nodal z coordinates

*DIM,uxnodes,ARRAY,nnodes,1 			!Array for nodal displacements in x
*MOPER,uxnodes,nodes,MAP,uxknown,upos,3,,0,100		!Interpolate x displacement to nodes

*DIM,uynodes,ARRAY,nnodes,1 			!Array for nodal displacements in y
*MOPER,uynodes,nodes,MAP,uyknown,upos,3,,0,100		!Interpolate y displacement to nodes

*DIM,uznodes,ARRAY,nnodes,1 			!Array for nodal displacements in z
*MOPER,uznodes,nodes,MAP,uzknown,upos,3,,0,100		!Interpolate z displacement to nodes

!Apply displacements to nodes
*DO,I,1,nnodes					
   D,I,UX,uxnodes(I)				!Displacement in x
   D,I,UY,uynodes(I)				!Displacement in y
   D,I,UZ,uznodes(I)				!Displacement in z
*ENDDO

!Apply DOF constraints (fixed substrate bottom)
!*********************************************************************************************
NSEL,S,LOC,Z,rng(3,1)		!Select nodes on the substrate bottom
D,ALL,UX,0					!Fix selected nodes in x direction
D,ALL,UY,0					!Fix selected nodes in y direction
D,ALL,UZ,0					!Fix selected nodes in z direction 

NLGEOM,OFF

!SOLVE
!*************************************************************************************************
FINISH						!Exit from processor
/SOLU						!Enter solver
ANTYPE,STATIC,NEW				!New static analysis
ALLSEL						!Select all entities
SOLVE						!Starts a solution

!Post-processing
!*************************************************************************************************
FINISH						!Exit from processor
/POST1						!Enter general postprocessor

*DIM,ns,ARRAY,nnodes,10				!Define array for stresses
*DO,I,1,nnodes
   ns(I,1)=I 					!Node number
   ns(I,2)=nodes(I,1)				!Nodal x coordinates
   ns(I,3)=nodes(I,2)				!Nodal y coordinates
   ns(I,4)=nodes(I,3)				!Nodal z coordinates
   *GET,ns(I,5),NODE,I,S,X			!Normal stress in x direction
   *GET,ns(I,6),NODE,I,S,Y			!Normal stress in y direction
   *GET,ns(I,7),NODE,I,S,Z			!Normal stress in z direction
   *GET,ns(I,8),NODE,I,S,XY			!Shear stress in xy plane
   *GET,ns(I,9),NODE,I,S,XZ			!Shear stress in xz plane
   *GET,ns(I,10),NODE,I,S,YZ			!Shear stress in yz plane
*ENDDO

*CFOPEN,stress,txt 				!Create file
*CFWRITE,Node,N,x,y,z,Tx,Ty,Tz,Txy,Txz,Tyz
*DO,I,1,nnodes					!Write data to file
   *CFWRITE,Node,ns(I,1),ns(I,2),ns(I,3),ns(I,4),ns(I,5),ns(I,6),ns(I,7),ns(I,8),ns(I,9),ns(I,10)
*ENDDO
*CFCLOS

FINISH						!Exit from processor