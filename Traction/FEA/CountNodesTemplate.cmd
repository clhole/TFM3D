/BATCH
/FILNAME,CountNodes,0


!*******************************************************************************************
!**************      Check import and boolean operation of .sat file      ******************
!*******************************************************************************************
!CL


!Modeling
!*******************************************************************************************
/PREP7 						!Enter preprocessor

/CWD,'***working directory***'			!Change working directory

~SATIN,'***modelname***','sat','***path***',SOLIDS,0
						!Load cell model
BLOCK,'***xmin,xmax,ymin,ymax,zmin,zmax***'	!Create substrate bulk
VGEN,,1,,,'***dx,dy,dz***',,,1			!Shift cell to substrate center
VSBV,2,1					!Subtract cell from the substrate bulk
					
!Group substrate edges for mesh definition
LSEL,S,LOC,X,0		
LSEL,A,LOC,Y,0
LSEL,A,LOC,Y,8
LSEL,A,LOC,Z,0
LSEL,A,LOC,Z,8
CM,edges,LINE

!Meshing
ET,1,SOLID285					!Element type
LESIZE,edges,***size***				!Mesh spacing on substrate edges
ESIZE,***size***				!Define minimum mesh spacing
VMESH,3						!Mesh the substrate (volume 3)

NSEL,ALL					!Select all nodes
*GET,nnodes,NODE,,COUNT				!Get number of nodes

!Create output file for total number of nodes
*CFOPEN,CountNodes,txt
*CFWRITE,Nodes,nnodes
*CFCLOS

FINISH