/BATCH
/FILNAME,sat2ansTest,0


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
BLOCK,-200,200,-200,200,-200,200		!Create substrate bulk
VSBV,2,1					!Subtract cell from the substrate bulk

!Create output file
*CFOPEN,sat2ansTest,txt
*CFWRITE,STATUS,1
*CFCLOS

FINISH