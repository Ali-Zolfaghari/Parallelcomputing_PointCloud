
! ==== Parameters Definition 
Module param

Integer          ,Parameter    ::      dp = 8

Real(kind=dp)    ,Parameter    :: epsilon = 1.e-12_dp
Real(kind=dp)    ,Parameter    ::      pi = 4._dp*atan(1._dp)

End Module param

! ==== Matrices & Vectors Definition 
Module arrays

Use param        ,Only         : dp

Real(kind=dp)    ,Dimension( :,: )    ,Allocatable    :: Point,Sphere,Distance

End Module arrays

! ==== Variables Definition 
Module var

Use param        ,Only         : dp

Integer          :: i,j,k,ii,jj,kk
Integer          :: NDNum,SphNum,IN_Point,IN_SPH,IOS
Real             :: TempData
Real(kind=dp)    :: ST_Time,FN_Time,MAX_DIS,MIN_DIS

End Module var

!==============================================================================

Program Main

Use param
Use arrays
Use var

Implicit None

Call Readonly

Call CPU_TIME(ST_Time)

IN_Point=0
Do i = 1,NDNum
    Do j = 1,SphNum
        IF((((Point(i,1)-Sphere(j,1))**2+(Point(i,2)-Sphere(j,2))**2+  &
		&    (Point(i,3)-Sphere(j,3))**2)<=(Sphere(j,4))**2) .AND. Point(i,4)== -1.0_dp)Then
        Point(i,4)= 0.0_dp
        IN_Point=IN_Point+1
        End IF
    End Do
End Do

Allocate(Distance(NDNum+1+IN_SPH+1,IN_Point))
Distance = 0.0_dp

j=1
Do i = 1,NDNum
    IF(Point(i,4)== 0.0_dp)Then
    Distance(NDNum+1,j)=i
    j=j+1
    End IF
End Do

Do i = 1,NDNum
    Do j = 1,IN_Point
    Distance(i,j)=SQRT((Point(i,1)-Point(INT(Distance(NDNum+1,j)),1))**2  &
				 &    +(Point(i,2)-Point(INT(Distance(NDNum+1,j)),2))**2  &
				 &    +(Point(i,3)-Point(INT(Distance(NDNum+1,j)),3))**2)
    End Do
End Do

Do i = 1,IN_Point
MAX_DIS = MAXVAL(Distance(1:NDNum,i))
    Do j = 1,IN_SPH
    MIN_DIS = MINVAL(Distance(1:NDNum,i))
	Distance(NDNum+1+IN_SPH+1,i)=Distance(NDNum+1+IN_SPH+1,i)+MIN_DIS
        Do k = 1,NDNum
            IF(Distance(k,i)==MIN_DIS)Then
            Distance(NDNum+1+j,i)=k
            Distance(k,i)=MAX_DIS
            EXIT
            End IF
        End Do
    End Do
Distance(NDNum+1+IN_SPH+1,i)=Distance(NDNum+1+IN_SPH+1,i)/(IN_SPH-1)
End Do

Open(40,File='OUTPUT_SERIAL.txt')

Write(40,*),IN_Point
Do i=1,IN_Point
Write(40,"(F15.5,F15.5,F15.5,F15.5)"),Point(INT(Distance(NDNum+1,i)),1)  &
                                   & ,Point(INT(Distance(NDNum+1,i)),2)  &
		                           & ,Point(INT(Distance(NDNum+1,i)),3)  &
		                           & ,Distance(NDNum+1+IN_SPH+1,i)
End Do


Call CPU_TIME(FN_Time)

Write(40,*),'    FINAL CALCULATION TIME  : ',FN_Time-ST_Time
Write(40,*),'|===========================================================================|'
Write(40,*),'|===========================================================================|'
Write(40,*),'|                                                                           |'
Write(40,*),'|    ***** ***** *   * ***** *     ***** ***** ***** *****    ***** *   *   |'
Write(40,*),'|     *  * *     *   * *     *     *   * *   * *      *  *    *   * *   *   |'
Write(40,*),'|     *  * ***** *   * ***** *     *   * ***** *****  *  *    ***** *****   |'
Write(40,*),'|     *  * *      * *  *     *     *   * *     *      *  *    *   *     *   |'
Write(40,*),'|    ***** *****   *   ***** ***** ***** *     ***** *****    ***** *****   |'
Write(40,*),'|                                                                           |'
Write(40,*),'|    *****    ***** ***** *     ***** ***** ***** *   * ***** ***** *****   |'
Write(40,*),'|    *   *       *  *   * *     *     *   * *     *   * *   * *   *   *     |'
Write(40,*),'|    *****      *   *   * *     ***** ***** ***** ***** ***** *****   *     |'
Write(40,*),'|    *   * ##  *    *   * *     *     *   * *   * *   * *   * * *     *     |'
Write(40,*),'|    *   * ## ***** ***** ***** *     *   * ***** *   * *   * *   * *****   |'
Write(40,*),'|                                                                           |'
Write(40,*),'|===========================================================================|'
Write(40,*),'|===========================================================================|' 
Close(40)

Print*,'    FINAL CALCULATION TIME  : ',FN_Time-ST_Time

Deallocate(Point,Sphere,Distance)

End Program Main

!==============================================================================

SubRoutine ReadOnly

Use param
Use arrays
Use var

Implicit None


NDNum=0
IOS=0
Open(10,File='INPUT_Pxyz.txt')
Do While (IOS==0)
Read (10,*,iostat=IOS),TempData    
NDNum=NDNum+1
End Do    
Close(10)
NDNum=NDNum-1
Print*,'    Node Number  : ',NDNum
Allocate(Point(NDNum,4))
Point(1:NDNum,4)= -1.0_dp

SphNum=0
IOS=0
Open(20,File='INPUT_Sxyz.txt')
Do While (IOS==0)
Read (20,*,iostat=IOS),TempData    
SphNum=SphNum+1
End Do    
Close(20)
SphNum=SphNum-1
Print*,'    Sphere Number  : ',SphNum
Allocate(Sphere(SphNum,4))


Open(10,File='INPUT_Pxyz.txt')
Do i=1,NDNum
Read(10,*),Point(i,1),Point(i,2),Point(i,3)
End Do
Close(10)

Open(20,File='INPUT_Sxyz.txt')
Do i=1,SphNum
Read(20,*),Sphere(i,1),Sphere(i,2),Sphere(i,3),Sphere(i,4)
End Do
Close(20)

Open(30,File='INPUT_MPI.txt')
Read(30,*),IN_SPH
IN_SPH=IN_SPH+1
Close(30)

End SubRoutine ReadOnly

