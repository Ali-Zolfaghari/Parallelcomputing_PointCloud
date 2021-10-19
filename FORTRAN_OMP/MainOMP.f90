
! ==== Parameters Definition 
Module param

Integer          ,Parameter    ::      dp = 8

Real(kind=dp)    ,Parameter    :: epsilon = 1.e-12_dp
Real(kind=dp)    ,Parameter    ::      pi = 4._dp*atan(1._dp)

End Module param

! ==== Matrices & Vectors Definition 
Module arrays

Use param        ,Only         : dp

Integer          ,Dimension( : )      ,Allocatable    :: Counter,CHECK

Real(kind=dp)    ,Dimension( :,: )    ,Allocatable    :: Point,Sphere,Distance

End Module arrays

! ==== Variables Definition 
Module var

Use param        ,Only         : dp

Integer          :: i,j,k,ii,jj,kk
Integer          :: NDNum,SphNum,NTID,CHUNK,TID,IN_Point,PT_PRT,PT_RES,IN_SPH,IOS
Real             :: TempData
Real(kind=dp)    :: ST_Time,FN_Time,MAX_DIS,MIN_DIS

End Module var

!==============================================================================

Program Main

Use param
Use arrays
Use var
Use omp_lib

Implicit None

CHARACTER  FC*15,FR*40

Call Readonly

ST_Time = OMP_GET_WTIME()

Call OMP_SET_NUM_THREADS(NTID)

!$OMP PARALLEL SHARED(Point,Sphere,NDNum,SphNum,Counter,CHUNK) PRIVATE(I,J,TID)
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
Do i = 1,NDNum
    Do j = 1,SphNum
        IF((((Point(i,1)-Sphere(j,1))**2+(Point(i,2)-Sphere(j,2))**2+  &
		&    (Point(i,3)-Sphere(j,3))**2)<=(Sphere(j,4))**2) .AND. Point(i,4)== -1.0_dp)Then
        Point(i,4)= 0.0_dp
        TID = OMP_GET_THREAD_NUM()
        Counter(TID+1)=Counter(TID+1)+1
        End IF
    End Do
End Do
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

IN_Point=0
Do i=1,NTID
IN_Point=IN_Point+Counter(i)
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

!$OMP PARALLEL SHARED(Point,Distance,NDNum,IN_Point) PRIVATE(I,J)
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
Do i = 1,NDNum
    Do j = 1,IN_Point
    Distance(i,j)=SQRT((Point(i,1)-Point(INT(Distance(NDNum+1,j)),1))**2  &
				 &    +(Point(i,2)-Point(INT(Distance(NDNum+1,j)),2))**2  &
				 &    +(Point(i,3)-Point(INT(Distance(NDNum+1,j)),3))**2)
    End Do
End Do
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

PT_PRT = INT(IN_Point/NTID)
PT_RES = MOD(IN_Point,NTID)

CHECK(1)=1
Do TID=0,NTID-1
IF( TID < PT_RES ) Then
	CHECK(TID+2) = CHECK(TID+1)+PT_PRT+1
Else
	CHECK(TID+2) = CHECK(TID+1)+PT_PRT
End IF
End Do

!$OMP PARALLEL SHARED(Distance,IN_SPH,NDNum,CHECK) PRIVATE(I,J,K,II,TID,MAX_DIS,MIN_DIS)
TID = OMP_GET_THREAD_NUM()
Do i = CHECK(TID+1),CHECK(TID+2)-1
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
!$OMP END PARALLEL 

Write(FC,'(I5)'),NTID
FR='OUTPUT'//FC//'.txt'
Open(40,File=FR)

Write(40,*),IN_Point
Do i=1,IN_Point
Write(40,"(F15.5,F15.5,F15.5,F15.5)"),Point(INT(Distance(NDNum+1,i)),1)  &
                                   & ,Point(INT(Distance(NDNum+1,i)),2)  &
		                           & ,Point(INT(Distance(NDNum+1,i)),3)  &
								   & ,Distance(NDNum+1+IN_SPH+1,i)
End Do

FN_Time = OMP_GET_WTIME()

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

Deallocate(Counter,CHECK)
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

Open(30,File='INPUT_OMP.txt')
Read(30,*),NTID
Allocate(Counter(NTID))
Allocate(CHECK(NTID+1))
Counter(1:NTID)= 0
Read(30,*),CHUNK
Read(30,*),IN_SPH
IN_SPH=IN_SPH+1
Close(30)

End SubRoutine ReadOnly

