
! ==== Parameters Definition 

Module param

Integer          ,Parameter    ::      dp = 8
Integer          ,Parameter    ::  MASTER = 0
Integer          ,Parameter    ::     M2W = 1
Integer          ,Parameter    ::     W2M = 2

Real(kind=dp)    ,Parameter    :: epsilon = 1.e-12_dp
Real(kind=dp)    ,Parameter    ::      pi = 4._dp*atan(1._dp)

End Module param

! ==== Matrices & Vectors Definition 
Module arrays

Use param        ,Only         : dp

Real(kind=dp)    ,Dimension( : )      ,Allocatable    :: IN_P
Real(kind=dp)    ,Dimension( :,: )    ,Allocatable    :: Point,Sphere,Distance

End Module arrays

! ==== Variables Definition 
Module var

Use param        ,Only         : dp

Integer          :: i,j,k,ii,jj,kk
Integer          :: ierr,INPT_Type,ArryID,Dest,Src,COL
Integer          :: RID,NTP,Counter
Integer          :: NDNum,SphNum,IN_Point,IN_SPH,PT_PRT,PT_RES,IOS
Real             :: TempData

Real(kind=dp)    :: ST_Time,FN_Time,MAX_DIS,MIN_DIS

End Module var
!==============================================================================

Program Main

Use param
Use arrays
Use var
Use mpi

Implicit None

CHARACTER  FC*15,FR*40
Integer    status(MPI_STATUS_SIZE)

Call MPI_INIT( ierr )

Call MPI_COMM_RANK( MPI_COMM_WORLD , RID , ierr )
Call MPI_COMM_SIZE( MPI_COMM_WORLD , NTP , ierr )

IF(RID == MASTER)Then

ST_Time = MPI_Wtime()

Call Readonly

PT_PRT = INT(NDNum/(NTP-1))
PT_RES = MOD(NDNum,(NTP-1))
ArryID = 1

Do Dest = 1,NTP-1
    IF(Dest <= PT_RES) Then
    COL = PT_PRT + 1
    Else
    COL = PT_PRT
    End IF

Call MPI_SEND( NDNum           , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( SphNum          , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( ArryID          , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( COL             , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( Sphere          , 4*SphNum , MPI_DOUBLE_PRECISION , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( Point(1,ArryID) , 4*COL    , MPI_DOUBLE_PRECISION , Dest , M2W , MPI_COMM_WORLD , ierr )

ArryID = ArryID + COL
End Do


IN_Point = 0
Do Src = 1, NTP-1

Call MPI_RECV( ArryID          , 1     , MPI_INTEGER          , Src , W2M , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( COL             , 1     , MPI_INTEGER          , Src , W2M , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( Counter         , 1     , MPI_INTEGER          , Src , W2M , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( Point(1,ArryID) , 4*COL , MPI_DOUBLE_PRECISION , Src , W2M , MPI_COMM_WORLD , status , ierr )

IN_Point = IN_Point+Counter
End Do

Allocate(Distance(NDNum+1+IN_SPH+1,IN_Point))
Allocate(IN_P(IN_Point))

j=1
Do i = 1,NDNum
    IF(Point(4,i)== 0.0_dp)Then
    Distance(NDNum+1,j)=i
	IN_P(j)=i
    j=j+1
    End IF
End Do

PT_PRT = INT(IN_Point/(NTP-1))
PT_RES = MOD(IN_Point,(NTP-1))
ArryID = 1

Do Dest = 1,NTP-1
    IF(Dest <= PT_RES) Then
    COL = PT_PRT + 1
    Else
    COL = PT_PRT
    End IF

Call MPI_SEND( IN_SPH             , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( IN_Point           , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( ArryID             , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( COL                , 1        , MPI_INTEGER          , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( Point              , 4*NDNum  , MPI_DOUBLE_PRECISION , Dest , M2W , MPI_COMM_WORLD , ierr )
Call MPI_SEND( IN_P(ArryID)       , COL      , MPI_DOUBLE_PRECISION , Dest , M2W , MPI_COMM_WORLD , ierr )

ArryID = ArryID + COL
End Do

Do Src = 1, NTP-1 

Call MPI_RECV( ArryID             , 1                      , MPI_INTEGER          , Src , W2M , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( COL                , 1                      , MPI_INTEGER          , Src , W2M , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( Distance(1,ArryID) , COL*(NDNum+1+IN_SPH+1) , MPI_DOUBLE_PRECISION , Src , W2M , MPI_COMM_WORLD , status , ierr )

End Do

Write(FC,'(I5)'),NTP
FR='OUTPUT'//FC//'.txt'
Open(40,File=FR)

Write(40,*),IN_Point
Do i=1,IN_Point
Write(40,"(F15.5,F15.5,F15.5,F15.5)"),Point(1,INT(Distance(NDNum+1,i)))  &
                                   & ,Point(2,INT(Distance(NDNum+1,i)))  &
		                           & ,Point(3,INT(Distance(NDNum+1,i)))  &
								   & ,Distance(NDNum+1+IN_SPH+1,i)
End Do

FN_Time = MPI_Wtime()

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

Deallocate(Point,Sphere,Distance,IN_P)

Else

Call MPI_RECV( NDNum  , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( SphNum , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( ArryID , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( COL    , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )

Allocate(Point(4,NDNum))
Allocate(Sphere(4,SphNum))

Call MPI_RECV( Sphere , 4*SphNum , MPI_DOUBLE_PRECISION , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( Point  , 4*COL    , MPI_DOUBLE_PRECISION , MASTER , M2W , MPI_COMM_WORLD , status , ierr )

Counter = 0
Do i = 1,COL
    Do j = 1,SphNum
        IF((((Point(1,i)-Sphere(1,j))**2+(Point(2,i)-Sphere(2,j))**2  &
		&   +(Point(3,i)-Sphere(3,j))**2)<=(Sphere(4,j))**2) .AND. Point(4,i)== -1.0_dp)Then
        Point(4,i)= 0.0_dp
        Counter = Counter+1
        End IF
    End Do
End Do

Call MPI_SEND( ArryID  , 1     , MPI_INTEGER          , MASTER , W2M , MPI_COMM_WORLD , ierr )
Call MPI_SEND( COL     , 1     , MPI_INTEGER          , MASTER , W2M , MPI_COMM_WORLD , ierr )
Call MPI_SEND( Counter , 1     , MPI_INTEGER          , MASTER , W2M , MPI_COMM_WORLD , ierr )
Call MPI_SEND( Point   , 4*COL , MPI_DOUBLE_PRECISION , MASTER , W2M , MPI_COMM_WORLD , ierr )


Call MPI_RECV( IN_SPH   , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( IN_Point , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( ArryID   , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( COL      , 1        , MPI_INTEGER          , MASTER , M2W , MPI_COMM_WORLD , status , ierr )
Call MPI_RECV( Point    , 4*NDNum  , MPI_DOUBLE_PRECISION , MASTER , M2W , MPI_COMM_WORLD , status , ierr )

Allocate(Distance(NDNum+1+IN_SPH+1,COL))
Allocate(IN_P(COL))

Call MPI_RECV( IN_P     , COL      , MPI_DOUBLE_PRECISION , MASTER , M2W , MPI_COMM_WORLD , status , ierr )

Distance(NDNum+1,1:COL) = IN_P(1:COL)

Do i = 1,COL
    Do j = 1,NDNum
	    Distance(j,i)=SQRT((Point(1,j)-Point(1,INT(Distance(NDNum+1,i))))**2  &
				     &    +(Point(2,j)-Point(2,INT(Distance(NDNum+1,i))))**2  &
				     &    +(Point(3,j)-Point(3,INT(Distance(NDNum+1,i))))**2)
    End Do
End Do

Do i = 1,COL
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

Call MPI_SEND( ArryID   , 1                      , MPI_INTEGER          , MASTER , W2M , MPI_COMM_WORLD , ierr )
Call MPI_SEND( COL      , 1                      , MPI_INTEGER          , MASTER , W2M , MPI_COMM_WORLD , ierr )
Call MPI_SEND( Distance , COL*(NDNum+1+IN_SPH+1) , MPI_DOUBLE_PRECISION , MASTER , W2M , MPI_COMM_WORLD , ierr )

Deallocate(Point,Sphere,Distance,IN_P)

End IF

Call MPI_FINALIZE(ierr)


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
Allocate(Point(4,NDNum))
Point(4,1:NDNum)= -1.0_dp

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
Allocate(Sphere(4,SphNum))


Open(10,File='INPUT_Pxyz.txt')
Do i=1,NDNum
Read(10,*),Point(1,i),Point(2,i),Point(3,i)
End Do
Close(10)

Open(20,File='INPUT_Sxyz.txt')
Do i=1,SphNum
Read(20,*),Sphere(1,i),Sphere(2,i),Sphere(3,i),Sphere(4,i)
End Do
Close(20)

Open(30,File='INPUT_MPI.txt')
Read(30,*),IN_SPH
IN_SPH=IN_SPH+1
Close(30)

End SubRoutine ReadOnly

