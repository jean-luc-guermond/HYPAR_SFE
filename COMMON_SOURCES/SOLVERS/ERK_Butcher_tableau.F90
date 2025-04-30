MODULE Butcher_tableau
  IMPLICIT NONE
  PRIVATE
  TYPE, PUBLIC :: BT
     INTEGER :: s=2, v=1, sv=21
     REAL(KIND=8), DIMENSION(:,:), POINTER :: A, inc_A, MatRK
     REAL(KIND=8), DIMENSION(:),   POINTER :: C, inc_C
     INTEGER,      DIMENSION(:),   POINTER :: lp_of_l
   CONTAINS
     PROCEDURE :: init => init_bt
  END TYPE BT
CONTAINS
  SUBROUTINE init_bt(this, sv)
    IMPLICIT NONE
    CLASS(BT), INTENT(inout) :: this
    INTEGER :: asv, sv
    REAL(KIND=8) :: gamma
    this%sv = sv
    asv= ABS(this%sv)
    this%s = (asv-MODULO(asv,10))/10
    IF (ASSOCIATED(this%C)) THEN
       NULLIFY   (this%C,this%inc_C,this%lp_of_l,this%A,this%inc_A,This%MatRK)
    END IF
    ALLOCATE(this%A(this%s+1,this%s))
    ALLOCATE(this%inc_A(this%s+1,this%s))
    ALLOCATE(this%MatRK(this%s+1,this%s))
    ALLOCATE(this%C(this%s+1))
    ALLOCATE(this%inc_C(this%s+1))
    ALLOCATE(this%lp_of_l(this%s+1))

    this%A =0.d0
    this%inc_A =0.d0
    this%C = 0.d0
    this%inc_C = 0.d0
    this%lp_of_l = 0
    this%C(this%s+1)=1.d0
    SELECT CASE(asv)

    CASE(11) !===Equi RK2, ERK(2,2;1) two-stage, 2nd order
       this%C=(/0.d0,1.d0/)
       this%lp_of_l=(/1,1/)
       this%A(2,1)=this%C(2)
    CASE(21) !===Equi RK2, ERK(2,2;1) two-stage, 2nd order
       this%C=(/0.d0,0.5d0,1.d0/)
       this%lp_of_l=(/1,1,2/)
       this%A(2,1)=this%C(2)
       this%A(3,2)=this%C(3)
    CASE(22) !===SSP RK2, ERK(2,2;1/2) two-stage, 2nd order
       this%C=(/0.d0,1.d0,1.d0/)
       this%lp_of_l=(/1,1,2/)
       this%A(2,1)=this%C(2)
       this%A(3,1)=0.5d0
       this%A(3,2)=0.5d0
    CASE(31) !===Equi RK3 ERK(3,3;1) three-stage, 3nd order
       this%C=(/0.d0,1.d0/3,2.d0/3,1.d0/)
       this%lp_of_l=(/1,1,2,3/)
       this%A(2,1)=this%C(2)
       this%A(3,2)=this%C(3)
       this%A(4,1)=1.d0/4
       this%A(4,3)=3.d0/4
    CASE (32) !===SSP SSPRK3 ERK(3,3;1/3) three-stage, 3nd order
       this%C=(/0.d0,1.d0,0.5d0,1.d0/)
       this%lp_of_l=(/1,1,1,2/)
       this%A(2,1)=this%C(2)
       this%A(3,1)=0.25d0
       this%A(3,2)=0.25d0
       this%A(4,1)=1.d0/6
       this%A(4,2)=1.d0/6
       this%A(4,3)=2.d0/3
    CASE (33) !===Crouzeix-Norsett three-stage, 3nd order
       gamma = 1.d0/2*(1.d0 + 1.d0/SQRT(3.d0))
       this%C=(/0.d0,gamma,1.d0-gamma,1.d0/)
       this%lp_of_l=(/1,1,1,2/)
       this%A(2,1)=this%C(2)
       this%A(3,1)=gamma-1.d0
       this%A(3,2)=2*(1.d0-gamma)
       this%A(4,1)=0.d0
       this%A(4,2)=1.d0/2
       this%A(4,3)=1.d0/2
    CASE (41) !===Equi RK4  ERK(4,3,1) (seems to be only third-order)
       this%C=(/0.d0,1.d0/4,0.5d0,3.d0/4,1.d0/)
       this%lp_of_l=(/1,1,2,3,4/)
       this%A(2,1)=this%C(2)
       this%A(3,2)=this%C(3)
       this%A(4,2)=1.d0/4
       this%A(4,3)=1.d0/2
       this%A(5,2)=2.d0/3
       this%A(5,3)=-1.d0/3
       this%A(5,4)=2.d0/3
    CASE (42) !===popular RK4 ERK(4,4;1/2)
       this%C=(/0.d0,1.d0/2,1.d0/2,1.d0,1.d0/)
       this%lp_of_l=(/1,1,2,3,4/)
       this%A(2,1)=this%C(2)
       this%A(3,2)=this%C(3)
       this%A(4,3)=this%C(4)
       this%A(5,1)=1.d0/6
       this%A(5,2)=2.d0/6
       this%A(5,3)=2.d0/6
       this%A(5,4)=1.d0/6
    CASE (43) !===3/8 rule RK4 ER4 ERK(4,4;3/4)
       this%C=(/0.d0,1.d0/3,2.d0/3,1.d0,1.d0/)
       this%lp_of_l=(/1,1,2,3,4/)
       this%A(2,1)=this%C(2)
       this%A(3,1)=-1.d0/3
       this%A(3,2)=1.d0
       this%A(4,1)=1.d0
       this%A(4,2)=-1.d0
       this%A(4,3)=1.d0
       this%A(5,1)=1.d0/8
       this%A(5,2)=3.d0/8
       this%A(5,3)=3.d0/8
       this%A(5,4)=1.d0/8
    CASE (44) !===4 stages, third-order (Optimal, AEJLG, ERK(4,3;1))
       this%C=(/0.d0,1.d0/4,2.d0/4,3.d0/4,1.d0/)
       this%lp_of_l=(/1,1,2,3,4/)
       this%A(2,:)=(/ 0.25d0,               0.d0,                0.d0,                  0.d0/)
       this%A(3,:)=(/-0.0221890115021647d0, 0.5221890115021647d0, 0.d0,                 0.d0/) 
       this%A(4,:)=(/ 0.d0,                 0.2515028538497361d0, 0.4984971461502638d0, 0.d0/)
       !!this%A(3,:)=(/0.2500000000000000d0, 0.2500000000000000d0, 0.d0, 0.d0/) 
       !!this%A(4,:)=(/ 0.0000000000000000d0, 0.3595354020227154d0, 0.3904645979772846d0, 0.d0/) 
       this%A(5,:)=(/ 0.0264037560338726d0, 0.5874553985650487d0,-0.2541220652317153d0, 0.6402629106327939d0/)
    CASE (51) !===SSPRK(5,4;0.51).
       this%C(2)= 0.39175222657188905833d0
       this%C(3)= 0.58607968931154123111d0
       this%C(4)= 0.47454236312139913362d0
       this%C(5)= 0.93501063096765159845d0
       this%lp_of_l=(/1,1,2,2,3,5/)
       this%A(2,1)= this%C(2)
       this%A(3,1)= 0.21766909626116921036d0
       this%A(3,2)= 0.36841059305037202075d0
       this%A(4,1)= 0.08269208665781075441d0
       this%A(4,2)= 0.13995850219189573938d0
       this%A(4,3)= 0.25189177427169263984d0
       this%A(5,1)= 0.06796628363711496324d0
       this%A(5,2)= 0.11503469850463199467d0
       this%A(5,3)= 0.20703489859738471851d0
       this%A(5,4)= 0.54497475022851992204d0
       this%A(6,1)= 0.14681187608478644956d0
       this%A(6,2)= 0.24848290944497614757d0
       this%A(6,3)= 0.10425883033198029567d0
       this%A(6,4)= 0.27443890090134945681d0
       this%A(6,5)= 0.22600748323690765039d0
   CASE(52) !=== Five stages, fourth-order, optimal, (AEJLG, ERK(5,4;1))
       this%C=(/0.d0,1.d0/5,2.d0/5,3.d0/5,4.d0/5,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5/)
       this%A(2,:)=(/ 0.2d0,                 0.d0,                                         0.d0,                  0.d0/)
       this%A(3,:)=(/ 0.26075582269554909d0, 0.13924417730445096d0, 0.d0,                  0.d0,                  0.d0/)
       this%A(4,:)=(/-0.25856517872570289d0, 0.91136274166280729d0,-0.05279756293710430d0, 0.d0,                  0.d0/)
       this%A(5,:)=(/ 0.21623276431503774d0, 0.51534223099602405d0,-0.81662794199265554d0, 0.88505294668159373d0, 0.d0/)
       this%A(6,:)=(/-0.10511678454691901d0, 0.87880047152100838d0,-0.58903404061484477d0, 0.46213380485434047d0, &
            0.35321654878641495d0/)      
    CASE (61) !===Equi 1/5, ERK(6,5;5/6)
       this%C=(/0.d0,1.d0/5,2.d0/5,3.d0/5,4.d0/5,1.d0,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5,6/)
       this%A(2,1)=this%C(2)
       this%A(3,2)=this%C(3)
       this%A(4,1)=3.d0/20
       this%A(4,3)=9.d0/20
       this%A(5,1)= 4.d0/5
       this%A(5,2)=-8.d0/5
       this%A(5,3)= 8.d0/5
       this%A(6,1)=-71.d0/4
       this%A(6,2)= 40.d0
       this%A(6,3)=-75.d0/4
       this%A(6,4)=-10.d0
       this%A(6,5)= 15.d0/2
       this%A(7,1)=17.d0/144
       this%A(7,3)=25.d0/36
       this%A(7,4)=-25.d0/72
       this%A(7,5)=25.d0/48
       this%A(7,6)=1.d0/72
    CASE (62) !===Equi 1/6 Six stages, fourth-order (AEJLG, ERK(6,4;1))
       this%C=(/0.d0,1.d0/6,2.d0/6,3.d0/6,4.d0/6,5.d0/6,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5,6,7/)
       !===This is nonlinear 4th order and linear 5th order and IMEX compatibble 
       this%A(2,:) = (/ 0.1666666666666667d0, 0.d0,                 0.d0,                 0.d0,0.d0,0.d0/)
       this%A(3,:) = (/-0.4447518666866271d0, 0.7780852000199604d0, 0.d0,                 0.d0,0.d0,0.d0/)
       this%A(4,:) = (/ 0.0893971199002535d0, 0.1913734465774781d0, 0.2192294335222684d0, 0.d0,0.d0,0.d0/)
       this%A(5,:) = (/ 0.0635170175925083d0, 0.1428758587504731d0, 0.1359933602040204d0, 0.3242804301196647d0, 0.d0, 0.d0/)
       this%A(6,:) = (/ 0.0727304753901270d0, 0.2698992458411764d0,-0.0619049508228309d0, 0.2187862524098531d0, &
            0.3338223105150078d0, 0.d0/)
       this%A(7,:) = (/ 0.0830000000000030d0, 0.1349999999999911d0, 0.1300000000000078d0, 0.4700000000000000d0, &
            -0.2850000000000030d0, 0.4670000000000011d0/) 
    CASE (63) !===Butcher 1964, ERK(6,5;2/3)
       this%C=(/0.d0,1.d0/4,1.d0/4,1.d0/2,3.d0/4,1.d0,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5,6/)
       this%A(2,1)=this%C(2)
       this%A(3,1)=1.d0/8
       this%A(3,2)=1.d0/8
       this%A(4,2)=-1.d0/2
       this%A(4,3)= 1.d0
       this%A(5,1)=3.d0/16
       this%A(5,4)=9.d0/16
       this%A(6,1)= -3.d0/7
       this%A(6,2)=  2.d0/7
       this%A(6,3)= 12.d0/7
       this%A(6,4)=-12.d0/7
       this%A(6,5)=  8.d0/7
       this%A(7,1)= 7.d0/90
       this%A(7,3)=32.d0/90
       this%A(7,4)=12.d0/90
       this%A(7,5)=32.d0/90
       this%A(7,6)= 7.d0/90
    CASE(64) !====ERK(6,5;1) bug. Only fourth order. AEJLG
       this%C=(/0.d0,1.d0/6,2.d0/6,3.d0/6,4.d0/6,5.d0/6,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5,6/)
       this%lp_of_l=(/1,1,2,3,4,5,6/)
       this%A(2,:) = (/ 0.1666666666666667d0, 0.d0,                 0.d0,                 0.d0,0.d0,0.d0/)
       this%A(3,:) = (/-0.0005041509245550d0, 0.3338374842578883d0, 0.d0,                 0.d0,0.d0,0.d0/)
       this%A(4,:) = (/ 0.4949827381621411d0,-0.7382930557116373d0, 0.7433103175494962d0, 0.d0,0.d0,0.d0/)
       this%A(5,:) = (/ 0.1145718006265044d0, 0.6464301143951824d0,-0.9682665500489782d0, 0.8739313016939579d0, &
            0.d0, 0.d0/)
       this%A(6,:) = (/ 0.0195678784369044d0, 0.2638733780720444d0, 0.3395907234578278d0,-0.2989195535788929d0, &
            0.5092209069454497d0, 0.d0/)
       this%A(7,:) = (/ 0.1098931622638753d0, 0.0005341886806206d0, 0.3989316226387609d0, 0.2010683773612388d0,&
            -0.1505341886806198d0, 0.4401068377361241d0/)
    CASE(71) !====ERK(7,5;1) nonlinear 5th order, AEJLG
       this%C=(/0.d0,1.d0/7,2.d0/7,3.d0/7,4.d0/7,5.d0/7,6.d0/7,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5,6,7/)
       this%A(2,:) = (/ 0.1428571428571428d0, 0.d0,                 0.d0,                 0.d0,                 0.d0,0.d0,0.d0/)
       this%A(3,:) = (/ 0.0107112392440190d0, 0.2750030464702667d0, 0.d0,                 0.d0,                 0.d0,0.d0,0.d0/)
       this%A(4,:) = (/ 0.4812641640978049d0,-0.9634955610241845d0, 0.9108028254978082d0, 0.d0,                 0.d0,0.d0,0.d0/)
       this%A(5,:) = (/ 0.3718168921588643d0,-0.5615016072645083d0, 0.5590150320678782d0, 0.2020982544663372d0, 0.d0,0.d0,0.d0/)
       this%A(6,:) = (/ 0.2210152091353196d0, 0.3526985345185882d0,-0.8940286416538530d0, 0.8097519357353253d0, &
             0.2248486765503341d0, 0.d0,                 0.d0/)
       this%A(7,:) = (/ 0.2038005573304989d0,-0.4759394836773749d0, 1.0938423462713556d0,-0.2853403360393202d0, &
            -0.1249739792585170d0, 0.4457537525162146d0, 0.d0/)
       this%A(8,:) = (/ 0.0979996468518429d0,-0.0044680013474942d0, 0.3592897484042626d0, 0.0225280828210270d0, &
             0.2680292384753062d0,-0.1064595934043306d0, 0.3630808781993861d0/)    
    CASE default
       WRITE(*,*) ' BUG in init_bt: wrong sv'
       STOP 'fail'
    END SELECT
    IF (this%sv>0) THEN
       this%lp_of_l=1
       this%MatRK=>this%A
    ELSE
       CALL inc_bt(this)
       this%MatRK=>this%inc_A
    END IF
  END SUBROUTINE init_bt

  SUBROUTINE inc_bt(this)
    IMPLICIT NONE
    CLASS(BT), INTENT(inout) :: this
    INTEGER :: l
    this%inc_A = 0.d0
    this%inc_C = 0.d0
    this%inc_C = this%C - this%C(this%lp_of_l) !===c(l)-c(l')
    DO l = 2, this%s+1
       this%inc_A(l,:) = this%A(l,:)-this%A(this%lp_of_l(l),:) !===a(l,:) - a(l',:)
    END DO
  END SUBROUTINE inc_bt
END MODULE Butcher_tableau
