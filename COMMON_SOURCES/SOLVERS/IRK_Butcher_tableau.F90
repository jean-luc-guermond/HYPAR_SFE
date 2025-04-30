MODULE Implicit_Butcher_tableau
  IMPLICIT NONE
  PRIVATE
  TYPE, PUBLIC :: IBT
     INTEGER :: s=2, v=1, sv=21
     REAL(KIND=8), DIMENSION(:,:), POINTER :: A, inc_A, MatRK
     REAL(KIND=8), DIMENSION(:),   POINTER :: C, inc_C
     INTEGER,      DIMENSION(:),   POINTER :: lp_of_l
   CONTAINS
     PROCEDURE :: init => init_bt
  END TYPE IBT
CONTAINS
  SUBROUTINE init_bt(this, sv)
    IMPLICIT NONE
    CLASS(IBT), INTENT(inout) :: this
    INTEGER :: asv, sv
    REAL(KIND=8) :: gamma
    this%sv =  sv
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

    CASE(21) !===Equi IRK2, ERK(2,2;1) two-stage, 2nd order
       this%C=(/0.d0,0.5d0,1.d0/)
       this%lp_of_l=(/1,1,2/)
       this%A(2,:)=(/0.d0,0.5d0/)
       this%A(3,:)=(/0.d0,1.d0/)
      
    CASE(22) !===SSP IRK2, ERK(2,2;1/2) two-stage, 2nd order
       this%C=(/0.d0,1.d0,1.d0/)
       this%lp_of_l=(/1,1,2/)
       this%A(2,:)=(/0.5d0,0.5d0/)
       this%A(3,:)=(/0.5d0,0.5d0/)

    CASE(31) !===Equi IRK3 ERK(3,3;1) three-stage, 3nd order
       gamma=(1.d0+1.d0/sqrt(3.d0))/2
       this%C=(/0.d0,1.d0/3,2.d0/3,1.d0/)
       this%lp_of_l=(/1,1,2,3/)
       this%A(2,:)=(/1.d0/3-gamma, gamma, 0.d0 /)
       this%A(3,:)=(/gamma, 2.d0/3-2.d0*gamma, gamma /)
       this%A(4,:)=(/0.25d0,0.d0,0.75d0/)
       
    CASE (33) !===Crouzeix-Norsett three-stage, 3nd order
       gamma = 1.d0/2*(1.d0 + 1.d0/SQRT(3.d0))
       this%C=(/0.d0,gamma,1.d0-gamma,1.d0/)
       this%lp_of_l=(/1,1,1,2/)
       this%A(2,:)=(/0.d0, gamma, 0.d0 /)
       this%A(3,:)=(/0.d0, 1.d0-2.d0*gamma, gamma /)
       this%A(4,:)=(/0.d0,1.d0/2,1.d0/2/)
         
    CASE(41) !===Equi IRK4  ERK(4,3;1) 
       this%C=(/0.d0,1.d0/4,0.5d0,3.d0/4,1.d0/)
       this%lp_of_l=(/1,1,2,3,4/)
       this%A(2,:)=(/-0.1858665215084591d0, 0.4358665215084591d0, 0.d0, 0.d0/)
       this%A(3,:)=(/-0.4367256409878701d0, 0.5008591194794110d0, 0.4358665215084591d0, 0.d0/)
       this%A(4,:)=(/-0.0423391342724147d0, 0.7701152303135821d0,-0.4136426175496265d0, 0.4358665215084591d0/)
       this%A(5,:)=(/ 0.d0,                 2.d0/3,              -1.d0/3,               2.d0/3/)
       
   CASE(52) !=== Equi IRK(5,4;1) AEJLG
       this%C=(/0.d0,1.d0/5,2.d0/5,3.d0/5,4.d0/5,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5/)
       this%A(2,:)=(/-0.37281606248213511d0, 0.57281606248213512d0, 0.d0,                  0.d0,                  0.d0/)
       this%A(3,:)=(/-0.66007935107985416d0, 0.48726328859771911d0, 0.57281606248213512d0, 0.d0,                  0.d0/)
       this%A(4,:)=(/-0.69934543274239502d0, 1.82596107935553742d0,-1.09943170909527743d0, 0.57281606248213512d0, 0.d0/)
       this%A(5,:)=(/ 0.00000000000000000d0,-0.05144383172900784d0, 1.17898889035791732d0,-0.90036112111104449d0, &
            0.57281606248213512d0/)
       this%A(6,:)=(/-0.10511678454691901d0, 0.87880047152100838d0,-0.58903404061484477d0, 0.46213380485434047d0, &
            0.35321654878641495d0/)      
 
    CASE (62) !===Equi IRK(6,4;1) AEJLG
       this%C=(/0.d0,1.d0/6,2.d0/6,3.d0/6,4.d0/6,5.d0/6,1.d0/)
       this%lp_of_l=(/1,1,2,3,4,5,6,7/)
       this%A(2,:) = (/-0.1113871744697862d0, 0.2780538411364528d0, 0.d0,0.d0,0.d0,0.d0/)
       this%A(3,:) = (/-0.7193507615705692d0, 0.7746302537674498d0, 0.2780538411364528d0, 0.d0,0.d0,0.d0/)
       this%A(4,:) = (/ 0.5518029866688972d0, 0.1104050865166429d0,-0.4402619143219927d0, 0.2780538411364528d0, 0.d0,0.d0/)
       this%A(5,:) = (/ 0.2044212940947437d0, 0.7369116313032833d0,-0.6137248254193539d0, 0.0610047255515406d0, &
            0.2780538411364528d0, 0.d0/)
       this%A(6,:) = (/ 0.0660767687645300d0, 0.0489052670268613d0, 0.2501367454670004d0, 0.5829521002593755d0, &
            -0.3927913893208868d0, 0.2780538411364528d0/)
       this%A(7,:) = (/ 0.0830000000000000d0, 0.1349999999999999d0, 0.1300000000000000d0, 0.4700000000000000d0, &
            -0.28500000000000d0, 0.4670000000000000d0/)
      
    CASE default
       WRITE(*,*) ' BUG in init_bt: wrong sv'
       STOP
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
    CLASS(IBT), INTENT(inout) :: this
    INTEGER :: l
    this%inc_A = 0.d0
    this%inc_C = 0.d0
    this%inc_C = this%C - this%C(this%lp_of_l) !===c(l)-c(l')
    DO l = 2, this%s+1
       this%inc_A(l,:) = this%A(l,:)-this%A(this%lp_of_l(l),:) !===a(l,:) - a(l',:)
    END DO
  END SUBROUTINE inc_bt
  
END MODULE Implicit_Butcher_tableau
