module pntst_struct
!-----------------------------------------------------------------------
! defining pointset data structure
!-----------------------------------------------------------------------
implicit none
private
public pointset

type pointset
    real(8), dimension(:,:), allocatable :: pnts   !points
    integer, dimension(:,:), allocatable :: surfs  !2d surface elements
    integer, dimension(:,:), allocatable :: vols   !3d volume elements
    integer :: dim        !spatial dimension 2d or 3d 
    integer :: totpnts    !total points
    integer :: order      !element order
    integer :: totsurfs   !total surface elements
    integer :: surfshape  !type of surface shape (tri/quad)
    integer :: surfpnts   !total points per surface element
    integer :: totvols    !total 3d volume elements
    integer :: volshape   !type of volume shape (tet/hex)
    integer :: volpnts    !total points per volume element
    real(8) :: dx         !approximated average point spacing
end type

end module pntst_struct

module mat_struct
!-----------------------------------------------------------------------
! defining materials data structure
!-----------------------------------------------------------------------
implicit none
private
public materials

type materials
    real(8), dimension(:), allocatable   ::  d   !diffusivity
    real(8), dimension(:), allocatable   ::  nu  !viscosity or poisson's ratio
    real(8), dimension(:), allocatable   ::  ro  !density
    real(8), dimension(:), allocatable   ::  moe  !modulaous of elasiticity
    integer :: total !total phases
end type

end module mat_struct

module slvr_prmtrs_struct
!-----------------------------------------------------------------------
! defining solver parameters data structure
!-----------------------------------------------------------------------
implicit none
private
public slvr_prmtrs

type slvr_prmtrs
    real(8) :: dt, h, tol, av, p
    real(8) :: pr, ra, t, thot, tcold
    real(8) :: shift_prmtr, shift_surf
    real(8), dimension(3) :: bg_b, bg_l
    real(8), dimension(:), allocatable :: g
    integer :: itrs, nt, prnt_frq, bg_nx
    real(8) :: rbf_alpha !RBF parameter
    logical :: rbf_polyex !augment RBF with polynomials
    integer :: krnl, rbf, mls, wls, sph, gfd, order
    logical :: upwind, shift, lgr, sharp_intrf, vtk, periodic
    real(8) :: segma, a_2, delta, eps_4, theta0, tau0, alpha, gamma, teq, lhf, lamda, cp    !allen cahn solidification
    real(8) :: si_beta, si_eps, si_dt                                      !intrf sharpening
    real(8) :: m, w                                                     !phase field mobility
    integer :: si_nt                                                    !intrf sharpening
    integer :: intrf_slvr                                              !type of intrf solver
    integer :: invert_phi
    !integer :: i_adv_v, i_nrm_v
    !real(8) :: i_adv_v_scl, i_vx_scl, i_vy_scl

end type

end module slvr_prmtrs_struct

module krnl_struct
!-----------------------------------------------------------------------
! defining kernel data structure
!-----------------------------------------------------------------------
implicit none
private
public kernel
type kernel
    integer :: totnbrs
    integer, dimension(:), allocatable   :: nbrs
    real(8), dimension(:), allocatable   :: n, nx, ny, nz, nxx, nyy, nzz, nxy, nxz, nyz, nabla2
    real(8), dimension(:,:), allocatable :: gfd_inv_a, gfd_trm
end type

end module krnl_struct

module bc_struct
!-----------------------------------------------------------------------
! defining boundary conditions data structure
!-----------------------------------------------------------------------
implicit none
private
public bc

type bc
    integer, dimension(:), allocatable   ::  pnts
    real(8), dimension(:), allocatable   ::  vals
    integer                              ::  tot
end type

end module bc_struct
