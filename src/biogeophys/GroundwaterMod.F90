module GroundwaterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates prognostic groundwater and Water Table Depth (WTD) considering pumping
  ! and groundwater lateral flow
  !
  ! Usage:
  !
  !   - Call ... in order to compute ... 
  !     . This should be called once per
  !     timestep.
  ! 
  !   - Call ... in order to calculate ... This should be called
  !     exactly ... per time step, 
  !
  !   - Access the timestep's ...
  !
  ! Design notes:
  !
  !   In principle, Theim equation is applied between r = r_e and r = dx (center to center of four
  !   surrounding cells) to take into account the GW pumping.
  !
  !   Assumptions: The head across a given radius is constant (i.e., cone of depression is circular)
  !
  !   First Version: Farshid Felfelani 2019-03-01
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type, get_proc_global
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog
  use clm_varcon       , only : isecspday, degpsec, denh2o, spval, namec
  use clm_varpar       , only : nlevsoi, nlevgrnd
  use clm_time_manager , only : get_step_size
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  use GridcellType     , only : grc                
  use ColumnType       , only : col                
  use PatchType        , only : patch                
  use subgridAveMod    , only : p2c, c2g
  use filterColMod     , only : filter_col_type, col_filter_from_logical_array
  !
  implicit none
  private

  ! !PUBLIC TYPES:
  
  ! This type is public (and its components are public, too) to aid unit testing
  type, public :: irrigation_params_type
     ! Minimum LAI for irrigation
     real(r8) :: 
     integer  :: 
     integer  :: 
     real(r8) :: 
     real(r8) :: 
     real(r8) :: 
     real(r8) :: 
     logical  :: 

  end type irrigation_params_type


  type, public :: irrigation_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: 
     real(r8), pointer, public :: 

     ! Private data members; set in initialization:
     type(irrigation_params_type) :: params
     integer :: dtime                ! land model time step (sec)
     real(r8), pointer :: (:,:) 
     real(r8), pointer :: (:,:) 

     ! Private data members; time-varying:
     real(r8), pointer ::      (:) 
     real(r8), pointer ::      (:) 
     integer , pointer ::      (:) 
     real(r8), pointer ::      (:) 

   contains
     ! Public routines
     procedure, public :: Init => IrrigationInit
     procedure, public :: Restart
     procedure, public :: ApplyIrrigation

     ! Private routines
     procedure, private :: ReadNamelist
     procedure, private :: CheckNamelistValidity   ! Check for validity of input parameters
  end type irrigation_type

  interface irrigation_params_type
     module procedure irrigation_params_constructor
  end interface irrigation_params_type

  
contains

  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================
  
  !------------------------------------------------------------------------
  subroutine UpdateGWTheim(this, bounds, NLFilename, &

    ! !DESCRIPTION:
    !   In principle, Theim equation is applied between 
	!   r = r_e and r = dx (center to center of four
    !   surrounding cells) to take into account the GW pumping.
    !
    !   Assumptions: The head across a given radius is constant 
	!  (i.e., cone of depression is circular) 

    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom
    use decompMod       , only : ldecomp, get_proc_global
	! get_proc_global: total gridcells, landunits, columns, patchs across all processors
    !
    ! !ARGUMENTS:

    !
    ! !LOCAL VARIABLES:	
  
    real(r8), pointer :: neighbors_count(:)  ! complete grid cell array of count
    real(r8), pointer :: B_id_glob(:)        ! complete grid cell array of B_id
    real(r8), pointer :: B_od_long(:)        ! B_od array for all grid cells
    real(r8), pointer :: B_od_glob(:)        ! complete grid cell array of B_od
    real(r8), pointer :: lodgepole_wtgcell_long(:)  ! same pair for...
    real(r8), pointer :: lodgepole_wtgcell_glob(:)  ! ...lodgepole_wtgcell
    integer :: ng                                   ! total number of grid cells
    integer :: ng          ! total number of grid cells
    integer :: p, g        ! patch, gridcell indices
    integer :: g_id, g_od  ! gridcell indices in/out-dispersing cells
    integer :: ier         ! error code
    integer :: yr          ! year
    integer :: mon         ! month
    integer :: day         ! day
    integer :: tod         ! seconds



       ! Initialize for the mpi_allreduce located between the out-dispersal and
       ! in-dispersal loops
       call get_proc_global(ng=ng)
       ! Variables to gather from all PEs, while between the out-dispersal and
       ! in-dispersal loops
       allocate(B_id_glob(ng))
       allocate(B_od_long(ng))
       allocate(B_od_glob(ng))
       allocate(lodgepole_wtgcell_long(ng))
       allocate(lodgepole_wtgcell_glob(ng))

	   
       allocate(neighbors_count(ng))
       neighbors_count = 0._r8  ! initialize counter vector
	   
	   
       call mpi_allreduce(B_od_long, B_od_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)
       call mpi_allreduce(lodgepole_wtgcell_long, lodgepole_wtgcell_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)




       ! Loop to find all out-dispersing beetles' neighbors and
       ! sum each out-dispersing cell's neighbors' lodgepole pine weights so
       ! that neighbors_count is the weighted sum of grid cells with lodgepole
       ! pine > 0
       do g_od = 1, ng
          if (B_od_glob(g_od) > 0._r8) then
             do g_id = 1, ng
                ! identify neighbors with the ixy, jxy indices of grid cells
                if (ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1) then
                   neighbors_count(g_od) = neighbors_count(g_od) + &
                                           lodgepole_wtgcell_glob(g_id)
                end if  ! find surrounding neighbors
             end do  ! g_id loop
          end if  ! B_od > 0
          ! B_id gets beetles from its neighbors
          ! Disperse beetles from out-dispersing cells to in-dispersing cells
          ! Weight out-dispersing beetles by lodgepole pine weights divided by
          ! the neighbors_count sum
          if (neighbors_count(g_od) > 0._r8) then
             do g_id = 1, ng
                if (ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1) then
                   B_id_glob(g_id) = B_id_glob(g_id) + &
                                     lodgepole_wtgcell_glob(g_id) * &
                                     B_od_glob(g_od) / neighbors_count(g_od)
                end if  ! find surrounding neighbors
             end do  ! g_id
          end if  ! neighbors_count > 0
       end do  ! g_od loop
						  


						  
  end subroutine IrrigationInit

  !-----------------------------------------------------------------------


end module GroundwaterMod
