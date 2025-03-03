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
  ! HISTORY:
  !
  !   The communication between the neighbouring cells is adopted from a FORTRAN code 
  !   simulating the tree mortality from GW wrote by Sam Levis   
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use decompMod         , only : bounds_type, get_proc_global
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use clm_varctl        , only : iulog, use_pumping
  use clm_varcon        , only : isecspday, degpsec, denh2o, spval, namec, rpi, aquifer_water_baseline, watmin, pondmx
  use shr_const_mod     , only : SHR_CONST_PI
  use clm_varpar        , only : nlevsoi, nlevgrnd
  use GridcellType      , only : grc
  use LandunitType      , only : lun    
  use ColumnType        , only : col
  use PatchType         , only : patch                
  use subgridAveMod     , only : p2c, c2g
  use filterColMod      , only : filter_col_type, col_filter_from_logical_array
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use IrrigationMod     , only : irrigation_type
  use spmdMod           , only : iam, masterproc  ! FFelfelani: to get processor number
  use decompMod         , only : get_proc_global, get_proc_bounds, get_clump_bounds,get_proc_clumps  ! FFelfelani: to get number of gridcells
  use clm_time_manager  , only : get_curr_date, get_nstep
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

  ! !PUBLIC TYPES:
  implicit none
  private
  
   type, public :: groundwater_type
     private
    contains
   
     ! Public routines
     procedure, public :: GWFanLatSpinup
     procedure, public :: UpdateGWFanLatPump
     procedure, public :: UpdateGWFanLatTheimPump
     procedure, public :: UpdateGWDefaultPump
     ! Private routines
     procedure, private :: TransmissivityFromFan
     procedure, private :: TheimLateral
     ! procedure, private :: FanLatVal
   end type groundwater_type 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
contains

  ! ======================================================================== 
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================

  !------------------------------------------------------------------------
  subroutine GWFanLatSpinup(this, bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst)

    ! !DESCRIPTION:
    !   In principle, Theim equation is applied between 
    !   r = r_e and r = dx (center to center of four
    !   surrounding cells) to take into account the GW pumping.
    !
    !   Assumptions: The head across a given radius is constant 
    !  (i.e., cone of depression is circular) 

    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom, MPI_INTEGER
    use decompMod       , only : ldecomp, get_proc_global
    use shr_const_mod   , only : SHR_CONST_PI
    use GridcellType    , only : grc
    use clm_time_manager, only : get_step_size, get_curr_date, get_nstep
    use landunit_varcon , only : istwet, istsoil, istice_mec, istcrop

    ! !ARGUMENTS:
    class(groundwater_type)  , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst

    ! !LOCAL VARIABLES:	
 
    character(len=32) :: subname = 'GroundwaterMod' ! subroutine name 
    real(r8), pointer :: zwt_long(:)        , zwt_glob(:)            ! GW_out array for all grid cells; complete grid cell array of GW_out
    real(r8), pointer :: Qn_glob(:)                                  ! Lateral flow in (mm3/sec)
    real(r8), pointer :: QnConserved_long(:), QnConserved_glob(:)	
    real(r8), pointer :: slopelev_long(:), slopelev_glob(:)
    real(r8), pointer :: hgt_long(:), hgt_glob(:)
    real(r8), pointer :: g_totCweight(:)
    real(r8), pointer :: g_watsat_long(:)
    real(r8), pointer :: g_sucst_long(:)
    real(r8), pointer :: g_bsw_long(:)
    real(r8), pointer :: g_cellarea_long(:) , g_cellarea_glob(:)
    real(r8), pointer :: GW_ratio_long(:)
    real(r8), pointer :: AqTransmiss_long(:), AqTransmiss_glob(:)

    integer,  pointer :: ZeroHydroCell(:)          ! cells that have no hydrological columns
    integer,  pointer :: ZeroHydroCell_glob(:)     ! global cells that have no hydrological columns

    real(r8) :: rous, aRatio, colArea              ! aquifer yield (-); area ratio of center/neighbor
    real(r8) :: AqTransmissMean                    ! mean aquifer transmissivity of two adjucent cells
    real(r8) :: widMean, lenMean, deltaxMean       ! mean contact width and lenght of two adjucent cells
    real(r8) :: s_y, dummysum, dummysumabs, dummysum2
    real(r8) :: slopeHeadTop,    slopeHeadBot,    slopeHeadLft,    slopeHeadRgt
    real(r8) :: slopeHeadTopLft, slopeHeadTopRgt, slopeHeadBotLft, slopeHeadBotRgt
    real(r8) :: TransmissMean
    real(r8) :: pump_tot, pump_layer
    real(r8) :: Qgw_lateral_tot, Qgw_lateral_layer
    integer  :: jwt(bounds%begc:bounds%endc)       ! index of the soil layer right above the water table (-)
    real(r8) :: dtime                              ! land model time step (sec)
    real(r8) :: QLateral, QLateralVal
    real(r8) :: l_edge,r_edge,t_edge,b_edge        ! GW lateral on left, right, top, and bottom edge of the cell
	
    integer :: ng, nl, nc, np, nCohorts            ! total number of grid cells,landunits,columns,patches
    integer :: g, c                                ! patch, gridcell, column indices
    integer :: ier                                 ! error code
    integer :: j,fc,i

    integer :: begg, endg       ! beginning and ending gridcell index of current proc
    integer :: begl, endl       ! beginning and ending landunit index of current proc
    integer :: begc, endc       ! beginning and ending column index of current proc
    integer :: begp, endp       ! beginning and ending pft index of current proc

    integer :: year       ! year (0, ...) for nstep
    integer :: month      ! month (1, ..., 12) for nstep
    integer :: day        ! day of month (1, ..., 31) for nstep
    integer :: secs       ! seconds into current date for nstep
    integer :: nstep, cnt

    ! Conversion factors

    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km_to_m     = 1.e3_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: mm_to_m     = 1.e-3_r8
    real(r8), parameter :: m_to_mm     = 1.e3_r8
    real(r8), parameter :: HydroThresh = 0.1_r8
    !-----------------------------------------------------------------------
     associate(&                               
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
		  
          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)

          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input and Output: [real(r8) (:)   ]  water table depth (m)                                        
          wa                 =>    soilhydrology_inst%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
          Qgw_lateral        =>    soilhydrology_inst%Qgw_lateral_col    , & ! Output: [real(r8) (:)   ]  GW lateral flow (mm/s)
          AqTransmiss        =>    soilhydrology_inst%AqTransmiss_col    , & ! Output: [real(r8) (:)   ]  Aquifer Transmissivity(mm2/s)

          Qgw_lateral_g      =>    soilhydrology_inst%Qgw_lateral_grc    & ! Output: [real(r8) (:)   ]  GW lateral flow at gridcell level (mm)

          )
       !-----------------------------------------------
       dtime = get_step_size()
       nstep = get_nstep()
       call get_curr_date (year, month, day, secs)

       ! Initialize for the mpi_allreduce located between the out and 
       ! in loops
       call get_proc_global(ng=ng, nl=nl, nc=nc, np=np, nCohorts=nCohorts)
       ! Variables to gather from all PEs, while between the out and  
       ! in loops 

       allocate(Qn_glob(ng))
       allocate(ZeroHydroCell(ng))
       allocate(ZeroHydroCell_glob(ng))
       allocate(zwt_long(ng))
       allocate(zwt_glob(ng))
       allocate(g_totCweight(ng))
       allocate(g_watsat_long(ng))
       allocate(g_sucst_long(ng))
       allocate(g_bsw_long(ng))
       allocate(g_cellarea_long(ng))
       allocate(g_cellarea_glob(ng))
       allocate(AqTransmiss_long(ng))
       allocate(AqTransmiss_glob(ng))

       allocate(slopelev_long(ng))
       allocate(slopelev_glob(ng))
       allocate(hgt_long(ng))
       allocate(hgt_glob(ng))

       ! Initialize to 0 so as to MPI_SUM zeros in all PEs but one per grid cell
       ! FFELFELANI: the length is the total number of gridcell in the entire domain
       ! but only those cells taken care of each processor get value and the rest
       ! remain zero, finally mpi_allreduce would sum them up
       zwt_long(:)                = 0._r8
       Qn_glob(:)                 = 0._r8
       g_totCweight(:)            = 0._r8
       g_watsat_long(:)           = 0._r8
       g_sucst_long(:)            = 0._r8
       g_bsw_long(:)              = 0._r8
       ZeroHydroCell(:)           = 0
       ZeroHydroCell_glob(:)      = 10000
       g_cellarea_long(:)         = 0._r8
       AqTransmiss_long(:)        = 0._r8
       slopelev_long(:)           = 0._r8
       hgt_long(:)                = 0._r8

       !Initialize to 1e36 to help make errors stand out
       zwt_glob(:)           = 1.e36_r8
       g_cellarea_glob(:)    = 1.e36_r8
       AqTransmiss_glob(:)   = 1.e36_r8
       slopelev_glob(:)      = 1.e36_r8
       hgt_glob(:)           = 1.e36_r8

       AqTransmiss(bounds%begc:bounds%endc) = this%TransmissivityFromFan(bounds, num_hydrologyc, filter_hydrologyc, soilstate_inst, soilhydrology_inst)

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc) 
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer  
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1 
                exit
             end if
          enddo
       end do

       ! To get the weighted average ZWT accross all columns within each cell  
       ! loop over the columns that each processor handles
       ! do c = bounds%begc,bounds%endc
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)
         ! Saving for mpi_allreduce coming up
          g_totCweight(g)    = g_totCweight(g)     + col%wtgcell(c)
          g_cellarea_long(g) = g_cellarea_long(g)  + col%wtgcell(c) * grc%area(g)

          zwt_long(g)          = zwt_long(g)         + col%wtgcell(c) * zwt(c)  ! GW_od goes to history
          AqTransmiss_long(g)  = AqTransmiss_long(g) + col%wtgcell(c) * AqTransmiss(c)

          g_watsat_long(g)  = g_watsat_long(g) + col%wtgcell(c) * watsat(c,nlevsoi) 
          g_sucst_long(g)   = g_sucst_long(g)  + col%wtgcell(c) * sucsat(c,nlevsoi) 
          g_bsw_long(g)     = g_bsw_long(g)    + col%wtgcell(c) * bsw(c,nlevsoi)
       end do  

       do g = bounds%begg,bounds%endg

        slopelev_long(g) = grc%slopelev(g)
        hgt_long(g) = grc%HGT_M(g)

          if (g_totCweight(g) <= 1.000001_r8 .and. g_totCweight(g) > HydroThresh) then

              zwt_long(g)         = zwt_long(g) / g_totCweight(g)
              AqTransmiss_long(g) = AqTransmiss_long(g)/ g_totCweight(g)
              g_watsat_long(g)    = g_watsat_long(g) / g_totCweight(g)
              g_sucst_long(g)     = g_sucst_long(g) / g_totCweight(g)
              g_bsw_long(g)       = g_bsw_long(g) / g_totCweight(g)
          else
              ZeroHydroCell(g) = 1
          end if
       end do


       ! Need action: what happens to the zwt_long of cells g_totCweight(g) < 0.01???????????????
	
       ! ------------------------------------------------
       ! Between the out and in loops
       ! all processors share all the relevant data needed to come up with GW_id
       ! Could the same be done in if (iam == 0) followed by call mpi_bcast?
       ! ------------------------------------------------

       ! Gather zwt_glob from zwt_long, ultimately
       !                  from GW_od

       ! I do MPI_SUM because each processor only has values for the patch/clump its dealing  
       ! with and the rest of the array is zero   
       call mpi_allreduce(zwt_long, zwt_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(AqTransmiss_long, AqTransmiss_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(g_cellarea_long, g_cellarea_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)

       ! we need to exclude those cells that have no hydrologically active columns  
       call mpi_allreduce(ZeroHydroCell, ZeroHydroCell_glob, ng, &
                          MPI_INTEGER, MPI_SUM, mpicom, ier)

       call mpi_allreduce(slopelev_long, slopelev_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(hgt_long, hgt_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)
						  
       call mpi_barrier(mpicom,ier)

       ! gathering the information from the neigboring cells. 
       ! do  g = 1, ng
       do g = bounds%begg,bounds%endg
          l_edge = 0._r8
          r_edge = 0._r8
          t_edge = 0._r8
          b_edge = 0._r8
          ! The GW lateral flow is ruled by Darcy's
          if (ZeroHydroCell_glob(g)== 0) then

             if (ldecomp%gtoplft(g) <= ng .and. ldecomp%gtoplft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoplft(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gtoplft(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtoplft(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gtoplft(g)), zwt_glob(g), zwt_glob(ldecomp%gtoplft(g)), ldecomp%gtoplftUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if

             if (ldecomp%gtop(g) <= ng .and. ldecomp%gtop(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtop(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gtop(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtop(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gtop(g)), zwt_glob(g), zwt_glob(ldecomp%gtop(g)), ldecomp%gtopUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if

             if (ldecomp%gtoprgt(g) <= ng .and. ldecomp%gtoprgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoprgt(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gtoprgt(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtoprgt(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gtoprgt(g)), zwt_glob(g), zwt_glob(ldecomp%gtoprgt(g)), ldecomp%gtoprgtUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if

             if (ldecomp%grgt(g) <= ng .and. ldecomp%grgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%grgt(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%grgt(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%grgt(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%grgt(g)), zwt_glob(g), zwt_glob(ldecomp%grgt(g)), ldecomp%grgtUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if	

             if (ldecomp%gbotrgt(g) <= ng .and. ldecomp%gbotrgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotrgt(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gbotrgt(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbotrgt(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gbotrgt(g)), zwt_glob(g), zwt_glob(ldecomp%gbotrgt(g)), ldecomp%gbotrgtUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if

             if (ldecomp%gbot(g) <= ng .and. ldecomp%gbot(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbot(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gbot(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbot(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gbot(g)), zwt_glob(g), zwt_glob(ldecomp%gbot(g)), ldecomp%gbotUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if

             if (ldecomp%gbotlft(g) <= ng .and. ldecomp%gbotlft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotlft(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gbotlft(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbotlft(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gbotlft(g)), zwt_glob(g), zwt_glob(ldecomp%gbotlft(g)), ldecomp%gbotlftUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if

             if (ldecomp%glft(g) <= ng .and. ldecomp%glft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%glft(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%glft(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%glft(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%glft(g)), zwt_glob(g), zwt_glob(ldecomp%glft(g)), ldecomp%glftUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

             end if

          ! IF there is pumping, the GW lateral flow is ruled by Combination of Darcy's and Theim
          ! else if (ZeroHydroCell_glob(g)== 0 .and. GW_ratio_long(g) * qirrig_long(g) > 0._r8) then

             Qgw_lateral_g(g) = (Qn_glob(g) / (grc%area(g) * km2_to_mm2)) * dtime
          end if
 
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)

          aRatio  = col%wtgcell(c) / g_totCweight(g)
          colArea = col%wtgcell(c) * grc%area(g) * km2_to_mm2
          Qgw_lateral(c) = Qn_glob(g) * aRatio / colArea  !unit is converted to mm/sec 

          rous = watsat(c,nlevsoi) &
               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
          rous=max(rous,0.02_r8)
          ! if (nstep  > 400 .and. c == 838456) write(*,*) 'c,', c, Qgw_lateral(c)
          if (Qgw_lateral(c) > 0._r8) then

                Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
                if(jwt(c) == nlevsoi) then             
                   ! wa(c)  = wa(c) + Qgw_lateral_tot
                   zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                else   
                   do j = jwt(c)+1, 1,-1
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=min(Qgw_lateral_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                       Qgw_lateral_layer=max(Qgw_lateral_layer,0._r8)
				
                       if(s_y > 0._r8) zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer
                       if (Qgw_lateral_tot <= 0._r8) exit
                   enddo
                end if

          else if (Qgw_lateral(c) < 0._r8) then

              Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
              !! --  water table is below the soil column  -------------------------------------- 
              if(jwt(c) == nlevsoi) then             
                 ! wa(c)  = wa(c) + Qgw_lateral_tot
                 zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
              else                                
                 !! -- water table within soil layers 1-9  -------------------------------------
                 !! ============================== Qgw_lateral_tot ========================================= 
                 !! --  Now remove water via Qgw_lateral_tot

                 !! should never be positive... but include for completeness 
                 if(Qgw_lateral_tot > 0.) then !rising water table

                    call endrun(msg="Qgw_lateral_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                 else ! deepening water table
                    do j = jwt(c)+1, nlevsoi
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=max(Qgw_lateral_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                       Qgw_lateral_layer=min(Qgw_lateral_layer,0._r8)

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer

                       if (Qgw_lateral_tot >= 0.) then 
                         zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8
                          exit
                       else
                          zwt(c) = zi(c,j)
                       endif
                    enddo

                    !! --  remove residual Qgw_lateral_tot  ---------------------------------------------
                    zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                    ! wa(c) = wa(c) + Qgw_lateral_tot
                 endif

                 !! -- recompute jwt  ---------------------------------------------------------
                 !! allow jwt to equal zero when zwt is in top layer
                 jwt(c) = nlevsoi
                 do j = 1,nlevsoi
                    if(zwt(c) <= zi(c,j)) then
                       jwt(c) = j-1
                       exit
                    end if
                 enddo
              end if! end of jwt if construct

              zwt(c) = max(0.0_r8,zwt(c))
              ! zwt(c) = min(80._r8,zwt(c))
          end if
       end do

       deallocate(Qn_glob)
       deallocate(ZeroHydroCell)
       deallocate(ZeroHydroCell_glob)
       deallocate(zwt_long)
       deallocate(zwt_glob)
       deallocate(g_totCweight)
       deallocate(g_watsat_long)
       deallocate(g_sucst_long)
       deallocate(g_bsw_long)
       deallocate(g_cellarea_long)
       deallocate(g_cellarea_glob)
       deallocate(AqTransmiss_long)
       deallocate(AqTransmiss_glob)
       deallocate(slopelev_long)
       deallocate(slopelev_glob) 
       deallocate(hgt_long)
       deallocate(hgt_glob) 
     end associate
  end subroutine GWFanLatSpinup

  !------------------------------------------------------------------------
  subroutine UpdateGWFanLatPump(this, bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst,waterstate_inst, irrigation_inst, waterflux_inst)

    ! !DESCRIPTION:
    !   In principle, Theim equation is applied between 
    !   r = r_e and r = dx (center to center of four
    !   surrounding cells) to take into account the GW pumping.
    !
    !   Assumptions: The head across a given radius is constant 
    !  (i.e., cone of depression is circular) 

    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom, MPI_INTEGER
    use decompMod       , only : ldecomp, get_proc_global
    use shr_const_mod   , only : SHR_CONST_PI
    use GridcellType    , only : grc
    use clm_time_manager, only : get_step_size, get_curr_date, get_nstep
    use landunit_varcon , only : istwet, istsoil, istice_mec, istcrop

    ! !ARGUMENTS:
    class(groundwater_type)  , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(waterstate_type)    , intent(inout) :: waterstate_inst
    type(irrigation_type)    , intent(in)    :: irrigation_inst	
    type(waterflux_type)     , intent(inout) :: waterflux_inst

    ! !LOCAL VARIABLES:	
 
    character(len=32) :: subname = 'GroundwaterMod' ! subroutine name 
    real(r8), pointer :: zwt_long(:)        , zwt_glob(:)            ! GW_out array for all grid cells; complete grid cell array of GW_out
    real(r8), pointer :: Qn_glob(:)                                  ! Lateral flow in (mm3/sec)
    real(r8), pointer :: QnConserved_long(:), QnConserved_glob(:)	
    real(r8), pointer :: slopelev_long(:), slopelev_glob(:)
    real(r8), pointer :: hgt_long(:), hgt_glob(:)
    real(r8), pointer :: g_totCweight(:)
    real(r8), pointer :: g_watsat_long(:)
    real(r8), pointer :: g_sucst_long(:)
    real(r8), pointer :: g_bsw_long(:)
    real(r8), pointer :: g_cellarea_long(:) , g_cellarea_glob(:)
    real(r8), pointer :: GW_ratio_long(:)
    real(r8), pointer :: AqTransmiss_long(:), AqTransmiss_glob(:)

    integer,  pointer :: ZeroHydroCell(:)          ! cells that have no hydrological columns
    integer,  pointer :: ZeroHydroCell_glob(:)     ! global cells that have no hydrological columns

    real(r8) :: rous, aRatio, colArea              ! aquifer yield (-); area ratio of center/neighbor
    real(r8) :: AqTransmissMean                    ! mean aquifer transmissivity of two adjucent cells
    real(r8) :: widMean, lenMean, deltaxMean       ! mean contact width and lenght of two adjucent cells
    real(r8) :: s_y, dummysum, dummysumabs, dummysum2
    real(r8) :: slopeHeadTop,    slopeHeadBot,    slopeHeadLft,    slopeHeadRgt
    real(r8) :: slopeHeadTopLft, slopeHeadTopRgt, slopeHeadBotLft, slopeHeadBotRgt
    real(r8) :: TransmissMean
    real(r8) :: pump_tot, pump_layer
    real(r8) :: Qgw_lateral_tot, Qgw_lateral_layer
    integer  :: jwt(bounds%begc:bounds%endc)       ! index of the soil layer right above the water table (-)
    real(r8) :: dtime                              ! land model time step (sec)
    real(r8) :: QLateral, QLateralVal
    real(r8) :: l_edge,r_edge,t_edge,b_edge        ! GW lateral on left, right, top, and bottom edge of the cell
	
    integer :: ng, nl, nc, np, nCohorts            ! total number of grid cells,landunits,columns,patches
    integer :: g, c                                ! patch, gridcell, column indices
    integer :: ier                                 ! error code
    integer :: j,fc,i

    integer :: begg, endg       ! beginning and ending gridcell index of current proc
    integer :: begl, endl       ! beginning and ending landunit index of current proc
    integer :: begc, endc       ! beginning and ending column index of current proc
    integer :: begp, endp       ! beginning and ending pft index of current proc

    integer :: year       ! year (0, ...) for nstep
    integer :: month      ! month (1, ..., 12) for nstep
    integer :: day        ! day of month (1, ..., 31) for nstep
    integer :: secs       ! seconds into current date for nstep
    integer :: nstep, cnt

    ! Conversion factors

    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km_to_m     = 1.e3_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: mm_to_m     = 1.e-3_r8
    real(r8), parameter :: m_to_mm     = 1.e3_r8
    real(r8), parameter :: HydroThresh = 0.1_r8
    !-----------------------------------------------------------------------
     associate(&                               
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
          GW_ratio           =>    col%GW_ratio                          , & ! Input:  [real(r8) (:)   ]  USGS GW ratio as irrigation source                                                
                                              
          qflx_irrig         =>    irrigation_inst%qflx_irrig_col        , & ! irrigation flux (mm H2O /s)

          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)

          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input and Output: [real(r8) (:)   ]  water table depth (m)                                        
          wa                 =>    soilhydrology_inst%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
          qcharge            =>    soilhydrology_inst%qcharge_col        , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
          Qgw_lateral        =>    soilhydrology_inst%Qgw_lateral_col    , & ! Output: [real(r8) (:)   ]  GW lateral flow (mm/s)
          AqTransmiss        =>    soilhydrology_inst%AqTransmiss_col    , & ! Output: [real(r8) (:)   ]  Aquifer Transmissivity(mm2/s)
          Pump_wa            =>    soilhydrology_inst%Pump_wa_col        , & ! Output: [real(r8) (:)   ]  Pumped Water from the aquifer(mm/s)
          QlatField_north    =>    soilhydrology_inst%QlatField_northing_grc , & !  Output: [real(r8) (:)   ] Northward lateral GW flow
          QlatField_east     =>    soilhydrology_inst%QlatField_easting_grc  , & !  Output: [real(r8) (:)   ] Eastward lateral GW flow

          Qgw_lateral_g      =>    soilhydrology_inst%Qgw_lateral_grc    , & ! Output: [real(r8) (:)   ]  GW lateral flow at gridcell level (mm)

          qflx_drain         =>    waterflux_inst%qflx_drain_col         , & ! Input and Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)                    

          h2osoi_liq         =>    waterstate_inst%h2osoi_liq_col        & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          )
       !-----------------------------------------------
       dtime = get_step_size()
       nstep = get_nstep()
       call get_curr_date (year, month, day, secs)		  
 
       ! Initialize for the mpi_allreduce located between the out and 
       ! in loops
       call get_proc_global(ng=ng, nl=nl, nc=nc, np=np, nCohorts=nCohorts)
       ! Variables to gather from all PEs, while between the out and  
       ! in loops 

       allocate(Qn_glob(ng))
       allocate(QnConserved_long(ng))
       allocate(QnConserved_glob(ng))
       allocate(ZeroHydroCell(ng))
       allocate(ZeroHydroCell_glob(ng))
       allocate(zwt_long(ng))
       allocate(zwt_glob(ng))
       allocate(g_totCweight(ng))
       allocate(g_watsat_long(ng))
       allocate(g_sucst_long(ng))
       allocate(g_bsw_long(ng))
       allocate(g_cellarea_long(ng))
       allocate(g_cellarea_glob(ng))
       allocate(GW_ratio_long(ng))
       allocate(AqTransmiss_long(ng))
       allocate(AqTransmiss_glob(ng))

       allocate(slopelev_long(ng))
       allocate(slopelev_glob(ng))
       allocate(hgt_long(ng))
       allocate(hgt_glob(ng))

       ! Initialize to 0 so as to MPI_SUM zeros in all PEs but one per grid cell
       ! FFELFELANI: the length is the total number of gridcell in the entire domain
       ! but only those cells taken care of each processor get value and the rest
       ! remain zero, finally mpi_allreduce would sum them up
       zwt_long(:)                = 0._r8
       Qn_glob(:)                 = 0._r8
       QnConserved_long(:)        = 0._r8
       QnConserved_glob(:)        = 0._r8
       g_totCweight(:)            = 0._r8
       g_watsat_long(:)           = 0._r8
       g_sucst_long(:)            = 0._r8
       g_bsw_long(:)              = 0._r8
       ZeroHydroCell(:)           = 0
       ZeroHydroCell_glob(:)      = 10000
       g_cellarea_long(:)         = 0._r8
       GW_ratio_long(:)           = 0._r8
       AqTransmiss_long(:)        = 0._r8
       slopelev_long(:)           = 0._r8
       hgt_long(:)                = 0._r8
	   
       !Initialize to 1e36 to help make errors stand out
       zwt_glob(:)           = 1.e36_r8
       g_cellarea_glob(:)    = 1.e36_r8
       AqTransmiss_glob(:)   = 1.e36_r8
       slopelev_glob(:)      = 1.e36_r8
       hgt_glob(:)           = 1.e36_r8

       AqTransmiss(bounds%begc:bounds%endc) = this%TransmissivityFromFan(bounds, num_hydrologyc, filter_hydrologyc, soilstate_inst, soilhydrology_inst)

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc) 
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer  
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1 
                exit
             end if
          enddo
       end do

      ! Removing the pumped water from the soil column	   
       if (use_pumping == .true.) then
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             g = col%gridcell(c)
                 !!! use analytical expression for aquifer specific yield
                 rous = watsat(c,nlevsoi) &
                      * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
                 rous=max(rous,0.02_r8)

                 pump_tot = - GW_ratio(c) * qflx_irrig(c) * dtime
                 Pump_wa(c) = GW_ratio(c) * qflx_irrig(c)
                 ! if (nstep > 400 .and. c == 838456) write(*,*) 'c,', c, pump_tot, Pump_wa(c) 
                 !!!--  water table is below the soil column  --------------------------------------
                 if(jwt(c) == nlevsoi) then             
                    wa(c)  = wa(c) + pump_tot
                    zwt(c) = zwt(c) - pump_tot/1000._r8/rous
 
                 else                                
                    !!!-- water table within soil layers 1-9  --------------------------------------
                    !!!============================== pump_tot ========================================= 
                    !!!--  Now remove water via pump_tot

                    !!!should never be positive... but include for completeness
                    if(pump_tot > 0.) then !rising water table

                       call endrun(msg="pump_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                    else ! deepening water table
                       do j = jwt(c)+1, nlevsoi
                          !!! use analytical expression for specific yield
                          s_y = watsat(c,j) &
                               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                          s_y=max(s_y,0.02_r8)

                          pump_layer=max(pump_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                          pump_layer=min(pump_layer,0._r8)
                          h2osoi_liq(c,j) = h2osoi_liq(c,j) + pump_layer

                          pump_tot = pump_tot - pump_layer

                          if (pump_tot >= 0.) then 
                             zwt(c) = zwt(c) - pump_layer/s_y/1000._r8
                             exit
                          else
                             zwt(c) = zi(c,j)
                          endif
                       enddo

                       !!!--  remove residual pump_tot  ---------------------------------------------
                       zwt(c) = zwt(c) - pump_tot/1000._r8/rous
                       wa(c)  = wa(c) + pump_tot
                    endif

                    !!!-- recompute jwt  ---------------------------------------------------------
                    !!! allow jwt to equal zero when zwt is in top layer
                    jwt(c) = nlevsoi
                    do j = 1,nlevsoi
                       if(zwt(c) <= zi(c,j)) then
                          jwt(c) = j-1
                          exit
                       end if
                    enddo
                 end if! end of jwt if construct

                 zwt(c) = max(0.0_r8,zwt(c))
                 ! zwt(c) = min(80._r8,zwt(c))    
          end do
       end if
       ! To get the weighted average ZWT accross all columns within each cell  
       ! loop over the columns that each processor handles
       ! do c = bounds%begc,bounds%endc
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)
         ! Saving for mpi_allreduce coming up
          g_totCweight(g)    = g_totCweight(g)     + col%wtgcell(c)
          g_cellarea_long(g) = g_cellarea_long(g)  + col%wtgcell(c) * grc%area(g)

          zwt_long(g)          = zwt_long(g)         + col%wtgcell(c) * zwt(c)  ! GW_od goes to history
          AqTransmiss_long(g)  = AqTransmiss_long(g) + col%wtgcell(c) * AqTransmiss(c)

          GW_ratio_long(g) = GW_ratio_long(g) + col%wtgcell(c) * GW_ratio(c)

          g_watsat_long(g)  = g_watsat_long(g) + col%wtgcell(c) * watsat(c,nlevsoi) 
          g_sucst_long(g)   = g_sucst_long(g)  + col%wtgcell(c) * sucsat(c,nlevsoi) 
          g_bsw_long(g)     = g_bsw_long(g)    + col%wtgcell(c) * bsw(c,nlevsoi)
       end do  

       do g = bounds%begg,bounds%endg

        slopelev_long(g) = grc%slopelev(g)
        hgt_long(g) = grc%HGT_M(g)

          if (g_totCweight(g) <= 1.000001_r8 .and. g_totCweight(g) > HydroThresh) then

              zwt_long(g)         = zwt_long(g) / g_totCweight(g)
              AqTransmiss_long(g) = AqTransmiss_long(g)/ g_totCweight(g)
              GW_ratio_long(g)    = GW_ratio_long(g) / g_totCweight(g)
              g_watsat_long(g)    = g_watsat_long(g) / g_totCweight(g)
              g_sucst_long(g)     = g_sucst_long(g) / g_totCweight(g)
              g_bsw_long(g)       = g_bsw_long(g) / g_totCweight(g)
          else
              ZeroHydroCell(g) = 1
          end if
       end do


       ! Need action: what happens to the zwt_long of cells g_totCweight(g) < 0.01???????????????
	
       ! ------------------------------------------------
       ! Between the out and in loops
       ! all processors share all the relevant data needed to come up with GW_id
       ! Could the same be done in if (iam == 0) followed by call mpi_bcast?
       ! ------------------------------------------------

       ! Gather zwt_glob from zwt_long, ultimately
       !                  from GW_od

       ! I do MPI_SUM because each processor only has values for the patch/clump its dealing  
       ! with and the rest of the array is zero   
       call mpi_allreduce(zwt_long, zwt_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(AqTransmiss_long, AqTransmiss_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(g_cellarea_long, g_cellarea_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)

       ! we need to exclude those cells that have no hydrologically active columns  
       call mpi_allreduce(ZeroHydroCell, ZeroHydroCell_glob, ng, &
                          MPI_INTEGER, MPI_SUM, mpicom, ier)

       call mpi_allreduce(slopelev_long, slopelev_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 
       call mpi_allreduce(hgt_long, hgt_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 
						  
       call mpi_barrier(mpicom,ier)

       ! gathering the information from the neigboring cells. 
       do g = bounds%begg,bounds%endg
          l_edge = 0._r8
          r_edge = 0._r8
          t_edge = 0._r8
          b_edge = 0._r8
          ! The GW lateral flow is ruled by Darcy's
          if (ZeroHydroCell_glob(g)== 0) then

             if (ldecomp%gtoplft(g) <= ng .and. ldecomp%gtoplft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoplft(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gtoplft(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtoplft(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gtoplft(g)), zwt_glob(g), zwt_glob(ldecomp%gtoplft(g)), ldecomp%gtoplftUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 l_edge = l_edge + QLateralVal * sqrt(2._r8)/2._r8  !why positive: because the + direction in eastward
                 t_edge = t_edge - QLateralVal * sqrt(2._r8)/2._r8  !why negative: because the + direction in upward

             end if

             if (ldecomp%gtop(g) <= ng .and. ldecomp%gtop(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtop(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gtop(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtop(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gtop(g)), zwt_glob(g), zwt_glob(ldecomp%gtop(g)), ldecomp%gtopUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 t_edge = t_edge - QLateralVal !why negative: because the + direction in upward

             end if

             if (ldecomp%gtoprgt(g) <= ng .and. ldecomp%gtoprgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoprgt(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gtoprgt(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtoprgt(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gtoprgt(g)), zwt_glob(g), zwt_glob(ldecomp%gtoprgt(g)), ldecomp%gtoprgtUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 r_edge = r_edge - QLateralVal * sqrt(2._r8)/2._r8
                 t_edge = t_edge - QLateralVal * sqrt(2._r8)/2._r8

             end if

             if (ldecomp%grgt(g) <= ng .and. ldecomp%grgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%grgt(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%grgt(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%grgt(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%grgt(g)), zwt_glob(g), zwt_glob(ldecomp%grgt(g)), ldecomp%grgtUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 r_edge = r_edge - QLateralVal

             end if	

             if (ldecomp%gbotrgt(g) <= ng .and. ldecomp%gbotrgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotrgt(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gbotrgt(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbotrgt(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gbotrgt(g)), zwt_glob(g), zwt_glob(ldecomp%gbotrgt(g)), ldecomp%gbotrgtUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 r_edge = r_edge - QLateralVal * sqrt(2._r8)/2._r8
                 b_edge = b_edge + QLateralVal * sqrt(2._r8)/2._r8

             end if

             if (ldecomp%gbot(g) <= ng .and. ldecomp%gbot(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbot(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gbot(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbot(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gbot(g)), zwt_glob(g), zwt_glob(ldecomp%gbot(g)), ldecomp%gbotUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 b_edge = b_edge + QLateralVal

             end if

             if (ldecomp%gbotlft(g) <= ng .and. ldecomp%gbotlft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotlft(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%gbotlft(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbotlft(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%gbotlft(g)), zwt_glob(g), zwt_glob(ldecomp%gbotlft(g)), ldecomp%gbotlftUP(g), '___Diag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 l_edge = l_edge + QLateralVal * sqrt(2._r8)/2._r8
                 b_edge = b_edge + QLateralVal * sqrt(2._r8)/2._r8

             end if

             if (ldecomp%glft(g) <= ng .and. ldecomp%glft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%glft(g)) == 0) then

                 QLateralVal = FanLatVal_hgt(AqTransmiss_glob(g), AqTransmiss_glob(ldecomp%glft(g)), g_cellarea_glob(g), g_cellarea_glob(ldecomp%glft(g)), &
                           hgt_glob(g), hgt_glob(ldecomp%glft(g)), zwt_glob(g), zwt_glob(ldecomp%glft(g)), ldecomp%glftUP(g), 'nonDiag')

                 Qn_glob(g) = Qn_glob(g) + QLateralVal

                 l_edge = l_edge + QLateralVal

             end if

          ! IF there is pumping, the GW lateral flow is ruled by Combination of Darcy's and Theim
          ! else if (ZeroHydroCell_glob(g)== 0 .and. GW_ratio_long(g) * qirrig_long(g) > 0._r8) then

             QlatField_north(g) = (t_edge + b_edge)/2._r8
             QlatField_east(g)  = (r_edge + l_edge)/2._r8

             Qgw_lateral_g(g) = (Qn_glob(g) / (grc%area(g) * km2_to_mm2)) * dtime
          end if
 
       end do

       do g = bounds%begg,bounds%endg
        QnConserved_long(g) = Qn_glob(g)
       end do
       call mpi_allreduce(QnConserved_long, QnConserved_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)

          aRatio  = col%wtgcell(c) / g_totCweight(g)
          colArea = col%wtgcell(c) * grc%area(g) * km2_to_mm2
          Qgw_lateral(c) = Qn_glob(g) * aRatio / colArea  !unit is converted to mm/sec 

          rous = watsat(c,nlevsoi) &
               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
          rous=max(rous,0.02_r8)
          if (Qgw_lateral(c) > 0._r8) then

                Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
                if(jwt(c) == nlevsoi) then             
                   wa(c)  = wa(c) + Qgw_lateral_tot
                   zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                else   
                   do j = jwt(c)+1, 1,-1
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=min(Qgw_lateral_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                       Qgw_lateral_layer=max(Qgw_lateral_layer,0._r8)

                       h2osoi_liq(c,j) = h2osoi_liq(c,j) + Qgw_lateral_layer
				
                       if(s_y > 0._r8) zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer
                       if (Qgw_lateral_tot <= 0._r8) exit
                   enddo

                   ! add residual to the sub-surface runoff (both lateral flow and runoff are positive here)
				   qflx_drain(c) = qflx_drain(c) + Qgw_lateral_tot / dtime

                end if

          else if (Qgw_lateral(c) < 0._r8) then

              Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
              !! --  water table is below the soil column  -------------------------------------- 
              if(jwt(c) == nlevsoi) then             
                 wa(c)  = wa(c) + Qgw_lateral_tot
                 zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
              else                                
                 !! -- water table within soil layers 1-9  -------------------------------------
                 !! ============================== Qgw_lateral_tot ========================================= 
                 !! --  Now remove water via Qgw_lateral_tot

                 !! should never be positive... but include for completeness 
                 if(Qgw_lateral_tot > 0.) then !rising water table

                    call endrun(msg="Qgw_lateral_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                 else ! deepening water table
                    do j = jwt(c)+1, nlevsoi
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=max(Qgw_lateral_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                       Qgw_lateral_layer=min(Qgw_lateral_layer,0._r8)
                       h2osoi_liq(c,j) = h2osoi_liq(c,j) + Qgw_lateral_layer

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer

                       if (Qgw_lateral_tot >= 0.) then 
                         zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8
                          exit
                       else
                          zwt(c) = zi(c,j)
                       endif
                    enddo

                    !! --  remove residual Qgw_lateral_tot  ---------------------------------------------
                    zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                    wa(c) = wa(c) + Qgw_lateral_tot
                 endif

                 !! -- recompute jwt  ---------------------------------------------------------
                 !! allow jwt to equal zero when zwt is in top layer
                 jwt(c) = nlevsoi
                 do j = 1,nlevsoi
                    if(zwt(c) <= zi(c,j)) then
                       jwt(c) = j-1
                       exit
                    end if
                 enddo
              end if! end of jwt if construct

              zwt(c) = max(0.0_r8,zwt(c))
              ! zwt(c) = min(80._r8,zwt(c))
          end if
       end do

       deallocate(Qn_glob)
       deallocate(QnConserved_long)
       deallocate(QnConserved_glob)
       deallocate(ZeroHydroCell)
       deallocate(ZeroHydroCell_glob)
       deallocate(zwt_long)
       deallocate(zwt_glob)
       deallocate(g_totCweight)
       deallocate(g_watsat_long)
       deallocate(g_sucst_long)
       deallocate(g_bsw_long)
       deallocate(g_cellarea_long)
       deallocate(g_cellarea_glob)
       deallocate(GW_ratio_long)
       deallocate(AqTransmiss_long)
       deallocate(AqTransmiss_glob)
       deallocate(slopelev_long)
       deallocate(slopelev_glob)
       deallocate(hgt_long)
       deallocate(hgt_glob)
     end associate
  end subroutine UpdateGWFanLatPump

  !------------------------------------------------------------------------
  subroutine UpdateGWFanLatTheimPump(this, bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst,waterstate_inst, irrigation_inst)

    ! !DESCRIPTION:
    !   In principle, Theim equation is applied between 
    !   r = r_e and r = dx (center to center of four
    !   surrounding cells) to take into account the GW pumping.
    !
    !   Assumptions: The head across a given radius is constant 
    !  (i.e., cone of depression is circular) 

    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom, MPI_INTEGER
    use decompMod       , only : ldecomp, get_proc_global
    use shr_const_mod   , only : SHR_CONST_PI
    use GridcellType    , only : grc
    use clm_time_manager, only : get_step_size, get_curr_date, get_nstep
    use landunit_varcon , only : istwet, istsoil, istice_mec, istcrop

    ! !ARGUMENTS:
    class(groundwater_type)  , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(waterstate_type)    , intent(inout) :: waterstate_inst
    type(irrigation_type)    , intent(in)    :: irrigation_inst	

    ! !LOCAL VARIABLES:	
 
    character(len=32) :: subname = 'GroundwaterMod' ! subroutine name 
    real(r8), pointer :: zwt_long(:)         , zwt_glob(:)            ! GW_out array for all grid cells; complete grid cell array of GW_out
    real(r8), pointer :: Qn_glob(:)                                   ! Lateral flow in (mm3/sec)
    real(r8), pointer :: g_totCweight(:)
    real(r8), pointer :: g_watsat_long(:)
    real(r8), pointer :: g_sucst_long(:)
    real(r8), pointer :: g_bsw_long(:)
    real(r8), pointer :: g_cellarea_long(:)  , g_cellarea_glob(:)
    real(r8), pointer :: Pump_wa_long(:)     , Pump_wa_glob(:)
    real(r8), pointer :: GW_ratio_long(:)
    real(r8), pointer :: AqTransmiss_long(:) , AqTransmiss_glob(:)
    real(r8), pointer :: qcharge_long(:)     , qcharge_glob(:)

    integer,  pointer :: ZeroHydroCell(:)          ! cells that have no hydrological columns
    integer,  pointer :: ZeroHydroCell_glob(:)     ! global cells that have no hydrological columns

    real(r8) :: rous, aRatio, colArea              ! aquifer yield (-); area ratio of center/neighbor
    real(r8) :: AqTransmissMean                    ! mean aquifer transmissivity of two adjucent cells
    real(r8) :: widMean, lenMean, deltaxMean       ! mean contact width and lenght of two adjucent cells
    real(r8) :: s_y, dummysum, dummysum2
    real(r8) :: TransmissMean
    real(r8) :: pump_tot, pump_layer
    real(r8) :: Qgw_lateral_tot, Qgw_lateral_layer
    real(r8) :: QLateral    
    integer  :: jwt(bounds%begc:bounds%endc)       ! index of the soil layer right above the water table (-)
    real(r8) :: dtime                              ! land model time step (sec)
    real(r8) :: l_edge,r_edge,t_edge,b_edge        ! GW lateral on left, right, top, and bottom edge of the cell
	
    integer :: ng, nl, nc, np, nCohorts            ! total number of grid cells,landunits,columns,patches
    integer :: g, c                                ! patch, gridcell, column indices
    integer :: ier                                 ! error code
    integer :: j,fc,i

    integer :: begg, endg       ! beginning and ending gridcell index of current proc
    integer :: begl, endl       ! beginning and ending landunit index of current proc
    integer :: begc, endc       ! beginning and ending column index of current proc
    integer :: begp, endp       ! beginning and ending pft index of current proc

    integer :: year       ! year (0, ...) for nstep
    integer :: month      ! month (1, ..., 12) for nstep
    integer :: day        ! day of month (1, ..., 31) for nstep
    integer :: secs       ! seconds into current date for nstep
    integer :: nstep

    ! Conversion factors

    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: mm_to_m     = 1.e-3_r8
    real(r8), parameter :: m_to_mm     = 1.e3_r8
    real(r8), parameter :: HydroThresh = 0.1_r8
    !-----------------------------------------------------------------------
     associate(&                               
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
          GW_ratio           =>    col%GW_ratio                          , & ! Input:  [real(r8) (:)   ]  USGS GW ratio as irrigation source                                                
                                              
          qflx_irrig         =>    irrigation_inst%qflx_irrig_col        , & ! irrigation flux (mm H2O /s)

          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)

          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input and Output: [real(r8) (:)   ]  water table depth (m)                                        
          wa                 =>    soilhydrology_inst%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
          qcharge            =>    soilhydrology_inst%qcharge_col        , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
          Qgw_lateral        =>    soilhydrology_inst%Qgw_lateral_col    , & ! Output: [real(r8) (:)   ]  GW lateral flow (mm/s)
          AqTransmiss        =>    soilhydrology_inst%AqTransmiss_col    , & ! Output: [real(r8) (:)   ]  Aquifer Transmissivity(mm2/s)
          Pump_wa            =>    soilhydrology_inst%Pump_wa_col        , & ! Output: [real(r8) (:)   ]  Pumped Water from the aquifer(mm/s)
          QlatField_north    =>    soilhydrology_inst%QlatField_northing_grc , & !  Output: [real(r8) (:)   ] Northward lateral GW flow
          QlatField_east     =>    soilhydrology_inst%QlatField_easting_grc  , & !  Output: [real(r8) (:)   ] Eastward lateral GW flow
		  
          h2osoi_liq         =>    waterstate_inst%h2osoi_liq_col        & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          )
       !-----------------------------------------------
       dtime = get_step_size()
       nstep = get_nstep()
       call get_curr_date (year, month, day, secs)
       ! if (masterproc) then
       !if (iam == 200) then
       !   write(*,*) 'year, month, day, secs, dtime, ', year, month, day, secs, dtime
       !end if    		  
 
       ! Initialize for the mpi_allreduce located between the out and 
       ! in loops
       call get_proc_global(ng=ng, nl=nl, nc=nc, np=np, nCohorts=nCohorts)
       ! Variables to gather from all PEs, while between the out and  
       ! in loops 

       allocate(Qn_glob(ng))
       allocate(ZeroHydroCell(ng))
       allocate(ZeroHydroCell_glob(ng))
       allocate(zwt_long(ng))
       allocate(zwt_glob(ng))
       allocate(g_totCweight(ng))
       allocate(g_watsat_long(ng))
       allocate(g_sucst_long(ng))
       allocate(g_bsw_long(ng))
       allocate(g_cellarea_long(ng))
       allocate(g_cellarea_glob(ng))
       allocate(Pump_wa_long(ng))
       allocate(Pump_wa_glob(ng))
       allocate(GW_ratio_long(ng))
       allocate(AqTransmiss_long(ng))
       allocate(AqTransmiss_glob(ng))
       allocate(qcharge_long(ng))
       allocate(qcharge_glob(ng))


       ! Initialize to 0 so as to MPI_SUM zeros in all PEs but one per grid cell
       ! FFELFELANI: the length is the total number of gridcell in the entire domain
       ! but only those cells taken care of each processor get value and the rest
       ! remain zero, finally mpi_allreduce would sum them up
       zwt_long(:)                = 0._r8
       Qn_glob(:)                 = 0._r8
       g_totCweight(:)            = 0._r8
       g_watsat_long(:)           = 0._r8
       g_sucst_long(:)            = 0._r8
       g_bsw_long(:)              = 0._r8
       ZeroHydroCell(:)           = 0
       ZeroHydroCell_glob(:)      = 10000
       g_cellarea_long(:)         = 0._r8
       Pump_wa_long(:)            = 0._r8
       GW_ratio_long(:)           = 0._r8
       AqTransmiss_long(:)        = 0._r8
       qcharge_long(:)            = 0._r8

       !Initialize to 1e36 to help make errors stand out
       zwt_glob(:)           = 1.e36_r8
       g_cellarea_glob(:)    = 1.e36_r8
       AqTransmiss_glob(:)   = 1.e36_r8
       Pump_wa_glob(:)       = 1.e36_r8
       qcharge_glob(:)       = 1.e36_r8
	   
       AqTransmiss(bounds%begc:bounds%endc) = this%TransmissivityFromFan(bounds, num_hydrologyc, filter_hydrologyc, soilstate_inst, soilhydrology_inst)

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc) 
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer  
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1 
                exit
             end if
          enddo
       end do

      ! Removing the pumped water from the soil column	   
       if (use_pumping == .true.) then
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             g = col%gridcell(c)
                 !!! use analytical expression for aquifer specific yield
                 rous = watsat(c,nlevsoi) &
                      * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
                 rous=max(rous,0.02_r8)

                 pump_tot = - GW_ratio(c) * qflx_irrig(c) * dtime
                 Pump_wa(c) = GW_ratio(c) * qflx_irrig(c)
                 ! if (nstep > 400 .and. c == 838456) write(*,*) 'c,', c, pump_tot, Pump_wa(c) 
                 !!!--  water table is below the soil column  --------------------------------------
                 if(jwt(c) == nlevsoi) then             
                    wa(c)  = wa(c) + pump_tot
                    zwt(c) = zwt(c) - pump_tot/1000._r8/rous
 
                 else                                
                    !!!-- water table within soil layers 1-9  --------------------------------------
                    !!!============================== pump_tot ========================================= 
                    !!!--  Now remove water via pump_tot

                    !!!should never be positive... but include for completeness
                    if(pump_tot > 0.) then !rising water table

                       call endrun(msg="pump_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                    else ! deepening water table
                       do j = jwt(c)+1, nlevsoi
                          !!! use analytical expression for specific yield
                          s_y = watsat(c,j) &
                               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                          s_y=max(s_y,0.02_r8)

                          pump_layer=max(pump_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                          pump_layer=min(pump_layer,0._r8)
                          h2osoi_liq(c,j) = h2osoi_liq(c,j) + pump_layer

                          pump_tot = pump_tot - pump_layer

                          if (pump_tot >= 0.) then 
                             zwt(c) = zwt(c) - pump_layer/s_y/1000._r8
                             exit
                          else
                             zwt(c) = zi(c,j)
                          endif
                       enddo

                       !!!--  remove residual pump_tot  ---------------------------------------------
                       zwt(c) = zwt(c) - pump_tot/1000._r8/rous
                       wa(c)  = wa(c) + pump_tot
                    endif

                    !!!-- recompute jwt  ---------------------------------------------------------
                    !!! allow jwt to equal zero when zwt is in top layer
                    jwt(c) = nlevsoi
                    do j = 1,nlevsoi
                       if(zwt(c) <= zi(c,j)) then
                          jwt(c) = j-1
                          exit
                       end if
                    enddo
                 end if! end of jwt if construct

                 zwt(c) = max(0.0_r8,zwt(c))
                 ! zwt(c) = min(80._r8,zwt(c))    
          end do
       end if
	
       ! To get the weighted average ZWT accross all columns within each cell  
       ! loop over the columns that each processor handles
       ! do c = bounds%begc,bounds%endc
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)
         ! Saving for mpi_allreduce coming up
          g_totCweight(g)    = g_totCweight(g)     + col%wtgcell(c)
          g_cellarea_long(g) = g_cellarea_long(g)  + col%wtgcell(c) * grc%area(g)

          zwt_long(g)          = zwt_long(g)         + col%wtgcell(c) * zwt(c)  ! GW_od goes to history
          AqTransmiss_long(g)  = AqTransmiss_long(g) + col%wtgcell(c) * AqTransmiss(c)

          Pump_wa(c)        = GW_ratio(c) * qflx_irrig(c)
          Pump_wa_long(g)   = Pump_wa_long(g) + col%wtgcell(c) * GW_ratio(c) * qflx_irrig(c)
		  
          GW_ratio_long(g) = GW_ratio_long(g) + col%wtgcell(c) * GW_ratio(c)
          qcharge_long(g)  = qcharge_long(g)  + col%wtgcell(c) * qcharge(c)

          g_watsat_long(g)  = g_watsat_long(g) + col%wtgcell(c) * watsat(c,nlevsoi) 
          g_sucst_long(g)   = g_sucst_long(g)  + col%wtgcell(c) * sucsat(c,nlevsoi) 
          g_bsw_long(g)     = g_bsw_long(g)    + col%wtgcell(c) * bsw(c,nlevsoi)
       end do  

       do g = bounds%begg,bounds%endg
          if (g_totCweight(g) <= 1.000001_r8 .and. g_totCweight(g) > HydroThresh) then

              zwt_long(g)         = zwt_long(g) / g_totCweight(g)
              AqTransmiss_long(g) = AqTransmiss_long(g)/ g_totCweight(g)
              GW_ratio_long(g)    = GW_ratio_long(g) / g_totCweight(g)
              qcharge_long(g)     = qcharge_long(g)/ g_totCweight(g)
              g_watsat_long(g)    = g_watsat_long(g) / g_totCweight(g)
              g_sucst_long(g)     = g_sucst_long(g) / g_totCweight(g)
              g_bsw_long(g)       = g_bsw_long(g) / g_totCweight(g)
          else
              ZeroHydroCell(g) = 1
          end if
       end do


       ! Need action: what happens to the zwt_long of cells g_totCweight(g) < 0.01???????????????
	
       ! ------------------------------------------------
       ! Between the out and in loops
       ! all processors share all the relevant data needed to come up with GW_id
       ! Could the same be done in if (iam == 0) followed by call mpi_bcast?
       ! ------------------------------------------------

       ! Gather zwt_glob from zwt_long, ultimately
       !                  from GW_od

       ! I do MPI_SUM because each processor only has values for the patch/clump its dealing  
       ! with and the rest of the array is zero   
       call mpi_allreduce(zwt_long, zwt_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(AqTransmiss_long, AqTransmiss_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(g_cellarea_long, g_cellarea_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)

       call mpi_allreduce(Pump_wa_long, Pump_wa_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(qcharge_long, qcharge_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 						  

       ! we need to exclude those cells that have no hydrologically active columns  
       call mpi_allreduce(ZeroHydroCell, ZeroHydroCell_glob, ng, &
                          MPI_INTEGER, MPI_SUM, mpicom, ier)

       call mpi_barrier(mpicom,ier)
	   
       ! gathering the information from the neigboring cells.
       ! do  g = 1, ng
       do g = bounds%begg,bounds%endg
          l_edge = 0._r8
          r_edge = 0._r8
          t_edge = 0._r8
          b_edge = 0._r8

          ! The GW lateral flow is ruled by Darcy's
          if (ZeroHydroCell_glob(g)== 0) then

             if (ldecomp%gtoplft(g) <= ng .and. ldecomp%gtoplft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoplft(g)) == 0) then
			 
                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%gtoplft(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtoplft(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%gtoplft(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%gtoplft(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%gtoplft(g)),'___Diag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral

                 l_edge = l_edge + QLateral * sqrt(2._r8)/2._r8  !why positive: because the + direction in eastward
                 t_edge = t_edge - QLateral * sqrt(2._r8)/2._r8  !why negative: because the + direction in upward

             end if

             if (ldecomp%gtop(g) <= ng .and. ldecomp%gtop(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtop(g)) == 0) then

                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%gtop(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtop(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%gtop(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%gtop(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%gtop(g)),'nonDiag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral

                 t_edge = t_edge - QLateral !why negative: because the + direction in upward 

             end if

             if (ldecomp%gtoprgt(g) <= ng .and. ldecomp%gtoprgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoprgt(g)) == 0) then

                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%gtoprgt(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%gtoprgt(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%gtoprgt(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%gtoprgt(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%gtoprgt(g)),'___Diag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral

                 r_edge = r_edge - QLateral * sqrt(2._r8)/2._r8
                 t_edge = t_edge - QLateral * sqrt(2._r8)/2._r8
				 
             end if

             if (ldecomp%grgt(g) <= ng .and. ldecomp%grgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%grgt(g)) == 0) then

                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%grgt(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%grgt(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%grgt(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%grgt(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%grgt(g)),'nonDiag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral				 

                 r_edge = r_edge - QLateral

             end if	

             if (ldecomp%gbotrgt(g) <= ng .and. ldecomp%gbotrgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotrgt(g)) == 0) then

                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%gbotrgt(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbotrgt(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%gbotrgt(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%gbotrgt(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%gbotrgt(g)),'___Diag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral

                 r_edge = r_edge - QLateral * sqrt(2._r8)/2._r8
                 b_edge = b_edge + QLateral * sqrt(2._r8)/2._r8

             end if

             if (ldecomp%gbot(g) <= ng .and. ldecomp%gbot(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbot(g)) == 0) then

                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%gbot(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbot(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%gbot(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%gbot(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%gbot(g)),'nonDiag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral		

                 b_edge = b_edge + QLateral

             end if

             if (ldecomp%gbotlft(g) <= ng .and. ldecomp%gbotlft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotlft(g)) == 0) then

                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%gbotlft(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%gbotlft(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%gbotlft(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%gbotlft(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%gbotlft(g)),'___Diag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral		

                 l_edge = l_edge + QLateral * sqrt(2._r8)/2._r8
                 b_edge = b_edge + QLateral * sqrt(2._r8)/2._r8

             end if

             if (ldecomp%glft(g) <= ng .and. ldecomp%glft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%glft(g)) == 0) then


                 QLateral = this%TheimLateral(AqTransmiss_glob(g),AqTransmiss_glob(ldecomp%glft(g)), &
                                              g_cellarea_glob(g), g_cellarea_glob(ldecomp%glft(g)), & 
                                              qcharge_glob(g), qcharge_glob(ldecomp%glft(g)), &
                                              Pump_wa_glob(g), Pump_wa_glob(ldecomp%glft(g)), &
                                              zwt_glob(g), zwt_glob(ldecomp%glft(g)),'nonDiag')
										  
                 Qn_glob(g) = Qn_glob(g) + QLateral	

                 l_edge = l_edge + QLateral

             end if

             QlatField_north(g) = (t_edge + b_edge)/2._r8
             QlatField_east(g)  = (r_edge + l_edge)/2._r8 

          end if
       end do
   
       ! Checking the lateral water balance (in terms of volume)
       !if (iam == 200) then
       !   dummysum = 0._r8
       !   do  g = 1, ng
       !       dummysum = dummysum + Qn_glob(g)
       !   end do
       !   write(*,*) 'Felfelani: this is the sum of the Latera GW; ', dummysum
       !end if


       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)

          ! going back from Grid level to Column Level
          aRatio  = col%wtgcell(c) / g_totCweight(g)
          colArea = col%wtgcell(c) * grc%area(g) * km2_to_mm2
          Qgw_lateral(c) = Qn_glob(g) * aRatio / colArea  !unit is converted to mm/sec 

          rous = watsat(c,nlevsoi) &
               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
          rous=max(rous,0.02_r8)
          ! if (nstep  > 400 .and. c == 838456) write(*,*) 'c,', c, Qgw_lateral(c)
          if (Qgw_lateral(c) > 0._r8) then

                Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
                if(jwt(c) == nlevsoi) then             
                   wa(c)  = wa(c) + Qgw_lateral_tot
                   zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                else   
                   do j = jwt(c)+1, 1,-1
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=min(Qgw_lateral_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                       Qgw_lateral_layer=max(Qgw_lateral_layer,0._r8)

                       h2osoi_liq(c,j) = h2osoi_liq(c,j) + Qgw_lateral_layer
				
                       if(s_y > 0._r8) zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer
                       if (Qgw_lateral_tot <= 0.) exit
                   enddo
                end if

          else if (Qgw_lateral(c) < 0._r8) then

              Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
              !! --  water table is below the soil column  -------------------------------------- 
              if(jwt(c) == nlevsoi) then             
                 wa(c)  = wa(c) + Qgw_lateral_tot
                 zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
              else                                
                 !! -- water table within soil layers 1-9  -------------------------------------
                 !! ============================== Qgw_lateral_tot ========================================= 
                 !! --  Now remove water via Qgw_lateral_tot

                 !! should never be positive... but include for completeness 
                 if(Qgw_lateral_tot > 0.) then !rising water table

                    call endrun(msg="Qgw_lateral_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                 else ! deepening water table
                    do j = jwt(c)+1, nlevsoi
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=max(Qgw_lateral_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                       Qgw_lateral_layer=min(Qgw_lateral_layer,0._r8)
                       h2osoi_liq(c,j) = h2osoi_liq(c,j) + Qgw_lateral_layer

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer

                       if (Qgw_lateral_tot >= 0.) then 
                         zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8
                          exit
                       else
                          zwt(c) = zi(c,j)
                       endif
                    enddo

                    !! --  remove residual Qgw_lateral_tot  ---------------------------------------------
                    zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                    wa(c) = wa(c) + Qgw_lateral_tot
                 endif

                 !! -- recompute jwt  ---------------------------------------------------------
                 !! allow jwt to equal zero when zwt is in top layer
                 jwt(c) = nlevsoi
                 do j = 1,nlevsoi
                    if(zwt(c) <= zi(c,j)) then
                       jwt(c) = j-1
                       exit
                    end if
                 enddo
              end if! end of jwt if construct

              zwt(c) = max(0.0_r8,zwt(c))
              ! zwt(c) = min(80._r8,zwt(c))
          end if
       end do

       deallocate(Qn_glob)
       deallocate(ZeroHydroCell)
       deallocate(ZeroHydroCell_glob)
       deallocate(zwt_long)
       deallocate(zwt_glob)
       deallocate(g_totCweight)
       deallocate(g_watsat_long)
       deallocate(g_sucst_long)
       deallocate(g_bsw_long)
       deallocate(g_cellarea_long)
       deallocate(g_cellarea_glob)
       deallocate(Pump_wa_long)
       deallocate(Pump_wa_glob)
       deallocate(GW_ratio_long)
       deallocate(AqTransmiss_long)
       deallocate(AqTransmiss_glob)
       deallocate(qcharge_long)
       deallocate(qcharge_glob)

     end associate
  end subroutine UpdateGWFanLatTheimPump

  !-----------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine UpdateGWDefaultPump(this, bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst,waterstate_inst, irrigation_inst)

    ! !DESCRIPTION:
    !   In principle, Theim equation is applied between 
    !   r = r_e and r = dx (center to center of four
    !   surrounding cells) to take into account the GW pumping.
    !
    !   Assumptions: The head across a given radius is constant 
    !  (i.e., cone of depression is circular) 

    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom, MPI_INTEGER
    use decompMod       , only : ldecomp, get_proc_global
    use shr_const_mod   , only : SHR_CONST_PI
    use GridcellType    , only : grc
    use clm_time_manager, only : get_step_size, get_curr_date, get_nstep
    use landunit_varcon , only : istwet, istsoil, istice_mec, istcrop

    ! !ARGUMENTS:
    class(groundwater_type)  , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(waterstate_type)    , intent(inout) :: waterstate_inst
    type(irrigation_type)    , intent(in)    :: irrigation_inst	

    ! !LOCAL VARIABLES:	
 
    character(len=32) :: subname = 'GroundwaterMod' ! subroutine name 

    real(r8) :: rous, aRatio, colArea              ! aquifer yield (-); area ratio of center/neighbor
    real(r8) :: s_y, dummysum, dummysum2
    real(r8) :: pump_tot, pump_layer
    integer  :: jwt(bounds%begc:bounds%endc)       ! index of the soil layer right above the water table (-)
    real(r8) :: dtime                              ! land model time step (sec)

	
    integer :: ng, nl, nc, np, nCohorts            ! total number of grid cells,landunits,columns,patches
    integer :: g, c                                ! patch, gridcell, column indices
    integer :: ier                                 ! error code
    integer :: j,fc,i

    integer :: begg, endg       ! beginning and ending gridcell index of current proc
    integer :: begl, endl       ! beginning and ending landunit index of current proc
    integer :: begc, endc       ! beginning and ending column index of current proc
    integer :: begp, endp       ! beginning and ending pft index of current proc

    integer :: year       ! year (0, ...) for nstep
    integer :: month      ! month (1, ..., 12) for nstep
    integer :: day        ! day of month (1, ..., 31) for nstep
    integer :: secs       ! seconds into current date for nstep
    integer :: nstep

    ! Conversion factors

    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: mm_to_m     = 1.e-3_r8
    real(r8), parameter :: m_to_mm     = 1.e3_r8
    real(r8), parameter :: HydroThresh = 0.1_r8
    !-----------------------------------------------------------------------
     associate(&                               
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
          GW_ratio           =>    col%GW_ratio                          , & ! Input:  [real(r8) (:)   ]  USGS GW ratio as irrigation source                                                
                                              
          qflx_irrig         =>    irrigation_inst%qflx_irrig_col        , & ! irrigation flux (mm H2O /s)

          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)

          zwt                =>    soilhydrology_inst%zwt_col                , & ! Input and Output: [real(r8) (:)   ]  water table depth (m)
          wa                 =>    soilhydrology_inst%wa_col                 , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)
          qcharge            =>    soilhydrology_inst%qcharge_col            , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
          Qgw_lateral        =>    soilhydrology_inst%Qgw_lateral_col        , & ! Output: [real(r8) (:)   ]  GW lateral flow (mm/s)
          AqTransmiss        =>    soilhydrology_inst%AqTransmiss_col        , & ! Output: [real(r8) (:)   ]  Aquifer Transmissivity(mm2/s)
          Pump_wa            =>    soilhydrology_inst%Pump_wa_col            , & ! Output: [real(r8) (:)   ]  Pumped Water from the aquifer(mm/s)
          QlatField_north    =>    soilhydrology_inst%QlatField_northing_grc , & !  Output: [real(r8) (:)   ] Northward lateral GW flow
          QlatField_east     =>    soilhydrology_inst%QlatField_easting_grc  , & !  Output: [real(r8) (:)   ] Eastward lateral GW flow

          h2osoi_liq         =>    waterstate_inst%h2osoi_liq_col              & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          )
       !-----------------------------------------------
       dtime = get_step_size()
       nstep = get_nstep()
       call get_curr_date (year, month, day, secs)
       call get_proc_global(ng=ng, nl=nl, nc=nc, np=np, nCohorts=nCohorts)


       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc) 
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer  
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1 
                exit
             end if
          enddo
       end do

      ! Removing the pumped water from the soil column	   
       if (use_pumping == .true.) then
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             g = col%gridcell(c)
                 !!! use analytical expression for aquifer specific yield
                 rous = watsat(c,nlevsoi) &
                      * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
                 rous=max(rous,0.02_r8)

                 pump_tot = - GW_ratio(c) * qflx_irrig(c) * dtime
                 Pump_wa(c) = GW_ratio(c) * qflx_irrig(c)
                 ! if (nstep > 400 .and. c == 838456) write(*,*) 'c,', c, pump_tot, Pump_wa(c) 
                 !!!--  water table is below the soil column  --------------------------------------
                 if(jwt(c) == nlevsoi) then             
                    wa(c)  = wa(c) + pump_tot
                    zwt(c) = zwt(c) - pump_tot/1000._r8/rous
 
                 else                                
                    !!!-- water table within soil layers 1-9  --------------------------------------
                    !!!============================== pump_tot ========================================= 
                    !!!--  Now remove water via pump_tot

                    !!!should never be positive... but include for completeness
                    if(pump_tot > 0.) then !rising water table

                       call endrun(msg="pump_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                    else ! deepening water table
                       do j = jwt(c)+1, nlevsoi
                          !!! use analytical expression for specific yield
                          s_y = watsat(c,j) &
                               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                          s_y=max(s_y,0.02_r8)

                          pump_layer=max(pump_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                          pump_layer=min(pump_layer,0._r8)
                          h2osoi_liq(c,j) = h2osoi_liq(c,j) + pump_layer

                          pump_tot = pump_tot - pump_layer

                          if (pump_tot >= 0.) then 
                             zwt(c) = zwt(c) - pump_layer/s_y/1000._r8
                             exit
                          else
                             zwt(c) = zi(c,j)
                          endif
                       enddo

                       !!!--  remove residual pump_tot  ---------------------------------------------
                       zwt(c) = zwt(c) - pump_tot/1000._r8/rous
                       wa(c)  = wa(c) + pump_tot
                    endif

                    !!!-- recompute jwt  ---------------------------------------------------------
                    !!! allow jwt to equal zero when zwt is in top layer
                    jwt(c) = nlevsoi
                    do j = 1,nlevsoi
                       if(zwt(c) <= zi(c,j)) then
                          jwt(c) = j-1
                          exit
                       end if
                    enddo
                 end if! end of jwt if construct

                 zwt(c) = max(0.0_r8,zwt(c))
                 ! zwt(c) = min(80._r8,zwt(c))    
          end do
       end if

     end associate
  end subroutine UpdateGWDefaultPump
  !-----------------------------------------------------------------------

  function TransmissivityFromFan(this, bounds, num_hydrologyc, filter_hydrologyc, &
                                 soilstate_inst, soilhydrology_inst) &
    result(Transmiss)

    ! !DESCRIPTION:
    !  Calculating the transmissivity based on the Fan et al. (2007)
    !  and Zeng et al. (2016)

    use decompMod       , only : ldecomp

    ! !ARGUMENTS:
    real(r8),allocatable                     :: Transmiss(:) ! Func Result: Transmissivity of the aquifer (mm^2/s)
    class(groundwater_type)  , intent(in)    :: this
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc                     ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:)               ! column filter for soil points
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst

    ! !LOCAL VARIABLES:
    integer  :: c,j,fc,i,g                                   ! indices
    real(r8) :: e_folding_length
    real(r8) :: beta_rad                                     ! terrain slope (rad)
    integer  :: jwt2(bounds%begc:bounds%endc)                ! index of the soil layer right above the water table (-)

    real(r8), parameter :: m_to_mm = 1.e3_r8	     
    character(len=8), parameter :: bedrockORregolit = 'regolith'  ! Choose between 'regolith' and 'bedrock_'


    associate(                                                             &                              
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)  

          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          hk_l               =>    soilstate_inst%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)
          clayP              =>    soilstate_inst%cellclay_col           , & ! Input:  [real(r8) (:,:) ] percent clay (0 < ~ < 100)-----> It is not a fraction!!!
          zwt                =>    soilhydrology_inst%zwt_col              & ! Input: [real(r8) (:)   ]  water table depth (m)
          )
    !-----------------------------------------------------------------------
    allocate(Transmiss(bounds%begc:bounds%endc))
    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt2(c) = 100 ! arbitrarily 100 for water table below the soil column 
       ! allow jwt2 to equal zero when zwt is in top layer
       do j = 1,nlevsoi
          if(zwt(c) <= zi(c,j)) then
             jwt2(c) = j
             exit
          end if
       enddo
    end do

    Transmiss(:) = 0._r8
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       g = col%gridcell(c)


       if (bedrockORregolit == 'bedrock_') then 
          !!!!!!! Zeng et al.(2018): e_folding_length (m) calculations for bedrock	 
          beta_rad = (rpi/180._r8) * col%topo_slope(c)
          if (beta_rad <= 0.16) then
              e_folding_length = 20._r8/(1._r8 + 125._r8*beta_rad)
          else if (beta_rad > 0.16) then
              e_folding_length = 1._r8
          end if
       else if (bedrockORregolit == 'regolith') then 
          !!!!!!! Fan et al.(2007): e_folding_length (m) calculations for regolith
          beta_rad = (rpi/180._r8) * col%topo_slope(c)
          if (beta_rad <= 0.16) then
              e_folding_length = 120._r8/(1._r8 + 150._r8*beta_rad)
          else if (beta_rad > 0.16) then
              e_folding_length = 5._r8
          end if

       end if

       !!!! 1.5-m Trannsmissible Tickness -- start
       if (jwt2(c) < 10) then
          Transmiss(c) = clayP(c,jwt2(c)) * hksat(c,jwt2(c)) * (zi(c,jwt2(c))-zwt(c)) * m_to_mm
          do j = jwt2(c)+1,10
             Transmiss(c) = Transmiss(c) + clayP(c,j) * hksat(c,j) * dz(c,j) * m_to_mm
          end do
          Transmiss(c) = Transmiss(c) + clayP(c,10) * hksat(c,10) * col%FEDEPTH(c) * m_to_mm

       else if (jwt2(c) .eq. 10) then
          Transmiss(c) = clayP(c,10) * hksat(c,10) * (zi(c,10)-zwt(c)) * m_to_mm &
                         + clayP(c,10) * hksat(c,10) * col%FEDEPTH(c) * m_to_mm

       else if (jwt2(c) > 10) then
          Transmiss(c) = clayP(c,10) * hksat(c,10) * col%FEDEPTH(c) * m_to_mm &
                         * exp((zi(c,10)-zwt(c))/col%FEDEPTH(c))
       end if
	   !Transmiss(c) = 0.1_r8 * Transmiss(c)
       !!!! 1.5-m Trannsmissible Tickness -- end

       !!!! Original formulation -- start
       !if (jwt2(c) < nlevsoi) then
       !   Transmiss(c) = clayP(c,jwt2(c)) * hksat(c,jwt2(c)) * (zi(c,jwt2(c))-zwt(c)) * m_to_mm
       !   do j = jwt2(c)+1,nlevsoi
       !      Transmiss(c) = Transmiss(c) + clayP(c,j) * hksat(c,j) * dz(c,j) * m_to_mm
       !   end do
       !   Transmiss(c) = Transmiss(c) + clayP(c,nlevsoi) * hksat(c,nlevsoi) * e_folding_length * m_to_mm

       !else if (jwt2(c) .eq. nlevsoi) then
       !   Transmiss(c) = clayP(c,nlevsoi) * hksat(c,nlevsoi) * (zi(c,nlevsoi)-zwt(c)) * m_to_mm &
       !                  + clayP(c,nlevsoi) * hksat(c,nlevsoi) * e_folding_length * m_to_mm

       !else if (jwt2(c) .eq. 100) then
       !   Transmiss(c) = clayP(c,nlevsoi) * hksat(c,nlevsoi) * e_folding_length * m_to_mm &
       !                  * exp((zi(c,nlevsoi)-zwt(c))/e_folding_length)
       !end if
       !!!! Original formulation -- end

       !!!! If considering hydraulic conductivity instead of saturated hydraulic conductivity

       ! if (jwt2(c) < nlevsoi) then
          ! Transmiss(c) = clayP(c,jwt2(c)) * hk_l(c,jwt2(c)) * (zi(c,jwt2(c))-zwt(c)) * m_to_mm
          ! do j = jwt2(c)+1,nlevsoi
             ! Transmiss(c) = Transmiss(c) + clayP(c,j) * hk_l(c,j) * dz(c,j) * m_to_mm
          ! end do
          ! Transmiss(c) = Transmiss(c) + clayP(c,nlevsoi) * hk_l(c,nlevsoi) * e_folding_length * m_to_mm

       ! else if (jwt2(c) .eq. nlevsoi) then
          ! Transmiss(c) = clayP(c,nlevsoi) * hk_l(c,nlevsoi) * (zi(c,nlevsoi)-zwt(c)) * m_to_mm &
                         ! + clayP(c,nlevsoi) * hk_l(c,nlevsoi) * e_folding_length * m_to_mm

       ! else if (jwt2(c) .eq. 100) then
          ! Transmiss(c) = clayP(c,nlevsoi) * hk_l(c,nlevsoi) * e_folding_length * m_to_mm &
                         ! * exp((zi(c,nlevsoi)-zwt(c))/e_folding_length)
       ! end if

	   ! if (grc%latdeg(g) > 37.270 .and. grc%latdeg(g) < 37.280 .and. grc%londeg(g) > 257.270 .and. grc%londeg(g) < 257.280) then
	       ! write(*,*)  '--------------------------------------------------------------'   
           ! write(*, '(a60, i8, i8, i8, f20.6, f20.6, f20.6)')        'Line 1: iam, g, c, Nneighbor, grc%latdeg, grc%londeg  ', iam, g, c, ldecomp%gneighbors(g), grc%latdeg(g), grc%londeg(g)
		   ! write(*, '(a60, i8, i8, ES15.5, ES15.5, ES15.5, ES15.5, ES15.5)') 'Line 2: g, c, clayP, hksat, e_fold, zwt, trans        ', g, c, clayP(c,nlevsoi), hksat(c,nlevsoi), e_folding_length, zwt(c), Transmiss(c)
		   ! write(*, '(a60, i8, i8, ES15.5, ES15.5, ES15.5, ES15.5, ES15.5)') 'Line 3: g, c, beta_rad, col%slop, grc%slop, rpi, brad ', g, c, beta_rad, col%topo_slope(c), grc%slopelev(g), rpi, (rpi/180._r8) * col%topo_slope(c)

	   ! end if
    end do

    end associate
  end function TransmissivityFromFan

  !-----------------------------------------------------------------------

  function FanLatVal(trans_c, trans_n, a_c, a_n, &
                           slope_c, slope_n, zwt_c, zwt_n, gUP, Diagonal) &
    result(latval)

    ! !DESCRIPTION:
    !  Calculating the lateral flow value based on the Fan et al. (2007)
    !  and Zeng et al. (2016)

    ! !ARGUMENTS:
    real(r8)                                 :: latval ! Func Result: lateral flow value [mm^3/s]
    ! class(groundwater_type)  , intent(in)    :: this
    real(r8)                 , intent(in)    :: trans_c, trans_n  ! Transmissivity of the center and neighbor cell [mm^2/s]
    real(r8)                 , intent(in)    :: a_c, a_n          ! Area of the center and neighbor cell [km^2]
    real(r8)                 , intent(in)    :: slope_c, slope_n  ! Slope of the center and neighbor cell [degrees]
    real(r8)                 , intent(in)    :: zwt_c, zwt_n      ! Depth to WT of the center and neighbor cell [m]
    real(r8)                 , intent(in)    :: gUP               ! Uphill/downhill
    character(len=7)         , intent(in)    :: Diagonal

    ! !LOCAL VARIABLES:
    integer             :: c,j,fc,i,g, nstep                      ! indices
    real(r8)            :: TransMean_
    real(r8)            :: deltaxMean_
    real(r8)            :: widMean_
    real(r8)            :: lenMean_
    real(r8)            :: slopHMean_

    real(r8), parameter :: m_to_mm     = 1.e3_r8	     
    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: km_to_m     = 1.e3_r8

    !-----------------------------------------------------------------------
    TransMean_ = (trans_c + trans_n)/2._r8

    if (Diagonal == '___Diag') then

       lenMean_ = (sqrt(a_c) + sqrt(a_n)) * km_to_mm * sqrt(2._r8) / 2._r8
       deltaxMean_ = (sqrt(a_c) + sqrt(a_n)) * km_to_mm / 2._r8
       widMean_ = deltaxMean_ * sqrt(0.5_r8 * tan(rpi/8._r8))
       slopHMean_ = (sqrt(a_c)*tan(rpi/180._r8*slope_c) + &
                     sqrt(a_n)*tan(rpi/180._r8*slope_n)) * km_to_m * sqrt(2._r8) / 2._r8

    else if (Diagonal == 'nonDiag') then

       lenMean_ = (sqrt(a_c) + sqrt(a_n)) * km_to_mm / 2._r8
       widMean_ = lenMean_ * sqrt(0.5_r8 * tan(rpi/8._r8))
       slopHMean_ = (sqrt(a_c)*tan(rpi/180._r8*slope_c) + &
                     sqrt(a_n)*tan(rpi/180._r8*slope_n)) * km_to_m / 2._r8

    end if

    latval = widMean_ * TransMean_ * (abs(gUP) * (zwt_c - zwt_n) + gUP * slopHMean_) * m_to_mm / lenMean_

  end function FanLatVal

  !-----------------------------------------------------------------------

  function FanLatVal_hgt(trans_c, trans_n, a_c, a_n, &
                           hgt_c, hgt_n, zwt_c, zwt_n, gUP, Diagonal) &
    result(latval)

    ! !DESCRIPTION:
    !  Calculating the lateral flow value based on the Fan et al. (2007)
    !  and Zeng et al. (2016)

    ! !ARGUMENTS:
    real(r8)                                 :: latval ! Func Result: lateral flow value [mm^3/s]
    ! class(groundwater_type)  , intent(in)    :: this
    real(r8)                 , intent(in)    :: trans_c, trans_n  ! Transmissivity of the center and neighbor cell [mm^2/s]
    real(r8)                 , intent(in)    :: a_c, a_n          ! Area of the center and neighbor cell [km^2]
    real(r8)                 , intent(in)    :: hgt_c, hgt_n      ! Hight of the center and neighbor cell [m]
    real(r8)                 , intent(in)    :: zwt_c, zwt_n      ! Depth to WT of the center and neighbor cell [m]
    real(r8)                 , intent(in)    :: gUP               ! Uphill/downhill
    character(len=7)         , intent(in)    :: Diagonal

    ! !LOCAL VARIABLES:
    integer             :: c,j,fc,i,g, nstep                      ! indices
    real(r8)            :: TransMean_
    real(r8)            :: deltaxMean_
    real(r8)            :: widMean_
    real(r8)            :: lenMean_
    real(r8)            :: slopHMean_

    real(r8), parameter :: m_to_mm     = 1.e3_r8	     
    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: km_to_m     = 1.e3_r8

    !-----------------------------------------------------------------------
    TransMean_ = (trans_c + trans_n)/2._r8

    if (Diagonal == '___Diag') then

       lenMean_ = (sqrt(a_c) + sqrt(a_n)) * km_to_mm * sqrt(2._r8) / 2._r8
       deltaxMean_ = (sqrt(a_c) + sqrt(a_n)) * km_to_mm / 2._r8
       widMean_ = deltaxMean_ * sqrt(0.5_r8 * tan(rpi/8._r8))


    else if (Diagonal == 'nonDiag') then

       lenMean_ = (sqrt(a_c) + sqrt(a_n)) * km_to_mm / 2._r8
       widMean_ = lenMean_ * sqrt(0.5_r8 * tan(rpi/8._r8))


    end if

    latval = widMean_ * TransMean_ * (abs(gUP) * ((hgt_n - zwt_n)-(hgt_c - zwt_c))) * m_to_mm / lenMean_

  end function FanLatVal_hgt

  !-----------------------------------------------------------------------

  function TheimLateral(this, AqTransmiss_cent, AqTransmiss_neig, gcellarea_cent, gcellarea_neig, &
                              qcharge_cent, qcharge_neig, PumpWa_cent, PumpWa_neig, zwt_cent, zwt_neig, Diagonal) &
    result(QLateral)

    ! !DESCRIPTION:
    !  Calculating the lateral flow based on the Fan et al. (2007)
    ! and also Theim Theory

    ! !ARGUMENTS:
    real(r8)                                 :: QLateral                            ! Fuction output: Lateral Flow to Cell i,j (mm3/sec)
    class(groundwater_type)  , intent(in)    :: this
    real(r8)                 , intent(in)    :: AqTransmiss_cent, AqTransmiss_neig  ! number of column soil points in column filter
    real(r8)                 , intent(in)    :: gcellarea_cent, gcellarea_neig      ! column filter for soil points
    real(r8)                 , intent(in)    :: qcharge_cent, qcharge_neig
    real(r8)                 , intent(in)    :: PumpWa_cent, PumpWa_neig
    real(r8)                 , intent(in)    :: zwt_cent, zwt_neig 
    character(len=7)         , intent(in)    :: Diagonal

    ! !LOCAL VARIABLES:
    integer             :: c,j,fc,i,g, nstep                                               ! indices
    real(r8)            :: AqTransmissMean
    real(r8)            :: qchargeMean
    real(r8)            :: deltaxMean
    real(r8)            :: widMean
    real(r8)            :: cellareaMean
    real(r8)            :: lenMean
    real(r8)            :: QLatDummy, QLatDummy2
    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: mm_to_m     = 1.e-3_r8
    real(r8), parameter :: m_to_mm     = 1.e3_r8
    real(r8), parameter :: HydroThresh = 0.1_r8

    !----------------------------------------------------------------------- 
    nstep = get_nstep()
    AqTransmissMean = (AqTransmiss_cent + AqTransmiss_neig)/2._r8
    qchargeMean     = (qcharge_cent + qcharge_neig)/2._r8
    cellareaMean    = (gcellarea_cent + gcellarea_neig)/2._r8

    if (Diagonal == '___Diag') then

       lenMean = (sqrt(gcellarea_cent) + sqrt(gcellarea_neig)) * km_to_mm * sqrt(2._r8) / 2._r8
       deltaxMean = (sqrt(gcellarea_cent) + sqrt(gcellarea_neig)) * km_to_mm / 2._r8
       widMean = deltaxMean * sqrt(0.5_r8 * tan(rpi/8._r8))

    else if (Diagonal == 'nonDiag') then

       lenMean = (sqrt(gcellarea_cent) + sqrt(gcellarea_neig)) * km_to_mm / 2._r8
       widMean = lenMean * sqrt(0.5_r8 * tan(rpi/8._r8))

    end if

    ! Theim for the center cell; Fan lateral for the neighbor
    if (PumpWa_cent > 0._r8 .and. PumpWa_neig == 0._r8 .and. zwt_cent > zwt_neig) then

        !QLateral = SHR_CONST_PI * AqTransmissMean * (zwt_cent - zwt_neig) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8) + &
        !           qchargeMean * SHR_CONST_PI * cellareaMean * km2_to_mm2 * (1._r8 - 0.17803_r8**2) / 16._r8 / log(1._r8 / 0.17803_r8) - &
        !           qchargeMean * SHR_CONST_PI * (0.17803_r8 * sqrt(cellareaMean) * km_to_mm)**2 / 8._r8

        QLateral = SHR_CONST_PI * AqTransmissMean * (zwt_cent - zwt_neig) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8)

    ! Theim for the neighbor cell; Fan lateral for the center
    else if (PumpWa_cent == 0._r8 .and. PumpWa_neig > 0._r8 .and. zwt_cent < zwt_neig) then 

        !QLateral = -(SHR_CONST_PI * AqTransmissMean * (zwt_neig - zwt_cent) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8) + &
        !           qchargeMean * SHR_CONST_PI * cellareaMean * km2_to_mm2 * (1._r8 - 0.17803_r8**2) / 16._r8 / log(1._r8 / 0.17803_r8) - &
        !           qchargeMean * SHR_CONST_PI * (0.17803_r8 * sqrt(cellareaMean) * km_to_mm)**2 / 8._r8) 

        QLateral = -(SHR_CONST_PI * AqTransmissMean * (zwt_neig - zwt_cent) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8))

    ! Theim for the neighbor cell
    else if (PumpWa_cent > 0._r8 .and. PumpWa_neig > 0._r8 .and. zwt_cent < zwt_neig) then

        !QLateral = -(SHR_CONST_PI * AqTransmissMean * (zwt_neig - zwt_cent) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8) + &
        !           qchargeMean * SHR_CONST_PI * cellareaMean * km2_to_mm2 * (1._r8 - 0.17803_r8**2) / 16._r8 / log(1._r8 / 0.17803_r8) - &
        !           qchargeMean * SHR_CONST_PI * (0.17803_r8 * sqrt(cellareaMean) * km_to_mm)**2 / 8._r8)  


        QLateral =-(SHR_CONST_PI * AqTransmissMean * (zwt_neig - zwt_cent) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8))

    ! Theim for the center cell
    else if (PumpWa_cent > 0._r8 .and. PumpWa_neig > 0._r8 .and. zwt_cent > zwt_neig) then

        !QLateral = SHR_CONST_PI * AqTransmissMean * (zwt_cent - zwt_neig) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8) + &
        !           qchargeMean * SHR_CONST_PI * cellareaMean * km2_to_mm2 * (1._r8 - 0.17803_r8**2) / 16._r8 / log(1._r8 / 0.17803_r8) - &
        !           qchargeMean * SHR_CONST_PI * (0.17803_r8 * sqrt(cellareaMean) * km_to_mm)**2 / 8._r8 


        QLateral = SHR_CONST_PI * AqTransmissMean * (zwt_cent - zwt_neig) * m_to_mm / 4._r8 / log(1._r8 / 0.17803_r8)

    else

        QLateral = widMean * AqTransmissMean * (zwt_cent - zwt_neig) * m_to_mm / lenMean

    end if
    QLatDummy =  widMean * AqTransmissMean * (zwt_cent - zwt_neig) * m_to_mm / lenMean

    !if (QLatDummy .ne. QLateral) write(*,*) "QLateral_Fan, QLateral_Theim", QLatDummy, QLateral
    !if (nstep == 300 .and. (QLatDummy .ne. QLateral)) write(*,*) "QLateral_Fan, QLateral_Theim w/wo E", QLatDummy, QLateral, QLatDummy2
    ! write(*,*) 'QLatDummy, QLateral',QLatDummy, QLateral
    ! write(*,*) 'Trans_cent, Trans_neig',AqTransmiss_cent, AqTransmiss_neig 
    ! write(*,*) 'gcellarea_cent, gcellarea_neig',gcellarea_cent, gcellarea_neig
    ! write(*,*) 'qcharge_cent, qcharge_neig', qcharge_cent, qcharge_neig
    ! write(*,*) 'PumpWa_cent, PumpWa_neig', PumpWa_cent, PumpWa_neig
    ! write(*,*) 'zwt_cent, zwt_neig', zwt_cent, zwt_neig
    ! write(*,*) 'lenMean,widMean',lenMean,widMean
    ! write(*,*) ''
  end function TheimLateral
  
  
end module GroundwaterMod