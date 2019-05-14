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
  use clm_varctl        , only : iulog
  use clm_varcon        , only : isecspday, degpsec, denh2o, spval, namec, rpi
  use clm_varpar        , only : nlevsoi, nlevgrnd
  use clm_time_manager  , only : get_step_size
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  use GridcellType      , only : grc                
  use ColumnType        , only : col                
  use PatchType         , only : patch                
  use subgridAveMod     , only : p2c, c2g
  use filterColMod      , only : filter_col_type, col_filter_from_logical_array
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use IrrigationMod     , only : irrigation_type
  use spmdMod           , only : iam  ! FFelfelani: to get processor number
  ! !PUBLIC TYPES:
  implicit none
  private
  


  ! Public routines
  public :: UpdateGWTheim
  ! public :: UpdateGWFanLateral

  ! Private routines
  private :: TransmissivityFromFan
  
contains

  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================
  
  !------------------------------------------------------------------------
  subroutine UpdateGWTheim(bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst,irrigation_inst)

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
    use shr_const_mod   , only : SHR_CONST_PI
	! get_proc_global: total gridcells, landunits, columns, patchs across all processors
    ! get_proc_global(ng, nl, nc, np, nCohorts) --->  all the arguments are intent(out)
	
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(irrigation_type)    , intent(in)    :: irrigation_inst	
	
    ! !LOCAL VARIABLES:	
  
    real(r8), pointer :: neighbors_count(:)  ! complete grid cell array of count
    real(r8), pointer :: GW_sum_glob(:)        ! First term of Theim Theory
    real(r8), pointer :: GW_out_long(:)        ! GW_out array for all grid cells
    real(r8), pointer :: GW_out_glob(:)        ! complete grid cell array of GW_out
    real(r8), pointer :: lodgepole_wtgcell_long(:)  ! same pair for...
    real(r8), pointer :: lodgepole_wtgcell_glob(:)  ! ...lodgepole_wtgcell
    real(r8), pointer :: AqTransmiss(:)  ! Transmissivity of the aquifer (mm^2/s)
    real(r8)          :: dtime           ! land model time step (sec)
	
    integer :: ng          ! total number of grid cells
    integer :: p, g ,c     ! patch, gridcell, column indices
    integer :: g_in, g_out  ! gridcell indices in/out cells
    integer :: ier         ! error code
    integer :: yr          ! year
    integer :: mon         ! month
    integer :: day         ! day
    integer :: tod         ! seconds
    integer :: j,fc,i                            
                           ! id ---> in-dispersing
						   ! od ---> out-dispersing

    ! Conversion factors
    real(r8), parameter :: km_to_mm   = 1.e6_r8
    real(r8), parameter :: km2_to_mm2 = 1.e12_r8						   
    real(r8), parameter :: mm_to_m = 1.e-3_r8	
	!-----------------------------------------------------------------------
     associate(& 
          !z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          !zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
          GW_ratio           =>    col%GW_ratio                          , & ! Input:  [real(r8) (:)   ]  USGS GW ratio as irrigation source                     
                                               
          qflx_irrig         =>    irrigation_inst%qflx_irrig_col        , & ! irrigation flux (mm H2O /s)
		  
          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
          eff_porosity       =>    soilstate_inst%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         

          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input and Output: [real(r8) (:)   ]  water table depth (m)                             
          zwt_perched        =>    soilhydrology_inst%zwt_perched_col    , & ! Output: [real(r8) (:)   ]  perched water table depth (m)                     
          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Output: [real(r8) (:)   ]  frost table depth (m)                             
          wa                 =>    soilhydrology_inst%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
          qcharge            =>    soilhydrology_inst%qcharge_col          & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
          )
       !-----------------------------------------------
	   
       allocate(AqTransmiss(bounds%begg:bounds%endg))       		  
 
       ! Initialize for the mpi_allreduce located between the out and
       ! in loops
       call get_proc_global(ng=ng)
       ! Variables to gather from all PEs, while between the out and
       ! in loops
       allocate(GW_sum_glob(ng))
       allocate(GW_out_long(ng))
       allocate(GW_out_glob(ng))  

       ! Initialize to 0 so as to MPI_SUM zeros in all PEs but one per grid cell
       GW_out_long(:) = 0._r8

       ! Initialize to 1e36 to help make errors stand out
       GW_out_glob(:) = 1.e36_r8

       ! Initialize to 0 to sum correctly
       GW_sum_glob(:) = 0._r8

       do c = bounds%begc,bounds%endc  ! first p loop: out
          g = col%gridcell(c)
		  
		  ! Saving for mpi_allreduce coming up between the do p loops
          GW_out_long(g) = zwt(c)  ! GW_od goes to history

       end do  ! first p loop: out	   

	   
       ! ------------------------------------------------
       ! Between the out and in loops
       ! all processors share all the relevant data needed to come up with GW_id
       ! Could the same be done in if (iam == 0) followed by call mpi_bcast?
       ! ------------------------------------------------

       ! Gather GW_out_glob from GW_out_long, ultimately
       !                  from GW_od

       ! I do MPI_SUM because each processor only has values for the patch/clump its dealing
       ! with and the rest of the array is zero   
       call mpi_allreduce(GW_out_long, GW_out_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 


       ! From out GW (GW_out_glob) get in GW
       ! (GW_id) by finding nearest neighbors with the help of ixy and jxy

       allocate(neighbors_count(ng))
       neighbors_count = 0._r8  ! initialize counter vector

       ! Loop to find all out GW' neighbors and
       ! sum each out cell's neighbors' lodgepole pine weights so
       ! that neighbors_count is the weighted sum of grid cells with lodgepole
       ! pine > 0
       do g_out = 1, ng
          if (GW_out_glob(g_out) > 0._r8) then
             do g_in = 1, ng
                ! identify neighbors with the ixy, jxy indices of grid cells
                if (ldecomp%ixy(g_out) == ldecomp%ixy(g_in) - 1 .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in)     .or. &

                    ldecomp%ixy(g_out) == ldecomp%ixy(g_in)     .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in) + 1 .or. &

                    ldecomp%ixy(g_out) == ldecomp%ixy(g_in) + 1 .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in)     .or. &

                    ldecomp%ixy(g_out) == ldecomp%ixy(g_in)     .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in) - 1) then


                    neighbors_count(g_out) = neighbors_count(g_out) + 1
                end if  ! find surrounding neighbors
             end do  ! g_in loop
          end if  ! GW_od > 0
          ! GW_id gets GW from its neighbors
          ! Disperse GW from out cells to in cells
          ! Weight out GW by lodgepole pine weights divided by
          ! the neighbors_count sum
          if (neighbors_count(g_out) > 0._r8) then
             do g_in = 1, ng
                if (ldecomp%ixy(g_out) == ldecomp%ixy(g_in) - 1 .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in)     .or. &

                    ldecomp%ixy(g_out) == ldecomp%ixy(g_in)     .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in) + 1 .or. &


                    ldecomp%ixy(g_out) == ldecomp%ixy(g_in) + 1 .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in)     .or. &


                    ldecomp%ixy(g_out) == ldecomp%ixy(g_in)     .and.  &
                    ldecomp%jxy(g_out) == ldecomp%jxy(g_in) - 1) then

                    
                    GW_sum_glob(g_out) = GW_sum_glob(g_out) + GW_out_glob(g_in) / neighbors_count(g_out)				

                end if  ! find surrounding neighbors
             end do  ! g_in
          end if  ! neighbors_count > 0      
       end do  ! g_out loop

       AqTransmiss = TransmissivityFromFan(bounds, num_hydrologyc, filter_hydrologyc, soilstate_inst, soilhydrology_inst)	   

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)

          if (grc%latdeg(g) < 34.03 .and. grc%latdeg(g) > 34.02 .and. grc%londeg(g) < 257.03 .and. grc%londeg(g) > 257.02) then
              write(*,*) 'Felfelani: Processor Num: g, lat(g), lon(g), zwt(c)', iam, g, grc%latdeg(j), grc%londeg(j), zwt(c)
          end if

 
          zwt(c) = GW_sum_glob(g) + &
                  ((GW_ratio(c) * qflx_irrig(c) * grc%area(g) * km2_to_mm2) / (2 * SHR_CONST_PI * AqTransmiss(c)) + &
                  (qcharge(c) * (0.208_r8 * sqrt(grc%area(g)) * km_to_mm)**2/(2 * AqTransmiss(c)))) * log(1/0.208_r8) * mm_to_m - &
                  (qcharge(c) * grc%area(g) * km2_to_mm2 * (1-0.208_r8**2)/(4 * AqTransmiss(c))) * mm_to_m	

          wa(c)  = wa(c) - GW_ratio(c) * qflx_irrig(c) * dtime

          if (grc%latdeg(g) < 34.03 .and. grc%latdeg(g) > 34.02 .and. grc%londeg(g) < 257.03 .and. grc%londeg(g) > 257.02) then
              write(*,*) 'Felfelani: Processor Num: g, lat(g), lon(g), zwt(c)', iam, g, grc%latdeg(j), grc%londeg(j), zwt(c)
          end if

       end do       
	   
	   
       deallocate(neighbors_count)
       deallocate(GW_sum_glob)
       deallocate(GW_out_long)
       deallocate(GW_out_glob)
       deallocate(AqTransmiss)

     end associate
  end subroutine UpdateGWTheim

  !-----------------------------------------------------------------------

  !subroutine UpdateGWFanLateral(this, bounds, &)

  !end subroutine UpdateGWFanLateral

  !-----------------------------------------------------------------------

  function TransmissivityFromFan(bounds, num_hydrologyc, filter_hydrologyc, &
                                 soilstate_inst, soilhydrology_inst) &
    result(Transmiss)

    ! !DESCRIPTION:
    !  Calculating the transmissivity based on the Fan et al. (2007)
    !  and Zeng et al. (2016)

    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst

    ! !LOCAL VARIABLES:
    integer  :: c,j,fc,i,g,l                                 ! indices
    real(r8) :: e_folding_length
    real(r8) :: Transmiss(bounds%begc:bounds%endc)           ! Transmissivity of the aquifer (mm^2/s)
    real(r8) :: beta_rad                                     ! terrain slope (rad)
    integer  :: jwt2(bounds%begc:bounds%endc)                ! index of the soil layer right above the water table (-)
     
    associate(                                                            & 	
          !z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)  

          hk_l               =>    soilstate_inst%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)                    
          clayP              =>    soilstate_inst%cellclay_col           , & ! Input:  [real(r8) (:,:) ] percent clay (0 < ~ < 100)-----> Make sure it is not a fraction!!!
          zwt                =>    soilhydrology_inst%zwt_col              & ! Input: [real(r8) (:)   ]  water table depth (m)
          )
    !-----------------------------------------------------------------------

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
       l = col%landunit(c)

       ! e_folding_length (m) calculations	 
       beta_rad = (rpi/180.) * col%topo_slope(c)
       if (beta_rad <= 0.16) then
           e_folding_length = 20._r8/(1._r8 + 125._r8*beta_rad)
       else if (beta_rad > 0.16) then
           e_folding_length = 1._r8
       end if
 
       if (jwt2(c) < nlevsoi) then
          Transmiss(c) = clayP(c,jwt2(c)) * hk_l(c,jwt2(c)) * (zi(c,jwt2(c))-zwt(c)) * 1.e3_r8
          do j = jwt2(c)+1,nlevsoi
             Transmiss(c) = Transmiss(c) + clayP(c,j) * hk_l(c,j) * dz(c,j) * 1.e3_r8
          end do
          Transmiss(c) = Transmiss(c) + clayP(c,nlevsoi) * hk_l(c,nlevsoi) * e_folding_length * 1.e3_r8

       else if (jwt2(c) .eq. nlevsoi) then
          Transmiss(c) = clayP(c,nlevsoi) * hk_l(c,nlevsoi) * (zi(c,nlevsoi)-zwt(c)) * 1.e3_r8 &
                         + clayP(c,nlevsoi) * hk_l(c,nlevsoi) * e_folding_length * 1.e3_r8

       else if (jwt2(c) .eq. 100) then
          Transmiss(c) = clayP(c,nlevsoi) * hk_l(c,nlevsoi) * e_folding_length * 1.e3_r8 &
                         * exp((zi(c,nlevsoi)-zwt(c))/e_folding_length)
       end if 
    end do
    end associate
  end function TransmissivityFromFan



end module GroundwaterMod