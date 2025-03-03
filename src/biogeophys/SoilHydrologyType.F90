Module SoilHydrologyType

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use decompMod             , only : bounds_type
  use clm_varpar            , only : nlevgrnd, nlayer, nlayert, nlevsoi
  use clm_time_manager      , only : get_curr_date, get_start_date
  use clm_varcon            , only : spval
  use clm_varctl            , only : iulog
  use LandunitType          , only : lun                
  use ColumnType            , only : col                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  type, public :: soilhydrology_type

     integer :: h2osfcflag              ! true => surface water is active (namelist)        
     integer :: origflag                ! used to control soil hydrology properties (namelist)

     real(r8), pointer :: num_substeps_col   (:)    ! col adaptive timestep counter     
     ! NON-VIC
     real(r8), pointer :: frost_table_col   (:)     ! col frost table depth                    
     real(r8), pointer :: zwt_col           (:)     ! col water table depth
     real(r8), pointer :: rechclim_col           (:)     ! col rechclim
     real(r8), pointer :: zwts_col          (:)     ! col water table depth, the shallower of the two water depths
     real(r8), pointer :: zwt_perched_col   (:)     ! col perched water table depth
     real(r8), pointer :: wa_col            (:)     ! col water in the unconfined aquifer (mm)
     real(r8), pointer :: Qgw_lateral_col   (:)     ! col Groundwater lateral flow (mm/s)
     real(r8), pointer :: AqTransmiss_col   (:)     ! col Aquifer Transmissivity (mm)
     real(r8), pointer :: Pump_wa_col       (:)     ! col pumped water (mm/s)
     real(r8), pointer :: QlatField_northing_grc (:)! grc Groundwater lateral flow towards north (+)(mm)
     real(r8), pointer :: QlatField_easting_grc  (:)! grc Groundwater lateral flow towards east (+)(mm)
     real(r8), pointer :: Qgw_lateral_grc   (:)     ! grc Groundwater lateral budget (mm)
     ! real(r8), pointer :: Nneighbors             (:)
     real(r8), pointer :: qcharge_col       (:)     ! col aquifer recharge rate (mm/s) 
     real(r8), pointer :: qcharge_org_col   (:)     ! col original aquifer recharge rate (mm/s) 
     real(r8), pointer :: fracice_col       (:,:)   ! col fractional impermeability (-)
     real(r8), pointer :: icefrac_col       (:,:)   ! col fraction of ice       
     real(r8), pointer :: fcov_col          (:)     ! col fractional impermeable area
     real(r8), pointer :: fsat_col          (:)     ! col fractional area with water table at surface
     real(r8), pointer :: h2osfc_thresh_col (:)     ! col level at which h2osfc "percolates"   (time constant)

     ! VIC 
     real(r8), pointer :: hkdepth_col       (:)     ! col VIC decay factor (m) (time constant)                    
     real(r8), pointer :: b_infil_col       (:)     ! col VIC b infiltration parameter (time constant)                    
     real(r8), pointer :: ds_col            (:)     ! col VIC fracton of Dsmax where non-linear baseflow begins (time constant)                    
     real(r8), pointer :: dsmax_col         (:)     ! col VIC max. velocity of baseflow (mm/day) (time constant)
     real(r8), pointer :: Wsvic_col         (:)     ! col VIC fraction of maximum soil moisutre where non-liear base flow occurs (time constant)
     real(r8), pointer :: porosity_col      (:,:)   ! col VIC porosity (1-bulk_density/soil_density)
     real(r8), pointer :: vic_clm_fract_col (:,:,:) ! col VIC fraction of VIC layers in CLM layers 
     real(r8), pointer :: depth_col         (:,:)   ! col VIC layer depth of upper layer  
     real(r8), pointer :: c_param_col       (:)     ! col VIC baseflow exponent (Qb) 
     real(r8), pointer :: expt_col          (:,:)   ! col VIC pore-size distribution related paramter(Q12) 
     real(r8), pointer :: ksat_col          (:,:)   ! col VIC Saturated hydrologic conductivity 
     real(r8), pointer :: phi_s_col         (:,:)   ! col VIC soil moisture dissusion parameter 
     real(r8), pointer :: moist_col         (:,:)   ! col VIC soil moisture (kg/m2) for VIC soil layers 
     real(r8), pointer :: moist_vol_col     (:,:)   ! col VIC volumetric soil moisture for VIC soil layers 
     real(r8), pointer :: max_moist_col     (:,:)   ! col VIC max layer moist + ice (mm) 
     real(r8), pointer :: max_infil_col     (:)     ! col VIC maximum infiltration rate calculated in VIC
     real(r8), pointer :: i_0_col           (:)     ! col VIC average saturation in top soil layers 
     real(r8), pointer :: ice_col           (:,:)   ! col VIC soil ice (kg/m2) for VIC soil layers

   contains

     ! Public routines
     procedure, public  :: Init
     procedure, public  :: Restart

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private :: ReadNL

  end type soilhydrology_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)

    class(soilhydrology_type) :: this
    type(bounds_type), intent(in)    :: bounds  
    character(len=*), intent(in) :: NLFilename

    call this%ReadNL(NLFilename)
    call this%InitAllocate(bounds) 
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%num_substeps_col  (begc:endc))                 ; this%num_substeps_col  (:)     = nan
    allocate(this%frost_table_col   (begc:endc))                 ; this%frost_table_col   (:)     = nan
    allocate(this%zwt_col           (begc:endc))                 ; this%zwt_col           (:)     = nan
    allocate(this%rechclim_col      (begc:endc))                 ; this%rechclim_col      (:)     = nan
    allocate(this%zwt_perched_col   (begc:endc))                 ; this%zwt_perched_col   (:)     = nan
    allocate(this%zwts_col          (begc:endc))                 ; this%zwts_col          (:)     = nan

    allocate(this%wa_col            (begc:endc))                 ; this%wa_col            (:)     = nan
    allocate(this%Qgw_lateral_col   (begc:endc))                 ; this%Qgw_lateral_col   (:)     = nan
    allocate(this%AqTransmiss_col   (begc:endc))                 ; this%AqTransmiss_col   (:)     = nan
    allocate(this%Pump_wa_col       (begc:endc))                 ; this%Pump_wa_col       (:)     = nan
    allocate(this%QlatField_northing_grc (begg:endg))            ; this%QlatField_northing_grc(:) = nan
    allocate(this%QlatField_easting_grc  (begg:endg))            ; this%QlatField_easting_grc (:) = nan
    allocate(this%Qgw_lateral_grc   (begg:endg))                 ; this%Qgw_lateral_grc   (:)     = nan
    ! allocate(this%Nneighbors (begg:endg))                        ; this%Nneighbors(:) = nan
	
    allocate(this%qcharge_col       (begc:endc))                 ; this%qcharge_col       (:)     = nan
    allocate(this%qcharge_org_col   (begc:endc))                 ; this%qcharge_org_col   (:)     = nan
    allocate(this%fracice_col       (begc:endc,nlevgrnd))        ; this%fracice_col       (:,:)   = nan
    allocate(this%icefrac_col       (begc:endc,nlevgrnd))        ; this%icefrac_col       (:,:)   = nan
    allocate(this%fcov_col          (begc:endc))                 ; this%fcov_col          (:)     = nan   
    allocate(this%fsat_col          (begc:endc))                 ; this%fsat_col          (:)     = nan
    allocate(this%h2osfc_thresh_col (begc:endc))                 ; this%h2osfc_thresh_col (:)     = nan

    allocate(this%hkdepth_col       (begc:endc))                 ; this%hkdepth_col       (:)     = nan
    allocate(this%b_infil_col       (begc:endc))                 ; this%b_infil_col       (:)     = nan
    allocate(this%ds_col            (begc:endc))                 ; this%ds_col            (:)     = nan
    allocate(this%dsmax_col         (begc:endc))                 ; this%dsmax_col         (:)     = nan
    allocate(this%Wsvic_col         (begc:endc))                 ; this%Wsvic_col         (:)     = nan
    allocate(this%depth_col         (begc:endc,nlayert))         ; this%depth_col         (:,:)   = nan
    allocate(this%porosity_col      (begc:endc,nlayer))          ; this%porosity_col      (:,:)   = nan
    allocate(this%vic_clm_fract_col (begc:endc,nlayer, nlevsoi)) ; this%vic_clm_fract_col (:,:,:) = nan
    allocate(this%c_param_col       (begc:endc))                 ; this%c_param_col       (:)     = nan
    allocate(this%expt_col          (begc:endc,nlayer))          ; this%expt_col          (:,:)   = nan
    allocate(this%ksat_col          (begc:endc,nlayer))          ; this%ksat_col          (:,:)   = nan
    allocate(this%phi_s_col         (begc:endc,nlayer))          ; this%phi_s_col         (:,:)   = nan
    allocate(this%moist_col         (begc:endc,nlayert))         ; this%moist_col         (:,:)   = nan
    allocate(this%moist_vol_col     (begc:endc,nlayert))         ; this%moist_vol_col     (:,:)   = nan
    allocate(this%max_moist_col     (begc:endc,nlayer))          ; this%max_moist_col     (:,:)   = nan
    allocate(this%max_infil_col     (begc:endc))                 ; this%max_infil_col     (:)     = nan
    allocate(this%i_0_col           (begc:endc))                 ; this%i_0_col           (:)     = nan
    allocate(this%ice_col           (begc:endc,nlayert))         ; this%ice_col           (:,:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    use decompMod       , only : ldecomp
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begc, endc
    integer           :: begg, endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%wa_col(begc:endc) = spval
    call hist_addfld1d (fname='WA',  units='mm',  &
         avgflag='A', long_name='water in the unconfined aquifer (vegetated landunits only)', &
         ptr_col=this%wa_col, l2g_scale_type='veg')

    this%Qgw_lateral_col(begc:endc) = spval
    call hist_addfld1d (fname='Qgw_lateral',  units='mm/s',  &
         avgflag='A', long_name='Groundwater lateral water in the unconfined aquifer (vegetated landunits only)', &
         ptr_col=this%Qgw_lateral_col, l2g_scale_type='veg')

    this%Qgw_lateral_grc(begg:endg) = spval
    call hist_addfld1d (fname='Qgw_lateral_grc',  units='mm',  &
         avgflag='SUM', long_name='Groundwater lateral water in the unconfined aquifer (gridcell level)', &
         ptr_lnd=this%Qgw_lateral_grc, l2g_scale_type='veg')

    this%AqTransmiss_col(begc:endc) = spval
    call hist_addfld1d (fname='Aq_Transmissivity',  units='mm2/s',  &
         avgflag='A', long_name='Transmissivity of the unconfined aquifer (vegetated landunits only)', &
         ptr_col=this%AqTransmiss_col, l2g_scale_type='veg')

    this%Pump_wa_col(begc:endc) = spval
    call hist_addfld1d (fname='Pumped_Wa',  units='mm/s',  &
         avgflag='A', long_name='Pumped Water from the unconfined aquifer (vegetated landunits only)', &
         ptr_col=this%Pump_wa_col, l2g_scale_type='veg')

    this%QlatField_northing_grc(begg:endg) = spval
    call hist_addfld1d (fname='QlatField_northing_grc',  units='mm',  &
         avgflag='A', long_name='Northward groundwater lateral flow', &
         ptr_lnd=this%QlatField_northing_grc, l2g_scale_type='veg')

    this%QlatField_easting_grc(begg:endg) = spval
    call hist_addfld1d (fname='QlatField_easting_grc',  units='mm',  &
         avgflag='A', long_name='Eastward groundwater lateral flow', &
         ptr_lnd=this%QlatField_easting_grc, l2g_scale_type='veg')

    ! this%Nneighbors(begg:endg) = spval
    ! call hist_addfld1d (fname='Nneighbors',  units='unitless',  &
         ! avgflag='A', long_name='Number of Neighbors', &
         ! ptr_lnd=ldecomp%gneighbors, l2g_scale_type='veg')

    this%qcharge_col(begc:endc) = spval
    call hist_addfld1d (fname='QCHARGE',  units='mm/s',  &
         avgflag='A', long_name='aquifer recharge rate (vegetated landunits only)', &
         ptr_col=this%qcharge_col, l2g_scale_type='veg')

    this%qcharge_org_col(begc:endc) = spval
    call hist_addfld1d (fname='QCHARGE_ORG',  units='mm/s',  &
         avgflag='A', long_name='original aquifer recharge rate (vegetated landunits only)', &
         ptr_col=this%qcharge_org_col, l2g_scale_type='veg')

    this%fcov_col(begc:endc) = spval
    call hist_addfld1d (fname='FCOV',  units='unitless',  &
         avgflag='A', long_name='fractional impermeable area', &
         ptr_col=this%fcov_col, l2g_scale_type='veg')

    this%fsat_col(begc:endc) = spval
    call hist_addfld1d (fname='FSAT',  units='unitless',  &
         avgflag='A', long_name='fractional area with water table at surface', &
         ptr_col=this%fsat_col, l2g_scale_type='veg')

    this%num_substeps_col(begc:endc) = spval
    call hist_addfld1d (fname='NSUBSTEPS',  units='unitless',  &
         avgflag='A', long_name='number of adaptive timesteps in CLM timestep', &
         ptr_col=this%num_substeps_col, l2g_scale_type='veg', &
         default='inactive')

    this%frost_table_col(begc:endc) = spval
    call hist_addfld1d (fname='FROST_TABLE',  units='m',  &
         avgflag='A', long_name='frost table depth (vegetated landunits only)', &
         ptr_col=this%frost_table_col, l2g_scale_type='veg', default='inactive')

    this%zwt_col(begc:endc) = spval
    call hist_addfld1d (fname='ZWT',  units='m',  &
         avgflag='A', long_name='water table depth (vegetated landunits only)', &
         ptr_col=this%zwt_col, l2g_scale_type='veg')

    this%zwt_perched_col(begc:endc) = spval
    call hist_addfld1d (fname='ZWT_PERCH',  units='m',  &
         avgflag='A', long_name='perched water table depth (vegetated landunits only)', &
         ptr_col=this%zwt_perched_col, l2g_scale_type='veg')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type) , intent(in)    :: bounds
    ! !LOCAL VARIABLES:
    integer :: c ! indices

    !-----------------------------------------------------------------------

    ! Nothing for now

    ! needs to be initialized to spval to avoid problems when 
    ! averaging for the accum field
    do c = bounds%begc, bounds%endc
       this%num_substeps_col(c) = spval
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_io, ncd_double, ncd_pio_openfile
    use clm_varctl , only : gwFanInit, fsurdat
	use clm_varcon , only : grlnd
	use fileutils  , only : getfil
	use abortutils , only : endrun
	use spmdMod        , only : masterproc
	use column_varcon  , only : icol_road_perv
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
	integer :: p,c,j,l,g,lev,nlevs ! indices
    logical :: readvar      ! determine if variable is on initial file
    type(file_desc_t)  :: ncid_srf 
	character(len=256) :: locfn 
    real(r8) ,pointer  :: wtd_Fan      (:)   ! read in - WTD	
    real(r8) ,pointer  :: rechclim_Fan (:)   ! read in - Climatologic Recharge	

    integer :: year, rsyr      ! year (0, ...) for nstep
    integer :: month, rsmon    ! month (1, ..., 12) for nstep
	integer :: day, rsday
    integer :: secs, tod       ! seconds into current date for nstep

    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='FROST_TABLE', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='frost table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%frost_table_col)
    if (flag == 'read' .and. .not. readvar) then
       this%frost_table_col(bounds%begc:bounds%endc) = col%zi(bounds%begc:bounds%endc,nlevsoi)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='WA', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='water in the unconfined aquifer', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%wa_col)

    call restartvar(ncid=ncid, flag=flag, varname='Qgw_lateral', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='Groundwater lateral water in the unconfined aquifer', units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%Qgw_lateral_col)

    call restartvar(ncid=ncid, flag=flag, varname='Aq_Transmissivity', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='Transmissivity of the unconfined aquifer', units='mm2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%AqTransmiss_col)

    call restartvar(ncid=ncid, flag=flag, varname='Pumped_Wa', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='Pumped water from the unconfined aquifer', units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%Pump_wa_col)

    call restartvar(ncid=ncid, flag=flag, varname='QlatField_northing_grc', xtype=ncd_double,  & 
         dim1name='gridcell', &
         long_name='Northward groundwater lateral flow', units='mm', &
         interpinic_flag='skip', readvar=readvar, data=this%QlatField_northing_grc)		 
		 
    call restartvar(ncid=ncid, flag=flag, varname='QlatField_easting_grc', xtype=ncd_double,  & 
         dim1name='gridcell', &
         long_name='Eastward groundwater lateral flow', units='mm', &
         interpinic_flag='skip', readvar=readvar, data=this%QlatField_easting_grc)		 

    call get_curr_date (year, month, day, secs)
	call get_start_date(rsyr, rsmon, rsday, tod)

	if (flag == 'read' .and. gwFanInit == .true. .and. year == rsyr) then
        if (masterproc) write(iulog,*) '                                        '
        if (masterproc) write(iulog,*) '****************************************'
        if (masterproc) write(iulog,*) '****************************************'
        if (masterproc) write(iulog,*) '****************************************'
        if (masterproc) write(iulog,*) '                                        '
        if (masterproc) write(iulog,*) 'GW is initialized by Fan et al. data: Flag ', flag
		if (masterproc) write(iulog,*) 'year and rsyr: ', year, rsyr
        if (masterproc) write(iulog,*) '                                        '

         call restartvar(ncid=ncid, flag=flag, varname='ZWT', xtype=ncd_double,  & 
              dim1name='column', &
              long_name='water table depth', units='m', &
              interpinic_flag='interp', readvar=readvar, data=this%zwt_col)

        this%rechclim_col(bounds%begc:bounds%endc) = 0._r8
		this%zwt_col(bounds%begc:bounds%endc) = 0._r8
	
        allocate(wtd_Fan(bounds%begg:bounds%endg))
        allocate(rechclim_Fan(bounds%begg:bounds%endg))

        call getfil (fsurdat, locfn, 0)
        call ncd_pio_openfile (ncid_srf, locfn, 0)
	
        call ncd_io(ncid=ncid_srf, varname='EQZWT', flag='read', data=wtd_Fan, dim1name=grlnd, readvar=readvar)

		if (.not. readvar) then
           call endrun(msg=' ERROR: EQZWT NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
        end if

        call ncd_io(ncid=ncid_srf, varname='RECHCLIM', flag='read', data=rechclim_Fan, dim1name=grlnd, readvar=readvar)

		if (.not. readvar) then
           call endrun(msg=' ERROR: RECHCLIM NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
        end if

        do c = bounds%begc,bounds%endc
            g = col%gridcell(c)
            l = col%landunit(c)
            if (.not. lun%lakpoi(l)) then  !not lake
                if (lun%urbpoi(l)) then
                    if (col%itype(c) == icol_road_perv) then
                        ! Note that the following hard-coded constants (on the next two lines)
                        ! seem implicitly related to aquifer_water_baseline
                        this%zwt_col(c) = wtd_Fan(g)
					    this%rechclim_col(c) = rechclim_Fan(g)
                    else
                        this%zwt_col(c) = spval
					    this%rechclim_col(c) = spval
                    end if

                else
                    ! Note that the following hard-coded constants (on the next two lines) seem
                    ! implicitly related to aquifer_water_baseline
                    this%zwt_col(c) = wtd_Fan(g)
				    this%rechclim_col(c) = rechclim_Fan(g)
                end if
            end if
        end do

        deallocate(wtd_Fan)
        deallocate(rechclim_Fan)

    else
        if (masterproc) write(iulog,*) '                                        '
        if (masterproc) write(iulog,*) '****************************************'
        if (masterproc) write(iulog,*) '****************************************'
        if (masterproc) write(iulog,*) '****************************************'
        if (masterproc) write(iulog,*) '                                        '
        if (masterproc) write(iulog,*) 'GW is defined/read/written to the restart file: Flag ', flag
		if (masterproc) write(iulog,*) 'year and rsyr: ', year, rsyr
        if (masterproc) write(iulog,*) '                                        '

         call restartvar(ncid=ncid, flag=flag, varname='ZWT', xtype=ncd_double,  & 
              dim1name='column', &
              long_name='water table depth', units='m', &
              interpinic_flag='interp', readvar=readvar, data=this%zwt_col)

    end if

    call restartvar(ncid=ncid, flag=flag, varname='ZWT_PERCH', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='perched water table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%zwt_perched_col)
    if (flag == 'read' .and. .not. readvar) then
       this%zwt_perched_col(bounds%begc:bounds%endc) = col%zi(bounds%begc:bounds%endc,nlevsoi)
    end if

  end subroutine Restart

   !-----------------------------------------------------------------------
   subroutine ReadNL( this, NLFilename )
     !
     ! !DESCRIPTION:
     ! Read namelist for SoilHydrology
     !
     ! !USES:
     use shr_mpi_mod    , only : shr_mpi_bcast
     use shr_log_mod    , only : errMsg => shr_log_errMsg
     use spmdMod        , only : masterproc, mpicom
     use fileutils      , only : getavu, relavu, opnfil
     use clm_nlUtilsMod , only : find_nlgroup_name
     use clm_varctl     , only : iulog 
     use abortutils     , only : endrun
     !
     ! !ARGUMENTS:
     class(soilhydrology_type) :: this
     character(len=*), intent(IN) :: NLFilename ! Namelist filename
     !
     ! !LOCAL VARIABLES:
     integer :: ierr                 ! error code
     integer :: unitn                ! unit for namelist file
     integer :: origflag=0            !use to control soil hydraulic properties
     integer :: h2osfcflag=1          !If surface water is active or not
     character(len=32) :: subname = 'SoilHydrology_readnl'  ! subroutine name
     !-----------------------------------------------------------------------

     namelist / clm_soilhydrology_inparm / h2osfcflag, origflag

     ! preset values

     origflag = 0          
     h2osfcflag = 1        

     if ( masterproc )then

        unitn = getavu()
        write(iulog,*) 'Read in clm_soilhydrology_inparm  namelist'
        call opnfil (NLFilename, unitn, 'F')
        call find_nlgroup_name(unitn, 'clm_soilhydrology_inparm', status=ierr)
        if (ierr == 0) then
           read(unitn, clm_soilhydrology_inparm, iostat=ierr)
           if (ierr /= 0) then
              call endrun(msg="ERROR reading clm_soilhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
           end if
        else
           call endrun(msg="ERROR finding clm_soilhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
        end if
        call relavu( unitn )

     end if

     call shr_mpi_bcast(h2osfcflag, mpicom)
     call shr_mpi_bcast(origflag,   mpicom)

     this%h2osfcflag = h2osfcflag
     this%origflag   = origflag

   end subroutine ReadNL

end Module SoilHydrologyType
