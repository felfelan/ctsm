module mkUSGSGWratioMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkUSGSGWratioMod
!
! !DESCRIPTION:
! make USGS GW contribution ratio to the total irrigation
!
! !REVISION HISTORY:
! Author: Farshid Felfelani
!
!-----------------------------------------------------------------------
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public mkUSGSGWratio            ! regrid USGS GW data
!
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkUSGSGWratio
!
! !INTERFACE:
subroutine mkUSGSGWratio(ldomain, mapfname, datfname, ndiag, USGS_o)
!
! !DESCRIPTION:
! make USGS GW from input USGS GW data
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkncdio
  use mkdiagnosticsMod, only : output_diagnostics_continuous
  use mkchecksMod, only : min_bad
!
! !ARGUMENTS:
  
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  real(r8)          , intent(out):: USGS_o(:)  ! output grid: USGS GW
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Farshid Felfelani
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain            ! local domain
  real(r8), allocatable :: data_i(:)          ! data on input grid
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status

  real(r8), parameter :: min_valid = 0._r8    ! minimum valid value

  character(len=32) :: subname = 'mkUSGSGWratio'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make USGS GW.....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read domain and mapping information, check for consistency
  ! -----------------------------------------------------------------

  call domain_read(tdomain,datfname)
  
  call gridmap_mapread(tgridmap, mapfname )
  call gridmap_check( tgridmap, subname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! -----------------------------------------------------------------
  ! Open input file, allocate memory for input data
  ! -----------------------------------------------------------------

  write(6,*)'Open USGS GW file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(data_i(tdomain%ns), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid USGS GW
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'USGS_mean', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, USGS_o, nodata=0._r8)

  ! Check validity of output data
  if (min_bad(USGS_o, min_valid, 'USGS_mean')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, USGS_o, tgridmap, "USGS_mean", "USGS contribution of GW to irrigation(fraction)", ndiag)

  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (data_i)

  write (6,*) 'Successfully made USGS GW'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkUSGSGWratio

end module mkUSGSGWratioMod
