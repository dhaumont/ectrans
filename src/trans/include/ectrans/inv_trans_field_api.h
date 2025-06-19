interface

subroutine inv_trans_field_api(ydfspvor,ydfspdiv,ydfspscalar, &
                             & ydfu, ydfv, ydfvor,ydfdiv,ydfscalar, &
                             & ydfudm, ydfvdm, ydfscalardm, ydfscalardl,& 
                             & kspec, kproma, kgpblks, kgptot, kflevg, kflevl,&
                             & ldacc, ldverbose, &
                             & fspgl_proc)


type(field_basic_ptr),intent(in), optional  :: ydfspvor(:), ydfspdiv(:)        ! spectral vector fields : vorticity and divergence fields (in)
type(field_basic_ptr),intent(in), optional  :: ydfspscalar(:)                  ! spectral scalar fields (in)

type(field_basic_ptr),intent(in), optional  :: ydfu(:),ydfv(:)                 ! grid vector fields     (out)
type(field_basic_ptr),intent(in), optional  :: ydfvor(:),ydfdiv(:)             ! grid vector fields :vorticity and divergence     (out)
type(field_basic_ptr),intent(in), optional  :: ydfscalar(:)                    ! grid scalar fields     (out)

type(field_basic_ptr),intent(in), optional  :: ydfudm(:),ydfvdm(:)             ! grid vector fields derivatives ew (out)
type(field_basic_ptr),intent(in), optional  :: ydfscalardm(:), ydfscalardl(:)  ! grid scalar fields derivatives ew and ns (out)

integer(kind=jpim),   intent(in)            :: kspec
integer(kind=jpim),   intent(in)            :: kproma
integer(kind=jpim),   intent(in)            :: kgpblks
integer(kind=jpim),   intent(in)            :: kgptot
integer(kind=jpim),   intent(in)            :: kflevg
integer(kind=jpim),   intent(in)            :: kflevl
logical,              intent(in), optional  :: ldacc
logical,              intent(in), optional  :: ldverbose
procedure (fspgl_intf),           optional  :: fspgl_proc

#include "inv_trans.h"
#include "abor1.intfb.h"

end subroutine inv_trans_field_api

end interface
