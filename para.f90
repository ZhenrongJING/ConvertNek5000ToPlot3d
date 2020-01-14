Module Nekdata
    implicit none
    character(128):: a
    character(20) :: filename = 'airfoil0.f00001'
    integer, parameter :: odr = 4, npts = odr*odr*odr
    integer, parameter :: nz = 3, nnormal = 20, nchord = 50, nwake = 30, nleading = 18 
    integer, parameter :: netot = 10680
    integer :: ndex(netot), nout, nbeg
End module
