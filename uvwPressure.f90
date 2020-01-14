include "para.f90"

Module Param_one_block
    use Nekdata
    implicit none

    character(50) ::mshfile = 'mesh/pressureSide_block.xyz'    
    character(50) :: infile = 'uvw/pressureSide_block.fun'
    character(50) ::outfile = 'Combined_uvw/pressure.fun'
    integer, parameter :: nex = nnormal, ney = nz, nez = nchord
    
    integer, parameter :: npx = nex*(odr-1)+1, npy = ney*(odr-1)+1, npz = nez*(odr-1)+1
End module
! _______________________________________________________________________________________

Program main
    implicit none
    call combine
end


include "uvwBlock.f90"
