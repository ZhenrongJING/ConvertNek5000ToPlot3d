include "para.f90"

Module Share
use Nekdata
    implicit none
    integer::  ic1(4), jc1(4), kc1(4)
    integer::  ic2(4), jc2(4), kc2(4)
    real, allocatable :: x(:,:,:,:,:,:)
    real, allocatable :: y(:,:,:,:,:,:)
    real, allocatable :: z(:,:,:,:,:,:)
End module

Module Param_one_block
    use Nekdata
    implicit none
    
    character(50) :: infile = 'mesh/pressureSide_block.xyz'
    character(50) ::outfile = 'Combined_mesh/Combined_PressureSide_Block.xyz'
    integer, parameter :: nex = nnormal, ney = nz, nez = nchord
    integer, parameter :: npx = nex*(odr-1)+1, npy = ney*(odr-1)+1, npz = nez*(odr-1)+1
End module
! _______________________________________________________________________________________

Program main
    implicit none
    call combine
end

include "meshBlock.f90"
