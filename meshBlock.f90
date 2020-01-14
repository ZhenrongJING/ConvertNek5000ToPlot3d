subroutine combine
use Nekdata
use Share
use Param_one_block
implicit none
    
    integer:: ncount, n
    
    integer::  i,  j, k
    integer:: ii, jj, kk
    integer:: ni, nj, nk

    real:: tol
    
    allocate(x(odr, odr, odr, nex, ney, nez))
    allocate(y(odr, odr, odr, nex, ney, nez))
    allocate(z(odr, odr, odr, nex, ney, nez))


    open(101,file = infile,form='unformatted')
    read(101) ncount
    read(101)(ii, jj, kk, i = 1, ncount)
    
    do nk = 1, nez
    do nj = 1, ney
    do ni = 1, nex
    read(101)   (((x(i, j, k, ni, nj, nk),i=1,odr),j=1,odr),k=1,odr),&
                (((y(i, j, k, ni, nj, nk),i=1,odr),j=1,odr),k=1,odr),&
                (((z(i, j, k, ni, nj, nk),i=1,odr),j=1,odr),k=1,odr)
    end do
    end do
    end do
    close(101)
    
    call determin_i
    call determin_j
    call determin_k
    call first_line
    call first_plane
    call first_block
    
    call output

endsubroutine

subroutine output
use Share
use Param_one_block
implicit none
integer:: i, j, k
integer::ni,nj,nk
integer::ei,ej,ek
real, allocatable :: xo(:,:,:)
real, allocatable :: yo(:,:,:)
real, allocatable :: zo(:,:,:)
    
    allocate(xo(npx, npy, npz))
    allocate(yo(npx, npy, npz))
    allocate(zo(npx, npy, npz))
    
    do nk = 1, nez
    do nj = 1, ney
    do ni = 1, nex
        do ek = 1, odr
        do ej = 1, odr
        do ei = 1, odr
            i = (ni-1)*(odr-1) + ei
            j = (nj-1)*(odr-1) + ej
            k = (nk-1)*(odr-1) + ek
            
            xo(i,j,k) = x(ei,ej,ek,ni,nj,nk)
            yo(i,j,k) = y(ei,ej,ek,ni,nj,nk)
            zo(i,j,k) = z(ei,ej,ek,ni,nj,nk)
        end do
        end do
        end do
    end do
    end do
    end do
    
    
    open(101,file = outfile, form = 'unformatted')
    write(101)1
    write(101)npx, npy, npz
    write(101)  (((xo(i, j, k),i=1,npx),j=1,npy),k=1,npz),&
                (((yo(i, j, k),i=1,npx),j=1,npy),k=1,npz),&
                (((zo(i, j, k),i=1,npx),j=1,npy),k=1,npz)
    close(101)
    
!     open(101,file = '../Mesh/test.xyz', form = 'unformatted')
!     write(101)netot
!     write(101)(odr, odr, odr, i = 1, netot)
!     do nk = 1, nez
!     do nj = 1, ney
!     do ni = 1, nex
!     write(101)   (((x(i, j, k, ni, nj, nk),i=1,odr),j=1,odr),k=1,odr),&
!                  (((y(i, j, k, ni, nj, nk),i=1,odr),j=1,odr),k=1,odr),&
!                  (((z(i, j, k, ni, nj, nk),i=1,odr),j=1,odr),k=1,odr)
!     end do
!     end do
!     end do
!     close(101)

deallocate(x, y, z, xo, yo, zo)
endsubroutine

subroutine first_block
use Share
use Param_one_block
implicit none
integer::i, j, k
integer::iflag

do k = 2, nez
    do j = 1, ney
        do i = 1, nex
            call find_connection(i, j, k-1, i, j, k)
            if(ic2(1) == ic2(2) .and. ic2(1) == ic2(3)) call rotate_ki(i,j,k)
            if(jc2(1) == jc2(2) .and. jc2(1) == jc2(3)) call rotate_jk(i,j,k)
            
            call find_connection(i, j, k-1, i, j, k)
            if(ic1(1) == ic1(2) .and. ic2(1) /= ic2(2)) call rotate_ij(i,j,k)
            if(ic1(1) /= ic1(2) .and. ic2(1) == ic2(2)) call rotate_ij(i,j,k)
            
            
            call find_connection(i, j, k-1, i, j, k)
            if(ic1(1) /= ic2(1) .and. ic1(2) /= ic2(2)) call mirror_i(i,j,k)            
            call find_connection(i, j, k-1, i, j, k)
            if(jc1(1) /= jc2(1) .and. jc1(2) /= jc2(2)) call mirror_j(i,j,k)
            call find_connection(i, j, k-1, i, j, k)
            if(kc2(1) == 6) call mirror_k(i,j,k)

            
!             do iflag = 1, 4
!                 print*, ic1(iflag),jc1(iflag),kc1(iflag),' ici', ic2(iflag),jc2(iflag),kc2(iflag)
!             end do
!             stop
        end do
    end do
end do

endsubroutine

subroutine first_plane
use Share
use Param_one_block
implicit none
integer::i, j
integer::iflag

do j = 2, ney
    do i = 1, nex
        call find_connection(i, j-1, 1, i, j, 1)
        if(ic2(1) == ic2(2) .and. ic2(1) == ic2(3)) call rotate_ij(i,j,1)
        if(kc2(1) == kc2(2) .and. kc2(1) == kc2(3)) call rotate_jk(i,j,1)
        
        call find_connection(i, j-1, 1, i, j, 1)
        if(ic1(1) == ic1(2) .and. ic2(1) /= ic2(2)) call rotate_ki(i,j,1)
        if(ic1(1) /= ic1(2) .and. ic2(1) == ic2(2)) call rotate_ki(i,j,1)

        call find_connection(i, j-1, 1, i, j, 1)
        if(kc1(1) /= kc2(1) .and. kc1(2) /= kc2(2)) call mirror_k(i,j,1)

        call find_connection(i, j-1, 1, i, j, 1)
        if(ic1(1) /= ic2(1) .and. ic1(2) /= ic2(2)) call mirror_i(i,j,1)
        
        call find_connection(i, j-1, 1, i, j, 1)
        if(jc2(1) == odr) call mirror_j(i,j,1)

    end do
end do
endsubroutine



subroutine first_line
use Share
use Param_one_block
implicit none
integer::i
integer::iflag

do i = 2, nex
    
    call find_connection(i-1, 1, 1, i, 1, 1)
    if(jc2(1) == jc2(2).and. jc2(1)==jc2(3)) call rotate_ij(i, 1, 1)    
    if(kc2(1) == kc2(2).and. kc2(1)==kc2(3)) call rotate_ki(i, 1, 1)
    
    call find_connection(i-1, 1, 1, i, 1, 1)
    if(jc1(1) == jc1(2) .and. jc2(1) /= jc2(2)) call rotate_jk(i,1,1)
    if(jc1(1) /= jc1(2) .and. jc2(1) == jc2(2)) call rotate_jk(i,1,1)

!    do iflag = 1, 4
!        print*, ic1(iflag),jc1(iflag),kc1(iflag),' ici', ic2(iflag),jc2(iflag),kc2(iflag)
!    end do
!    stop
    call find_connection(i-1, 1, 1, i, 1, 1)
    if (jc1(1).ne.jc2(1))  call mirror_j(i, 1, 1)
    
    call find_connection(i-1, 1, 1, i, 1, 1)
    if (kc1(1).ne.kc2(1))  call mirror_k(i, 1, 1)
    
    call find_connection(i-1, 1, 1, i, 1, 1)
    if(ic2(1) == odr) call mirror_i(i,1,1)

end do
endsubroutine


subroutine determin_i
use Share
implicit none
integer:: iflag
    call find_connection(1,1,1,2,1,1)
    if(jc1(1) == jc1(2).and. jc1(1)==jc1(3)) call rotate_ij(1, 1, 1)    
    if(kc1(1) == kc1(2).and. kc1(1)==kc1(3)) call rotate_ki(1, 1, 1)    
    
    call find_connection(1,1,1,2,1,1)
    if(ic1(1) == ic1(2).and. ic1(1)==ic1(3) .and. ic1(1) == 1) call mirror_i(1, 1, 1)
    
endsubroutine


subroutine determin_j
use Share
implicit none
integer:: iflag
    call find_connection(1,1,1,1,2,1)
    if(kc1(1) == kc1(2).and. kc1(1)==kc1(3)) call rotate_jk(1, 1, 1)    
    
    call find_connection(1,1,1,1,2,1)
    if(jc1(1) == jc1(2).and. jc1(1)==jc1(3) .and. jc1(1) == 1) call mirror_j(1, 1, 1)
    
!     do iflag = 1, 4
!         print*, ic1(iflag),jc1(iflag),kc1(iflag),' ici', ic2(iflag),jc2(iflag),kc2(iflag)
!     end do
!     
endsubroutine

subroutine determin_k
use Share
implicit none
integer:: iflag
    
    call find_connection(1,1,1,1,1,2)
    if(kc1(1) == kc1(2).and. kc1(1)==kc1(3) .and. kc1(1) == 1) call mirror_k(1, 1, 1)
!     do iflag = 1, 4
!         print*, ic1(iflag),jc1(iflag),kc1(iflag),' ici', ic2(iflag),jc2(iflag),kc2(iflag)
!     end do
endsubroutine

subroutine find_connection(i_fst, j_fst, k_fst, i_scd, j_scd, k_scd)
use Share
implicit none
integer:: i, j, k
integer::ii,jj,kk
integer:: i_fst, j_fst, k_fst, i_scd, j_scd, k_scd
integer:: iflag
real:: tol
    iflag = 0
    do i = 1, odr, odr - 1
    do j = 1, odr, odr - 1
    do k = 1, odr, odr - 1
        do ii = 1, odr, odr - 1
        do jj = 1, odr, odr - 1
        do kk = 1, odr, odr - 1
            tol =       abs(x(i, j, k, i_fst, j_fst, k_fst) - x(ii, jj, kk, i_scd, j_scd, k_scd)) 
            tol = tol + abs(y(i, j, k, i_fst, j_fst, k_fst) - y(ii, jj, kk, i_scd, j_scd, k_scd))
            tol = tol + abs(z(i, j, k, i_fst, j_fst, k_fst) - z(ii, jj, kk, i_scd, j_scd, k_scd))
            if(tol.le.0.000001) then
                iflag = iflag + 1
                ic1(iflag) = i
                jc1(iflag) = j
                kc1(iflag) = k
                ic2(iflag) = ii
                jc2(iflag) = jj
                kc2(iflag) = kk
            end if
        end do
        end do
        end do
    end do
    end do
    end do
endsubroutine

subroutine mirror_i(ni,nj,nk)
use Share
implicit none
    integer::  i,  j, k
    integer:: ni, nj, nk
    
    real:: x_s(odr, odr, odr)
    real:: y_s(odr, odr, odr)
    real:: z_s(odr, odr, odr)
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x_s(i,j,k) = x(i, j, k, ni, nj, nk)
        y_s(i,j,k) = y(i, j, k, ni, nj, nk)
        z_s(i,j,k) = z(i, j, k, ni, nj, nk)
    end do
    end do
    end do
    
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x(i, j, k, ni, nj, nk) = x_s(odr+1-i, j, k)
        y(i, j, k, ni, nj, nk) = y_s(odr+1-i, j, k)
        z(i, j, k, ni, nj, nk) = z_s(odr+1-i, j, k)
    end do
    end do
    end do

endsubroutine

subroutine mirror_j(ni,nj,nk)
use Share
implicit none
    integer::  i,  j, k
    integer:: ni, nj, nk
    
    real:: x_s(odr, odr, odr)
    real:: y_s(odr, odr, odr)
    real:: z_s(odr, odr, odr)
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x_s(i,j,k) = x(i, j, k, ni, nj, nk)
        y_s(i,j,k) = y(i, j, k, ni, nj, nk)
        z_s(i,j,k) = z(i, j, k, ni, nj, nk)
    end do
    end do
    end do
    
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x(i, j, k, ni, nj, nk) = x_s(i, odr-j+1, k)
        y(i, j, k, ni, nj, nk) = y_s(i, odr-j+1, k)
        z(i, j, k, ni, nj, nk) = z_s(i, odr-j+1, k)
    end do
    end do
    end do

endsubroutine

subroutine mirror_k(ni,nj,nk)
use Share
implicit none
    integer::  i,  j, k
    integer:: ni, nj, nk
    
    real:: x_s(odr, odr, odr)
    real:: y_s(odr, odr, odr)
    real:: z_s(odr, odr, odr)
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x_s(i,j,k) = x(i, j, k, ni, nj, nk)
        y_s(i,j,k) = y(i, j, k, ni, nj, nk)
        z_s(i,j,k) = z(i, j, k, ni, nj, nk)
    end do
    end do
    end do
    
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x(i, j, k, ni, nj, nk) = x_s(i, j, odr-k+1)
        y(i, j, k, ni, nj, nk) = y_s(i, j, odr-k+1)
        z(i, j, k, ni, nj, nk) = z_s(i, j, odr-k+1)
    end do
    end do
    end do

endsubroutine

subroutine rotate_ij(ni,nj,nk)
use Share
implicit none
    integer::  i,  j, k
    integer:: ni, nj, nk
    
    real:: x_s(odr, odr, odr)
    real:: y_s(odr, odr, odr)
    real:: z_s(odr, odr, odr)
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x_s(i,j,k) = x(i, j, k, ni, nj, nk)
        y_s(i,j,k) = y(i, j, k, ni, nj, nk)
        z_s(i,j,k) = z(i, j, k, ni, nj, nk)
    end do
    end do
    end do

    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x(i, j, k, ni, nj, nk) = x_s(j, i, k)
        y(i, j, k, ni, nj, nk) = y_s(j, i, k)
        z(i, j, k, ni, nj, nk) = z_s(j, i, k)
    end do
    end do
    end do

endsubroutine


subroutine rotate_jk(ni,nj,nk)
use Share
implicit none
    integer::  i,  j, k
    integer:: ni, nj, nk
    
    real:: x_s(odr, odr, odr)
    real:: y_s(odr, odr, odr)
    real:: z_s(odr, odr, odr)
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x_s(i,j,k) = x(i, j, k, ni, nj, nk)
        y_s(i,j,k) = y(i, j, k, ni, nj, nk)
        z_s(i,j,k) = z(i, j, k, ni, nj, nk)
    end do
    end do
    end do

    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x(i, j, k, ni, nj, nk) = x_s(i,k,j)
        y(i, j, k, ni, nj, nk) = y_s(i,k,j)
        z(i, j, k, ni, nj, nk) = z_s(i,k,j)
    end do
    end do
    end do

endsubroutine


subroutine rotate_ki(ni,nj,nk)
use Share
implicit none
    integer::  i,  j, k
    integer:: ni, nj, nk
    
    real:: x_s(odr, odr, odr)
    real:: y_s(odr, odr, odr)
    real:: z_s(odr, odr, odr)
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x_s(i,j,k) = x(i, j, k, ni, nj, nk)
        y_s(i,j,k) = y(i, j, k, ni, nj, nk)
        z_s(i,j,k) = z(i, j, k, ni, nj, nk)
    end do
    end do
    end do
    
    do k = 1, odr
    do j = 1, odr
    do i = 1, odr    
        x(i, j, k, ni, nj, nk) = x_s(k, j, i)
        y(i, j, k, ni, nj, nk) = y_s(k, j, i)
        z(i, j, k, ni, nj, nk) = z_s(k, j, i)
    end do
    end do
    end do
endsubroutine
