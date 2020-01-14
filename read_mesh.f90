include "para.f90"

Program main
    use Nekdata
    implicit none
    integer:: i
    call readnek
end

subroutine readnek
use Nekdata
implicit none
    integer:: i, j, k, n, ncount, imax = 0
    integer:: ii, jj, kk
    integer:: ni, nj, nk
    integer:: nex, ney, nez
    integer:: nstart

    real:: time
    
    real, allocatable :: x(:,:,:,:)
    real, allocatable :: y(:,:,:,:)
    real, allocatable :: z(:,:,:,:)
    
    allocate(x(odr,odr,odr,netot))
    allocate(y(odr,odr,odr,netot))
    allocate(z(odr,odr,odr,netot))
    
    open(101,file=filename,form='unformatted')
    read(101)a
    write(*,*)a
    close(101)
    
    ncount = 34
    open(101,file = filename,form='unformatted',recl =4, access='direct', convert='little_endian')
    do i = 1, netot
        ncount = ncount + 1
        read(101,rec = ncount) ndex(i)
        if(imax.le.ndex(i)) imax = ndex(i)
    enddo
    close(101)
    
    if(imax.ne.netot) then
        print*, imax, netot
        stop
    end if
    
    open(101,file = filename,form='unformatted',recl =4, access='direct', convert='little_endian')
    do n = 1, netot
        if(mod(n, 10000).eq.0) print*, n
        do i = 1, npts
            ncount = ncount + 1
            read(101,rec = ncount) x(i,1,1,ndex(n))
        enddo
        do i = 1, npts
            ncount = ncount + 1
            read(101,rec = ncount) y(i,1,1,ndex(n))
        enddo
        do i = 1, npts
            ncount = ncount + 1
            read(101,rec = ncount) z(i,1,1,ndex(n))
        enddo
    enddo
    close(101)

! output! output! output! output! output! output! output! output! output
    nout = nz*nleading*nnormal
    open(103,file = 'mesh/leadingEdge_block.xyz', form = 'unformatted')
    write(103)nout
    write(103)(odr, odr, odr,i=1, nout)
    do n = nstart + 1, nstart + nout
        write(103)  (((x(i,j,k,n),i=1,odr),j=1,odr),k=1,odr),&
                    (((y(i,j,k,n),i=1,odr),j=1,odr),k=1,odr),&
                    (((z(i,j,k,n),i=1,odr),j=1,odr),k=1,odr)
    enddo
    close(103)
    nstart = nstart + nout

    nout = nz*nchord*nnormal
    open(103,file = 'mesh/suctionSide_block.xyz', form = 'unformatted')
    write(103)nout
    write(103)(odr, odr, odr,i=1, nout)
    do n = nstart + 1, nstart + nout
        write(103)  (((x(i,j,k,n),i=1,odr),j=1,odr),k=1,odr),&
                    (((y(i,j,k,n),i=1,odr),j=1,odr),k=1,odr),&
                    (((z(i,j,k,n),i=1,odr),j=1,odr),k=1,odr)
    enddo
    close(103)
    nstart = nstart + nout

    nout = nz*nchord*nnormal
    open(103,file = 'mesh/pressureSide_block.xyz', form = 'unformatted')
    write(103)nout
    write(103)(odr, odr, odr,i=1, nout)
    do n = nstart + 1, nstart + nout
        write(103)  (((x(i,j,k,n),i=1,odr),j=1,odr),k=1,odr),&
                    (((y(i,j,k,n),i=1,odr),j=1,odr),k=1,odr),&
                    (((z(i,j,k,n),i=1,odr),j=1,odr),k=1,odr)
    enddo
    close(103)
    nstart = nstart + nout
    print*, nstart
endsubroutine
