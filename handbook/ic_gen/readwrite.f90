      integer nx,ny,nz,i,j,k
      double precision u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
      double precision ubin(nx*ny*nz),vbin(nx*ny*nz),wbin(nx*ny*nz)

C     Open the binary flow fields (in this case the initial field)

      open(unit=530,file='CBC_64_qr_U.bin',access='stream')
      open(unit=531,file='CBC_64_qr_V.bin',access='stream')
      open(unit=532,file='CBC_64_qr_W.bin',access='stream')

C     Store the data files in arrays

      read (530) (ubin(ibin), ibin=1,nx*ny*nz)
      read (531) (vbin(ibin), ibin=1,nx*ny*nz)
      read (532) (wbin(ibin), ibin=1,nx*ny*nz)

C     Close the files

      close(530)
      close(531)
      close(532)

C     The initial values of u, v, w and p are set
 
      do k=1,nz
      do j=1,ny
      do i=1,nx
         u(i,j,k) = ubin((k-1)*ny*nx + (j-1)*nx + i)
         v(i,j,k) = vbin((k-1)*ny*nx + (j-1)*nx + i)
         w(i,j,k) = wbin((k-1)*ny*nx + (j-1)*nx + i)
         p(i,j,k) = 0.0D0
      enddo
      enddo
      enddo
        </pre>
        <p>
            and written using
        </p>
        <pre style="background-color: LightGrey">
      integer i,j,k,ibin
      double precision arru(nx*ny*nz),arrv(nx*ny*nz),arrw(nx*ny*nz)

C     Make the data arrays

      do k=1,nz
      do j=1,ny
      do i=1,nx
         arru((k-1)*nx*ny + (j-1)*nx + i) = u(i,j,k)
         arrv((k-1)*nx*ny + (j-1)*nx + i) = v(i,j,k)
         arrw((k-1)*nx*ny + (j-1)*nx + i) = w(i,j,k)
      enddo
      enddo
      enddo

C     Write the data array (in this case at time 98)

      open(unit=678,file='CBC_64_dsm_98_U.bin',access='stream')
      open(unit=679,file='CBC_64_dsm_98_V.bin',access='stream')
      open(unit=680,file='CBC_64_dsm_98_W.bin',access='stream')

C     Write the binary flow fields

      write (678) (arru(ibin),ibin=1,nx*ny*nz)
      write (679) (arrv(ibin),ibin=1,nx*ny*nz)
      write (680) (arrw(ibin),ibin=1,nx*ny*nz)

C     close the binary data files

      close(678)
      close(679)
      close(680)
