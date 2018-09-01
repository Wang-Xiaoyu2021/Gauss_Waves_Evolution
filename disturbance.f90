subroutine inv_mat(A,invA)
   implicit none
   integer(kind = 8)        :: i,j,k,l,is(0:100),js(0:100)
   complex(kind = 8)        :: A(0:100,0:100),invA(0:100,0:100),A_old(0:100,0:100)
   complex(kind = 8)        :: t
   integer(kind = 8)        :: d
   print*,"test3"
   A_old=A
   l=1
   do k=0,100
      d=1D-3
      do i=k,100
         do j=k,100
            if (abs(a(i,j)).gt.d) then
               d=abs(a(i,j))
               is(k)=i
               js(k)=j
            endif
         enddo
      enddo
      if ((d+1.0).eq.1.0) then
         l=0
      endif
      do j=0,100
         t=a(k,j)
         a(k,j)=a(is(k),j)
         a(is(k),j)=t
      enddo
      do i=0,100
         t=a(i,k)
         a(i,k)=a(i,js(k))
         a(i,js(k))=t
      enddo
      a(k,k)=1/a(k,k)
      do j=0,100
         if(j.ne.k) then
            a(k,j)=a(k,j)*a(k,k)
         endif
      enddo
      do i=0,100
         if (i.ne.k) then
            do j=0,100
               if (j.ne.k) then
                  a(i,j)=a(i,j)-a(i,k)*a(k,j)
               endif
            enddo
         endif
      enddo
      do i=0,100
         if(i.ne.k) then
            a(i,k)=-a(i,k)*a(k,k)
         endif
      enddo
   enddo
   do k=100,0,-1
      do j=0,100
         t=a(k,j)
         a(k,j)=a(js(k),j)
         a(js(k),j)=t
      enddo
      do i=0,100
         t=a(i,k)
         a(i,k)=a(i,is(k))
         a(i,is(k))=t
      enddo
   enddo
   invA=A
   A=A_old
   end subroutine
