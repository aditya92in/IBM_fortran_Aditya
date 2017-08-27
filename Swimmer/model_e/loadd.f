      subroutine loadd(in)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer reclengte,in
      reclengte=8
      if(in.eq.0) then
      open(15,file='data_start.txt',form='unformatted')
      read(15)Unew,Vnew,Uold,Vold,p,pold,Cnew,Cold
      close(15)
      endif
      if(in.eq.1) then
      open(15,file='data_start.txt',form='unformatted')
      write(15)Unew,Vnew,Uold,Vold,p,pold,Cnew,Cold
      close(15)
      endif
      return
      end
