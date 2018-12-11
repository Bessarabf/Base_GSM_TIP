       subroutine intut(uts,tut,sut,dut,hut,n,nut,
     *                  et,ed,eh,ieisc,isat)

      dimension tut(isat,100),sut(isat,100),dut(isat,100),hut(isat,100)
      ut=uts/3600.
      if(ut.lt.tut(1,n).and.(ut+24.).gt.tut(nut,n))then
        ieisc=0
        return
      end if
      i=1
    1 continue
      ip=i+1
      if(ip.gt.nut)then
        ieisc=0
      else
        if(ut.ge.tut(i,n).and.ut.le.tut(ip,n))then
          ieisc=1
          if(nut.eq.2)then
            et=sut(1,n)
            ed=dut(1,n)
            eh=hut(1,n)
          else
            a=tut(i,n)
            b=tut(ip,n)
            del=(ut-a)/(b-a)
            a=sut(i,n)
            b=sut(ip,n)
            et=a+(b-a)*del
            a=dut(i,n)
            b=dut(ip,n)
            ed=a+(b-a)*del
            a=hut(i,n)
            b=hut(ip,n)
            eh=a+(b-a)*del
          end if
        else
          if((ut+24.).ge.tut(i,n).and.(ut+24.).le.tut(ip,n))then
            ieisc=1
            if(nut.eq.2)then
              et=sut(1,n)
              ed=dut(1,n)
              eh=hut(1,n)
            else
              a=tut(i,n)
              b=tut(ip,n)
              del=(ut+24.-a)/(b-a)
              a=sut(i,n)
              b=sut(ip,n)
              et=a+(b-a)*del
              a=dut(i,n)
              b=dut(ip,n)
              ed=a+(b-a)*del
              a=hut(i,n)
              b=hut(ip,n)
              eh=a+(b-a)*del
            end if
          else
            i=ip
            goto1
          end if
        end if
      end if
      return
      end
                                                                                                                                                                                                                   
