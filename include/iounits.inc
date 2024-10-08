      integer iunit,ounit,logunit,datunit,keyunit,consunit
      parameter (iunit=1,ounit=2,logunit=3,datunit=4,keyunit=5,
     *           consunit=6)

      character*4 short(-2:2,-2:2)
      data short / 
     *  '-2-2','-1-2','#0-2','+1-2','+2-2',
     *  '-2-1','-1-1','#0-1','+1-1','+2-1',
     *  '-2#0','-1#0','#0#0','+1#0','+2#0',
     *  '-2+1','-1+1','#0+1','+1+1','+2+1',
     *  '-2+2','-1+2','#0+2','+1+2','+2+2' /

      character*6 parname(9)
      data parname /
     *  'betad ','dprper','pratio',
     *  'dcgper','cratio','betacc',
     *  'cu(0) ','cu(2) ','rex   ' /
