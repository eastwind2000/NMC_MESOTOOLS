'reinit'
'open grid.ctl'
'open  dbz.ctl'
'enable print ./gmf/hf_radar_2007061200.33A_ne01.gmf'
'set parea 1 9.5 1 8'
'set xlopts 1 4 0.16'
'set ylopts 1 4 0.16'
'set mpdset q_cnworld'
'set lat  20 28  '
'set lon 105 114'
'set grads off'
'set grid  off'
'set ylint 1'
'set xlint 1'
 
'define a=oacres(gd(t=1),dbz.2,10,7,3,2,1)'
'set cmin 20'
'set cint 5'
'set gxout shaded'
'set rgb 30 51 51 255'
'set rgb 31 51 153 255'
'set rgb 32 0 255 255'
'set rgb 33 0 255 0'
'set rgb 34 51 204 51'
'set rgb 35 0 128 0'
'set rgb 36 255 255 0'
'set rgb 37 255 204 0'
'set rgb 38 255 153 51'
'set rgb 39 255 0 0'
'set rgb 40 204 51 0'
'set rgb 41 153 0 0'
'set rgb 42 255 0 255'
'set rgb 43 153 102 255'
'set rgb 44 255 255 255'
'set clevs        20 25 30 35 40 45 50 55 60 65 70'
'set ccols       0  34 35 36 37 38 39 40 41 42 43 44'
'd a'
'run cbarn 1 1'
 

****** ra1ch%x0=109.4561 ; ra1ch%y0=24.3569 ! Liuzhou 
 
 
x1=117.2578 ;   y1=31.8669
x2=115.8167  ;  y2=32.9167
'q w2xy 'x1' 'y1''
x01=subwrd(result,3)
y01=subwrd(result,6)
'q w2xy 'x2' 'y2''
x02=subwrd(result,3)
y02=subwrd(result,6)
'draw mark 5 'x01' 'y01' 0.1'
'draw mark 5 'x02' 'y02' 0.1'
L=178.696;s1=50;s2=100;s3=150;s4=200;s5=250
s0=3.7
m1=s0*s1/L;m2=s0*s2/L;m3=s0*s3/L;m4=s0*s4/L;m5=s0*s5/L
'set line 4 1 6'
'draw mark 2 'x01' 'y01' 'm1''
'draw mark 2 'x01' 'y01' 'm2''
'draw mark 2 'x01' 'y01' 'm3''
'draw mark 2 'x01' 'y01' 'm4''
'draw mark 2 'x01' 'y01' 'm5''
 
'set mpdset cnworld cnriver'
'set map  8  1  8'
'draw map'
'run cbarn 1 1'
'draw title dBz 2007061200.33A'
'print'
'printim ./gif/hf_radar_2007061200.33A_ne01.gif white '

 
