;; This plots a long version meant to be displayed in landscape mode.
;;
;; Carry out Gaussian smoothing with reflective boundaries along x-y
;;dimensions. Also write out the smoothed version for further visualization.

;; convert ra, dec, and z into comoving Mpc in the same coordinate
;; system as the reconstructions. REFER TO GEN_DACH_INPUT FOR CORRECT
;; X-Y ZERO POINTS!!
;;
;; Also generate little transverse maps along each slice, indicating
;; position on the sky.
;;
;; 2017-09-01 - Extended redshift range to z=2.05-2.55
;; 2017-09-09 - Read galaxies from a single text file. Different
;;              symbols for each survey (no error bars)
;; 2017-09-12 - Plot color bar twice, so that the figure can be split
;;              in two in the paper.
;; 2017-09-29 - Added horizontal colorbar at top, label RA range at
;;              each mini-sightline plot

Om = 0.31
Ol = 1.-Om

ra0 = 149.95
dec0 = 2.15
ra1 = 150.39
dec1 = 2.501

zmin = 2.05
zmax = 2.55
zmid = avg([zmin, zmax])
comdist0 = 2998. * comdis(zmin, Om, Ol)
dcomdist_dz = dcomdisdz(zmid, Om, Ol) * 2998.
mpc_amin = comdis(zmid, 0.31, 0.69)*2998. * !pi/180./60.

r_sm = 2.
mapfil = 'map_2017_v2.bin'
mapfil_sm = 'map_2017_v2_sm2.0.bin'
outsuf='slice_yz_Lpar2.0_Lperp2.5_v2_sm2.0_skewers_skymap'

;; Read skewer xy positions
readcol, '/Users/kheegan/lya/3d_recon/data/cl2017_redux/' + $
         'cl2017_valueadded_20170426_widez.txt', tomoid, gmag, conf, $
         ra_sk, dec_sk, sn1, sn2, sn3,f='x,l,f,f,x,f,f,f,f,f',/silent

skewerall = where((sn1 GE 1.2 OR sn2 GE 1.2 OR sn3 GE 1.2) AND conf GE 3.)
ra_sk = ra_sk[skewerall]
ra_sk = ra_sk - 150.
dec_sk = dec_sk[skewerall]


;; positions of known protoclusters
ra_pc =[150.09,150.00, 150.13, 150.24]
dec_pc = [2.25, 2.24, 2.37, 2.34]
z_pc = [2.095,2.45,2.47,2.506]

ra_npc =150.06
dec_npc = 2.20
z_npc = 2.30

ddeg_dchi = 180./!dpi / 2997./ comdis((zmin+zmax)/2., Om, Ol)
;; This factor is to account for the binning of the map grid
binfac = 2.

catdir = '/Users/kheegan/lya/3d_recon/ts/ts_cats/'
cosmos_specz_fil = '/Users/kheegan/scratch/cosmos/spec/' + $
                   'OBSERVED_TARGETS_15April2015_withHeader.dat'
mosdef2016_fil = catdir+'mosdef_zcat.16aug2016.fits'


nx = 60L
ny = 48L
nz = 876L

if file_test(mapfil_sm) then begin
   map_in = dblarr(nx*ny*nz)

   openr, 12, mapfil_sm
   readu, 12, map_in
   close, 12
   
   map_arr = dblarr(nx,ny,nz)
   
   ctr = 0L
   for ii = 0, nx-1 do begin
      for jj=0, ny-1 do begin
         for kk=0, nz-1 do begin
            map_arr[ii,jj,kk] = map_in[ctr]
            ctr++
         endfor
      endfor
   endfor
endif else begin
   
   map_in = dblarr(nx*ny*nz)
   
   openr, 12, mapfil
   readu, 12, map_in
   close, 12
   
   map_arr = dblarr(nx,ny,nz)
   
   ctr = 0L
   for ii = 0, nx-1 do begin
      for jj=0, ny-1 do begin
         for kk=0, nz-1 do begin
            map_arr[ii,jj,kk] = map_in[ctr]
            ctr++
         endfor
      endfor
   endfor
   
   
;; Carry out smoothing on a rebinned map
   maptmp = congrid(map_arr, nx/2., ny/2., nz/2.)
   
;; This function will take the array, pad it, and return a nz^3 array 
   map_pad = pad_reflect_map(maptmp,438)
   map_sm = smooth_3d(map_pad, r_sm)
   
   map_arr = map_sm[219-nx/4:219+nx/4-1, 219-ny/4:219+ny/4-1, *]

   map_arr = double(rebin(map_arr, nx, ny, nz))

   openw, 15, mapfil_sm
   writeu, 15, transpose(map_arr)
   close, 15
   
endelse 

;; Done with smoothing ---------------------------------------------

xpc = round(comdist0 * (ra_pc - ra0)/cos(0.5*abs(dec_pc+dec0)*!pi/180.) *!pi / 180. *binfac)
ypc = comdist0 * (dec_pc - dec0)*!pi / 180. * binfac
zpc = (z_pc - zmin) * dcomdist_dz * binfac

readcol, 'cat_tomoxyz_cl2017_uniq_v2_no3dhst.dat', xpos, ypos, zpos, $
         ori_id, f='f,f,f,a', skipline=4


;; DONE READING IN GALAXIES ---------------------------------------------

dz_gal = 3.3  ;; Delta comoving distance corresponding to 300km/s at z=2.3
dz_nir = 0.66 ;; Delta comoving distance corresponding to 60km/s at z=2.3

;; Read in skewer positions  ----------------------------------------
readcol, '/Users/kheegan/lya/3d_recon/map2017/list_tomo_input_2017.txt', $
         objid, xpos_sk, ypos_sk, zpos0_sk, zpos1_sk, f='l, x, x, x, x, f, f, f, f'


set_plot, 'ps'

device, file=outsuf+'.ps', /color, $
        /encap, ysize=73.5, xsize=51.08, /inch
!p.font=0
nslice=15
!p.multi = [0,1,nslice]
 

deltamin = -0.30;min(sliceyz)
deltamax = 0.20 ;max(sliceyz)

erase

;; Generate redshift axis
zmid = avg([zmin, zmax])
dcomdist_dz = dcomdisdz(zmid, Om, Ol)*2998.

y0 = 1./(2. * nslice) + (0.95/float(nslice)) * findgen(nslice)
dy_plot = 1./float(2.*nslice)

plotsym, 8, 1.3,/fill

;------------------------------------------------------------
for ii=0, nslice-1  do begin
   ;; These are the indices for the map array, NOT the correct
   ;; x-boundaries for discrete objects
   x0 = 4 * ii
   x1 = 4 * ii + 3
   sliceyz = transpose(reform(avg(map_arr[x0:x1,*,*],0)))
   fgridimg1=bytscl(sliceyz , top=240, min=deltamin, max=deltamax )
   fgridimg1 = -(fgridimg1 - 240.)

   x0obj = float(4*ii)
   x1obj = float(4*ii) + 4.
   
   galcut = where(xpos GE x0obj AND xpos LT x1obj, ngalhere)
   if ngalhere GT 0 then begin
      ygal = ypos[galcut]
      zgal = zpos[galcut]
      gal_id = ori_id[galcut]
   endif

   pccut = where(xpc GE x0obj AND xpc LT x1obj, npchere)
   if npchere GT 0 then begin
      ypctmp = ypc[pccut]
      zpctmp = zpc[pccut]
   endif

   position_win = [0.13, y0[ii]-dy_plot, 0.92, y0[ii]+dy_plot]

   xyratio = 48./876.
   xypos = cgAspect(xyratio, position=position_win)
   
   plot, findgen(876)/2., findgen(24), /nodata, xticks=0, yticks=0, $
         xsty=13, ysty=13, position=xypos,/norm, charsize=1.7

; These sequence of commands is to get the y-position of a few of the
; slices in order to determine the colorbar positions.
; First colorbar will be between the 4th and 5th slice, while the
; other one will be adjacent to the 12th slice... this is assuming
; that we want to split the plot into 8-panel and 7-panel 
if ii EQ 3 then ywindow3 = !y.window
if ii EQ 10 then ywindow10 = !y.window
if ii EQ 11 then ywindow11 = !y.window

   ywin0 = !y.window[0]
   ywin1 = !y.window[1]
;xyouts, 98, 22,textoidl('0 < x_{perp} (h^{-1} Mpc) < 2') , $
;        color=djs_icolor('red'), charthick=3, $
;        charsize=1.5, /data
   
   loadct, 33, /silent
   
   TV, fgridimg1, $
       !X.WINDOW(0), !Y.WINDOW(0), $
       XSIZE = (!X.WINDOW(1) - !X.WINDOW(0)), $
       YSIZE = !Y.WINDOW(1) - !Y.WINDOW(0), /NORM

   ra_min = ra0 + ii * 2.*ddeg_dchi
   ra_max = ra_min + 2.*ddeg_dchi
   
   ;xyouts, 133, 32, 'Slice #'+strtrim(ii+1,2)+': '+ $
   ;        string(ra_min,'(f7.3)')+' < RA (deg) < '+string(ra_max, '(f7.3)'), $
   ;     color=djs_icolor('red'), charthick=3, $
   ;     charsize=2., /data

      
   loadct, 0, /silent
   if ngalhere GT 0 then begin
      getzDeep = where(strmatch(gal_id, 'zDeep'),nzd)
      getVUDS = where(strmatch(gal_id, 'VUDS'),nvu)
      getMOSDEF = where(strmatch(gal_id, 'MOSDEF'),nmd)
      getZFIRE = where(strmatch(gal_id, 'ZFIRE'),nzf)
      getCLAMATO = where(strmatch(gal_id, 'CLAMATO'),ncl)
      ;; Squares for zCOSMOS
      if nzd GT 0 then begin
         oplot, zgal[getzDeep]/2., ygal[getzDeep]/2., $
                  psym=cgsymcat(15), thick=5, color=150, symsize=2.2
         oplot, zgal[getzDeep]/2., ygal[getzDeep]/2., $
                psym=cgsymcat('OPENSQUARE'), thick=5, color=0, symsize=2.2
      endif
      ;; Diamonds for VUDS
      if nvu GT 0 then begin
         oplot, zgal[getVUDS]/2., ygal[getVUDS]/2., $
                psym=cgsymcat(14), color=150, symsize=3.3
         oplot, zgal[getVUDS]/2., ygal[getVUDS]/2., $
                psym=cgsymcat('OPENDIAMOND'), color=0, symsize=3.3, thick=3
      endif
      ;; Circles for CLAMATO
      if ncl GT 0 then begin
         oplot, zgal[getCLAMATO]/2., ygal[getCLAMATO]/2., $
                psym=cgsymcat('FILLEDCIRCLE'), color=150, symsize=2.7
         oplot, zgal[getCLAMATO]/2., ygal[getCLAMATO]/2., $
                psym=cgsymcat('OPENCIRCLE'), color=0, symsize=2.7, thick=3
      endif
      ;; Up Triangles for ZFIRE
      if nzf GT 0 then begin
         oplot, zgal[getZFIRE]/2., ygal[getZFIRE]/2., $
             psym=cgsymcat('FILLEDUPTRIANGLE'), color=150, symsize=2.7
         oplot, zgal[getZFIRE]/2., ygal[getZFIRE]/2., $
             psym=cgsymcat('OPENUPTRIANGLE'), color=0, symsize=2.7, thick=3
      endif
      ;; Down triangles for MOSDEF
      if nmd GT 0 then begin
         oplot, zgal[getMOSDEF]/2., ygal[getMOSDEF]/2., $
             psym=cgsymcat('FILLEDDOWNTRIANGLE'), color=150, symsize=2.7
         oplot, zgal[getMOSDEF]/2., ygal[getMOSDEF]/2., $
             psym=cgsymcat('OPENDOWNTRIANGLE'), color=0, symsize=2.7, thick=3
      endif
   endif

   if npchere GT 0 then begin
      oplot, zpctmp/2., ypctmp/2., color=0, psym=cgsymcat(46),symsize=9.5
      oplot, zpctmp/2., ypctmp/2., color=255, psym=cgsymcat(45),symsize=9.5,thick=8
   endif

   ;; Plot skewer positions
   skewercut = where(xpos_sk*binfac GE x0obj AND xpos_sk*binfac LT x1obj AND $
                     ypos_sk GE 0. AND ypos_sk LT 24., nskewer)

   if nskewer GT 0 then begin
      print, 'Skewers in slice #',ii+1
      for iskewer=0, nskewer-1 do begin
         id_skw = objid[skewercut[iskewer]]
         yskw = ypos_sk[skewercut[iskewer]]
         zskw0 = zpos0_sk[skewercut[iskewer]]
         zskw1 = zpos1_sk[skewercut[iskewer]]
         print, id_skw, xpos_sk[skewercut[iskewer]], yskw 
         oplot, [zskw0, zskw1], [yskw,yskw], thick=4, color=djs_icolor('white')
      endfor
   endif

   ;; Plot circle around z=2.3 galaxies
   if ra_min LE ra_npc AND ra_max GT ra_npc then begin
      print, 'plotting z=2.30 non-protocluster'
      ynpc = comdist0 * (dec_npc - dec0)*!pi / 180. * binfac
      znpc = (z_npc - zmin) * dcomdist_dz * binfac
;      tvellipse, 4., 4., znpc/2,ynpc/2., /data, thick=5, color=djs_icolor('pink')
   endif
      
   axis, yaxis=1, charsize=5., charthick=4, yran=[0,24],ysty=1, $
         ytit=textoidl('y_{perp} (h^{-1} Mpc)'),ytickv=[0,5,10,15,20], $
         yticks=4, yminor=5

   decran = dec0 + [0.,24.]/(2998.*comdis((zmin+zmax)/2., Om, Ol))*180./!dpi
   axis, yaxis=0, charsize=5., charthick=4, yran=decran, /ysty,ytickformat='(A1)';, $
  ;       ytit='Dec (deg)'
   
   comdis0 = comdis(zmin, Om, Ol)*2998.
   comdis1 = comdis0 +438.
   axis, xaxis=1, xran=[comdis0, comdis1], charthick=4, charsize=5., $
         xsty=1, xtit=textoidl('Comoving Distance (h^{-1} Mpc)')
   
   z0 = zmin
   z1 = z0 + 438./dcomdist_dz
   axis, xaxis=0, xran=[z0, z1], charthick=4, charsize=5., $ 
         xsty=1, xtit='Redshift'

;; Plot little sky map
   position_sky=[0.005,ywin0, 0.18, ywin1]

   skyratio = 24./30.
   skypos = cgAspect(skyratio, position=position_sky,align='center')
   ychsize = !d.y_ch_size
   yvsize  = !d.y_vsize
   
   dslice_ra = 2./ mpc_amin / 60.
   raslice0 = ra0 + (ii)*dslice_ra - 150.
   titstr = string(ra0+(ii)*dslice_ra , '(f7.3)')+textoidl('^\circ')+'<RA<'+ $
            string(ra0+(ii+1)*dslice_ra, '(f7.3)')+textoidl('^\circ')
   plot, [avg(ra0+ra1)/2.], [avg(dec0+dec1)/2.], /nodata, /norm, $
         position=skypos,xsty=13, ysty=13,/noerase,xran=[ra1-150.,ra0-150.], $
         yran=[dec0,dec1], title=titstr,charsize=3., charthick=4
   ;; Plot slice footprint
   x_vert = [raslice0, raslice0+dslice_ra, raslice0+dslice_ra, raslice0]
   y_vert = [dec0,    dec0,    dec1, dec1]
   polyfill, x_vert, y_vert, color=150, /data
   ;; plot sightlines
   plotsym, 0, 0.65,/fill,color=djs_icolor('red')
   oplot, ra_sk, dec_sk, psym=8
   
   xran=[ra1-150., ra0-150.]
   axis, yaxis=0, charsize=4., charthick=4, yran=decran, /ysty, $
         ytit='Dec (deg)'
   axis, yaxis=1, charsize=4., charthick=4, yran=decran, /ysty, $
         YTICKFORMAT="(A1)"
   axis, xaxis=0, charsize=4., charthick=4, xran=xran, /xsty, $
         xtit='RA+150.0 (deg)', xtickinterval=0.1, xtickformat='(f4.1)' ;$
   axis, xaxis=1, charsize=4., charthick=4, xran=xran, /xsty, $
         xtickformat='(A1)'
   
endfor

;-------------------------------------------------------------

   ;; Now plot color bar(s)
   multiplot, /reset
   !p.multi = replicate(0,5)
   
   loadct, 33, /silent
   tickpts =deltamin + (deltamax - deltamin)/ 4.*findgen(5)
   ticknames = string(tickpts, '(f5.2)')

   ; The two vertical bars
   ycen_tmp = (ywindow10[1]+ywindow11[0])/2.
   barpos = [0.955, ycen_tmp - 0.06, 0.97, ycen_tmp + 0.06]
   colorbar, /vertical,position=barpos, $
             divisions=4, ticknames=ticknames, ncolors=240, $
             charsize=2.5, charthick=4,/right,/norm, /invert
   xyouts, 0.96, ycen_tmp+0.065, textoidl('\delta^{rec}_F'), charsize=3., $
           charthick=4, /norm

   ycen_tmp = (ywindow3[0]+ywindow3[1])/2.
   barpos = [0.955, ycen_tmp - 0.06, 0.97, ycen_tmp + 0.06]
   colorbar, /vertical,position=barpos, $
             divisions=4, ticknames=ticknames, ncolors=240, $
             charsize=2.5, charthick=4,/right,/norm, /invert
   xyouts, 0.96, ycen_tmp+0.065, textoidl('\delta^{rec}_F'), charsize=3., $
           charthick=4, /norm

   ; One horizontal bar at top
   xcen_tmp = (0.14+0.92)/2.
   ybottom_tmp = 0.975
   ythick = 0.008
   barpos = [xcen_tmp-0.08, ybottom_tmp, xcen_tmp+0.08, ybottom_tmp+ythick]
   colorbar, /horizontal, position=barpos, $
             divisions=4, ticknames=ticknames, ncolors=240, $
             charsize=2.5, charthick=4,/norm, /invert
   xyouts, xcen_tmp-0.003, ybottom_tmp+ythick+0.005, $
           textoidl('\delta^{rec}_F'), charsize=3., $
           charthick=4, /norm

device, /close

set_plot, 'x'

stop

end
