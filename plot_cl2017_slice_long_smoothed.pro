;; This plots a long version meant to be displayed in landscape mode.
;;
;; Carry out Gaussian smoothing with reflective boundaries along x-y
;;dimensions. Also write out the smoothed version for further visualization.

;; convert ra, dec, and z into comoving Mpc in the same coordinate
;; system as the reconstructions. REFER TO GEN_DACH_INPUT FOR CORRECT
;; X-Y ZERO POINTS!!

Om = 0.31
Ol = 1.-Om

ra0 = 149.95
dec0 = 2.15

zmin = 2.15
zmax = 2.55
zmid = avg([zmin, zmax])
comdist0 = 2997. * comdis(zmin, Om, Ol)
dcomdist_dz = dcomdisdz(zmid, Om, Ol) * 2997.

r_sm = 4.
mapfil_sm = 'map_2017_sm4.0.bin'

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

outsuf='slice_yz_Lpar2.2_Lperp2.0_v0_sm2.0_skewers'
mapfil = 'map.bin'


nx = 60L
ny = 48L
nz = 680L

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
   map_pad = pad_reflect_map2016(maptmp)
   map_sm = smooth_3d(map_pad, r_sm)
   
   map_arr = map_sm[170-nx/4:170+nx/4-1, 170-ny/4:170+ny/4-1, *]

   map_arr = double(rebin(map_arr, nx, ny, nz))

   openw, 15, mapfil_sm
   writeu, 15, transpose(map_arr)
   close, 15
   
endelse 

;; Done with smoothing ---------------------------------------------

print, 'Reading in MOSDEF redshifts'
cat_md = mrdfits(mosdef2016_fil,1)

;; We only want redshifts within COSMOS, then match in ra+dec 
getcosmos = where(strmatch(cat_md.field,'*COSMOS*') AND $
                  cat_md.z_mosfire GE 2.,nmosdef)
cat_md = cat_md[getcosmos]

;; Read ZFIRE catalog --------------------------------------------
readcol, catdir+'KG_ZFIRE_catalog.txt', ra_zf, dec_zf, z_zf, conf_zf, $
         f='x,f,f,f,l'

;; Find and remove duplicates between MOSDEF and ZFIRE that are
;; within ~0.4" of each other
nmatch = djs_angle_nmatch(cat_md.ra, cat_md.dec, ra_zf, dec_zf, $
                          1.e-4)
remove, where(nmatch GT 0), ra_zf, dec_zf, z_zf

;; Combine MOSDEF and ZFIRE redshifts
ra_mosf = [cat_md.ra, double(ra_zf)]
dec_mosf = [cat_md.dec, double(dec_zf)]
z_mosf = [cat_md.z_mosfire, double(z_zf)]

zcut_tmp = where(z_mosf LE 2.15 OR z_mosf GT 2.54)
remove, zcut_tmp, ra_mosf, dec_mosf, z_mosf

;; Now read in COSMOS redshifts and throw out low-z stuff and
;; low-confidence 
readcol, cosmos_specz_fil, ori_id,ra_gal, dec_gal,  $
         zspec, qflag, f='a,f,f,x,f,f', skip=120
zcut = where(zspec LE 2.15 OR zspec GT 2.54)
remove, zcut, ra_gal, dec_gal, zspec, qflag, ori_id
;; Remove low-quality redshifts
qualcut = where(qflag LT 3 OR qflag GE 10)
remove, qualcut, ra_gal, dec_gal, zspec, qflag, ori_id

;; Read in CLAMATO galaxies with similar cuts. 
;readcol, '/Users/kheegan/lya/3d_recon/data/cl2016_redux/cpilot_insp_20160711_v2.txt', $
;         fname, catnum, zcl, qual_cl, $
;         f='a,l, f,f'
;zcut = where(zcl LE 2.2 OR zcl GT 2.5 OR qual_cl LE 2)
;remove, zcut, catnum, zcl, qual_cl, fname
;ncl = n_elements(catnum)
;; Read in CLAMATO catalog
;cat_fil = '/Users/kheegan/lya/3d_recon/ts/pilot/' + $
;          'cosmos_ts_pilot_mastercat.fits'
;cat = mrdfits(cat_fil, 1, /silent)
;ra_cl = cat[catnum].ra
;dec_cl = cat[catnum].dec

;dupl_list = []
;for ii=0, ncl-1 do begin
;   ratmp = (ra_cl[ii])[0]
;   dectmp = (dec_cl[ii])[0]
;   deltapos = sqrt( (ratmp - ra_gal)^2 + (dectmp - dec_gal)^2)
;   matchpos = where(deltapos LE 2.e-4, nmat)
;   if nmat GT 0 then dupl_list = [dupl_list, ii]
;endfor
;remove, dupl_list, ra_cl, dec_cl, zcl, catnum
;; Finally, check for cases where a source is targeted in 2 separate
;; masks 
;catsort = sort(catnum)
;catnum = catnum[catsort]
;ra_cl = ra_cl[catsort]
;dec_cl = dec_cl[catsort]
;zcl = zcl[catsort]

;uniqcat = uniq(catnum)

;ra_gal = [ra_gal, ra_cl[uniqcat]]
;dec_gal = [dec_gal, dec_cl[uniqcat]]
;zspec = [zspec, zcl[uniqcat]]

xpos = comdist0 * (ra_gal - ra0)/cos(0.5*abs(dec_gal+dec0)*!pi/180.) *!pi / 180. *binfac
ypos = comdist0 * (dec_gal - dec0)*!pi / 180. * binfac
zpos = (zspec - zmin) * dcomdist_dz * binfac

xpc = round(comdist0 * (ra_pc - ra0)/cos(0.5*abs(dec_pc+dec0)*!pi/180.) *!pi / 180. *binfac)
ypc = comdist0 * (dec_pc - dec0)*!pi / 180. * binfac
zpc = (z_pc - zmin) * dcomdist_dz * binfac

xmos= comdist0 * (ra_mosf - ra0)/cos(0.5*abs(dec_mosf+dec0)*!pi/180.) *!pi / 180. *binfac
ymos= comdist0 * (dec_mosf - dec0)*!pi / 180. * binfac
zmos= (z_mosf - zmin) * dcomdist_dz * binfac

in_vol = where(xpos GE 1. AND xpos LE 59. AND ypos GE 1. AND ypos LE 47., n_invol)
print, n_invol, ' galaxies within map volume'

ra_gal = ra_gal[in_vol]
dec_gal = dec_gal[in_vol]
zspec = zspec[in_vol]
ori_id = ori_id[in_vol]

xpos = round(xpos[in_vol])
ypos = ypos[in_vol]
zpos = zpos[in_vol]

in_vol_mosf = where(xmos GE 1. AND xmos LE 47. AND ymos GE 1. $
                    AND ymos LE 47., n_invol_mosf)
print, n_invol_mosf, ' MOSFIRE galaxies within map volume'
xmos = round(xmos[in_vol_mosf])
ymos = ymos[in_vol_mosf]
zmos = zmos[in_vol_mosf]


;; DONE READING IN GALAXIES ---------------------------------------------

dz_gal = 3.3  ;; Delta comoving distance corresponding to 300km/s at z=2.3
dz_nir = 0.66 ;; Delta comoving distance corresponding to 60km/s at z=2.3

;; Read in skewer positions  ----------------------------------------
readcol, '/Users/kheegan/lya/3d_recon/map2017/list_tomo_input_2017.txt', $
         objid, xpos_sk, ypos_sk, zpos0_sk, zpos1_sk, f='l, x, x, x, x, f, f, f, f'


set_plot, 'ps'

device, file=outsuf+'.ps', /color, $
        /encap, ysize=64.28, xsize=34.08, /inch
!p.font=0
nslice=15
!p.multi = [0,1,nslice]
 

deltamin = -0.25;min(sliceyz)
deltamax = 0.15 ;max(sliceyz)

erase

;; Generate redshift axis
zmid = avg([zmin, zmax])
dcomdist_dz = dcomdisdz(zmid, Om, Ol)*2997.

y0 = 1./(2. * nslice) + (0.98/float(nslice)) * findgen(nslice)
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
   
   galcut = where(xpos GE x0obj AND xpos LT x1obj AND $
                  strmatch(ori_id,'MOSFIRE') EQ 0, ngalhere)
   if ngalhere GT 0 then begin
      ygal = ypos[galcut]
      zgal = zpos[galcut]
   endif

   nircut = where(xmos GE x0obj AND xmos LT x1obj, nnirhere)
   if nnirhere GT 0 then begin
      ynir = ymos[nircut]
      znir = zmos[nircut]
   endif

   pccut = where(xpc GE x0obj AND xpc LT x1obj, npchere)
   if npchere GT 0 then begin
      ypctmp = ypc[pccut]
      zpctmp = zpc[pccut]
   endif

   position_win = [0.04, y0[ii]-dy_plot, 0.9, y0[ii]+dy_plot]

   xyratio = 48./680.
   xypos = cgAspect(xyratio, position=position_win)
   
   plot, findgen(680)/2., findgen(24), /nodata, xticks=0, yticks=0, $
         xsty=13, ysty=13, position=xypos,/norm, charsize=1.7
   
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
   
   xyouts, 133, 32, 'Slice #'+strtrim(ii+1,2)+': '+ $
           string(ra_min,'(f7.3)')+' < RA (deg) < '+string(ra_max, '(f7.3)'), $
        color=djs_icolor('red'), charthick=3, $
        charsize=2., /data

      
   loadct, 0, /silent
   if ngalhere GT 0 then begin
      oploterror, zgal/2., ygal/2., $
                  replicate(dz_gal, ngalhere), replicate(0., ngalhere), $
                  psym=cgsymcat(15), thick=5, color=150, symsize=2.
   endif

   if nnirhere Gt 0 then begin
      oplot, znir/2., ynir/2., color=150,  psym=cgsymcat(15), symsize=2.
   endif

   if npchere GT 0 then begin
      oplot, zpctmp/2., ypctmp/2., color=0, psym=cgsymcat(46),symsize=6.
      oplot, zpctmp/2., ypctmp/2., color=255, psym=cgsymcat(45),symsize=6.,thick=8
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
      tvellipse, 4., 4., znpc/2,ynpc/2., /data, thick=5, color=djs_icolor('pink')
   endif
      
   axis, yaxis=0, charsize=3.8, charthick=3, yran=[0,24],ysty=1, $
         ytit=textoidl('y_{perp} (h^{-1} Mpc)'),ytickv=[0,4,8,12,16,20,24], $
         yticks=6, yminor=4

   decran = dec0 + [0.,24.]/(2997.*comdis((zmin+zmax)/2., Om, Ol))*180./!dpi
   axis, yaxis=1, charsize=3.8, charthick=3, yran=decran, /ysty, $
         ytit='Dec (deg)'
   
   comdis0 = comdis(zmin, Om, Ol)*2997.
   comdis1 = comdis0 + 340.
   axis, xaxis=1, xran=[comdis0, comdis1], charthick=3, charsize=3.8, $
         xsty=1, xtit=textoidl('Comoving Distance (h^{-1} Mpc)')
   
   z0 = zmin
   z1 = z0 + 340./dcomdist_dz
   axis, xaxis=0, xran=[z0, z1], charthick=3, charsize=3.8, $ 
         xsty=1, xtit='z'

endfor

;-------------------------------------------------------------

   ;; Now plot color bar(s)
   multiplot, /reset
   !p.multi = replicate(0,5)
   
   loadct, 33, /silent
   tickpts =deltamin + (deltamax - deltamin)/ 4.*findgen(5)

   ticknames = string(tickpts, '(f5.2)')
   colorbar, /vertical,position=[0.95,0.37,0.965,0.495] $
             ,divisions=4, ticknames=ticknames, ncolors=240, $
             charsize=2.3, charthick=3,/right,/norm, /invert
   xyouts, 0.954, 0.50, textoidl('\delta^{rec}_F'), charsize=3., $
           charthick=3, /norm

device, /close

set_plot, 'x'

stop

end
