;; Generate input files from CLAMATO data, suitable for input to
;; dachshund for reconstruction.
;;
;; New for 2017A:
;; - New map zero-points
;; - No longer manually discard spectra; instead flag with negative
;;   numbers in catalog
;  - reduced restframe range to >1185A to avoid SiII lines at ~1192A  

Om = 0.31
Olamb = 1.-Om

datadir = '~/lya/3d_recon/data/cl2017_redux/'
specdir = datadir+'spec_v0/'

;bad_id = [12960, 15190]

;; This defines the Cartesian zeros for the output map
ra0 = 149.95  
dec0 = 2.15   

zmin = 2.15
zmax = 2.55
waveminlya = 1215.67 * (1.+zmin)
wavemaxlya = 1215.67 * (1.+zmax)
;zmid = avg([zmin, zmax])
zmid = 2.35


; transverse comoving distance = angular separation(in rad) *
; comoving distance. So evaluate comoving distance
comdist = 2997. * comdis(zmid, Om, Olamb)
; For LOS distances we use a fixed dD/dz evaluated at zmid
dcomdist_dz = dcomdisdz(zmid, Om, Olamb) * 2997.

insp_fil =  datadir+'cl2017_valueadded_20170426.txt'

;; Read in CLAMATO catalog
;cat_fil = '/Users/kheegan/lya/3d_recon/ts/pilot/' + $
;          'cosmos_ts_pilot_mastercat.fits'
;cat = mrdfits(cat_fil, 1, /silent)

readcol, insp_fil, specfil, catnum, mag, zconf, zsp, ra, dec, $
         snrlya1, snrlya2, f='a, l, f, f, f,f ,f , f,f', /silent

qualcut = where((snrlya1 GE 1.2 OR snrlya2 GE 1.2) AND zconf GE 3., $
                nsel, complement=failcut) 

remove, failcut, specfil, catnum, zsp, zconf, snrlya1,snrlya2, ra, dec, mag

;; Manually remove objects that are deemed bad
;match, catnum, bad_id, sub1, sub2
;if sub1 NE [-1] then remove, sub1, specfil, catnum, mag, zconf, $
;                             zsp, ra, dec, snrlya1, snrlya2

nsel = n_elements(specfil)

;; Read in continuum template. We also define a finite spectrum up to
;; 1400A because mean-flux regulation might act funny if it's
;; just the LyaF segment
contfil = 'contlya_musyc_allstack.dat'
readcol, contfil, wavecont1, cont1, /silent
wavecont = 1000. + findgen(400)
cont = replicate(avg(cont1), 400)
wc0 = min(wavecont1)
wc1 = max(wavecont1)
intrange = where(wavecont GE wc0 AND wavecont LE wc1)
cont[intrange] = interpol(cont1, wavecont1, wavecont[intrange])

;; prepare output vectors
x_out = []
y_out = []
z_out = []
delta_out = []
noise_out = []

xlist = []
ylist = []
z0list = []
z1list = []
i0list = []
i1list = []
weightlist = []

ctr = 0L

close, 19
openw, 19, 'list_tomo_input_2017.txt'
printf, 19, '# TomoID, zspec, gmag, ra, dec, delta_x (Mpc/h), ' + $
        'delta_y (Mpc/h), z0 (Mpc/h), z1 (Mpc/h), i_start, i_end'

for ii=0, nsel-1 do begin
   
   infil_tmp = specdir+specfil[ii]
      
   ;; Parse object name
   namesplit = strsplit(specfil[ii], '_', /extract)
   if strmatch(specfil[ii], 'c16_*') OR strmatch(specfil[ii], 'c17_*') then $
      objname = namesplit[2] + '_'+ namesplit[3] $
                else objname = namesplit[1] + '_'+ namesplit[2]

   print, objname
   
   ztmp = (zsp[ii])[0]
   
   catnum_tmp = (catnum[ii])[0]
   ra_tmp = ra[ii]
   dec_tmp = dec[ii]
   gmag_tmp = mag[ii]

   ;; Is this object a QSO?
   if zconf[ii] GT 10. then qsoflag = 1 else qsoflag = 0

   print, infil_tmp
   ;; Read in spectrum
   wave = mrdfits(infil_tmp, 2, /silent)
   flux = mrdfits(infil_tmp, 0, /silent)
   sig  = mrdfits(infil_tmp, 1, /silent)
   ivar = (sig GT 0)/(sig^2 + (sig EQ 0))
   
   ;; Check if mask files exist
   maskfil = file_search(datadir+'masks/mask_'+objname+'_F.fits', count=nmask)
   if nmask EQ 1 then begin
      mask = mrdfits(maskfil, 0, /silent)
      ivar = ivar * float(mask)
      sig = sig * float(mask)
   endif

   ;; Ditch the first few pixels, which seem problematic in some
   ;; spectra
   remove, lindgen(4), wave, flux, sig, ivar
   
   if qsoflag EQ 0 then begin
      ;; Mask regions near absorption lines
      maskvec = replicate(1., n_elements(wave))
      masklines = where((wave GT (1175.7*(1.+ztmp)-7.5) AND $
                         wave LT (1175.7*(1.+ztmp)+7.5)) OR $
                        (wave GT (1134.4*(1.+ztmp)-7.5) AND $
                         wave LT (1134.4*(1.+ztmp)+7.5))OR $
                        (wave GT (1085.7*(1.+ztmp)-7.5) AND $
                         wave LT (1085.7*(1.+ztmp)+7.5)) OR $
                        (wave GT (1192*(1.+ztmp)-7.5) AND $
                         wave LT (1192*(1.+ztmp)+7.5)) )
      
      maskvec[masklines] = 0.
      
      ivar =ivar * maskvec
      sig = sig * maskvec
   endif

   ;; Skip objects if insufficient # of LyaF pix present
   wavemintmp = 1041.*(1.+ztmp) > waveminlya
   wavemaxtmp = 1185.*(1.+ztmp) < wavemaxlya
   forestrange = where(wave GE wavemintmp AND wave LT wavemaxtmp AND $
                       ivar GT 0., npix_forest)
   if npix_forest LT 20 then begin
      print, 'Skipping '+objname+' due to insufficient pixels'
      continue
   endif
   
   ;; Fit continuum... first check if there's a hand-fitted continuum 
   contfil = file_search(datadir+'cont/cont_'+objname+'_F.fits', $
                         count=ncont)
   if ncont EQ 1 then begin
      cont_tmp = mrdfits(contfil, 0)
      wavemintmp = 1041.*(1.+ztmp) > waveminlya
      wavemaxtmp = 1185.*(1.+ztmp) < wavemaxlya
      
      
   endif else begin
;;for quasars we do MF-PCA
      if qsoflag EQ 1 then begin
         wavemintmp = 1041.*(1.+ztmp) > waveminlya
         wavemaxtmp = 1185.*(1.+ztmp) < wavemaxlya
         cont_tmp = pcaspec_func(wave/(1.+ztmp), flux, ivar=ivar, /dr7eigen, zqso=ztmp)
      endif else begin
         
         ;; For galaxies, we do MF-regulation if there is enough forest
         ;; pathlength (~250A)
         forestrange = where(wave GE 1040.*(1.+ztmp) AND wave LT 1185.*(1.+ztmp))
         waveforest = wave[forestrange]
         
         dwave_forest = max(waveforest) - (min(waveforest) > 3600.)
         
         if dwave_forest GT 250. then begin
            cont_tmp = interpol(cont, wavecont*(1.+ztmp), wave)
            normrange= where(wave/(1.+ztmp) GE 1050 AND $
                             wave/(1.+ztmp) LE 1190)
            cont_tmp = cont_tmp * avg(flux[normrange])/avg(cont_tmp[normrange])
            cont_tmp = mfreg(wave, flux, cont_tmp, ivar=ivar, $
                             lyaf_range=[1040, 1185.], zqso=ztmp, wavemin=3800.)
         endif else begin
            
            ;; CNR_FOREST is for estimating SNR but also estimates a
            ;; power-law for galaxies
            void = cnr_forest(wave, flux, sig, zem=ztmp, cont_out=cont_tmp)
            cont_tmp[where(wave/(1.+ztmp) LT 1200)] *= 0.95
         endelse
         
         wavemintmp = 1040.*(1.+ztmp) > waveminlya
         wavemaxtmp = 1185.*(1.+ztmp) < wavemaxlya
      endelse
   endelse 
      
;   xran=(1.+ztmp)*[1030., 1350.]
;   yvals = flux[sort(flux)]
;   yran=[yvals[round(0.05*n_elements(yvals))], yvals[round(0.95*n_elements(yvals))]]
;   plot, wave, smooth(flux,3), xran=xran, tit=objname, /xsty, yran=yran, /ysty
;   oplot, wave, cont_tmp, color=djs_icolor('green')
;   oplot, wave, sig, color=djs_icolor('red')
   
;   vline, 1040 * (1.+ztmp)
;   vline, 1185 * (1.+ztmp)
;   vline, 1085 * (1.+ztmp);

;   vline, wavemintmp, linesty=2
;   vline, wavemaxtmp, linesty=2
   
;   dummy = ''
;   read, dummy

   forestcut = where(wave GE wavemintmp AND wave LE wavemaxtmp, npix_tmp)

   zvec = wave[forestcut]/ 1215.67 - 1.
   losvec = (zvec - zmin) * dcomdist_dz
   
   gcirc, 2, ra_tmp, dec_tmp, ra0, dec_tmp , delta_ra
   gcirc, 2, ra_tmp, dec_tmp, ra_tmp, dec0, delta_dec

   
   xtmp = comdist * delta_ra * !pi/180./3600.
   ytmp = comdist * delta_dec* !pi/180./3600.

   ; GCIRC only gives absolute differences, so...
   if ra_tmp LT ra0 then xtmp = -xtmp
   if dec_tmp LT dec0 THEN ytmp = -ytmp

   f_vec = exp(-taueff_evo(zvec))
   delta_tmp = flux[forestcut]/cont_tmp[forestcut]/f_vec - 1
   wtmp = f_vec^2 * cont_tmp[forestcut]^2 * ivar[forestcut]
   noise_tmp =  (wtmp NE 0)/(sqrt(wtmp) + (wtmp EQ 0))
   avgweight = avg(f_vec^2 * cont_tmp[forestcut]^2 * ivar[forestcut])

   zerow = where(wtmp EQ 0, nzero)
   if nzero GT 0 then remove, zerow, delta_tmp, noise_tmp, losvec
   
   npix_tmp = n_elements(noise_tmp)

   delta_out = [delta_out, delta_tmp]
   noise_out =[noise_out , noise_tmp]
   weightlist = [weightlist, avgweight]

   x_out = [x_out, replicate(xtmp, npix_tmp)]
   y_out = [y_out, replicate(ytmp, npix_tmp)]
   z_out = [z_out, losvec]

   xlist = [xlist, xtmp]
   ylist = [ylist, ytmp]
   z0list = [z0list, min(losvec)]
   z1list = [z1list, max(losvec)]

   i0list = [i0list, ctr]
   i1list = [i1list, ctr + npix_tmp-1]

   printf, 19,catnum_tmp, ztmp, gmag_tmp,  ra_tmp, dec_tmp, xtmp, ytmp, $
           min(losvec), max(losvec), ctr, ctr+npix_tmp-1, $
           format='(I05,2x, f7.4,2x,f5.2,3x, f7.3,2x,f8.5,3x, f8.3,2x,' + $
           'f8.3,2(2x,f8.3), 2(2x,i6))' 
  
   ctr += npix_tmp

endfor
close, 19

;zerow = where(noise_out EQ 0.)
;remove, zerow, x_out, y_out, z_out, noise_out, delta_out

   ;; set a noise floor
   noise_out = noise_out > 0.2

print, n_elements(x_out), ' pixels'

pixel_arr = double([[x_out], [y_out], [z_out], [noise_out], [delta_out]])
pixel_arr = transpose(pixel_arr)

openw, 11, 'pixel_data.bin'
writeu, 11, pixel_arr
close, 11

stop


end
