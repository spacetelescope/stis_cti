
dir = '/home/lockwood/cti_testing/test_data/science/'

root = 'obr101010'
a = mrdfits(dir + 'orig/' + root + '_crj.fits.gz', 1)
b = mrdfits(dir + root + '_crc.fits', 1)

zmin = -30
zmax =  20

xr = [300, 800]
yr = [  0, 500]

psopen, root + '_comparison.eps', /encapsulated, xsize=8, ysize=4, /inches

!p.multi = [0,2,1]
charsize = 1  ; !d.name eq 'PS' ? 1.5 : 1
thick    = !d.name eq 'PS' ? 3 : 1


display, a[xr[0]:xr[1], yr[0]:yr[1]], $
    findgen(xr[1]-xr[0]+1)+xr[0], findgen(yr[1]-yr[0]+1)+yr[0], $
    min=zmin, max=zmax, title=root + '_crj', $
    xtitle='X', ytitle='Y', $
    xthick=thick, ythick=thick, charthick=thick, charsize=charsize
display, b[xr[0]:xr[1], yr[0]:yr[1]], min=zmin, max=zmax, title=root + '_crc', /noerase, $
    findgen(xr[1]-xr[0]+1)+xr[0], findgen(yr[1]-yr[0]+1)+yr[0], $
    xtitle='X', ytitle='Y', $
    xthick=thick, ythick=thick, charthick=thick, charsize=charsize


if !d.name eq 'PS' then psclose
END
