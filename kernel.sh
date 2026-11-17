#!/bin/bash
#
# Perform the 2D Gaussian kernel estimator
#
# usage: ./kernel.sh colvars_final.txt [output]
#   input : x y w
#   output: x y F dF/dx dF/dy   (F in kcal/mol)

DATA=${1:?usage: $0 data.txt [out.txt]}
OUT=${2:-kernel.txt}
KERNEL=${KERNEL:-./kernel}

# settings
X0=-2.5 ; XM=1.5     # grid range in x (must enclose the data)
Y0=-1.0 ; YM=3.0     # grid range in y
N=200                # grid points per dimension (nx = ny = N)
R=1.0                # bandwidth scale
                     # Start at R = 1; Scan a few values (0.5,1,2)
		     # Pick the smallest R where those numbers stop changing 
		     # and spurious minima are gone; Keep h at least about 0.5-1 grid spacing
DH=10                # integer: kernel cut-off in units of h
TEMP=300             # temp

CAP=4.0              # cap for FE

# R<1: sharper wells and barriers, but noisier;
# R>1: smooth, but wells come out shallower and broader and barriers lower

# clean up (get rid of seperators)
awk '$1<999 && $2<999 && $3>0 {print $1,$2,$3}' "$DATA" > data.txt

# write kernel.in
awk -v x0=$X0 -v xm=$XM -v y0=$Y0 -v ym=$YM -v N=$N -v r=$R -v dh=$DH -v temp=$TEMP '
NR==1 {minx=maxx=$1; miny=maxy=$2}
{
    x=$1; y=$2; n++
    sx+=x; sy+=y; sx2+=x*x; sy2+=y*y
    if (x<minx) minx=x; if (x>maxx) maxx=x
    if (y<miny) miny=y; if (y>maxy) maxy=y
}
END {
    sdx=sqrt((n*sx2-sx*sx)/(n*(n-1)))
    sdy=sqrt((n*sy2-sy*sy)/(n*(n-1)))
    print x0; print y0; print xm; print ym
    print N;  print N
    print r*0.9*exp(-log(n)/5)*sdx      # h_x
    print r*0.9*exp(-log(n)/5)*sdy      # h_y
    print dh
    print temp
    print "comments:"
    print "data x range", minx, maxx, " y range", miny, maxy
    print "average", sx/n, sy/n
    print "s.d.", sdx, sdy
    print n
    if (minx<=x0 || maxx>=xm || miny<=y0 || maxy>=ym)
        print "WARNING: data extend beyond the grid; those points are ignored" > "/dev/stderr"
}' data.txt > kernel.in

# run the kernel program
awk -v n=$(tail -1 kernel.in) 'BEGIN{print n}{print}' data.txt | $KERNEL > "$OUT"

awk -v cap=$CAP 'NR==FNR{if(NR==1||$3<m)m=$3;next} {f=$3-m; if(f>cap)f=cap; print $1,$2,f}' kernel.txt kernel.txt > fe2d.txt
