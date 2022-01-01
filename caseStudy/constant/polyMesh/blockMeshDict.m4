/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           | Copyright (C) 2016 Ehsan Madadi-Kandjani        |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General macros to create cylinder mesh

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(btQuad, ($1b $2b $2t $1t))
define(topQuad, ($1t $4t $3t $2t))
define(bottomQuad, ($1b $2b $3b $4b))

//Mathematical constants:
define(pi, 3.1415926536)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

// Base z
define(Zb, 0)

// shaft radious
define(rs, 49.9)

// bearing radius
define(rb, 50)

// Height of cylinder
define(zh, 100)

// Outlet z
define(Zt, calc(Zb + zh))

// Number of cells at shaft
define(Ns, 45)

// Number of cells between shaft and bearing
define(Nb, 5)

// Number of cells in the bearing height
define(Nz, 10)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// esessentric x and y magnitude 

// x 
define(kx, 7.8e-2)

// y
define(ky, 0)

// bearing x and y position

// x
define(xb0, calc(rb*cos( (pi/180)*0 ) ) ) 
define(xb1, calc(rb*cos( (pi/180)*90 ) ) ) 
define(xb2, calc(rb*cos( (pi/180)*180 ) ) ) 
define(xb3, calc(rb*cos( (pi/180)*270 ) ) ) 

// y
define(yb0, calc(rb*sin( (pi/180)*0 ) ) ) 
define(yb1, calc(rb*sin( (pi/180)*90 ) ) ) 
define(yb2, calc(rb*sin( (pi/180)*180 ) ) )
define(yb3, calc(rb*sin( (pi/180)*270 ) ) )

//  x and y positions beetween arcs bearing

// x
define(bx01b, calc(rb*cos((pi/180)*45) ) )
define(bx12b, calc(rb*cos((pi/180)*135) ) )
define(bx23b, calc(rb*cos((pi/180)*225) ) )
define(bx30b, calc(rb*cos((pi/180)*315) ) )

// y
define(by01b, calc(rb*sin((pi/180)*45) ) )
define(by12b, calc(rb*sin((pi/180)*135) ) )
define(by23b, calc(rb*sin((pi/180)*225) ) )
define(by30b, calc(rb*sin((pi/180)*315) ) )

// shaft x and y position

// x
define(xs0, calc( rs*cos((pi/180)*0) + kx ) )
define(xs1, calc( rs*cos((pi/180)*90) + kx ) )
define(xs2, calc( rs*cos((pi/180)*180) + kx ) )
define(xs3, calc( rs*cos((pi/180)*270) + kx ) )

// y
define(ys0, calc( rs*sin((pi/180)*0) + ky ) )
define(ys1, calc( rs*sin((pi/180)*90) + ky ) )
define(ys2, calc( rs*sin((pi/180)*180) + ky ) )
define(ys3, calc( rs*sin((pi/180)*270) + ky ) )

//  x and y positions beetween arcs shaft

// x
define(bx01s, calc( rs*cos((pi/180)*45) + kx ) )
define(bx12s, calc( rs*cos((pi/180)*135) + kx ) )
define(bx23s, calc( rs*cos((pi/180)*225) + kx ) )
define(bx30s, calc( rs*cos((pi/180)*315) + kx ) ) 

// y
define(by01s, calc( rs*sin((pi/180)*45) + ky ) ) 
define(by12s, calc( rs*sin((pi/180)*135) + ky ) )
define(by23s, calc( rs*sin((pi/180)*225) + ky ) )
define(by30s, calc( rs*sin((pi/180)*315) + ky ) ) 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    (xb0 yb0 Zb) vlabel(b0b)
    (xb1 yb1 Zb) vlabel(b1b)
    (xb2 yb2 Zb) vlabel(b2b)
    (xb3 yb3 Zb) vlabel(b3b)
    
    (xs0 ys0 Zb) vlabel(s0b)
    (xs1 ys1 Zb) vlabel(s1b)
    (xs2 ys2 Zb) vlabel(s2b)
    (xs3 ys3 Zb) vlabel(s3b)
    
    (xb0 yb0 Zt) vlabel(b0t)
    (xb1 yb1 Zt) vlabel(b1t)
    (xb2 yb2 Zt) vlabel(b2t)
    (xb3 yb3 Zt) vlabel(b3t)
    
    (xs0 ys0 Zt) vlabel(s0t)
    (xs1 ys1 Zt) vlabel(s1t)
    (xs2 ys2 Zt) vlabel(s2t)
    (xs3 ys3 Zt) vlabel(s3t)
);

blocks
(
    //block0
    hex2D(s0, b0, b1, s1)
    (Nb Ns Nz)
    simpleGrading (1 1 1)
    
    //block1
    hex2D(s1, b1, b2, s2)
    (Nb Ns Nz)
    simpleGrading (1 1 1)
    
    //block2
    hex2D(s2, b2, b3, s3)
    (Nb Ns Nz)
    simpleGrading (1 1 1)
    
    //block3
    hex2D(s3, b3, b0, s0)
    (Nb Ns Nz)
    simpleGrading (1 1 1)   
    

);

edges
(
    //bearing edges
    arc b0b b1b (bx01b by01b Zb)
    arc b1b b2b (bx12b by12b Zb)
    arc b2b b3b (bx23b by23b Zb)
    arc b3b b0b (bx30b by30b Zb)

    arc b0t b1t (bx01b by01b Zt)
    arc b1t b2t (bx12b by12b Zt)
    arc b2t b3t (bx23b by23b Zt)
    arc b3t b0t (bx30b by30b Zt)
    
    //shaft edges
    arc s0b s1b (bx01s by01s Zb)
    arc s1b s2b (bx12s by12s Zb)
    arc s2b s3b (bx23s by23s Zb)
    arc s3b s0b (bx30s by30s Zb)

    arc s0t s1t (bx01s by01s Zt)
    arc s1t s2t (bx12s by12s Zt)
    arc s2t s3t (bx23s by23s Zt)
    arc s3t s0t (bx30s by30s Zt)      
    
);

patches
(
    wall walls
    (
        btQuad(b0, b1)
        btQuad(b1, b2)
        btQuad(b2, b3)
        btQuad(b3, b0)
    )

    wall movingwalls
    (
        btQuad(s0, s1)
        btQuad(s1, s2)
        btQuad(s2, s3)
        btQuad(s3, s0)
     )
    
    patch inlet
    (
        bottomQuad(b0, b1, s1, s0)
        bottomQuad(b1, b2, s2, s1)
        bottomQuad(b2, b3, s3, s2)
        bottomQuad(b3, b0, s0, s3)
    )
    
    patch outlet
    (
        topQuad(b0, b1, s1, s0)
        topQuad(b1, b2, s2, s1)
        topQuad(b2, b3, s3, s2)
        topQuad(b3, b0, s0, s3)
    )
);

mergePatchPairs
(
);
