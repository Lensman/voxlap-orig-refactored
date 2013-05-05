#include <math.h>

#include "../../include/sysmain.h"
#include "../../include/porthacks.h"
#include "../../include/voxlap5.h"
#include "../../include/kglobals.h"
#include "../../include/cpu_detect.h"
#include "../../include/ksnippits.h"

#define MAXCLIPIT (VSID*4) //VSID*2+4 is not a power of 2!
static lpoint2d clipit[MAXCLIPIT];

double findmaxcr (double px, double py, double pz, double cr)
{
	double f, g, maxcr, thresh2;
	long x, y, z, i0, i1, ix, y0, y1, z0, z1;
	char *v;

	thresh2 = cr+1.7321+1; thresh2 *= thresh2;
	maxcr = cr*cr;

		//Find closest point of all nearby cubes to (px,py,pz)
	x = (long)px; y = (long)py; z = (long)pz; i0 = i1 = 0; ix = x; y0 = y1 = y;
	while (1)
	{
		f = max(fabs((double)x+.5-px)-.5,0);
		g = max(fabs((double)y+.5-py)-.5,0);
		f = f*f + g*g;
		if (f < maxcr)
		{
			if (((unsigned long)x >= VSID) || ((unsigned long)y >= VSID))
				{ z0 = z1 = 0; }
			else
			{
				v = sptr[y*VSID+x];
				if (z >= v[1])
				{
					while (1)
					{
						if (!v[0]) { z0 = z1 = 0; break; }
						v += v[0]*4;
						if (z < v[1]) { z0 = v[3]; z1 = v[1]; break; }
					}
				}
				else { z0 = MAXZDIM-2048; z1 = v[1]; }
			}

			if ((pz <= z0) || (pz >= z1))
				maxcr = f;
			else
			{
				g = min(pz-(double)z0,(double)z1-pz);
				f += g*g; if (f < maxcr) maxcr = f;
			}
		}

		if ((x-px)*(x-px)+(y-py)*(y-py) < thresh2)
		{
			if ((x <= ix) && (x > 0))
				{ clipit[i1].x = x-1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
			if ((x >= ix) && (x < VSID-1))
				{ clipit[i1].x = x+1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
			if ((y <= y0) && (y > 0))
				{ clipit[i1].x = x; clipit[i1].y = y-1; i1 = ((i1+1)&(MAXCLIPIT-1)); y0--; }
			if ((y >= y1) && (y < VSID-1))
				{ clipit[i1].x = x; clipit[i1].y = y+1; i1 = ((i1+1)&(MAXCLIPIT-1)); y1++; }
		}
		if (i0 == i1) break;
		x = clipit[i0].x; y = clipit[i0].y; i0 = ((i0+1)&(MAXCLIPIT-1));
	}
	return(sqrt(maxcr));
}

#if 0

	//Point: (x,y), line segment: (px,py)-(px+vx,py+vy)
	//Returns 1 if point is closer than sqrt(cr2) to line
long dist2linept2d (double x, double y, double px, double py, double vx, double vy, double cr2)
{
	double f, g;
	x -= px; y -= py; f = x*vx + y*vy; if (f <= 0) return(x*x + y*y <= cr2);
	g = vx*vx + vy*vy; if (f >= g) { x -= vx; y -= vy; return(x*x + y*y <= cr2); }
	x = x*g-vx*f; y = y*g-vy*f; return(x*x + y*y <= cr2*g*g);
}

static char clipbuf[MAXZDIM+16]; //(8 extra on each side)
long sphtraceo (double px, double py, double pz,    //start pt
					double vx, double vy, double vz,    //move vector
					double *nx, double *ny, double *nz, //new pt after collision
					double *fx, double *fy, double *fz, //pt that caused collision
					double cr, double acr)
{
	double t, u, ex, ey, ez, Za, Zb, Zc, thresh2;
	double vxyz, vyz, vxz, vxy, rvxyz, rvyz, rvxz, rvxy, rvx, rvy, rvz, cr2;
	long i, i0, i1, x, y, z, xx, yy, zz, v, vv, ix, y0, y1, z0, z1;
	char *vp;

	t = 1;
	(*nx) = px + vx;
	(*ny) = py + vy;
	(*nz) = pz + vz;

	z0 = max((long)(min(pz,*nz)-cr)-2,-1);
	z1 = min((long)(max(pz,*nz)+cr)+2,MAXZDIM);

	thresh2 = cr+1.7321+1; thresh2 *= thresh2;

	vyz = vz*vz; vxz = vx*vx; vxy = vy*vy; cr2 = cr*cr;
	vyz += vxy; vxy += vxz; vxyz = vyz + vxz; vxz += vz*vz;
	rvx = 1.0 / vx; rvy = 1.0 / vy; rvz = 1.0 / vz;
	rvyz = 1.0 / vyz; rvxz = 1.0 / vxz; rvxy = 1.0 / vxy;
	rvxyz = 1.0 / vxyz;

		//Algorithm fails (stops short) if cr < 2 :(
	i0 = i1 = 0; ix = x = (long)px; y = y0 = y1 = (long)py;
	while (1)
	{
		for(z=z0;z<=z1;z++) clipbuf[z+8] = 0;
		i = 16;
		for(yy=y;yy<y+2;yy++)
			for(xx=x;xx<x+2;xx++,i<<=1)
			{
				z = z0;
				if ((unsigned long)(xx|yy) < VSID)
				{
					vp = sptr[yy*VSID+xx];
					while (1)
					{
						if (vp[1] > z) z = vp[1];
						if (!vp[0]) break;
						vp += vp[0]*4;
						zz = vp[3]; if (zz > z1) zz = z1;
						while (z < zz) clipbuf[(z++)+8] |= i;
					}
				}
				while (z <= z1) clipbuf[(z++)+8] |= i;
			}

		xx = x+1; yy = y+1; v = clipbuf[z0+8];
		for(z=z0;z<z1;z++)
		{
			zz = z+1; v = (v>>4)|clipbuf[zz+8];
			if ((!v) || (v == 255)) continue;

//---------------Check 1(8) corners of cube (sphere intersection)-------------

			//if (((v-1)^v) >= v)  //True if v is: {1,2,4,8,16,32,64,128}
			if (!(v&(v-1)))      //Same as above, but {0,1,2,4,...} (v's never 0)
			{
				ex = xx-px; ey = yy-py; ez = zz-pz;
				Zb = ex*vx + ey*vy + ez*vz;
				Zc = ex*ex + ey*ey + ez*ez - cr2;
				u = Zb*Zb - vxyz*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
						//   //Proposed compare optimization:
						//f = Zb*Zb-u; g = vxyz*t; h = (Zb*2-g)*g;
						//if ((unsigned __int64 *)&f < (unsigned __int64 *)&h)
					u = (Zb - sqrt(u)) * rvxyz;
					if ((u >= 0) && (u < t))
					{
						*fx = xx; *fy = yy; *fz = zz; t = u;
						*nx = vx*u + px; *ny = vy*u + py; *nz = vz*u + pz;
					}
				}
			}

//---------------Check 3(12) edges of cube (cylinder intersection)-----------

			vv = v&0x55; if (((vv-1)^vv) >= vv)  //True if (v&0x55)={1,4,16,64}
			{
				ey = yy-py; ez = zz-pz;
				Zb = ey*vy + ez*vz;
				Zc = ey*ey + ez*ez - cr2;
				u = Zb*Zb - vyz*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
					u = (Zb - sqrt(u)) * rvyz;
					if ((u >= 0) && (u < t))
					{
						ex = vx*u + px;
						if ((ex >= x) && (ex <= xx))
						{
							*fx = ex; *fy = yy; *fz = zz; t = u;
							*nx = ex; *ny = vy*u + py; *nz = vz*u + pz;
						}
					}
				}
			}
			vv = v&0x33; if (((vv-1)^vv) >= vv) //True if (v&0x33)={1,2,16,32}
			{
				ex = xx-px; ez = zz-pz;
				Zb = ex*vx + ez*vz;
				Zc = ex*ex + ez*ez - cr2;
				u = Zb*Zb - vxz*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
					u = (Zb - sqrt(u)) * rvxz;
					if ((u >= 0) && (u < t))
					{
						ey = vy*u + py;
						if ((ey >= y) && (ey <= yy))
						{
							*fx = xx; *fy = ey; *fz = zz; t = u;
							*nx = vx*u + px; *ny = ey; *nz = vz*u + pz;
						}
					}
				}
			}
			vv = v&0x0f; if (((vv-1)^vv) >= vv) //True if (v&0x0f)={1,2,4,8}
			{
				ex = xx-px; ey = yy-py;
				Zb = ex*vx + ey*vy;
				Zc = ex*ex + ey*ey - cr2;
				u = Zb*Zb - vxy*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
					u = (Zb - sqrt(u)) * rvxy;
					if ((u >= 0) && (u < t))
					{
						ez = vz*u + pz;
						if ((ez >= z) && (ez <= zz))
						{
							*fx = xx; *fy = yy; *fz = ez; t = u;
							*nx = vx*u + px; *ny = vy*u + py; *nz = ez;
						}
					}
				}
			}

//---------------Check 3(6) faces of cube (plane intersection)---------------

			if (vx)
			{
				switch(v&0x03)
				{
					case 0x01: ex = xx+cr; if ((vx > 0) || (px < ex)) goto skipfacex; break;
					case 0x02: ex = xx-cr; if ((vx < 0) || (px > ex)) goto skipfacex; break;
					default: goto skipfacex;
				}
				u = (ex - px) * rvx;
				if ((u >= 0) && (u < t))
				{
					ey = vy*u + py;
					ez = vz*u + pz;
					if ((ey >= y) && (ey <= yy) && (ez >= z) && (ez <= zz))
					{
						*fx = xx; *fy = ey; *fz = ez; t = u;
						*nx = ex; *ny = ey; *nz = ez;
					}
				}
			}
skipfacex:;
			if (vy)
			{
				switch(v&0x05)
				{
					case 0x01: ey = yy+cr; if ((vy > 0) || (py < ey)) goto skipfacey; break;
					case 0x04: ey = yy-cr; if ((vy < 0) || (py > ey)) goto skipfacey; break;
					default: goto skipfacey;
				}
				u = (ey - py) * rvy;
				if ((u >= 0) && (u < t))
				{
					ex = vx*u + px;
					ez = vz*u + pz;
					if ((ex >= x) && (ex <= xx) && (ez >= z) && (ez <= zz))
					{
						*fx = ex; *fy = yy; *fz = ez; t = u;
						*nx = ex; *ny = ey; *nz = ez;
					}
				}
			}
skipfacey:;
			if (vz)
			{
				switch(v&0x11)
				{
					case 0x01: ez = zz+cr; if ((vz > 0) || (pz < ez)) goto skipfacez; break;
					case 0x10: ez = zz-cr; if ((vz < 0) || (pz > ez)) goto skipfacez; break;
					default: goto skipfacez;
				}
				u = (ez - pz) * rvz;
				if ((u >= 0) && (u < t))
				{
					ex = vx*u + px;
					ey = vy*u + py;
					if ((ex >= x) && (ex <= xx) && (ey >= y) && (ey <= yy))
					{
						*fx = ex; *fy = ey; *fz = zz; t = u;
						*nx = ex; *ny = ey; *nz = ez;
					}
				}
			}
skipfacez:;
		}

		if ((x <= ix) && (x > 0) && (dist2linept2d(x-1,y,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x-1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((x >= ix) && (x < VSID-1) && (dist2linept2d(x+1,y,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x+1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((y <= y0) && (y > 0) && (dist2linept2d(x,y-1,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x; clipit[i1].y = y-1; i1 = ((i1+1)&(MAXCLIPIT-1)); y0--; }
		if ((y >= y1) && (y < VSID-1) && (dist2linept2d(x,y+1,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x; clipit[i1].y = y+1; i1 = ((i1+1)&(MAXCLIPIT-1)); y1++; }
		if (i0 == i1) break;
		x = clipit[i0].x; y = clipit[i0].y; i0 = ((i0+1)&(MAXCLIPIT-1));
	}

	if ((*nx) < acr) (*nx) = acr;
	if ((*ny) < acr) (*ny) = acr;
	if ((*nx) > VSID-acr) (*nx) = VSID-acr;
	if ((*ny) > VSID-acr) (*ny) = VSID-acr;
	if ((*nz) > MAXZDIM-1-acr) (*nz) = MAXZDIM-1-acr;
	if ((*nz) < MAXZDIM-2048) (*nz) = MAXZDIM-2048;

	return (t == 1);
}

#endif

static double gx0, gy0, gcrf2, grdst, gendt, gux, guy;
static long gdist2square (double x, double y)
{
	double t;
	x -= gx0; y -= gy0; t = x*gux + y*guy; if (t <= 0) t = gcrf2;
	else if (t*grdst >= gendt) { x -= gux*gendt; y -= guy*gendt; t = gcrf2; }
	else t = t*t*grdst + gcrf2;
	return(x*x + y*y <= t);
}

long sphtrace (double x0, double y0, double z0,          //start pt
					double vx, double vy, double vz,          //move vector
					double *hitx, double *hity, double *hitz, //new pt after collision
					double *clpx, double *clpy, double *clpz, //pt causing collision
					double cr, double acr)
{
	double f, t, dax, day, daz, vyx, vxy, vxz, vyz, rvz, cr2, fz, fc;
	double dx, dy, dx1, dy1;
	double nx, ny, intx, inty, intz, dxy, dxz, dyz, dxyz, rxy, rxz, ryz, rxyz;
	long i, j, x, y, ix, iy0, iy1, i0, i1, iz[2], cz0, cz1;
	char *v;

		 //Precalculate global constants for ins & getval functions
	if ((vx == 0) && (vy == 0) && (vz == 0))
		{ (*hitx) = x0; (*hity) = y0; (*hitz) = z0; return(1); }
	gux = vx; guy = vy; gx0 = x0; gy0 = y0; dxy = vx*vx + vy*vy;
	if (dxy != 0) rxy = 1.0 / dxy; else rxy = 0;
	grdst = rxy; gendt = 1; cr2 = cr*cr; t = cr + 0.7072; gcrf2 = t*t;

	if (((long *)&vz)[1] >= 0) { dtol(   z0-cr-.5,&cz0); dtol(vz+z0+cr-.5,&cz1); }
								 else { dtol(vz+z0-cr-.5,&cz0); dtol(   z0+cr-.5,&cz1); }

		//Precalculate stuff for closest point on cube finder
	dax = 0; day = 0; vyx = 0; vxy = 0; rvz = 0; vxz = 0; vyz = 0;
	if (vx != 0) { vyx = vy/vx; if (((long *)&vx)[1] >= 0) dax = x0+cr; else dax = x0-cr-1; }
	if (vy != 0) { vxy = vx/vy; if (((long *)&vy)[1] >= 0) day = y0+cr; else day = y0-cr-1; }
	if (vz != 0)
	{
		rvz = 1.0/vz; vxz = vx*rvz; vyz = vy*rvz;
		if (((long *)&vz)[1] >= 0) daz = z0+cr; else daz = z0-cr;
	}

	dxyz = vz*vz;
	dxz = vx*vx+dxyz; if (dxz != 0) rxz = 1.0 / dxz;
	dyz = vy*vy+dxyz; if (dyz != 0) ryz = 1.0 / dyz;
	dxyz += dxy; rxyz = 1.0 / dxyz;

	dtol(x0-.5,&x); dtol(y0-.5,&y);
	ix = x; iy0 = iy1 = y;
	i0 = 0; clipit[0].x = x; clipit[0].y = y; i1 = 1;
	do
	{
		x = clipit[i0].x; y = clipit[i0].y; i0 = ((i0+1)&(MAXCLIPIT-1));

		dx = (double)x; dx1 = (double)(x+1);
		dy = (double)y; dy1 = (double)(y+1);

			//closest point on cube finder
			//Plane intersection (both vertical planes)
#if 0
		intx = dbound((dy-day)*vxy + x0,dx,dx1);
		inty = dbound((dx-dax)*vyx + y0,dy,dy1);
#else
		intx = (dy-day)*vxy + x0;
		inty = (dx-dax)*vyx + y0;
		if (((long *)&intx)[1] < ((long *)&dx)[1]) intx = dx;
		if (((long *)&inty)[1] < ((long *)&dy)[1]) inty = dy;
		if (((long *)&intx)[1] >= ((long *)&dx1)[1]) intx = dx1;
		if (((long *)&inty)[1] >= ((long *)&dy1)[1]) inty = dy1;
		//if (intx < (double)x) intx = (double)x;
		//if (inty < (double)y) inty = (double)y;
		//if (intx > (double)(x+1)) intx = (double)(x+1);
		//if (inty > (double)(y+1)) inty = (double)(y+1);
#endif

		do
		{
			if (((long *)&dxy)[1] == 0) { t = -1.0; continue; }
			nx = intx-x0; ny = inty-y0; t = vx*nx + vy*ny; if (((long *)&t)[1] < 0) continue;
			f = cr2 - nx*nx - ny*ny; if (((long *)&f)[1] >= 0) { t = -1.0; continue; }
			f = f*dxy + t*t; if (((long *)&f)[1] < 0) { t = -1.0; continue; }
			t = (t-sqrt(f))*rxy;
		} while (0);
		if (t >= gendt) goto sphtracecont;
		if (((long *)&t)[1] < 0) intz = z0; else intz = vz*t + z0;

			//Find closest ceil(iz[0]) & flor(iz[1]) in (x,y) column
		dtol(intz-.5,&i);
		if ((unsigned long)(x|y) < VSID)
		{
			v = sptr[y*VSID+x]; iz[0] = MAXZDIM-2048; iz[1] = v[1];
			while (i >= iz[1])
			{
				if (!v[0]) { iz[1] = -1; break; }
				v += v[0]*4;
				iz[0] = v[3]; if (i < iz[0]) { iz[1] = -1; break; }
				iz[1] = v[1];
			}
		}
		else iz[1] = -1;

			//hit xz plane, yz plane or z-axis edge?
		if (iz[1] < 0) //Treat whole column as solid
		{
			if (((long *)&t)[1] >= 0) { gendt = t; (*clpx) = intx; (*clpy) = inty; (*clpz) = intz; goto sphtracecont; }
		}

			//Must check tops & bottoms of slab
		for(i=1;i>=0;i--)
		{
				//Ceil/flor outside of quick&dirty bounding box
			if ((iz[i] < cz0) || (iz[i] > cz1)) continue;

				//Plane intersection (parallel to ground)
			intz = (double)iz[i]; t = intz-daz;
			intx = t*vxz + x0;
			inty = t*vyz + y0;

			j = 0;                         // A ³ 8 ³ 9
			//     if (intx < dx)  j |= 2; //ÄÄÄÅÄÄÄÅÄÄÄ
			//else if (intx > dx1) j |= 1; // 2 ³ 0 ³ 1
			//     if (inty < dy)  j |= 8; //ÄÄÄÅÄÄÄÅÄÄÄ
			//else if (inty > dy1) j |= 4; // 6 ³ 4 ³ 5
				  if (((long *)&intx)[1] <  ((long *)&dx)[1])  j |= 2;
			else if (((long *)&intx)[1] >= ((long *)&dx1)[1]) j |= 1;
				  if (((long *)&inty)[1] <  ((long *)&dy)[1])  j |= 8;
			else if (((long *)&inty)[1] >= ((long *)&dy1)[1]) j |= 4;

				//NOTE: only need to check once per "for"!
			if ((!j) && (vz != 0)) //hit xy plane?
			{
				t *= rvz;
				if ((((long *)&t)[1] >= 0) && (t < gendt)) { gendt = t; (*clpx) = intx; (*clpy) = inty; (*clpz) = intz; }
				continue;
			}

				//common calculations used for rest of checks...
			fz = intz-z0; fc = cr2-fz*fz; fz *= vz;

			if (j&3)
			{
				nx = (double)((j&1)+x);
				if (((long *)&dxz)[1] != 0) //hit y-axis edge?
				{
					f = nx-x0; t = vx*f + fz; f = (fc - f*f)*dxz + t*t;
					if (((long *)&f)[1] >= 0) t = (t-sqrt(f))*rxz; else t = -1.0;
				} else t = -1.0;
				ny = vy*t + y0;
					  if (((long *)&ny)[1] > ((long *)&dy1)[1]) j |= 0x10;
				else if (((long *)&ny)[1] >= ((long *)&dy)[1])
				{
					if ((((long *)&t)[1] >= 0) && (t < gendt)) { gendt = t; (*clpx) = nx; (*clpy) = ny; (*clpz) = intz; }
					continue;
				}
				inty = (double)(((j>>4)&1)+y);
			}
			else inty = (double)(((j>>2)&1)+y);

			if (j&12)
			{
				ny = (double)(((j>>2)&1)+y);
				if (((long *)&dyz)[1] != 0) //hit x-axis edge?
				{
					f = ny-y0; t = vy*f + fz; f = (fc - f*f)*dyz + t*t;
					if (((long *)&f)[1] >= 0) t = (t-sqrt(f))*ryz; else t = -1.0;
				} else t = -1.0;
				nx = vx*t + x0;
					  if (((long *)&nx)[1] > ((long *)&dx1)[1]) j |= 0x20;
				else if (((long *)&nx)[1] >= ((long *)&dx)[1])
				{
					if ((((long *)&t)[1] >= 0) && (t < gendt)) { gendt = t; (*clpx) = nx; (*clpy) = ny; (*clpz) = intz; }
					continue;
				}
				intx = (double)(((j>>5)&1)+x);
			}
			else intx = (double)((j&1)+x);

				//hit corner?
			nx = intx-x0; ny = inty-y0;
			t = vx*nx + vy*ny + fz; if (((long *)&t)[1] < 0) continue;
			f = fc - nx*nx - ny*ny; if (((long *)&f)[1] >= 0) continue;
			f = f*dxyz + t*t; if (((long *)&f)[1] < 0) continue;
			t = (t-sqrt(f))*rxyz;
			if (t < gendt) { gendt = t; (*clpx) = intx; (*clpy) = inty; (*clpz) = intz; }
		}
sphtracecont:;
		if ((x <= ix)  && (x >      0) && (gdist2square(dx- .5,dy+ .5))) { clipit[i1].x = x-1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((x >= ix)  && (x < VSID-1) && (gdist2square(dx+1.5,dy+ .5))) { clipit[i1].x = x+1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((y <= iy0) && (y >      0) && (gdist2square(dx+ .5,dy- .5))) { clipit[i1].x = x; clipit[i1].y = y-1; i1 = ((i1+1)&(MAXCLIPIT-1)); iy0 = y-1; }
		if ((y >= iy1) && (y < VSID-1) && (gdist2square(dx+ .5,dy+1.5))) { clipit[i1].x = x; clipit[i1].y = y+1; i1 = ((i1+1)&(MAXCLIPIT-1)); iy1 = y+1; }
	} while (i0 != i1);
#if 1
	(*hitx) = dbound(vx*gendt + x0,acr,VSID-acr);
	(*hity) = dbound(vy*gendt + y0,acr,VSID-acr);
	(*hitz) = dbound(vz*gendt + z0,MAXZDIM-2048+acr,MAXZDIM-1-acr);
#else
	(*hitx) = min(max(vx*gendt + x0,acr),VSID-acr);
	(*hity) = min(max(vy*gendt + y0,acr),VSID-acr);
	(*hitz) = min(max(vz*gendt + z0,MAXZDIM-2048+acr),MAXZDIM-1-acr);
#endif
	return(gendt == 1);
}

void clipmove (dpoint3d *p, dpoint3d *v, double acr)
{
	double f, gx, gy, gz, nx, ny, nz, ex, ey, ez, hitx, hity, hitz, cr;
	//double nx2, ny2, nz2, ex2, ey2, ez2; //double ox, oy, oz;
	long i, j, k;

	//ox = p->x; oy = p->y; oz = p->z;
	gx = p->x+v->x; gy = p->y+v->y; gz = p->z+v->z;

	cr = findmaxcr(p->x,p->y,p->z,acr);
	vx5.clipmaxcr = cr;

	vx5.cliphitnum = 0;
	for(i=0;i<3;i++)
	{
		if ((v->x == 0) && (v->y == 0) && (v->z == 0)) break;

		cr -= 1e-7;  //Shrinking radius error control hack

		//j = sphtraceo(p->x,p->y,p->z,v->x,v->y,v->z,&nx,&ny,&nz,&ex,&ey,&ez,cr,acr);
		//k = sphtraceo(p->x,p->y,p->z,v->x,v->y,v->z,&nx2,&ny2,&nz2,&ex2,&ey2,&ez2,cr,acr);

		j = sphtrace(p->x,p->y,p->z,v->x,v->y,v->z,&nx,&ny,&nz,&ex,&ey,&ez,cr,acr);

		//if ((j != k) || (fabs(nx-nx2) > .000001) || (fabs(ny-ny2) > .000001) || (fabs(nz-nz2) > .000001) ||
		//   ((j == 0) && ((fabs(ex-ex2) > .000001) || (fabs(ey-ey2) > .000001) || (fabs(ez-ez2) > .000001))))
		//{
		//   printf("%d %f %f %f %f %f %f\n",i,p->x,p->y,p->z,v->x,v->y,v->z);
		//   printf("%f %f %f ",nx,ny,nz); if (!j) printf("%f %f %f\n",ex,ey,ez); else printf("\n");
		//   printf("%f %f %f ",nx2,ny2,nz2); if (!k) printf("%f %f %f\n",ex2,ey2,ez2); else printf("\n");
		//   printf("\n");
		//}
		if (j) { p->x = nx; p->y = ny; p->z = nz; break; }

		vx5.cliphit[i].x = ex; vx5.cliphit[i].y = ey; vx5.cliphit[i].z = ez;
		vx5.cliphitnum = i+1;
		p->x = nx; p->y = ny; p->z = nz;

			//Calculate slide vector
		v->x = gx-nx; v->y = gy-ny; v->z = gz-nz;
		switch(i)
		{
			case 0:
				hitx = ex-nx; hity = ey-ny; hitz = ez-nz;
				f = (v->x*hitx + v->y*hity + v->z*hitz) / (cr * cr);
				v->x -= hitx*f; v->y -= hity*f; v->z -= hitz*f;
				break;
			case 1:
				nx -= ex; ny -= ey; nz -= ez;
				ex = hitz*ny - hity*nz;
				ey = hitx*nz - hitz*nx;
				ez = hity*nx - hitx*ny;
				f = ex*ex + ey*ey + ez*ez; if (f <= 0) break;
				f = (v->x*ex + v->y*ey + v->z*ez) / f;
				v->x = ex*f; v->y = ey*f; v->z = ez*f;
				break;
			default: break;
		}
	}

		//If you didn't move much, then don't move at all. This helps prevents
		//cliprad from shrinking, but you get stuck too much :(
	//if ((p->x-ox)*(p->x-ox) + (p->y-oy)*(p->y-oy) + (p->z-oz)*(p->z-oz) < 1e-12)
	//   { p->x = ox; p->y = oy; p->z = oz; }
}

long cansee (point3d *p0, point3d *p1, lpoint3d *hit)
{
	lpoint3d a, c, d, p, i;
	point3d f, g;
	long cnt;

	ftol(p0->x-.5,&a.x); ftol(p0->y-.5,&a.y); ftol(p0->z-.5,&a.z);
	if (isvoxelsolid(a.x,a.y,a.z)) { hit->x = a.x; hit->y = a.y; hit->z = a.z; return(0); }
	ftol(p1->x-.5,&c.x); ftol(p1->y-.5,&c.y); ftol(p1->z-.5,&c.z);
	cnt = 0;

		  if (c.x <  a.x) { d.x = -1; f.x = p0->x-a.x;   g.x = (p0->x-p1->x)*1024; cnt += a.x-c.x; }
	else if (c.x != a.x) { d.x =  1; f.x = a.x+1-p0->x; g.x = (p1->x-p0->x)*1024; cnt += c.x-a.x; }
	else f.x = g.x = 0;
		  if (c.y <  a.y) { d.y = -1; f.y = p0->y-a.y;   g.y = (p0->y-p1->y)*1024; cnt += a.y-c.y; }
	else if (c.y != a.y) { d.y =  1; f.y = a.y+1-p0->y; g.y = (p1->y-p0->y)*1024; cnt += c.y-a.y; }
	else f.y = g.y = 0;
		  if (c.z <  a.z) { d.z = -1; f.z = p0->z-a.z;   g.z = (p0->z-p1->z)*1024; cnt += a.z-c.z; }
	else if (c.z != a.z) { d.z =  1; f.z = a.z+1-p0->z; g.z = (p1->z-p0->z)*1024; cnt += c.z-a.z; }
	else f.z = g.z = 0;

	ftol(f.x*g.z - f.z*g.x,&p.x); ftol(g.x,&i.x);
	ftol(f.y*g.z - f.z*g.y,&p.y); ftol(g.y,&i.y);
	ftol(f.y*g.x - f.x*g.y,&p.z); ftol(g.z,&i.z);

		//NOTE: GIGO! This can happen if p0,p1 (cansee input) is NaN, Inf, etc...
	if ((unsigned long)cnt > (VSID+VSID+2048)*2) cnt = (VSID+VSID+2048)*2;
	while (cnt > 0)
	{
		if (((p.x|p.y) >= 0) && (a.z != c.z)) { a.z += d.z; p.x -= i.x; p.y -= i.y; }
		else if ((p.z >= 0) && (a.x != c.x))  { a.x += d.x; p.x += i.z; p.z -= i.y; }
		else                                  { a.y += d.y; p.y += i.z; p.z += i.x; }
		if (isvoxelsolid(a.x,a.y,a.z)) break;
		cnt--;
	}
	hit->x = a.x; hit->y = a.y; hit->z = a.z; return(!cnt);
}

/**
 *  @param p0: start position
 *  @param v0: direction
 *  @param spr: pointer of sprite to test collision with
 *  @param h: coordinate of voxel hit in sprite coordinates (if any)
 *  @param ind: pointer to voxel hit (kv6voxtype) (0 if none hit)
 *  @param vsc:  input: max multiple/fraction of v0's length to scan (1.0 for |v0|)
 *  @return multiple/fraction of v0's length of hit point
 */
void sprhitscan (dpoint3d *p0, dpoint3d *v0, vx5sprite *spr, lpoint3d *h, kv6voxtype **ind, float *vsc)
{
	kv6voxtype *vx[4];
	kv6data *kv;
	point3d t, u, v;
	lpoint3d a, d, p, q;
	float f, g;
	long i, x, y, xup, ix0, ix1;

	(*ind) = 0;
	if (spr->flags&2)
	{
		kfatype *kf = spr->kfaptr;
			//This sets the sprite pointer to be the parent sprite (voxnum
			//   of the main sprite is invalid for KFA sprites!)
		spr = &kf->spr[kf->hingesort[(kf->numhin)-1]];
	}
	kv = spr->voxnum; if (!kv) return;

		//d transformed to spr space (0,0,0,kv->xsiz,kv->ysiz,kv->zsiz)
	v.x = v0->x*spr->s.x + v0->y*spr->s.y + v0->z*spr->s.z;
	v.y = v0->x*spr->h.x + v0->y*spr->h.y + v0->z*spr->h.z;
	v.z = v0->x*spr->f.x + v0->y*spr->f.y + v0->z*spr->f.z;

		//p transformed to spr space (0,0,0,kv->xsiz,kv->ysiz,kv->zsiz)
	t.x = p0->x-spr->p.x;
	t.y = p0->y-spr->p.y;
	t.z = p0->z-spr->p.z;
	u.x = t.x*spr->s.x + t.y*spr->s.y + t.z*spr->s.z;
	u.y = t.x*spr->h.x + t.y*spr->h.y + t.z*spr->h.z;
	u.z = t.x*spr->f.x + t.y*spr->f.y + t.z*spr->f.z;
	u.x /= (spr->s.x*spr->s.x + spr->s.y*spr->s.y + spr->s.z*spr->s.z);
	u.y /= (spr->h.x*spr->h.x + spr->h.y*spr->h.y + spr->h.z*spr->h.z);
	u.z /= (spr->f.x*spr->f.x + spr->f.y*spr->f.y + spr->f.z*spr->f.z);
	u.x += kv->xpiv; u.y += kv->ypiv; u.z += kv->zpiv;

	ix0 = max(vx5.xplanemin,0);
	ix1 = min(vx5.xplanemax,kv->xsiz);

		//Increment ray until it hits bounding box
		// (ix0,0,0,ix1-1ulp,kv->ysiz-1ulp,kv->zsiz-1ulp)
	g = (float)ix0;
	t.x = (float)ix1;      (*(long *)&t.x)--;
	t.y = (float)kv->ysiz; (*(long *)&t.y)--;
	t.z = (float)kv->zsiz; (*(long *)&t.z)--;
		  if (u.x <   g) { if (v.x <= 0) return; f = (  g-u.x)/v.x; u.x =   g; u.y += v.y*f; u.z += v.z*f; }
	else if (u.x > t.x) { if (v.x >= 0) return; f = (t.x-u.x)/v.x; u.x = t.x; u.y += v.y*f; u.z += v.z*f; }
		  if (u.y <   0) { if (v.y <= 0) return; f = (  0-u.y)/v.y; u.y =   0; u.x += v.x*f; u.z += v.z*f; }
	else if (u.y > t.y) { if (v.y >= 0) return; f = (t.y-u.y)/v.y; u.y = t.y; u.x += v.x*f; u.z += v.z*f; }
		  if (u.z <   0) { if (v.z <= 0) return; f = (  0-u.z)/v.z; u.z =   0; u.x += v.x*f; u.y += v.y*f; }
	else if (u.z > t.z) { if (v.z >= 0) return; f = (t.z-u.z)/v.z; u.z = t.z; u.x += v.x*f; u.y += v.y*f; }

	ix1 -= ix0;

	g = 262144.0 / sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

		//Note: (a.x,a.y,a.z) MUST be rounded towards -inf
	ftol(u.x-.5,&a.x); if ((unsigned long)(a.x-ix0) >= ix1) return;
	ftol(u.y-.5,&a.y); if ((unsigned long)a.y >= kv->ysiz) return;
	ftol(u.z-.5,&a.z); if ((unsigned long)a.z >= kv->zsiz) return;
	if (*(long *)&v.x < 0) { d.x = -1; u.x -= a.x;      v.x *= -g; }
							else { d.x =  1; u.x = a.x+1-u.x; v.x *=  g; }
	if (*(long *)&v.y < 0) { d.y = -1; u.y -= a.y;      v.y *= -g; }
							else { d.y =  1; u.y = a.y+1-u.y; v.y *=  g; }
	if (*(long *)&v.z < 0) { d.z = -1; u.z -= a.z;      v.z *= -g; }
							else { d.z =  1; u.z = a.z+1-u.z; v.z *=  g; }
	ftol(u.x*v.z - u.z*v.x,&p.x); ftol(v.x,&q.x);
	ftol(u.y*v.z - u.z*v.y,&p.y); ftol(v.y,&q.y);
	ftol(u.y*v.x - u.x*v.y,&p.z); ftol(v.z,&q.z);

		//Check if voxel at: (a.x,a.y,a.z) is solid
	vx[0] = kv->vox;
	for(x=0;x<a.x;x++) vx[0] += kv->xlen[x];
	vx[1] = vx[0]; xup = x*kv->ysiz;
	for(y=0;y<a.y;y++) vx[1] += kv->ylen[xup+y];
	vx[2] = vx[1]; vx[3] = &vx[1][kv->ylen[xup+y]];

	while (1)
	{
		//vs = kv->vox; //Brute force: remove all vx[?] code to enable this
		//for(x=0;x<a.x;x++) vs += kv->xlen[x];
		//for(y=0;y<a.y;y++) vs += kv->ylen[x*kv->ysiz+y];
		//for(ve=&vs[kv->ylen[x+y]];vs<ve;vs++) if (vs->z == a.z) break;

			//Check if voxel at: (a.x,a.y,a.z) is solid
		if (vx[1] < vx[3])
		{
			while ((a.z < vx[2]->z) && (vx[2] > vx[1]  )) vx[2]--;
			while ((a.z > vx[2]->z) && (vx[2] < vx[3]-1)) vx[2]++;
			if (a.z == vx[2]->z) break;
		}

		if ((p.x|p.y) >= 0)
		{
			a.z += d.z; if ((unsigned long)a.z >= kv->zsiz) return;
			p.x -= q.x; p.y -= q.y;
		}
		else if (p.z < 0)
		{
			a.y += d.y; if ((unsigned long)a.y >= kv->ysiz) return;
			p.y += q.z; p.z += q.x;

			if (a.y < y) { y--; vx[1] -= kv->ylen[xup+y];      }
			if (a.y > y) {      vx[1] += kv->ylen[xup+y]; y++; }
			vx[2] = vx[1]; vx[3] = &vx[1][kv->ylen[xup+y]];
		}
		else
		{
			a.x += d.x; if ((unsigned long)(a.x-ix0) >= ix1) return;
			p.x += q.z; p.z -= q.y;

			if (a.x < x) { x--; vx[0] -= kv->xlen[x];      xup -= kv->ysiz; }
			if (a.x > x) {      vx[0] += kv->xlen[x]; x++; xup += kv->ysiz; }
			if ((a.y<<1) < kv->ysiz) //Start y-slice search from closer side
			{
				vx[1] = vx[0];
				for(y=0;y<a.y;y++) vx[1] += kv->ylen[xup+y];
			}
			else
			{
				vx[1] = &vx[0][kv->xlen[x]];
				for(y=kv->ysiz;y>a.y;y--) vx[1] -= kv->ylen[xup+y-1];
			}
			vx[2] = vx[1]; vx[3] = &vx[1][kv->ylen[xup+y]];
		}
	}

		//given: a = kv6 coordinate, find: v = vxl coordinate
	u.x = (float)a.x-kv->xpiv;
	u.y = (float)a.y-kv->ypiv;
	u.z = (float)a.z-kv->zpiv;
	v.x = u.x*spr->s.x + u.y*spr->h.x + u.z*spr->f.x + spr->p.x;
	v.y = u.x*spr->s.y + u.y*spr->h.y + u.z*spr->f.y + spr->p.y;
	v.z = u.x*spr->s.z + u.y*spr->h.z + u.z*spr->f.z + spr->p.z;

		//Stupid dot product stuff...
	f = ((v.x-p0->x)*v0->x + (v.y-p0->y)*v0->y + (v.z-p0->z)*v0->z) /
		  (v0->x*v0->x + v0->y*v0->y + v0->z*v0->z);
	if (f >= (*vsc)) return;
	{ (*vsc) = f; (*h) = a; (*ind) = vx[2]; (*vsc) = f; }
}

unsigned long calcglobalmass ()
{
	unsigned long i, j;
	char *v;

	j = VSID*VSID*256;
	for(i=0;i<VSID*VSID;i++)
	{
		v = sptr[i]; j -= v[1];
		while (v[0]) { v += v[0]*4; j += v[3]-v[1]; }
	}
	return(j);
}

void orthonormalize (point3d *v0, point3d *v1, point3d *v2)
{
	float t;

	t = 1.0 / sqrt((v0->x)*(v0->x) + (v0->y)*(v0->y) + (v0->z)*(v0->z));
	(v0->x) *= t; (v0->y) *= t; (v0->z) *= t;
	t = (v1->x)*(v0->x) + (v1->y)*(v0->y) + (v1->z)*(v0->z);
	(v1->x) -= t*(v0->x); (v1->y) -= t*(v0->y); (v1->z) -= t*(v0->z);
	t = 1.0 / sqrt((v1->x)*(v1->x) + (v1->y)*(v1->y) + (v1->z)*(v1->z));
	(v1->x) *= t; (v1->y) *= t; (v1->z) *= t;
	(v2->x) = (v0->y)*(v1->z) - (v0->z)*(v1->y);
	(v2->y) = (v0->z)*(v1->x) - (v0->x)*(v1->z);
	(v2->z) = (v0->x)*(v1->y) - (v0->y)*(v1->x);
}

void dorthonormalize (dpoint3d *v0, dpoint3d *v1, dpoint3d *v2)
{
	double t;

	t = 1.0 / sqrt((v0->x)*(v0->x) + (v0->y)*(v0->y) + (v0->z)*(v0->z));
	(v0->x) *= t; (v0->y) *= t; (v0->z) *= t;
	t = (v1->x)*(v0->x) + (v1->y)*(v0->y) + (v1->z)*(v0->z);
	(v1->x) -= t*(v0->x); (v1->y) -= t*(v0->y); (v1->z) -= t*(v0->z);
	t = 1.0 / sqrt((v1->x)*(v1->x) + (v1->y)*(v1->y) + (v1->z)*(v1->z));
	(v1->x) *= t; (v1->y) *= t; (v1->z) *= t;
	(v2->x) = (v0->y)*(v1->z) - (v0->z)*(v1->y);
	(v2->y) = (v0->z)*(v1->x) - (v0->x)*(v1->z);
	(v2->z) = (v0->x)*(v1->y) - (v0->y)*(v1->x);
}

void orthorotate (float ox, float oy, float oz, point3d *ist, point3d *ihe, point3d *ifo)
{
	float f, t, dx, dy, dz, rr[9];

	fcossin(ox,&ox,&dx);
	fcossin(oy,&oy,&dy);
	fcossin(oz,&oz,&dz);
	f = ox*oz; t = dx*dz; rr[0] =  t*dy + f; rr[7] = -f*dy - t;
	f = ox*dz; t = dx*oz; rr[1] = -f*dy + t; rr[6] =  t*dy - f;
	rr[2] = dz*oy; rr[3] = -dx*oy; rr[4] = ox*oy; rr[8] = oz*oy; rr[5] = dy;
	ox = ist->x; oy = ihe->x; oz = ifo->x;
	ist->x = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->x = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->x = ox*rr[2] + oy*rr[5] + oz*rr[8];
	ox = ist->y; oy = ihe->y; oz = ifo->y;
	ist->y = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->y = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->y = ox*rr[2] + oy*rr[5] + oz*rr[8];
	ox = ist->z; oy = ihe->z; oz = ifo->z;
	ist->z = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->z = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->z = ox*rr[2] + oy*rr[5] + oz*rr[8];
	//orthonormalize(ist,ihe,ifo);
}

void dorthorotate (double ox, double oy, double oz, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	double f, t, dx, dy, dz, rr[9];

	dcossin(ox,&ox,&dx);
	dcossin(oy,&oy,&dy);
	dcossin(oz,&oz,&dz);
	f = ox*oz; t = dx*dz; rr[0] =  t*dy + f; rr[7] = -f*dy - t;
	f = ox*dz; t = dx*oz; rr[1] = -f*dy + t; rr[6] =  t*dy - f;
	rr[2] = dz*oy; rr[3] = -dx*oy; rr[4] = ox*oy; rr[8] = oz*oy; rr[5] = dy;
	ox = ist->x; oy = ihe->x; oz = ifo->x;
	ist->x = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->x = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->x = ox*rr[2] + oy*rr[5] + oz*rr[8];
	ox = ist->y; oy = ihe->y; oz = ifo->y;
	ist->y = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->y = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->y = ox*rr[2] + oy*rr[5] + oz*rr[8];
	ox = ist->z; oy = ihe->z; oz = ifo->z;
	ist->z = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->z = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->z = ox*rr[2] + oy*rr[5] + oz*rr[8];
	//dorthonormalize(ist,ihe,ifo);
}

void axisrotate (point3d *p, point3d *axis, float w)
{
	point3d ax;
	float t, c, s, ox, oy, oz, k[9];

	fcossin(w,&c,&s);
	t = axis->x*axis->x + axis->y*axis->y + axis->z*axis->z; if (t == 0) return;
	t = 1.0 / sqrt(t); ax.x = axis->x*t; ax.y = axis->y*t; ax.z = axis->z*t;

	t = 1.0-c;
	k[0] = ax.x*t; k[7] = ax.x*s; oz = ax.y*k[0];
	k[4] = ax.y*t; k[2] = ax.y*s; oy = ax.z*k[0];
	k[8] = ax.z*t; k[3] = ax.z*s; ox = ax.z*k[4];
	k[0] = ax.x*k[0] + c; k[5] = ox - k[7]; k[7] += ox;
	k[4] = ax.y*k[4] + c; k[6] = oy - k[2]; k[2] += oy;
	k[8] = ax.z*k[8] + c; k[1] = oz - k[3]; k[3] += oz;

	ox = p->x; oy = p->y; oz = p->z;
	p->x = ox*k[0] + oy*k[1] + oz*k[2];
	p->y = ox*k[3] + oy*k[4] + oz*k[5];
	p->z = ox*k[6] + oy*k[7] + oz*k[8];
}

void slerp (point3d *istr, point3d *ihei, point3d *ifor,
				point3d *istr2, point3d *ihei2, point3d *ifor2,
				point3d *ist, point3d *ihe, point3d *ifo, float rat)
{
	point3d ax;
	float c, s, t, ox, oy, oz, k[9];

	ist->x = istr->x; ist->y = istr->y; ist->z = istr->z;
	ihe->x = ihei->x; ihe->y = ihei->y; ihe->z = ihei->z;
	ifo->x = ifor->x; ifo->y = ifor->y; ifo->z = ifor->z;

	ax.x = istr->y*istr2->z - istr->z*istr2->y + ihei->y*ihei2->z - ihei->z*ihei2->y + ifor->y*ifor2->z - ifor->z*ifor2->y;
	ax.y = istr->z*istr2->x - istr->x*istr2->z + ihei->z*ihei2->x - ihei->x*ihei2->z + ifor->z*ifor2->x - ifor->x*ifor2->z;
	ax.z = istr->x*istr2->y - istr->y*istr2->x + ihei->x*ihei2->y - ihei->y*ihei2->x + ifor->x*ifor2->y - ifor->y*ifor2->x;
	t = ax.x*ax.x + ax.y*ax.y + ax.z*ax.z; if (t == 0) return;

		//Based on the vector suck-out method (see ROTATE2.BAS)
	ox = istr->x*ax.x + istr->y*ax.y + istr->z*ax.z;
	oy = ihei->x*ax.x + ihei->y*ax.y + ihei->z*ax.z;
	if (fabs(ox) < fabs(oy))
		{ c = istr->x*istr2->x + istr->y*istr2->y + istr->z*istr2->z; s = ox*ox; }
	else
		{ c = ihei->x*ihei2->x + ihei->y*ihei2->y + ihei->z*ihei2->z; s = oy*oy; }
	if (t == s) return;
	c = (c*t - s) / (t-s);
	if (c < -1) c = -1;
	if (c > 1) c = 1;
	fcossin(acos(c)*rat,&c,&s);

	t = 1.0 / sqrt(t); ax.x *= t; ax.y *= t; ax.z *= t;

	t = 1.0f-c;
	k[0] = ax.x*t; k[7] = ax.x*s; oz = ax.y*k[0];
	k[4] = ax.y*t; k[2] = ax.y*s; oy = ax.z*k[0];
	k[8] = ax.z*t; k[3] = ax.z*s; ox = ax.z*k[4];
	k[0] = ax.x*k[0] + c; k[5] = ox - k[7]; k[7] += ox;
	k[4] = ax.y*k[4] + c; k[6] = oy - k[2]; k[2] += oy;
	k[8] = ax.z*k[8] + c; k[1] = oz - k[3]; k[3] += oz;

	ox = ist->x; oy = ist->y; oz = ist->z;
	ist->x = ox*k[0] + oy*k[1] + oz*k[2];
	ist->y = ox*k[3] + oy*k[4] + oz*k[5];
	ist->z = ox*k[6] + oy*k[7] + oz*k[8];

	ox = ihe->x; oy = ihe->y; oz = ihe->z;
	ihe->x = ox*k[0] + oy*k[1] + oz*k[2];
	ihe->y = ox*k[3] + oy*k[4] + oz*k[5];
	ihe->z = ox*k[6] + oy*k[7] + oz*k[8];

	ox = ifo->x; oy = ifo->y; oz = ifo->z;
	ifo->x = ox*k[0] + oy*k[1] + oz*k[2];
	ifo->y = ox*k[3] + oy*k[4] + oz*k[5];
	ifo->z = ox*k[6] + oy*k[7] + oz*k[8];
}