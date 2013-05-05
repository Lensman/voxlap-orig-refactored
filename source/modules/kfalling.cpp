#include <math.h>
#include "../../include/sysmain.h"
#include "../../include/porthacks.h"
#include "../../include/voxlap5.h"
#include "../../include/kglobals.h"
#include "../../include/kmodelling.h"

extern long *vbuf, *vbit;
extern long bbuf[], p2c[], p2m[]; 

	//float detection & falling code variables...
	//WARNING: VLSTSIZ,FSTKSIZ,FLCHKSIZ can all have bounds errors! :(
#define VLSTSIZ 65536 //Theoretically should be at least: VOXSIZ\8
#define LOGHASHEAD 12
#define FSTKSIZ 8192

extern double logint[];      //was static, now in kmodelling
extern float tempfloatbuf[]; //was static, now in kmodelling
extern long factr[][2];      //was static, now in kmodelling

typedef struct { long v, b; } vlstyp;
vlstyp vlst[VLSTSIZ];
long hhead[1<<LOGHASHEAD], vlstcnt = 0x7fffffff;
lpoint3d fstk[FSTKSIZ]; //Note .z is actually used as a pointer, not z!
#define FLCHKSIZ 4096
lpoint3d flchk[FLCHKSIZ]; long flchkcnt = 0;

//float detection & falling code begins --------------------------------------
//How to use this section of code:
//Step 1: Call checkfloatinbox after every "deleting" set* call
//Step 2: Call dofalls(); at a constant rate in movement code

	//Adds all slabs inside box (inclusive) to "float check" list
void checkfloatinbox (long x0, long y0, long z0, long x1, long y1, long z1)
{
	long x, y;
	char *ov, *v;

	if (flchkcnt >= FLCHKSIZ) return;

		//Make all off-by-1 hacks in other code unnecessary
	x0 = max(x0-1,0); x1 = min(x1+1,VSID);
	y0 = max(y0-1,0); y1 = min(y1+1,VSID);
	z0 = max(z0-1,0); z1 = min(z1+1,MAXZDIM);

		//Add local box's slabs to flchk list - checked in next dofalls()
	for(y=y0;y<y1;y++)
		for(x=x0;x<x1;x++)
		{
			v = sptr[y*VSID+x];
			while (1)
			{
				ov = v; if ((z1 <= v[1]) || (!v[0])) break;
				v += v[0]*4; if (z0 >= v[3]) continue;
				flchk[flchkcnt].x = x;
				flchk[flchkcnt].y = y;
				flchk[flchkcnt].z = ov[1];
				flchkcnt++; if (flchkcnt >= FLCHKSIZ) return;
			}
		}
}

void isnewfloatingadd (const long f)
{
	long v = (((f>>(LOGHASHEAD+3))-(f>>3)) & ((1<<LOGHASHEAD)-1));
	vlst[vlstcnt].b = hhead[v]; hhead[v] = vlstcnt;
	vlst[vlstcnt].v = f; vlstcnt++;
}

long isnewfloatingot (const long f)
{
	long v = hhead[((f>>(LOGHASHEAD+3))-(f>>3)) & ((1<<LOGHASHEAD)-1)];
	while (1)
	{
		if (v < 0) return(-1);
		if (vlst[v].v == f) return(v);
		v = vlst[v].b;
	}
}

	//removes a & adds b while preserving index; used only by meltfall(...)
	//Must do nothing if 'a' not in hash

void isnewfloatingchg (long a, long b)
{
	long ov, v, i, j;

	i = (((a>>(LOGHASHEAD+3))-(a>>3)) & ((1<<LOGHASHEAD)-1));
	j = (((b>>(LOGHASHEAD+3))-(b>>3)) & ((1<<LOGHASHEAD)-1));

	v = hhead[i]; ov = -1;
	while (v >= 0)
	{
		if (vlst[v].v == a)
		{
			vlst[v].v = b; if (i == j) return;
			if (ov < 0) hhead[i] = vlst[v].b; else vlst[ov].b = vlst[v].b;
			vlst[v].b = hhead[j]; hhead[j] = v;
			return;
		}
		ov = v; v = vlst[v].b;
	}
}

long isnewfloating (flstboxtype *flb)
{
	float f;
	lpoint3d p, cen;
	long i, j, nx, ny, z0, z1, fend, ovlstcnt, mass;
	char *v, *ov;

	p.x = flb->chk.x; p.y = flb->chk.y; p.z = flb->chk.z;
	v = sptr[p.y*VSID+p.x];
	while (1)
	{
		ov = v;
		if ((p.z < v[1]) || (!v[0])) return(0);
		v += v[0]*4;
		if (p.z < v[3]) break;
	}

	if (isnewfloatingot((long)ov) >= 0) return(0);
	ovlstcnt = vlstcnt;
	isnewfloatingadd((long)ov);
	if (vlstcnt >= VLSTSIZ) return(0); //EVIL HACK TO PREVENT CRASH!

		//Init: centroid, mass, bounding box
	cen = p; mass = 1;
	flb->x0 = p.x-1; flb->y0 = p.y-1; flb->z0 = p.z-1;
	flb->x1 = p.x+1; flb->y1 = p.y+1; flb->z1 = p.z+1;

	fend = 0;
	while (1)
	{
		z0 = ov[1];         if (z0 < flb->z0) flb->z0 = z0;
		z1 = ov[ov[0]*4+3]; if (z1 > flb->z1) flb->z1 = z1;

		i = z1-z0;
		cen.x += p.x*i;
		cen.y += p.y*i;
		cen.z += (((z0+z1)*i)>>1); //sum(z0 to z1-1)
		mass += i;

		for(i=0;i<8;i++) //26-connectivity
		{
			switch(i)
			{
				case 0: nx = p.x-1; ny = p.y  ; if (nx < flb->x0) flb->x0 = nx; break;
				case 1: nx = p.x  ; ny = p.y-1; if (ny < flb->y0) flb->y0 = ny; break;
				case 2: nx = p.x+1; ny = p.y  ; if (nx > flb->x1) flb->x1 = nx; break;
				case 3: nx = p.x  ; ny = p.y+1; if (ny > flb->y1) flb->y1 = ny; break;
				case 4: nx = p.x-1; ny = p.y-1; break;
				case 5: nx = p.x+1; ny = p.y-1; break;
				case 6: nx = p.x-1; ny = p.y+1; break;
				case 7: nx = p.x+1; ny = p.y+1; break;
				default: _gtfo(); //tells MSVC default can't be reached
			}
			if ((unsigned long)(nx|ny) >= VSID) continue;

			v = sptr[ny*VSID+nx];
			while (1)
			{
				if (!v[0])
				{
					if (v[1] <= z1) return(0);  //This MUST be <=, (not <) !!!
					break;
				}
				ov = v; v += v[0]*4; //NOTE: this is a 'different' ov
				if ((ov[1] > z1) || (z0 > v[3])) continue; //26-connectivity
				j = isnewfloatingot((long)ov);
				if (j < 0)
				{
					isnewfloatingadd((long)ov);
					if (vlstcnt >= VLSTSIZ) return(0); //EVIL HACK TO PREVENT CRASH!
					fstk[fend].x = nx; fstk[fend].y = ny; fstk[fend].z = (long)ov;
					fend++; if (fend >= FSTKSIZ) return(0); //EVIL HACK TO PREVENT CRASH!
					continue;
				}
				if ((unsigned long)j < ovlstcnt) return(0);
			}
		}

		if (!fend)
		{
			flb->i0 = ovlstcnt;
			flb->i1 = vlstcnt;
			flb->mass = mass; f = 1.0 / (float)mass;
			flb->centroid.x = (float)cen.x*f;
			flb->centroid.y = (float)cen.y*f;
			flb->centroid.z = (float)cen.z*f;
			return(1);
		}
		fend--;
		p.x = fstk[fend].x; p.y = fstk[fend].y; ov = (char *)fstk[fend].z;
	}
}

void startfalls ()
{
	long i, z;

		//This allows clear to be MUCH faster when there isn't much falling
	if (vlstcnt < ((1<<LOGHASHEAD)>>1))
	{
		for(i=vlstcnt-1;i>=0;i--)
		{
			z = vlst[i].v;
			hhead[((z>>(LOGHASHEAD+3))-(z>>3)) & ((1<<LOGHASHEAD)-1)] = -1;
		}
	}
	else { for(z=0;z<(1<<LOGHASHEAD);z++) hhead[z] = -1; }

		//Float detection...
		//flstcnt[].i0/i1 tell which parts of vlst are floating
	vlstcnt = 0;

		//Remove any current pieces that are no longer floating
	for(i=vx5.flstnum-1;i>=0;i--)
		if (!isnewfloating(&vx5.flstcnt[i])) //Modifies flstcnt,vlst[],vlstcnt
			vx5.flstcnt[i] = vx5.flstcnt[--vx5.flstnum]; //onground, so delete flstcnt[i]

		//Add new floating pieces (while space is left on flstcnt)
	if (vx5.flstnum < FLPIECES)
		for(i=flchkcnt-1;i>=0;i--)
		{
			vx5.flstcnt[vx5.flstnum].chk = flchk[i];
			if (isnewfloating(&vx5.flstcnt[vx5.flstnum])) //Modifies flstcnt,vlst[],vlstcnt
			{
				vx5.flstcnt[vx5.flstnum].userval = -1; //New piece: let game programmer know
				vx5.flstnum++; if (vx5.flstnum >= FLPIECES) break;
			}
		}
	flchkcnt = 0;
}

	//Call 0 or 1 times (per flstcnt) between startfalls&finishfalls
void dofall (long i)
{
	long j, z;
	char *v;

		//Falling code... call this function once per piece
	vx5.flstcnt[i].chk.z++;
	for(z=vx5.flstcnt[i].i1-1;z>=vx5.flstcnt[i].i0;z--)
	{
		v = (char *)vlst[z].v; v[1]++; v[2]++;
		v = &v[v[0]*4];
		v[3]++;
		if ((v[3] == v[1]) && (vx5.flstcnt[i].i1 >= 0))
		{
			j = isnewfloatingot((long)v);
				//Make sure it's not part of the same floating object
			if ((j < vx5.flstcnt[i].i0) || (j >= vx5.flstcnt[i].i1))
				vx5.flstcnt[i].i1 = -1; //Mark flstcnt[i] for scum2 fixup
		}
	}

	if (vx5.vxlmipuse > 1)
	{
		long x0, y0, x1, y1;
		x0 = max(vx5.flstcnt[i].x0,0); x1 = min(vx5.flstcnt[i].x1+1,VSID);
		y0 = max(vx5.flstcnt[i].y0,0); y1 = min(vx5.flstcnt[i].y1+1,VSID);
		//FIX ME!!!
		//if ((x1 > x0) && (y1 > y0)) genmipvxl(x0,y0,x1,y1); //Don't replace with bbox!
	}
}

	//Sprite structure is already allocated
	//kv6, vox, xlen, ylen are all malloced in here!
long meltfall (vx5sprite *spr, const long fi, const long delvxl)
{
	long i, j, k, x, y, z, xs, ys, zs, xe, ye, ze;
	long oxvoxs, oyvoxs, numvoxs;
	char *v, *ov, *nv;
	kv6data *kv;
	kv6voxtype *voxptr;
	unsigned long *xlenptr;
	unsigned short *ylenptr;

	if (vx5.flstcnt[fi].i1 < 0) return(0);

	xs = max(vx5.flstcnt[fi].x0,0); xe = min(vx5.flstcnt[fi].x1,VSID-1);
	ys = max(vx5.flstcnt[fi].y0,0); ye = min(vx5.flstcnt[fi].y1,VSID-1);
	zs = max(vx5.flstcnt[fi].z0,0); ze = min(vx5.flstcnt[fi].z1,MAXZDIM-1);
	if ((xs > xe) || (ys > ye) || (zs > ze)) return(0);

		//Need to know how many voxels to allocate... SLOW :(
	numvoxs = vx5.flstcnt[fi].i0-vx5.flstcnt[fi].i1;
	for(i=vx5.flstcnt[fi].i0;i<vx5.flstcnt[fi].i1;i++)
		numvoxs += ((char *)vlst[i].v)[0];
	if (numvoxs <= 0) return(0); //No voxels found!

	spr->p = vx5.flstcnt[fi].centroid;
	spr->s.x = 1.f; spr->h.x = 0.f; spr->f.x = 0.f;
	spr->s.y = 0.f; spr->h.y = 1.f; spr->f.y = 0.f;
	spr->s.z = 0.f; spr->h.z = 0.f; spr->f.z = 1.f;

	x = xe-xs+1; y = ye-ys+1; z = ze-zs+1;

	j = sizeof(kv6data) + numvoxs*sizeof(kv6voxtype) + x*4 + x*y*2;
	i = (long)malloc(j); if (!i) return(0); if (i&3) { free((void *)i); return(0); }
	spr->voxnum = kv = (kv6data *)i; spr->flags = 0;
	kv->leng = j;
	kv->xsiz = x;
	kv->ysiz = y;
	kv->zsiz = z;
	kv->xpiv = spr->p.x - xs;
	kv->ypiv = spr->p.y - ys;
	kv->zpiv = spr->p.z - zs;
	kv->numvoxs = numvoxs;
	kv->namoff = 0;
	kv->lowermip = 0;
	kv->vox = (kv6voxtype *)((long)spr->voxnum+sizeof(kv6data));
	kv->xlen = (unsigned long *)(((long)kv->vox)+numvoxs*sizeof(kv6voxtype));
	kv->ylen = (unsigned short *)(((long)kv->xlen) + kv->xsiz*4);

	voxptr = kv->vox; numvoxs = 0;
	xlenptr = kv->xlen; oxvoxs = 0;
	ylenptr = kv->ylen; oyvoxs = 0;

	for(x=xs;x<=xe;x++)
	{
		for(y=ys;y<=ye;y++)
		{
			for(v=sptr[y*VSID+x];v[0];v=nv)
			{
				nv = v+v[0]*4;

				i = isnewfloatingot((long)v);
				if (((unsigned long)i >= vx5.flstcnt[fi].i1) || (i < vx5.flstcnt[fi].i0))
					continue;

				for(z=v[1];z<=v[2];z++)
				{
					voxptr[numvoxs].col = lightvox(*(long *)&v[((z-v[1])<<2)+4]);
					voxptr[numvoxs].z = z-zs;

					voxptr[numvoxs].vis = 0; //OPTIMIZE THIS!!!
					if (!isvoxelsolid(x-1,y,z)) voxptr[numvoxs].vis |= 1;
					if (!isvoxelsolid(x+1,y,z)) voxptr[numvoxs].vis |= 2;
					if (!isvoxelsolid(x,y-1,z)) voxptr[numvoxs].vis |= 4;
					if (!isvoxelsolid(x,y+1,z)) voxptr[numvoxs].vis |= 8;
					//if (z == v[1]) voxptr[numvoxs].vis |= 16;
					//if (z == nv[3]-1) voxptr[numvoxs].vis |= 32;
					if (!isvoxelsolid(x,y,z-1)) voxptr[numvoxs].vis |= 16;
					if (!isvoxelsolid(x,y,z+1)) voxptr[numvoxs].vis |= 32;

					voxptr[numvoxs].dir = 0; //FIX THIS!!!
					numvoxs++;
				}
				for(z=nv[3]+v[2]-v[1]-v[0]+2;z<nv[3];z++)
				{
					voxptr[numvoxs].col = lightvox(*(long *)&nv[(z-nv[3])<<2]);
					voxptr[numvoxs].z = z-zs;

					voxptr[numvoxs].vis = 0; //OPTIMIZE THIS!!!
					if (!isvoxelsolid(x-1,y,z)) voxptr[numvoxs].vis |= 1;
					if (!isvoxelsolid(x+1,y,z)) voxptr[numvoxs].vis |= 2;
					if (!isvoxelsolid(x,y-1,z)) voxptr[numvoxs].vis |= 4;
					if (!isvoxelsolid(x,y+1,z)) voxptr[numvoxs].vis |= 8;
					//if (z == v[1]) voxptr[numvoxs].vis |= 16;
					//if (z == nv[3]-1) voxptr[numvoxs].vis |= 32;
					if (!isvoxelsolid(x,y,z-1)) voxptr[numvoxs].vis |= 16;
					if (!isvoxelsolid(x,y,z+1)) voxptr[numvoxs].vis |= 32;

					voxptr[numvoxs].dir = 0; //FIX THIS!!!
					numvoxs++;
				}

#if 0
				if (delvxl) //Quick&dirty dealloc from VXL (bad for holes!)
				{
						//invalidate current vptr safely
					isnewfloatingchg((long)v,0);

					k = nv-v; //perform slng(nv) and adjust vlst at same time
					for(ov=nv;ov[0];ov+=ov[0]*4)
						isnewfloatingchg((long)ov,((long)ov)-k);

					j = (long)ov-(long)nv+(ov[2]-ov[1]+1)*4+4;

						//shift end of RLE column up
					v[0] = nv[0]; v[1] = nv[1]; v[2] = nv[2];
					for(i=4;i<j;i+=4) *(long *)&v[i] = *(long *)&nv[i];

						//remove end of RLE column from vbit
					i = ((((long)(&v[i]))-(long)vbuf)>>2); j = (k>>2)+i;
#if 0
					while (i < j) { vbit[i>>5] &= ~(1<<i); i++; }
#else
					if (!((j^i)&~31))
						vbit[i>>5] &= ~(p2m[j&31]^p2m[i&31]);
					else
					{
						vbit[i>>5] &=   p2m[i&31];  i >>= 5;
						vbit[j>>5] &= (~p2m[j&31]); j >>= 5;
						for(j--;j>i;j--) vbit[j] = 0;
					}
#endif
					nv = v;
				}
#endif
			}
			*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs;
		}
		*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
	}

	if (delvxl)
		for(x=xs;x<=xe;x++)
			for(y=ys;y<=ye;y++)
				for(v=sptr[y*VSID+x];v[0];v=nv)
				{
					nv = v+v[0]*4;

					i = isnewfloatingot((long)v);
					if (((unsigned long)i >= vx5.flstcnt[fi].i1) || (i < vx5.flstcnt[fi].i0))
						continue;

						//Quick&dirty dealloc from VXL (bad for holes!)

						//invalidate current vptr safely
					isnewfloatingchg((long)v,0);

					k = nv-v; //perform slng(nv) and adjust vlst at same time
					for(ov=nv;ov[0];ov+=ov[0]*4)
						isnewfloatingchg((long)ov,((long)ov)-k);

					j = (long)ov-(long)nv+(ov[2]-ov[1]+1)*4+4;

						//shift end of RLE column up
					v[0] = nv[0]; v[1] = nv[1]; v[2] = nv[2];
					for(i=4;i<j;i+=4) *(long *)&v[i] = *(long *)&nv[i];

						//remove end of RLE column from vbit
					i = ((((long)(&v[i]))-(long)vbuf)>>2); j = (k>>2)+i;
#if 0
					while (i < j) { vbit[i>>5] &= ~(1<<i); i++; }
#else
					if (!((j^i)&~31))
						vbit[i>>5] &= ~(p2m[j&31]^p2m[i&31]);
					else
					{
						vbit[i>>5] &=   p2m[i&31];  i >>= 5;
						vbit[j>>5] &= (~p2m[j&31]); j >>= 5;
						for(j--;j>i;j--) vbit[j] = 0;
					}
#endif
					nv = v;
				}

	vx5.flstcnt[fi].i1 = -2; //Mark flstcnt[i] invalid; no scum2 fixup

	if (vx5.vxlmipuse > 1) genmipvxl(xs,ys,xe+1,ye+1);

	return(vx5.flstcnt[fi].mass);
}

void finishfalls ()
{
	long i, x, y;

		//Scum2 box fixup: refreshes rle voxel data inside a bounding rectangle
	for(i=vx5.flstnum-1;i>=0;i--)
		if (vx5.flstcnt[i].i1 < 0)
		{
			if (vx5.flstcnt[i].i1 == -1)
			{
				for(y=vx5.flstcnt[i].y0;y<=vx5.flstcnt[i].y1;y++)
					for(x=vx5.flstcnt[i].x0;x<=vx5.flstcnt[i].x1;x++)
						scum2(x,y);
				scum2finish();
				updatebbox(vx5.flstcnt[i].x0,vx5.flstcnt[i].y0,vx5.flstcnt[i].z0,vx5.flstcnt[i].x1,vx5.flstcnt[i].y1,vx5.flstcnt[i].z1,0);
			}
			vx5.flstcnt[i] = vx5.flstcnt[--vx5.flstnum]; //onground, so delete flstcnt[i]
		}
}

/**
 * This converts a spherical cut-out of the VXL map into a .KV6 sprite in
 * memory. This function can be used to make walls fall over (with full
 * rotation). It allocates a new vx5sprite sprite structure and you are
 * responsible for freeing the memory using "free" in your own code.
 *
 * @param spr new vx5sprite structure. Position & orientation are initialized
 *        so when you call drawsprite, it exactly matches the VXL map.
 * @param hit center of sphere
 * @param hitrad radius of sphere
 * @param returns 0:bad, >0:mass of captured object (# of voxels)
 */
long meltsphere (vx5sprite *spr, lpoint3d *hit, long hitrad)
{
	long i, j, x, y, z, xs, ys, zs, xe, ye, ze, sq, z0, z1;
	long oxvoxs, oyvoxs, numvoxs, cx, cy, cz, cw;
	float f, ff;
	kv6data *kv;
	kv6voxtype *voxptr;
	unsigned long *xlenptr;
	unsigned short *ylenptr;

	xs = max(hit->x-hitrad,0); xe = min(hit->x+hitrad,VSID-1);
	ys = max(hit->y-hitrad,0); ye = min(hit->y+hitrad,VSID-1);
	zs = max(hit->z-hitrad,0); ze = min(hit->z+hitrad,MAXZDIM-1);
	if ((xs > xe) || (ys > ye) || (zs > ze)) return(0);

	if (hitrad >= SETSPHMAXRAD-1) hitrad = SETSPHMAXRAD-2;

	tempfloatbuf[0] = 0.0f;
#if 0
		//Totally unoptimized
	for(i=1;i<=hitrad;i++) tempfloatbuf[i] = pow((float)i,vx5.curpow);
#else
	tempfloatbuf[1] = 1.0f;
	for(i=2;i<=hitrad;i++)
	{
		if (!factr[i][0]) tempfloatbuf[i] = exp(logint[i]*vx5.curpow);
		else tempfloatbuf[i] = tempfloatbuf[factr[i][0]]*tempfloatbuf[factr[i][1]];
	}
#endif
	*(long *)&tempfloatbuf[hitrad+1] = 0x7f7fffff; //3.4028235e38f; //Highest float

// ---------- Need to know how many voxels to allocate... SLOW!!! :( ----------
	cx = cy = cz = 0; //Centroid
	cw = 0;       //Weight (1 unit / voxel)
	numvoxs = 0;
	sq = 0; //pow(fabs(x-hit->x),vx5.curpow) + "y + "z < pow(vx5.currad,vx5.curpow)
	for(x=xs;x<=xe;x++)
	{
		ff = tempfloatbuf[hitrad]-tempfloatbuf[labs(x-hit->x)];
		for(y=ys;y<=ye;y++)
		{
			f = ff-tempfloatbuf[labs(y-hit->y)];
			if (*(long *)&f > 0) //WARNING: make sure to always write ylenptr!
			{
				while (*(long *)&tempfloatbuf[sq] <  *(long *)&f) sq++;
				while (*(long *)&tempfloatbuf[sq] >= *(long *)&f) sq--;
				z0 = max(hit->z-sq,zs); z1 = min(hit->z+sq+1,ze);
				for(z=z0;z<z1;z++)
				{
					i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
					if (i)
					{
						cx += (x-hit->x); cy += (y-hit->y); cz += (z-hit->z); cw++;
					}
					if ((i == 0) || ((i == 1) && (1))) continue; //not_on_border))) continue; //FIX THIS!!!
					numvoxs++;
				}
			}
		}
	}
	if (numvoxs <= 0) return(0); //No voxels found!
// ---------------------------------------------------------------------------

	f = 1.0 / (float)cw; //Make center of sprite the centroid
	spr->p.x = (float)hit->x + (float)cx*f;
	spr->p.y = (float)hit->y + (float)cy*f;
	spr->p.z = (float)hit->z + (float)cz*f;
	spr->s.x = 1.f; spr->h.x = 0.f; spr->f.x = 0.f;
	spr->s.y = 0.f; spr->h.y = 1.f; spr->f.y = 0.f;
	spr->s.z = 0.f; spr->h.z = 0.f; spr->f.z = 1.f;

	x = xe-xs+1; y = ye-ys+1; z = ze-zs+1;

	j = sizeof(kv6data) + numvoxs*sizeof(kv6voxtype) + x*4 + x*y*2;
	i = (long)malloc(j); if (!i) return(0); if (i&3) { free((void *)i); return(0); }
	spr->voxnum = kv = (kv6data *)i; spr->flags = 0;
	kv->leng = j;
	kv->xsiz = x;
	kv->ysiz = y;
	kv->zsiz = z;
	kv->xpiv = spr->p.x - xs;
	kv->ypiv = spr->p.y - ys;
	kv->zpiv = spr->p.z - zs;
	kv->numvoxs = numvoxs;
	kv->namoff = 0;
	kv->lowermip = 0;
	kv->vox = (kv6voxtype *)((long)spr->voxnum+sizeof(kv6data));
	kv->xlen = (unsigned long *)(((long)kv->vox)+numvoxs*sizeof(kv6voxtype));
	kv->ylen = (unsigned short *)(((long)kv->xlen) + kv->xsiz*4);

	voxptr = kv->vox; numvoxs = 0;
	xlenptr = kv->xlen; oxvoxs = 0;
	ylenptr = kv->ylen; oyvoxs = 0;

	sq = 0; //pow(fabs(x-hit->x),vx5.curpow) + "y + "z < pow(vx5.currad,vx5.curpow)
	for(x=xs;x<=xe;x++)
	{
		ff = tempfloatbuf[hitrad]-tempfloatbuf[labs(x-hit->x)];
		for(y=ys;y<=ye;y++)
		{
			f = ff-tempfloatbuf[labs(y-hit->y)];
			if (*(long *)&f > 0) //WARNING: make sure to always write ylenptr!
			{
				while (*(long *)&tempfloatbuf[sq] <  *(long *)&f) sq++;
				while (*(long *)&tempfloatbuf[sq] >= *(long *)&f) sq--;
				z0 = max(hit->z-sq,zs); z1 = min(hit->z+sq+1,ze);
				for(z=z0;z<z1;z++)
				{
					i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
					if ((i == 0) || ((i == 1) && (1))) continue; //not_on_border))) continue; //FIX THIS!!!
					voxptr[numvoxs].col = lightvox(*(long *)i);
					voxptr[numvoxs].z = z-zs;
					voxptr[numvoxs].vis = 63; //FIX THIS!!!
					voxptr[numvoxs].dir = 0; //FIX THIS!!!
					numvoxs++;
				}
			}
			*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs;
		}
		*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
	}
	return(cw);
}

/**
 * This function is similar to meltsphere, except you can use any user-
 * defined shape (with some size limits). The user-defined shape is
 * described by a list of vertical columns in the "vspans" format:
 *
 * typedef struct { char z1, z0, x, y; } vspans;
 *
 * The list MUST be ordered first in increasing Y, then in increasing X
 * or else the function will crash! Fortunately, the structure is
 * arranged in a way that the data can be sorted quite easily using a
 * simple trick: if you use a typecast from vspans to "unsigned long",
 * you can use a generic sort code on 32-bit integers to achieve a
 * correct sort. The vspans members are all treated as unsigned chars,
 * so it's usually a good idea to bias your columns by 128, and then
 * reverse-bias them in the "offs" offset.
 *
 * @param spr new vx5sprite structure. Position & orientation are initialized
 *            so when you call drawsprite, it exactly matches the VXL map.
 * @param lst list in "vspans" format
 * @param lstnum number of columns on list
 * @param offs offset of top-left corner in VXL coordinates
 * @return mass (in voxel units), returns 0 if error (or no voxels)
 */
long meltspans (vx5sprite *spr, vspans *lst, long lstnum, lpoint3d *offs)
{
	float f;
	long i, j, x, y, z, xs, ys, zs, xe, ye, ze, z0, z1;
	long ox, oy, oxvoxs, oyvoxs, numvoxs, cx, cy, cz, cw;
	kv6data *kv;
	kv6voxtype *voxptr;
	unsigned long *xlenptr;
	unsigned short *ylenptr;

	if (lstnum <= 0) return(0);
// ---------- Need to know how many voxels to allocate... SLOW!!! :( ----------
	cx = cy = cz = 0; //Centroid
	cw = 0;       //Weight (1 unit / voxel)
	numvoxs = 0;
	xs = xe = ((long)lst[0].x)+offs->x;
	ys = ((long)lst[       0].y)+offs->y;
	ye = ((long)lst[lstnum-1].y)+offs->y;
	zs = ze = ((long)lst[0].z0)+offs->z;
	for(j=0;j<lstnum;j++)
	{
		x = ((long)lst[j].x)+offs->x;
		y = ((long)lst[j].y)+offs->y; if ((x|y)&(~(VSID-1))) continue;
			  if (x < xs) xs = x;
		else if (x > xe) xe = x;
		z0 = ((long)lst[j].z0)+offs->z;   if (z0 < 0) z0 = 0;
		z1 = ((long)lst[j].z1)+offs->z+1; if (z1 > MAXZDIM) z1 = MAXZDIM;
		if (z0 < zs) zs = z0;
		if (z1 > ze) ze = z1;
		for(z=z0;z<z1;z++) //getcube too SLOW... FIX THIS!!!
		{
			i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
			if (i) { cx += x-offs->x; cy += y-offs->y; cz += z-offs->z; cw++; }
			if (i&~1) numvoxs++;
		}
	}
	if (numvoxs <= 0) return(0); //No voxels found!
// ---------------------------------------------------------------------------

	f = 1.0 / (float)cw; //Make center of sprite the centroid
	spr->p.x = (float)offs->x + (float)cx*f;
	spr->p.y = (float)offs->y + (float)cy*f;
	spr->p.z = (float)offs->z + (float)cz*f;
	spr->x.x = 0.f; spr->y.x = 1.f; spr->z.x = 0.f;
	spr->x.y = 1.f; spr->y.y = 0.f; spr->z.y = 0.f;
	spr->x.z = 0.f; spr->y.z = 0.f; spr->z.z = 1.f;

	x = xe-xs+1; y = ye-ys+1; z = ze-zs;

	j = sizeof(kv6data) + numvoxs*sizeof(kv6voxtype) + y*4 + x*y*2;
	i = (long)malloc(j); if (!i) return(0); if (i&3) { free((void *)i); return(0); }
	spr->voxnum = kv = (kv6data *)i; spr->flags = 0;
	kv->leng = j;
	kv->xsiz = y;
	kv->ysiz = x;
	kv->zsiz = z;
	kv->xpiv = spr->p.y - ys;
	kv->ypiv = spr->p.x - xs;
	kv->zpiv = spr->p.z - zs;
	kv->numvoxs = numvoxs;
	kv->namoff = 0;
	kv->lowermip = 0;
	kv->vox = (kv6voxtype *)((long)spr->voxnum+sizeof(kv6data));
	kv->xlen = (unsigned long *)(((long)kv->vox)+numvoxs*sizeof(kv6voxtype));
	kv->ylen = (unsigned short *)(((long)kv->xlen) + kv->xsiz*4);

	voxptr = kv->vox; numvoxs = 0;
	xlenptr = kv->xlen; oxvoxs = 0;
	ylenptr = kv->ylen; oyvoxs = 0;
	ox = xs; oy = ys;
	for(j=0;j<lstnum;j++)
	{
		x = ((long)lst[j].x)+offs->x;
		y = ((long)lst[j].y)+offs->y; if ((x|y)&(~(VSID-1))) continue;
		while ((ox != x) || (oy != y))
		{
			*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs; ox++;
			if (ox > xe)
			{
				*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
				ox = xs; oy++;
			}
		}
		z0 = ((long)lst[j].z0)+offs->z;   if (z0 < 0) z0 = 0;
		z1 = ((long)lst[j].z1)+offs->z+1; if (z1 > MAXZDIM) z1 = MAXZDIM;
		for(z=z0;z<z1;z++) //getcube TOO SLOW... FIX THIS!!!
		{
			i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
			if (!(i&~1)) continue;
			voxptr[numvoxs].col = lightvox(*(long *)i);
			voxptr[numvoxs].z = z-zs;

			voxptr[numvoxs].vis = 63; //FIX THIS!!!
			//if (!isvoxelsolid(x-1,y,z)) voxptr[numvoxs].vis |= 1;
			//if (!isvoxelsolid(x+1,y,z)) voxptr[numvoxs].vis |= 2;
			//if (!isvoxelsolid(x,y-1,z)) voxptr[numvoxs].vis |= 4;
			//if (!isvoxelsolid(x,y+1,z)) voxptr[numvoxs].vis |= 8;
			//if (!isvoxelsolid(x,y,z-1)) voxptr[numvoxs].vis |= 16;
			//if (!isvoxelsolid(x,y,z+1)) voxptr[numvoxs].vis |= 32;

			voxptr[numvoxs].dir = 0; //FIX THIS!!!
			numvoxs++;
		}
	}
	while (1)
	{
		*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs; ox++;
		if (ox > xe)
		{
			*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
			ox = xs; oy++; if (oy > ye) break;
		}
	}
	return(cw);
}
