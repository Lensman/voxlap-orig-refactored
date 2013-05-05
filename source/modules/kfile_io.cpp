#include <math.h>
#include <stdio.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <stdarg.h>
#include <string.h>
#include <conio.h>
#include <dos.h>
#define MAX_PATH 260
#endif
#include <stdlib.h>
#define USE_MEMMAPPED_IO = 1
#include "../../include/kglobals.h"
#include "../../include/kfile_io.h"
#include "../../include/kplib.h"
#include "../../include/ksnippits.h"


long totclk;
long numsprites;
char sxlfilnam[MAX_PATH+1] = "";
char vxlfilnam[MAX_PATH+1] = "";
char skyfilnam[MAX_PATH+1] = "";

	//Displayed text message:
char message[256] = {0};
long messagetimeout = 0;
long sxlparspos, sxlparslen;

extern long *vbuf, *vbit, vbiti;
extern char *khashbuf;
extern long khashead[256], khashpos, khashsiz;
extern vx5sprite spr[];

extern long numsprites;
extern char message[256];
extern long messagetimeout;
extern long totclk;

extern void evilquit (const char *);
extern unsigned long calcglobalmass ();
extern long colorjit (long, long);
extern void scum2finish();
extern inline void delslab (long *, const long, long);
extern inline long *scum2 (const long &, const long &);
extern inline long lbound0 (long, long);
extern inline long lbound (long, long, long);
extern long slng (const char *);

extern void floodsucksprite (vx5sprite *, kv6data *, long, long, kv6voxtype *, kv6voxtype *);
extern void kfasorthinge (hingetype *, long, long *);

#if defined(USE_MEMMAPPED_IO)
	mem_mapped_file_info vbuf_file_info, vbit_file_info;
	extern mem_mapped_file_info span_file_info;
#endif


/** Given a name associated with the kv6data/kfatype return a pointer
 * Notice that each structure has a "namoff" member. Since I
 * use remalloc(), I have to make these offsets, not pointers. Use this
 * function to convert the offsets into pointers.
 *
 * @param namoff offset to the name
 * @return a pointer to the filename associated with the kv6data/kfatype
 */
char *getkfilname ( const long &namoff) { return(&khashbuf[namoff]); }

/**
 * Loads and begins parsing of an .SXL file. Always call this first before
 * using parspr().
 * @param sxlnam .SXL filename
 * @param vxlnam pointer to .VXL filename (written by loadsxl)
 * @param vxlnam pointer to .SKY filename (written by loadsxl)
 * @param globst pointer to global user string. You parse this yourself!
 *               You can edit this in Voxed by pressing F6.
 * @return 0: loadsxl failed (file not found or malloc failed)
 *         1: loadsxl successful; call parspr()!
 */
long loadsxl (const char *sxlnam, char **vxlnam, char **skynam, char **globst)
{
	long j, k, m, n;

		//NOTE: MUST buffer file because insertsprite uses kz file code :/
	if (!kzopen(sxlnam)) return(0);
	sxlparslen = kzfilelength();
	if (sxlbuf) { free(sxlbuf); sxlbuf = 0; }
	if (!(sxlbuf = (char *)malloc(sxlparslen))) return(0);
	kzread(sxlbuf,sxlparslen);
	kzclose();

	j = n = 0;

		//parse vxlnam
	(*vxlnam) = &sxlbuf[j];
	while ((sxlbuf[j]!=13)&&(sxlbuf[j]!=10) && (j < sxlparslen)) j++; sxlbuf[j++] = 0;
	while (((sxlbuf[j]==13)||(sxlbuf[j]==10)) && (j < sxlparslen)) j++;

		//parse skynam
	(*skynam) = &sxlbuf[j];
	while ((sxlbuf[j]!=13)&&(sxlbuf[j]!=10) && (j < sxlparslen)) j++; sxlbuf[j++] = 0;
	while (((sxlbuf[j]==13)||(sxlbuf[j]==10)) && (j < sxlparslen)) j++;

		//parse globst
	m = n = j; (*globst) = &sxlbuf[n];
	while (((sxlbuf[j] == ' ') || (sxlbuf[j] == 9)) && (j < sxlparslen))
	{
		j++;
		while ((sxlbuf[j]!=13) && (sxlbuf[j]!=10) && (j < sxlparslen)) sxlbuf[n++] = sxlbuf[j++];
		sxlbuf[n++] = 13; j++;
		while (((sxlbuf[j]==13) || (sxlbuf[j]==10)) && (j < sxlparslen)) j++;
	}
	if (n > m) sxlbuf[n-1] = 0; else (*globst) = &nullst;

		//Count number of sprites in .SXL file (helpful for GAME)
	sxlparspos = j;
	return(1);
}

/**
 * If loadsxl returns a 1, then you should call parspr with a while loop
 * that terminates when the return value is 0.
 * @param spr pointer to sprite structure (written by parspr) You allocate
 *            the vx5sprite, & parspr fills in the position&orientation.
 * @param userst pointer to user string associated with the given sprite.
 *               You can edit this in Voxed by right-clicking the sprite.
 * @return pointer to .KV6 filename OR NULL if no more sprites left.
 *         You must load the .KV6 to memory yourself by doing:
 *         char *kv6filename = parspr(...)
 *         if (kv6filename) spr->voxnum = getkv6(kv6filename);
 */
char *parspr (vx5sprite *spr, char **userst)
{
	float f;
	long j, k, m, n;
	char *namptr;

	j = sxlparspos; //unnecessary temp variable (to shorten code)

		//Automatically free temp sxlbuf when done reading sprites
	if (((j+2 < sxlparslen) && (sxlbuf[j] == 'e') && (sxlbuf[j+1] == 'n') && (sxlbuf[j+2] == 'd') &&
		((j+3 == sxlparslen) || (sxlbuf[j+3] == 13) || (sxlbuf[j+3] == 10))) || (j > sxlparslen))
		return(0);

		//parse kv6name
	for(k=j;(sxlbuf[k]!=',') && (k < sxlparslen);k++); sxlbuf[k] = 0;
	namptr = &sxlbuf[j]; j = k+1;

		//parse 12 floats
	for(m=0;m<12;m++)
	{
		if (m < 11) { for(k=j;(sxlbuf[k]!=',') && (k < sxlparslen);k++); }
		else { for(k=j;(sxlbuf[k]!=13) && (sxlbuf[k]!=10) && (k < sxlparslen);k++); }

		sxlbuf[k] = 0; f = atof(&sxlbuf[j]); j = k+1;
		switch(m)
		{
			case  0: spr->p.x = f; break;
			case  1: spr->p.y = f; break;
			case  2: spr->p.z = f; break;
			case  3: spr->s.x = f; break;
			case  4: spr->s.y = f; break;
			case  5: spr->s.z = f; break;
			case  6: spr->h.x = f; break;
			case  7: spr->h.y = f; break;
			case  8: spr->h.z = f; break;
			case  9: spr->f.x = f; break;
			case 10: spr->f.y = f; break;
			case 11: spr->f.z = f; break;
			default: _gtfo(); //tells MSVC default can't be reached
		}
	}
	while (((sxlbuf[j]==13) || (sxlbuf[j]==10)) && (j < sxlparslen)) j++;

	spr->flags = 0;

		//parse userst
	m = n = j; (*userst) = &sxlbuf[n];
	while (((sxlbuf[j] == ' ') || (sxlbuf[j] == 9)) && (j < sxlparslen))
	{
		j++;
		while ((sxlbuf[j]!=13) && (sxlbuf[j]!=10) && (j < sxlparslen)) sxlbuf[n++] = sxlbuf[j++];
		sxlbuf[n++] = 13; j++;
		while (((sxlbuf[j]==13) || (sxlbuf[j]==10)) && (j < sxlparslen)) j++;
	}
	if (n > m) sxlbuf[n-1] = 0; else (*userst) = &nullst;

	sxlparspos = j; //unnecessary temp variable (for short code)
	return(namptr);
}
/**
 * Loads a sky into memory. Sky must be PNG,JPG,TGA,GIF,BMP,PCX formatted as
 * a Mercator projection on its side. This means x-coordinate is latitude
 * and y-coordinate is longitude. Loadsky() can be called at any time.
 *
 * If for some reason you don't want to load a textured sky, you call call
 * loadsky with these 2 built-in skies:
 * loadsky("BLACK");  //pitch black
 * loadsky("BLUE");   //a cool ramp of bright blue to blue to greenish
 *
 * @param skyfilnam the name of the image to load
 * @return -1:bad, 0:good
 */

long loadsky (const char *skyfilnam)
{
	long x, y, xoff, yoff;
	float ang, f;

	if (skypic) { free((void *)skypic); skypic = skyoff = 0; }
	xoff = yoff = 0;

	if (!strcasecmp(skyfilnam,"BLACK")) return(0);
	if (!strcasecmp(skyfilnam,"BLUE")) goto loadbluesky;

	kpzload(skyfilnam,(int*)&skypic,(int*)&skybpl,(int*)&skyxsiz,(int*)&skyysiz);
	if (!skypic)
	{
		long r, g, b, *p;
loadbluesky:;
			//Load default sky
		skyxsiz = 512; skyysiz = 1; skybpl = skyxsiz*4;
		if (!(skypic = (long)malloc(skyysiz*skybpl))) return(-1);

		p = (long *)skypic; y = skyxsiz*skyxsiz;
		for(x=0;x<=(skyxsiz>>1);x++)
		{
			p[x] = ((((x*1081 - skyxsiz*252)*x)/y + 35)<<16)+
					 ((((x* 950 - skyxsiz*198)*x)/y + 53)<<8)+
					  (((x* 439 - skyxsiz* 21)*x)/y + 98);
		}
		p[skyxsiz-1] = 0x50903c;
		r = ((p[skyxsiz>>1]>>16)&255);
		g = ((p[skyxsiz>>1]>>8)&255);
		b = ((p[skyxsiz>>1])&255);
		for(x=(skyxsiz>>1)+1;x<skyxsiz;x++)
		{
			p[x] = ((((0x50-r)*(x-(skyxsiz>>1)))/(skyxsiz-1-(skyxsiz>>1))+r)<<16)+
					 ((((0x90-g)*(x-(skyxsiz>>1)))/(skyxsiz-1-(skyxsiz>>1))+g)<<8)+
					 ((((0x3c-b)*(x-(skyxsiz>>1)))/(skyxsiz-1-(skyxsiz>>1))+b));
		}
		y = skyxsiz*skyysiz;
		for(x=skyxsiz;x<y;x++) p[x] = p[x-skyxsiz];
	}

		//Initialize look-up table for longitudes
	if (skylng) free((void *)skylng);
	if (!(skylng = (point2d *)malloc(skyysiz*8))) return(-1);
	f = PI*2.0 / ((float)skyysiz);
	for(y=skyysiz-1;y>=0;y--)
		fcossin((float)y*f+PI,&skylng[y].x,&skylng[y].y);
	skylngmul = (float)skyysiz/(PI*2);
		//This makes those while loops in gline() not lockup when skyysiz==1
	if (skyysiz == 1) { skylng[0].x = 0; skylng[0].y = 0; }

		//Initialize look-up table for latitudes
	if (skylat) free((void *)skylat);
	if (!(skylat = (long *)malloc(skyxsiz*4))) return(-1);
	f = PI*.5 / ((float)skyxsiz);
	for(x=skyxsiz-1;x;x--)
	{
		ang = (float)((x<<1)-skyxsiz)*f;
		ftol(cos(ang)*32767.0,&xoff);
		ftol(sin(ang)*32767.0,&yoff);
		skylat[x] = (xoff<<16)+((-yoff)&65535);
	}
	skylat[0] = 0; //Hack to make sure assembly index never goes < 0
	skyxsiz--; //Hack for assembly code

	return(0);
}

/**
 * Saves a native Voxlap5 .VXL file & specified position to disk
 * @param filnam .VXL map formatted like this: "UNTITLED.VXL"
 * @param ipo default starting camera position
 * @param ist RIGHT unit vector
 * @param ihe DOWN unit vector
 * @param ifo FORWARD unit vector
 * @return 0:bad, 1:good
 */
long savevxl (const char *savfilnam, dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	FILE *fil;
	long i;

	if (!(fil = fopen(savfilnam,"wb"))) return(0);
	i = 0x09072000; fwrite(&i,4,1,fil);  //Version
	i = VSID; fwrite(&i,4,1,fil);
	i = VSID; fwrite(&i,4,1,fil);
	fwrite(ipo,24,1,fil);
	fwrite(ist,24,1,fil);
	fwrite(ihe,24,1,fil);
	fwrite(ifo,24,1,fil);
	for(i=0;i<VSID*VSID;i++) fwrite(sptr[i],slng(sptr[i]),1,fil);
	fclose(fil);
	return(1);

}

char relpathbase[MAX_PATH];
void relpathinit (char *st)
{
	long i;

	for(i=0;st[i];i++) relpathbase[i] = st[i];
	if ((i) && (relpathbase[i-1] != '/') && (relpathbase[i-1] != '\\'))
		relpathbase[i++] = '\\';
	relpathbase[i] = 0;
}

	//Makes path relative to voxed directory (for sxl).
	//(relpathbase is "C:\kwin\voxlap\" on my machine.)
	//Examples:
	//   "C:\KWIN\VOXLAP\KV6\ANVIL.KV6" -> "KV6\ANVIL.KV6"
	//   "KV6\ANVIL.KV6" -> "KV6\ANVIL.KV6"
static char *relpath (char *st)
{
	long i;
	char ch0, ch1;

	for(i=0;relpathbase[i];i++)
	{
		ch0 = st[i];
		if ((ch0 >= 'a') && (ch0 <= 'z')) ch0 -= 32;
		if (ch0 == '/') ch0 = '\\';

		ch1 = relpathbase[i];
		if ((ch1 >= 'a') && (ch1 <= 'z')) ch1 -= 32;
		if (ch1 == '/') ch1 = '\\';

		if (ch0 != ch1) return(st);
	}
	return(&st[i]);
}

void floatprint (float f, FILE *fil)
{
	char buf[32];

	sprintf(buf,"%g",f);
	if ((buf[0] == '0') && (buf[1] == '.'))
		fprintf(fil,"%s",&buf[1]); //"0." -> "."
	else if ((buf[0] == '-') && (buf[1] == '0') && (buf[2] == '.'))
		{ buf[1] = '-'; fprintf(fil,"%s",&buf[1]); } //"-0." -> "-."
	else fprintf(fil,"%s",buf);
}

void savesxl (char *sxlnam, vx5sprite spr[], long sxlind[] )
{
	FILE *fil;
	long i, j, k;

	if (!(fil = fopen(sxlnam,"wb"))) return;
	fprintf(fil,"%s\r\n",relpath(vxlfilnam));
	fprintf(fil,"%s\r\n",relpath(skyfilnam));

	i = 0; goto savesxlskip;
	do
	{
		if (spr[i].voxnum)
		{
			if (spr[i].flags&2) j = spr[i].kfaptr->namoff;
								else j = spr[i].voxnum->namoff;
			fprintf(fil,"%s,",relpath(getkfilname(j)));
		} else fprintf(fil,"DOT,");
		floatprint(spr[i].p.x,fil); fputc(',',fil);
		floatprint(spr[i].p.y,fil); fputc(',',fil);
		floatprint(spr[i].p.z,fil); fputc(',',fil);
		floatprint(spr[i].s.x,fil); fputc(',',fil);
		floatprint(spr[i].s.y,fil); fputc(',',fil);
		floatprint(spr[i].s.z,fil); fputc(',',fil);
		floatprint(spr[i].h.x,fil); fputc(',',fil);
		floatprint(spr[i].h.y,fil); fputc(',',fil);
		floatprint(spr[i].h.z,fil); fputc(',',fil);
		floatprint(spr[i].f.x,fil); fputc(',',fil);
		floatprint(spr[i].f.y,fil); fputc(',',fil);
		floatprint(spr[i].f.z,fil); fputc(13,fil); fputc(10,fil);
savesxlskip:;
			//User string
		for(j=sxlind[i],k=0;j<sxlind[i+1];j++)
		{
			if (!sxlbuf[j]) { fprintf(fil,"\r\n"); k = 0; continue; }
			if (!k) { fputc(32,fil); k = 1; }
			fputc(sxlbuf[j],fil);
		}

		i++;
	} while (i < numsprites);
	fprintf(fil,"end\r\n");
	fclose(fil);

	sprintf(message,"Saved %s",sxlnam); messagetimeout = totclk+4000;
}


/**
 * Loads a Comanche format map into memory.
 * @param filename Should be formatted like this: "C1.DTA"
 *                 It replaces the first letter with C&D to get both height&color
 * @return 0:bad, 1:good
*/
long loaddta (const char *filename, dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	long i, j, p, leng, minz = 255, maxz = 0, h[5], longpal[256];
	char dat, *dtahei, *dtacol, *v, dafilename[MAX_PATH];
	float f;
	FILE *fp;

	if (!vbuf) { vbuf = (long *)malloc((VOXSIZ>>2)<<2); if (!vbuf) evilquit("loaddta: vbuf malloc failed"); }
	if (!vbit) { vbit = (long *)malloc((VOXSIZ>>7)<<2); if (!vbit) evilquit("loaddta: vbit malloc failed"); }

	if (VSID != 1024) return(0);
	v = (char *)(&vbuf[1]); //1st dword for voxalloc compare logic optimization

	strcpy(dafilename,filename);

	dtahei = (char *)(&vbuf[(VOXSIZ-2097152)>>2]);
	dtacol = (char *)(&vbuf[(VOXSIZ-1048576)>>2]);

	dafilename[0] = 'd';
	if (!kzopen(dafilename)) return(0);
	kzseek(128,SEEK_SET); p = 0;
	while (p < 1024*1024)
	{
		dat = kzgetc();
		if (dat >= 192) { leng = dat-192; dat = kzgetc(); }
					  else { leng = 1; }
		dat = 255-dat;
		if (dat < minz) minz = dat;
		if (dat > maxz) maxz = dat;
		while (leng-- > 0) dtahei[p++] = dat;
	}
	kzclose();

	dafilename[0] = 'c';
	if (!kzopen(dafilename)) return(0);
	kzseek(128,SEEK_SET);
	p = 0;
	while (p < 1024*1024)
	{
		dat = kzgetc();
		if (dat >= 192) { leng = dat-192; dat = kzgetc(); }
					  else { leng = 1; }
		while (leng-- > 0) dtacol[p++] = dat;
	}

	dat = kzgetc();
	if (dat == 0xc)
		for(i=0;i<256;i++)
		{
			longpal[i] = kzgetc();
			longpal[i] = (longpal[i]<<8)+kzgetc();
			longpal[i] = (longpal[i]<<8)+kzgetc() + 0x80000000;
		}

	kzclose();

		//Fill board data
	minz = lbound(128-((minz+maxz)>>1),-minz,255-maxz);
	for(p=0;p<1024*1024;p++)
	{
		h[0] = (long)dtahei[p];
		h[1] = (long)dtahei[((p-1)&0x3ff)+((p     )&0xffc00)];
		h[2] = (long)dtahei[((p+1)&0x3ff)+((p     )&0xffc00)];
		h[3] = (long)dtahei[((p  )&0x3ff)+((p-1024)&0xffc00)];
		h[4] = (long)dtahei[((p  )&0x3ff)+((p+1024)&0xffc00)];

		j = 1;
		for(i=4;i>0;i--) if (h[i]-h[0] > j) j = h[i]-h[0];

		sptr[p] = v;
		v[0] = 0;
		v[1] = dtahei[p]+minz;
		v[2] = dtahei[p]+minz+j-1;
		v[3] = 0; //dummy (z top)
		v += 4;
		for(;j;j--) { *(long *)v = colorjit(longpal[dtacol[p]],0x70707); v += 4; }
	}

	//memset(&sptr[VSID*VSID],0,sizeof(sptr)-VSID*VSID*4);
	vbiti = (((long)v-(long)vbuf)>>2); //# vbuf longs/vbit bits allocated
	clearbuf((void *)vbit,vbiti>>5,-1);
	clearbuf((void *)&vbit[vbiti>>5],(VOXSIZ>>7)-(vbiti>>5),0);
	vbit[vbiti>>5] = (1<<vbiti)-1;

	vx5.globalmass = calcglobalmass();

	ipo->x = VSID*.5; ipo->y = VSID*.5; ipo->z = 128;
	f = 0.0*PI/180.0;
	ist->x = cos(f); ist->y = sin(f); ist->z = 0;
	ihe->x = 0; ihe->y = 0; ihe->z = 1;
	ifo->x = sin(f); ifo->y = -cos(f); ifo->z = 0;

	gmipnum = 1; vx5.flstnum = 0;
	updatebbox(0,0,0,VSID,VSID,MAXZDIM,0);
	return(1);

}

/**
 * Loads a heightmap from PNG (or TGA) into memory; alpha channel is height.
 * @param filename Any 1024x1024 PNG or TGA file with alpha channel.
 * @return 0:bad, 1:good
 */

long loadpng (const char *filename, dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{

	unsigned long *pngdat, dat[5];
	long i, j, k, l, p, leng, minz = 255, maxz = 0;
	char *v, *buf;
	float f;
	FILE *fp;

	if (!vbuf) { vbuf = (long *)malloc((VOXSIZ>>2)<<2); if (!vbuf) evilquit("loadpng: vbuf malloc failed"); }
	if (!vbit) { vbit = (long *)malloc((VOXSIZ>>7)<<2); if (!vbit) evilquit("loadpng: vbuf malloc failed"); }

	if (VSID != 1024) return(0);
	v = (char *)(&vbuf[1]); //1st dword for voxalloc compare logic optimization

	if (!kzopen(filename)) return(0);
	leng = kzfilelength();
	buf = (char *)malloc(leng); if (!buf) { kzclose(); return(0); }
	kzread(buf,leng);
	kzclose();

	kpgetdim(buf,leng,(int*)&i,(int*)&j); if ((i != VSID) && (j != VSID)) { free(buf); return(0); }
	pngdat = (unsigned long *)(&vbuf[(VOXSIZ-VSID*VSID*4)>>2]);
	if (kprender(buf,leng,(long)pngdat,VSID<<2,VSID,VSID,0,0) < 0) return(0);
	free(buf);

	for(i=0;i<VSID*VSID;i++)
	{
		if ((pngdat[i]>>24) < minz) minz = (pngdat[i]>>24);
		if ((pngdat[i]>>24) > maxz) maxz = (pngdat[i]>>24);
	}

		//Fill board data
	minz = lbound(128-((minz+maxz)>>1),-minz,255-maxz);
	for(p=0;p<VSID*VSID;p++)
	{
		dat[0] = pngdat[p];
		dat[1] = pngdat[((p-1)&(VSID-1))+((p     )&((VSID-1)*VSID))];
		dat[2] = pngdat[((p+1)&(VSID-1))+((p     )&((VSID-1)*VSID))];
		dat[3] = pngdat[((p  )&(VSID-1))+((p-VSID)&((VSID-1)*VSID))];
		dat[4] = pngdat[((p  )&(VSID-1))+((p+VSID)&((VSID-1)*VSID))];

		j = 1; l = dat[0];
		for(i=4;i>0;i--)
			if (((signed long)((dat[i]>>24)-(dat[0]>>24))) > j)
				{ j = (dat[i]>>24)-(dat[0]>>24); l = dat[i]; }

		sptr[p] = v;
		v[0] = 0;
		v[1] = (pngdat[p]>>24)+minz;
		v[2] = (pngdat[p]>>24)+minz+j-1;
		v[3] = 0; //dummy (z top)
		v += 4;
		k = (pngdat[p]&0xffffff)|0x80000000;
		if (j == 2)
		{
			l = (((  l     &255)-( k     &255))>>1)      +
				 (((((l>> 8)&255)-((k>> 8)&255))>>1)<< 8) +
				 (((((l>>16)&255)-((k>>16)&255))>>1)<<16);
		}
		else if (j > 2)
		{
			l = (((  l     &255)-( k     &255))/j)      +
				 (((((l>> 8)&255)-((k>> 8)&255))/j)<< 8) +
				 (((((l>>16)&255)-((k>>16)&255))/j)<<16);
		}
		*(long *)v = k; v += 4; j--;
		while (j) { k += l; *(long *)v = colorjit(k,0x30303); v += 4; j--; }
	}

	//memset(&sptr[VSID*VSID],0,sizeof(sptr)-VSID*VSID*4);
	vbiti = (((long)v-(long)vbuf)>>2); //# vbuf longs/vbit bits allocated
	clearbuf((void *)vbit,vbiti>>5,-1);
	clearbuf((void *)&vbit[vbiti>>5],(VOXSIZ>>7)-(vbiti>>5),0);
	vbit[vbiti>>5] = (1<<vbiti)-1;

	vx5.globalmass = calcglobalmass();

	ipo->x = VSID*.5; ipo->y = VSID*.5; ipo->z = 128;
	f = 0.0*PI/180.0;
	ist->x = cos(f); ist->y = sin(f); ist->z = 0;
	ihe->x = 0; ihe->y = 0; ihe->z = 1;
	ifo->x = sin(f); ifo->y = -cos(f); ifo->z = 0;

	gmipnum = 1; vx5.flstnum = 0;
	updatebbox(0,0,0,VSID,VSID,MAXZDIM,0);

	return(1);
}
void * allocate_mem_mapped_file( mem_mapped_file_info file_info, char * file_name, const size_t min_file_size)
{
	int windows_open_mode;
    int open_mode = e_open_mode::if_exists_truncate_if_not_exists_create;

    switch (open_mode)
    {
    case if_exists_fail_if_not_exists_create:
        windows_open_mode = CREATE_NEW;
        break;
    case if_exists_keep_if_dont_exists_fail:
        windows_open_mode = OPEN_EXISTING;
        break;
    case if_exists_keep_if_dont_exists_create:
        windows_open_mode = OPEN_ALWAYS;
        break;
    case if_exists_truncate_if_not_exists_fail:
        windows_open_mode = TRUNCATE_EXISTING;
        break;
    case if_exists_truncate_if_not_exists_create:
        windows_open_mode = CREATE_ALWAYS;
        break;
    default: return NULL;
    }

    //file_info.file_handle_ = CreateFile( file_name, GENERIC_READ | GENERIC_WRITE,
    //    0, 0, windows_open_mode, FILE_ATTRIBUTE_NORMAL, 0);
    //if (file_info.file_handle_ == INVALID_HANDLE_VALUE) return NULL;
    size_t initial_file_size = GetFileSize(file_info.file_handle_, 0);
    size_t adjusted_file_size = initial_file_size == 0 ? min_file_size : initial_file_size;
    file_info.file_mapping_handle_ = CreateFileMapping( INVALID_HANDLE_VALUE, 0, PAGE_WRITECOPY, 0, min_file_size, 0 );
    
    //file_info.file_mapping_handle_ = CreateFileMapping( INVALID_HANDLE_VALUE, 0, PAGE_WRITECOPY, 0, min_file_size, 0 );
    if (file_info.file_mapping_handle_ == INVALID_HANDLE_VALUE) return NULL;
    file_info.size_ = initial_file_size;
    file_info.capacity_ = adjusted_file_size;
    return (MapViewOfFile( file_info.file_mapping_handle_, FILE_MAP_COPY, 0, 0, 0) );
}
void loadnul (dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	lpoint3d lp0, lp1;
	long i, x, y;
	char *v;
	float f;
	#if defined(_WIN32) && defined(USE_MEMMAPPED_IO)

    //file_handle_ = CreateFile("mappedfile", GENERIC_READ,
    //    FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
    //if (file_handle_ == INVALID_HANDLE_VALUE) return;
    //file_mapping_handle_ = CreateFileMapping( file_handle_, 0, PAGE_READONLY, 0, 0, 0 );
    //if (file_mapping_handle_ == INVALID_HANDLE_VALUE) return;
    //vbuf = static_cast<long*>(MapViewOfFile( file_mapping_handle_, FILE_MAP_READ, 0, 0, 0) );
    //if (vbuf) size_ = GetFileSize(file_handle_, 0);
    
    //Resize
    //UnmapViewOfFile(vbuf);
    //CloseHandle(file_mapping_handle_);
    //file_mapping_handle_ = CreateFileMapping( file_handle_, 0, PAGE_READWRITE, 0, ((VOXSIZ>>2)<<2), 0 );
    //capacity_ = ((VOXSIZ>>2)<<2);
    //vbuf = static_cast<long*>(MapViewOfFile( file_mapping_handle_, FILE_MAP_WRITE, 0, 0, 0) );
//#if defined(__unix__)
//    ::munmap(const_cast<char*>(data_), size_);
//    ::close(file_handle_);
//#elif defined(_WIN32)
//    ::UnmapViewOfFile(data_);
//    ::CloseHandle(file_mapping_handle_);
//    ::CloseHandle(file_handle_);
//#endif
	vbuf = static_cast<long*>(allocate_mem_mapped_file( vbuf_file_info, "vbuf_mapped_memory", (((VOXSIZ)>>2)<<2) ));
	if (!vbuf) 
		evilquit("loadnul: vbuf malloc failed");

	vbit = static_cast<long*>(allocate_mem_mapped_file( vbit_file_info, "vbit_mapped_memory", (((VOXSIZ)>>7)<<2) ));
	if (!vbuf) 
		evilquit("loadnul: vbit malloc failed");

	#else
	if (!vbuf) { vbuf = (long *)malloc((VOXSIZ>>2)<<2); if (!vbuf) evilquit("loadnul: vbuf malloc failed"); }
	if (!vbit) { vbit = (long *)malloc((VOXSIZ>>7)<<2); if (!vbit) evilquit("loadnul: vbit malloc failed"); }
	#endif
	
	

	v = (char *)(&vbuf[1]); //1st dword for voxalloc compare logic optimization

		//Completely re-compile vbuf
	for(x=0;x<VSID;x++)
		for(y=0;y<VSID;y++)
		{
			sptr[y*VSID+x] = v;
			i = 0; // + (rand()&1);  //i = default height of plain
			v[0] = 0;
			v[1] = i;
			v[2] = i;
			v[3] = 0;  //z0 (Dummy filler)
			//i = ((((x+y)>>3) + ((x^y)>>4)) % 231) + 16;
			//i = (i<<16)+(i<<8)+i;
			v += 4;
			(*(long *)v) = ((x^y)&15)*0x10101+0x807c7c7c; //colorjit(i,0x70707)|0x80000000;
			v += 4;
		}

	//memset(&sptr[VSID*VSID],0,sizeof(sptr)-VSID*VSID*4);
	vbiti = (((long)v-(long)vbuf)>>2); //# vbuf longs/vbit bits allocated
	clearbuf((void *)vbit,vbiti>>5,-1);
	clearbuf((void *)&vbit[vbiti>>5],(VOXSIZ>>7)-(vbiti>>5),0);
	vbit[vbiti>>5] = (1<<vbiti)-1;

		
	vx5.colfunc = jitcolfunc; vx5.curcol = 0x80704030;

	lp0.x = VSID*.5-90; lp0.y = VSID*.5-90; lp0.z = MAXZDIM*.5-45;
	lp1.x = VSID*.5+90; lp1.y = VSID*.5+90; lp1.z = MAXZDIM*.5+45;
	setrect(&lp0,&lp1,-1);
		//Blow out sphere and stick you inside map
	//lp0.x = VSID*.5; lp0.y = VSID*.5; lp0.z = MAXZDIM*.5; setsphere(&lp0,64,-1);

	vx5.globalmass = calcglobalmass();

	ipo->x = VSID*.5; ipo->y = VSID*.5; ipo->z = MAXZDIM*.5; //ipo->z = -16;
	f = 0.0*PI/180.0;
	ist->x = cos(f); ist->y = sin(f); ist->z = 0;
	ihe->x = 0; ihe->y = 0; ihe->z = 1;
	ifo->x = sin(f); ifo->y = -cos(f); ifo->z = 0;

	gmipnum = 1; vx5.flstnum = 0;
	updatebbox(0,0,0,VSID,VSID,MAXZDIM,0);
}

//Quake3 .BSP loading code begins --------------------------------------------
typedef struct { long c, i; float z, z1; } vlinerectyp;
static point3d q3pln[5250];
static float q3pld[5250], q3vz[256];
static long q3nod[4850][3], q3lf[4850];

long vlinebsp (float x, float y, float z0, float z1, float *dvz)
{
	vlinerectyp vlrec[64];
	float z, t;
	long i, j, vcnt, vlcnt;
	char vt[256];

	vcnt = 1; i = 0; vlcnt = 0; vt[0] = 17;
	while (1)
	{
		if (i < 0)
		{
			if (vt[vcnt-1] != (q3lf[~i]&255))
				{ dvz[vcnt] = z0; vt[vcnt] = (q3lf[~i]&255); vcnt++; }
		}
		else
		{
			j = q3nod[i][0]; z = q3pld[j] - q3pln[j].x*x - q3pln[j].y*y;
			t = q3pln[j].z*z0-z;
			if ((t < 0) == (q3pln[j].z*z1 < z))
				{ vlrec[vlcnt].c = 0; i = q3nod[i][(t<0)+1]; }
			else
			{
				z /= q3pln[j].z; j = (q3pln[j].z<0)+1;
				vlrec[vlcnt].c = 1; vlrec[vlcnt].i = q3nod[i][j];
				vlrec[vlcnt].z = z; vlrec[vlcnt].z1 = z1;
				i = q3nod[i][3-j]; z1 = z;
			}
			vlcnt++; continue;
		}
		do { vlcnt--; if (vlcnt < 0) return(vcnt); } while (!vlrec[vlcnt].c);
		vlrec[vlcnt].c = 0; i = vlrec[vlcnt].i;
		z0 = vlrec[vlcnt].z; z1 = vlrec[vlcnt].z1;
		vlcnt++;
	}
	return(0);
}

/**
 * Loads a Quake 3 Arena .BSP format map into memory. First extract the map
 * from the .PAK file. NOTE: only tested with Q3DM(1,7,17) and Q3TOURNEY
 * @param filnam .BSP map formatted like this: "Q3DM17.BSP"
 * @param ipo default starting camera position
 * @param ist RIGHT unit vector
 * @param ihe DOWN unit vector
 * @param ifo FORWARD unit vector
 */

void loadbsp (const char *filnam, dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	FILE *fp;
	dpoint3d dp;
	float f, xof, yof, zof, sc, rsc;
	long numplanes, numnodes, numleafs, fpos[17], flng[17];
	long i, x, y, z, z0, z1, vcnt, *lptr, minx, miny, minz, maxx, maxy, maxz;
	char *v;

	if (!vbuf) { vbuf = (long *)malloc((VOXSIZ>>2)<<2); if (!vbuf) evilquit("loadbsp: vbuf malloc failed"); }
	if (!vbit) { vbit = (long *)malloc((VOXSIZ>>7)<<2); if (!vbit) evilquit("loadbsp: vbit malloc failed"); }

		//Completely re-compile vbuf
	v = (char *)(&vbuf[1]); //1st dword for voxalloc compare logic optimization
	for(x=0;x<VSID;x++)
		for(y=0;y<VSID;y++)
		{
			(sptr)[y*VSID+x] = v; v[0] = 0; v[1] = 0; v[2] = 0; v[3] = 0; v += 4;
			(*(long *)v) = ((x^y)&15)*0x10101+0x807c7c7c; v += 4;
		}

	//memset((void *)&sptr[VSID*VSID],0,sizeof(sptr)-VSID*VSID*4);
	vbiti = (((long)v-(long)vbuf)>>2); //# vbuf longs/vbit bits allocated
	clearbuf((void *)vbit,vbiti>>5,-1);
	clearbuf((void *)&vbit[vbiti>>5],(VOXSIZ>>7)-(vbiti>>5),0);
	vbit[vbiti>>5] = (1<<vbiti)-1;

	if (!kzopen(filnam)) return;
	kzread(&i,4); if (i != 0x50534249) { kzclose(); return; }
	kzread(&i,4); if (i != 0x2e) { kzclose(); return; }
	for(i=0;i<17;i++) { kzread(&fpos[i],4); kzread(&flng[i],4); }
	kzseek(fpos[2],SEEK_SET); numplanes = flng[2]/16;
	for(i=0;i<numplanes;i++) { kzread(&q3pln[i].x,12); kzread(&q3pld[i],4); }
	kzseek(fpos[3],SEEK_SET); numnodes = flng[3]/36;
	minx = 0x7fffffff; miny = 0x7fffffff; minz = 0x7fffffff;
	maxx = 0x80000000; maxy = 0x80000000; maxz = 0x80000000;
	for(i=0;i<numnodes;i++)
	{
		kzread(&q3nod[i][0],12);
		kzread(&x,4); if (x < minx) minx = x;
		kzread(&x,4); if (x < miny) miny = x;
		kzread(&x,4); if (x < minz) minz = x;
		kzread(&x,4); if (x > maxx) maxx = x;
		kzread(&x,4); if (x > maxy) maxy = x;
		kzread(&x,4); if (x > maxz) maxz = x;
	}
	kzseek(fpos[4]+4,SEEK_SET); numleafs = flng[4]/48;
	for(i=0;i<numleafs;i++) { kzread(&q3lf[i],4); kzseek(44,SEEK_CUR); }
	kzclose();

	sc = (float)(VSID-2)/(float)(maxx-minx);
	rsc = (float)(VSID-2)/(float)(maxy-miny); if (rsc < sc) sc = rsc;
	rsc = (float)(MAXZDIM-2)/(float)(maxz-minz); if (rsc < sc) sc = rsc;
	//i = *(long *)sc; i &= 0xff800000; sc = *(float *)i;
	xof = (-(float)(minx+maxx)*sc + VSID   )*.5;
	yof = (+(float)(miny+maxy)*sc + VSID   )*.5;
	zof = (+(float)(minz+maxz)*sc + MAXZDIM)*.5;

	rsc = 1.0 / sc;
	vx5.colfunc = curcolfunc; //0<x0<x1<VSID, 0<y0<y1<VSID, 0<z0<z1<256,
	for(y=0;y<VSID;y++)
		for(x=0;x<VSID;x++)
		{
			lptr = scum2(x,y);

				//voxx = q3x*+sc + xof;
				//voxy = q3y*-sc + yof;
				//voxz = q3z*-sc + zof;
			vcnt = vlinebsp(((float)x-xof)*rsc,((float)y-yof)*-rsc,-65536.0,65536.0,q3vz);
			for(i=vcnt-2;i>0;i-=2)
			{
				ftol(-q3vz[i+1]*sc+zof,&z0); if (z0 < 0) z0 = 0;
				ftol(-q3vz[i  ]*sc+zof,&z1); if (z1 > MAXZDIM) z1 = MAXZDIM;
				delslab(lptr,z0,z1);
			}
		}
	scum2finish();

	vx5.globalmass = calcglobalmass();

		//Find a spot that isn't too close to a wall
	sc = -1; ipo->x = VSID*.5; ipo->y = VSID*.5; ipo->z = -16;
	for(i=4096;i>=0;i--)
	{
		x = (rand()%VSID); y = (rand()%VSID); z = (rand()%MAXZDIM);
		if (!isvoxelsolid(x,y,z))
		{
			rsc = findmaxcr((double)x+.5,(double)y+.5,(double)z+.5,5.0);
			if (rsc <= sc) continue;
			ipo->x = (double)x+.5; ipo->y = (double)x+.5; ipo->z = (double)x+.5;
			sc = rsc; if (sc >= 5.0) break;
		}
	}
	f = 0.0*PI/180.0;
	ist->x = cos(f); ist->y = sin(f); ist->z = 0;
	ihe->x = 0; ihe->y = 0; ihe->z = 1;
	ifo->x = sin(f); ifo->y = -cos(f); ifo->z = 0;

	gmipnum = 1; vx5.flstnum = 0;
	updatebbox(0,0,0,VSID,VSID,MAXZDIM,0);
}

/** Save .KV6 sprites to disk.
 * This could be a handy function for debugging I suppose. Use it to save
 * .KV6 sprites to disk.
 *
 * @param filnam filename of .KV6 to save to disk. It's your responsibility to
 *               make sure it doesn't overwrite a file of the same name.
 * @param kv pointer to .KV6 object to save to disk.
 */
void savekv6 (const char *filnam, kv6data *kv)
{
	FILE *fil;
	long i;

	if (fil = fopen(filnam,"wb"))
	{
		i = 0x6c78764b; fwrite(&i,4,1,fil); //Kvxl
		fwrite(&kv->xsiz,4,1,fil); fwrite(&kv->ysiz,4,1,fil); fwrite(&kv->zsiz,4,1,fil);
		fwrite(&kv->xpiv,4,1,fil); fwrite(&kv->ypiv,4,1,fil); fwrite(&kv->zpiv,4,1,fil);
		fwrite(&kv->numvoxs,4,1,fil);
		fwrite(kv->vox,kv->numvoxs*sizeof(kv6voxtype),1,fil);
		fwrite(kv->xlen,kv->xsiz*sizeof(long),1,fil);
		fwrite(kv->ylen,kv->xsiz*kv->ysiz*sizeof(short),1,fil);
		fclose(fil);
	}
}

/** @note should make this inline to getkv6!
    \Warning  This can evilquit when misaligned malloc happens.
*/
static kv6data *loadkv6 (const char *filnam)
{
	FILE *fil;
	kv6data tk, *newkv6;
	long i;

	if (!kzopen(filnam))
	{
			//File not found, but allocate a structure anyway
			//   so it can keep track of the filename
		if (!(newkv6 = (kv6data *)malloc(sizeof(kv6data)))) return(0);
		newkv6->leng = sizeof(kv6data);
		newkv6->xsiz = newkv6->ysiz = newkv6->zsiz = 0;
		newkv6->xpiv = newkv6->ypiv = newkv6->zpiv = 0;
		newkv6->numvoxs = 0;
		newkv6->namoff = 0;
		newkv6->lowermip = 0;
		newkv6->vox = (kv6voxtype *)(((long)newkv6)+sizeof(kv6data));
		newkv6->xlen = (unsigned long *)newkv6->vox;
		newkv6->ylen = (unsigned short *)newkv6->xlen;
		return(newkv6);
	}

	kzread((void *)&tk,32);

	i = tk.numvoxs*sizeof(kv6voxtype) + tk.xsiz*4 + tk.xsiz*tk.ysiz*2;
	newkv6 = (kv6data *)malloc(i+sizeof(kv6data));
	if (!newkv6) { kzclose(); return(0); }
	if (((long)newkv6)&3) evilquit("getkv6 malloc not 32-bit aligned!");

	newkv6->leng = i+sizeof(kv6data);
	memcpy(&newkv6->xsiz,&tk.xsiz,28);
	newkv6->namoff = 0;
	newkv6->lowermip = 0;
	newkv6->vox = (kv6voxtype *)(((long)newkv6)+sizeof(kv6data));
	newkv6->xlen = (unsigned long *)(((long)newkv6->vox)+tk.numvoxs*sizeof(kv6voxtype));
	newkv6->ylen = (unsigned short *)(((long)newkv6->xlen) + tk.xsiz*4);

	kzread((void *)newkv6->vox,i);
	kzclose();
	return(newkv6);
}

/**
 * Loads a .KV6 voxel sprite into memory. It malloc's the array for you and
 * returns the pointer to the loaded vx5sprite. If the same filename was
 * passed before to this function, it will return the pointer to the
 * previous instance of the .KV6 buffer in memory (It will NOT load the
 * same file twice). Uninitvoxlap() de-allocates all .KV6 sprites for
 * you.
 *
 * Other advanced info: Uses a 256-entry hash table to compare filenames, so
 * it should be fast. If you want to modify a .KV6 without affecting all
 * instances, you must allocate&de-allocate your own kv6data structure,
 * and use memcpy. The buffer is kv6data.leng bytes long (inclusive).
 *
 * Cover-up function for LOADKV6: returns a pointer to the loaded kv6data
 * structure. Loads file only if not already loaded before with getkv6.
 *
 * @param kv6nam .KV6 filename
 * @return pointer to malloc'ed kv6data structure. Do NOT free this buffer
 *         yourself! Returns 0 if there's an error - such as bad filename.
 */
kv6data *getkv6 (const char *filnam)
{
	kv6data *kv6ptr;
	long i;

	if (inkhash(filnam,&i)) return(*(kv6data **)&khashbuf[i+4]);
	if (i == -1) return(0);

	if (kv6ptr = loadkv6((char *)&khashbuf[i+9]))
		kv6ptr->namoff = i+9; //Must use offset for ptr->name conversion

	*(kv6data **)&khashbuf[i+4] = kv6ptr;
	*(char *)&khashbuf[i+8] = 0; //0 for KV6
	return(kv6ptr);
}

long loadvxl (const char *lodfilnam, dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	FILE *fil;
	long i, j, fsiz;
	char *v, *v2;
	#if defined(_WIN32) && defined(USE_MEMMAPPED_IO)

	vbuf = static_cast<long*>(allocate_mem_mapped_file( vbuf_file_info, "loadvxl_vbuf_mapped_memory", (((VOXSIZ)>>2)<<2) ));
	if (!vbuf) 
		evilquit("loadvxl: vbuf malloc failed");

	vbit = static_cast<long*>(allocate_mem_mapped_file( vbit_file_info, "loadvxl_vbit_mapped_memory", (((VOXSIZ)>>7)<<2) ));
	if (!vbuf) 
		evilquit("loadvxl: vbit malloc failed");
	#else
	if (!vbuf) { vbuf = (long *)malloc((VOXSIZ>>2)<<2); if (!vbuf) evilquit("loadvxl: vbuf malloc failed"); }
	if (!vbit) { vbit = (long *)malloc((VOXSIZ>>7)<<2); if (!vbit) evilquit("loadvxl: vbit malloc failed"); }
	#endif

	if (!kzopen(lodfilnam)) return(0);
	fsiz = kzfilelength();

	kzread(&i,4); if (i != 0x09072000) return(0);
	kzread(&i,4); if (i != VSID) return(0);
	kzread(&i,4); if (i != VSID) return(0);
	kzread(ipo,24);
	kzread(ist,24);
	kzread(ihe,24);
	kzread(ifo,24);

	v = (char *)(&vbuf[1]); //1st dword for voxalloc compare logic optimization
	kzread((void *)v,fsiz-kztell());

	for(i=0;i<VSID*VSID;i++)
	{
		sptr[i] = v;
		while (v[0]) v += (((long)v[0])<<2);
		v += ((((long)v[2])-((long)v[1])+2)<<2);
	}
	kzclose();

	//memset(&sptr[VSID*VSID],0,sizeof(sptr)-VSID*VSID*4);
	vbiti = (((long)v-(long)vbuf)>>2); //# vbuf longs/vbit bits allocated
	clearbuf((void *)vbit,vbiti>>5,-1);
	clearbuf((void *)&vbit[vbiti>>5],(VOXSIZ>>7)-(vbiti>>5),0);
	vbit[vbiti>>5] = (1<<vbiti)-1;

	vx5.globalmass = calcglobalmass();
	backedup = -1;

	gmipnum = 1; vx5.flstnum = 0;
	updatebbox(0,0,0,VSID,VSID,MAXZDIM,0);
	return(1);
}

/** Uses a 256-entry hash to compare names very quickly.
 *  @return 2 values packed into a long
 *  0,retptr=-1: Error! (bad filename or out of memory)
 *	0,retptr>=0: Not in hash; new name allocated, valid index
 *	1,retptr>=0: Already in hash, valid index
*/
long inkhash (const char *filnam, long *retind)
{
	long i, j, hashind;

	(*retind) = -1;

	if (!filnam) return(0);
	j = strlen(filnam); if (!j) return(0);
	j += 10;
	if (khashpos+j > khashsiz) //Make sure string fits in khashbuf
	{
		i = khashsiz; do { i <<= 1; } while (khashpos+j > i);
		if (!(khashbuf = (char *)realloc(khashbuf,i))) return(0);
		khashsiz = i;
	}

		//Copy filename to avoid destroying original string
		//Also, calculate hash index (which hopefully is uniformly random :)
	strcpy(&khashbuf[khashpos+9],filnam);
	for(i=khashpos+9,hashind=0;khashbuf[i];i++)
	{
		if ((khashbuf[i] >= 'a') && (khashbuf[i] <= 'z')) khashbuf[i] -= 32;
		if (khashbuf[i] == '/') khashbuf[i] = '\\';
		hashind = (khashbuf[i] - hashind*3);
	}
	hashind %= (sizeof(khashead)/sizeof(khashead[0]));

		//Find if string is already in hash...
	for(i=khashead[hashind];i>=0;i=(*(long *)&khashbuf[i]))
		if (!strcmp(&khashbuf[i+9],&khashbuf[khashpos+9]))
			{ (*retind) = i; return(1); } //Early out: already in hash

	(*retind) = khashpos;
	*(long *)&khashbuf[khashpos] = khashead[hashind];
	*(long *)&khashbuf[khashpos+4] = 0; //Set to 0 just in case load fails
	khashead[hashind] = khashpos; khashpos += j;
	return(0);
}

/** Skip past single forward slash or backslash in filename.
 *  @param filnam Pointer to string containing filname
 *  @param filnam Pointer to string containing filname after stripping leading / or double backslash
*/
char *stripdir (char *filnam)
{
	long i, j;
	for(i=0,j=-1;filnam[i];i++)
		if ((filnam[i] == '/') || (filnam[i] == '\\')) j = i;
	return(&filnam[j+1]);
}

/**
 * Loads a .KFA file and its associated .KV6 voxel sprite into memory. Works
 * just like getkv6() for for .KFA files: (Returns a pointer to the loaded
 * kfatype structure. Loads data only if not already loaded before with
 * getkfa.)
 *
 * @param kfanam .KFA filename
 * @return pointer to malloc'ed kfatype structure. Do NOT free this buffer
 *         yourself! Returns 0 if there's an error - such as bad filename.
 */
kfatype *getkfa (const char *kfanam)
{
	kfatype *kfa;
	kv6voxtype *v, *ov, *ve;
	kv6data *kv;
	long i, j, x, y;
	char *cptr, snotbuf[MAX_PATH];

	if (inkhash(kfanam,&i)) return(*(kfatype **)&khashbuf[i+4]);
	if (i == -1) return(0);

	if (!(kfa = (kfatype *)malloc(sizeof(kfatype)))) return(0);
	memset(kfa,0,sizeof(kfatype));

	kfa->namoff = i+9; //Must use offset for ptr->name conversion
	*(kfatype **)&khashbuf[i+4] = kfa;
	*(char *)&khashbuf[i+8] = 1; //1 for KFA

	if (!kzopen(kfanam)) return(0);
	kzread(&i,4); if (i != 0x6b6c774b) { kzclose(); return(0); } //Kwlk
	kzread(&i,4); strcpy(snotbuf,kfanam); cptr = stripdir(snotbuf);
	kzread(cptr,i); cptr[i] = 0;
	kzread(&kfa->numhin,4);

		//Actual allocation for ->spr is numspr, which is <= numhin!
	if (!(kfa->spr = (vx5sprite *)malloc(kfa->numhin*sizeof(vx5sprite)))) { kzclose(); return(0); }

	if (!(kfa->hinge = (hingetype *)malloc(kfa->numhin*sizeof(hingetype)))) { kzclose(); return(0); }
	kzread(kfa->hinge,kfa->numhin*sizeof(hingetype));

	kzread(&kfa->numfrm,4);
	if (!(kfa->frmval = (short *)malloc(kfa->numfrm*kfa->numhin*2))) { kzclose(); return(0); }
	kzread(kfa->frmval,kfa->numfrm*kfa->numhin*2);

	kzread(&kfa->seqnum,4);
	if (!(kfa->seq = (seqtyp *)malloc(kfa->seqnum*sizeof(seqtyp)))) { kzclose(); return(0); }
	kzread(kfa->seq,kfa->seqnum*sizeof(seqtyp));

	kzclose();

		//MUST load the associated KV6 AFTER the kzclose :/
	kfa->numspr = 0;
	kv = getkv6(snotbuf); if (!kv) return(0);
	kfa->basekv6 = kv;
	if (!kv->numvoxs) return(0);
	v = kv->vox;
	for(x=kv->numvoxs-1;x>=0;x--) v[x].vis |= ((v[x].vis&32)<<1);
	for(x=0;x<kv->xsiz;x++)
		for(y=0;y<kv->ysiz;y++)
			for(ve=&v[kv->ylen[x*kv->ysiz+y]];v<ve;v++)
			{
				if (v->vis&16) ov = v;
				if ((v->vis&(64+32)) == 64+32)
				{
					floodsucksprite(&kfa->spr[kfa->numspr],kv,x,y,ov,v);
					kfa->numspr++;
				}
			}

	kfa->hingesort = (long *)malloc(kfa->numhin*4);
	kfasorthinge(kfa->hinge,kfa->numhin,kfa->hingesort);

		//Remember position offsets of limbs with no parent in hinge[?].p[0]
	for(i=(kfa->numhin)-1;i>=0;i--)
	{
		j = kfa->hingesort[i]; if (j >= kfa->numspr) continue;
		if (kfa->hinge[j].parent < 0)
		{
			kfa->hinge[j].p[0].x = -kfa->spr[j].p.x;
			kfa->hinge[j].p[0].y = -kfa->spr[j].p.y;
			kfa->hinge[j].p[0].z = -kfa->spr[j].p.z;
		}
	}

	return(kfa);
}

/**
 * Cover-up function to handle both .KV6 and .KFA files. It looks at the
 * filename extension and uses the appropriate function (either getkv6
 * or getkfa) and sets the sprite flags depending on the type of file.
 * The file must have either .KV6 or .KFA as the filename extension. If
 * you want to use weird filenames, then use getkv6/getkfa instead.
 *
 * @param spr Pointer to sprite structure that you provide. getspr() writes:
 *            only to the kv6data/voxtype, kfatim, and flags members.
 * @param filnam filename of either a .KV6 or .KFA file.
 */
void getspr (vx5sprite *s, const char *filnam)
{
	long i;

	if (!filnam) return;
	i = strlen(filnam); if (!i) return;

	if ((filnam[i-1] == 'A') || (filnam[i-1] == 'a'))
		{ s->kfaptr = getkfa(filnam); s->flags = 2; s->kfatim = 0; }
	else if (filnam[i-1] == '6')
		{ s->voxnum = getkv6(filnam); s->flags = 0; }
}
