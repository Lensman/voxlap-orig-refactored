#pragma once
#include <stdarg.h>
#include <stdio.h>

#include "kglobals.h"
extern "C" tiletype gdd;
const long font4x6[] = //256 DOS chars, from Ken's Build SMALLFNT
{
	0x000000,0x6f9f60,0x69f960,0xaffe40,0x4efe40,0x6ff6f0,0x66f6f0,0x000000,
	0xeeaee0,0x000000,0x000000,0x000000,0x000000,0x000000,0x7755c0,0x96f690,
	0x8cec80,0x26e620,0x4e4e40,0xaaa0a0,0x7dd550,0x7ca6c0,0x000ee0,0x4e4ee0,
	0x4e4440,0x444e40,0x02f200,0x04f400,0x000000,0x000000,0x000000,0x000000,
	0x000000,0x444040,0xaa0000,0xafafa0,0x6c46c0,0xa248a0,0x4a4ce0,0x240000,
	0x488840,0x422240,0x0a4a00,0x04e400,0x000224,0x00e000,0x000040,0x224480,
	0xeaaae0,0x444440,0xe2e8e0,0xe2e2e0,0xaae220,0xe8e2e0,0xe8eae0,0xe22220,
	0xeaeae0,0xeae220,0x040400,0x040480,0x248420,0x0e0e00,0x842480,0xc24040,
	0xeaece0,0x4aeaa0,0xcacac0,0x688860,0xcaaac0,0xe8c8e0,0xe8c880,0xe8aae0,
	0xaaeaa0,0xe444e0,0xe22a60,0xaacaa0,0x8888e0,0xaeeaa0,0xaeeea0,0xeaaae0,
	0xeae880,0xeaae60,0xeacaa0,0xe8e2e0,0xe44440,0xaaaae0,0xaaa440,0xaaeea0,
	0xaa4aa0,0xaae440,0xe248e0,0xc888c0,0x844220,0x622260,0x4a0000,0x0000e0,
	0x420000,0x006a60,0x88eae0,0x00e8e0,0x22eae0,0x006e60,0x24e440,0x06a62c,
	0x88eaa0,0x040440,0x040448,0x88aca0,0x444440,0x08eee0,0x00caa0,0x00eae0,
	0x00eae8,0x00eae2,0x00e880,0x0064c0,0x04e440,0x00aa60,0x00aa40,0x00eee0,
	0x00a4a0,0x00aa6c,0x00c460,0x648460,0x440440,0xc424c0,0x6c0000,0x04ae00,
	0x68886c,0xa0aa60,0x606e60,0xe06a60,0xa06a60,0xc06a60,0x046a60,0x00e8e4,
	0xe06e60,0xa06e60,0xc06e60,0x0a0440,0x0e0440,0x0c0440,0xa4aea0,0x404ea0,
	0x60ece0,0x007a70,0x7afab0,0xe0eae0,0xa0eae0,0xc0eae0,0xe0aa60,0xc0aa60,
	0xa0aa6c,0xa0eae0,0xa0aae0,0x4e8e40,0x65c4f0,0xa4ee40,0xcafab0,0x64e4c0,
	0x606a60,0x060440,0x60eae0,0x60aa60,0xc0caa0,0xe0aea0,0x6a60e0,0xeae0e0,
	0x404860,0x007400,0x00c400,0x8a4e60,0x8a6e20,0x404440,0x05a500,0x0a5a00,
	0x282828,0x5a5a5a,0xd7d7d7,0x444444,0x44c444,0x44cc44,0x66e666,0x00e666,
	0x00cc44,0x66ee66,0x666666,0x00ee66,0x66ee00,0x66e000,0x44cc00,0x00c444,
	0x447000,0x44f000,0x00f444,0x447444,0x00f000,0x44f444,0x447744,0x667666,
	0x667700,0x007766,0x66ff00,0x00ff66,0x667766,0x00ff00,0x66ff66,0x44ff00,
	0x66f000,0x00ff44,0x00f666,0x667000,0x447700,0x007744,0x007666,0x66f666,
	0x44ff44,0x44c000,0x007444,0xffffff,0x000fff,0xcccccc,0x333333,0xfff000,
	0x00dad0,0xcacac8,0xea8880,0x00f660,0xe848e0,0x007a40,0x0aac80,0x05a220,
	0xe4a4e0,0x4aea40,0x6996f0,0x646ae0,0x06f600,0x16f680,0x68c860,0x4aaaa0,
	0xe0e0e0,0x4e40e0,0x4240e0,0x4840e0,0x4a8880,0x222a40,0x40e040,0x06ec00,
	0xeae000,0x0cc000,0x00c000,0x644c40,0xcaa000,0xc4e000,0x0eee00,0x000000,
};
	//NOTE: font is stored vertically first! (like .ART files)
const int64_t font6x8[] = //256 DOS chars, from: DOSAPP.FON (tab blank)
{
	0x3E00000000000000,0x6F6B3E003E455145,0x1C3E7C3E1C003E6B,0x3000183C7E3C1800,
	0x7E5C180030367F36,0x000018180000185C,0x0000FFFFE7E7FFFF,0xDBDBC3FF00000000,
	0x0E364A483000FFC3,0x6000062979290600,0x0A7E600004023F70,0x2A1C361C2A003F35,
	0x0800081C3E7F0000,0x7F361400007F3E1C,0x005F005F00001436,0x22007F017F090600,
	0x606060002259554D,0x14B6FFB614000060,0x100004067F060400,0x3E08080010307F30,
	0x08083E1C0800081C,0x0800404040407800,0x3F3C3000083E083E,0x030F3F0F0300303C,
	0x0000000000000000,0x0003070000065F06,0x247E247E24000307,0x630000126A2B2400,
	0x5649360063640813,0x0000030700005020,0x00000000413E0000,0x1C3E080000003E41,
	0x08083E080800083E,0x0800000060E00000,0x6060000008080808,0x0204081020000000,
	0x00003E4549513E00,0x4951620000407F42,0x3649494922004649,0x2F00107F12141800,
	0x494A3C0031494949,0x0305097101003049,0x0600364949493600,0x6C6C00001E294949,
	0x00006CEC00000000,0x2400004122140800,0x2241000024242424,0x0609590102000814,
	0x7E001E555D413E00,0x49497F007E111111,0x224141413E003649,0x7F003E4141417F00,
	0x09097F0041494949,0x7A4949413E000109,0x00007F0808087F00,0x4040300000417F41,
	0x412214087F003F40,0x7F00404040407F00,0x04027F007F020402,0x3E4141413E007F08,
	0x3E00060909097F00,0x09097F005E215141,0x3249494926006619,0x3F0001017F010100,
	0x40201F003F404040,0x3F403C403F001F20,0x0700631408146300,0x4549710007087008,
	0x0041417F00000043,0x0000201008040200,0x01020400007F4141,0x8080808080800402,
	0x2000000007030000,0x44447F0078545454,0x2844444438003844,0x38007F4444443800,
	0x097E080008545454,0x7CA4A4A418000009,0x0000007804047F00,0x8480400000407D00,
	0x004428107F00007D,0x7C0000407F000000,0x04047C0078041804,0x3844444438000078,
	0x380038444444FC00,0x44784400FC444444,0x2054545408000804,0x3C000024443E0400,
	0x40201C00007C2040,0x3C6030603C001C20,0x9C00006C10106C00,0x54546400003C60A0,
	0x0041413E0800004C,0x0000000077000000,0x02010200083E4141,0x3C2623263C000001,
	0x3D001221E1A11E00,0x54543800007D2040,0x7855555520000955,0x2000785554552000,
	0x5557200078545555,0x1422E2A21C007857,0x3800085555553800,0x5555380008555455,
	0x00417C0100000854,0x0000004279020000,0x2429700000407C01,0x782F252F78007029,
	0x3400455554547C00,0x7F097E0058547C54,0x0039454538004949,0x3900003944453800,
	0x21413C0000384445,0x007C20413D00007D,0x3D00003D60A19C00,0x40413C00003D4242,
	0x002466241800003D,0x29006249493E4800,0x16097F00292A7C2A,0x02097E8840001078,
	0x0000785555542000,0x4544380000417D00,0x007D21403C000039,0x7A0000710A097A00,
	0x5555080000792211,0x004E51514E005E55,0x3C0020404D483000,0x0404040404040404,
	0x506A4C0817001C04,0x0000782A34081700,0x0014080000307D30,0x0814000814001408,
	0x55AA114411441144,0xEEBBEEBB55AA55AA,0x0000FF000000EEBB,0x0A0A0000FF080808,
	0xFF00FF080000FF0A,0x0000F808F8080000,0xFB0A0000FE0A0A0A,0xFF00FF000000FF00,
	0x0000FE02FA0A0000,0x0F0800000F080B0A,0x0F0A0A0A00000F08,0x0000F80808080000,
	0x080808080F000000,0xF808080808080F08,0x0808FF0000000808,0x0808080808080808,
	0xFF0000000808FF08,0x0808FF00FF000A0A,0xFE000A0A0B080F00,0x0B080B0A0A0AFA02,
	0x0A0AFA02FA0A0A0A,0x0A0A0A0AFB00FF00,0xFB00FB0A0A0A0A0A,0x0A0A0B0A0A0A0A0A,
	0x0A0A08080F080F08,0xF808F8080A0AFA0A,0x08080F080F000808,0x00000A0A0F000000,
	0xF808F8000A0AFE00,0x0808FF00FF080808,0x08080A0AFB0A0A0A,0xF800000000000F08,
	0xFFFFFFFFFFFF0808,0xFFFFF0F0F0F0F0F0,0xFF000000000000FF,0x0F0F0F0F0F0FFFFF,
	0xFE00241824241800,0x01017F0000344A4A,0x027E027E02000003,0x1800006349556300,
	0x2020FC00041C2424,0x000478040800001C,0x3E00085577550800,0x02724C00003E4949,
	0x0030595522004C72,0x1800182418241800,0x2A2A1C0018247E24,0x003C02023C00002A,
	0x0000002A2A2A2A00,0x4A4A510000242E24,0x00514A4A44000044,0x20000402FC000000,
	0x2A08080000003F40,0x0012241224000808,0x0000000609090600,0x0008000000001818,
	0x02023E4030000000,0x0900000E010E0100,0x3C3C3C0000000A0D,0x000000000000003C,
};

/** Draws 4x6 font on screen (very fast!)
 *
 *  @param x x of top-left corner
 *  @param y y of top-left corner
 *  @param fcol foreground color (32-bit RGB format)
 *  @param bcol background color (32-bit RGB format) or -1 for transparent
 *  @param fmt string - same syntax as printf
 */
void print4x6 (long x, long y, long fcol, const long bcol, const char *fmt, ...)
{
	va_list arglist;
	unsigned char st[280], *c;
	long i, j;

	if (!fmt) return;
	va_start(arglist,fmt);
	vsprintf((char*)st,fmt,arglist);
	va_end(arglist);

	y = y*bytesperline+(x<<2)+frameplace;
	if (bcol < 0)
	{
		for(j=20;j>=0;y+=bytesperline,j-=4)
			for(c=st,x=y;*c;c++,x+=16)
			{
				i = (font4x6[*c]>>j);
				if (i&8) *(long *)(x   ) = fcol;
				if (i&4) *(long *)(x+ 4) = fcol;
				if (i&2) *(long *)(x+ 8) = fcol;
				if (i&1) *(long *)(x+12) = fcol;
				if ((*c) == 9) x += 32;
			}
		return;
	}
	fcol -= bcol;
	for(j=20;j>=0;y+=bytesperline,j-=4)
		for(c=st,x=y;*c;c++,x+=16)
		{
			i = (font4x6[*c]>>j);
			*(long *)(x   ) = (((i<<28)>>31)&fcol)+bcol;
			*(long *)(x+ 4) = (((i<<29)>>31)&fcol)+bcol;
			*(long *)(x+ 8) = (((i<<30)>>31)&fcol)+bcol;
			*(long *)(x+12) = (((i<<31)>>31)&fcol)+bcol;
			if ((*c) == 9) { for(i=16;i<48;i+=4) *(long *)(x+i) = bcol; x += 32; }
		}
}

/** Draws 6x8 font on screen (very fast!)
 *
 *  @param x x of top-left corner
 *  @param y y of top-left corner
 *  @param fcol foreground color (32-bit RGB format)
 *  @param bcol background color (32-bit RGB format) or -1 for transparent
 *  @param fmt string - same syntax as printf
 */
void print6x8 (long x, long y, long fcol, const long bcol, char *fmt, ...)
{
	va_list arglist;
	char st[280], *c, *v;
	long i, j;

	if (!fmt) return;
	va_start(arglist,fmt);
	vsprintf(st,fmt,arglist);
	va_end(arglist);

	y = y*bytesperline+(x<<2)+frameplace;
	if (bcol < 0)
	{
		for(j=1;j<256;y+=bytesperline,j<<=1)
			for(c=st,x=y;*c;c++,x+=24)
			{
				v = (char *)(((long)font6x8) + ((long)c[0])*6);
				if (v[0]&j) *(long *)(x   ) = fcol;
				if (v[1]&j) *(long *)(x+ 4) = fcol;
				if (v[2]&j) *(long *)(x+ 8) = fcol;
				if (v[3]&j) *(long *)(x+12) = fcol;
				if (v[4]&j) *(long *)(x+16) = fcol;
				if (v[5]&j) *(long *)(x+20) = fcol;
				if ((*c) == 9) x += ((2*6)<<2);
			}
		return;
	}
	fcol -= bcol;
	for(j=1;j<256;y+=bytesperline,j<<=1)
		for(c=st,x=y;*c;c++,x+=24)
		{
			v = (char *)(((long)font6x8) + ((long)c[0])*6);
			*(long *)(x   ) = (((-(v[0]&j))>>31)&fcol)+bcol;
			*(long *)(x+ 4) = (((-(v[1]&j))>>31)&fcol)+bcol;
			*(long *)(x+ 8) = (((-(v[2]&j))>>31)&fcol)+bcol;
			*(long *)(x+12) = (((-(v[3]&j))>>31)&fcol)+bcol;
			*(long *)(x+16) = (((-(v[4]&j))>>31)&fcol)+bcol;
			*(long *)(x+20) = (((-(v[5]&j))>>31)&fcol)+bcol;
			if ((*c) == 9) { for(i=24;i<72;i+=4) *(long *)(x+i) = bcol; x += ((2*6)<<2); }
		}
}

/** Draws 6x8 font on screen (very fast!) specific to slab6
 *
 *  @param ox x of top-left corner
 *  @param y y of top-left corner
 *  @param fcol foreground color (32-bit RGB format)
 *  @param bcol background color (32-bit RGB format) or -1 for transparent
 *  @param fmt string - same syntax as printf
 */
void print6x8_s6 (long ox, long y, long fcol, const long bcol, const char *fmt, ...)
{
	va_list arglist;
	char st[1024], *c, *v, *cp, *cpx;
	long i, j, ie, x;

	if (!fmt) return;
	va_start(arglist,fmt);
	if (_vsnprintf((char *)&st,sizeof(st)-1,fmt,arglist)) st[sizeof(st)-1] = 0;
	va_end(arglist);

	cp = (char *)(y*gdd.p+gdd.f);
	for(j=1;j<256;y++,cp=(char *)(((long)cp)+gdd.p),j+=j)
		if ((unsigned)y < (unsigned)gdd.y)
			for(c=st,x=ox;*c;c++,x+=6)
			{
				v = (char *)(((long)*c)*6 + (long)font6x8);
				ie = MIN(gdd.x-x,6); cpx = &cp[x];
				for(i=MAX(-x,0);i<ie;i++) { if (v[i]&j) cpx[i] = fcol; else if (bcol >= 0) cpx[i] = bcol; }
				if ((*c) == 9) { if (bcol >= 0) for(i=6;i<18;i++) cpx[i] = bcol; x += 2*6; }
			}
}
