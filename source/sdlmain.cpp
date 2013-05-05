/**************************************************************************************************
 * SDLMAIN.CPP & SYSMAIN.H                                                                        *
 * By Jonathon Fowler                                                                             *
 **************************************************************************************************/

// This file has been modified from Ken Silverman's original release
//min & max failsafe
#undef max
#undef min

#define GL_GLEXT_PROTOTYPES
//#define NOGL_PBO // Define to turn off pixel buffer objects

#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>

#include <stdio.h>
#include <stdlib.h>

//need this for unused inline function in ksnippits
#include <math.h>

//Ericson2314's dirty porting tricks
#include "porthacks.h"

//Ken's short, general-purpose to-be-inlined functions mainly consisting of inline assembly are now here
#include "ksnippits.h"
#include "cpu_detect.h"
#include "ktimer.h"

#ifndef __i386__
    #error i386 targets only.
#endif

//for code_rwx_unlock only
#ifndef _WIN32
    #include <unistd.h>
    #include <sys/mman.h>
    typedef void * PVOID;
    typedef PVOID HANDLE;
    typedef HANDLE HWND;
#elif defined(_WIN32)
    #include "SDL/SDL_syswm.h"
#endif

#ifdef __APPLE__
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif

#if !defined( NOGL_PBO )
    #include "../include/glext.h"
    GLuint pboIds[2];                                                    //Use 2 buffers to prevent stalls
    GLuint textureId;
    GLubyte* imageData = 0;                                              // pointer to texture buffer
    #ifdef _WIN32
        PFNGLGENBUFFERSARBPROC pglGenBuffersARB = 0;                     // VBO Name Generation Procedure
        PFNGLBINDBUFFERARBPROC pglBindBufferARB = 0;                     // VBO Bind Procedure
        PFNGLBUFFERDATAARBPROC pglBufferDataARB = 0;                     // VBO Data Loading Procedure
        PFNGLBUFFERSUBDATAARBPROC pglBufferSubDataARB = 0;               // VBO Sub Data Loading Procedure
        PFNGLDELETEBUFFERSARBPROC pglDeleteBuffersARB = 0;               // VBO Deletion Procedure
        PFNGLGETBUFFERPARAMETERIVARBPROC pglGetBufferParameterivARB = 0; // return various parameters of VBO
        PFNGLMAPBUFFERARBPROC pglMapBufferARB = 0;                       // map VBO procedure
        PFNGLUNMAPBUFFERARBPROC pglUnmapBufferARB = 0;                   // unmap VBO procedure
        #define glGenBuffersARB           pglGenBuffersARB
        #define glBindBufferARB           pglBindBufferARB
        #define glBufferDataARB           pglBufferDataARB
        #define glBufferSubDataARB        pglBufferSubDataARB
        #define glDeleteBuffersARB        pglDeleteBuffersARB
        #define glGetBufferParameterivARB pglGetBufferParameterivARB
        #define glMapBufferARB            pglMapBufferARB
        #define glUnmapBufferARB          pglUnmapBufferARB
    #endif

    bool has_pbo = false;
    bool gl_init = false;
    const GLenum PIXEL_FORMAT   = GL_BGRA;
#endif

#ifndef NOSOUND
    #define ENTERMUTX SDL_LockAudio()
    #define LEAVEMUTX SDL_UnlockAudio()
	#include "ksound.h"
	#ifdef USEKZ
		//#define KPLIB_C  //if kplib is compiled as C
		#include "kplib.h"
	#endif
#endif

#define SYSMAIN
//#undef _WIN32
//We never want to define C bindings if this is compiled as C++.
#undef SYSMAIN_C	//Putting this here just in case.
#include "sysmain.h"

#define evalmacro(x) evalmacrox(x)
#define evalmacrox(x) #x

char *prognam = "SDLMain App (by KS/TD/JE)";
long progresiz = 1;

extern long initapp (long argc, char **argv);
extern void uninitapp ();
extern void doframe ();

static long quitprogram = 0, quitparam;
void breath();

long ActiveApp = 1, alwaysactive = 0;
long xres = 640, yres = 480, colbits = 8, fullscreen = 1, maxpages = 8;

#if !defined( NOGL_PBO )
    const int CHANNEL_COUNT  = 4;
#endif
long ddrawuseemulation = 0;

void* gethwnd()
{
#ifdef WIN32
    SDL_SysWMinfo wmi;
    SDL_VERSION(&wmi.version);
    if(!SDL_GetWMInfo(&wmi)) return NULL;
    return wmi.window;
#else
    return NULL;
#endif
}

//================== Fast & accurate TIMER FUNCTIONS begins ==================
void initklock ()
{
    perfrq = 1000.0;
    rdtimbase = rdtsc64();
    pertimbase = SDL_GetTicks();
    nextimstep = 4194304;
    klockmul = 0.000000001;
    klockadd = 0.0;
}

void readklock (double *tim)
{
    uint64_t q = rdtsc64()-rdtimbase;
    if (q > nextimstep)
    {
        uint32_t p;
        double d;
        p = SDL_GetTicks();
        d = klockmul;
        klockmul = ((double)(p-pertimbase))/(((double)q)*perfrq);
        klockadd += (d-klockmul)*((double)q);
        do
        {
            nextimstep <<= 1;
        }
        while (q > nextimstep);
    }
    (*tim) = ((double)q)*klockmul + klockadd;
}

//=================== Fast & accurate TIMER FUNCTIONS ends ===================

//SDL Video VARIABLES & CODE--------------------------------------------------
#ifndef NODRAW

#define BASICSURFBITS (SDL_HWSURFACE|SDL_DOUBLEBUF)

PALETTEENTRY pal[256];

static SDL_Surface *sdlsurf = NULL;
static int surflocked = 0;

#define MAXVALIDMODES 256
//typedef struct { long x, y; char c, r0, g0, b0, a0, rn, gn, bn, an; } validmodetype;
static validmodetype validmodelist[MAXVALIDMODES];
static long validmodecnt = 0;
validmodetype curvidmodeinfo;

/** Flip Surface */
void ddflip2gdi(void)
{
    if ( SDL_Flip(sdlsurf)  != 0)
    {
        // Flip failed
    }
    else
    {
        // Everything is fine
    }
}

static void gvmladd(long x, long y, char c)
{
    if (validmodecnt == MAXVALIDMODES) return;
    memset(&validmodelist[validmodecnt], 0, sizeof(validmodetype));
    validmodelist[validmodecnt].x = x;
    validmodelist[validmodecnt].y = y;
    validmodelist[validmodecnt].c = c;
    validmodecnt++;
}

long getvalidmodelist (validmodetype **davalidmodelist)
{
    int cd[5] = { 8,15,16,24,32 }, i, j;
    SDL_Rect **modes;
    SDL_PixelFormat pf = { NULL, 8, 1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0 };

    if (!validmodecnt)
    {
        for (i=0; i<5; i++)
        {
            pf.BitsPerPixel = cd[i];
            pf.BytesPerPixel = cd[i] >> 3;
            modes = SDL_ListModes(&pf, BASICSURFBITS|SDL_FULLSCREEN);
            if (modes == (SDL_Rect**)0)
            {
                // none available
            } if (modes == (SDL_Rect**)-1)
            {
                // add some defaults
                gvmladd(640, 480, cd[i]);
                gvmladd(800, 600, cd[i]);
                gvmladd(1024, 768, cd[i]);
                gvmladd(1152, 864, cd[i]);
                gvmladd(1280, 960, cd[i]);
                gvmladd(1280, 1024, cd[i]);
                gvmladd(1600, 1200, cd[i]);
            }
            else
            {
                for (j=0; modes[j]; j++)
                    gvmladd(modes[j]->w, modes[j]->h, cd[i]);
            }
        }
    }
    (*davalidmodelist) = validmodelist;
    return(validmodecnt);
}

void updatepalette (long start, long danum)
{
    SDL_Color pe[256];
    int i,j;

    if (!sdlsurf) return;

    for (i=0,j=danum; j>0; i++,j--)
    {
        pe[i].r = pal[start+i].peRed;
        pe[i].g = pal[start+i].peGreen;
        pe[i].b = pal[start+i].peBlue;
        pe[i].unused = 0;
    }
    SDL_SetColors(sdlsurf, pe, start, danum);
}

void stopdirectdraw()
{
    if (!surflocked) return;

#if !defined( NOGL_PBO )
    if (has_pbo)
    {
        glUnmapBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB); // release pointer to mapping buffer
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        glBindTexture(GL_TEXTURE_2D, textureId);
        glColor4f(1, 1, 1, 1);
        glBegin(GL_QUADS);
        glNormal3f(0, 0, -1);
        glTexCoord2f(0.0f, 0.0f);
        glVertex3f( -1.0f, 1.0f, 0.0f);
        glTexCoord2f(1.0f, 0.0f);
        glVertex3f( 1.0f, 1.0f, 0.0f);
        glTexCoord2f(1.0f, 1.0f);
        glVertex3f(1.0f,  -1.0f, 0.0f);
        glTexCoord2f(0.0f, 1.0f);
        glVertex3f( -1.0f,  -1.0f, 0.0f);
        glEnd();

        // unbind texture
        glBindTexture(GL_TEXTURE_2D, 0);

    }
    else
#endif
        if (sdlsurf && SDL_MUSTLOCK(sdlsurf)) SDL_UnlockSurface(sdlsurf);

    surflocked = 0;
}
#if !defined(NOGL_PBO)
void updatePixels(void* dst, int size)
{
    static int color = 0x00000000;

    if(!dst)
        return;

    int* ptr = (int*)(GLubyte*)dst;

    // copy 4 bytes at once
    for(int i = 0; i < yres; ++i)
    {
        for(int j = 0; j < xres; ++j)
        {
            *ptr = color;
            ++ptr;
        }
        color += 23;   // add an arbitary number (no meaning)
    }
    ++color;            // scroll down
}
#endif

long startdirectdraw(long *vidplc, long *dabpl, long *daxres, long *dayres)
{
    if (surflocked) stopdirectdraw();

    static int index = 0;
    int nextIndex = 0;

#if !defined(NOGL_PBO)
    if ( has_pbo )
    {
        *daxres = sdlsurf->w;
        *dayres = sdlsurf->h;
        *dabpl = sdlsurf->pitch;

        index = (index + 1) % 2;
        nextIndex = (index + 1) % 2;
        // bind the texture and PBO
        glBindTexture(GL_TEXTURE_2D, textureId );
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pboIds[index]);

        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, xres, yres, PIXEL_FORMAT, GL_UNSIGNED_BYTE, 0);
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pboIds[nextIndex]);
        glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, (xres) * (yres) *CHANNEL_COUNT, 0, GL_STREAM_DRAW_ARB);

        GLubyte* ptr = (GLubyte*)glMapBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, GL_WRITE_ONLY_ARB);
        *vidplc = (long)ptr;
        //memset(ptr, 0, xres*yres*CHANNEL_COUNT);
        //*vidplc = (long)sdlsurf->pixels;
    }
    else
    {
#endif
        if (SDL_LockSurface(sdlsurf) < 0) return 0;
        *dabpl = sdlsurf->pitch;
        *daxres = sdlsurf->w;
        *dayres = sdlsurf->h;
        *vidplc = (long)sdlsurf->pixels;
#if !defined(NOGL_PBO)
    }
#endif

    surflocked = 1;
    return 1;
}

void nextpage()
{
    if (!sdlsurf) return;
    //if (surflocked) stopdirectdraw();
#if !defined(NOGL_PBO)
    if (has_pbo)
    {
        
        SDL_GL_SwapBuffers();
    }
    else
#endif
        SDL_Flip(sdlsurf);
}

long clearscreen(long fillcolor)
{
#if !defined(NOGL_PBO)
    if (has_pbo){
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        return 0;
    }
    else{
#endif
        SDL_FillRect(sdlsurf, NULL, fillcolor);
        return 0;
#if !defined(NOGL_PBO)
    }
#endif
}

#if !defined(NOGL_PBO)
void initopengl_withpbo(uint32_t * surfbits)
{
    *surfbits |= SDL_OPENGL;
    *surfbits |= SDL_HWACCEL;

    SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 8 );
    SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 8 );
    SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 8 );
    SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 0 );
    SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 0 );
    SDL_GL_SetAttribute( SDL_GL_ACCELERATED_VISUAL, 1 );
    SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
    SDL_GL_SetAttribute( SDL_GL_SWAP_CONTROL, 0 );
}

bool initSharedMem()
{
    // allocate texture buffer
    imageData = new GLubyte[xres*yres*CHANNEL_COUNT];
    memset(imageData, 0, xres*yres*CHANNEL_COUNT);

    glShadeModel(GL_FLAT);                      // shading mathod: GL_SMOOTH or GL_FLAT
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);      // 4-byte pixel alignment

    // enable /disable features
    //@glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    //glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    //@glEnable(GL_LIGHTING);
    glDisable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_CULL_FACE);

    glGenTextures(1, &textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8_EXT
                 , xres
                 , yres
                 , 0
                 , PIXEL_FORMAT
                 , GL_UNSIGNED_BYTE
                 , (GLvoid*)imageData);
    glBindTexture(GL_TEXTURE_2D, 0);

    return true;
}

void clearSharedMem()
{
    // deallocate texture buffer
    delete [] imageData;
    imageData = 0;

    // clean up texture
    glDeleteTextures(1, &textureId);

    // clean up PBOs
    if( has_pbo )
    {
        glDeleteBuffersARB(2, pboIds);
    }
}

void init_opengl_extensions()
{
    initSharedMem();

 

    /* Get function pointer */
#ifdef _WIN32
    glGenBuffersARB =           (PFNGLGENBUFFERSARBPROC)            wglGetProcAddress("glGenBuffersARB");
    glBindBufferARB =           (PFNGLBINDBUFFERARBPROC)            wglGetProcAddress("glBindBufferARB");
    glBufferDataARB =           (PFNGLBUFFERDATAARBPROC)            wglGetProcAddress("glBufferDataARB");
    glBufferSubDataARB =        (PFNGLBUFFERSUBDATAARBPROC)         wglGetProcAddress("glBufferSubDataARB");
    glDeleteBuffersARB =        (PFNGLDELETEBUFFERSARBPROC)         wglGetProcAddress("glDeleteBuffersARB");
    glGetBufferParameterivARB = (PFNGLGETBUFFERPARAMETERIVARBPROC)  wglGetProcAddress("glGetBufferParameterivARB");
    glMapBufferARB =            (PFNGLMAPBUFFERARBPROC)             wglGetProcAddress("glMapBufferARB");
    glUnmapBufferARB =          (PFNGLUNMAPBUFFERARBPROC)           wglGetProcAddress("glUnmapBufferARB");

    // check once again PBO extension
    if(glGenBuffersARB && glBindBufferARB && glBufferDataARB && glBufferSubDataARB &&
            glMapBufferARB && glUnmapBufferARB && glDeleteBuffersARB && glGetBufferParameterivARB)
    {
        has_pbo = true;
    }
    else
    {
        has_pbo = false;
    }
#endif
    fprintf(stdout, "Screen BPP: %d\n", SDL_GetVideoSurface()->format->BitsPerPixel);

    fprintf( stdout, "hw_available: %d\n", SDL_GetVideoInfo()->hw_available );
    fprintf( stdout, "blit_hw: %d\n", SDL_GetVideoInfo()->blit_hw );
    fprintf( stdout, "blit_hw_A: %d\n", SDL_GetVideoInfo()->blit_hw_A );
    fprintf( stdout, "video_mem: %d\n", SDL_GetVideoInfo()->video_mem );


    fprintf( stdout, "Vendor     : %s\n", glGetString( GL_VENDOR ) );
    fprintf( stdout, "Renderer   : %s\n", glGetString( GL_RENDERER ) );
    fprintf( stdout, "Version    : %s\n", glGetString( GL_VERSION ) );

    if (has_pbo)
    {
        fprintf(stdout, "Video card supports GL_ARB_pixel_buffer_object.\n");
    }
    else
    {
        fprintf(stderr, "****NO video card support for GL_ARB_pixel_buffer_object.****\n");
    }

    fprintf( stdout, "Extensions : %s\n", glGetString( GL_EXTENSIONS ) );
    fprintf( stdout, "\n");
    if (has_pbo)
    {
        glGenBuffersARB(2, pboIds);
        fprintf(stdout, "Created 2 pixel buffer objects: %d and %d of size: %d \n", pboIds[0],pboIds[1], xres*yres*CHANNEL_COUNT );
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pboIds[0]);
        glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, xres*yres*CHANNEL_COUNT, 0, GL_STREAM_DRAW_ARB);
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pboIds[1]);
        glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, xres*yres*CHANNEL_COUNT, 0, GL_STREAM_DRAW_ARB);
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
    }
    gl_init = true;
}
#endif

long initdirectdraw(long daxres, long dayres, long dacolbits)
{
    uint32_t surfbits = 0;

#ifndef NOINPUT
    extern long mouse_acquire;
    if (mouse_acquire)
    {
        SDL_WM_GrabInput(SDL_GRAB_OFF);
        SDL_ShowCursor(SDL_ENABLE);
    }
#endif

    xres = daxres;
    yres = dayres;
    colbits = dacolbits;

    if (fullscreen) surfbits |= SDL_FULLSCREEN;
    else if (progresiz) surfbits |= SDL_RESIZABLE; // Needs fix

#ifndef NOGL_PBO
    // Map an opengl PBO to use as a write buffer, instead of SDL framebuffer
    initopengl_withpbo( &surfbits );
    sdlsurf = SDL_SetVideoMode(daxres, dayres, 0, surfbits);
#else
    surfbits = BASICSURFBITS;
    surfbits |= SDL_HWPALETTE;
    surfbits |= SDL_HWSURFACE;
    sdlsurf = SDL_SetVideoMode(daxres, dayres, dacolbits, surfbits);
#endif

    if (!sdlsurf)
    {
        fputs("video init failed!",stderr);
        return 0;
    }

#ifndef NOGL_PBO
    // For PBO extension
    if (gl_init)
    {
        clearSharedMem();
    }
    init_opengl_extensions();
#endif

    // FIXME: should I?
    xres = sdlsurf->w;
    yres = sdlsurf->h;
    colbits = sdlsurf->format->BitsPerPixel;
    fullscreen = (sdlsurf->flags & SDL_FULLSCREEN) == SDL_FULLSCREEN;

    memset(&curvidmodeinfo, 0, sizeof(validmodetype));
    curvidmodeinfo.x = sdlsurf->w;
    curvidmodeinfo.y = sdlsurf->h;
    curvidmodeinfo.c = sdlsurf->format->BitsPerPixel;
    curvidmodeinfo.r0 = sdlsurf->format->Rshift;
    curvidmodeinfo.rn = 8 - sdlsurf->format->Rloss;
    curvidmodeinfo.g0 = sdlsurf->format->Gshift;
    curvidmodeinfo.gn = 8 - sdlsurf->format->Gloss;
    curvidmodeinfo.b0 = sdlsurf->format->Bshift;
    curvidmodeinfo.bn = 8 - sdlsurf->format->Bloss;
    curvidmodeinfo.a0 = sdlsurf->format->Ashift;
    curvidmodeinfo.an = 8 - sdlsurf->format->Aloss;

    if (colbits == 8) updatepalette(0,256);

#ifndef NOINPUT
    if (mouse_acquire)
    {
        SDL_WM_GrabInput(SDL_GRAB_ON);
        SDL_ShowCursor(SDL_DISABLE);
    }
#endif

    return 1;
}

long changeres(long daxres, long dayres, long dacolbits, long dafullscreen)
{
    fullscreen = dafullscreen;
    return initdirectdraw(daxres, dayres, dacolbits);
}

#endif

//SDL Input VARIABLES & CODE--------------------------------------------------
char keystatus[256];
static long shkeystatus = 0;
#define KEYBUFSIZ 256
static long keybuf[KEYBUFSIZ], keybufr = 0, keybufw = 0, keybufw2 = 0;

#ifndef NOINPUT

char ext_keystatus[256]; // +TD
char ext_mbstatus[8] = {0}; // +TD extended mouse button status
long mouse_acquire = 1, kbd_acquire = 1;
static void (*setmousein)(long, long) = NULL;
static long mousex = 0, mousey = 0, gbstatus = 0;
long mousmoth=0;
float mousper=0.0;

void readkeyboard ()
{
    return; // all gets done in sdlmsgloop now
}

void readmouse (float *fmousx, float *fmousy, long *bstatus)
{
    if (!mouse_acquire)
    {
        *fmousx = *fmousy = 0;
        *bstatus = 0;
        return;
    }

    *fmousx = (float)mousex;
    *fmousy = (float)mousey;
    *bstatus = gbstatus;

    mousex = mousey = 0;
}

#endif

// Kensound code begins -------------------------------------------------------

#ifndef NOSOUND

#include <math.h>

//Same as: stricmp(st0,st1) except: '/' == '\'

#define PI 3.14159265358979323

void sdlmixcallback(void *userdata, uint8_t *stream, int len)
{
    long amt;

    if (!streambuf) return;

    // copy over the normal kensound mixing buffer
    amt = MIN(streambuf->buflen - streambuf->writecur, len);
    memcpy(stream, &((uint8_t*)streambuf->buf)[ streambuf->writecur ], amt);
    streambuf->writecur += amt;
    if (streambuf->writecur == streambuf->buflen) streambuf->writecur = 0;
    stream += amt;
    len -= amt;
    if (len > 0)
    {
        memcpy(stream, &((uint8_t*)streambuf->buf)[ streambuf->writecur ], len);
        streambuf->writecur += len;
    }

    // mix in the umixer buffers
}

static void kensoundclose ()
{
    numrendersnd = 0;

    if (streambuf) SDL_CloseAudio();

    if (audhashbuf)
    {
        long i;
        for(i=audlastloadedwav; i>=0; i=(*(long *)&audhashbuf[i+4]))
            free((void *)(*(long *)&audhashbuf[i+8]));
        free(audhashbuf);
        audhashbuf = 0;
    }
    audhashpos = audhashsiz = 0;

    if (streambuf)
    {
        KSoundBuffer_Release(streambuf);
        streambuf = 0;
    }
}

static int kensoundinit (int samprate, int numchannels, int bytespersamp)
{
    SDL_AudioSpec audwant,audgot;
    WAVEFORMATEX wft;
    KSBufferDesc bufdesc;
    void *w0, *w1;
    unsigned long l0, l1;

    if (streambuf) kensoundclose();

    audwant.freq = samprate;
    audwant.format = bytespersamp == 1 ? AUDIO_U8 : AUDIO_S16;
    audwant.channels = numchannels;
    audwant.samples = MINBREATHREND;
    audwant.callback = sdlmixcallback;
    audwant.userdata = NULL;

    if (SDL_OpenAudio(&audwant, &audgot) < 0) return -1;
    samprate = audgot.freq;
    switch (audgot.format)
    {
    case AUDIO_U8:
        bytespersamp = 1;
        break;
    case AUDIO_S16:
        bytespersamp = 2;
        break;
    default:
        SDL_CloseAudio();
        return -1;
    }
    numchannels = audgot.channels;

    samplerate = samprate;
    rsamplerate = 1.0/(float)samplerate;
    numspeakers = numchannels;
    bytespersample = bytespersamp;
    gshiftval = (numspeakers+bytespersample-2);

    wft.wFormatTag = WAVE_FORMAT_PCM;
    wft.nSamplesPerSec = samprate;
    wft.wBitsPerSample = (bytespersamp<<3);
    wft.nChannels = numchannels;
    wft.nBlockAlign = bytespersamp*wft.nChannels;
    wft.nAvgBytesPerSec = wft.nSamplesPerSec * wft.nBlockAlign;
    wft.cbSize = 0;

    bufdesc.dwBufferBytes = SNDSTREAMSIZ;
    bufdesc.lpwfxFormat = &wft;
    if (KSound_CreateSoundBuffer(&bufdesc,&streambuf) != KS_OK)
    {
        fprintf(stderr, "kensoundinit(): failed creating sound buffer\n");
        SDL_CloseAudio();
        return(-1);
    }

    //streambuf->SetVolume(curmusicvolume);
    //streambuf->SetFrequency(curmusicfrequency);
    //streambuf->SetPan(curmusicpan);

    //Zero out streaming buffer before beginning play
    KSoundBuffer_Lock(streambuf, 0,SNDSTREAMSIZ,&w0,&l0,&w1,&l1,0);
    if (w0) memset(w0,(bytespersample-2)&128,l0);
    if (w1) memset(w1,(bytespersample-2)&128,l1);
    //streambuf->Unlock(w0,l0,w1,l1);

    initfilters();

    //streambuf->SetCurrentPosition(0);
    //streambuf->Play(0,0,DSBPLAY_LOOPING);

    SDL_PauseAudio(0);

    return(0);
}

#endif

void setvolume (long percentmax)
{
    //Dummy function to get rid of linking error for game.cpp
}

void quitloop ()
{
    SDL_Event ev;
    ev.type = SDL_QUIT;
    SDL_PushEvent(&ev);
}

void evilquit (const char *str) //Evil because this function makes awful assumptions!!!
{
#ifndef NOSOUND
    kensoundclose();
#endif
#if !defined(NOGL_PBO)
    if ( gl_init )
        clearSharedMem();
#endif
    SDL_Quit();
    if (str) fprintf(stderr,"fatal error! %s\n",str);
    uninitapp();
    exit(0);
}

#ifndef NOINPUT
void setacquire (long mouse, long kbd)
{
    SDL_GrabMode g;

    // SDL doesn't let the mouse and keyboard be grabbed independantly
    if ((mouse || kbd) != mouse_acquire)
    {
        g = SDL_WM_GrabInput( (mouse || kbd) ? SDL_GRAB_ON : SDL_GRAB_OFF );
        mouse_acquire = kbd_acquire = (g == SDL_GRAB_ON);

        SDL_ShowCursor(mouse_acquire ? SDL_DISABLE : SDL_ENABLE);
    }
}

void setmouseout (void (*in)(long,long), long x, long y)
{
    if (fullscreen) return;
    setmousein = in;
    setacquire(0, kbd_acquire);

    SDL_WarpMouse(x,y);
}
#endif

static unsigned char keytranslation[SDLK_LAST] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 14, 15, 0, 0, 0, 28, 0, 0, 0, 0, 0, 89, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 57, 2, 40, 4, 5, 6, 8, 40, 10, 11, 9, 13, 51,
    12, 52, 53, 11, 2, 3, 4, 5, 6, 7, 8, 9, 10, 39, 39, 51, 13, 52, 53, 3, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    26, 43, 27, 7, 12, 41, 30, 48, 46, 32, 18, 33, 34, 35, 23, 36, 37, 38, 50,
    49, 24, 25, 16, 19, 31, 20, 22, 47, 17, 45, 21, 44, 0, 0, 0, 0, 211, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 82, 79, 80, 81, 75, 76, 77, 71, 72, 73, 83, 181, 55, 74, 78, 156, 0,
    200, 208, 205, 203, 210, 199, 207, 201, 209, 59, 60, 61, 62, 63, 64, 65,
    66, 67, 68, 87, 88, 0, 0, 0, 0, 0, 0, 69, 58, 70, 54, 42, 157, 29, 184, 56,
    0, 0, 219, 220, 0, 0, 0, /*-2*/0, 84, 183, 221, 0, 0, 0
};

static void sdlmsgloop(void)
{
    SDL_Event ev;
    unsigned char sc;

    while (SDL_PollEvent(&ev))
    {
        switch (ev.type)
        {
        case SDL_ACTIVEEVENT:
            if (ev.active.state & SDL_APPINPUTFOCUS)
            {
                shkeystatus = 0;
                gbstatus = 0;
                ActiveApp = ev.active.gain;
                if (mouse_acquire)
                {
                    if (ActiveApp)
                    {
                        SDL_WM_GrabInput(SDL_GRAB_ON);
                        SDL_ShowCursor(SDL_DISABLE);
                    }
                    else
                    {
                        SDL_WM_GrabInput(SDL_GRAB_OFF);
                        SDL_ShowCursor(SDL_ENABLE);
                    }
                }
            }
            break;
#ifndef NODRAW
        case SDL_VIDEORESIZE:
            ActiveApp = 1;
            initdirectdraw(ev.resize.w, ev.resize.h, colbits);
            break;
#endif
        case SDL_KEYDOWN:
            sc = keytranslation[ev.key.keysym.sym];
            keystatus[ sc ] = 1;   // FIXME: verify this is kosher
            ext_keystatus[ sc ] = 1|2;
            switch (ev.key.keysym.sym)
            {
            case SDLK_LSHIFT:
                shkeystatus |= (1<<16);
                break;
            case SDLK_RSHIFT:
                shkeystatus |= (1<<17);
                break;
            case SDLK_LCTRL:
                shkeystatus |= (1<<18);
                break;
            case SDLK_RCTRL:
                shkeystatus |= (1<<19);
                break;
            case SDLK_LALT:
                shkeystatus |= (1<<20);
                break;
            case SDLK_RALT:
                shkeystatus |= (1<<21);
                break;
            default:
            {
                long i = ((keybufw2+1)&(KEYBUFSIZ-1));
                keybuf[keybufw2&(KEYBUFSIZ-1)] = ((long)sc << 8)+shkeystatus;   // FIXME: verify
                if (i != keybufr) keybufw2 = i; //prevent fifo overlap
            }
            break;
            }
            // The WM_CHAR should be in here somewhere
            if (ev.key.keysym.unicode != 0 && (ev.key.keysym.unicode & 0xff80) == 0)
            {
                if (keybufw2 != keybufr) //stick ASCII code in last FIFO value
                    keybuf[(keybufw2-1)&(KEYBUFSIZ-1)] |= (ev.key.keysym.unicode & 0x7f);
            }
            break;
        case SDL_KEYUP:
            sc = keytranslation[ev.key.keysym.sym];
            keystatus[ sc ] = 0;   // FIXME: verify this is kosher
            ext_keystatus[ sc ] &= ~1; // preserve bit 2 only
            switch (ev.key.keysym.sym)
            {
            case SDLK_LSHIFT:
                shkeystatus &= ~(1<<16);
                break;
            case SDLK_RSHIFT:
                shkeystatus &= ~(1<<17);
                break;
            case SDLK_LCTRL:
                shkeystatus &= ~(1<<18);
                break;
            case SDLK_RCTRL:
                shkeystatus &= ~(1<<19);
                break;
            case SDLK_LALT:
                shkeystatus &= ~(1<<20);
                break;
            case SDLK_RALT:
                shkeystatus &= ~(1<<21);
                break;
            default:
                break;
            }
            break;
#ifndef NOINPUT
            // FIXME: mousewheel?
            //if (zdelta > 0) ext_mbstatus[6] = 2;
            //else if (zdelta < 0) ext_mbstatus[7] = 2;
        case SDL_MOUSEMOTION:
            //      if (gpMouse && ActiveApp && setmousein) {
            //         mouse_acquire = 1;
            //         gpMouse->Acquire(); gbstatus = 0;
            //         setmousein(mx, my);
            //         setmousein = NULL;
            //      }
            if (ActiveApp)
            {
                mousex += ev.motion.xrel;
                mousey += ev.motion.yrel;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
        case SDL_MOUSEBUTTONUP:
        {
            int i=-1;
            switch (ev.button.button)
            {
            case SDL_BUTTON_LEFT:
                i = 0;
                break;
            case SDL_BUTTON_RIGHT:
                i = 1;
                break;
            case SDL_BUTTON_MIDDLE:
                i = 2;
                break;
            }
            if (i<0) break;

            if (ev.type == SDL_MOUSEBUTTONDOWN)
            {
                ext_mbstatus[i] = 1|2;
                gbstatus |= 1<<i;
            }
            else
            {
                ext_mbstatus[i] &= 2;
                gbstatus &= ~(1<<i);
            }
        }
        break;
#endif
        case SDL_QUIT:
#ifndef NOSOUND
            kensoundclose();
#endif
            SDL_Quit();
            uninitapp();
            quitprogram = 1;
            break;
        }
    }
}

long keyread ()
{
    long i;

    if (keybufr == keybufw) return(0);
    i = keybuf[keybufr];
    keybufr = ((keybufr+1)&(KEYBUFSIZ-1));
    return(i);
}

void breath ()
{
    sdlmsgloop();
#ifndef NOSOUND
    ENTERMUTX;
    kensoundbreath(MINBREATHREND);
    LEAVEMUTX;
#endif
}

#if defined(__GNUC__) && !( defined(_WIN32) || defined(_WINDOWS_) || (defined(__MINGW32__) || defined(__MINGW64__)) ) //glibc doesn't have this
static void strlwr(char *s)
{
    for ( ; *s; s++)
    {
        if ((*s) >= 'A' || (*s) <= 'Z') *s &= ~0x20;
    }
}
#endif
//Call like this: arg2filename(argv[1],".ksm",curfilename);
//Make sure curfilename length is lon(wParam&255)g enough!
void arg2filename (const char *oarg, const char *ext, char *narg)
{
    long i;

    //Chop off quotes at beginning and end of filename
    for(i=strlen(oarg); i>0; i--)
        if (oarg[i-1] == '\\') break;
    if ((!i) && (oarg[0] == '\"')) i = 1;
    strcpy(narg,&oarg[i]);
    if (narg[0] == 0) return;
    if (narg[strlen(narg)-1] == '\"') narg[strlen(narg)-1] = 0;
    strlwr(narg);
    if (!strstr(narg,ext)) strcat(narg,ext);
}

int main(int argc, char **argv)
{
    uint32_t sdlinitflags;
    int i;

    cputype = getcputype();
    if ((cputype&((1<<0)+(1<<4))) != ((1<<0)+(1<<4)))
    {
        fputs("Sorry, this program requires FPU&RDTSC support (>=Pentium)", stderr);
        return(-1);
    }

    sdlinitflags = SDL_INIT_TIMER;
#ifndef NOINPUT
    for(i=0; i<256; i++) keystatus[i] = 0;
    for(i=0; i<256; i++) ext_keystatus[i] = 0;
#endif
#ifndef NODRAW
    sdlinitflags |= SDL_INIT_VIDEO;
#endif
#ifndef NOSOUND
    sdlinitflags |= SDL_INIT_AUDIO;
#endif

    if (SDL_Init(sdlinitflags) < 0)
    {
        fputs("Failure initialising SDL.", stderr);
        return -1;
    }
    atexit(SDL_Quit);   // In case we exit() somewhere

    initklock();

    if (initapp(argc, argv) < 0) return -1;

#ifndef NODRAW
    if (!initdirectdraw(xres,yres,colbits))
    {
        SDL_Quit();
        return 0;
    }
#endif

#ifndef NOSOUND
    kensoundinit(44100,2,2);
#endif

    breath();
    while (!quitprogram)
    {
        if (alwaysactive || ActiveApp)
        {
            fpuinit(0x2);
            doframe();
        }
        else SDL_WaitEvent(NULL);
        breath();
    }

#if !defined( NOGL_PBO)
    if ( gl_init )
        clearSharedMem();
#endif

    return quitparam;
}

/*
 * vim:ts=4:
 */
