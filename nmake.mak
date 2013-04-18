# sub-directory macros
locSRC                 =$(MAKEDIR)/source
locINC                 =$(MAKEDIR)/include
locLIB                 =$(MAKEDIR)/libraries
locBINroot             =$(MAKEDIR)/binaries

# adding local equivalents to environment variables
INCLUDE                =$(locINC);$(INCLUDE)
LIB                    =$(locLIB);$(LIB)
#PATH                  =$(locBIN);$(PATH)

# -----------------------------------
# Importing Environment Variables
!IFNDEF USEV5ASM
USEV5ASM               =1
!ENDIF

!IFDEF _ASMNAME
AsmName                =$(_ASMNAME)
!ELSE
AsmName                =masm
!ENDIF

!IFNDEF HAVE_OPENGL
HAVE_OPENGL =1
!ENDIF

!IFNDEF BUILD
BUILD                  =Debug
!ENDIF

!IFNDEF _GFX
GFX                    =win
GFX                    =$(_GFX)
!ELSE
GFX                    =$(_GFX)
!ENDIF


!IF "$(BUILD)"=="Release"
CXX_MODE               =$(CXX_Release)
LNK_MODE               =$(LNK_Release)
locBIN                 =$(locBINroot)/release_$(GFX)
!ELSE IF "$(BUILD)"=="Debug"
CXX_MODE               =$(CXX_Debug)
LNK_MODE               =$(LNK_Debug)
locBIN                 =$(locBINroot)/debug_$(GFX)
!ENDIF

# END Importing Environment Variables
# -----------------------------------

# -----------------------------------
# Micrososft Visual C Macros

CXX                    =cl   #for Micrsoft Compiler
LNK                    =link #for Microsfoft Linker

# Flags
CXXFLAGS               =/Fo$(@R) /fp:fast /arch:SSE2 /c /J $(C_TYPE) $(CXX_MODE) $(GFX_CFLAGS) $(Random_Macros) /I $(locINC) # for Micrsoft Compiler(cl)
CXX_Debug              =$(GFX_CXX_Debug) /ZI /Fx /Gy /GF /Gs /MD /Fdbinaries\ /RTCsuc /Od
CXX_Release            =$(GFX_CXX_Release) /ZI /Fx /Gy /GF /Gs /MD /Fdbinaries\

# can use to fix some issues : /NODEFAULTLIB:msvcrt.lib
LNKFLAGS               =/out:$(@) /SUBSYSTEM:WINDOWS $(LNK_MODE)# for Microsfoft Linker (link)
LNK_Debug              =/DEBUG
LNK_Release            =

# END Micrososft Visual C Macros
# -----------------------------------

# -----------------------------------
# Assembler Macros
# Japheth (Open) Watcomm Assembler (jwasm)
!IF "$(AsmName)"=="jwasm"
AS                     =jwasm
AFLAGS                 =-Fo$(@R) -c -coff -8 -DWIN32
!ENDIF

# Netwide Assembler (nasm)
!IF "$(AsmName)"=="nasm"
AS                     =nasm
AFLAGS                 =-o $(@) -f win32 -DWIN32 --prefix _
!ENDIF

# Micrsoft Macro Assembler (masm)
!IF "$(AsmName)"=="masm"
AS                     =ml
AFLAGS                 =/Fo$(@R) /c /coff
!ENDIF
# END Assembler Macros
# -----------------------------------

# -----------------------------------
# Graphics Backend





!IF "$(GFX)"=="win"
C_TYPE                 =/TP
GFX_LIBS               =ddraw.lib dinput.lib dxguid.lib


GFX_CXX_Debug          =/MD
GFX_CXX_Release        =/MD

!ELSE IF "$(GFX)"=="sdl"
C_TYPE                 =/TP
GFX_LIBS               =SDL.lib SDLmain.lib opengl32.lib

GFX_CXX_Debug          =/MDd
GFX_CXX_Release        =/MD
!ENDIF

# END Importing Environment Variables
# -----------------------------------

# -----------------------------------
# Toggle Random Macros
!IF "$(USEV5ASM)"=="1"
if_USEV5ASM            =$(locBIN)/v5$(OBJSuf)
!ENDIF
Random_Macros          =/DUSEV5ASM=$(USEV5ASM)
# END Toggle Random Macros
# -----------------------------------

# -----------------------------------
# generalObjs : The general set of object files to link, for a program with graphics
# List is in order
# -----------------------------------

GEN_CFLAGS             =/DSYSMAIN_C /DKPLIB_C /DUSEKZ /DUSESCREEN /DFILEIO

generalObjs            =$(locBIN)/cpu_detect$(OBJSuf) \
                        $(locBIN)/kscreen$(OBJSuf) \
                        $(locBIN)/kfile_io$(OBJSuf) \
                        $(locBIN)/voxlap5$(OBJSuf) \
                        $(locBIN)/kfastvoxel$(OBJSuf) \
                        $(locBIN)/kvoxel_modelling$(OBJSuf) \
                        $(if_USEV5ASM) \
                        $(locBIN)/kcolors$(OBJSuf) \
                        $(locBIN)/kplib$(OBJSuf)

# -----------------------------------
# Libraries to link
# -----------------------------------
simpleLIBs             =$(GFX_LIBS) user32.lib gdi32.lib ole32.lib
gameLIBs               =$(GFX_LIBS) user32.lib gdi32.lib ole32.lib
voxedLIBs              =$(GFX_LIBS) user32.lib gdi32.lib           comdlg32.lib
kwalkLIBs              =$(GFX_LIBS) user32.lib gdi32.lib ole32.lib comdlg32.lib

OBJSuf                 =.obj
EXESuf                 =.exe

Phony:                        all
all:                          voxlap slab6
voxlap:                       $(locBIN)/simple$(EXESuf) \
                              $(locBIN)/game$(EXESuf) \
                              $(locBIN)/voxed$(EXESuf) \
                              $(locBIN)/kwalk$(EXESuf)

# executable ($(EXESuf)) (meta)targets
simple:                       $(locBIN)/simple$(EXESuf)
$(locBIN)/simple$(EXESuf):    $(locBIN)/simple$(OBJSuf) $(locBIN)/$(GFX)main1$(OBJSuf)
	$(LNK) $(LNKFLAGS)        $(locBIN)/simple$(OBJSuf) $(generalObjs) $(locBIN)/$(GFX)main1$(OBJSuf) $(simpleLIBs)

game:                         $(locBIN)/game$(EXESuf)
$(locBIN)/game$(EXESuf):      $(locBIN)/game$(OBJSuf) $(locBIN)/$(GFX)main1$(OBJSuf)
	$(LNK) $(LNKFLAGS)        $(locBIN)/game$(OBJSuf) $(generalObjs) $(locBIN)/$(GFX)main1$(OBJSuf) $(gameLIBs)

voxed:                        $(locBIN)/voxed$(EXESuf)
$(locBIN)/voxed$(EXESuf):     $(locBIN)/voxed$(OBJSuf) $(locBIN)/$(GFX)main2$(OBJSuf)
	$(LNK) $(LNKFLAGS)        $(locBIN)/voxed$(OBJSuf) $(generalObjs) $(locBIN)/$(GFX)main2$(OBJSuf) $(voxedLIBs)

kwalk:                        $(locBIN)/kwalk$(EXESuf)
$(locBIN)/kwalk$(EXESuf):     $(locBIN)/kwalk$(OBJSuf) $(locBIN)/$(GFX)main2$(OBJSuf)
	$(LNK) $(LNKFLAGS)        $(locBIN)/kwalk$(OBJSuf) $(generalObjs) $(locBIN)/$(GFX)main2$(OBJSuf) $(kwalkLIBs)

slab6:                        $(locBIN)/slab6$(EXESuf)
$(locBIN)/slab6$(EXESuf):     $(locBIN)/slab6$(OBJSuf) $(locBIN)/s6$(OBJSuf) $(locBIN)/slab6.res                          $(locBIN)/$(GFX)main2$(OBJSuf)
	$(LNK) $(LNKFLAGS)        $(locBIN)/slab6$(OBJSuf) $(generalObjs) $(locBIN)/s6$(OBJSuf) $(locBIN)/slab6.res                          $(locBIN)/$(GFX)main2$(OBJSuf) $(kwalkLIBs)

# binary object ($(OBJSuf)) targets


# Primary Object

simple_o:                          $(locSRC)/simple.cpp
$(locBIN)/simple$(OBJSuf):         $(locSRC)/simple.cpp $(if_USEV5ASM) $(locINC)/voxlap5.h $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/simple.cpp /DUSEKZ /DUSESCREEN /DFILEIO
#   used to use /QIfist

# binary object ($(OBJSuf)) targets
game_o:                            $(locSRC)/game.cpp
$(locBIN)/game$(OBJSuf):           $(locSRC)/game.cpp   $(locINC)/voxlap5.h $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/game.cpp /DUSEKZ /DUSESCREEN /DFILEIO
#   used to use /QIfist

# Primary Object
simple_o:                          $(locSRC)/simple.cpp
$(locBIN)/simple$(OBJSuf):         $(locSRC)/simple.cpp $(locINC)/voxlap5.h $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/simple.cpp /DUSEKZ /DUSESCREEN /DFILEIO
#   used to use /QIfist

voxed_o:                           $(locSRC)/voxed.cpp
$(locBIN)/voxed$(OBJSuf):          $(locSRC)/voxed.cpp  $(locINC)/voxlap5.h $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/voxed.cpp /DUSEKZ /DUSESCREEN /DFILEIO /DKVOXEL_MODELLING /DUSEPHYSICS

kwalk_o:                           $(locSRC)/kwalk.cpp
$(locBIN)/kwalk$(OBJSuf):          $(locSRC)/kwalk.cpp  $(locINC)/voxlap5.h $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/kwalk.cpp /DUSEKZ /DUSESCREEN /DFILEIO /DKVOXEL_MODELLING /DUSEPHYSICS

slab6_o:                           $(locSRC)/slab6.cpp
$(locBIN)/slab6$(OBJSuf):          $(locSRC)/slab6.cpp  $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/slab6.cpp /DUSEKZ /DUSESCREEN /DFILEIO

slab6_res:
$(locBIN)/slab6.res:               $(locSRC)/slab6.rc $(locSRC)/slab6.ico
	rc -r /fo $(locBIN)/slab6.res  $(locSRC)/slab6.rc

# Secondary Objects
voxlap:                            $(locSRC)/modules/voxlap5.cpp
$(locBIN)/voxlap5$(OBJSuf):        $(locSRC)/modules/voxlap5.cpp  $(if_USEV5ASM)  $(locINC)/voxlap5.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/modules/voxlap5.cpp  /DUSEKZ /DUSESCREEN

v5:                                $(locBIN)/v5$(OBJSuf)
$(locBIN)/v5$(OBJSuf):             $(locSRC)/modules/v5.$(AsmName)
	$(AS)  $(AFLAGS)               $(locSRC)/modules/v5.$(AsmName)

s6:                                $(locBIN)/v5$(OBJSuf)
$(locBIN)/s6$(OBJSuf):             $(locSRC)/modules/s6.$(AsmName)
	$(AS)  $(AFLAGS)               $(locSRC)/modules/s6.$(AsmName)

main1:                             $(locBIN)/$(GFX)main1$(OBJSuf)
$(locBIN)/$(GFX)main1$(OBJSuf):    $(locSRC)/$(GFX)main.cpp $(locINC)/voxlap5.h $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/$(GFX)main.cpp /DUSEKZ /DUSESCREEN /DZOOM_TEST

main2:                             $(locBIN)/$(GFX)main2$(OBJSuf)
$(locBIN)/$(GFX)main2$(OBJSuf):    $(locSRC)/$(GFX)main.cpp $(locINC)/voxlap5.h $(locINC)/sysmain.h
	$(CXX) $(CXXFLAGS)             $(locSRC)/$(GFX)main.cpp /DUSEKZ /DUSESCREEN

#-------------------------------
# Clean Scripts
#-------------------------------

# cleanall cleans release and debug builds
cleanall: tidyall
    cd $(locBINroot)
	del /f /S "*$(OBJSuf)"
	del /f /S "*$(EXESuf)"

# Remove cruft for all builds, but leave executables ( so entire folder can be zipped for nightly builds )
tidyall:
    cd $(locBINroot)
    del /f /S "$(locBINroot)/*std*.txt"
    del /f /S "$(locBINroot)/*.*db"
    del /f /S "$(locBINroot)/*.ilk"
    del /f /S "$(locBINroot)/*.bak"
    @echo Directories cleared

# clean only cleans current build set
clean: tidy
    cd "$(locBIN)"
    del /f /S "*$(OBJSuf)"
    del /f /Q "*$(EXESuf)"

# Remove cruft from defined build, but leave executables ( so entire folder can be zipped for nightly builds )
tidy:
    cd "$(locBIN)"
	del /f /S "$(locBIN)/*std*.txt"
	del /f /S "$(locBIN)/*.*db"
    del /f /S "$(locBIN)/*.ilk"
    @echo Directories cleared
#-------------------------------
# Module build system
#-------------------------------
BUILD_CMD=nmake -y -f nmake.mak
modules:
    $(BUILD_CMD) module name=voxlap5
    $(BUILD_CMD) module name=cpu_detect
    $(BUILD_CMD) module name=kvoxel_modelling
    $(BUILD_CMD) module name=kfile_io
    $(BUILD_CMD) module name=kplib
    $(BUILD_CMD) module name=kscreen
    $(BUILD_CMD) module name=kfastvoxel
    $(BUILD_CMD) module name=kcolors


!IF "$(USEV5ASM)"=="1"
    $(BUILD_CMD) assembly name=s6
    $(BUILD_CMD) assembly name=v5
!ENDIF

module:                        $(locBIN)/$(name)$(OBJSuf)
$(locBIN)/$(name)$(OBJSuf):    $(locSRC)/modules/$(name).cpp $(locINC)/$(name).h $(locINC)/kglobals.h
	$(CXX) $(CXXFLAGS)         $(locSRC)/modules/$(name).cpp

assembly:
    $(AS)  $(AFLAGS)               $(locSRC)/modules/$(name).$(AsmName)
