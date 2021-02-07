#============================================================
#====================== MAKEFILE ============================
#============================================================

# Shell utilisé
SHELL = /bin/sh

# Compilateur FORTRAN utilisé
FF = gfortran

# Vérifie le type de compilateur
#COMPILER = $(shell $(FF) -v | cut -d' ' -f1)

#============================================================
#================ Options et variables ======================
#============================================================
FFLAGS           = -Wall -O3 -llapack
FFLAGS_OPT       = -O3 -fPIC -Ofast
FFLAGS_DEBUG     = -g -fPIC -debug
FFLAGS_INSPECTOR = -O0 -g -check none -fPIC  -warn all -warn errors
FFLIBS           = -L/usr/local/opt/lapack/lib -llapack


#============================================================
#================= Liste des dossiers =======================
#============================================================
DATADIR = Data
SRCDIR  = F90
OBJDIR  = Objet
BINDIR  = Bin


ifeq ($(COMPILER),ifort)
  MODULEFLAG := -module   # INTEL
else
  MODULEFLAG = -J     # GNU etc.
endif


#============================================================
#============= Définition objets à complier =================
#============================================================
TARGET = iso_diff.exe iso_diff_abs.exe rayl_diff.exe hg_diff.exe polyn_diff.exe rayl_RGB.exe stereo_RGB.exe test_bhmie.exe test_bhmie_rb.exe env_eff.exe mie_diff.exe mie_diff_CloudySky.exe mie_RGB.exe mie_stereo_RGB.exe

OBJFILES_ISO_DIFF = mod_param.o mod_random.o mie.o mod_diff.o mod_iso.o iso_diff.o
OBJFILES_ISO_ABS =  mod_param.o mod_random.o mie.o mod_diff mod_iso.o iso_diff_abs.o
OBJFILES_RAYL_DIFF = mod_param.o mod_random.o mie.o mod_diff.o mod_proj_stereo.o rayl_diff.o
OBJFILES_HG_DIFF = mod_param.o mod_random.o mie.o mod_diff.o hg_diff.o
OBJFILES_POLYN_DIFF = mod_param.o mod_random.o mie.o mod_diff.o polyn_diff.o
OBJFILES_RAYL_RGB = mod_param.o mod_random.o mie.o mod_diff.o rayl_RGB.o
OBJFILES_RAYL_STEREO_RGB = mod_param.o mod_proj_stereo.o rayl_stereo_RGB.o
OBJFILES_TEST_BHMIE = mod_param.o mie.o mod_diff.o test_bhmie.o
OBJFILES_TEST_BHMIE_RB = mod_param.o mie.o mod_proj_stereo.o test_bhmie_rb.o
OBJFILES_ENV_EFF = mod_param.o mod_random.o mie.o mod_diff.o mod_proj_stereo.o mod_enveloppe.o env_eff.o
OBJFILES_MIE_DIFF = mod_param.o mod_random.o mod_proj_stereo.o mod_fzero.o mie.o mod_diff.o mod_enveloppe.o mod_mie.o mie_diff.o 
OBJFILES_MIE_DIFF_CLOUDYSKY = mod_param.o mod_random.o mod_proj_stereo.o mod_fzero.o mie.o mod_diff.o mod_enveloppe.o mod_mie.o mie_diff_CloudySky.o
OBJFILES_MIE_RGB = mod_param.o mod_random.o mod_proj_stereo.o mod_fzero.o mie.o mod_diff.o mod_enveloppe.o mod_mie.o mie_RGB.o
OBJFILES_MIE_STEREO_RGB = mod_param.o mod_proj_stereo.o mie_stereo_RGB.o


#============================================================
#============ Complation des F90 en Objet ===================
#============================================================
FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(OBJDIR)

.SUFFIXES: .f90 .o
.f90.o: ; @mkdir -p $(BINDIR) $(OBJDIR) $(DATADIR)
	$(FF) -c $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.mod


#============================================================
#============== Complation en exécutables ===================
#============================================================
all: iso_diff.exe iso_diff_abs.exe rayl_diff.exe hg_diff.exe polyn_diff.exe rayl_RGB.exe rayl_stereo_RGB.exe test_bhmie.exe test_bhmie_rb.exe env_eff.exe mie_diff.exe mie_diff_CloudySky.exe mie_RGB.exe

iso_diff.exe: $(OBJFILES_ISO_DIFF)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_ISO_DIFF)) $(FFLIBS)

iso_diff_abs.exe: $(OBJFILES_ISO_ABS)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_ISO_ABS)) $(FFLIBS)

rayl_diff.exe: $(OBJFILES_RAYL_DIFF)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_RAYL_DIFF)) $(FFLIBS)

hg_diff.exe: $(OBJFILES_HG_DIFF)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_HG_DIFF)) $(FFLIBS)

polyn_diff.exe: $(OBJFILES_POLYN_DIFF)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_POLYN_DIFF)) $(FFLIBS)

rayl_RGB.exe: $(OBJFILES_RAYL_RGB)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_RAYL_RGB)) $(FFLIBS)

rayl_stereo_RGB.exe: $(OBJFILES_RAYL_STEREO_RGB)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_RAYL_STEREO_RGB)) $(FFLIBS)

test_bhmie.exe: $(OBJFILES_TEST_BHMIE)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_TEST_BHMIE)) $(FFLIBS)

test_bhmie_rb.exe: $(OBJFILES_TEST_BHMIE_RB)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_TEST_BHMIE_RB)) $(FFLIBS)

env_eff.exe: $(OBJFILES_ENV_EFF)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_ENV_EFF)) $(FFLIBS)

mie_diff.exe: $(OBJFILES_MIE_DIFF)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_MIE_DIFF)) $(FFLIBS)

mie_diff_CloudySky.exe: $(OBJFILES_MIE_DIFF_CLOUDYSKY)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_MIE_DIFF_CLOUDYSKY)) $(FFLIBS)

mie_RGB.exe: $(OBJFILES_MIE_RGB)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_MIE_RGB)) $(FFLIBS)

mie_stereo_RGB.exe: $(OBJFILES_MIE_STEREO_RGB)
	$(FF) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $(OBJFILES_MIE_STEREO_RGB)) $(FFLIBS)


#============================================================
#===================== Nettoyage ============================
#============================================================
.PHONEY: movedata
movedata:
	@mv *.dat $(DATADIR)

.PHONEY: deepclean
deepclean:
	@rm -rf *~ $(OBJDIR) $(BINDIR) $(DATADIR)

.PHONEY: clean
clean:
	@rm -rf *~ $(OBJDIR) $(BINDIR) 


#============================================================
#=============== Dépendances des modules ====================
#============================================================
mod_random.o: mod_param.o

mod_fzero.o: mod_param.o

mod_proj_stereo.o: mod_param.o

mod_iso.o: mod_param.o mod_random.o

mod_diff.o: mod_param.o mod_random.o

mod_mie.o: mod_param.o mod_random.o mie.o mod_fzero.o mod_diff.o

mod_enveloppe.o: mod_param.o mod_random.o mie.o mod_diff.o mod_proj_stereo.o


#============================================================
#============= Dépendances des programmes ===================
#============================================================
iso_diff.o: mod_param.o mod_random.o mod_diff.o mod_iso.o

iso_diff_abs.o: mod_param.o mod_random.o mod_iso.o

rayl_diff.o: mod_param.o mod_random.o mod_diff.o mod_proj_stereo.o

hg_diff.o: mod_param.o mod_random.o mod_diff.o

polyn_diff.o: mod_param.o mod_random.o mod_diff.o

rayl_RGB.o: mod_param.o mod_random.o mod_diff.o

rayl_stereo_RGB.o: mod_param.o mod_proj_stereo.o

test_bhmie.o: mod_param.o mie.o mod_diff.o

test_bhmie_rb.o: mod_param.o mie.o mod_proj_stereo.o

env_eff.o: mod_param.o mod_random.o mie.o mod_proj_stereo.o mod_enveloppe.o mod_diff.o

mie_diff.o: mod_param.o mod_random.o mod_proj_stereo.o mod_fzero.o mie.o mod_diff.o mod_enveloppe.o mod_mie.o

mie_diff_CloudySky.o: mod_param.o mod_random.o mod_proj_stereo.o mod_fzero.o mie.o mod_diff.o mod_enveloppe.o mod_mie.o

mie_RGB.o: mod_param.o mod_random.o mod_proj_stereo.o mod_fzero.o mie.o mod_diff.o mod_enveloppe.o mod_mie.o

mie_stereo_RGB.o: mod_param.o mod_proj_stereo.o
